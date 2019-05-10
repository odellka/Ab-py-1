"""General 3D topology optimization code by Zhihao Zuo and Yimin Xie. Note that the CAE file shall contain a model 'Model-1' with a dependent part 'Part-1' and a static step 'Step-1'."""
import math,customKernel,random
from abaqus import getInput,getInputs
from odbAccess import openOdb
## Function of formatting Abaqus model for stiffness optimisation
def fmtMdb(Mdb):
    mdl = Mdb.models['Model-1']
    part = mdl.parts['Part-1']
    # Build sections and assign solid section
    mdl.Material('Material01').Elastic(((1.0, 0.3), ))
    mdl.HomogeneousSolidSection('sldSec','Material01')
    mdl.Material('Material02').Elastic(((0.001**3, 0.3), ))
    mdl.HomogeneousSolidSection('voidSec','Material02')
    part.SectionAssignment(part.Set('ss',part.elements),'sldSec')
    # Define output request
    mdl.FieldOutputRequest('SEDensity','Step-1',variables=('ELEDEN', ))
    mdl.HistoryOutputRequest('ExtWork','Step-1',variables=('ALLWK', ))
## Function of running FEA for raw sensitivities and objective function
def FEA(Iter,Mdb,Xe,Ae):
    var=79
    Mdb.Job('Design_Job'+str(Iter+var),'Model-1').submit()
    Mdb.jobs['Design_Job'+str(Iter+var)].waitForCompletion()
    #status=Mdb.Job('Design_Job'+str(Iter+var),'Model-1').status
    #print status
    #opdb = openOdb('Design_Job'+str(Iter+var)+'.odb')
    #f1 = open('Design_Job'+str(Iter+var)+".msg", "r")
    #last_line = f1.readlines()[-1]
    #print (last_line)
    #seng = opdb.steps['Step-1'].frames[-1].fieldOutputs['ESEDEN'].values
    #for en in seng: 
        #print(Xe[en.elementLabel])
    #    Ae[en.elementLabel]=en.data/Xe[en.elementLabel]
    elsetName = 'SET-4'
    
    #if 'TERMINATED' in open('Design_Job'+str(Iter+var)+".msg", "r").read():
    #    obj = (1000, random.randint(1,len(elmts)+1))
    #    print ('ERRRRRRROR')
    #else:      
    #    obj=getMaxMises('Design_Job'+str(Iter+var)+'.odb',elsetName)
    obj=getMaxMises('Design_Job'+str(Iter+var)+'.odb',elsetName)    
    #opdb.steps['Step-1'].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLWK'].data[-1][1]
    #opdb.close()
    return obj

def getMaxMises(odbName,elsetName):
    """ Print max mises location and value given odbName
        and elset(optional)
    """
    elset = elemset = None
    region = "over the entire model"
    """ Open the output database """
    odb = openOdb(odbName)
    assembly = odb.rootAssembly

    """ Check to see if the element set exists
        in the assembly
    """
    if elsetName:
        try:
            elemset = assembly.elementSets[elsetName]
            region = " in the element set : " + elsetName;
        except KeyError:
            print 'An assembly level elset named %s does' \
                   'not exist in the output database %s' \
                   % (elsetName, odbName)
            odb.close()
            
    """ Initialize maximum values """
    maxMises = -0.1
    maxElem = 0
    maxStep = "Step-1"
    maxFrame = -1
    Stress = 'S'
    isStressPresent = 0
    for step in odb.steps.values():
        print 'Processing Step:', step.name
        for frame in step.frames:
            allFields = frame.fieldOutputs
            if (allFields.has_key(Stress)):
                isStressPresent = 1
                stressSet = allFields[Stress]
                if elemset:
                    stressSet = stressSet.getSubset(
                        region=elemset)      
                for stressValue in stressSet.values:                
                    if (stressValue.mises > maxMises):
                        maxMises = stressValue.mises
                        maxElem = stressValue.elementLabel
                        maxStep = step.name
                        maxFrame = frame.incrementNumber
    if(isStressPresent):
        print 'Maximum von Mises stress %s is %f in element %d'%(
            region, maxMises, maxElem)
        print 'Location: frame # %d  step:  %s '%(maxFrame,maxStep)
    else:
        print 'Stress output is not available in' \
              'the output database : %s\n' %(odb.name)
    
    """ Close the output database before exiting the program """
    obj=(maxMises, maxElem)
    return obj
    odb.close()
    
    
    ## Function of preparing filter map (Fm={elm1:[[el1,el2,...],[wf1,wf2,...]],...})
def preFlt(Rmin,Elmts,Nds,Fm):
    # Calculate element centre coordinates
    c0 = {}
    for el in Elmts:
        nds = el.connectivity
        c0[el.label]=[sum([Nds[nd].coordinates[i]/len(nds) for nd in nds]) for i in range(3)]
        print (c0[el.label])
    # Weighting factors
    for el in Elmts:
        Fm[el.label] = [[],[]]
        for em in Elmts:
            dis=math.sqrt(sum([(c0[el.label][i]-c0[em.label][i])**2 for i in range(3)]))
            if dis<Rmin:
                Fm[el.label][0].append(em.label)
                Fm[el.label][1].append(Rmin - dis)
        sm = sum(Fm[el.label][1])
        for i in range(len(Fm[el.label][0])): Fm[el.label][1][i] /= sm
## Function of filtering sensitivities
def fltAe(Ae,Fm):
    raw = Ae.copy()
    for el in Fm.keys():
        Ae[el] = 0.0
        for i in range(len(Fm[el][0])): Ae[el]+=raw[Fm[el][0][i]]*Fm[el][1][i]
## Function of optimality update for design variables and Abaqus model
def BESO(Vf,Xe,Ae,Part,Elmts):
    lo, hi = min(Ae.values()), max(Ae.values())
    tv = Vf*len(Elmts)
    while (hi-lo)/hi > e-5:
        th = (lo+hi)/2.0
        for key in Xe.keys(): Xe[key] = 1.0 if Ae[key]>th else 0.001
        if sum(Xe.values())-tv>0: lo = th
        else: hi = th
    # Label elements as solid or void
    vlb, slb = [], []
    for el in Elmts:
        
        if Xe[el.label] == 1.0: slb.append(el.label)
        else: vlb.append(el.label)
    # Assign solid and void elements to each section
    Part.SectionAssignment(Part.SetFromElementLabels('ss',slb),'sldSec')
    Part.SectionAssignment(Part.SetFromElementLabels('vs',vlb),'voidSec')
## ====== MAIN PROGRAM ======
if __name__ == '__main__':
    # Set parameters and inputs
    vf,rmin,ert = (0.5, 1, 0.02)
    mddb = openMdb('2D-2.cae')
    # Design initialization
    fmtMdb(mddb)
    part = mddb.models['Model-1'].parts['Part-1']
    elmts, nds = part.elements, part.nodes
    oh, vh = [], []


    xe, ae, oae, fm = {}, {}, {}, {}
    #for el in elmts: 
    #    xe[el.label] = 1.0
    #    print(0)
    #if rmin>0: preFlt(rmin,elmts,nds,fm)
    
    # Optimisation iteration
    
    min_stress=0;
    change, iter, obj = 1, -1, 0
    while change > 0.001 and iter<30:
        iter += 1
        if iter==0:
            vlb, slb = [], []
            slb=range(1,len(elmts)+1)
            vlb=[81, 20, 19, 18, 17, 16, 15, 14, 13, 102, 12, 11, 10, 9, 181, 182, 183, 184, 8, 185, 7, 6, 186, 200, 5, 199, 198, 187, 4, 197]
            for el in vlb:
                slb.remove(el)
            print slb
            #for el in elmts:
                #slb.append(el.label)
            part.SectionAssignment(part.SetFromElementLabels('ss',slb),'sldSec')
            (stress, maxElem)=FEA(iter,mddb,xe,ae);
            oh.append(stress)
            vlb_best=vlb[:]
            slb_best=slb[:]
            min_stress=stress
        else:
            #vlb=vlb_best[:]
            #slb=slb_best[:]

            eltoremove=maxElem#random.randint(1,len(elmts)+1)
            print(maxElem)
            if eltoremove in slb:
                slb.remove(eltoremove)
                vlb.append(eltoremove)
            elif eltoremove in vlb:
                vlb.remove(eltoremove)
                slb.append(eltoremove)
            print(vlb)

            if len(slb)>0: part.SectionAssignment(part.SetFromElementLabels('ss',slb),'sldSec')
            if len(vlb)>0: part.SectionAssignment(part.SetFromElementLabels('vs',vlb),'voidSec')
            #for el in elmts:
            #    #print(el.label)
            #    if el.label==randel: vlb.append(el.label)
            #    else: slb.append(el.label)
            #xe[el.label] = 1.0
            #Setting section materials

            # Run FEA
            (stress, maxElem)=FEA(iter,mddb,xe,ae);
            oh.append(stress)
            if stress <= min_stress:
                print('NEW RECORD')
                vlb_best=vlb[:]
                slb_best=slb[:]
                min_stress=stress
            print(oh)
        # Process sensitivities
        #if rmin>0: fltAe(ae,fm)
        #for el in Fm.keys():
            #Ae[el] = 0.0
        #if iter > 0: ae=dict([(k,(ae[k]+oae[k])/2.0) for k in ae.keys()])
        #oae = ae.copy()
        # BESO optimisation
        #for key in Xe.keys(): Xe[key] = 1.0 
        

            #else: vlb.append(el.label)
                

        # Assign solid and void elements to each section

        
        #vh.append(sum(xe.values())/len(xe))
        ##nv = max(vf,vh[-1]*(1.0-ert))
        #BESO(nv,xe,ae,part,elmts)
        #if iter>10: change=math.fabs((sum(oh[iter-4:iter+1])-sum(oh[iter-9:iter-4]))/sum(oh[iter-9:iter-4]))
    # Save results
    #mddb.customData.History = {'vol':vh,'obj':oh}
    #mddb.saveAs('Final_design.cae')
