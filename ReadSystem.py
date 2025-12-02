import os
def readsystem(foamdictionary,path):
    for filedir in os.listdir(path):
        if 'system' in filedir:
            print('Reading OpenFoam System file...')
            for file in os.listdir(path+'/'+filedir):
                if 'controlDict' in file:
                    with open(path+'/'+filedir+'/'+file,'r') as f:
                        lst = f.readlines()
                        for i,ele in enumerate(lst):
                            if not ele.strip():
                                lst=lst[i+1:]
                                break
                        for ele in lst:
                            if ele.strip():
                                if 'application' == ele.strip().split()[0]:
                                    foamdictionary.controlDict.application=ele.strip().strip(';').split()[1]
                                if 'startFrom' ==ele.strip().split()[0]:
                                    foamdictionary.controlDict.startFrom = ele.strip().strip(';').split()[1]
                                if 'startTime'==ele.strip().split()[0]:
                                    foamdictionary.controlDict.startTime =ele.strip().strip(';').split()[1]
                                if 'stopAt'==ele.strip().split()[0]:
                                    foamdictionary.controlDict.stopAt=ele.strip().strip(';').split()[1]
                                if 'endTime'==ele.strip().split()[0]:
                                    foamdictionary.controlDict.endTime = ele.strip().strip(';').split()[1]
                                if 'deltaT' ==ele.strip().split()[0]:
                                    foamdictionary.controlDict.deltaT = ele.strip().strip(';').split()[1]
                                if 'writeControl' ==ele.strip().split()[0]:
                                    foamdictionary.controlDict.writeControl=ele.strip().strip(';').split()[1]
                                if 'writeInterval' ==ele.strip().split()[0]:
                                    foamdictionary.controlDict.writeInterval = ele.strip().strip(';').split()[1]
                                if 'purgeWrite' ==ele.strip().split()[0]:
                                    foamdictionary.controlDict.purgeWrite = ele.strip().strip(';').split()[1]

                if 'fvSchemes' in file:
                    with open(path+'/'+filedir+'/'+file,'r') as f:
                        lst = f.readlines()
                        for i,ele in enumerate(lst):
                            if not ele.strip():
                                lst=lst[i+1:]
                                break
                        for i in range(len(lst)):
                            lst[i]=lst[i].strip().strip('{};')
                        lst=[x for x in lst if x!='']
                        for i in range(len(lst)):
                            if 'ddtSchemes' in lst[i]:
                                foamdictionary.fvSchemes.ddtSchemes.default=' '.join(lst[i+1].split()[1:])
                            if 'gradSchemes' in lst[i]:
                                foamdictionary.fvSchemes.gradSchemes.default=' '.join(lst[i+1].split()[1:])
                            if 'divSchemes' in lst[i]:
                                foamdictionary.fvSchemes.divSchemes.default=' '.join(lst[i+1].split()[1:])
                            if 'laplacianSchemes' in lst[i]:
                                foamdictionary.fvSchemes.laplacianSchemes.default=' '.join(lst[i+1].split()[1:])
                            if 'interpolationSchemes' in lst[i]:
                                foamdictionary.fvSchemes.interpolationSchemes.default=' '.join(lst[i+1].split()[1:])
                            if 'snGradSchemes' in lst[i]:
                                foamdictionary.fvSchemes.snGradSchemes.default=' '.join(lst[i+1].split()[1:])

                if 'fvSolution' in file:
                    with open(path + '/' + filedir + '/' + file, 'r') as f:
                        lst = f.readlines()
                        for i,ele in enumerate(lst):
                            if not ele.strip():
                                lst=lst[i+1:]
                                break
                        for i in range(len(lst)):
                            lst[i]=lst[i].strip().strip('{};')
                        lst=[x for x in lst if x!='']
                        for i in range(len(lst)):
                            if 'solvers' == lst[i]:
                                temp_sov=foamdictionary.fvSolution.solvers
                            if 'SIMPLE' == lst[i]:
                                temp_sim=foamdictionary.fvSolution.SIMPLE
                            if 'relaxationFactors' == lst[i]:
                                temp_rex=foamdictionary.fvSolution.relaxationFactors

                            #solvers
                            if 'p' == lst[i]:
                                temp=temp_sov.p
                            if 'U' == lst[i]:
                                temp=temp_sov.U
                            if 'T' == lst[i]:
                                temp=temp_sov.T

                            if lst[i].split()[0]=='solver':
                                temp.solver=lst[i].split()[1]
                            if lst[i].split()[0]=='preconditioner':
                                temp.preconditioner=lst[i].split()[1]
                            if lst[i].split()[0]=='smoother':
                                temp.smoother=lst[i].split()[1]
                            if lst[i].split()[0]=='tolerance':
                                temp.tolerance=float(lst[i].split()[1])
                            if lst[i].split()[0]=='relTol':
                                temp.relTol=float(lst[i].split()[1])
                            if lst[i].split()[0]=='nPreSweeps':
                                temp.nPreSweeps=int(lst[i].split()[1])
                            if lst[i].split()[0]=='nPreSweeps':
                                temp.nPreSweeps=int(lst[i].split()[1])
                            if lst[i].split()[0]=='nPostSweeps':
                                temp.nPostSweeps=int(lst[i].split()[1])
                            if lst[i].split()[0]=='maxIter':
                                temp.maxIter=int(lst[i].split()[1])
                            if lst[i].split()[0] == 'nFinestSweeps':
                                temp.nFinestSweeps = int(lst[i].split()[1])

                            #SIMPLE
                            if lst[i].split()[0] == 'nCorrectors':
                                temp_sim.nCorrectors = int(lst[i].split()[1])
                            if lst[i].split()[0] == 'pRefCell':
                                temp_sim.pRefCell = int(lst[i].split()[1])
                            if lst[i].split()[0] == 'pRefValue':
                                temp_sim.pRefValue = int(lst[i].split()[1])
                            if lst[i].split()[0] == 'residualControl':
                                temp=temp_sim.residualControl

                            if lst[i].split()[0] =='p' and len(lst[i].split())==2:
                                temp.p=float(lst[i].split()[1])
                            if lst[i].split()[0] =='U' and len(lst[i].split())==2:
                                temp.U=float(lst[i].split()[1])
                            if lst[i].split()[0] =='T' and len(lst[i].split())==2:
                                temp.T=float(lst[i].split()[1])

                            #relaxationFactors
                            if lst[i] == 'fields':
                                temp=temp_rex.fields
                            if lst[i] == 'equations':
                                temp=temp_rex.equations

                    #不存在取默认值
                    #solvers
                    for fieldName in dir(foamdictionary.fvSolution.solvers):
                        if not '__' in fieldName:
                            if getattr(foamdictionary.fvSolution.solvers,fieldName).maxIter==None:
                                getattr(foamdictionary.fvSolution.solvers, fieldName).maxIter=20

                            if getattr(foamdictionary.fvSolution.solvers,fieldName).solver=='GAMG':
                                if getattr(foamdictionary.fvSolution.solvers, fieldName).nPreSweeps==None:
                                    getattr(foamdictionary.fvSolution.solvers, fieldName).nPreSweeps=0
                                if getattr(foamdictionary.fvSolution.solvers, fieldName).nPostSweeps==None:
                                    getattr(foamdictionary.fvSolution.solvers, fieldName).nPostSweeps=2
                                if getattr(foamdictionary.fvSolution.solvers, fieldName).nFinestSweeps==None:
                                    getattr(foamdictionary.fvSolution.solvers, fieldName).nFinestSweeps=2

                    #SIMPLE
                    for fieldName in dir(foamdictionary.fvSolution.SIMPLE.residualControl):
                        if not '__' in fieldName:
                            if not fieldName in dir(foamdictionary.fvSolution.solvers):
                                setattr(foamdictionary.fvSolution.SIMPLE.residualControl,fieldName,1e-6)
                    if foamdictionary.fvSolution.SIMPLE.pRefCell==None:
                        foamdictionary.fvSolution.SIMPLE.pRefCell=0
                    if foamdictionary.fvSolution.SIMPLE.pRefValue == None:
                        foamdictionary.fvSolution.SIMPLE.pRefValue = 0

                    #Relaxation Factor
                    for fieldName in dir(foamdictionary.fvSolution.relaxationFactors.fields):
                        if not '__' in fieldName:
                            if not fieldName in dir(foamdictionary.fvSolution.solvers):
                                setattr(foamdictionary.fvSolution.relaxationFactors.fields,fieldName,1.0)
                    for fieldName in dir(foamdictionary.fvSolution.relaxationFactors.equations):
                        if not '__' in fieldName:
                            if not fieldName in dir(foamdictionary.fvSolution.solvers):
                                setattr(foamdictionary.fvSolution.relaxationFactors.equations,fieldName,1.0)




