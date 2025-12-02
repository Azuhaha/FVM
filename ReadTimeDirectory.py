import os
import numpy as np
from FVM.Class.FluidClass import BoundaryPatchRef
def readtimedirectory(reg):
    print('Reading OpenFoam time file...')
    if reg.foamDictionary.controlDict.startFrom=='startTime':
        timeDirectory=reg.foamDictionary.controlDict.startTime
    elif reg.foamDictionary.controlDict.startFrom=='firstTime':
        timeDirectory='0'
    elif reg.foamDictionary.controlDict.startFrom=='lastestTime':
        raise Exception("To be improved1")
    else:
        raise Exception("\nError in controlDict. startFrom is not valid\n")

    for file in os.listdir(reg.caseDirectoryPath+'/'+timeDirectory):
        if os.path.isfile(reg.caseDirectoryPath+'/'+timeDirectory+'/'+file):

            with open(reg.caseDirectoryPath+'/'+timeDirectory+'/'+file,'r') as f:
                lst=f.readlines()
                for i in range(len(lst)):
                    lst[i]=lst[i].strip().strip('{};')
                lst = [x for x in lst if x != '']

                for i in range(len(lst)):
                    if lst[i].split()[0]=='class':
                        fieldtype=lst[i].split()[1]

                    if lst[i].split()[0]=='object':
                        objectname=lst[i].split()[1]
                        fluidfield=getattr(reg.fluid,objectname)

                    if lst[i].split()[0]=='dimensions':
                        str1=' '.join(lst[i].split()[1:])
                        fielddimensions=[int(x) for x in str1.strip('[]').split()]

                    if lst[i].split()[0]=='internalField':
                        if lst[i].split()[1]=='uniform':
                            if len(lst[i].split())==3:
                                fieldint=np.array([float(lst[i].split()[2])])
                            else:
                                str2=' '.join(lst[i].split()[2:])
                                fieldint=np.array([float(x) for x in str2.strip('()').split()])
                        else:
                            raise Exception("To be improved2")#nonimiform
                    if 'boundaryField' ==lst[i]:
                        idx=-1
                        j=i+1
                        while len(lst[j].split())==1:
                            idx+=1
                            bpf=BoundaryPatchRef()
                            bpf.name=lst[j]
                            j+=1
                            while len(lst[j].split())>1 and '*' not in lst[j]:
                                if 'type' in lst[j].split()[0]:
                                    bpf.type=lst[j].split()[1]
                                if 'value' in lst[j].split()[0]:
                                    if lst[j].split()[1]=='uniform':
                                        bpf.valueType='uniform'
                                        if len(lst[j].split())==3:
                                            bpf.value=np.array([float(lst[j].split()[2])])
                                        else:
                                            str3 = ' '.join(lst[j].split()[2:])
                                            bpf.value = np.array([float(x) for x in str3.strip('()').split()])
                                    else:
                                        raise Exception("To be improved3")  # nonimiform
                                j+=1
                            #若无给定值则默认边界处的值为0
                            if not isinstance(bpf.value, np.ndarray):
                                if bpf.value==None:
                                    if 'Vector' in fieldtype:
                                        bpf.value=np.array([0.,0.,0.])
                                    else:
                                        bpf.value = np.array([0.])
                            fluidfield.boundaryPatchRef[idx]=bpf


                fluidfield.type=fieldtype
                fluidfield.dimensions=fielddimensions
                if 'vol' in fluidfield.type:
                    fluidfield.phi=np.tile(fieldint,(reg.mesh.numberOfElements,1))
                    lst=[[] for _ in range(reg.mesh.numberOfFaces-reg.mesh.numberOfInteriorFace)]
                    for idx in range(len(fluidfield.boundaryPatchRef)):
                        istart=reg.mesh.numberOfElements+reg.mesh.cfdBoundaryPatchesArray[idx].startFaceIndex-reg.mesh.numberOfInteriorFace
                        iend=istart+reg.mesh.cfdBoundaryPatchesArray[idx].nFaces
                        lst[istart-reg.mesh.numberOfElements:iend-reg.mesh.numberOfElements]=[fluidfield.boundaryPatchRef[idx].value]*(iend-istart)
                    fluidfield.phi=np.append(fluidfield.phi,lst,axis=0)

                else:
                    fluidfield.phi = np.tile(fieldint, (reg.mesh.numberOfInteriorFace, 1))
                    lst = [[] for _ in range(reg.mesh.numberOfFaces - reg.mesh.numberOfInteriorFace)]
                    for idx in range(len(fluidfield.boundaryPatchRef)):
                        istart=reg.mesh.cfdBoundaryPatchesArray[idx].startFaceIndex
                        iend=istart+reg.mesh.cfdBoundaryPatchesArray[idx].nFaces
                        lst[istart - reg.mesh.numberOfInteriorFace:iend - reg.mesh.numberOfInteriorFace] = [fluidfield.boundaryPatchRef[idx].value] * (iend - istart)

                    fluidfield.phi=np.append(fluidfield.phi,lst,axis=0)









