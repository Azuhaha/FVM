import os
import numpy as np
from Class.FluidClass import BoundaryPatchRef
from UpdateScales import updatescale

def readtransportproperties(reg):
    if os.path.exists(reg.caseDirectoryPath + '/constant/transportProperties'):
        print('Reading OpenFoam TransportProperties file...')
        with open(reg.caseDirectoryPath + '/constant/transportProperties','r' ) as f:
            lst = f.readlines()
            for i, ele in enumerate(lst):
                if not ele.strip():
                    lst = lst[i + 1:]
                    break
            for i, ele in enumerate(lst):
                lst[i]=ele.strip().strip(';/*').strip().strip('*')

            lst = [x for x in lst if x != '']
            for i, ele in enumerate(lst):
                name=ele.split()[0]
                val=np.array([float(ele.split()[-1])])
                temp=' '.join(ele.split()[1:-1])
                dimension=[int(x) for x in temp.strip('[]').split()]
                fluidfield=getattr(reg.fluid,name)

                fluidfield.phi=np.tile(val,(reg.mesh.numberOfElements+reg.mesh.numberOfBElements,1))
                fluidfield.previterphi=np.tile(np.array([0.]),(reg.mesh.numberOfElements+reg.mesh.numberOfBElements,1))
                fluidfield.prevtimestepphi=np.tile(np.array([0.]),(reg.mesh.numberOfElements+reg.mesh.numberOfBElements,1))

                fluidfield.dimensions=dimension

                for iBP in range(reg.mesh.numberOfBoundaryPatches):
                    bpr=BoundaryPatchRef()
                    bpr.type='zeroGradient'
                    bpr.value=val
                    fluidfield.boundaryPatchRef[iBP]=bpr

                transportpropertyval=getattr(reg.foamDictionary.transportProperties,name)
                transportpropertyval.name=name
                transportpropertyval.dimensions=dimension
                transportpropertyval.propertyValue=val

                updatescale(fluidfield,reg)

                if name=='rho':
                    reg.compressible=False


#后续完善

