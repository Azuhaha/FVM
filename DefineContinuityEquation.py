import os
import numpy as np
from FVM.InterpolateFromElementsToFaces import *
from FVM.Class.FluidClass import BoundaryPatchRef
def definecontinuityequation(reg):
    if not os.path.exists(os.path.join(reg.caseDirectoryPath,'0','U')):
        raise Exception("0文件夹中p文件不存在")
    else:
        #print('Defining Continuity Equation...')
        #初始化
        reg.model.p.residuals.rmsResidual=np.array([1.])
        reg.model.p.residuals.maxResidual = np.array([1.])
        reg.model.p.residuals.initResidual = np.array([1.])
        reg.model.p.residuals.finalResidual= np.array([1.])

        reg.model.p.terms=['massDivergenceTerm']

        if reg.compressible:
            if reg.STEADY_STATE_RUN:
                reg.model.p.terms.append('FalseTransient')
            else:
                reg.model.p.terms.append('Transient')

        #mdot_f初始化
        #print('Initializing mdot_f,DU1,DU2,DU3,DUT1,DUT2,DUT3,pp...')
        reg.fluid.mdot_f.type='surfaceScalarField'

        reg.fluid.mdot_f.phi=np.zeros((reg.mesh.numberOfInteriorFace+reg.mesh.numberOfBFaces,1))
        reg.fluid.mdot_f.previterphi = np.zeros((reg.mesh.numberOfInteriorFace+reg.mesh.numberOfBFaces,1))
        reg.fluid.mdot_f.prevtimestepphi=np.zeros((reg.mesh.numberOfInteriorFace+reg.mesh.numberOfBFaces,1))
        reg.fluid.mdot_f.dimensions=[0,0,0,0,0,0,0]

        #由U计算mdot_f
        U_phi_f=interpolateFromElementsToFaces('linear',reg.fluid.U.phi,reg)
        rho_phi_f = interpolateFromElementsToFaces('linear', reg.fluid.rho.phi, reg)

        U_Sf=np.sum(U_phi_f * reg.mesh.faceSf, axis=1)
        reg.fluid.mdot_f.phi=rho_phi_f*U_Sf[:,np.newaxis]

        #初始化DU1 DU2 DU3 DUT1 DUT2 DUT3 pp
        for fluidfleidname in ['DU1','DU2','DU3','DUT1','DUT2','DUT3','pp']:
            fluidfleid=getattr(reg.fluid,fluidfleidname)
            fluidfleid.type='volScalarField'

            fluidfleid.phi=np.zeros((reg.mesh.numberOfElements+reg.mesh.numberOfBElements,1))
            fluidfleid.previterphi = np.zeros((reg.mesh.numberOfElements+reg.mesh.numberOfBElements,1))
            fluidfleid.prevtimestepphi = np.zeros((reg.mesh.numberOfElements+reg.mesh.numberOfBElements,1))

            if fluidfleidname!='pp':
                for j in range(reg.mesh.numberOfBoundaryPatches):
                    bpr=BoundaryPatchRef()
                    bpr.type='zeroGradient'
                    fluidfleid.boundaryPatchRef[j] =bpr
            else: #pp
                fluidfleid.boundaryPatchRef=reg.fluid.p.boundaryPatchRef
                for j in range(reg.mesh.numberOfBoundaryPatches):
                    if fluidfleid.boundaryPatchRef[j].type=='fixedValue':
                        fluidfleid.boundaryPatchRef[j].value=np.array([0.])







