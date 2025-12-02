import numpy as np
from InterpolateGradientsFromElementsToInteriorFaces import *
from InterpolateFromElementsToFaces import *
from AssembleMomentumTransientTerm import *
from AssembleMomentumConvectionTerm import *
from AssembleMomentumStressTerm import *
from AssembleMomentumPressureGradientTerm import *
from AssembleMomentumFalseTransientTerm import *
from AssembleMomentumGravitationalForceTerm  import *

from ZeroElementFLUXCoefficients import *
from ZeroFaceFLUXCoefficients import *
from AssembleIntoGlobalMatrixFaceFluxes import *
from AssembleIntoGlobalMatrixElementFluxes import *

def assembleMomentumEquationTerms(iComponent, reg):
    equation=getattr(reg.model,'U')

    for termname in equation.terms:
        if termname=='Transient':
            zeroElementFLUXCoefficients(reg)
            assembleMomentumTransientTerm(iComponent,reg)
            assembleIntoGlobalMatrixElementFluxes(reg)
        elif termname=='Convection':
            zeroFaceFLUXCoefficients(reg)
            assembleMomentumConvectionTerm(iComponent,reg)
            assembleMomentumDCSchemeTerm(iComponent,reg) #DC延迟修正
            assembleIntoGlobalMatrixFaceFluxes(reg)
            zeroElementFLUXCoefficients(reg)
            assembleMomentumDivergenceCorrectionTerm(iComponent,reg) #单元通量修正，未理解？
            assembleIntoGlobalMatrixElementFluxes(reg)
        elif termname == 'Stress':
            zeroFaceFLUXCoefficients(reg)
            assembleMomentumStressTerm(iComponent,reg)
            assembleIntoGlobalMatrixFaceFluxes(reg)
        elif termname == 'PressureGradient':
            zeroElementFLUXCoefficients(reg)
            assembleMomentumPressureGradientTerm(iComponent,reg)
            assembleIntoGlobalMatrixElementFluxes(reg)
        elif termname == 'FalseTransient':
            zeroElementFLUXCoefficients(reg)
            assembleMomentumFalseTransientTerm(iComponent,reg)
            assembleIntoGlobalMatrixElementFluxes(reg)
        elif termname == 'GravitationalForce':
            zeroElementFLUXCoefficients(reg)
            assembleMomentumGravitationalForceTerm(iComponent,reg)
            assembleIntoGlobalMatrixElementFluxes(reg)
        else:
            raise Exception("term is not defined")

def assembleMomentumDCSchemeTerm(iComponent,reg):
    if reg.foamDictionary.fvSchemes.divSchemes=='Gauss upwind':
        return
    # 二阶迎风
    elif reg.foamDictionary.fvSchemes.divSchemes=='Gauss linear':
        #确定面的上游单元
        pos=np.where(reg.fluid.mdot_f.phi[0:reg.mesh.numberOfInteriorFace,:]>0,1,0)
        iupwind=pos*np.array(reg.mesh.owners[0:reg.mesh.numberOfInteriorFace])+(1-pos)*np.array(reg.mesh.neighbours)
        #内部面上游的phi梯度值
        phiGradC=reg.fluid.U.phiGradient[iupwind,:,iComponent]
        phiGrad_f = interpolateGradientsFromElementsToInteriorFaces('Gauss linear corrected', reg.fluid.U.phiGradient[:reg.mesh.numberOfElements,:,iComponent],reg.fluid.U.phi[:reg.mesh.numberOfElements,iComponent],reg)

        Cf=reg.mesh.faceCentroid[:reg.mesh.numberOfInteriorFace,:]-reg.mesh.elementCentroids[iupwind,:]
        #mdot_f*(2*Ugrad_upwind-Ugrad_f)*Cf
        dc_corr=reg.fluid.mdot_f.phi[0:reg.mesh.numberOfInteriorFace,:]*np.sum((2*phiGradC-phiGrad_f[:,:,iComponent])*Cf,axis=1)[:,np.newaxis]
        reg.fluxes.FluxTf[0:reg.mesh.numberOfInteriorFace]+=dc_corr

#避免散度误差
def assembleMomentumDivergenceCorrectionTerm(iComponent,reg):
    effdiv=np.zeros((reg.mesh.numberOfElements,1))
    for iFace in range(reg.mesh.numberOfInteriorFace):
        own=reg.mesh.owners[iFace]
        nei=reg.mesh.neighbours[iFace]

        effdiv[own,:]+=reg.fluid.mdot_f.phi[iFace,:]
        effdiv[nei,:]-=reg.fluid.mdot_f.phi[iFace,:]

    for iBFace in range(reg.mesh.numberOfInteriorFace,reg.mesh.numberOfFaces):
        own=reg.mesh.owners[iBFace]
        effdiv[own,:]+=reg.fluid.mdot_f.phi[iBFace,:]

    reg.fluxes.FluxC=np.maximum(effdiv,0.)-effdiv
    U_phi = reg.fluid.U.phi[:reg.mesh.numberOfElements, iComponent]
    reg.fluxes.FluxV=-1*np.maximum(effdiv,0.)*U_phi[:,np.newaxis]
    reg.fluxes.FluxT=reg.fluxes.FluxC*U_phi[:,np.newaxis]+reg.fluxes.FluxV