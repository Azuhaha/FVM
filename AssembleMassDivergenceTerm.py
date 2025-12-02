import numpy as np
from InterpolateFromElementsToFaces import interpolateFromElementsToInteriorFaces
from InterpolateGradientsFromElementsToInteriorFaces import *

def assembleMassDivergenceTerm(reg):
    ##assemble at interior faces
    local_FluxCf=np.zeros((reg.mesh.numberOfInteriorFace,1))
    local_FluxFf = np.zeros((reg.mesh.numberOfInteriorFace, 1))
    local_FluxVf = np.zeros((reg.mesh.numberOfInteriorFace, 1))

    Sf=reg.mesh.faceSf[:reg.mesh.numberOfInteriorFace,:]
    CF=reg.mesh.faceCf[:reg.mesh.numberOfInteriorFace,:]

    MagCF = np.linalg.norm(CF, axis=1)
    e=CF/MagCF[:,np.newaxis]

    DU1_f=interpolateFromElementsToInteriorFaces('linear',reg.fluid.DU1.phi,reg)
    DU2_f = interpolateFromElementsToInteriorFaces('linear', reg.fluid.DU2.phi, reg)
    DU3_f = interpolateFromElementsToInteriorFaces('linear', reg.fluid.DU3.phi, reg)

    U_bar_f=interpolateFromElementsToInteriorFaces('linear', reg.fluid.U.phi, reg)

    rho_f=interpolateFromElementsToInteriorFaces('linearUpwind', reg.fluid.rho.phi, reg)

    p_grad_bar_f=interpolateGradientsFromElementsToInteriorFaces('linear',reg.fluid.p.phiGradient,reg.fluid.p.phi,reg)
    p_grad_f=interpolateGradientsFromElementsToInteriorFaces('Gauss linear corrected',reg.fluid.p.phiGradient,reg.fluid.p.phi,reg)

    #assemble coefficient
    #term1
    U_bar_f=np.sum(U_bar_f*Sf,axis=1)
    local_FluxVf+=rho_f*U_bar_f[:,np.newaxis]

    #term2
    DUSf=np.hstack((DU1_f,DU2_f,DU3_f))*Sf
    MagDUSf=np.linalg.norm(DUSf,axis=1)
    eDUSf=DUSf/MagDUSf[:,np.newaxis]

    eDUSfe=np.sum(eDUSf*e,axis=1)
    DUEf = e*MagDUSf[:,np.newaxis]/eDUSfe[:,np.newaxis]
    geoDiff=np.linalg.norm(DUEf,axis=1)/MagCF

    DUTf=DUSf-DUEf

    local_FluxCf+=rho_f*geoDiff[:,np.newaxis]
    local_FluxFf -= rho_f * geoDiff[:, np.newaxis]

    local_FluxVf -= rho_f * np.sum(p_grad_f[:,:,0]*DUTf,axis=1)[:,np.newaxis]

    local_FluxVf += rho_f * np.sum(p_grad_bar_f[:,:,0] * DUSf, axis=1)[:, np.newaxis]

    #update global fluxes
    reg.fluxes.FluxCf=local_FluxCf
    reg.fluxes.FluxFf = local_FluxFf
    reg.fluxes.FluxVf = local_FluxVf



    ##assemble at boundary faces
    reg.fluxes.FluxCf=np.vstack((reg.fluxes.FluxCf,np.zeros((reg.mesh.numberOfBFaces,1))))
    reg.fluxes.FluxFf = np.vstack((reg.fluxes.FluxFf, np.zeros((reg.mesh.numberOfBFaces, 1))))
    reg.fluxes.FluxVf = np.vstack((reg.fluxes.FluxVf, np.zeros((reg.mesh.numberOfBFaces, 1))))

    for iBPatch in range(reg.mesh.numberOfBoundaryPatches):
        thePhysicalType=reg.mesh.cfdBoundaryPatchesArray[iBPatch].type
        theBCType=reg.fluid.p.boundaryPatchRef[iBPatch].type

        if thePhysicalType=='wall':
            if theBCType=='noSlip':
                assembleMassDivergenceTermWallNoslipBC(iBPatch,reg)
            elif theBCType=='slip':
                assembleMassDivergenceTermWallSlipBC(iBPatch,reg)
            elif theBCType=='zeroGradient':
                assembleMassDivergenceTermWallZeroGradientBC(iBPatch,reg)
            else:
                raise Exception('the BC Type not implemented')
        elif thePhysicalType=='inlet':
            if theBCType=='inlet' or theBCType=='zeroGradient':
                assembleMassDivergenceTermInletZeroGradientBC(iBPatch,reg)
            elif theBCType=='fixedValue':
                assembleMassDivergenceTermInletFixedValueBC(iBPatch,reg)
            else:
                raise Exception('the BC Type not implemented')
        elif thePhysicalType=='outlet':
            if theBCType=='outlet' or theBCType=='zeroGradient':
                assembleMassDivergenceTermOutletZeroGradientBC(iBPatch,reg)
            elif theBCType=='fixedValue':
                assembleMassDivergenceTermOutletFixedValueBC(iBPatch,reg)
            else:
                raise Exception('the BC Type not implemented')
        elif thePhysicalType in ['empty','symmetry','symmetryPlane']:
            continue
        else:
            raise Exception('the Physical Type not implemented')


def assembleMassDivergenceTermWallNoslipBC(iBPatch,reg):
    pass

def assembleMassDivergenceTermWallSlipBC(iBPatch,reg):
    pass

def assembleMassDivergenceTermWallZeroGradientBC(iBPatch,reg):
    pass

def assembleMassDivergenceTermInletZeroGradientBC(iBPatch,reg):
    local_FluxCb=np.zeros((reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces,1))
    local_FluxFb = np.zeros((reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces, 1))

    startBFace=reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    endBFace=startBFace+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    startBElement=reg.mesh.numberOfElements+startBFace-reg.mesh.numberOfInteriorFace
    endBElement=startBElement+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    Sf=reg.mesh.faceSf[startBFace:endBFace,:]

    U_b=reg.fluid.U.phi[startBElement:endBElement,:]
    rho_b=reg.fluid.rho.phi[startBElement:endBElement,:]

    #step1 assemble rhie-chow interpolation term1
    U_b=np.sum(U_b*Sf,axis=1)
    local_FluxVb=rho_b*U_b[:,np.newaxis]

    #update global fluxes
    reg.fluxes.FluxCf[startBFace:endBFace,:]=local_FluxCb
    reg.fluxes.FluxFf[startBFace:endBFace, :] = local_FluxFb
    reg.fluxes.FluxVf[startBFace:endBFace, :] = local_FluxVb


def assembleMassDivergenceTermOutletZeroGradientBC(iBPatch, reg):
    local_FluxCb = np.zeros((reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces, 1))
    local_FluxFb = np.zeros((reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces, 1))
    local_FluxVb = np.zeros((reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces, 1))

    startBFace = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    endBFace = startBFace + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    startBElement = reg.mesh.numberOfElements + startBFace - reg.mesh.numberOfInteriorFace
    endBElement = startBElement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    Sf = reg.mesh.faceSf[startBFace:endBFace, :]

    U_b = reg.fluid.U.phi[startBElement:endBElement, :]
    rho_b = reg.fluid.rho.phi[startBElement:endBElement, :]

    # step1 assemble rhie-chow interpolation term1
    U_b = np.sum(U_b * Sf, axis=1)
    local_FluxVb += rho_b * U_b[:, np.newaxis]

    # update global fluxes
    reg.fluxes.FluxCf[startBFace:endBFace, :] = local_FluxCb
    reg.fluxes.FluxFf[startBFace:endBFace, :] = local_FluxFb
    reg.fluxes.FluxVf[startBFace:endBFace, :] = local_FluxVb

def assembleMassDivergenceTermOutletFixedValueBC(iBPatch, reg):

    local_FluxVb = np.zeros((reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces, 1))

    startBFace = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    endBFace = startBFace + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    startBElement = reg.mesh.numberOfElements + startBFace - reg.mesh.numberOfInteriorFace
    endBElement = startBElement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    Sf = reg.mesh.faceSf[startBFace:endBFace, :]
    CF= reg.mesh.faceCF[startBFace:endBFace, :]
    magCF=np.linalg.norm(CF,axis=1)
    e=CF/magCF[:,np.newaxis]


    U_b = reg.fluid.U.phi[startBElement:endBElement, :]
    rho_b = reg.fluid.rho.phi[startBElement:endBElement, :]
    p_grad_b=reg.fluid.p.phiGradient[startBElement:endBElement, :]

    #assemble coefficients
    DU1_b=reg.fluid.DU1.phi[startBElement:endBElement, :]
    DU2_b = reg.fluid.DU2.phi[startBElement:endBElement, :]
    DU3_b = reg.fluid.DU3.phi[startBElement:endBElement, :]

    DUSb=np.hstack((DU1_b,DU2_b,DU3_b))*Sf
    magDUSb=np.linalg.norm(DUSb,axis=1)
    eDUSb=DUSb/magDUSb[:,np.newaxis]

    eDUSbe=np.sum(e*eDUSb,axis=1)
    DUEb=e*magDUSb[:,np.newaxis]/eDUSbe[:,np.newaxis]
    geoDiff=np.linalg.norm(DUEb,axis=1)/np.linalg.norm(CF,axis=1)


    #assemble term1
    U_bar_b = np.sum(U_b * Sf, axis=1)
    local_FluxVb += rho_b * U_bar_b[:, np.newaxis]
    #assemble term 2
    local_FluxCb=rho_b*geoDiff[:,np.newaxis]
    local_FluxFb = -1*rho_b * geoDiff[:, np.newaxis]

    #assemble term3
    local_FluxVb+=rho_b*np.sum(p_grad_b[:,:,0]*DUSb,axis=1)[:,np.newaxis]

    # update global fluxes
    reg.fluxes.FluxCf[startBFace:endBFace, :] = local_FluxCb
    reg.fluxes.FluxFf[startBFace:endBFace, :] = local_FluxFb
    reg.fluxes.FluxVf[startBFace:endBFace, :] = local_FluxVb


def assembleMassDivergenceTermInletFixedValueBC(iBPatch, reg):
    local_FluxCb = np.zeros((reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces, 1))
    local_FluxFb = np.zeros((reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces, 1))
    local_FluxVb = np.zeros((reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces, 1))


    startBFace = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    endBFace = startBFace + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    startBElement = reg.mesh.numberOfElements + startBFace - reg.mesh.numberOfInteriorFace
    endBElement = startBElement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    Sf = reg.mesh.faceSf[startBFace:endBFace, :]
    CF = reg.mesh.faceCF[startBFace:endBFace, :]
    magCF = np.linalg.norm(CF, axis=1)
    e = CF / magCF[:, np.newaxis]

    U_b = reg.fluid.U.phi[startBElement:endBElement, :]
    rho_b = reg.fluid.rho.phi[startBElement:endBElement, :]
    p_grad_b = reg.fluid.p.phiGradient[startBElement:endBElement, :]

    # assemble coefficients
    DU1_b = reg.fluid.DU1.phi[startBElement:endBElement, :]
    DU2_b = reg.fluid.DU2.phi[startBElement:endBElement, :]
    DU3_b = reg.fluid.DU3.phi[startBElement:endBElement, :]

    DUSb = np.hstack((DU1_b, DU2_b, DU3_b)) * Sf
    magDUSb = np.linalg.norm(DUSb, axis=1)
    eDUSb = DUSb / magDUSb[:, np.newaxis]

    eDUSbe = np.sum(e * eDUSb, axis=1)
    DUEb = e * magDUSb[:, np.newaxis] / eDUSbe[:, np.newaxis]
    geoDiff = np.linalg.norm(DUEb, axis=1) / np.linalg.norm(CF, axis=1)

    # assemble term1
    U_bar_b = np.sum(U_b * Sf, axis=1)
    local_FluxVb += rho_b * U_bar_b[:, np.newaxis]
    # assemble term 2
    local_FluxCb += rho_b * geoDiff[:, np.newaxis]

    # assemble term3
    local_FluxVb += rho_b * np.sum(p_grad_b * DUSb, axis=1)[:, np.newaxis]

    # update global fluxes
    reg.fluxes.FluxCf[startBFace:endBFace, :] = local_FluxCb
    reg.fluxes.FluxFf[startBFace:endBFace, :] = local_FluxFb
    reg.fluxes.FluxVf[startBFace:endBFace, :] = local_FluxVb