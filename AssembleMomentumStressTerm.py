import numpy as np
from InterpolateGradientsFromElementsToInteriorFaces import *
from InterpolateFromElementsToFaces import *

def assembleMomentumStressTerm(iComponent,reg):
    #/////////////////////////////////////////////////////内部面
    CF=reg.mesh.faceCF[:reg.mesh.numberOfInteriorFace,:]
    Sf = reg.mesh.faceSf[:reg.mesh.numberOfInteriorFace, :]

    n=Sf/np.linalg.norm(Sf,axis=1)[:,np.newaxis]
    e=CF/np.linalg.norm(CF,axis=1)[:,np.newaxis]

    CF_limited=np.maximum(np.sum(CF*n,axis=1),0.05*np.linalg.norm(CF,axis=1))

    mag_Sf=np.linalg.norm(Sf,axis=1)

    #正交分量和非正交分量 正交分量为Sf在CF的分量
    Ef=e*mag_Sf[:,np.newaxis] #正交修正法
    #Ef=e*np.sum(Sf*e,axis=1)[:,np.newaxis] #最小修正法
    #Ef=e*np.sum(Sf*Sf,axis=1)[:,np.newaxis]/np.sum(e*Sf,axis=1)[:,np.newaxis] #超松弛法
    Tf=Sf-Ef

    UGrad_f=interpolateGradientsFromElementsToInteriorFaces('linear',reg.fluid.U.phiGradient[:reg.mesh.numberOfElements,:,:],reg.fluid.U.phi[:reg.mesh.numberOfElements,:],reg)

    #正交分量上的几何比重 Ef
    geoDiff_f=np.linalg.norm(Ef,axis=1)/CF_limited #?

    mu_f=interpolateFromElementsToFaces('linear',reg.fluid.mu.phi,reg)

    #mu_f*gDiff
    reg.fluxes.FluxCf[:reg.mesh.numberOfInteriorFace,:]=mu_f[:reg.mesh.numberOfInteriorFace,:]*geoDiff_f[:,np.newaxis]
    #-mu_f*gDiff
    reg.fluxes.FluxFf[:reg.mesh.numberOfInteriorFace, :] = -1*mu_f[:reg.mesh.numberOfInteriorFace, :] * geoDiff_f[:,np.newaxis]
    #非正交部分 -mu_f*gradU_f*Tf
    reg.fluxes.FluxVf[:reg.mesh.numberOfInteriorFace, :]=-1*mu_f[:reg.mesh.numberOfInteriorFace, :]*np.sum(UGrad_f[:,:,iComponent]*Tf,axis=1)[:,np.newaxis]

    #添加速度梯度的转置项
    UGrad_f_T=UGrad_f.transpose(0, 2, 1)
    #减/加?
    reg.fluxes.FluxVf[:reg.mesh.numberOfInteriorFace, :]+=mu_f[:reg.mesh.numberOfInteriorFace, :]*np.sum(UGrad_f_T[:,:,iComponent]*Sf,axis=1)[:,np.newaxis]

    #添加可压缩的产生项
    if reg.compressible:
        UDiv_Sf=(UGrad_f[:, 0, 0] + UGrad_f[:, 1, 1] + UGrad_f[:, 2, 2]) * Sf[:, iComponent]
        reg.fluxes.FluxVf[:reg.mesh.numberOfInteriorFace, :] += 2./3.*mu_f[:reg.mesh.numberOfInteriorFace, :]*UDiv_Sf[:,np.newaxis]

    #mu_f*gDiff*U_C-mu_f*gDiff*U_F-mu_f*gradU_f*Tf+mu_f*UGrad_f_T*Sf
    reg.fluxes.FluxTf[:reg.mesh.numberOfInteriorFace, :]=(reg.fluxes.FluxCf[:reg.mesh.numberOfInteriorFace,:]*reg.fluid.U.phi[reg.mesh.owners[:reg.mesh.numberOfInteriorFace],iComponent][:,np.newaxis]
                                                          +reg.fluxes.FluxFf[:reg.mesh.numberOfInteriorFace, :]*reg.fluid.U.phi[reg.mesh.neighbours,iComponent][:,np.newaxis]
                                                          +reg.fluxes.FluxVf[:reg.mesh.numberOfInteriorFace, :])

    #/////////////////////////////////////////////////////边界面
    for iBPatch in range(reg.mesh.numberOfBoundaryPatches):
        if reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='wall':
            if reg.fluid.U.boundaryPatchRef[iBPatch].type=='slip':
                continue
            elif reg.fluid.U.boundaryPatchRef[iBPatch].type=='noSlip':
                assembleStressTermWallNoslipBC(iBPatch, iComponent,reg)
            elif reg.fluid.U.boundaryPatchRef[iBPatch].type=='fixedValue':
                assembleStressTermWallFixedValueBC(iBPatch, iComponent, reg)
            else:
                raise Exception(reg.fluid.U.boundaryPatchRef[iBPatch].type+' bc in not implemented')
        elif reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='inlet':
            if reg.fluid.U.boundaryPatchRef[iBPatch].type=='fixedValue':
                assembleStressTermInletFixedValueBC(iBPatch, iComponent, reg)
            elif reg.fluid.U.boundaryPatchRef[iBPatch].type=='zeroGradient' or reg.fluid.U.boundaryPatchRef[iBPatch].type=='inlet':
                assembleStressTermInletZeroGradientBC(iBPatch, iComponent,reg)
            else:
                raise Exception(reg.fluid.U.boundaryPatchRef[iBPatch].type + ' bc in not implemented')
        elif reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='outlet':
            if reg.fluid.U.boundaryPatchRef[iBPatch].type=='zeroGradient' or reg.fluid.U.boundaryPatchRef[iBPatch].type=='outlet':
                assembleStressTermOutletZeroGradientBC(iBPatch, iComponent, reg)
            elif reg.fluid.U.boundaryPatchRef[iBPatch].type=='fixedValue':
                assembleStressTermOutletFixedValueBC(iBPatch, iComponent,reg)
            else:
                raise Exception(reg.fluid.U.boundaryPatchRef[iBPatch].type + ' bc in not implemented')
        elif reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='symmetry' or reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='symmetryPlane':
            assembleStressTermSymmetryBC(iBPatch, iComponent, reg)
        elif reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='empty':
            assembleStressTermEmptyBC(iBPatch, iComponent, reg)
        else:
            raise Exception(reg.mesh.cfdBoundaryPatchesArray[iBPatch].type + ' physical condition is not implemented')


def assembleStressTermWallNoslipBC(iBPatch, iComponent,reg):
    istartface=reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iendface=istartface+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    istartelement=reg.mesh.numberOfElements+istartface-reg.mesh.numberOfInteriorFace
    iendelement=istartelement+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    mu_b=reg.fluid.mu.phi[istartelement:iendelement,:]
    Sf_b=reg.mesh.faceSf[istartface:iendface,:]
    magSf_b=np.linalg.norm(Sf_b,axis=1)
    wallDist_b=reg.mesh.wallDistLimited[istartface:iendface,:]

    n=Sf_b/magSf_b[:,np.newaxis]
    nx=n[:,0]
    ny=n[:,1]
    nz=n[:,2]
    nx2= nx*nx
    ny2 = ny*ny
    nz2 = nz*nz

    U_b=reg.fluid.U.phi[istartelement:iendelement,:]
    U_C=reg.fluid.U.phi[reg.mesh.owners[istartface:iendface],:]

    u_b=U_b[:,0]
    v_b=U_b[:,1]
    w_b=U_b[:,2]

    u_C=U_C[:,0]
    v_C=U_C[:,1]
    w_C=U_C[:,2]

    if iComponent==0:
        reg.fluxes.FluxCf[istartface:iendface,:]=mu_b*magSf_b[:,np.newaxis]/wallDist_b*(1-nx2[:,np.newaxis])
        reg.fluxes.FluxFf[istartface:iendface, :] =np.zeros((iendface-istartface,1))
        reg.fluxes.FluxVf[istartface:iendface, :]=-1*mu_b*magSf_b[:,np.newaxis]/wallDist_b*(u_b*(1-nx2)+(v_C-v_b)*ny*nx+(w_C-w_b)*nz*nx)[:,np.newaxis]
    elif iComponent==1:
        reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * magSf_b[:, np.newaxis] / wallDist_b * (1 - ny2[:, np.newaxis])
        reg.fluxes.FluxFf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
        reg.fluxes.FluxVf[istartface:iendface, :] = -1*mu_b * magSf_b[:, np.newaxis] / wallDist_b * ((u_C-u_b) * nx*ny +v_b* (1-ny2) + (w_C - w_b) * nz * ny)[:, np.newaxis]
    else:
        reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * magSf_b[:, np.newaxis] / wallDist_b * (1 - nz2[:, np.newaxis])
        reg.fluxes.FluxFf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
        reg.fluxes.FluxVf[istartface:iendface, :] = -1*mu_b * magSf_b[:, np.newaxis] / wallDist_b * ((u_C-u_b) * nx*nz +(v_C-v_b)*ny*nz +w_b * (1-nz2))[:, np.newaxis]

    reg.fluxes.FluxTf[istartface:iendface, :] = reg.fluxes.FluxCf[istartface:iendface, :]*U_C[:,[iComponent]]+reg.fluxes.FluxFf[istartface:iendface, :]*U_b[:,[iComponent]]+reg.fluxes.FluxVf[istartface:iendface, :]

#与assembleStressTermWallNoslipBC一致
def assembleStressTermWallFixedValueBC(iBPatch, iComponent, reg):
    istartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iendface = istartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    istartelement = reg.mesh.numberOfElements + istartface - reg.mesh.numberOfInteriorFace
    iendelement = istartelement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    mu_b = reg.fluid.mu.phi[istartelement, iendelement:]
    Sf_b = reg.mesh.faceSf[istartface:iendface, :]
    magSf_b = np.linalg.norm(Sf_b, axis=1)
    wallDist_b = reg.mesh.wallDistLimited[istartface:iendface, :]

    n = Sf_b / magSf_b[:, np.newaxis]
    nx = n[:, 0]
    ny = n[:, 1]
    nz = n[:, 2]
    nx2= nx*nx
    ny2 = ny*ny
    nz2 = nz*nz

    U_b = reg.fluid.U.phi[istartelement:iendelement, :]
    U_C = reg.fluid.U.phi[reg.mesh.owners[istartface:iendface], :]

    u_b = U_b[:, 0]
    v_b = U_b[:, 1]
    w_b = U_b[:, 2]

    u_C = U_C[:, 0]
    v_C = U_C[:, 1]
    w_C = U_C[:, 2]

    if iComponent == 0:
        reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * magSf_b[:, np.newaxis] / wallDist_b * (1 - nx2[:, np.newaxis])
        reg.fluxes.FluxFf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
        reg.fluxes.FluxVf[istartface:iendface, :] = -1 * mu_b * magSf_b[:, np.newaxis] / wallDist_b * (u_b * (1 - nx2) + (v_C - v_b) * ny * nx + (w_C - w_b) * nz * nx)[:, np.newaxis]
    elif iComponent == 1:
        reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * magSf_b[:, np.newaxis] / wallDist_b * (1 - ny2[:, np.newaxis])
        reg.fluxes.FluxFf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
        reg.fluxes.FluxVf[istartface:iendface, :] = -1 * mu_b * magSf_b[:, np.newaxis] / wallDist_b * ((u_C - u_b) * nx * ny + v_b * (1 - ny2) + (w_C - w_b) * nz * ny)[:, np.newaxis]
    else:
        reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * magSf_b[:, np.newaxis] / wallDist_b * (1 - nz2[:, np.newaxis])
        reg.fluxes.FluxFf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
        reg.fluxes.FluxVf[istartface:iendface, :] = -1 * mu_b * magSf_b[:, np.newaxis] / wallDist_b * ((u_C - u_b) * nx * nz + (v_C - v_b) * ny * nz + w_b * (1 - nz2))[:, np.newaxis]

    reg.fluxes.FluxTf[istartface:iendface, :] = reg.fluxes.FluxCf[istartface:iendface, :] * U_C[:, [iComponent]] + reg.fluxes.FluxFf[istartface:iendface,:] * U_b[:,[iComponent]] + reg.fluxes.FluxVf[istartface:iendface,:]

def assembleStressTermInletFixedValueBC(iBPatch, iComponent, reg):
    istartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iendface = istartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    istartelement = reg.mesh.numberOfElements + istartface - reg.mesh.numberOfInteriorFace
    iendelement = istartelement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    mu_b = reg.fluid.mu.phi[istartelement:iendelement,:]
    Sf_b = reg.mesh.faceSf[istartface:iendface, :]
    magSf_b = np.linalg.norm(Sf_b, axis=1)
    wallDist_b = reg.mesh.wallDistLimited[istartface:iendface, :]

    Ui_b = reg.fluid.U.phi[istartelement:iendelement, iComponent]
    Ui_C = reg.fluid.U.phi[reg.mesh.owners[istartface:iendface], iComponent]

    geoDiff_b=magSf_b[:,np.newaxis]/wallDist_b

    reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * geoDiff_b
    reg.fluxes.FluxFf[istartface:iendface, :] = -1*mu_b * geoDiff_b
    reg.fluxes.FluxVf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
    reg.fluxes.FluxTf[istartface:iendface, :] = reg.fluxes.FluxCf[istartface:iendface, :] * Ui_C[:,np.newaxis] + reg.fluxes.FluxFf[istartface:iendface,:] * Ui_b[:,np.newaxis] + reg.fluxes.FluxVf[istartface:iendface,:]

#assembleStressTermInletFixedValueBC一致
def assembleStressTermInletZeroGradientBC(iBPatch, iComponent,reg):
    istartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iendface = istartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    istartelement = reg.mesh.numberOfElements + istartface - reg.mesh.numberOfInteriorFace
    iendelement = istartelement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    mu_b = reg.fluid.mu.phi[istartelement, iendelement:]
    Sf_b = reg.mesh.faceSf[istartface:iendface, :]
    magSf_b = np.linalg.norm(Sf_b, axis=1)
    wallDist_b = reg.mesh.wallDistLimited[istartface:iendface, :]

    Ui_b = reg.fluid.U.phi[istartelement:iendelement, iComponent]
    Ui_C = reg.fluid.U.phi[reg.mesh.owners[istartface:iendface], iComponent]

    geoDiff_b=magSf_b[:,np.newaxis]/wallDist_b

    reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * geoDiff_b
    reg.fluxes.FluxFf[istartface:iendface, :] = -1*mu_b * geoDiff_b
    reg.fluxes.FluxVf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
    reg.fluxes.FluxTf[istartface:iendface, :] = reg.fluxes.FluxCf[istartface:iendface, :] * Ui_C[:,np.newaxis] + reg.fluxes.FluxFf[istartface:iendface,:] * Ui_b[:,np.newaxis] + reg.fluxes.FluxVf[istartface:iendface,:]


def assembleStressTermOutletZeroGradientBC(iBPatch, iComponent, reg):
    istartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iendface = istartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    istartelement = reg.mesh.numberOfElements + istartface - reg.mesh.numberOfInteriorFace
    iendelement = istartelement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    mu_b = reg.fluid.mu.phi[istartelement:iendelement,:]
    Sf_b = reg.mesh.faceSf[istartface:iendface, :]
    magSf_b = np.linalg.norm(Sf_b, axis=1)
    wallDist_b = reg.mesh.wallDistLimited[istartface:iendface, :]

    Ui_b = reg.fluid.U.phi[istartelement:iendelement, iComponent]
    Ui_C = reg.fluid.U.phi[reg.mesh.owners[istartface:iendface], iComponent]

    geoDiff_b=magSf_b[:,np.newaxis]/wallDist_b

    reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * geoDiff_b
    reg.fluxes.FluxFf[istartface:iendface, :] = -1*mu_b * geoDiff_b
    reg.fluxes.FluxVf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
    reg.fluxes.FluxTf[istartface:iendface, :] = reg.fluxes.FluxCf[istartface:iendface, :] * Ui_C[:,np.newaxis] + reg.fluxes.FluxFf[istartface:iendface,:] * Ui_b[:,np.newaxis] + reg.fluxes.FluxVf[istartface:iendface,:]



def assembleStressTermOutletFixedValueBC(iBPatch, iComponent,reg):
    istartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iendface = istartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    istartelement = reg.mesh.numberOfElements + istartface - reg.mesh.numberOfInteriorFace
    iendelement = istartelement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    mu_b = reg.fluid.mu.phi[istartelement:iendelement,:]
    Sf_b = reg.mesh.faceSf[istartface:iendface, :]
    magSf_b = np.linalg.norm(Sf_b, axis=1)
    wallDist_b = reg.mesh.wallDistLimited[istartface:iendface, :]

    Ui_b = reg.fluid.U.phi[istartelement:iendelement, iComponent]
    Ui_C = reg.fluid.U.phi[reg.mesh.owners[istartface:iendface], iComponent]

    geoDiff_b=magSf_b[:,np.newaxis]/wallDist_b

    reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * geoDiff_b
    reg.fluxes.FluxFf[istartface:iendface, :] = -1*mu_b * geoDiff_b
    reg.fluxes.FluxVf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
    reg.fluxes.FluxTf[istartface:iendface, :] = reg.fluxes.FluxCf[istartface:iendface, :] * Ui_C[:,np.newaxis] + reg.fluxes.FluxFf[istartface:iendface,:] * Ui_b[:,np.newaxis] + reg.fluxes.FluxVf[istartface:iendface,:]

def assembleStressTermSymmetryBC(iBPatch, iComponent, reg):
    istartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iendface = istartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    istartelement = reg.mesh.numberOfElements + istartface - reg.mesh.numberOfInteriorFace
    iendelement = istartelement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    mu_b = reg.fluid.mu.phi[istartelement:iendelement,:]
    Sf_b = reg.mesh.faceSf[istartface:iendface, :]
    magSf_b = np.linalg.norm(Sf_b, axis=1)
    wallDist_b = reg.mesh.wallDistLimited[istartface:iendface, :]

    n = Sf_b / magSf_b[:, np.newaxis]
    nx = n[:, 0]
    ny = n[:, 1]
    nz = n[:, 2]
    nx2= nx*nx
    ny2 = ny*ny
    nz2 = nz*nz

    U_b = reg.fluid.U.phi[istartelement:iendelement, :]
    U_C = reg.fluid.U.phi[reg.mesh.owners[istartface:iendface], :]

    u_b = U_b[:, 0]
    v_b = U_b[:, 1]
    w_b = U_b[:, 2]

    u_C = U_C[:, 0]
    v_C = U_C[:, 1]
    w_C = U_C[:, 2]

    if iComponent == 0:
        reg.fluxes.FluxCf[istartface:iendface, :] = 2*mu_b * magSf_b[:, np.newaxis] / wallDist_b * nx2[:, np.newaxis]
        reg.fluxes.FluxFf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
        reg.fluxes.FluxVf[istartface:iendface, :] = 2 * mu_b * magSf_b[:, np.newaxis] / wallDist_b * ((v_C * ny + w_C * nz )*nx)[:, np.newaxis]
    elif iComponent == 1:
        reg.fluxes.FluxCf[istartface:iendface, :] = 2*mu_b * magSf_b[:, np.newaxis] / wallDist_b * ny2[:, np.newaxis]
        reg.fluxes.FluxFf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
        reg.fluxes.FluxVf[istartface:iendface, :] = 2 * mu_b * magSf_b[:, np.newaxis] / wallDist_b * ((u_C * nx+w_C * nz)*ny)[:, np.newaxis]
    else:
        reg.fluxes.FluxCf[istartface:iendface, :] = 2*mu_b * magSf_b[:, np.newaxis] / wallDist_b * nz2[:, np.newaxis]
        reg.fluxes.FluxFf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
        reg.fluxes.FluxVf[istartface:iendface, :] = 2 * mu_b * magSf_b[:, np.newaxis] / wallDist_b * ((u_C * nx + v_C * ny)*nz)[:, np.newaxis]

    reg.fluxes.FluxTf[istartface:iendface, :] = reg.fluxes.FluxCf[istartface:iendface, :] * U_C[:, [iComponent]] + reg.fluxes.FluxFf[istartface:iendface,:] * U_b[:,[iComponent]] + reg.fluxes.FluxVf[istartface:iendface,:]

#与assembleStressTermWallNoslipBC一致
def assembleStressTermEmptyBC(iBPatch, iComponent, reg):
    istartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iendface = istartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    istartelement = reg.mesh.numberOfElements + istartface - reg.mesh.numberOfInteriorFace
    iendelement = istartelement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    mu_b = reg.fluid.mu.phi[istartelement:iendelement,:]
    Sf_b = reg.mesh.faceSf[istartface:iendface, :]
    magSf_b = np.linalg.norm(Sf_b, axis=1)
    wallDist_b = reg.mesh.wallDistLimited[istartface:iendface, :]

    n = Sf_b / magSf_b[:, np.newaxis]
    nx = n[:, 0]
    ny = n[:, 1]
    nz = n[:, 2]
    nx2= nx*nx
    ny2 = ny*ny
    nz2 = nz*nz

    U_b = reg.fluid.U.phi[istartelement:iendelement, :]
    U_C = reg.fluid.U.phi[reg.mesh.owners[istartface:iendface], :]

    u_b = U_b[:, 0]
    v_b = U_b[:, 1]
    w_b = U_b[:, 2]

    u_C = U_C[:, 0]
    v_C = U_C[:, 1]
    w_C = U_C[:, 2]

    if iComponent == 0:
        reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * magSf_b[:, np.newaxis] / wallDist_b * (1 - nx2[:, np.newaxis])
        reg.fluxes.FluxFf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
        reg.fluxes.FluxVf[istartface:iendface, :] = -1 * mu_b * magSf_b[:, np.newaxis] / wallDist_b * (u_b * (1 - nx2) + (v_C - v_b) * ny * nx + (w_C - w_b) * nz * nx)[:, np.newaxis]
    elif iComponent == 1:
        reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * magSf_b[:, np.newaxis] / wallDist_b * (1 - ny2[:, np.newaxis])
        reg.fluxes.FluxFf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
        reg.fluxes.FluxVf[istartface:iendface, :] = -1 * mu_b * magSf_b[:, np.newaxis] / wallDist_b * ((u_C - u_b) * nx * ny + v_b * (1 - ny2) + (w_C - w_b) * nz * ny)[:, np.newaxis]
    else:
        reg.fluxes.FluxCf[istartface:iendface, :] = mu_b * magSf_b[:, np.newaxis] / wallDist_b * (1 - nz2[:, np.newaxis])
        reg.fluxes.FluxFf[istartface:iendface, :] = np.zeros((iendface - istartface, 1))
        reg.fluxes.FluxVf[istartface:iendface, :] = -1 * mu_b * magSf_b[:, np.newaxis] / wallDist_b * ((u_C - u_b) * nx * nz + (v_C - v_b) * ny * nz + w_b * (1 - nz2))[:, np.newaxis]

    reg.fluxes.FluxTf[istartface:iendface, :] = reg.fluxes.FluxCf[istartface:iendface, :] * U_C[:, [iComponent]] + reg.fluxes.FluxFf[istartface:iendface,:] * U_b[:,[iComponent]] + reg.fluxes.FluxVf[istartface:iendface,:]