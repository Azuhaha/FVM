import numpy as np
def assembleMomentumConvectionTerm(iComponent,reg):
    ###内部面
    #||mdot_f,0||
    reg.fluxes.FluxCf[0:reg.mesh.numberOfInteriorFace,:]=np.maximum(reg.fluid.mdot_f.phi[:reg.mesh.numberOfInteriorFace,:],0.)
    # -||-mdot_f,0||
    reg.fluxes.FluxFf[0:reg.mesh.numberOfInteriorFace,:]=-1.*np.maximum(-1.*reg.fluid.mdot_f.phi[:reg.mesh.numberOfInteriorFace,:],0.)
    reg.fluxes.FluxVf[0:reg.mesh.numberOfInteriorFace,:]=np.zeros((reg.mesh.numberOfInteriorFace,1))
    #||mdot_f,0||*Uown-||-mdot_f,0||*Unei+0 每个内部面的通量
    U_phi_own=reg.fluid.U.phi[reg.mesh.owners[0:reg.mesh.numberOfInteriorFace], iComponent]
    U_phi_nei=reg.fluid.U.phi[reg.mesh.neighbours, iComponent]
    reg.fluxes.FluxTf[0:reg.mesh.numberOfInteriorFace,:]=reg.fluxes.FluxCf[0:reg.mesh.numberOfInteriorFace,:]*U_phi_own[:,np.newaxis]+reg.fluxes.FluxFf[0:reg.mesh.numberOfInteriorFace,:]*U_phi_nei[:,np.newaxis]+reg.fluxes.FluxVf[0:reg.mesh.numberOfInteriorFace,:]

    ###边界面
    for iBPatch in range(reg.mesh.numberOfBoundaryPatches):

        if reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='wall':
            continue
        elif reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='inlet':
            if reg.fluid.U.boundaryPatchRef[iBPatch].type=='fixedValue':
                assembleConvectionTermSpecifiedValue(iComponent,iBPatch,reg)
            elif reg.fluid.U.boundaryPatchRef[iBPatch].type=='zeroGradient':
                assembleConvectionTermZeroGradient(iComponent, iBPatch, reg)
            else:
                raise Exception("BCtype not implemented")

        elif reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='outlet':
            if reg.fluid.U.boundaryPatchRef[iBPatch].type=='fixedValue':
                assembleConvectionTermSpecifiedValue(iComponent,iBPatch,reg)
            elif reg.fluid.U.boundaryPatchRef[iBPatch].type=='zeroGradient':
                assembleConvectionTermZeroGradient(iComponent, iBPatch, reg)
            else:
                raise Exception("BCtype not implemented")
        elif reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='empty' or reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='symmetry':
            continue
        else:
            raise Exception("PhysicalType nor implemented")

def assembleConvectionTermSpecifiedValue(iComponent,iBPatch,reg):
    iBstartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iBendface = iBstartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    iBstartelement = reg.mesh.numberOfElements + iBstartface - reg.mesh.numberOfInteriorFace
    iBendelement = iBstartelement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    reg.fluxes.FluxCf[iBstartface:iBendface, :] = np.maximum(reg.fluid.mdot_f.phi[iBstartface:iBendface, :], 0.)
    reg.fluxes.FluxFf[iBstartface:iBendface, :] = -1. * np.maximum(-1. * reg.fluid.mdot_f.phi[iBstartface:iBendface, :], 0.)
    reg.fluxes.FluxVf[iBstartface:iBendface, :] = np.zeros((reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces, 1))
    U_phi_own=reg.fluid.U.phi[reg.mesh.owners[iBstartface:iBendface], iComponent]
    U_phi_B=reg.fluid.U.phi[iBstartelement:iBendelement, iComponent]
    reg.fluxes.FluxTf[iBstartface:iBendface, :] = reg.fluxes.FluxCf[iBstartface:iBendface, :] * U_phi_own[:,np.newaxis] + reg.fluxes.FluxFf[iBstartface:iBendface, :] *U_phi_B[:,np.newaxis]+reg.fluxes.FluxVf[iBstartface:iBendface, :]

def assembleConvectionTermZeroGradient(iComponent, iBPatch, reg)  :
    iBstartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iBendface = iBstartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    iBstartelement = reg.mesh.numberOfElements + iBstartface - reg.mesh.numberOfInteriorFace
    iBendelement = iBstartelement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    reg.fluxes.FluxCf[iBstartface:iBendface, :] = np.maximum(reg.fluid.mdot_f.phi[iBstartface:iBendface, :], 0.)
    reg.fluxes.FluxFf[iBstartface:iBendface, :] = -1. * np.maximum(-1. * reg.fluid.mdot_f.phi[iBstartface:iBendface, :], 0.)
    reg.fluxes.FluxVf[iBstartface:iBendface, :] = np.zeros((reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces, 1))
    U_phi_own = reg.fluid.U.phi[reg.mesh.owners[iBstartface:iBendface], iComponent]
    U_phi_B = reg.fluid.U.phi[iBstartelement:iBendelement, iComponent]
    reg.fluxes.FluxTf[iBstartface:iBendface, :] = reg.fluxes.FluxCf[iBstartface:iBendface, :] * U_phi_own[:,np.newaxis] + reg.fluxes.FluxFf[iBstartface:iBendface,:] * U_phi_B[:,np.newaxis] + reg.fluxes.FluxVf[iBstartface:iBendface,:]