import numpy as np
def correctMdot(reg):
    ###interior
    owner_f=reg.mesh.owners[:reg.mesh.numberOfInteriorFace]
    nei_f=reg.mesh.neighbours

    #correct
    #0.75系数由于未考虑非正交分量
    reg.fluid.mdot_f.phi[:reg.mesh.numberOfInteriorFace,:]+=0.75*(reg.fluxes.FluxCf[:reg.mesh.numberOfInteriorFace,:]*reg.fluid.pp.phi[owner_f,:]+reg.fluxes.FluxFf[:reg.mesh.numberOfInteriorFace,:]*reg.fluid.pp.phi[nei_f,:])

    ###boundary
    for iBPatch in range(reg.mesh.numberOfBoundaryPatches):
        thePhysicalType=reg.mesh.cfdBoundaryPatchesArray[iBPatch].type
        theBCType=reg.fluid.p.boundaryPatchRef[iBPatch].type

        if thePhysicalType=='wall':
            continue
        elif thePhysicalType=='inlet':
            if theBCType=='inlet' or theBCType=='zeroGradient':
                correctMdotZeroGradient(iBPatch,reg)
            elif theBCType=='fixedValue':
                correctMdotFixedValue(iBPatch,reg)
            else:
                raise Exception('theBCType not implemented')
        elif thePhysicalType=='outlet':
            if theBCType=='outlet' or theBCType=='zeroGradient':
                correctMdotZeroGradient(iBPatch,reg)
            elif theBCType=='fixedValue':
                correctMdotFixedValue(iBPatch,reg)
            else:
                raise Exception('theBCType not implemented')
        elif thePhysicalType=='symmetry' or thePhysicalType=='empty':
            continue
        else:
            raise Exception('thePhysicalType not implemented')


def correctMdotZeroGradient(iBPatch,reg):
    iBstartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iBendface = iBstartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    iBstartelement = reg.mesh.numberOfElements + iBstartface - reg.mesh.numberOfInteriorFace
    iBendelement = iBstartelement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    rho_b=reg.fluid.rho.phi[iBstartelement:iBendelement,:]
    Sb=reg.mesh.faceSf[iBstartface:iBendface,:]
    owner_b=reg.mesh.owners[iBstartface:iBendface]
    U_b=reg.fluid.U.phi[owner_b,:]

    #update and store mdot
    reg.fluid.mdot_f.phi[iBstartface:iBendface,:]=0.75*rho_b*np.sum(U_b*Sb,axis=1)[:,np.newaxis]+(1-0.75)*reg.fluid.mdot_f.phi[iBstartface:iBendface,:]

def correctMdotFixedValue(iBPatch, reg):
    iBstartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iBendface = iBstartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    owner_b = reg.mesh.owners[iBstartface:iBendface]
    # correct
    reg.fluid.mdot_f.phi[iBstartface:iBendface, :] += 0.75 * reg.fluxes.FluxCf[iBstartface:iBendface, :]*reg.fluid.pp.phi[owner_b,:]

