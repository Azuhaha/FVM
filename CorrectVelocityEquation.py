import numpy as np

def correctVelocityEquation(iComponent,reg):
    #interior
    theNumOfElement=reg.mesh.numberOfElements
    #根据dphi更新phi
    reg.fluid.U.phi[:theNumOfElement,[iComponent]]+=reg.coefficient[0].dphi

    reg.coefficient[0].dphi=np.zeros((theNumOfElement,1))

    #boundary
    for iBPatch in range(reg.mesh.numberOfBoundaryPatches):
        thePhysicalType=reg.mesh.cfdBoundaryPatchesArray[iBPatch].type
        theBCType=reg.fluid.U.boundaryPatchRef[iBPatch].type

        #WALL
        if thePhysicalType=='wall':
            if theBCType=='noSlip':
                continue
            elif theBCType=='slip':
                correctVelocitySlipWall(iBPatch,iComponent,reg)
            elif theBCType=='fixedValue':
                continue
            else:
                raise Exception('Wall Condition not implemented')

        elif thePhysicalType=='inlet':
            if theBCType=='fixedValue':
                continue
            elif theBCType=='zeroGradient' or theBCType=='inlet':
                correctVelocityZeroGradient(iBPatch, iComponent,reg)
            else:
                raise Exception('Inlet Condition not implemented')

        elif thePhysicalType=='outlet':
            if theBCType=='fixedValue':
                continue
            elif theBCType=='zeroGradient' or theBCType=='outlet':
                correctVelocityZeroGradient(iBPatch, iComponent,reg)
            else:
                raise Exception('Outlet Condition not implemented')

        elif thePhysicalType in ['symmetry','symmetryPlane','empty']:
            correctVelocitySymmetry(iBPatch, iComponent,reg)

        else:
            raise Exception('thePhysicalType Condition not implemented')

#边界面的速度值等于owner单元的值
def correctVelocityZeroGradient(iBPatch, iComponent,reg):
    startBElement=reg.mesh.numberOfElements+reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex-reg.mesh.numberOfInteriorFace
    endBElement=startBElement+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    startBFace=reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    endBFace=startBFace+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    iOwners=reg.mesh.owners[startBFace:endBFace]

    reg.fluid.U.phi[startBElement:endBElement,iComponent]=reg.fluid.U.phi[iOwners,iComponent]

#边界速度值=owner单元的值-owner单元的值在Sf方向的分量*n在icom的分量
def correctVelocitySlipWall(iBPatch, iComponent, reg):
    startBElement=reg.mesh.numberOfElements+reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex-reg.mesh.numberOfInteriorFace
    endBElement=startBElement+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    startBFace=reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    endBFace=startBFace+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    iOwners=reg.mesh.owners[startBFace:endBFace]

    Sb=reg.mesh.faceSf[startBFace:endBFace,:]
    n=Sb/np.linalg.norm(Sb,axis=1)[:,np.newaxis]

    U_normal=np.sum(reg.fluid.U.phi[iOwners,:]*n,axis=1)
    reg.fluid.U.phi[startBElement:endBElement,iComponent]=reg.fluid.U.phi[iOwners,iComponent]-U_normal*n[:,iComponent]

def correctVelocitySymmetry(iBPatch, iComponent, reg):
    startBElement=reg.mesh.numberOfElements+reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex-reg.mesh.numberOfInteriorFace
    endBElement=startBElement+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    startBFace=reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    endBFace=startBFace+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    iOwners=reg.mesh.owners[startBFace:endBFace]

    Sb=reg.mesh.faceSf[startBFace:endBFace,:]
    n=Sb/np.linalg.norm(Sb,axis=1)[:,np.newaxis]

    U_normal=np.sum(reg.fluid.U.phi[iOwners,:]*n,axis=1)
    reg.fluid.U.phi[startBElement:endBElement,iComponent]=reg.fluid.U.phi[iOwners,iComponent]-U_normal*n[:,iComponent]