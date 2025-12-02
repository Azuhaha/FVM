import numpy as np
def updateVectorFieldForAllBoundaryPatches(reg,fluidfieldname):
    for iBPatch in range(reg.mesh.numberOfBoundaryPatches):
        thePhysicalPatchType=reg.mesh.cfdBoundaryPatchesArray[iBPatch].type
        fluidfield=getattr(reg.fluid, fluidfieldname)
        theBCType=fluidfield.boundaryPatchRef[iBPatch].type

        #wall
        if thePhysicalPatchType=='wall':
            if theBCType=='fixedValue':
                updateFixedValue(reg,iBPatch,fluidfield)
            elif theBCType=='zeroGradient' or theBCType=='noSlip' or theBCType=='slip':
                updateZeroGradient(reg,iBPatch,fluidfield)
            else:
                raise ValueError("wall bc not defined")
        #inlet
        elif thePhysicalPatchType=='inlet':
            if theBCType=='fixedValue':
                updateFixedValue(reg,iBPatch,fluidfield)
            elif theBCType=='zeroGradient':
                updateZeroGradient(reg,iBPatch,fluidfield)
            else:
                raise ValueError("inlet bc not defined")

        #outlet
        elif thePhysicalPatchType=='outlet':
            if theBCType=='fixedValue':
                updateFixedValue(reg,iBPatch,fluidfield)
            elif theBCType=='zeroGradient' or theBCType=='outlet':
                updateZeroGradient(reg,iBPatch,fluidfield)
            else:
                raise ValueError("outlet bc not defined")

        #sym
        elif thePhysicalPatchType=='symmetry':
                updateSymmetry(reg,iBPatch,fluidfield)

        #empty
        elif thePhysicalPatchType=='empty':
                updateSymmetry(reg,iBPatch,fluidfield)

        else:
            raise ValueError("Physical Condition bc not defined")

def updateFixedValue(reg,iBPatch,fluidfield):
    startBElement=reg.mesh.numberOfElements+reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex-reg.mesh.numberOfInteriorFace
    endBElement=reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces+startBElement

    fluidfield.phi[startBElement:endBElement,:]=np.tile(fluidfield.boundaryPatchRef[iBPatch].value,(endBElement-startBElement,1))


def updateZeroGradient(reg, iBPatch, fluidfield):
    startBElement=reg.mesh.numberOfElements+reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex-reg.mesh.numberOfInteriorFace
    endBElement=reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces+startBElement
    startBFace=reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    endBFace=startBFace+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    fluidfield.phi[startBElement:endBElement,:]=fluidfield.phi[reg.mesh.owners[startBFace:endBFace],:]


def updateSymmetry(reg, iBPatch, fluidfield):
    startBElement=reg.mesh.numberOfElements+reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex-reg.mesh.numberOfInteriorFace
    endBElement=reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces+startBElement
    startBFace=reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    endBFace=startBFace+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    Sf=reg.mesh.faceSf[startBFace:endBFace, :]
    norm_Sf=np.linalg.norm(Sf,axis=1)
    n=Sf/norm_Sf[:,np.newaxis]

    phi_normal_mag=np.sum(fluidfield.phi[reg.mesh.owners[startBFace:endBFace],:]*n,axis=1)
    phi_normal=n*phi_normal_mag[:,np.newaxis]

    fluidfield.phi[startBElement:endBElement,:]=fluidfield.phi[reg.mesh.owners[startBFace:endBFace],:]-phi_normal





