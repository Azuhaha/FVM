import numpy as np
def updateGradients(reg):
    for fluidfieldname in dir(reg.fluid):
        if 'U' == fluidfieldname or 'p' == fluidfieldname or 'T' == fluidfieldname or 'rho' == fluidfieldname:
            updateGradient(reg,fluidfieldname)


def updateGradient(reg,fluidfieldname):
    if 'Gauss linear' == reg.foamDictionary.fvSchemes.gradSchemes.default:
        fluidfield = getattr(reg.fluid, fluidfieldname)
        computeGradientGaussLinear0(reg,fluidfield) #高斯线性计算梯度

    else:
        raise Exception(f"{reg.foamDictionary.fvSchemes.gradSchemes.default} is develpoed\n")

    #边界处梯度更新
    for iBPatch in range(reg.mesh.numberOfBoundaryPatches):
        if reg.mesh.cfdBoundaryPatchesArray[iBPatch].type=='wall':
            updateBoundaryGradients(iBPatch,fluidfield,reg)
        elif reg.mesh.cfdBoundaryPatchesArray[iBPatch].type == 'inlet':
            updateBoundaryGradients(iBPatch, fluidfield,reg)
        elif reg.mesh.cfdBoundaryPatchesArray[iBPatch].type == 'outlet':
            updateBoundaryGradients(iBPatch, fluidfield,reg)
        elif reg.mesh.cfdBoundaryPatchesArray[iBPatch].type == 'symmetry' or reg.mesh.cfdBoundaryPatchesArray[iBPatch].type == 'empty':
            pass
        else:
            raise Exception('Boundary condition not recognized')


def computeGradientGaussLinear0(reg,fluidfield):
    numOfComponent = fluidfield.phi.shape[1]
    fluidfield.phiGradient = np.zeros((reg.mesh.numberOfElements + reg.mesh.numberOfBElements, 3, numOfComponent))
    for iCom in range(numOfComponent):
        phi_f = reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace,0] * fluidfield.phi[reg.mesh.neighbours, iCom] + (1 - reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace,0]) *fluidfield.phi[reg.mesh.owners[:reg.mesh.numberOfInteriorFace], iCom]
        for iface in range(reg.mesh.numberOfInteriorFace):

            fluidfield.phiGradient[reg.mesh.owners[iface], :, iCom] += phi_f[iface] * reg.mesh.faceSf[iface, :]
            fluidfield.phiGradient[reg.mesh.neighbours[iface], :, iCom] -= phi_f[iface] * reg.mesh.faceSf[iface, :]
        for iface in range(reg.mesh.numberOfInteriorFace, reg.mesh.numberOfFaces):
            fluidfield.phiGradient[reg.mesh.owners[iface], :, iCom] += fluidfield.phi[iface-reg.mesh.numberOfInteriorFace+reg.mesh.numberOfElements] * reg.mesh.faceSf[iface, :]
        fluidfield.phiGradient[:reg.mesh.numberOfElements,:, iCom] /= reg.mesh.elementVolumes

    # 边界单元的梯度等于边界面的owner单元的梯度
    fluidfield.phiGradient[reg.mesh.numberOfElements:, :, :] = fluidfield.phiGradient[reg.mesh.owners[reg.mesh.numberOfInteriorFace:], :, :]


#边界单元的梯度=owner单元的梯度-owner单元的梯度在Cf的分量+通过C和f处的phi值计算的Cf方向的梯度
def updateBoundaryGradients(iBPatch,fluidfield,reg):
    if 'volVectorField' == fluidfield.type:
        numOfComponent=3
    else:
        numOfComponent=1

    startbface=reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    endbface=startbface+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    startbelement=startbface-reg.mesh.numberOfInteriorFace+reg.mesh.numberOfElements
    endbelement=startbelement+reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    for iCom in range(numOfComponent):
        Cf=reg.mesh.faceCentroid[startbface:endbface,:]
        C=reg.mesh.elementCentroids[reg.mesh.owners[startbface:endbface],:]
        dCf=Cf-C
        norm_dCf=np.linalg.norm(dCf,axis=1)
        e=dCf/norm_dCf[:,np.newaxis]
        grad_owner=np.sum(fluidfield.phiGradient[reg.mesh.owners[startbface:endbface], :, iCom] * e, axis=1)
        grad_Cf=(fluidfield.phi[startbelement:endbelement, iCom] - fluidfield.phi[reg.mesh.owners[startbface:endbface], iCom]) / norm_dCf
        fluidfield.phiGradient[startbelement:endbelement,:,iCom]=fluidfield.phiGradient[reg.mesh.owners[startbface:endbface],:,iCom]-grad_owner[:,np.newaxis]*e+ grad_Cf[:,np.newaxis]*e





