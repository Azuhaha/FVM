import numpy as np
def interpolateGradientsFromElementsToInteriorFaces(interpolationScheme,gradPhi,phi,reg):
    numOfComponent=phi.shape[1]
    grad_f=np.zeros((reg.mesh.numberOfInteriorFace,3,numOfComponent))
    if interpolationScheme=='linear' or interpolationScheme=='Gauss linear':
        for iCom in range(numOfComponent):
            grad_f[:,0,iCom]=(1.-reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace,0])*gradPhi[reg.mesh.neighbours,0,iCom]+reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace,0]*gradPhi[reg.mesh.owners[0:reg.mesh.numberOfInteriorFace],0,iCom]
            grad_f[:, 1, iCom] = (1. - reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace, 0]) * gradPhi[reg.mesh.neighbours, 1, iCom] + reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace, 0] * gradPhi[reg.mesh.owners[0:reg.mesh.numberOfInteriorFace], 1, iCom]
            grad_f[:, 2, iCom] = (1. - reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace, 0]) * gradPhi[reg.mesh.neighbours, 2, iCom] + reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace, 0] * gradPhi[reg.mesh.owners[0:reg.mesh.numberOfInteriorFace], 2, iCom]
    elif interpolationScheme=='Gauss linear corrected':
        for iCom in range(numOfComponent):
            grad_f[:,0,iCom]=(1.-reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace,0])*gradPhi[reg.mesh.neighbours,0,iCom]+reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace,0]*gradPhi[reg.mesh.owners[0:reg.mesh.numberOfInteriorFace],0,iCom]
            grad_f[:, 1, iCom] = (1. - reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace, 0]) * gradPhi[reg.mesh.neighbours, 1, iCom] + reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace,0] * gradPhi[reg.mesh.owners[0:reg.mesh.numberOfInteriorFace], 1, iCom]
            grad_f[:, 2, iCom] = (1. - reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace, 0]) * gradPhi[reg.mesh.neighbours, 2, iCom] + reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace, 0] * gradPhi[reg.mesh.owners[0:reg.mesh.numberOfInteriorFace], 2, iCom]

            CF_mag=np.linalg.norm(reg.mesh.faceCF[:reg.mesh.numberOfInteriorFace,:],axis=1)
            CF_e=reg.mesh.faceCF[:reg.mesh.numberOfInteriorFace,:]/CF_mag[:,np.newaxis]

            local_grad_CF=(phi[reg.mesh.neighbours,iCom]-phi[reg.mesh.owners[:reg.mesh.numberOfInteriorFace],iCom])/CF_mag
            local_grad=local_grad_CF[:,np.newaxis]*CF_e

            local_avg_grad_f=np.sum(grad_f[:,:,iCom]*CF_e,axis=1)
            local_avg_grad=local_avg_grad_f[:,np.newaxis]*CF_e

            #强行使界面处沿CF方向的梯度等于由点C和F上值所定义的局部梯度,p213式9.33
            grad_f[:,:,iCom]=grad_f[:,:,iCom]-local_avg_grad+local_grad

    elif  interpolationScheme=='Gauss upwind':
        pos = np.where(reg.fluid.mdot_f.phi[0:reg.mesh.numberOfInteriorFace, :] > 0, 1, 0)
        #pos=1表示上游为Owner单元，pos=0表示上游为neighbour单元，这里相反？
        grad_f[:,:,:]=(pos[:,np.newaxis,np.newaxis]*gradPhi[reg.mesh.neighbours,:,:]
                       +(1-pos[:,np.newaxis,np.newaxis])*gradPhi[reg.mesh.owners[:reg.mesh.numberOfInteriorFace],:,:])

    else:
        raise Exception("interpolation Scheme incorrect")

    return grad_f
