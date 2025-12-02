import numpy as np
def interpolateFromElementsToFaces(interpolationscheme,phi,reg):

    if interpolationscheme=='linear':

        phi_f1=phi[reg.mesh.neighbours,:]*reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace,:]+phi[reg.mesh.owners[:reg.mesh.numberOfInteriorFace],:]*(1-reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace,:])
        phi_f2=phi[reg.mesh.numberOfElements:reg.mesh.numberOfElements+reg.mesh.numberOfBElements,:]

        phi_f=np.vstack((phi_f1, phi_f2))

        return phi_f


    else:
        raise Exception('other interpolationscheme developed')

def interpolateFromElementsToInteriorFaces(interpolationscheme,phi,reg):
    if interpolationscheme=='linear':
        phi_f=phi[reg.mesh.neighbours,:]*reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace,:]+phi[reg.mesh.owners[:reg.mesh.numberOfInteriorFace],:]*(1-reg.mesh.faceWeights[:reg.mesh.numberOfInteriorFace,:])
    elif interpolationscheme=='vanLeerV':
        #(Vol_nei+Vol_own)*phi_own*phi_nei/(Vol_nei*phi_own+Vol_own*phi_nei)
        phi_f = (reg.mesh.elementVolumes[reg.mesh.neighbours, :]+reg.mesh.elementVolumes[reg.mesh.owners[:reg.mesh.numberOfInteriorFace],:])*phi[reg.mesh.owners[:reg.mesh.numberOfInteriorFace],:]*phi[reg.mesh.neighbours, :]/(reg.mesh.elementVolumes[reg.mesh.neighbours, :]*phi[reg.mesh.owners[:reg.mesh.numberOfInteriorFace],:]+reg.mesh.elementVolumes[reg.mesh.owners[:reg.mesh.numberOfInteriorFace]]*phi[reg.mesh.neighbours, :])
    elif interpolationscheme=='linearUpwind':
        pos = np.where(reg.fluid.mdot_f.phi[0:reg.mesh.numberOfInteriorFace, :] > 0, 1, 0)
        phi_f=phi[reg.mesh.neighbours,:]*(1-pos)+phi[reg.mesh.owners[:reg.mesh.numberOfInteriorFace],:]*pos
    else:
        raise Exception('interpolation scheme incorrect')

    return phi_f


