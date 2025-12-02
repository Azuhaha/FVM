import numpy as np
#根据coefficient中的值计算residual
def computeResidualsArray(theCoefficient):
    ac=theCoefficient.ac
    anb=theCoefficient.anb
    bc = theCoefficient.bc
    cconn = theCoefficient.cconn
    phi = theCoefficient.dphi
    residual=np.zeros(ac.shape)
    # bc-ac*dphi-Σ(anb[nei]*phi[nei])
    for iElement in range(ac.shape[0]):
        residual[iElement,0]=bc[iElement,0]-ac[iElement,0]*phi[iElement,0]
        for nNeighbour in range(len(cconn[iElement])):
            iNeighbour=cconn[iElement][nNeighbour]
            residual[iElement,0]-=anb[iElement][nNeighbour,0]*phi[iNeighbour,0]
    return residual