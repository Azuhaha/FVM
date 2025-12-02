import numpy as np
from ComputeResidualsArray import *

def solveAlgebraicSystem(gridlevel, reg, *args):
    if len(args) == 0:
        smoother = 'DILU'
        maxIter = 20
        tolerance = 1e-6
        relTol = 0.1
    elif len(args) == 2:
        smoother = args[0]
        maxIter = 20
        tolerance = 1e-6
        relTol = 0.1
    elif len(args) == 3:
        smoother = args[0]
        maxIter = args[1]
        tolerance = 1e-6
        relTol = 0.1
    elif len(args) == 4:
        smoother = args[0]
        maxIter = args[1]
        tolerance = args[2]
        relTol = 0.1
    elif len(args) == 5:
        smoother = args[0]
        maxIter = args[1]
        tolerance = args[2]
        relTol = args[3]

    theCoefficients = reg.coefficient[gridlevel]
    ac = theCoefficients.ac
    anb = theCoefficients.anb
    bc = theCoefficients.bc
    cconn = theCoefficients.cconn
    dphi = theCoefficients.dphi
    theNumberOfElements = theCoefficients.numOfElements

    residualsArray = computeResidualsArray(theCoefficients)
    initialResidual = np.sum(np.abs(residualsArray)) / theNumberOfElements
    finalResidual = initialResidual

    if maxIter == 0:
        return
    if smoother == 'DILU':
        # Factorize Ax=b (Apply incomplete upper lower decomposition)
        [dc, rc] = factorizeILU(ac, anb, bc, cconn)

        # solve system
        for iter in range(maxIter):
            dphi = solveILU(ac, anb, bc, dc, rc, cconn, dphi)

            theCoefficients.dphi=dphi

            #check if termination criterion satisfied
            residualsArray=computeResidualsArray(theCoefficients)

            finalResidual=np.sum(np.abs(residualsArray)) / theNumberOfElements

            if finalResidual<relTol*initialResidual and finalResidual<tolerance:
                break

    elif smoother == 'SOR' or smoother=='GaussSeidel':
        # solve system
        for iter in range(maxIter):
            dphi = solveSOR(ac, anb, bc, cconn, dphi)
            theCoefficients.dphi = dphi

            # check if termination criterion satisfied
            #res=b-A*phi*
            residualsArray = computeResidualsArray(theCoefficients)
            finalResidual = np.sum(np.abs(residualsArray)) / theNumberOfElements
            if finalResidual < relTol * initialResidual and finalResidual < tolerance:
                break

    theCoefficients.dphi=dphi
    reg.coefficient[gridlevel]=theCoefficients

    return [initialResidual,finalResidual]


def factorizeILU(ac, anb, bc, cconn):
    numOfElements = ac.shape[0]

    dc=np.zeros(ac.shape)
    rc = np.zeros(ac.shape)

    for i1 in range(numOfElements):
        dc[i1]=ac[i1]
    for i1 in range(numOfElements):
        dc[i1, 0] = 1.0 / dc[i1, 0]
        rc[i1]=bc[i1]
        i1NbList = cconn[i1]
        i1NNb = len(i1NbList)
        if i1 != numOfElements - 2:
            # loop over neighbours of iElement
            for j1_ in range(i1NNb):
                jj1 = i1NbList[j1_]  # i1的相邻元素jj1,编号为j1_
                # for all neighbour j>i do
                if jj1 > i1 and jj1 < numOfElements:
                    j1NbList = cconn[jj1]
                    j1NNb = len(j1NbList)
                    for i1_ in range(j1NNb):
                        if j1NbList[i1_] == i1:
                            # dc=dc-A[j][i]*dc[i]*A[i][j]
                            dc[jj1, 0] = dc[jj1, 0] - anb[jj1][i1_, 0] * dc[i1, 0] * anb[i1][j1_, 0]
                            break
    return [dc, rc]


def solveILU(ac, anb, bc, dc, rc, cconn, dphi):
    numOfElements = ac.shape[0]

    # update residuals array
    # -ac*dphi+bc-Σanb[j]*dphi[j]
    for iElement in range(numOfElements):
        conn = cconn[iElement]
        res = -1 * ac[iElement] * dphi[iElement] + bc[iElement]
        for iLocalNeighbour in range(len(conn)):
            j = conn[iLocalNeighbour]
            res -= anb[iElement][iLocalNeighbour] * dphi[j]
        rc[iElement] = res

    # forward substitution
    # rc[j]=rc-Σi(j>i)anb[j][i]*dc[i]*rc[i]
    for i1 in range(numOfElements):
        mat1 = dc[i1, 0] * rc[i1, 0]

        i1NbList = cconn[i1]
        i1NNb = len(i1NbList)

        # loop over neighbours of iElement
        for j1_ in range(i1NNb):
            j1 = i1NbList[j1_]  # i1的相邻元素jj1,编号为j1_
            # for all neighbour j>i do
            if j1 > i1 and j1 < numOfElements:
                j1NbList = cconn[j1]
                j1NNb = len(j1NbList)
                for i1_ in range(j1NNb):
                    if j1NbList[i1_] == i1:
                        mat2 = anb[j1][i1_, 0] * mat1
                        rc[j1, 0] -= mat2
                        break

    # backward substitution
    for i1 in range(numOfElements - 1, -1, -1):
        if i1 < numOfElements - 1:
            i1NbList = cconn[i1]
            i1NNb = len(i1NbList)
            # loop over neighbours of iElement
            for j1_ in range(i1NNb):
                j = i1NbList[j1_]
                if j > i1:
                    rc[i1] -= anb[i1][j1_] * rc[j]

        # compute product D*R
        mat1 = dc[i1] * rc[i1]
        rc[i1] = mat1

        # update dphi
        dphi[i1] += mat1

    return dphi

def solveSOR(ac, anb, bc, cconn, dphi):
    raise Exception("SOR is developed")