import numpy as np
from FVM.Class.CoefficientClass import *
def setupCoefficients(reg,iLevel):
    reg.coefficient[iLevel]=Coefficient()
    if iLevel==0:
        theNumOfElements=reg.mesh.numberOfElements
        reg.coefficient[iLevel].csize = np.array([[len(x)] for x in reg.mesh.elementNeighbours])  # 相邻单元的数量
        reg.coefficient[iLevel].cconn = reg.mesh.elementNeighbours
    else:
        raise Exception('developed')
    reg.coefficient[iLevel].numOfElements=theNumOfElements

    reg.coefficient[iLevel].ac=np.zeros((theNumOfElements,1))
    reg.coefficient[iLevel].ac_old=np.zeros((theNumOfElements,1))
    reg.coefficient[iLevel].bc=np.zeros((theNumOfElements,1))
    reg.coefficient[iLevel].dc=np.zeros((theNumOfElements,1))
    reg.coefficient[iLevel].rc=np.zeros((theNumOfElements,1))
    reg.coefficient[iLevel].dphi=np.zeros((theNumOfElements,1))

    reg.coefficient[iLevel].anb=[]
    for ielement in range(theNumOfElements):
        reg.coefficient[iLevel].anb.append(np.zeros((reg.coefficient[iLevel].csize[ielement,0],1)))