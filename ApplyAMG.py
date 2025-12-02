from Class.CoefficientClass import *
import numpy as np
from SetupCoefficients import *
from SolveAlgebraicSystem import *
from ComputeResidualsArray import *
def applyAMG(preconditioner, maxIter, tolerance, relTol, nPreSweeps, nPostSweeps, nFinestSweeps,reg):
    #默认为V循环
    cycleType='V-Cycle'

    maxCoarseLevels=10 #最大的粗化等级

    minNumOfParents=5 #最小的聚合后单元的数量

    #单元聚合并更新各级别的系数
    iLevel=0
    while (iLevel+1)<=maxCoarseLevels:
        iLevel+=1
        numOfParents=agglomerateLevel(iLevel,reg) #创建parent的coefficient并更新其中的numofelement、csize和cconn
        assembleAgglomeratedLHS(iLevel,reg) #更新parent的coefficient中ac,anb
        if numOfParents<=minNumOfParents:
            break
    maxLevels=iLevel+1

    #calculate initial residual
    theCoefficient=reg.coefficient[0]
    residualsArray=computeResidualsArray(theCoefficient)
    initialResidual = np.sum(np.abs(residualsArray))
    finalResidual = initialResidual

    #Cycle
    if maxLevels<=3 or cycleType=='V_Cycle':
        for iter in range(maxIter):
            finalResidual=applyVCycle(0,preconditioner,maxLevels, nPreSweeps, nPostSweeps,relTol, nFinestSweeps,reg)
            if finalResidual<relTol*initialResidual and finalResidual<tolerance:
                break

    return [initialResidual,finalResidual]

# 创建parent的coefficient并更新其中的numofelement、csize和cconn
def agglomerateLevel(iLevel,reg):
    theCoefficients=reg.coefficient[iLevel-1]
    theNumberOfFineElements = theCoefficients.ac.shape[0]


    anb = theCoefficients.anb
    cconn = theCoefficients.cconn
    csize = theCoefficients.csize

    parents = -1.*np.ones((theNumberOfFineElements, 1))
    maxAnb = np.zeros((theNumberOfFineElements, 1))

    for iElement in range(theNumberOfFineElements):

        maxAnb[iElement,0]=np.max(-1*anb[iElement])

    iParent=0

    #Step1:Agglomeration
    for iSeed in range(theNumberOfFineElements):
        if parents[iSeed,0]==-1:
            parents[iSeed,0]=iParent
            children=[iSeed]
            for iNB_local in range(int(csize[iSeed,0])):
                iNB=cconn[iSeed][iNB_local]
                #print (cconn[iSeed][iNB_local])
                if parents[iNB,0]==0:
                    if -1*anb[iSeed][iNB_local,0]/maxAnb[iSeed,0]>0.5:
                        parents[iNB,0]=iParent
                        children.append(iNB)

            theNumberOfChildren=len(children)
            children2=[]
            for iChild_local in range(1,theNumberOfChildren):
                iChild=children[iChild_local]
                for iChildNB_local in range(csize[iChild,0]):
                    iChildNB=cconn[iChild][iChildNB_local]
                    if parents[iChildNB]==0:
                        if -1*anb[iChild][iChildNB_local,0]/maxAnb[iChild,0]>0.5:
                            parents[iChildNB,0]=iParent
                            children2.append(iChildNB)

            theNumberOfChildren=len(children)+len(children2)

            if theNumberOfChildren==1:
                parents[iSeed,0]=-1
            else:
                iParent+=1

    # Step2:Last step agglomeration 找到孤立单元并归入相关性最大的parents
    for iOrphan in range(theNumberOfFineElements):
        if parents[iOrphan,0]==-1:
            strength=0
            for iNB_local in range(int(csize[iOrphan,0])):
                iNB=cconn[iOrphan][iNB_local]
                if parents[iNB,0]!=-1:
                    if strength<-1*anb[iOrphan][iNB_local,0]/maxAnb[iNB,0]:
                        strength=-1*anb[iOrphan][iNB_local,0]/maxAnb[iNB,0]
                        parents[iOrphan,0]=parents[iNB,0]
        #仍未加入parents则单独定为parent
        if parents[iOrphan,0]==-1:
            parents[iOrphan,0]=iParent
            iParent+=1

    theNumberOfParents=iParent
    theCoefficients.parents=parents
    reg.coefficient[iLevel - 1]=theCoefficients


    #setup connectivity and csize
    #theParentCConn和theParentCSize赋给下一级的coefficient
    theParentCConn=[[] for _ in range(theNumberOfParents)]
    for iElement in range(theNumberOfFineElements):
        for iNB_local in range(int(csize[iElement,0])):
            iNB=cconn[iElement][iNB_local]
            if parents[iElement,0]!=parents[iNB,0]:
                if not parents[iNB,0] in theParentCConn[int(parents[iElement,0])]:
                    theParentCConn[int(parents[iElement, 0])].append(int(parents[iNB,0]))

    theParentCSize=np.zeros((theNumberOfParents,1))
    for iCoarseElement in range(theNumberOfParents):
        theParentCSize[iCoarseElement,0]=len(theParentCConn[iCoarseElement])

    #setup coefficients for coarser level
    reg.coefficient[iLevel]=Coefficient()
    reg.coefficient[iLevel].numOfElements=theNumberOfParents
    reg.coefficient[iLevel].ac=np.zeros((theNumberOfParents,1))
    reg.coefficient[iLevel].ac_old=np.zeros((theNumberOfParents,1))
    reg.coefficient[iLevel].bc=np.zeros((theNumberOfParents,1))
    reg.coefficient[iLevel].dc=np.zeros((theNumberOfParents,1))
    reg.coefficient[iLevel].rc=np.zeros((theNumberOfParents,1))
    reg.coefficient[iLevel].dphi=np.zeros((theNumberOfParents,1))
    reg.coefficient[iLevel].csize=theParentCSize
    reg.coefficient[iLevel].cconn=theParentCConn
    reg.coefficient[iLevel].anb = []
    for ielement in range(theNumberOfParents):
        reg.coefficient[iLevel].anb.append(np.zeros((int(theParentCSize[ielement,0]), 1)))

    return theNumberOfParents

#更新parent的coefficient中ac,anb
def assembleAgglomeratedLHS(iLevel,reg):
    theCoefficients = reg.coefficient[iLevel - 1]
    parents = theCoefficients.parents
    ac = theCoefficients.ac
    anb = theCoefficients.anb
    cconn = theCoefficients.cconn
    csize = theCoefficients.csize

    theParentCoefficients = reg.coefficient[iLevel]
    AC = theParentCoefficients.ac
    ANB = theParentCoefficients.anb
    CCONN = theParentCoefficients.cconn

    for iElement in range(theCoefficients.numOfElements):
        iParent=int(parents[iElement,0])
        AC[iParent,0]+=ac[iElement,0] #parent的AC中包含所有的元素的ac和该元素的相邻的属于parent的元素anb相加
        for iNB_local in range(int(csize[iElement,0])):
            iNB=cconn[iElement][iNB_local]
            iNBParent=int(parents[iNB,0])
            if iNBParent==iParent:
                AC[iParent,0]+=anb[iElement][iNB_local,0]
            else:
                iNBParent_local=CCONN[iParent].index(iNBParent)
                ANB[iParent][iNBParent_local,0]+=anb[iElement][iNB_local,0] #parent的ANB中包含两个parent中各自的且相邻的元素之间的anb之和

    reg.coefficient[iLevel].ac =AC
    reg.coefficient[iLevel].anb = ANB
    reg.coefficient[iLevel].cconn =CCONN


#
def applyVCycle(gridlevel,preconditioner,maxLevels, nPreSweeps, nPostSweeps,relTol, nFinestSweeps,reg):
    #restriction phase
    while gridlevel<maxLevels-1:
        #pre-sweep
        solveAlgebraicSystem(gridlevel,reg,preconditioner,nPreSweeps)
        #restrict residuals 限制算子，低频误差传递到粗网格
        restrict(gridlevel,reg)
        #update level
        gridlevel+=1

    #smoothening the coarsest level
    solveAlgebraicSystem(gridlevel,reg,preconditioner,nPreSweeps)

    #prolongatoin phase
    while gridlevel>0:
        if gridlevel==1:
            #prolongate correction to finer solution
            prolongate(gridlevel,reg)#在粗网格上进一步处理低频误差，通过延拓算子（prolongation）传回细网格
            #post-sweep finest level
            solveAlgebraicSystem(gridlevel-1,reg,preconditioner,nFinestSweeps)
        else:
            #prolongate correction to finer solution
            prolongate(gridlevel,reg)
            #post-sweep finest level
            solveAlgebraicSystem(gridlevel-1,reg,preconditioner,nPostSweeps)
        gridlevel-=1

    theCoefficient=reg.coefficient[0]
    residual=computeResidualsArray(theCoefficient)
    finalResidual=np.sum(np.abs(residual))

    return finalResidual

def restrict(gridlevel,reg):
    thecoefficient=reg.coefficient[gridlevel]
    residual=computeResidualsArray(thecoefficient)

    theParents=thecoefficient.parents
    theNumOfElements=thecoefficient.numOfElements

    theNumOfCoarseElements=reg.coefficient[gridlevel+1].numOfElements

    #initial coarse level RHS coefficient bc
    BC=np.zeros((theNumOfCoarseElements,1))
    for iFineElement in range(theNumOfElements):
        iParent=theParents[iFineElement,0]
        BC[iParent]+=residual[iFineElement]

    reg.coefficient[gridlevel + 1].bc=BC

def prolongate(gridlevel,reg):
    thecoefficient=reg.coefficient[gridlevel]
    theFinercoefficient = reg.coefficient[gridlevel-1]

    iParent=theFinercoefficient.parents
    theFinercoefficient.dphi+=thecoefficient.dphi[iParent]

    reg.coefficient[gridlevel - 1]=theFinercoefficient


































































