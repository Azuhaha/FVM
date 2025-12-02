from ComputeScaledRMSResiduals import *


def postAssembleContinuityEquation(reg):
    #enforce diagonal dominance as this may not be ensured
    assembleDiagDominance(reg)

    #compute scaled residuals
    computeScaledRMSResiduals('p',0,reg)

def assembleDiagDominance(reg):
    csize=reg.coefficient[0].csize
    anb=reg.coefficient[0].anb
    ac = reg.coefficient[0].ac

    for iElement in range(reg.mesh.numberOfElements):
        SumAik=0
        for k  in range(csize[iElement,0]):
            if anb[iElement][k,0]>0:
                anb[iElement][k, 0]=0
            SumAik-=anb[iElement][k, 0]
        ac[iElement,0]=max(ac[iElement,0],SumAik)

    reg.coefficient[0].anb=anb
    reg.coefficient[0].ac=ac
