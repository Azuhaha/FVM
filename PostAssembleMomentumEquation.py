import numpy as np
import math
from UpdateScalarFieldForAllBoundaryPatches import *
from ComputeScaledRMSResiduals import *
def postAssembleMomentumEquation(iComponent, reg):
    #Apply under-relaxation
    assembleImplicitRelaxation('U',reg)
    #Store DU and DUT
    assembleDCoefficients(iComponent,reg)
    #Compute RMS and MAX Residuals
    computeScaledRMSResiduals('U', iComponent,reg)

def assembleImplicitRelaxation(equationname,reg):
    urf=getattr(reg.foamDictionary.fvSolution.relaxationFactors.equations,equationname)
    #ac/urf
    reg.coefficient[0].ac /= urf #式14.26


def assembleDCoefficients(iComponent,reg):
    volume=reg.mesh.elementVolumes

    DUfield=getattr(reg.fluid,'DU'+str(iComponent+1))
    DUTfield = getattr(reg.fluid, 'DUT' + str(iComponent + 1))

    if 'SIMPLE' in dir(reg.foamDictionary.fvSolution):
        #V/ac
        DUfield.phi[:reg.mesh.numberOfElements,:]=volume/reg.coefficient[0].ac
        #ac_old/ac
        DUTfield.phi[:reg.mesh.numberOfElements, :] = reg.coefficient[0].ac_old / reg.coefficient[0].ac
    elif 'PIMPLE' in dir(reg.foamDictionary.fvSolution):
        raise Exception('PIMPLE Algorithm is developed')

    updateScalarFieldForAllBoundaryPatches(reg,DUfield.name)
    updateScalarFieldForAllBoundaryPatches(reg, DUTfield.name)




















