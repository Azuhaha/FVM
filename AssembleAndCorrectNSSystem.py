from SetupCoefficients import *
from PreAssembleMomentumEquation import *
from AssembleMomentumEquationTerms import  *
from PostAssembleMomentumEquation import *
from SolveEquation import *
from CorrectVelocityEquation import *
from UpdateGradients import updateGradient
from UpdateScales import *
from SetupCoefficients import *
from PreAssembleContinuityEquation import *
from AssembleContinuityEquationTerms import *
from PostAssembleContinuityEquation import *
from CorrectNSSystemField import *

def assembleAndCorrectNSSystem(reg):
    #Momentum
    for iComponent in range(3):
        setupCoefficients(reg,0) #系数初始化

        #print('\tAssembling Momentum Equation...')
        assembleMomentumEquation(iComponent,reg)

        #print('\tSolving Momentum Equation...')
        solveMomentumEquation(iComponent, reg)

        setupCoefficients(reg, 0)

    #Continuity
    #print('\tAssembling Continuity Equation...')
    assembleContinuityEquation(reg)

    #print('\tASolving Continuity Equation...')
    solveContinuityEquation(reg)

    setupCoefficients(reg, 0)






def assembleMomentumEquation(iComponent,reg):
    preAssembleMomentumEquation(iComponent,reg)

    assembleMomentumEquationTerms(iComponent, reg)
    #
    postAssembleMomentumEquation(iComponent, reg)

def solveMomentumEquation(iComponent, reg):
    #solve equation
    #print('\t\tSolving...')
    solveEquation('U',iComponent,reg)

    #correct equation 根据dphi更新phi
    #print('\t\tCorrecting Velocity Equation...')
    correctVelocityEquation(iComponent,reg)

    #updates
    #print('\t\tUpdating Gradient and Scale...')
    #根据最新的phi更新梯度值
    updateGradient(reg,'U')

    updatescale(getattr(reg.fluid,'U'),reg)


def assembleContinuityEquation(reg):
    preAssembleContinuityEquation(reg)

    assembleContinuityEquationTerms(reg)
    #
    postAssembleContinuityEquation(reg)

def solveContinuityEquation(reg):
    #print('\t\tSolving...')
    solveEquation('p',0,reg)

    #print('\t\tCorrecting NS System Field...')
    correctNSSystemField(reg)


