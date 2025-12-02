import os
from SetupCoefficients import *
from SetupFluxes import *
from UpdateFieldsForAllBoundaryPatches import *
from UpdateGradients import *
from UpdateScales import *
from UpdatePrevIter import *
from UpdateProperties import *
from AssembleAndCorrectNSSystem import *
from AssembleAndCorrectEnergyEquation import *
from PrintResiduals import *
from PlotResiduals import *
import matplotlib.pyplot as plt

def runFalseTransientCase(reg):
    #print('Initializing Coefficients[0],Fluxes...')
    setupCoefficients(reg,0) #内部单元
    setupFluxes(reg)

    #Initialize cfdRun time
    #print('Initializing StartTime...')
    if reg.foamDictionary.controlDict.startFrom=='startTime':
        reg.time.currentTime = reg.foamDictionary.controlDict.startTime
    else:
        raise Exception('developed')

    if not os.path.exists(os.path.join(reg.caseDirectoryPath,'convergence')):
        os.mkdir(os.path.join(reg.caseDirectoryPath,'convergence'))
    with open(os.path.join(reg.caseDirectoryPath,'convergence','convergenceUp.out'),'w') as f:
        f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'.format('noIter','Time[s]','UxResInit','UyResInit','UzResInit','pResInit','kResInit','epsilonResInit','omegaResInit','TResInit'))

    #Pre - updates(necessary on startup)
    #print('Updating Fields for All BoundaryPatches(U,p,T,rho,mu,Cp,k)...')
    updateFieldsForAllBoundaryPatches(reg)

    #print('Updating Gradients(U,p,T,rho)...')
    updateGradients(reg)

    #print('Updating Scales(U,p,T,rho)...')
    updatescales(reg)

    #开始迭代计算
    print('\nStarting Iteration...')
    numOfIterations=0

    while reg.time.currentTime<reg.foamDictionary.controlDict.endTime:
        numOfIterations+=1
        print(f'\nStep{numOfIterations} of Iteration...')

        reg.time.currentTime+=reg.foamDictionary.controlDict.deltaT


        #print

        #upadate
        #print(f'\tUpdating Data of Previous Iteration and Properties...')
        updatePrevIter(reg)
        updateProperties(reg)

        #important N-S Pressure correction
        assembleAndCorrectNSSystem(reg)

        #energy
        assembleAndCorrectEnergyEquation(reg)

        #printresidual
        printResiduals(numOfIterations,reg)
        plotResiduals(numOfIterations,reg)

        if numOfIterations==2:
            break


    plt.show()












