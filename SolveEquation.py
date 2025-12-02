from ApplyAMG import *

def solveEquation(equationname,iComponent,reg):
    solver=getattr(reg.foamDictionary.fvSolution.solvers,equationname).solver
    maxIter = getattr(reg.foamDictionary.fvSolution.solvers, equationname).maxIter
    tolerance = getattr(reg.foamDictionary.fvSolution.solvers, equationname).tolerance
    relTol = getattr(reg.foamDictionary.fvSolution.solvers, equationname).relTol

    if solver=='GAMG':
        preconditioner = getattr(reg.foamDictionary.fvSolution.solvers, equationname).preconditioner
        nPreSweeps = getattr(reg.foamDictionary.fvSolution.solvers, equationname).nPreSweeps
        nPostSweeps = getattr(reg.foamDictionary.fvSolution.solvers, equationname).nPostSweeps
        nFinestSweeps = getattr(reg.foamDictionary.fvSolution.solvers, equationname).nFinestSweeps
        #print('\t\t\tApplying AMG...')
        [initRes,finalRes]=applyAMG(preconditioner, maxIter, tolerance, relTol, nPreSweeps, nPostSweeps, nFinestSweeps,reg)

    elif solver=='smoothSolver':
        smoother=getattr(reg.foamDictionary.fvSolution.solvers, equationname).smoother
        [initRes,finalRes]=solveAlgebraicSystem(0,reg,smoother,maxIter,tolerance,relTol)

    else:
        raise Exception("solver is not defined")

    theequation=getattr(reg.model, equationname)

    theequation.residuals.initResidual[iComponent]=initRes
    theequation.residuals.finalResidual[iComponent]=finalRes





