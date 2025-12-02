from CorrectVelocityEquation import correctVelocitySlipWall
from CorrectVelocityEquation import correctVelocitySymmetry
from CorrectVelocityEquation import correctVelocityZeroGradient
from UpdateScalarFieldForAllBoundaryPatches import *
from UpdateGradients import updateGradient
from UpdateScales import updatescale
from UpdateGradients import updateGradient
from CorrectMdot import *
from UpdateProperty import *


def correctNSSystemField(reg):
    #correct pressure correction field
    setupPressureCorrection(reg)

    correctPressureCorrectionForInterior(reg)

    updateScalarFieldForAllBoundaryPatches(reg,'pp')
    updateGradient(reg,'pp')

    #correct mdot_f
    correctMdot(reg)

    #correct pressure
    correctPressureForInterior(reg)
    correctPressureForBoundaryPatches(reg)

    #correct U
    correctVelocityForInterior(reg)
    correctVelocityForBoundaryPatches(0,reg)
    correctVelocityForBoundaryPatches(1, reg)
    correctVelocityForBoundaryPatches(2, reg)

    updateProperty('rho',reg)

    #update
    updateScalarFieldForAllBoundaryPatches(reg,'p')
    updateScalarFieldForAllBoundaryPatches(reg, 'rho')
    updateScalarFieldForAllBoundaryPatches(reg, 'U')

    #update scales
    updatescale(reg.fluid.p,reg)
    updatescale(reg.fluid.U, reg)
    updatescale(reg.fluid.rho, reg)

    #update gradients
    updateGradient(reg,'rho')
    updateGradient(reg, 'p')
    updateGradient(reg, 'U')


def setupPressureCorrection(reg):

    pRefValue=reg.foamDictionary.fvSolution.SIMPLE.pRefValue

    #check if need reference pressure
    if reg.mesh.closed:
        pRefCell=reg.foamDictionary.fvSolution.SIMPLE.pRefCell
        pRefValue+=reg.coefficient[0].dphi[pRefCell,0]

    reg.coefficient[0].dphi-=pRefValue

def correctPressureCorrectionForInterior(reg):
    #update pressure correction field
    reg.fluid.pp.phi[:reg.mesh.numberOfElements,:]=reg.coefficient[0].dphi[:reg.mesh.numberOfElements,:]
    #reset correction
    reg.coefficient[0].dphi=np.zeros((reg.mesh.numberOfElements,1))

def correctVelocityForInterior(reg):
    ppGrad=reg.fluid.pp.phiGradient[:reg.mesh.numberOfElements,:]
    DU1=reg.fluid.DU1.phi[:reg.mesh.numberOfElements,:]
    DU2 = reg.fluid.DU2.phi[:reg.mesh.numberOfElements, :]
    DU3 = reg.fluid.DU3.phi[:reg.mesh.numberOfElements, :]

    #Dc*gradP 式15.101
    DUPPGrad=np.hstack((DU1,DU2,DU3))*ppGrad[:,:,0]

    #correct velocity
    reg.fluid.U.phi[:reg.mesh.numberOfElements,:]-=DUPPGrad


def correctPressureForInterior(reg):
    #get limits
    if 'pmax' in dir(reg.foamDictionary.fvSolution.SIMPLE):
        theMaximunAcceptedValue=reg.foamDictionary.fvSolution.SIMPLE.pmax
    else:
        theMaximunAcceptedValue=1e7
    if 'pmin' in dir(reg.foamDictionary.fvSolution.SIMPLE):
        theMinimunAcceptedValue=reg.foamDictionary.fvSolution.SIMPLE.pmax
    else:
        theMinimunAcceptedValue=-1e7

    urf_p=reg.foamDictionary.fvSolution.relaxationFactors.fields.p

    if reg.compressible:
        for iElement in range(reg.mesh.numberOfElements):
            pcorr=urf_p*reg.fluid.pp.phi[iElement,0]
            p_new=reg.fluid.p.phi[iElement,0]*pcorr

            if p_new>theMaximunAcceptedValue:
                p_new=reg.fluid.p.phi[iElement,0]+0.75*(theMaximunAcceptedValue-reg.fluid.p.phi[iElement,0])
            elif p_new<theMinimunAcceptedValue:
                p_new = reg.fluid.p.phi[iElement, 0] + 0.75 * (theMinimunAcceptedValue - reg.fluid.p.phi[iElement, 0])

            reg.fluid.p.phi[iElement,0]=p_new
    else:
        #式15.101
        reg.fluid.p.phi[:reg.mesh.numberOfElements, :]+=urf_p*reg.fluid.pp.phi[:reg.mesh.numberOfElements,:]

def correctPressureForBoundaryPatches(reg):
    for iBPatch in range(reg.mesh.numberOfBoundaryPatches):
        thePhysicalType = reg.mesh.cfdBoundaryPatchesArray[iBPatch].type
        theBCType = reg.fluid.p.boundaryPatchRef[iBPatch].type

        if thePhysicalType == 'wall':
            correctPressureZeroGradient(iBPatch,reg)
        elif thePhysicalType == 'inlet':
            if theBCType == 'inlet' or theBCType == 'zeroGradient':
                correctPressureZeroGradient(iBPatch, reg)
            else:
                raise Exception('inlet not defined')
        elif thePhysicalType == 'outlet':
            if theBCType == 'outlet' or theBCType == 'zeroGradient':
                correctPressureZeroGradient(iBPatch, reg)
            elif theBCType == 'fixedValue':
                correctPressureSpecifiedValue(iBPatch, reg)
            else:
                raise Exception('oulet not defined')
        elif thePhysicalType == 'symmetry' or thePhysicalType == 'empty':
            correctPressureZeroGradient(iBPatch, reg)
        else:
            raise Exception('thePhysicalType not implemented')


def correctVelocityForBoundaryPatches(iComponent,reg):
    for iBPatch in range(reg.mesh.numberOfBoundaryPatches):
        thePhysicalType = reg.mesh.cfdBoundaryPatchesArray[iBPatch].type
        theBCType = reg.fluid.U.boundaryPatchRef[iBPatch].type

        if thePhysicalType == 'wall':
            if theBCType=='noSlip':
                continue
            elif theBCType=='slip':
                correctVelocitySlipWall(iBPatch,iComponent,reg)
            elif theBCType=='fixedValue':
                continue
            else:
                raise Exception('wall condition not implemented')

        elif thePhysicalType == 'inlet':
            if theBCType == 'inlet' or theBCType == 'zeroGradient':
                correctVelocityZeroGradient(iBPatch, iComponent,reg)
            elif theBCType == 'fixedValue':
                continue
            else:
                raise Exception('inlet not implemented')

        elif thePhysicalType == 'outlet':
            if theBCType == 'outlet' or theBCType == 'zeroGradient':
                correctVelocityZeroGradient(iBPatch,iComponent, reg)
            elif theBCType == 'fixedValue':
                continue
            else:
                raise Exception('oulet not implemented')

        elif thePhysicalType == 'symmetry' or thePhysicalType == 'empty' or thePhysicalType == 'symmetryPlane':
            correctVelocitySymmetry(iBPatch,iComponent, reg)
        else:
            raise Exception('thePhysicalType not implemented')

def correctPressureZeroGradient(iBPatch, reg):
    iBstartface = reg.mesh.cfdBoundaryPatchesArray[iBPatch].startFaceIndex
    iBendface = iBstartface + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces
    iBstartelement = reg.mesh.numberOfElements + iBstartface - reg.mesh.numberOfInteriorFace
    iBendelement = iBstartelement + reg.mesh.cfdBoundaryPatchesArray[iBPatch].nFaces

    iOwners=reg.mesh.owners[iBstartface:iBendface]
    reg.fluid.p.phi[iBstartelement:iBendelement,:]=reg.fluid.p.phi[iOwners,:]

def correctPressureSpecifiedValue(iBPatch, reg):
    pass