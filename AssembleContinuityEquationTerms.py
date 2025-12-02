from ZeroFaceFLUXCoefficients import *
from AssembleIntoGlobalMatrixFaceFluxes import *
from AssembleMassDivergenceTerm import *
from AssembleMassDivergenceAdvectionTerm import *

def assembleContinuityEquationTerms(reg):
    for termName in reg.model.p.terms:
        if termName=='massDivergenceTerm':
            zeroFaceFLUXCoefficients(reg)
            assembleMassDivergenceTerm(reg)

            if reg.compressible:
                assembleMassDivergenceAdvectionTerm(reg)

            storeMassFlowRate(reg)

            assembleIntoGlobalMatrixFaceFluxes(reg)

        elif termName=='Transient':
            pass
        elif termName=='FalseTransient':
            pass
        else:
            raise Exception('the TermName is not defined')


def storeMassFlowRate(reg):
    #interior
    #?
    reg.fluxes.FluxTf[:reg.mesh.numberOfInteriorFace,:]=reg.fluxes.FluxCf[:reg.mesh.numberOfInteriorFace,:]*reg.fluid.p.phi[reg.mesh.owners[:reg.mesh.numberOfInteriorFace]]+reg.fluxes.FluxFf[:reg.mesh.numberOfInteriorFace,:]*reg.fluid.p.phi[reg.mesh.neighbours]+reg.fluxes.FluxVf[:reg.mesh.numberOfInteriorFace,:]
    #boundary
    iBFaces=range(reg.mesh.numberOfInteriorFace,reg.mesh.numberOfFaces)
    iBElements=range(reg.mesh.numberOfElements,reg.mesh.numberOfElements+reg.mesh.numberOfBElements)
    reg.fluxes.FluxTf[iBFaces,:]=reg.fluxes.FluxCf[iBFaces,:]*reg.fluid.p.phi[reg.mesh.owners[reg.mesh.numberOfInteriorFace:reg.mesh.numberOfFaces],:]+reg.fluxes.FluxFf[iBFaces,:]*reg.fluid.p.phi[iBElements,:]+reg.fluxes.FluxVf[iBFaces,:]

    reg.fluid.mdot_f.phi=reg.fluxes.FluxTf