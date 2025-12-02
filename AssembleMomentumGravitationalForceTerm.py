
def assembleMomentumGravitationalForceTerm(iComponent,reg):
    volume = reg.mesh.elementVolumes

    rho = reg.fluid.rho.phi[:reg.mesh.numberOfElements, :]

    #待完善...
    #gi= reg.foamDictionary.g.value[iComponent]

    #reg.fluxes.FluxT=rho*gi*volume


