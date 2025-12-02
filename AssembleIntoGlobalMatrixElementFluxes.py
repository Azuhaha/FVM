
def assembleIntoGlobalMatrixElementFluxes(reg):
    reg.coefficient[0].ac += reg.fluxes.FluxC
    reg.coefficient[0].ac_old += reg.fluxes.FluxC_old
    reg.coefficient[0].bc -= reg.fluxes.FluxT