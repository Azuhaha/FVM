import numpy as np
def assembleMomentumTransientTerm(iComponent,reg):
    if reg.foamDictionary.fvSchemes.ddtSchemes.default=='steadyState':
        return
    elif reg.foamDictionary.fvSchemes.ddtSchemes.default=='Euler':
        assembleFirstOrderEulerTransientTerm(iComponent,reg)
    else:
        raise Exception("ddiScheme is incorrect")

#函数assembleMomentumTransientTerm引用
def assembleFirstOrderEulerTransientTerm(iComponent,reg):
    deltat=reg.foamDictionary.controlDict.deltaT
    #ρ*V/dt
    reg.fluxes.FluxC=reg.mesh.elementVolumes*reg.fluid.rho.phi[0:reg.mesh.numberOfElements,:]/deltat
    #-ρ'*V/dt
    reg.fluxes.FluxC_old=-1*reg.mesh.elementVolumes*reg.fluid.rho.prevtimestepphi[0:reg.mesh.numberOfElements,:]/deltat
    reg.fluxes.FluxV=np.zeros(reg.fluxes.FluxC.shape)
    #ρ*V*U/dt-ρ'*V*U'/dt
    reg.fluxes.FluxT=reg.fluxes.FluxC*reg.fluid.U.phi[0:reg.mesh.numberOfElements,iComponent]+reg.fluxes.FluxC_old*reg.fluid.U.prevtimestepphi[0:reg.mesh.numberOfElements,iComponent]