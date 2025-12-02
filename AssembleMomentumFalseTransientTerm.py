import numpy as np


def assembleMomentumFalseTransientTerm(iComponent,reg):
    volume = reg.mesh.elementVolumes

    Ui=reg.fluid.U.phi[:reg.mesh.numberOfElements,[iComponent]]
    Ui_old=reg.fluid.U.prevtimestepphi[:reg.mesh.numberOfElements,[iComponent]]

    rho=reg.fluid.rho.phi[:reg.mesh.numberOfElements,:]
    rho_old = reg.fluid.rho.prevtimestepphi[:reg.mesh.numberOfElements, :]

    #非瞬态计算的dt取很大
    fasleDeltaT=1e6

    #rho*V/dt
    reg.fluxes.FluxC=volume*rho/fasleDeltaT
    #-rho'*V/dt
    reg.fluxes.FluxC_old = -1*volume * rho_old / fasleDeltaT
    reg.fluxes.FluxV = np.zeros(reg.fluxes.FluxC.shape)
    #rho*V*U/dt-rho'*V*U'/dt
    reg.fluxes.FluxT=reg.fluxes.FluxC*Ui+reg.fluxes.FluxC_old*Ui_old




