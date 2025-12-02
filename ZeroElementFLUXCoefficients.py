import numpy as np
def zeroElementFLUXCoefficients(reg):
    reg.fluxes.FluxC = np.zeros((reg.mesh.numberOfElements,1))
    reg.fluxes.FluxV = np.zeros((reg.mesh.numberOfElements, 1))
    reg.fluxes.FluxT = np.zeros((reg.mesh.numberOfElements, 1))
    reg.fluxes.FluxC_old = np.zeros((reg.mesh.numberOfElements, 1))