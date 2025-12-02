import numpy as np
def zeroFaceFLUXCoefficients(reg):
    reg.fluxes.FluxCf = np.zeros((reg.mesh.numberOfFaces,1))
    reg.fluxes.FluxFf = np.zeros((reg.mesh.numberOfFaces, 1))
    reg.fluxes.FluxVf = np.zeros((reg.mesh.numberOfFaces, 1))
    reg.fluxes.FluxTf = np.zeros((reg.mesh.numberOfFaces, 1))