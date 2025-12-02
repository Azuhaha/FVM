import numpy as np
def setupFluxes(reg):
    reg.fluxes.FluxCf=np.zeros((reg.mesh.numberOfFaces,1))
    reg.fluxes.FluxFf=np.zeros((reg.mesh.numberOfFaces,1))
    reg.fluxes.FluxVf=np.zeros((reg.mesh.numberOfFaces,1))
    reg.fluxes.FluxTf=np.zeros((reg.mesh.numberOfFaces,1))

    reg.fluxes.FluxC=np.zeros((reg.mesh.numberOfElements,1))
    reg.fluxes.FluxV=np.zeros((reg.mesh.numberOfElements,1))
    reg.fluxes.FluxT=np.zeros((reg.mesh.numberOfElements,1))
    reg.fluxes.FluxC_old=np.zeros((reg.mesh.numberOfElements,1))
