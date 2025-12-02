import numpy as np
def updatescales(reg):
    if 'U' in dir(reg.fluid):
        updatescale(getattr(reg.fluid,'U'),reg)
    if 'p' in dir(reg.fluid):
        updatescale(getattr(reg.fluid,'p'),reg)
    if 'T' in dir(reg.fluid):
        updatescale(getattr(reg.fluid,'T'),reg)
    if 'rho' in dir(reg.fluid):
        updatescale(getattr(reg.fluid,'rho'),reg)

def updatescale(fluidfield,reg):

    phiMax= np.max(np.linalg.norm(fluidfield.phi, axis=1))
    phiMin =np.min(np.linalg.norm(fluidfield.phi, axis=1))

    if fluidfield.name=='p':
        vel_scale=getattr(reg.fluid,'U').scale
        rho_scale=getattr(reg.fluid,'rho').scale
        p_dyn=0.5*rho_scale*vel_scale**2
        fluidfield.scale=max(p_dyn,phiMax)
    elif fluidfield.name=='U':
        GeoLengthScale=np.sum(reg.mesh.elementVolumes)**(1/3)
        fluidfield.scale=max(GeoLengthScale,phiMax)
    else:
        fluidfield.max=phiMax
        fluidfield.min = phiMin
        fluidfield.scale = phiMax