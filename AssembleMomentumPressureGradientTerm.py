def assembleMomentumPressureGradientTerm(iComponent,reg):
    volume=reg.mesh.elementVolumes
    p_grad=reg.fluid.p.phiGradient[:reg.mesh.numberOfElements,:]
    #V*Grad_p
    reg.fluxes.FluxT=volume*p_grad[:,iComponent,:]