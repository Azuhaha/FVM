
def updatePrevIter(reg):
    for fluidfieldname in ['U','p','rho','T','mdot_f','Cp','mu']:
        if fluidfieldname in dir(reg.fluid):
            fluidfield=getattr(reg.fluid,fluidfieldname)
            fluidfield.previterphi=fluidfield.phi
            fluidfield.prevtimestepphi = fluidfield.phi
