
def updateProperties(reg):
    for fluidfieldname in ['rho','k', 'mu']:
        if fluidfieldname in dir(reg.fluid):
            fluidfield=getattr(reg.fluid,fluidfieldname)
            pass
        #待完善