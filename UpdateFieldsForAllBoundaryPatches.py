from UpdateVectorFieldForAllBoundaryPatches import *
from UpdateScalarFieldForAllBoundaryPatches import *

def updateFieldsForAllBoundaryPatches(reg):
    if 'U' in dir(reg.fluid):
        updateVectorFieldForAllBoundaryPatches(reg,'U')
    elif 'p' in dir(reg.fluid):
        updateScalarFieldForAllBoundaryPatches(reg,'p')
    elif 'T' in dir(reg.fluid):
        updateScalarFieldForAllBoundaryPatches(reg,'T')
    elif 'rho' in dir(reg.fluid):
        updateScalarFieldForAllBoundaryPatches(reg,'rho')
    elif 'mu' in dir(reg.fluid):
        updateScalarFieldForAllBoundaryPatches(reg,'mu')
    elif 'Cp' in dir(reg.fluid):
        updateScalarFieldForAllBoundaryPatches(reg,'Cp')
    elif 'k' in dir(reg.fluid):
        updateScalarFieldForAllBoundaryPatches(reg,'k')

