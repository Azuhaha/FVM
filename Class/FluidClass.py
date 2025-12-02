class Fluid:
    def __init__(self):
        self.U = FluidVar('U')
        self.p = FluidVar('p')
        self.mu = FluidVar('mu')
        self.rho = FluidVar('rho')
        self.k = FluidVar('k')
        self.Cp = FluidVar('Cp')
        self.mdot_f = FluidVar('mdot_f')
        self.DU1 = FluidVar('DU1')
        self.DU2 = FluidVar('DU2')
        self.DU3 = FluidVar('DU3')
        self.DUT1 = FluidVar('DUT1')
        self.DUT2 = FluidVar('DUT2')
        self.DUT3 = FluidVar('DUT3')
        self.pp = FluidVar('pp')

class FluidVar:
    def __init__(self,name):
        self.name=name
        self.type = None
        self.phi = None
        self.previterphi=None
        self.prevtimestepphi = None
        self.dimensions = None
        self.boundaryPatchRef = {}
        self.phiGradient=None
        self.max=None
        self.min=None
        self.scale=None

class BoundaryPatchRef:
    def __init__(self):
        self.name=None
        self.type=None
        self.valueType=None
        self.value=None



