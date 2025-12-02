class FoamDictionary:
    def __init__(self):
        self.controlDict=ControlDict()
        self.fvSchemes=FvSchemes()
        self.fvSolution=FvSolution()
        self.transportProperties=TransportProperties()
        self.turbulenceProperties=TurbulenceProperties()
        self.g=Grav()

class ControlDict:
    def __init__(self):
        self.application=None
        self.startFrom=None
        self.startTime=None
        self.stopAt=None
        self.endTime=None
        self.deltaT=None
        self.writeControl = None
        self.writeInterval = None
        self.purgeWrite = None

class FvSchemes:
    def __init__(self):
        self.ddtSchemes=Default()
        self.gradSchemes=Default()
        self.divSchemes=Default()
        self.laplacianSchemes = Default()
        self.interpolationSchemes = Default()
        self.snGradSchemes = Default()

class FvSolution:
    def __init__(self):
        self.solvers=SolverVar()
        self.SIMPLE=SimpleControl()
        self.relaxationFactors=RelaxationFactors()

class TransportProperties:
    def __init__(self):
        self.mu=TransportPropertiesVar()
        self.rho=TransportPropertiesVar()
        self.k=TransportPropertiesVar()
        self.Cp = TransportPropertiesVar()

class TransportPropertiesVar:
    def __init__(self):
        self.name=None
        self.dimensions=[]
        self.propertyValue=None

class TurbulenceProperties:
    def __init__(self):
        self.RASModel=None
        self.turbulence=None
        self.printCoeffs=None

class Grav:
    def __init__(self):
        self.value=None

class Default:
    def __init__(self):
        self.default=None

class SolverVar:
    def __init__(self):
        self.p=Variable()
        self.U = Variable()
        self.T = Variable()

class Variable:
    def __init__(self):
        self.solver=None
        self.preconditioner=None
        self.smoother=None
        self.tolerance=None
        self.relTol=None
        self.nPreSweeps=None
        self.nPostSweeps=None
        self.maxIter=None
        self.nFinestSweeps=None

class SimpleControl:
    def __init__(self):
        self.residualControl=ResidualControl()
        self.nCorrectors=None
        self.pRefCell=None
        self.pRefValue=None

class ResidualControl:
    def __init__(self):
        self.p=None
        self.U=None
        self.T=None

class RelaxationFactors:
    def __init__(self):
        self.equations=RelaxationFactorsEquation()
        self.fields = RelaxationFactorsField()

class RelaxationFactorsEquation:
    def __init__(self):
        self.U=None
        self.T =None

class RelaxationFactorsField:
    def __init__(self):
        self.p=None