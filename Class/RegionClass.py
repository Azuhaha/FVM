from Class.MeshClass import Mesh
from Class.FoamDictionaryClass import FoamDictionary
from Class.FluidClass import Fluid
from Class.ModelClass import *
from Class.FluxClass import *
from Class.TimeClass import *

class Region:
    def __init__(self):
        self.caseDirectoryPath=None
        self.STEADY_STATE_RUN=True
        self.mesh=Mesh()
        self.foamDictionary=FoamDictionary()
        self.fluid=Fluid()
        self.compressible=False
        self.model=Model()
        self.coefficient={}
        self.fluxes=Flux()
        self.time=Time()





