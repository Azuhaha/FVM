from FVM.Class.MeshClass import Mesh
from FVM.Class.FoamDictionaryClass import FoamDictionary
from FVM.Class.FluidClass import Fluid
from FVM.Class.ModelClass import *
from FVM.Class.FluxClass import *
from FVM.Class.TimeClass import *

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





