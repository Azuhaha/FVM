class Model:
    def __init__(self):
        self.U = ModelVar('U')
        self.p = ModelVar('p')

class ModelVar:
    def __init__(self,name):
        self.name=name
        self.residuals=Residuals()
        self.terms=None

class Residuals:
    def __init__(self):
        self.rmsResidual=None
        self.maxResidual=None
        self.initResidual = None
        self.finalResidual = None