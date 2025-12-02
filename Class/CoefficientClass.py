class Coefficient:
    def __init__(self):
        self.ac=None
        self.bc=None
        self.ac_old=None
        self.anb=None
        self.dc=None
        self.rc=None
        self.dphi=None
        self.ccon=None #相邻网格的编号
        self.csize=None #相邻网格的数量
        self.numOfElements=None
        self.parents=None