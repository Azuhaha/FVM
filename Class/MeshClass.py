class Mesh:
    def __init__(self):
        self.nodeCentroids=None
        self.numberOfNodes=None
        self.faceNodes=None
        self.numberOfFaces=None
        self.owners=None
        self.neighbours=None
        self.numberOfInteriorFace=None
        self.numberOfNeighbour = None
        self.numberOfBFaces=None
        self.numberOfElements=None
        self.numberOfBElements=None
        self.numberOfBoundaryPatches=None
        self.cfdBoundaryPatchesArray={}

        self.elementNeighbours=None
        self.elementFaces=None
        self.elementNodes=None
        self.upperAnbCoeffIndex=None#某个内部面对应的owner单元的单元面局部索引
        self.lowerAnbCoeffIndex=None #某个内部面对应的neighbour单元的单元面局部索引

        self.nodeElements=None
        self.nodeFaces=None

        self.elementCentroids=None
        self.elementVolumes=None
        self.faceCentroid=None
        self.faceSf=None
        self.faceAreas=None
        self.faceWeights=None
        self.faceCF=None
        self.faceCf=None
        self.faceFf=None
        self.wallDist=None
        self.wallDistLimited=None

        self.closed=True


class Boundary:
    def __init__(self):
        self.name=None
        self.index=None
        self.nFaces=None
        self.type=None
        self.startFaceIndex = None