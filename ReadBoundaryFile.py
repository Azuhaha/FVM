from FVM.Class.MeshClass import *

def readboundaryfile(mesh,path_boundaryfile):
    with open(path_boundaryfile,"r") as f:
        i=0
        for line in f:
            if not line.strip():
                data = f.readline().strip()
                while not data:
                    data = f.readline().strip()
                mesh.numberOfBoundaryPatches=int(data)

                data=f.readline().strip()
                while '(' not in data:
                    data = f.readline().strip()

                boundary = Boundary()
                boundary.name = f.readline().strip()
                index=0
                for ln in f:
                    if 'type' in ln:
                        boundary.type=ln.strip().strip(';').split()[1]
                    if 'nFaces' in ln:
                        boundary.nFaces=int(ln.strip().strip(';').split()[1])
                    if 'startFace' in ln:
                        boundary.startFaceIndex=int(ln.strip().strip(';').split()[1])
                    if '}' in ln:
                        boundary.index=index
                        mesh.cfdBoundaryPatchesArray[index]=boundary
                        index+=1
                        boundary = Boundary()
                        boundary.name= f.readline().strip()
                break
