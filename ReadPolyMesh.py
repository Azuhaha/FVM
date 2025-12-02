import os
import numpy as np
from ReadPointsFile import *
from ReadFacesFile import *
from ReadOwnerFile import *
from ReadNeighbourFile import *
from ReadBoundaryFile import *
def readpolymesh(mesh,path):
    for file0 in os.listdir(path):
        if "constant" in file0:
            for file1 in os.listdir(path+'/'+file0):
                if "polyMesh" in file1:
                    print('Reading OpenFoam polymesh files...')
                    for file2 in os.listdir(path+'/'+file0+'/'+file1):
                        if "points" in file2:
                            path_points=path+'/'+file0+'/'+file1+'/'+file2
                            readpointsfile(mesh,path_points)
                        if "faces" in file2:
                            path_faces=path+'/'+file0+'/'+file1+'/'+file2
                            readfacesfile(mesh,path_faces)
                        if "owner" in file2:
                            path_owner=path+'/'+file0+'/'+file1+'/'+file2
                            readownerfile(mesh,path_owner)
                        if "neighbour" in file2:
                            path_neighbour=path+'/'+file0+'/'+file1+'/'+file2
                            readneighbourfile(mesh,path_neighbour)
                        if "boundary" in file2:
                            path_boundary=path+'/'+file0+'/'+file1+'/'+file2
                            readboundaryfile(mesh,path_boundary)


    # 输出信息
    print("节点个数：{}".format(mesh.numberOfNodes))
    print("面个数：{}".format(mesh.numberOfFaces))
    print("单元个数：{}".format(mesh.numberOfElements))
    boundaryinfo=""
    for i in range(mesh.numberOfBoundaryPatches):
        boundaryinfo+=mesh.cfdBoundaryPatchesArray[i].name+" "
    print("边界条件：{}".format(boundaryinfo))






