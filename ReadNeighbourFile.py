import numpy as np
def readneighbourfile(mesh,path_neighbourfile):
    with open(path_neighbourfile,"r") as f:
        lst = f.readlines()
        for iele,ele in enumerate(lst):
            if not ele.strip():
                lst=lst[iele+1:]
                break

        lst = [x for x in lst if x.strip().strip('()/* ')!='']

        mesh.numberOfNeighbour=int(lst[0])

        data_lst=[int(x) for x in lst[1:]]

        mesh.neighbours=data_lst


    mesh.numberOfInteriorFace=mesh.numberOfNeighbour
    mesh.numberOfBFaces=mesh.numberOfFaces-mesh.numberOfInteriorFace
    mesh.numberOfBElements=mesh.numberOfBFaces