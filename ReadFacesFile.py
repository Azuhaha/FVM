import numpy as np
def readfacesfile(mesh,path_facefile):
    with open(path_facefile,"r") as f:
        lst = f.readlines()
        for iele,ele in enumerate(lst):
            if not ele.strip():
                lst=lst[iele+1:]
                break

        lst = [x for x in lst if x.strip().strip('()/* ')!='']

        mesh.numberOfFaces=int(lst[0])

        data_lst=[]
        for ele in lst[1:]:
            ele=ele[1:].strip().strip('()')
            data_lst.append([int(x) for x in ele.split()])

        mesh.faceNodes=data_lst