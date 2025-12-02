import numpy as np
def readownerfile(mesh,path_ownerfile):
    with open(path_ownerfile,"r") as f:
        lst = f.readlines()
        for iele,ele in enumerate(lst):
            if not ele.strip():
                lst=lst[iele+1:]
                break

        lst = [x for x in lst if x.strip().strip('()/* ')!='']

        data_lst=[int(x) for x in lst[1:]]

        mesh.owners=data_lst
        mesh.numberOfElements=max(data_lst)+1