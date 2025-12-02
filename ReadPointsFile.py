import numpy as np

def readpointsfile(mesh,path_pointfile):
    with open(path_pointfile,"r") as f:
        lst = f.readlines()
        for iele,ele in enumerate(lst):
            if not ele.strip():
                lst=lst[iele+1:]
                break

        lst = [x.strip().strip('()') for x in lst if x.strip().strip('()/* ')!='']

        mesh.numberOfNodes=int(lst[0])

        data_lst=[]
        for ele in lst[1:]:
            data_lst.append([float(x) for x in ele.split()])

        mesh.nodeCentroids=np.array(data_lst)





