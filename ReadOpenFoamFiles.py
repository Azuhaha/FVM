from ReadPolyMesh import *
from ProcessTopology import *
from ReadSystem import *
from ReadTimeDirectory import *
from ReadTransportProperties import *
from ReadThermophysicalProperties import *
from ReadTurbulenceProperties import *
from ReadGravity import *

def readopenfoamfiles(region):

    # 读取polymesh文件
    readpolymesh(region.mesh, region.caseDirectoryPath)

    # 拓扑
    processtopology(region.mesh)


    # 读取system内容
    readsystem(region.foamDictionary, region.caseDirectoryPath)


    # 读取初始化场
    readtimedirectory(region)


    # 读取transportproperties
    readtransportproperties(region)


    # 读取thermophysicalproperties
    readthermophysicalproperties(region)


    # 读取turbulenceproperties
    readturbulenceproperties(region)

    # 读取gravity
    readgravity(region)