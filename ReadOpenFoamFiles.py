from FVM.ReadPolyMesh import *
from FVM.ProcessTopology import *
from FVM.ReadSystem import *
from FVM.ReadTimeDirectory import *
from FVM.ReadTransportProperties import *
from FVM.ReadThermophysicalProperties import *
from FVM.ReadTurbulenceProperties import *
from FVM.ReadGravity import *

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