
from Class.RegionClass import Region
from ReadOpenFoamFiles import *
from DefineMomentumEquation import *
from DefineContinuityEquation import *
from RunCase import *
from PlotMesh import *

casepath="./elbow"

isPlotMesh=False

region=Region()
region.caseDirectoryPath=casepath

#读取openfoam文件
readopenfoamfiles(region)

#显示网格
if isPlotMesh:
    plotMesh(region)

#定义region.model.U
definemomentumequation(region)

#定义region.model.p; 初始化mdot_f=rho_f*U_f*Sf; 初始化DU1 DU2 DU3 DUT1 DUT2 DUT3 pp为0值
definecontinuityequation(region)

#开始运行
runcase(region)

print('finished')


