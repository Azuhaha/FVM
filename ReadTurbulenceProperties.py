import os
def readturbulenceproperties(reg):
    if os.path.exists(reg.caseDirectoryPath + '/constant/turbulenceProperties'):
        print('Reading OpenFoam TurbulenceProperties file...')
        with open(reg.caseDirectoryPath + '/constant/turbulenceProperties','r' ) as f:
            lst = f.readlines()
            for i, ele in enumerate(lst):
                if not ele.strip():
                    lst = lst[i + 1:]
                    break
            for i, ele in enumerate(lst):
                lst[i]=ele.strip().strip(';/*{}').strip().strip('*/{}')
            lst = [x for x in lst if x != '']
            for ele in lst:
                if ele.split()[0]=='RASModel':
                    reg.foamDictionary.turbulenceProperties.RASModel=ele.split()[1]
                if ele.split()[0]=='turbulence':
                    reg.foamDictionary.turbulenceProperties.turbulence = ele.split()[1]
                if ele.split()[0]=='printCoeffs':
                    reg.foamDictionary.turbulenceProperties.printCoeffs = ele.split()[1]


    else:
        reg.foamDictionary.turbulenceProperties.turbulence='off'
        reg.foamDictionary.turbulenceProperties.RASModel='laminar'
        return

