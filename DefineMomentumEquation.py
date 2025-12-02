import os
import numpy as np
def definemomentumequation(reg):

    if not os.path.exists(os.path.join(reg.caseDirectoryPath,'0','U')):
        raise Exception("0文件夹中U文件不存在")
    else:
        #print('Defining Momentum Equation...')
        #初始化
        reg.model.U.residuals.rmsResidual=np.array([1.,1.,1.])
        reg.model.U.residuals.maxResidual = np.array([1.,1.,1.])
        reg.model.U.residuals.initResidual = np.array([1.,1.,1.])
        reg.model.U.residuals.finalResidual= np.array([1.,1.,1.])

        #terms增加项
        reg.model.U.terms=['Convection','Stress','PressureGradient']

        if reg.STEADY_STATE_RUN:
            reg.model.U.terms.append('FalseTransient')
        else:
            reg.model.U.terms.append('Transient')

        if reg.foamDictionary.g.value!=None:
            reg.model.U.terms.append('GravitationalForce')
