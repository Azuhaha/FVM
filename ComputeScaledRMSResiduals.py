import numpy as np
import math
def computeScaledRMSResiduals(equationname, iComponent,reg):
    scale=getattr(reg.fluid,equationname).scale

    if equationname=='p':
        #For pressure correction equation, the divergence of the mass flow rate is the residual.
        #Scale with scale value (max value)
        p_scale=getattr(reg.fluid,'p').scale
        maxScaledResidual=np.max(np.abs(reg.coefficient[0].bc)/(reg.coefficient[0].ac*p_scale))
        temp=np.abs(reg.coefficient[0].bc) / (reg.coefficient[0].ac * p_scale)
        rmsScaledResidual=np.sqrt(np.mean(temp**2))

    else:
        if not reg.STEADY_STATE_RUN:
            rho=reg.fluid.rho.phi
            if equationname=='T':
                Cp=reg.fluid.Cp.phi
                rho*=Cp
            deltaT=reg.foamDictionary.controlDict.deltaT

            volume=reg.mesh.elementVolumes

            at=volume*rho[:reg.mesh.numberOfElements]/deltaT
            local_ac=reg.coefficient[0].ac-at

            mask = local_ac < (1e-6*at)
            local_ac[mask]=at[mask]

            local_residual=reg.coefficient[0].bc/(local_ac*scale)

            maxScaledResidual=np.max(np.abs(local_residual))
            idx=np.argmax(np.abs(local_residual))
            maxResidualSquared=maxScaledResidual+local_residual[idx,0]**2

            rmsScaledResidual=math.sqrt(maxResidualSquared/reg.mesh.numberOfElements)

        else:
            local_residual = reg.coefficient[0].bc / (reg.coefficient[0].ac * scale)

            maxScaledResidual=np.max(np.abs(local_residual))
            idx=np.argmax(np.abs(local_residual))
            maxResidualSquared=maxScaledResidual+local_residual[idx,0]**2

            rmsScaledResidual=math.sqrt(maxResidualSquared/reg.mesh.numberOfElements)

    getattr(reg.model,equationname).residuals.rmsResidual[iComponent]=rmsScaledResidual
    getattr(reg.model, equationname).residuals.maxResidual[iComponent] = maxScaledResidual