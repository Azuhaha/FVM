
def assembleIntoGlobalMatrixFaceFluxes(reg):
    #内部面
    for iFace in range(reg.mesh.numberOfInteriorFace):
        own=reg.mesh.owners[iFace]
        nei=reg.mesh.neighbours[iFace]

        iOwnerNeighbourCoef=reg.mesh.upperAnbCoeffIndex[iFace]
        iNeighbourOwnerCoef=reg.mesh.lowerAnbCoeffIndex[iFace]

        reg.coefficient[0].ac[own,:]+=reg.fluxes.FluxCf[iFace,:]
        reg.coefficient[0].anb[own][iOwnerNeighbourCoef,:]+=reg.fluxes.FluxFf[iFace,:]
        #注意bc为残差形式 式14.22,14.23
        reg.coefficient[0].bc[own,:]-=reg.fluxes.FluxTf[iFace,:]

        reg.coefficient[0].ac[nei,:]-=reg.fluxes.FluxFf[iFace,:]
        reg.coefficient[0].anb[nei][iNeighbourOwnerCoef,:]-=reg.fluxes.FluxCf[iFace,:]
        reg.coefficient[0].bc[nei,:]+=reg.fluxes.FluxTf[iFace,:]

    #边界面
    for iFace in range(reg.mesh.numberOfInteriorFace,reg.mesh.numberOfFaces):
        own = reg.mesh.owners[iFace]

        reg.coefficient[0].ac[own, :] += reg.fluxes.FluxCf[iFace, :]
        reg.coefficient[0].bc[own, :] -= reg.fluxes.FluxTf[iFace, :]