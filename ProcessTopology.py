import numpy as np
def processtopology(mesh):
    #初始化
    #print('Processing Topology...')
    mesh.elementNeighbours=[[] for _ in range(mesh.numberOfElements)]
    mesh.elementFaces=[[] for _ in range(mesh.numberOfElements)]
    mesh.elementNodes=[[] for _ in range(mesh.numberOfElements)]

    mesh.upperAnbCoeffIndex=[[] for _ in range(mesh.numberOfInteriorFace)]
    mesh.lowerAnbCoeffIndex = [[] for _ in range(mesh.numberOfInteriorFace)]

    mesh.nodeElements=[[] for _ in range(mesh.numberOfNodes)]
    mesh.nodeFaces=[[] for _ in range(mesh.numberOfNodes)]



    for iface in range(mesh.numberOfInteriorFace):

        own=mesh.owners[iface]
        nei=mesh.neighbours[iface]

        mesh.elementNeighbours[own].append(nei)
        mesh.elementNeighbours[nei].append(own)

        mesh.elementFaces[own].append(iface)
        mesh.elementFaces[nei].append(iface)

    for iface in range(mesh.numberOfInteriorFace,mesh.numberOfFaces):
        own=mesh.owners[iface]
        mesh.elementFaces[own].append(iface)

    for ielement in range(mesh.numberOfElements):
        for iface in mesh.elementFaces[ielement]:
            mesh.elementNodes[ielement].extend(mesh.faceNodes[iface])
    for ielement,val_lst in enumerate(mesh.elementNodes):
        lst=list(set(val_lst))
        mesh.elementNodes[ielement]=lst

    for ielement in range(mesh.numberOfElements):
        iNb=0
        for iface in mesh.elementFaces[ielement]:
            if iface>=mesh.numberOfInteriorFace:
                continue
            own = mesh.owners[iface]
            nei = mesh.neighbours[iface]
            if ielement==own:
                mesh.upperAnbCoeffIndex[iface].clear()
                mesh.upperAnbCoeffIndex[iface].append(iNb)
            elif ielement==nei:
                mesh.lowerAnbCoeffIndex[iface].clear()
                mesh.lowerAnbCoeffIndex[iface].append(iNb)
            iNb+=1
    mesh.upperAnbCoeffIndex=[x[0] for x in mesh.upperAnbCoeffIndex]
    mesh.lowerAnbCoeffIndex = [x[0] for x in mesh.lowerAnbCoeffIndex]


    for ielement,lst in enumerate(mesh.elementNodes):
        for inode in lst:
            mesh.nodeElements[inode].append(ielement)

    for iface,lst in enumerate(mesh.faceNodes):
        for inode in lst:
            mesh.nodeFaces[inode].append(iface)

    #GeometryProcess
    #面形心
    lst_faceCentroid=[]
    lst_faceSf=[]
    lst_faceAreas=[]
    for iface in range(mesh.numberOfFaces):
        local_centre=np.array([0.,0.,0.])
        for inode in mesh.faceNodes[iface]:
            local_centre+=mesh.nodeCentroids[inode,:]

        local_centre=local_centre/len(mesh.faceNodes[iface])
        centroid=np.array([0.,0.,0.])
        Sf = np.array([0., 0., 0.])
        area=0
        for i in range(len(mesh.faceNodes[iface])):
            point1=local_centre
            point2=mesh.nodeCentroids[mesh.faceNodes[iface][i],:]
            if i<len(mesh.faceNodes[iface])-1:
                point3=mesh.nodeCentroids[mesh.faceNodes[iface][i+1],:]
            else:
                point3 = mesh.nodeCentroids[mesh.faceNodes[iface][0],:]
            local_centroid=(point1+point2+point3)/3.
            local_Sf=0.5*np.cross(point2-point1,point3-point1)
            local_area=np.linalg.norm(local_Sf)

            centroid += local_area * local_centroid
            Sf += local_Sf
            area += local_area
        centroid=centroid/area
        lst_faceCentroid.append(centroid.tolist())
        lst_faceSf.append(Sf.tolist())
        lst_faceAreas.append([area])
    mesh.faceCentroid=np.array(lst_faceCentroid)
    mesh.faceSf=np.array(lst_faceSf)
    mesh.faceAreas=np.array(lst_faceAreas)

    lst_elementCentroids=[]
    lst_elementVolumes=[]
    for ielement in range(mesh.numberOfElements):
        local_centre=np.array([0.,0.,0.])
        for iface in mesh.elementFaces[ielement]:
            local_centre+=mesh.faceCentroid[iface,:]
        local_centre=local_centre/len(mesh.elementFaces[ielement])

        localVolumeCentroidSum=np.array([0.,0.,0.])
        localVolumeSum=0
        for iface in mesh.elementFaces[ielement]:
            Cf=mesh.faceCentroid[iface,:]-local_centre

            faceSign = -1
            if ielement == mesh.owners[iface]:
                faceSign = 1
            local_Sf = faceSign * mesh.faceSf[iface]

            localVolume = np.dot(local_Sf, Cf)/3.

            localCentroid = 0.75 * mesh.faceCentroid[iface,:] + 0.25 * local_centre

            localVolumeCentroidSum += localCentroid * localVolume

            localVolumeSum += localVolume
        lst_elementCentroids.append(localVolumeCentroidSum / localVolumeSum)
        lst_elementVolumes.append([localVolumeSum])
    mesh.elementCentroids=np.array(lst_elementCentroids)
    mesh.elementVolumes=np.array(lst_elementVolumes)

    lst_faceCF=[]
    lst_faceCf = []
    lst_faceFf=[]
    lst_faceWeights=[]
    lst_wallDist=[]
    lst_wallDistLimited=[]
    for iface in range(mesh.numberOfInteriorFace):
        own=mesh.owners[iface]
        nei=mesh.neighbours[iface]
        lst_faceCF.append(mesh.elementCentroids[nei,:]-mesh.elementCentroids[own,:])
        lst_faceCf.append(mesh.faceCentroid[iface,:]-mesh.elementCentroids[own,:])
        lst_faceFf.append(mesh.faceCentroid[iface,:]-mesh.elementCentroids[nei,:])
        n=mesh.faceSf[iface,:]/np.linalg.norm(mesh.faceSf[iface,:])
        lst_faceWeights.append([np.dot(lst_faceCf[iface],n)/(np.dot(lst_faceCf[iface],n)-np.dot(lst_faceFf[iface],n))])
        lst_wallDist.append([0.])
        lst_wallDistLimited.append([0.])
    mesh.faceCF=np.array(lst_faceCF)
    mesh.faceCf=np.array(lst_faceCf)
    mesh.faceFf=np.array(lst_faceFf)
    mesh.faceWeights=np.array(lst_faceWeights)
    mesh.wallDist=np.array(lst_wallDist)
    mesh.wallDistLimited=np.array(lst_wallDistLimited)



    lst_faceCF=[]
    lst_faceCf = []
    lst_faceFf=[]
    lst_faceWeights=[]
    lst_wallDist=[]
    lst_wallDistLimited=[]
    for iface in range(mesh.numberOfInteriorFace,mesh.numberOfFaces):
        own=mesh.owners[iface]
        lst_faceCF.append(mesh.faceCentroid[iface,:]-mesh.elementCentroids[own,:])
        lst_faceCf.append(mesh.faceCentroid[iface,:]-mesh.elementCentroids[own,:])
        lst_faceWeights.append([1.])
        n = mesh.faceSf[iface,:] / np.linalg.norm(mesh.faceSf[iface,:])
        lst_wallDist.append([max(np.dot(lst_faceCf[iface-mesh.numberOfInteriorFace],n),1e-24)])
        lst_wallDistLimited.append([max(lst_wallDist[iface-mesh.numberOfInteriorFace][0],0.05*np.linalg.norm(lst_faceCf[iface-mesh.numberOfInteriorFace]))])
        lst_faceFf.append([0.,0.,0.])

    mesh.faceCF=np.append(mesh.faceCF,lst_faceCF,axis=0)
    mesh.faceCf = np.append(mesh.faceCf, lst_faceCf, axis=0)
    mesh.faceFf = np.append(mesh.faceFf, lst_faceFf, axis=0)
    mesh.faceWeights = np.append(mesh.faceWeights,lst_faceWeights, axis=0)
    mesh.wallDist=np.vstack((mesh.wallDist, np.array(lst_wallDist)))
    mesh.wallDistLimited=np.vstack((mesh.wallDistLimited, np.array(lst_wallDistLimited)))


