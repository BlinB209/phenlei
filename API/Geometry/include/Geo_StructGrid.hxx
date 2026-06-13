inline int StructGrid::GetMultigridStepI() const
{
    return multigridInfo->GetMultigridStepI();
}

inline int StructGrid::GetMultigridStepJ() const
{
    return multigridInfo->GetMultigridStepJ();
}

inline int StructGrid::GetMultigridStepK() const
{
    return multigridInfo->GetMultigridStepK();
}

inline void StructGrid::SetMultigridStepI(int istp)
{
    multigridInfo->SetMultigridStepI(istp);
}

inline void StructGrid::SetMultigridStepJ(int jstp)
{
    multigridInfo->SetMultigridStepJ(jstp);
}

inline void StructGrid::SetMultigridStepK(int kstp)
{
    multigridInfo->SetMultigridStepK(kstp);
}

inline RDouble3D * StructGrid::GetWallDist() const
{
    return walldist;
}

inline RDouble3D * StructGrid::GetNearestWallFaceNormalX() const
{
    return nearestwallfacenormalx;
}

inline RDouble3D * StructGrid::GetNearestWallFaceNormalY() const
{
    return nearestwallfacenormaly;
}

inline RDouble3D * StructGrid::GetNearestWallFaceNormalZ() const
{
    return nearestwallfacenormalz;
}

inline Int3D *StructGrid::GetCellBoundaryType()
{
    return cellBoundaryType;
}

inline RDouble3D * StructGrid::GetCellCenterX() const
{
    return cellMetrics->GetCellCenterX();
}

inline RDouble3D * StructGrid::GetCellCenterY() const
{
    return cellMetrics->GetCellCenterY();
}

inline RDouble3D * StructGrid::GetCellCenterZ() const
{
    return cellMetrics->GetCellCenterZ();
}

inline RDouble4D * StructGrid::GetFaceNormalX() const
{
    return faceMetrics->GetFaceNormalX();
}

inline RDouble4D * StructGrid::GetFaceNormalY() const
{
    return faceMetrics->GetFaceNormalY();
}

inline RDouble4D * StructGrid::GetFaceNormalZ() const
{
    return faceMetrics->GetFaceNormalZ();
}

inline RDouble4D * StructGrid::GetFaceArea() const
{
    return faceMetrics->GetFaceArea();
}

inline RDouble4D * StructGrid::GetFaceVectorX() const
{
    return faceMetrics->GetFaceVectorX();
}

inline RDouble4D * StructGrid::GetFaceVectorY() const
{
    return faceMetrics->GetFaceVectorY();
}

inline RDouble4D * StructGrid::GetFaceVectorZ() const
{
    return faceMetrics->GetFaceVectorZ();
}

inline RDouble3D * StructGrid::GetCellVolume() const
{
    return cellMetrics->GetCellVolume();
}

inline RDouble3D * StructGrid::GetCellLengthX() const
{
    return cellMetrics->GetCellLengthX();
}

inline RDouble3D * StructGrid::GetCellLengthY() const
{
    return cellMetrics->GetCellLengthY();
}

inline RDouble3D * StructGrid::GetCellLengthZ() const
{
    return cellMetrics->GetCellLengthZ();
}

inline RDouble3D * StructGrid::GetCellVolumeOld() const
{
    return dynamicGridMetrics->GetCellVolumeOld();
}

inline RDouble3D *StructGrid::GetCellVelocityX() const
{
    return dynamicGridMetrics->GetCellVelocityX();
}

inline RDouble3D *StructGrid::GetCellVelocityY() const
{
    return dynamicGridMetrics->GetCellVelocityY();
}

inline RDouble3D *StructGrid::GetCellVelocityZ() const
{
    return dynamicGridMetrics->GetCellVelocityZ();
}

inline RDouble5D * StructGrid::GetFaceVector_FD() const
{
    return faceMetrics->GetFaceVector_FD();
}

inline RDouble5D * StructGrid::GetFaceNormal_FD() const
{
    return faceMetrics->GetFaceNormal_FD();
}

inline RDouble3D * StructGrid::GetCellJacobian() const
{
    return cellMetrics->GetCellJacobian();
}

inline RDouble4D * StructGrid::GetFaceNormalVelocity() const
{
    return dynamicGridMetrics->GetFaceNormalVelocity();
}

inline void StructGrid::SetNI(int ni)
{
    nodeTopology->SetNI(ni);
}

inline void StructGrid::SetNJ(int nj)
{
    nodeTopology->SetNJ(nj);
}

inline void StructGrid::SetNK(int nk)
{
    nodeTopology->SetNK(nk);
}

inline int StructGrid::GetNI() const
{
    return nodeTopology->GetNI();
}

inline int StructGrid::GetNJ() const
{
    return nodeTopology->GetNJ();
}

inline int StructGrid::GetNK() const
{
    return nodeTopology->GetNK();
}

inline void StructGrid::SetZoneProbesCellNI(vector <int> zoneProbesCellNI)
{
    this->zoneProbesCellNI = zoneProbesCellNI;
}

inline void StructGrid::SetZoneProbesCellNJ(vector <int> zoneProbesCellNJ)
{
    this->zoneProbesCellNJ = zoneProbesCellNJ;
}

inline void StructGrid::SetZoneProbesCellNK(vector <int> zoneProbesCellNK)
{
    this->zoneProbesCellNK = zoneProbesCellNK;
}

inline vector <int> StructGrid::GetZoneProbesCellNI() const
{
    return zoneProbesCellNI;
}

inline vector <int> StructGrid::GetZoneProbesCellNJ() const
{
    return zoneProbesCellNJ;
}

inline vector <int> StructGrid::GetZoneProbesCellNK() const
{
    return zoneProbesCellNK;
}

inline void StructGrid::GetCellIterationIndex(int &iCellStart, int &iCellEnd, int &jCellStart, int &jCellEnd, int &kCellStart, int &kCellEnd, int layer)
{
    int ni = nodeTopology->GetNI();
    int nj = nodeTopology->GetNJ();
    int nk = nodeTopology->GetNK();

    iCellStart = 1 - layer;
    iCellEnd = ni - 1 + layer;
    jCellStart = 1 - layer;
    jCellEnd = nj - 1 + layer;

    if (nk == 1)
    {
        kCellStart = 1;
        kCellEnd = 1;
    }
    else
    {
        kCellStart = 1 - layer;
        kCellEnd = nk -1 + layer;
    }
}

inline void StructGrid::GetNodeIterationIndex(int &iNodeStart, int &iNodeEnd, int &jNodeStart, int &jNodeEnd, int &kNodeStart, int &kNodeEnd)
{
    int ni = nodeTopology->GetNI();
    int nj = nodeTopology->GetNJ();
    int nk = nodeTopology->GetNK();

    iNodeStart = 1;
    iNodeEnd = ni;
    jNodeStart = 1;
    jNodeEnd = nj;

    if (nk == 1)
    {
        kNodeStart = 1;
        kNodeEnd = 1;
    }
    else
    {
        kNodeStart = 1;
        kNodeEnd = nk;
    }
}

inline void StructGrid::GetFaceIterationIndex(int &iFaceStart, int &iFaceEnd, int &jFaceStart, int &jFaceEnd, int &kFaceStart, int &kFaceEnd, int iSurface)
{
    int ni = nodeTopology->GetNI();
    int nj = nodeTopology->GetNJ();
    int nk = nodeTopology->GetNK();

    if (iSurface == 1)
    {
        iFaceStart = 0;
        iFaceEnd = ni-1;
        jFaceStart = 1;
        jFaceEnd = nj-1;

        if (nk == 1)
        {
            kFaceStart = 1;
            kFaceEnd = 1;
        }
        else
        {
            kFaceStart = 1;
            kFaceEnd = nk-1;
        }
    }
    else if (iSurface == 2)
    {
        iFaceStart = 1;
        iFaceEnd = ni-1;
        jFaceStart = 0;
        jFaceEnd = nj-1;

        if (nk == 1)
        {
            kFaceStart = 1;
            kFaceEnd = 1;
        }
        else
        {
            kFaceStart = 1;
            kFaceEnd = nk-1;
        }
    }
    else
    {
        iFaceStart = 1;
        iFaceEnd = ni-1;
        jFaceStart = 1;
        jFaceEnd = nj-1;
        kFaceStart = 0;
        kFaceEnd = nk-1;
    }
}

inline void StructGrid::GetNsurfIndex(int &il1, int &jl1, int &kl1, int iDimension)
{
    if(iDimension == 1)
    {
        il1 = 1;
        jl1 = 0;
        kl1 = 0;
    }
    else if(iDimension == 2)
    {
        il1 = 0;
        jl1 = 1;
        kl1 = 0;
    }
    else if(iDimension == 3)
    {
        il1 = 0;
        jl1 = 0;
        kl1 = 1;
    }
}

inline void StructGrid::GetLeftCellOfFace(int &i, int &j, int &k, int &il1, int &jl1, int &kl1, int &ile, int &jle, int &kle)
{
	ile  = i - il1;
	jle  = j - jl1;
	kle  = k - kl1;
}

inline void StructGrid::GetRightCellOfFace(int &i, int &j, int &k, int &il1, int &jl1, int &kl1, int &ire, int &jre, int &kre)
{
	ire  = i + il1;
	jre  = j + jl1;
	kre  = k + kl1;
}

inline void StructGrid::GetND(int &si, int &sj, int &sk) const
{
    si = nodeTopology->GetNI();
    sj = nodeTopology->GetNJ();
    sk = nodeTopology->GetNK();
}

inline void StructGrid::GetLocalCenter(int gp, int &ip, int &jp, int &kp) const
{
    int ni = nodeTopology->GetNI();
    int nj = nodeTopology->GetNJ();
    kp = gp / ((ni - 1) * (nj - 1));

    int other = gp % ((ni - 1) * (nj - 1));

    jp = other / (ni - 1);

    ip = other % (ni - 1);

    return;
}

inline RDouble3D * StructGrid::GetStructX() const
{
    return nodeTopology->GetStructX();
}

inline RDouble3D * StructGrid::GetStructY() const
{
    return nodeTopology->GetStructY();
}

inline RDouble3D * StructGrid::GetStructZ() const
{
    return nodeTopology->GetStructZ();
}

inline StructBCSet * StructGrid::GetStructBCSet() const
{
    return structBCSet;
}

inline RDouble4D * StructGrid::GetFaceVelocityX() const
{
    return dynamicGridMetrics->GetFaceVelocityX();
}

inline RDouble4D * StructGrid::GetFaceVelocityY() const
{
    return dynamicGridMetrics->GetFaceVelocityY();
}

inline RDouble4D * StructGrid::GetFaceVelocityZ() const
{
    return dynamicGridMetrics->GetFaceVelocityZ();
}

inline int StructGrid::GetZoneStartPointLabel() const
{
    return oversetGridTopology->GetZoneStartPointLabel();
}

inline void StructGrid::SetZoneStartPointLabel(int zoneStartPointLabel)
{
    oversetGridTopology->SetZoneStartPointLabel(zoneStartPointLabel);
}

inline int * StructGrid::GetHingedPointContainer() const
{
    return oversetGridTopology->GetHingedPointContainer();
}

inline RDouble * StructGrid::GetXCoreContainer() const
{
    return oversetGridTopology->GetXCoreContainer();
}

inline RDouble * StructGrid::GetYCoreContainer() const
{
    return oversetGridTopology->GetYCoreContainer();
}

inline RDouble * StructGrid::GetZCoreContainer() const
{
    return oversetGridTopology->GetZCoreContainer();
}

inline int StructGrid::GetNumberOfCores() const
{
    return oversetGridTopology->GetNumberOfCores();
}

inline int * StructGrid::GetILinkPointLabel() const
{
    return oversetGridTopology->GetILinkPointLabel();
}

inline int * StructGrid::GetJLinkPointLabel() const
{
    return oversetGridTopology->GetJLinkPointLabel();
}

inline int * StructGrid::GetKLinkPointLabel() const
{
    return oversetGridTopology->GetKLinkPointLabel();
}

inline void StructGrid::SetZoneStartCenterLabel(int zoneStartCenterLabel)
{
    oversetGridTopology->SetZoneStartCenterLabel(zoneStartCenterLabel);
}