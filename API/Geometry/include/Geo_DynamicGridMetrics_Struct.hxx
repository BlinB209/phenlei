inline void Geo_DynamicGridMetrics_Struct::DeallocateAll()
{
    if (voln != NULL) { delete voln;   voln = NULL; }
    if (xfv  != NULL) { delete xfv;    xfv  = NULL; }
    if (yfv  != NULL) { delete yfv;    yfv  = NULL; }
    if (zfv  != NULL) { delete zfv;    zfv  = NULL; }
    if (xcv  != NULL) { delete xcv;    xcv  = NULL; }
    if (ycv  != NULL) { delete ycv;    ycv  = NULL; }
    if (zcv  != NULL) { delete zcv;    zcv  = NULL; }
    if (vgn  != NULL) { delete vgn;    vgn  = NULL; }
}

inline RDouble3D * Geo_DynamicGridMetrics_Struct::GetCellVolumeOld() const
{
    return voln;
}

inline RDouble4D * Geo_DynamicGridMetrics_Struct::GetFaceNormalVelocity() const
{
    return vgn;
}

inline RDouble4D * Geo_DynamicGridMetrics_Struct::GetFaceVelocityX() const
{
    return xfv;
}

inline RDouble4D * Geo_DynamicGridMetrics_Struct::GetFaceVelocityY() const
{
    return yfv;
}

inline RDouble4D * Geo_DynamicGridMetrics_Struct::GetFaceVelocityZ() const
{
    return zfv;
}

inline RDouble3D *Geo_DynamicGridMetrics_Struct::GetCellVelocityX() const
{
    return xcv;
}

inline RDouble3D *Geo_DynamicGridMetrics_Struct::GetCellVelocityY() const
{
    return ycv;
}

inline RDouble3D *Geo_DynamicGridMetrics_Struct::GetCellVelocityZ() const
{
    return zcv;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceVelocityX(RDouble4D *xfv)
{
    if (this->xfv && this->xfv != xfv)
    {
        delete this->xfv;
    }
    this->xfv = xfv;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceVelocityY(RDouble4D *yfv)
{
    if (this->yfv && this->yfv != yfv)
    {
        delete this->yfv;
    }
    this->yfv = yfv;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceVelocityZ(RDouble4D *zfv)
{
    if (this->zfv && this->zfv != zfv)
    {
        delete this->zfv;
    }
    this->zfv = zfv;
}

inline void Geo_DynamicGridMetrics_Struct::SetCellVelocityX(RDouble3D *xcv)
{
    if (this->xcv && this->xcv != xcv)
    {
        delete this->xcv;
    }
    this->xcv = xcv;
}

inline void Geo_DynamicGridMetrics_Struct::SetCellVelocityY(RDouble3D *ycv)
{
    if (this->ycv && this->ycv != ycv)
    {
        delete this->ycv;
    }
    this->ycv = ycv;
}

inline void Geo_DynamicGridMetrics_Struct::SetCellVelocityZ(RDouble3D *zcv)
{
    if (this->zcv && this->zcv != zcv)
    {
        delete this->zcv;
    }
    this->zcv = zcv;
}

inline void Geo_DynamicGridMetrics_Struct::SetCellVolumeOld(RDouble3D *voln)
{
    if (this->voln && this->voln != voln)
    {
        delete this->voln;
    }
    this->voln = voln;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceNormalVelocity(RDouble4D *vgn)
{
    if (this->vgn && this->vgn != vgn)
    {
        delete this->vgn;
    }
    this->vgn = vgn;
}