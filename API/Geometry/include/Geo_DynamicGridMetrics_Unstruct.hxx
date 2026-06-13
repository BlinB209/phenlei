inline void Geo_DynamicGridMetrics_Unstruct::SetFaceVelocityX(RDouble *xfv)
{
    this->xfv = xfv;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetFaceVelocityY(RDouble *yfv)
{
    this->yfv = yfv;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetFaceVelocityZ(RDouble *zfv)
{
    this->zfv = zfv;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetCellVelocityX(RDouble *xcv)
{
    this->xcv = xcv;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetCellVelocityY(RDouble *ycv)
{
    this->ycv = ycv;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetCellVelocityZ(RDouble *zcv)
{
    this->zcv = zcv;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetFaceNormalVelocity(RDouble *vgn)
{
    this->vgn = vgn;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetCellVolumeOld(RDouble *voln)
{
    this->voln = voln;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetFaceVelocityX() const
{
    return xfv;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetFaceVelocityY() const
{
    return yfv;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetFaceVelocityZ() const
{
    return zfv;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetCellVelocityX() const
{
    return xcv;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetCellVelocityY() const
{
    return ycv;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetCellVelocityZ() const
{
    return zcv;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetFaceNormalVelocity() const
{
    return vgn;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetCellVolumeOld() const
{
    return voln;
}