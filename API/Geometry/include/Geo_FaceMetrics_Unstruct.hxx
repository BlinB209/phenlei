inline void Geo_FaceMetrics_Unstruct::SetFaceCenterX(RDouble *xfc)
{
    if (this->xfc && this->xfc != xfc)
    {
        delete [] this->xfc;
    }
    this->xfc = xfc;
}

inline void Geo_FaceMetrics_Unstruct::SetFaceCenterY(RDouble *yfc)
{
    if (this->yfc && this->yfc != yfc)
    {
        delete [] this->yfc;
    }
    this->yfc = yfc;
}

inline void Geo_FaceMetrics_Unstruct::SetFaceCenterZ(RDouble *zfc)
{
    if (this->zfc && this->zfc != zfc)
    {
        delete [] this->zfc;
    }
    this->zfc = zfc;
}

inline void Geo_FaceMetrics_Unstruct::SetFaceNormalX(RDouble *xfn)
{
    if (this->xfn && this->xfn != xfn)
    {
        delete [] this->xfn;
    }
    this->xfn = xfn;
}

inline void Geo_FaceMetrics_Unstruct::SetFaceNormalY(RDouble *yfn)
{
    if (this->yfn && this->yfn != yfn)
    {
        delete [] this->yfn;
    }
    this->yfn = yfn;
}

inline void Geo_FaceMetrics_Unstruct::SetFaceNormalZ(RDouble *zfn)
{
    if (this->zfn && this->zfn != zfn)
    {
        delete [] this->zfn;
    }
    this->zfn = zfn;
}

inline void Geo_FaceMetrics_Unstruct::SetFaceArea(RDouble *area)
{
    if (this->area && this->area != area)
    {
        delete [] this->area;
    }
    this->area = area;
}

inline RDouble * Geo_FaceMetrics_Unstruct::GetFaceCenterX() const
{
    return this->xfc;
}

inline RDouble * Geo_FaceMetrics_Unstruct::GetFaceCenterY() const
{
    return this->yfc;
}

inline RDouble * Geo_FaceMetrics_Unstruct::GetFaceCenterZ() const
{
    return this->zfc;
}

inline RDouble * Geo_FaceMetrics_Unstruct::GetFaceNormalX() const
{
    return this->xfn;
}

inline RDouble * Geo_FaceMetrics_Unstruct::GetFaceNormalY() const
{
    return this->yfn;
}

inline RDouble * Geo_FaceMetrics_Unstruct::GetFaceNormalZ() const
{
    return this->zfn;
}

inline RDouble * Geo_FaceMetrics_Unstruct::GetFaceArea() const
{
    return this->area;
}