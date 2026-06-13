inline RDouble4D * Geo_FaceMetrics_Struct::GetFaceNormalX() const
{
    return xfn;
}

inline RDouble4D * Geo_FaceMetrics_Struct::GetFaceNormalY() const
{
    return yfn;
}

inline RDouble4D * Geo_FaceMetrics_Struct::GetFaceNormalZ() const
{
    return zfn;
}

inline RDouble4D * Geo_FaceMetrics_Struct::GetFaceArea() const
{
    return area;
}

inline RDouble4D * Geo_FaceMetrics_Struct::GetFaceVectorX() const
{
    return xFaceVector;
}

inline RDouble4D * Geo_FaceMetrics_Struct::GetFaceVectorY() const
{
    return yFaceVector;
}

inline RDouble4D * Geo_FaceMetrics_Struct::GetFaceVectorZ() const
{
    return zFaceVector;
}

inline RDouble5D * Geo_FaceMetrics_Struct::GetFaceVector_FD() const
{
    return xyzFaceVector_FD;
}

inline RDouble5D * Geo_FaceMetrics_Struct::GetFaceNormal_FD() const
{ 
    return xyzFaceNormal_FD;
}

inline void Geo_FaceMetrics_Struct::DeallocateAll()
{
    if (xfn)  { delete xfn;      xfn = 0; }
    if (yfn)  { delete yfn;      yfn = 0; }
    if (zfn)  { delete zfn;      zfn = 0; }
    if (area) { delete area;    area = 0; }

    if (xFaceVector) { delete xFaceVector;    xFaceVector = 0; }
    if (yFaceVector) { delete yFaceVector;    yFaceVector = 0; }
    if (zFaceVector) { delete zFaceVector;    zFaceVector = 0; }

    if (xyzFaceVector_FD) { delete xyzFaceVector_FD;    xyzFaceVector_FD = 0; }
    if (xyzFaceNormal_FD) { delete xyzFaceNormal_FD;    xyzFaceNormal_FD = 0; }
}