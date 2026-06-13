inline void Geo_LSQWeight_Unstruct::SetLeastSquareIWT(RDouble *iwt) 
{
    this->iwt = iwt;
}

inline void Geo_LSQWeight_Unstruct::SetLeastSquareIXX(RDouble *ixx)
{
    this->ixx = ixx;
}

inline void Geo_LSQWeight_Unstruct::SetLeastSquareIYY(RDouble *iyy) 
{
    this->iyy = iyy;
}

inline void Geo_LSQWeight_Unstruct::SetLeastSquareIZZ(RDouble *izz)
{
    this->izz = izz;
}

inline void Geo_LSQWeight_Unstruct::SetLeastSquareIXY(RDouble *ixy) 
{
    this->ixy = ixy;
}

inline void Geo_LSQWeight_Unstruct::SetLeastSquareIXZ(RDouble *ixz)
{
    this->ixz = ixz;
}

inline void Geo_LSQWeight_Unstruct::SetLeastSquareIYZ(RDouble *iyz) 
{
    this->iyz = iyz;
}

inline void Geo_LSQWeight_Unstruct::SetFaceMark(char *fMark)
{
    this->fMark = fMark;
}

inline RDouble * Geo_LSQWeight_Unstruct::GetLeastSquareIWT() const
{
    return iwt;
}

inline RDouble * Geo_LSQWeight_Unstruct::GetLeastSquareIXX() const
{
    return ixx;
}

inline RDouble * Geo_LSQWeight_Unstruct::GetLeastSquareIYY() const
{
    return iyy;
}

inline RDouble * Geo_LSQWeight_Unstruct::GetLeastSquareIZZ() const
{
    return izz;
}

inline RDouble * Geo_LSQWeight_Unstruct::GetLeastSquareIXY() const
{
    return ixy;
}

inline RDouble * Geo_LSQWeight_Unstruct::GetLeastSquareIXZ() const
{
    return ixz;
}

inline RDouble * Geo_LSQWeight_Unstruct::GetLeastSquareIYZ() const
{
    return iyz;
}

inline char * Geo_LSQWeight_Unstruct::GetFaceMark() const
{
    return fMark;
}
