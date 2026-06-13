inline int Geo_MultiGridInfo_Struct::GetMultigridStepI() const
{
    return istp;
}

inline int Geo_MultiGridInfo_Struct::GetMultigridStepJ() const
{
    return jstp;
}

inline int Geo_MultiGridInfo_Struct::GetMultigridStepK() const
{
    return kstp;
}

inline void Geo_MultiGridInfo_Struct::SetMultigridStepI(int istp)
{
    this->istp = istp;
}

inline void Geo_MultiGridInfo_Struct::SetMultigridStepJ(int jstp)
{
    this->jstp = jstp;
}

inline void Geo_MultiGridInfo_Struct::SetMultigridStepK(int kstp)
{
    this->kstp = kstp;
}