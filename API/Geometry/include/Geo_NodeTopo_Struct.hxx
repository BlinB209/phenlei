inline void Geo_NodeTopo_Struct::SetNI(int ni)
{
    this->ni = ni;
}

inline void Geo_NodeTopo_Struct::SetNJ(int nj)
{
    this->nj = nj;
}

inline void Geo_NodeTopo_Struct::SetNK(int nk)
{
    this->nk = nk;
}

inline int Geo_NodeTopo_Struct::GetNI() const
{
    return ni;
}

inline int Geo_NodeTopo_Struct::GetNJ() const
{
    return nj;
}

inline int Geo_NodeTopo_Struct::GetNK() const
{
    return nk;
}

inline RDouble3D * Geo_NodeTopo_Struct::GetStructX() const
{
    return structx;
}

inline RDouble3D * Geo_NodeTopo_Struct::GetStructY() const
{
    return structy;
}

inline RDouble3D * Geo_NodeTopo_Struct::GetStructZ() const
{
    return structz;
}

inline void Geo_NodeTopo_Struct::DeallocateAll()
{
    if (structx) { delete structx;    structx = 0; }
    if (structy) { delete structy;    structy = 0; }
    if (structz) { delete structz;    structz = 0; }
}