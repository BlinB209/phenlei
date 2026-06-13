inline void HoleGrid::SetZoneStart(int iZone)
{
    this->zone_st = iZone;
}

inline void HoleGrid::SetZoneEnd(int iZone)
{
    this->zone_ed = iZone;
}

inline void HoleGrid::SetZoneAttached(int iZone)
{
    this->zone_at = iZone;
}

inline void HoleGrid::SetNumberOfSurface(int nSurf)
{
    this->ns = nSurf;
}

inline void HoleGrid::SetNumberOfPoints(int np)
{
    this->np = np;
}

inline int HoleGrid::GetZoneStart()
{
    return this->zone_st;
}

inline int HoleGrid::GetZoneEnd()
{
    return this->zone_ed;
}

inline int HoleGrid::GetZoneAttached()
{
    return this->zone_at;
}

inline int HoleGrid::GetNumberOfSurface()
{
    return this->ns;
}

inline int HoleGrid::GetNumberOfPoints()
{
    return this->np;
}

inline int * HoleGrid::GetNSI()
{
    return & (*nsi)[0];
}

inline int * HoleGrid::GetNSJ() 
{
    return & (*nsj)[0];
}

inline int * HoleGrid::GetNSP()
{
    return & (*nsp)[0];
}

inline RDouble * HoleGrid::GetX()
{
    return & (*x)[0];
}

inline RDouble * HoleGrid::GetY()
{
    return & (*y)[0];
}

inline RDouble * HoleGrid::GetZ()
{
    return & (*z)[0];
}

inline RDouble * HoleGrid::GetBox()
{
    return & (*box)[0];
}

inline void BackgroundTree::SetZoneViewRange(int st, int ed)
{
    this->zone_st = st;
    this->zone_ed = ed;
}

inline void BackgroundTree::SetUnitCellLocation(int iUnitCell, int iZone, int iCore, int jCore, int kCore)
{
    zone_location[iUnitCell] = iZone;
    i_location[iUnitCell]    = iCore;
    j_location[iUnitCell]    = jCore;
    k_location[iUnitCell]    = kCore;
}

inline bool BackgroundTree::IsEmpty()
{
    return isEmptyBinaryTree;
}

inline void OversetCellCollector::SetZoneViewRange(int st, int ed)
{
    this->zone_st = st;
    this->zone_ed = ed;
}

inline int OversetStructGrid::GetZoneInverseIndex(int i_view)
{
    return zoneInverseIndex[i_view];
}