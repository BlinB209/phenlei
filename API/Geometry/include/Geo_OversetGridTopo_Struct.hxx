inline int Geo_OversetGridTopo_Struct::GetNumberOfCores() const
{
    return numberOfCores;
}

inline int Geo_OversetGridTopo_Struct::GetZoneStartPointLabel() const
{
    return zoneStartPointLabel;
}

inline int * Geo_OversetGridTopo_Struct::GetILinkPointLabel() const
{
    return iLinkPointLabel;
}

inline int * Geo_OversetGridTopo_Struct::GetJLinkPointLabel() const
{
    return jLinkPointLabel;
}

inline int * Geo_OversetGridTopo_Struct::GetKLinkPointLabel() const
{
    return kLinkPointLabel;
}

inline int * Geo_OversetGridTopo_Struct::GetCellTypeContainer() const
{
    return iBlank;
}

inline void Geo_OversetGridTopo_Struct::SetNumberOfCores(int numberOfCores)
{
    this->numberOfCores = numberOfCores;
}

inline void Geo_OversetGridTopo_Struct::SetZoneStartCenterLabel(int zoneStartCenterLabel)
{
    this->zoneStartCenterLabel = zoneStartCenterLabel;
}

inline void Geo_OversetGridTopo_Struct::SetZoneStartPointLabel(int zoneStartPointLabel)
{
    this->zoneStartPointLabel = zoneStartPointLabel;
}

inline int * Geo_OversetGridTopo_Struct::GetHingedPointContainer() const
{
    return hingedPointContainer;
}

inline RDouble * Geo_OversetGridTopo_Struct::GetXCoreContainer() const
{
    return xCoreContainer;
}

inline RDouble * Geo_OversetGridTopo_Struct::GetYCoreContainer() const
{
    return yCoreContainer;
}

inline RDouble * Geo_OversetGridTopo_Struct::GetZCoreContainer() const
{
    return zCoreContainer;
}
