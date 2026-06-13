inline void Geo_CellMetrics_Struct::DeallocateAll()
{
    if (xcc) { delete xcc;    xcc = 0; }
    if (ycc) { delete ycc;    ycc = 0; }
    if (zcc) { delete zcc;    zcc = 0; }
    if (vol) { delete vol;    vol = 0; }

    if (xlen) { delete xlen;    xlen = 0; }
    if (ylen) { delete ylen;    ylen = 0; }
    if (zlen) { delete zlen;    zlen = 0; }

    if (jacobian) { delete jacobian;    jacobian = 0; }
}

inline RDouble3D * Geo_CellMetrics_Struct::GetCellVolume() const
{
    return vol;
}

inline RDouble3D * Geo_CellMetrics_Struct::GetCellLengthX() const
{
    return xlen;
}

inline RDouble3D * Geo_CellMetrics_Struct::GetCellLengthY() const
{
    return ylen;
}

inline RDouble3D * Geo_CellMetrics_Struct::GetCellLengthZ() const
{
    return zlen;
}

inline RDouble3D * Geo_CellMetrics_Struct::GetCellJacobian() const
{
    return jacobian;
}

inline RDouble3D * Geo_CellMetrics_Struct::GetCellCenterX() const
{
    return xcc;
}

inline RDouble3D * Geo_CellMetrics_Struct::GetCellCenterY() const
{
    return ycc;
}

inline RDouble3D * Geo_CellMetrics_Struct::GetCellCenterZ() const
{
    return zcc;
}

inline RDouble3D * Geo_CellMetrics_Struct::GetLargestLocalGridLength()
{
    return largestLocalGridLength;
}

inline RDouble3D * Geo_CellMetrics_Struct::GetSmallestLocalGridLength()
{
    return smallestLocalGridLength;
}

inline RDouble3D * Geo_CellMetrics_Struct::GetSubgridLength()
{
    return subgridLength;
}

inline void Geo_CellMetrics_Struct::SetLargestLocalGridLength(RDouble3D *length)
{
    if (this->largestLocalGridLength && this->largestLocalGridLength != length)
    {
        delete this->largestLocalGridLength;
    }
    this->largestLocalGridLength = length;
}

inline void Geo_CellMetrics_Struct::SetSmallestLocalGridLength(RDouble3D *length)
{
    if (this->smallestLocalGridLength && this->smallestLocalGridLength != length)
    {
        delete this->smallestLocalGridLength;
    }
    this->smallestLocalGridLength = length;
}

inline void Geo_CellMetrics_Struct::SetSubgridLength(RDouble3D *length)
{
    if (this->subgridLength && this->subgridLength != length)
    {
        delete this->subgridLength;
    }
    this->subgridLength = length;
}