inline RDouble * Geo_CellMetrics_Unstruct::GetCellCenterX() const
{
    return xcc;
}

inline RDouble * Geo_CellMetrics_Unstruct::GetCellCenterY() const
{
    return ycc;
}

inline RDouble * Geo_CellMetrics_Unstruct::GetCellCenterZ() const
{
    return zcc;
}

inline RDouble * Geo_CellMetrics_Unstruct::GetCellVolume()  const
{
    return vol;
}

inline RDouble * Geo_CellMetrics_Unstruct::GetCellSkewness()  const
{
    return cellSkewness;
}

inline void Geo_CellMetrics_Unstruct::SetCellCenterX(RDouble *xcc)
{
    if (this->xcc && this->xcc != xcc)
    {
        delete [] this->xcc;
    }
    this->xcc = xcc;
}

inline void Geo_CellMetrics_Unstruct::SetCellCenterY(RDouble *ycc)
{
    if (this->ycc && this->ycc != ycc)
    {
        delete [] this->ycc;
    }
    this->ycc = ycc;
}

inline void Geo_CellMetrics_Unstruct::SetCellCenterZ(RDouble *zcc)
{
    if (this->zcc && this->zcc != zcc)
    {
        delete [] this->zcc;
    }
    this->zcc = zcc;
}

inline void Geo_CellMetrics_Unstruct::SetCellVolume(RDouble *vol)
{
    if (this->vol && this->vol != vol)
    {
        delete [] this->vol;
    }
    this->vol = vol;
}

inline void Geo_CellMetrics_Unstruct::SetCellSkewness(RDouble *skewness)
{
    if (this->cellSkewness && this->cellSkewness != skewness)
    {
        delete [] this->cellSkewness;
    }
    this->cellSkewness = skewness;
}

inline void Geo_CellMetrics_Unstruct::SetLargestLocalGridLength(RDouble *length)
{
    if (this->largestLocalGridLength && this->largestLocalGridLength != length)
    {
        delete [] this->largestLocalGridLength;
    }
    this->largestLocalGridLength = length;
}

inline void Geo_CellMetrics_Unstruct::SetSubgridLength(RDouble *length)
{
    if (this->subgridLength && this->subgridLength != length)
    {
        delete [] this->subgridLength;
    }
    this->subgridLength = length;
}

inline RDouble * Geo_CellMetrics_Unstruct::GetLargestLocalGridLength()
{
    return largestLocalGridLength;
}

inline RDouble * Geo_CellMetrics_Unstruct::GetSubgridLength()
{
    return subgridLength;
}