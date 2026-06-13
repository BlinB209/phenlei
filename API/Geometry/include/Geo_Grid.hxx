inline bool Grid::IsFinestGrid()
{
    if (this->GetFineGrid() == 0)
        return true;
    else
        return false;
}

inline bool Grid::IsCoarsestGrid()
{
    if (GetCoarseGrid() == 0)
        return true;
    else
        return false;
}

inline SimpleGrid * Grid::GetOldGrid()
{
    return oldGrid; 
}

inline void Grid::SetOldGrid(SimpleGrid *oldGrid)
{
    this->oldGrid = oldGrid;
}

inline void Grid::SetNTotalFace(int nTotalFace)
{
    this->nTotalFace = nTotalFace;
}

inline void Grid::SetNIFace(int nIFace)
{
    this->nIFace = nIFace;
}

inline void Grid::SetNTotalCell(int nTotalCell)
{
    this->nTotalCell = nTotalCell;
}

inline void Grid::SetNBoundFace(int nBoundFace)
{
    this->nBoundFace = nBoundFace;
}

inline int Grid::GetNTotalFace() const
{
    return nTotalFace;
}

inline int Grid::GetNTotalCell() const
{
    return nTotalCell;
}

inline int Grid::GetNTotalLine() const
{
    return nTotalLine;
}

inline int Grid::GetNBoundFace() const
{ 
    return nBoundFace;
}

inline int Grid::GetNIFace() const
{
    if (interfaceInfo)
    {
        return interfaceInfo->GetNIFace();
    }
    return 0;
}

inline void Grid::SetDim(int dim) 
{
    this->dimension = dim;
}

inline int Grid::GetDim() const
{
    return dimension;
}

inline int Grid::Type() const
{
    return type;
}

inline GridID * Grid::GetGridID() const
{
    return index;
}

inline void Grid::SetZoneID(int index)
{
    this->index->SetIndex(index);
}

inline void Grid::SetZoneLocalID(int index)
{
    this->index->SetLocalIndex(index);
}

inline int Grid::GetZoneID() const
{
    return index->GetIndex();
}

inline int Grid::GetZoneLocalID() const
{
    return index->GetLocalIndex();
}

inline SimpleVC * Grid::GetVolumeConditionIn()
{
    return this->volumeCondition;
}

inline void Grid::SetInterfaceInfo(InterfaceInfo *interfaceInfo)
{
    if (!interfaceInfo) return;

    delete this->interfaceInfo;
    this->interfaceInfo = interfaceInfo;

    if (interfaceFields) delete interfaceFields;
    interfaceFields = new InterfaceFields(interfaceInfo);
}

inline void Grid::SetInterpointInfo(InterpointInformation *interpointInformation)
{
    if (!interpointInformation)
    {
        return;
    }

    delete this->interpointInformation;
    this->interpointInformation = interpointInformation;
    
    delete interpointFields;
    interpointFields = new InterpointFields(interpointInformation);
}

inline InterfaceInfo * Grid::GetInterfaceInfo() const
{
    return interfaceInfo;
}

inline InterpointInformation *Grid::GetInterpointInfo() const
{
    return interpointInformation;
}

inline InterfaceFields *Grid::GetInterfaceFields() 
{
    return interfaceFields;
}

inline InterpointFields * Grid::GetInterpointFields()
{
    return interpointFields;
}

inline void Grid::SetOversetInfoProxy(OversetInfoProxy *oversetInfoProxy)
{
    this->oversetInfoProxy = oversetInfoProxy;
}

inline void Grid::SetOversetInformationProxy(OversetInformationProxy *oversetInformationProxyIn)
{
    this->oversetInformationProxy = oversetInformationProxyIn;
}

inline OversetInformationProxy * Grid::GetOversetInformationProxy() const
{
    return oversetInformationProxy;
}

inline OversetInfoProxy * Grid::GetOversetInfoProxy() const
{
    return oversetInfoProxy; 
}

inline int Grid::GetLevel() const
{
    return level;
}

inline void Grid::SetLevel(int level)
{ 
    this->level = level;
}

inline Grid * Grid::GetCoarseGrid() const
{
    return cGrid;
}

inline Grid * Grid::GetFineGrid() const
{
    return fGrid;
}

inline void Grid::SetCoarseGrid(Grid *cGrid)
{
    this->cGrid = cGrid;
}

inline void Grid::SetFineGrid(Grid *fGrid)
{
    this->fGrid = fGrid;
}

inline void Grid::CopyPara(Data_Param *gPara)
{
    this->gPara = gPara;
}

inline void Grid::CopyField(Data_Field *gField)
{
    this->gField = gField;
}

inline set < Data_IndexData > * Grid::GetDataSet()
{
    return gField->GetDataSet();
}

inline void Grid::UpdateDataPtr(const string &name, void *data)
{
    gField->UpdateDataPtr(name, data);
};

inline void * Grid::GetDataPtr(const string &name) const
{
    return gField->GetDataPtr(name);
};

inline void Grid::DeleteDataPtr(const string &name) const
{
    return gField->DeleteDataPtr(name);
};

inline void Grid::UpdateData(const string &name, void *data, int type, int size)
{
    gPara->UpdateData(name, data, type, size);
}

inline void Grid::GetData(const string &name, void *data, int type, int size)
{
    gPara->GetData(name, data, type, size);
}

inline int Grid::GetIBlock()
{
    return iBlock;
}

inline void Grid::SetIBlock(int iBlockIn)
{
    this->iBlock = iBlockIn;
}

inline int Grid::GetZoneProbesNumber() const
{
    return zoneProbesNumber;
}

inline void Grid::SetZoneProbesNumber(int zoneProbesNumberIn)
{
    this->zoneProbesNumber = zoneProbesNumberIn;
}

inline vector<int> Grid::GetZoneProbesGlobalID() const
{
    return zoneProbesGlobalID;
}

inline void Grid::SetZoneProbesGlobalID(vector<int> &zoneProbesGlobalIDIn)
{
    this->zoneProbesGlobalID = zoneProbesGlobalIDIn;
}

inline vector<int> Grid::GetZoneProbesLineID() const
{
    return zoneProbesLineID;
}

inline void Grid::SetZoneProbesLineID(vector<int> &zoneProbesLineIDIn)
{
    this->zoneProbesLineID = zoneProbesLineIDIn;
}

inline vector<int> Grid::GetZoneProbesSurfaceID() const
{
    return zoneProbesSurfaceID;
}

inline void Grid::SetZoneProbesSurfaceID(vector<int> &zoneProbesSurfaceIDIn)
{
    this->zoneProbesSurfaceID = zoneProbesSurfaceIDIn;
}

inline vector<vector<RDouble> > Grid::GetZoneProbesCoordinates() const
{
    return zoneProbesCoordinates;
}

inline void Grid::SetZoneProbesCoordinates(vector<vector<RDouble> > &zoneProbesCoordinatesIn)
{
    this->zoneProbesCoordinates = zoneProbesCoordinatesIn;
}

inline void Grid::SetVolumeCondition(SimpleVC *volumeConditionIn)
{
    if (volumeCondition)
    {
        delete volumeCondition;
    }
    this->volumeCondition = volumeConditionIn;
}

inline void Grid::SetFileIndex(int fileIndex)
{
    this->fileIndexCurrentGridBelong = fileIndex;
}

inline int Grid::GetFileIndex()
{
    return this->fileIndexCurrentGridBelong;
}
