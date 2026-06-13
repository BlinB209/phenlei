inline int * InterfacePatch::GetFaceIndexForRecv() const
{
    return faceIndexForRecv;
}

inline void InterfacePatch::SetFaceIndexForRecv(int *faceIndexForRecv)
{
    this->faceIndexForRecv = faceIndexForRecv;
}

inline int * InterfacePatch::GetFaceIndexForSend() const
{
    return faceIndexForSend;
}

inline void InterfacePatch::SetFaceIndexForSend(int *faceIndexForSend)
{
    this->faceIndexForSend = faceIndexForSend;
}

inline int InterfacePatch::GetNumberOfFace() const
{
    return faceNumberOfPatch;
}

inline void InterfacePatch::SetNumberOfFace(int faceNumberOfPatch)
{
    this->faceNumberOfPatch = faceNumberOfPatch;
}

inline int InterfacePatch::GetZoneIndexOfNeighbor() const
{
    return zoneIndexOfNeighbor;
}

inline void InterfacePatch::SetZoneIndexOfNeighbor(int zoneIndexOfNeighbor)
{
    this->zoneIndexOfNeighbor = zoneIndexOfNeighbor;
}

inline int InterfacePatch::GetZoneID() const
{
    return zoneid;
}

//! GMRESParallel
inline void InterfacePatch::SetCellIndexOfLocal(int *cellIndexOfLocal)
{
    this->cellIndexOfLocal = cellIndexOfLocal;
}

inline int* InterfacePatch::GetCellIndexOfLocal() const
{
    return this->cellIndexOfLocal;
}

inline void InterfacePatch::SetCellIndexOfLocalGhost(int *cellIndexLocalGhost)
{
    this->cellIndexOfLocalGhost = cellIndexLocalGhost;
}

inline int* InterfacePatch::GetCellIndexOfLocalGhost() const
{
    return this->cellIndexOfLocalGhost;
}

inline void InterfacePatch::SetCellIndexOfNeighbor(int *cellIndexOfNeighbor)
{
    this->cellIndexOfNeighbor = cellIndexOfNeighbor;
}

inline int* InterfacePatch::GetCellIndexOfNeighbor() const
{
    return this->cellIndexOfNeighbor;
}

inline void InterfaceInfo::SetNIFace(int nIFace)
{
    this->nIFace = nIFace;
}

inline int InterfaceInfo::GetNIFace() const
{
    return nIFace;
}

inline void InterfaceInfo::SetNumberOfNeighbor(int numberOfNeighbor)
{
    this->numberOfNeighbor = numberOfNeighbor;
}

inline int InterfaceInfo::GetNumberOfNeighbor() const
{
    return numberOfNeighbor;
}

inline vector <InterfacePatch *> InterfaceInfo::GetInterfacePatch() const
{
    return interFacePatch;
}

inline InterfacePatch * InterfaceInfo::GetInterfacePatch(int ineighbor) const
{
    return interFacePatch[ineighbor];
}

inline int * InterfaceInfo::GetInterFace2ZoneID() const
{
    return interFace2ZoneID;
}

inline int * InterfaceInfo::GetInterFace2CellID() const 
{ 
    return interFace2CellID; 
}

inline int * InterfaceInfo::GetInterFace2InterFaceID() const
{
    return interFace2InterFaceID;
}

inline int * InterfaceInfo::GetInterFace2BoundaryFace() const
{
    return interFace2BoundaryFace;
}

inline int * InterfaceInfo::GetInterFaceDirection() const
{
    return interFaceDirection;
}

inline void InterfaceInfo::SetIsNeighborInterfaceFound(bool *data)
{
    this->isNeighborInterfaceFound = data;
}

inline bool * InterfaceInfo::GetIsNeighborInterfaceFound()
{
    return isNeighborInterfaceFound;
}

inline int InterfaceInfo::GetZoneIndexOfNeighbor(int ineighbor) const
{
    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    return interFacePatch->GetZoneIndexOfNeighbor();
}

inline int * InterfaceInfo::GetFaceIndexForRecv(int ineighbor) const
{
    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    return interFacePatch->GetFaceIndexForRecv();
}

inline int * InterfaceInfo::GetFaceIndexForSend(int ineighbor) const
{
    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    return interFacePatch->GetFaceIndexForSend();
}

//! GMRESParallel
inline int * InterfaceInfo::GetCellIndexOfNeighbor(int ineighbor) const
{
    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    return interFacePatch->GetCellIndexOfNeighbor();
}

inline int * InterfaceInfo::GetGlobalNeighborCellIndex() const
{
    return this->globalNeighborCellIndex;
}

inline void InterfaceInfo::SetGlobalNeighborCellIndex(int* globalNeighborCellIndex)
{
    this->globalNeighborCellIndex = globalNeighborCellIndex;
}

inline int * InterfaceInfo::GetLocalInterfaceCellIndex() const
{
    return this->localInterfaceCellIndex;
}

inline void InterfaceInfo::SetLocalInterfaceCellIndex(int* localInterfaceCellIndex)
{
    this->localInterfaceCellIndex = localInterfaceCellIndex;
}

inline int * InterfaceInfo::GetLocalInterfacePhysicalCellIndex() const
{
    return this->localInterfacePhysicalCellIndex;
}

inline void InterfaceInfo::SetLocalInterfacePhysicalCellIndex(int *localInterfacePhysicalCellIndex)
{
    this->localInterfacePhysicalCellIndex = localInterfacePhysicalCellIndex;
}

/**
 * @brief return the ghost cell index with respect to the global mesh for interface
 * 
 * @param cell local ghost cell index
 * @return int 
 */
inline int InterfaceInfo::MatchedGlobalNeighborCellIndex(int cell)
{
    for (int i = 0; i < nIFace; ++i)
    {
        if(cell == localInterfaceCellIndex[i])
        {
            return globalNeighborCellIndex[i];
        }
    }
    return -1;
}

/**
 * @brief return the physical cell index with respect to the local mesh for interface
 * 
 * @param cell local ghost cell index
 * @return int 
 */
inline int InterfaceInfo::MatchedLocalPhysicalCellIndex(int cell)
{
    for (int i = 0; i < nIFace; ++i)
    {
        if(cell == localInterfaceCellIndex[i])
        {
            return localInterfacePhysicalCellIndex[i];
        }
    }
    return -1;
}

inline Data_ParamFieldSuite * InterfaceInfo::GetSendDataStorage(int istore) const
{
    return &dsend[istore];
}

inline Data_ParamFieldSuite * InterfaceInfo::GetRecvDataStorage(int istore) const
{
    return &drecv[istore];
}

inline int InterfaceInfo::GetNIFaceOfNeighbor(int ineighbor) const
{
    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    return interFacePatch->GetNumberOfFace();
}

inline uint_t InterfaceFields::Size()
{
    return variableNames.size();
}

inline int InterfaceFields::GetDim(const int iData)
{
    return variableDimensions[iData];
}

inline string & InterfaceFields::GetName(const int iData)
{
    return variableNames[iData];
}

inline int InterfaceFields::GetSolverIndex(const int iData)
{
    return solverIndex[iData];
}
