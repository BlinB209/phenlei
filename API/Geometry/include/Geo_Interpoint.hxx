inline int * InterpointPatch::GetPointIndexForReceive() const
{
    return pointIndexForReceive; 
}

inline void InterpointPatch::SetPointIndexForReceive(int *pointIndexForReceive)
{
    this->pointIndexForReceive = pointIndexForReceive; 
}

inline int * InterpointPatch::GetPointIndexForSend() const
{
    return pointIndexForSend;
}

inline void InterpointPatch::SetPointIndexForSend(int *pointIndexForSend)
{
    this->pointIndexForSend = pointIndexForSend; 
}

inline int InterpointPatch::GetNumberOfPoint() const
{
    return pointNumberOfPatch;
}

inline void InterpointPatch::SetNumberOfPoint(int pointNumberOfPatch) 
{ 
    this->pointNumberOfPatch = pointNumberOfPatch;
}

inline int InterpointPatch::GetZoneIndexOfNeighbor() const 
{ 
    return zoneIndexOfNeighbor; 
}

inline void InterpointPatch::SetZoneIndexOfNeighbor(int zoneIndexOfNeighbor) 
{
    this->zoneIndexOfNeighbor = zoneIndexOfNeighbor; 
}

inline int InterpointPatch::GetZoneID() const 
{ 
    return zoneID; 
}

inline void InterpointInformation::SetNumberOfInterpoints(int numberOfInterpoints)
{ 
    this->numberOfInterpoints = numberOfInterpoints; 
}

inline int InterpointInformation::GetNumberOfInterpoints() const
{ 
    return numberOfInterpoints;
}

inline void InterpointInformation::SetNumberOfNeighbor(int numberOfNeighbor) 
{ 
    this->numberOfNeighbor = numberOfNeighbor;
}

inline int InterpointInformation::GetNumberOfNeighbor() const 
{
    return numberOfNeighbor;
}

inline InterpointPatch ** InterpointInformation::GetInterpointPatch() const
{
    return interPointPatch;
}

inline InterpointPatch * InterpointInformation::GetInterpointPatch(int iNeighbor) const
{
    return interPointPatch[iNeighbor];
}

inline int * InterpointInformation::GetInterPoint2ZoneID() const 
{ 
    return interPoint2ZoneID; 
}

inline int * InterpointInformation::GetInterPoint2InterPointID() const 
{ 
    return interPoint2InterPointID;
}

inline int * InterpointInformation::GetInterPoint2GlobalPoint() const 
{
    return interPoint2GlobalPoint;
}

inline int * InterpointInformation::GetCellNumberOfInterPoint() const 
{
    return cellNumberOfInterPoint;
}

inline int * InterpointInformation::GetTotalZonesOfInterPoint() const 
{
    return totalZonesOfInterPoint;
}
inline int * InterpointInformation::GetLabelOfInterPoint() const
{
    return labelOfInterPoint;
}

inline void InterpointInformation::SetIsNeighborInterpointFound(bool *data)
{
    this->isNeighborInterpointFound = data;
}

inline bool * InterpointInformation::GetIsNeighborInterpointFound()
{
    return isNeighborInterpointFound;
}

inline int InterpointInformation::GetZoneIndexOfNeighbor(int iNeighbor) const 
{
    InterpointPatch *interPointPatch = this->GetInterpointPatch(iNeighbor);
    return interPointPatch->GetZoneIndexOfNeighbor();
}

inline int * InterpointInformation::GetPointIndexForReceive(int iNeighbor) const
{
    InterpointPatch *interPointPatch = this->GetInterpointPatch(iNeighbor);
    return interPointPatch->GetPointIndexForReceive();
}

inline int * InterpointInformation::GetPointIndexForSend(int iNeighbor) const
{ 
    InterpointPatch *interPointPatch = this->GetInterpointPatch(iNeighbor);
    return interPointPatch->GetPointIndexForSend();
}

inline Data_ParamFieldSuite * InterpointInformation::GetSendDataStorage(int iStore) const
{ 
    return &dataForSend[iStore]; 
}

inline Data_ParamFieldSuite * InterpointInformation::GetReceiveDataStorage(int iStore) const
{ 
    return &dataForReceive[iStore];
}

inline int InterpointInformation::GetNumberOfInterpointsForNeighbor(int iNeighbor) const
{ 
    InterpointPatch *interPointPatch = this->GetInterpointPatch(iNeighbor);
    return interPointPatch->GetNumberOfPoint();
}

inline uint_t InterpointFields::Size()
{
    return variableNames.size();
}

inline int InterpointFields::GetDimension(const int iData)
{
    return variableDimensions[iData];
}

inline string & InterpointFields::GetName(const int iData)
{
    return variableNames[iData];
}

inline int InterpointFields::GetSolverIndex(const int iData)
{
    return solverIndex[iData];
}