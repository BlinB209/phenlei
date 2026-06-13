inline void UnstructGrid::SetLeftCellOfFace(int *left_cell_of_face_in)
{
    this->faceTopology->SetLeftCellOfFace(left_cell_of_face_in);
}

inline void UnstructGrid::SetRightCellOfFace(int *right_cell_of_face_in)
{
    this->faceTopology->SetRightCellOfFace(right_cell_of_face_in);
}

inline void UnstructGrid::SetNodeNumberOfEachFace(int *node_number_of_each_face)
{
    this->faceTopology->SetNodeNumberOfEachFace(node_number_of_each_face);
}

inline void UnstructGrid::SetFace2Node(int *face2node)
{
    this->faceTopology->SetFace2Node(face2node);
}

inline void UnstructGrid::SetFace2NodeSubscript( long long int *face2nodeSubscript)
{
    this->faceTopology->SetFace2NodeSubscript(face2nodeSubscript);
}

inline void UnstructGrid::SetCell2CoarseGridCell(int *cell2coarsegridcell)
{
    this->cell2coarsegridcell = cell2coarsegridcell;
}

inline void UnstructGrid::SetLI_nLine(int nLine)
{
    this->lineTopology->SetLI_nLine(nLine);
}

inline void UnstructGrid::SetLI_nCellsInLine(int *nCellsInLines)
{
    this->lineTopology->SetLI_nCellsInLine(nCellsInLines);
}

inline void UnstructGrid::SetLI_CellOfLine(int **CellOfLine)
{
    this->lineTopology->SetLI_CellOfLine(CellOfLine);
}

inline void UnstructGrid::SetLI_FaceOfLine(int **FaceOfLine)
{
    this->lineTopology->SetLI_FaceOfLine(FaceOfLine);
}

inline void UnstructGrid::SetLI_LineOfCell(int *LineOfCell)
{
    this->lineTopology->SetLI_LineOfCell(LineOfCell);
}

inline void UnstructGrid::Setfc2cL(int *fc2cL)
{
    this->fc2cL = fc2cL;
}

inline void UnstructGrid::Setfc2cR(int *fc2cR)
{
    this->fc2cR = fc2cR;
}

inline void UnstructGrid::SetBCRecord(UnstructBCSet **bcr)
{
    this->bcr = bcr;
}

inline void UnstructGrid::SetUnstructBCSet(UnstructBCSet *unstructBCSet)
{
    this->unstructBCSet = unstructBCSet;
}

inline void UnstructGrid::SetBlankIndex(int *iBlank)
{
    this->iBlank = iBlank;
}

inline void UnstructGrid::SetZoneProbesCellID(vector <int> zoneProbesCellID)
{
    this->zoneProbesCellID = zoneProbesCellID;
}

inline int * UnstructGrid::GetLeftCellOfFace() const
{
    return this->faceTopology->GetLeftCellOfFace();
}

inline int * UnstructGrid::GetRightCellOfFace() const
{
    return this->faceTopology->GetRightCellOfFace();
}

inline int *UnstructGrid::GetNodeNumberOfEachFace() const
{
    return this->faceTopology->GetNodeNumberOfEachFace();
}

inline int *UnstructGrid::GetFace2Node() const
{
    return this->faceTopology->GetFace2Node();
}

inline RDouble *UnstructGrid::GetFaceCenterX() const
{
    return this->faceMetrics->GetFaceCenterX();
}

inline RDouble *UnstructGrid::GetFaceCenterY() const
{
    return this->faceMetrics->GetFaceCenterY();
}

inline RDouble *UnstructGrid::GetFaceCenterZ() const
{
    return this->faceMetrics->GetFaceCenterZ();
}

inline RDouble *UnstructGrid::GetFaceNormalX() const
{
    return this->faceMetrics->GetFaceNormalX();
}

inline RDouble *UnstructGrid::GetFaceNormalY() const
{
    return this->faceMetrics->GetFaceNormalY();
}

inline RDouble *UnstructGrid::GetFaceNormalZ() const
{
    return this->faceMetrics->GetFaceNormalZ();
}

inline RDouble *UnstructGrid::GetFaceArea() const
{
    return this->faceMetrics->GetFaceArea();
}

inline RDouble *UnstructGrid::GetCellCenterX() const
{
    return cellMetrics->GetCellCenterX();
}

inline RDouble *UnstructGrid::GetCellCenterY() const
{
    return cellMetrics->GetCellCenterY();
}

inline RDouble * UnstructGrid::GetCellCenterZ() const
{
    return cellMetrics->GetCellCenterZ();
}

inline RDouble * UnstructGrid::GetCellVolume() const
{
    return cellMetrics->GetCellVolume();
}

inline RDouble * UnstructGrid::GetCellSkewness() const
{
    return cellMetrics->GetCellSkewness();
}

inline RDouble * UnstructGrid::GetFaceVelocityX() const
{
    return dynamicGridMetrics->GetFaceVelocityX();
}

inline RDouble *UnstructGrid::GetFaceVelocityY() const
{
    return dynamicGridMetrics->GetFaceVelocityY();
}

inline RDouble * UnstructGrid::GetFaceVelocityZ() const
{
    return dynamicGridMetrics->GetFaceVelocityZ();
}

inline RDouble * UnstructGrid::GetCellVelocityX() const
{
    return dynamicGridMetrics->GetCellVelocityX();
}

inline RDouble * UnstructGrid::GetCellVelocityY() const
{
    return dynamicGridMetrics->GetCellVelocityY();
}

inline RDouble * UnstructGrid::GetCellVelocityZ() const
{
    return dynamicGridMetrics->GetCellVelocityZ();
}

inline RDouble * UnstructGrid::GetFaceNormalVelocity() const
{
    return dynamicGridMetrics->GetFaceNormalVelocity();
}

inline RDouble * UnstructGrid::GetCellVolumeOld() const
{
    return dynamicGridMetrics->GetCellVolumeOld();
}

inline RDouble *UnstructGrid::GetLeastSquareIWT() const
{
    return leastSquareWeights->GetLeastSquareIWT();
}

inline RDouble *UnstructGrid::GetLeastSquareIXX() const
{
    return leastSquareWeights->GetLeastSquareIXX();
}

inline RDouble *UnstructGrid::GetLeastSquareIYY() const
{
    return leastSquareWeights->GetLeastSquareIYY();
}

inline RDouble *UnstructGrid::GetLeastSquareIZZ() const
{
    return leastSquareWeights->GetLeastSquareIZZ();
}

inline RDouble *UnstructGrid::GetLeastSquareIXY() const
{
    return leastSquareWeights->GetLeastSquareIXY();
}

inline RDouble *UnstructGrid::GetLeastSquareIXZ() const
{
    return leastSquareWeights->GetLeastSquareIXZ();
}

inline RDouble *UnstructGrid::GetLeastSquareIYZ() const
{
    return leastSquareWeights->GetLeastSquareIYZ();
}

inline void UnstructGrid::SetLeastSquareIWT(RDouble *iwt)
{
    leastSquareWeights->SetLeastSquareIWT(iwt);
}

inline void UnstructGrid::SetLeastSquareIXX(RDouble *ixx)
{
    leastSquareWeights->SetLeastSquareIXX(ixx);
}

inline void UnstructGrid::SetLeastSquareIYY(RDouble *iyy)
{
    leastSquareWeights->SetLeastSquareIYY(iyy);
}

inline void UnstructGrid::SetLeastSquareIZZ(RDouble *izz)
{
    leastSquareWeights->SetLeastSquareIZZ(izz);
}

inline void UnstructGrid::SetLeastSquareIXY(RDouble *ixy)
{
    leastSquareWeights->SetLeastSquareIXY(ixy);
}

inline void UnstructGrid::SetLeastSquareIXZ(RDouble *ixz)
{
    leastSquareWeights->SetLeastSquareIXZ(ixz);
}

inline void UnstructGrid::SetLeastSquareIYZ(RDouble *iyz)
{
    leastSquareWeights->SetLeastSquareIYZ(iyz);
}

inline void UnstructGrid::SetFaceMark(char *fMark)
{
    leastSquareWeights->SetFaceMark(fMark);
}

inline char *UnstructGrid::GetFaceMark() const
{
    return leastSquareWeights->GetFaceMark();
}

inline int * UnstructGrid::GetCell2CoarseGridCell(void) const
{
    return cell2coarsegridcell;
}

inline int * UnstructGrid::Getfc2cL(void) const
{
    return fc2cL;
}

inline int *UnstructGrid::Getfc2cR(void) const
{
    return fc2cR;
}

inline UnstructBCSet ** UnstructGrid::GetBCRecord(void) const
{
    return bcr;
}

inline void UnstructGrid::CreateUnstructBCSet(int nBCRegion)
{
    if (!unstructBCSet)
    {
        unstructBCSet = new UnstructBCSet(this->GetZoneID());
    }

    unstructBCSet->CreatenBCRegion(nBCRegion);
}

inline RDouble * UnstructGrid::GetNormalDistanceOfC2C(void) const
{
    return normalDistanceC2C;
}

inline int * UnstructGrid::GetBlankIndex(void) const
{
    return iBlank;
}

inline vector <int> UnstructGrid::GetZoneProbesCellID() const
{
    return zoneProbesCellID;
}

inline RDouble *UnstructGrid::GetWallDist() const
{
    return walldist;
}

inline RDouble *UnstructGrid::GetNearestWallFaceNormalX() const
{
    return nearestwallfacenormalx;
}

inline RDouble *UnstructGrid::GetNearestWallFaceNormalY() const
{
    return nearestwallfacenormaly;
}

inline RDouble *UnstructGrid::GetNearestWallFaceNormalZ() const
{
    return nearestwallfacenormalz;
}

inline RDouble * UnstructGrid::GetWallDistNode() const
{
    return walldistNode;
}

inline RDouble * UnstructGrid::GetLamdax() const
{
    return lamdax;
}

inline RDouble *UnstructGrid::GetLamday() const
{
    return lamday;
}

inline RDouble *UnstructGrid::GetLamdaz() const
{
    return lamdaz;
}

inline int * UnstructGrid::Getknode() const
{
    return knode;
}

//! Set the cell to node: node index of each cell.
inline void UnstructGrid::SetCell2Node(int *cell2node)
{
    cellTopology->SetCell2Node(cell2node);
}

inline map<int, int> &UnstructGrid::GetBoundaryPointLabel()
{
    return boundaryPointLabel;
}

//! Set the number of nodes in each cell.
inline void UnstructGrid::SetNodeNumberOfEachCell(int *node_number_of_each_cell)
{
    cellTopology->SetNodeNumberOfEachCell(node_number_of_each_cell);
}

//! Set the number of nodes in each cell.
inline int *UnstructGrid::GetNodeNumberOfEachCell() const
{
    return cellTopology->GetNodeNumberOfEachCell();
}

inline void UnstructGrid::SetOversetConfig(OversetConfig *oversetConfigIn)
{
    this->oversetConfig = oversetConfigIn;
}

inline void UnstructGrid::SetAverageVolume(RDouble rIn)
{
    averageVolume = rIn;
}

inline int UnstructGrid::GetLI_nLine() const
{
    return this->lineTopology->GetLI_nLine();
}

inline int * UnstructGrid::GetLI_nCellsInLine() const
{
    return this->lineTopology->GetLI_nCellsInLine();
}

inline int ** UnstructGrid::GetLI_CellOfLine() const
{
    return this->lineTopology->GetLI_CellOfLine();
}

inline int ** UnstructGrid::GetLI_FaceOfLine() const
{
    return this->lineTopology->GetLI_FaceOfLine();
}

inline int * UnstructGrid::GetLI_LineOfCell() const
{
    return this->lineTopology->GetLI_LineOfCell();
}

#ifdef USE_INCOMSOLVER

inline Vector3D *UnstructGrid::GetCellCenterVector() const
{
    return cellCenterVector;
}

inline Vector3D *UnstructGrid::GetFaceNormalVector() const
{
    return faceNormalVector;
}

inline int *UnstructGrid::GetLocalCell2GlobalMap() const
{
    return localCell2GlobalMap;
}

inline int *UnstructGrid::GetLocalCell2LocalMap() const

{
    return localCell2LocalMap;
}


inline void UnstructGrid::SetLocalCell2GlobalMap(int *localCell2GlobalMap)
{
    this->localCell2GlobalMap = localCell2GlobalMap;
}

inline void UnstructGrid::SetLocalCell2LocalMap(int *localCell2LocalMap)
{
    this->localCell2LocalMap = localCell2LocalMap;
}

inline void UnstructGrid::SetCell2ZoneIDMap(int *cell2ZoneIDMap)
{
    this->cell2ZoneIDMap = cell2ZoneIDMap;
}

inline int* UnstructGrid::GetCell2ZoneIDMap() const
{
    return cell2ZoneIDMap;
}

inline void UnstructGrid::SetCellCenterVector(Vector3D *cellC)
{
    cellCenterVector = cellC;
}

inline void UnstructGrid::SetFaceNormalVector(Vector3D *faceN)
{
    faceNormalVector = faceN;
}

#endif

