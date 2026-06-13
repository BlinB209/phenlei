#include "Post_WriteVisualFile.h"
#include "Post_WriteTecplot.h"
#include "Post_WriteEnsight.h"
#include "Post_WriteParaview.h"
#include "Solver.h"
#include "Glb_Dimension.h"
#include "HyList.h"
#include "TK_Log.h"
#include "TK_Exit.h"

namespace PHSPACE
{

Post_WriteVisualFile *generateVisualFile = 0;

Post_WriteVisualFile::Post_WriteVisualFile()
{
    this->nodeMap.resize(0);
    this->gridMap.resize(0);
    this->VTKDataList.resize(0);
    this->boundaryName.resize(0);
    this->boundaryType.resize(0);
    this->boundaryGrid.resize(0);
    this->boundaryGridNeedConstruct = true;
    this->gridChanged = false;

    int nZones = PHMPI::GetNumberofGlobalZones();
    this->flowFieldDataOnNode.resize(nZones);
}

Post_WriteVisualFile::~Post_WriteVisualFile()
{
    for (int iGrid = 0; iGrid < this->boundaryGrid.size(); ++ iGrid)
    {
        delete this->boundaryGrid[iGrid];
        DelPointer2(this->nodeMap[iGrid]);
    }
}

void Post_WriteVisualFile::Run(int flowTypeIn)
{
    this->flowType = flowTypeIn;

    Initialize();
    ConstructBoundaryGrid();
    ComputeNodeData();
    StoreVTKData();
    StoreVisualizationData();
    WriteFile();
    ClearFieldData();
}

void Post_WriteVisualFile::BoundaryGridNeedConstruct()
{
    this->boundaryGridNeedConstruct = true;
}

void Post_WriteVisualFile::ConstructBoundaryGrid()
{
    if (this->boundaryGridNeedConstruct != true)
    {
        return;
    }

    for (int iGrid = 0; iGrid < this->boundaryGrid.size(); ++ iGrid)
    {
        delete this->boundaryGrid[iGrid];
        DelPointer2(this->nodeMap[iGrid]);
    }
    this->nodeMap.resize(0);
    this->gridMap.resize(0);
    this->boundaryName.resize(0);
    this->boundaryType.resize(0);
    this->boundaryGrid.resize(0);

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int myid = GetCurrentProcessorID();
    for (int iZone =  0; iZone < nZones; ++ iZone)
    {
        int iZoneProcess = GetZoneProcessorID(iZone);
        if (myid == iZoneProcess)
        {
            Grid *grid = GetGrid(iZone);
            int gridType = grid->Type();
            if (gridType == STRUCTGRID)
            {
                ConstructStrBoundaryGrid(grid);
            }
            else
            {
                ConstructUnsBoundaryGrid(grid);
            }
        }
    }

    this->boundaryGridNeedConstruct = false;
    this->gridChanged = true;
}

void Post_WriteVisualFile::ConstructStrBoundaryGrid(Grid *gridIn)
{
    StructGrid *StrGrid = StructGridCast(gridIn);
    int gridIndex = StrGrid->GetZoneID();

    StructBCSet *structBCSet = StrGrid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();
    int s_st[3] = {0}, s_ed[3] = {0};
    RDouble3D &xx = * StrGrid->GetStructX();
    RDouble3D &yy = * StrGrid->GetStructY();
    RDouble3D &zz = * StrGrid->GetStructZ();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int bctype = bcregion->GetBCType();
        if (IsInterface(bctype)) continue;

        string bcName = bcregion->GetBCName();

        int nDim = GetDim();
        GridID *index = new GridID(gridIndex);
        Grid *grid = new StructGrid();
        grid->InitGrid(index, 0, nDim, STRUCTGRID);
        StructGrid *strGrid = StructGridCast(grid);

        int *s_lr3d = bcregion->GetFaceDirectionIndex();
        for (int m = 0; m < nDim; ++ m)
        {
            s_st[m] = bcregion->GetStartPoint(m);
            s_ed[m] = bcregion->GetEndPoint(m);
            if (bcregion->GetFaceDirection() == m)
            {
                if (s_st[m] > 1 || s_lr3d[m] == 1)
                {
                    s_st[m] += 1;
                    s_ed[m] += 1;
                }
            }
            else
            {
                s_ed[m] += 1;
            }
        }

        int ist = s_st[0];
        int ied = s_ed[0];
        int jst = s_st[1];
        int jed = s_ed[1];
        int kst = s_st[2];
        int ked = s_ed[2];

        if (nDim == TWO_D)
        {
            kst = 1;
            ked = 1;
        }

        int ii = ied - ist + 1;
        int jj = jed - jst + 1;
        int kk = ked - kst + 1;
        strGrid->SetNI(ii);
        strGrid->SetNJ(jj);
        strGrid->SetNK(kk);

        int iCell = ii - 1 ? ii - 1 : 1;
        int jCell = jj - 1 ? jj - 1 : 1;
        int kCell = kk - 1 ? kk - 1 : 1;
        int nTotalCell = iCell * jCell * kCell;
        int nTotalNode = ii * jj * kk;
        strGrid->SetNTotalNode(nTotalNode);
        strGrid->SetNTotalCell(nTotalCell);

        RDouble *x = new RDouble [nTotalNode];
        RDouble *y = new RDouble [nTotalNode];
        RDouble *z = new RDouble [nTotalNode];

        int **currentNodeMap = NewPointer2<int>(3, nTotalNode);
        int nodeCount = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    x[nodeCount] = xx(i, j, k);
                    y[nodeCount] = yy(i, j, k);
                    z[nodeCount] = zz(i, j, k);

                    currentNodeMap[0][nodeCount] = i;
                    currentNodeMap[1][nodeCount] = j;
                    currentNodeMap[2][nodeCount] = k;
                    nodeCount ++;
                }
            }
        }
        strGrid->SetX(x);
        strGrid->SetY(y);
        strGrid->SetZ(z);

        boundaryGrid.push_back(grid);
        gridMap.push_back(gridIndex);
        nodeMap.push_back(currentNodeMap);
        boundaryName.push_back(bcName);
        boundaryType.push_back(bctype);
    }
}

void Post_WriteVisualFile::ConstructUnsBoundaryGrid(Grid *gridIn)
{
    UnstructGrid *UnsGrid = UnstructGridCast(gridIn);

    int gridIndex                 = UnsGrid->GetZoneID();
    int nTotalNode                = UnsGrid->GetNTotalNode();
    int nBoundFace                = UnsGrid->GetNBoundFace();
    int *node_number_of_each_face = UnsGrid->GetNodeNumberOfEachFace();

    RDouble *x = UnsGrid->GetX();
    RDouble *y = UnsGrid->GetY();
    RDouble *z = UnsGrid->GetZ();

    RDouble *xfc = UnsGrid->GetFaceCenterX();
    RDouble *yfc = UnsGrid->GetFaceCenterY();
    RDouble *zfc = UnsGrid->GetFaceCenterZ();

    set<pair<int, string> > bcnameMap;
    UnstructBCSet *unstructBCSet = UnsGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace      = unstructBCSet->GetBCFaceInfo();

    int nodeCount = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        string bcName = bcRegion->GetBCName();
        int bcType = bcRegion->GetBCType();

        bcnameMap.insert(pair<int, string>(bcType, bcName));
        nodeCount += node_number_of_each_face[iFace];
    }

    vector<vector<int> > face2nodelist(nBoundFace);
    vector<int> linkmap;
    for (set<pair<int, string> >::iterator iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter)
    {
        if ((*iter).first == PHENGLEI::INTERFACE || (*iter).first < 0) continue;

        face2nodelist.resize(0);
        linkmap.resize(0);

        GetFace2NodeList(UnsGrid, *iter, linkmap, face2nodelist);
        if (face2nodelist.size() == 0)
        {
            continue;
        }

        int bctype = (*iter).first;
        string bcName = (*iter).second;

        int nDim = GetDim();
        GridID *index = new GridID(gridIndex);
        Grid *grid = new UnstructGrid();
        grid->InitGrid(index, 0, nDim, UNSTRUCTGRID);
        UnstructGrid *unsGrid = UnstructGridCast(grid);

        int NumPts      = static_cast<int>(linkmap.size());
        int NumElements = static_cast<int>(face2nodelist.size());
        int NumFaces    = NumElements * 4;
        unsGrid->SetNTotalNode(NumPts);
        unsGrid->SetNTotalCell(NumElements);
        unsGrid->SetNTotalFace(NumFaces);

        int TotalNumFaceNodes_Rect = NumElements * 4;
        int *cell2node = new int [TotalNumFaceNodes_Rect];
        int *node_number_of_each_cell = new int [NumElements];

        int count = 0;
        for (size_t iFace = 0; iFace < face2nodelist.size(); ++ iFace)
        {
            node_number_of_each_cell[iFace] = 4;
            uint_t nodeNumber = face2nodelist[iFace].size();
            int index0 = face2nodelist[iFace][0];
            for (int iNode = 0; iNode < nodeNumber; ++ iNode)
            {
                cell2node[count] = face2nodelist[iFace][iNode];
                count ++;
            }

            for (uint_t iNode = nodeNumber; iNode < 4; ++ iNode)
            {
                cell2node[count] = index0;
                count ++;
            }
        }
        unsGrid->SetNodeNumberOfEachCell(node_number_of_each_cell);
        unsGrid->SetCell2Node(cell2node);

        RDouble *xCurrent = new RDouble[NumPts];
        RDouble *yCurrent = new RDouble[NumPts];
        RDouble *zCurrent = new RDouble[NumPts];
        int **currentNodeMap = NewPointer2<int>(1, NumPts);

        for (size_t iNode = 0; iNode < linkmap.size(); ++ iNode)
        {
            int nodeIndex = linkmap[iNode];
            currentNodeMap[0][iNode] = nodeIndex;

            if (nodeIndex < nTotalNode)
            {
                xCurrent[iNode] = x[nodeIndex];
                yCurrent[iNode] = y[nodeIndex];
                zCurrent[iNode] = z[nodeIndex];
            }
            else
            {
                int faceIndex = nodeIndex - nTotalNode;

                xCurrent[iNode] = xfc[faceIndex];
                yCurrent[iNode] = yfc[faceIndex];
                zCurrent[iNode] = zfc[faceIndex];
            }
        }
        unsGrid->SetX(xCurrent);
        unsGrid->SetY(yCurrent);
        unsGrid->SetZ(zCurrent);

        boundaryGrid.push_back(unsGrid);
        gridMap.push_back(gridIndex);
        nodeMap.push_back(currentNodeMap);
        boundaryName.push_back(bcName);
        boundaryType.push_back(bctype);
    }
}

void Post_WriteVisualFile::StoreVTKData()
{
    if (!VTKvisual)
    {
        return;
    }

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        DataContainer *cdata = new DataContainer();
        cdata->MoveToBegin();

        int send_proc = GetZoneProcessorID(iZone);
        int recv_proc = GetServerProcessorID();

        if (myid == send_proc)
        {
            StoreVTKData(iZone, cdata);
        }

        PH_Trade(cdata, send_proc, recv_proc, iZone);

        if (myid == recv_proc)
        {
            VTKDataList.push_back(cdata);
        }
        else
        {
            delete cdata;    cdata = nullptr;
        }
    }
}

void Post_WriteVisualFile::StoreVTKData(int zoneIndex, DataContainer *cdata)
{
    Grid *grid = flowFieldDataOnNode[zoneIndex]->GetGrid();
    if (grid->Type() == UNSTRUCTGRID)
    {
        int gridID = grid->GetZoneID();
        int dimension = grid->GetDim();
        int type = grid->Type();

        PHWrite(cdata, gridID);
        PHWrite(cdata, dimension);
        PHWrite(cdata, type);

        StoreUnsBoundaryVTKData(zoneIndex, cdata);

        if (WantVisualField(grid))
        {
            StoreUnsFieldVTKData(zoneIndex, cdata);
        }

        int TecioMission = Writecomplete;
        PHWrite(cdata, TecioMission);
    }
    else
    {
        int gridID = grid->GetZoneID();
        int dimension = grid->GetDim();
        int type = grid->Type();

        PHWrite(cdata, gridID);
        PHWrite(cdata, dimension);
        PHWrite(cdata, type);

        StoreStrBoundaryVTKData(zoneIndex, cdata);

        if (WantVisualField(grid))
        {
            StoreStrFieldVTKData(zoneIndex, cdata);
        }

        int TecioMission = Writecomplete;
        PHWrite(cdata, TecioMission);
    }
}

void Post_WriteVisualFile::StoreStrBoundaryVTKData(int zoneIndex, DataContainer *cdata)
{
    Post_Visual *FieldData = flowFieldDataOnNode[zoneIndex];
    int nvarplot = FieldData->GetNumberofVisualVariables();
    int nVisualVariables = nvarplot;
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nVisualVariables = nVisualVariables + 2 * numberOfSpecies;
    }
    RDouble4D **qn = new RDouble4D *[nVisualVariables];
    FieldData->GetAllVisualNodeVarPtr(qn);

    int TecioMission = WriteBoundary;
    for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
    {
        if (boundaryGrid[iBoundary]->GetZoneID() != zoneIndex)
        {
            continue;
        }

        int bctype    = boundaryType[iBoundary];
        string bcName = boundaryName[iBoundary];

        PHWrite(cdata, TecioMission);
        PHWrite(cdata, bctype);
        cdata->WriteString(bcName);

        int nTotalNode = boundaryGrid[iBoundary]->GetNTotalNode();
        int nTotalCell = boundaryGrid[iBoundary]->GetNTotalCell();
        PHWrite(cdata, nTotalNode);
        PHWrite(cdata, nTotalCell);

        RDouble *xBoundary = boundaryGrid[iBoundary]->GetX();
        RDouble *yBoundary = boundaryGrid[iBoundary]->GetY();
        RDouble *zBoundary = boundaryGrid[iBoundary]->GetZ();
        PHWrite(cdata, xBoundary, nTotalNode);
        PHWrite(cdata, yBoundary, nTotalNode);
        PHWrite(cdata, zBoundary, nTotalNode);

        for (int m = 0; m < nVisualVariables; ++ m)
        {
            RDouble *qBoundary = new RDouble[nTotalNode];

            for (int iNode = 0; iNode < nTotalNode; ++ iNode)
            {
                int originaI = nodeMap[iBoundary][0][iNode];
                int originaJ = nodeMap[iBoundary][1][iNode];
                int originaK = nodeMap[iBoundary][2][iNode];

                qBoundary[iNode] = (*qn[m])(originaI, originaJ, originaK, 0);
            }

            PHWrite(cdata, qBoundary, nTotalNode);
            delete [] qBoundary;
        }

        int ni = StructGridCast(boundaryGrid[iBoundary])->GetNI();
        int nj = StructGridCast(boundaryGrid[iBoundary])->GetNJ();
        int nk = StructGridCast(boundaryGrid[iBoundary])->GetNK();

        int iCellEnd = ni - 1 ? ni - 1 : 1;
        int jCellEnd = nj - 1 ? nj - 1 : 1;
        int kCellEnd = nk - 1 ? nk - 1 : 1;

        int ncount = 0;
        int *cell2node = new int[nTotalCell * 4];
        for (int kCell = 0; kCell < kCellEnd; ++ kCell)
        {
            for (int jCell = 0; jCell < jCellEnd; ++ jCell)
            {
                for (int iCell = 0; iCell < iCellEnd; ++ iCell)
                {
                    int nodeIndex1, nodeIndex2, nodeIndex3, nodeIndex4;
                    if (ni == 1)
                    {
                        nodeIndex1 = iCell + ni *  jCell      + ni * nj *  kCell;
                        nodeIndex2 = iCell + ni * (jCell + 1) + ni * nj *  kCell;
                        nodeIndex3 = iCell + ni * (jCell + 1) + ni * nj * (kCell + 1);
                        nodeIndex4 = iCell + ni *  jCell      + ni * nj * (kCell + 1);
                    }
                    else if (nj == 1)
                    {
                        nodeIndex1 =  iCell      + ni * jCell + ni * nj *  kCell;
                        nodeIndex2 = (iCell + 1) + ni * jCell + ni * nj *  kCell;
                        nodeIndex3 = (iCell + 1) + ni * jCell + ni * nj * (kCell + 1);
                        nodeIndex4 =  iCell      + ni * jCell + ni * nj * (kCell + 1);
                    }
                    else
                    {
                        nodeIndex1 =  iCell      + ni *  jCell      + ni * nj * kCell;
                        nodeIndex2 = (iCell + 1) + ni *  jCell      + ni * nj * kCell;
                        nodeIndex3 = (iCell + 1) + ni * (jCell + 1) + ni * nj * kCell;
                        nodeIndex4 =  iCell      + ni * (jCell + 1) + ni * nj * kCell;
                    }

                    cell2node[ncount++] = nodeIndex1 + 1;
                    cell2node[ncount++] = nodeIndex2 + 1;
                    cell2node[ncount++] = nodeIndex3 + 1;
                    cell2node[ncount++] = nodeIndex4 + 1;
                }
            }
        }

        PHWrite(cdata, cell2node, nTotalCell * 4);
        delete [] cell2node;
    }

    delete [] qn;
}

void Post_WriteVisualFile::StoreStrFieldVTKData(int zoneIndex, DataContainer *cdata)
{
    Grid *grid = flowFieldDataOnNode[zoneIndex]->GetGrid();
    StructGrid *strGrid = StructGridCast(grid);

    int ni = strGrid->GetNI();
    int nj = strGrid->GetNJ();
    int nk = strGrid->GetNK();

    RDouble3D &xx = *strGrid->GetStructX();
    RDouble3D &yy = *strGrid->GetStructY();
    RDouble3D &zz = *strGrid->GetStructZ();

    int nTotalCell = strGrid->GetNTotalCell();
    int nTotalNode = strGrid->GetNTotalNode();

    if (nk == 1)
    {
        int nVisualVariables = flowFieldDataOnNode[zoneIndex]->GetNumberofVisualVariables();
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        if (nchem)
        {
            int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
            nVisualVariables = nVisualVariables + 2 * numberOfSpecies;
        }
        RDouble4D **qn = new RDouble4D *[nVisualVariables];
        flowFieldDataOnNode[zoneIndex]->GetAllVisualNodeVarPtr(qn);

        int TecioMission = WriteBlock;
        PHWrite(cdata, TecioMission);

        PHWrite(cdata, nTotalNode);
        PHWrite(cdata, nTotalCell);

        RDouble *xBoundary = new RDouble[nTotalNode];
        RDouble *yBoundary = new RDouble[nTotalNode];
        RDouble *zBoundary = new RDouble[nTotalNode];

        RDouble **qBoundary = new RDouble *[nVisualVariables];
        for (int m = 0; m < nVisualVariables; ++ m)
        {
            qBoundary[m] = new RDouble[nTotalNode];
        }

        int ist, ied, jst, jed, kst, ked;
        strGrid->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

        int iNode = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    xBoundary[iNode] = xx(i, j, k);
                    yBoundary[iNode] = yy(i, j, k);
                    zBoundary[iNode] = zz(i, j, k);

                    for (int m = 0; m < nVisualVariables; ++ m)
                    {
                        qBoundary[m][iNode] = (*qn[m])(i, j, k, 0);
                    }

                    iNode ++ ;
                }
            }
        }

        PHWrite(cdata, xBoundary, nTotalNode);
        PHWrite(cdata, yBoundary, nTotalNode);
        PHWrite(cdata, zBoundary, nTotalNode);

        delete [] xBoundary;
        delete [] yBoundary;
        delete [] zBoundary;

        for (int m = 0; m < nVisualVariables; ++ m)
        {
            PHWrite(cdata, qBoundary[m], nTotalNode);
            delete [] qBoundary[m];
        }

        delete [] qBoundary;

        int *cell2node = new int[nTotalCell * 4];
        int ncount = 0;
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = 1; i <= ni - 1; ++ i)
            {
                cell2node[ncount++] = i     + (j - 1) * ni;
                cell2node[ncount++] = i + 1 + (j - 1) * ni;
                cell2node[ncount++] = i + 1 + (j   ) * ni;
                cell2node[ncount++] = i     + (j   ) * ni;
            }
        }

        PHWrite(cdata, cell2node, nTotalCell * 4);
        delete [] cell2node;
        delete [] qn;
    }
}

void Post_WriteVisualFile::StoreUnsBoundaryVTKData(int zoneIndex, DataContainer *cdata)
{
    Post_Visual *FieldData = flowFieldDataOnNode[zoneIndex];
    int nVisualVariables = FieldData->GetNumberofVisualVariables();
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nVisualVariables = nVisualVariables + 2 * numberOfSpecies;
    }

    RDouble **qn = new RDouble *[nVisualVariables];
    FieldData->GetAllVisualNodeVarPtr(qn);

    UnstructGrid *grid = UnstructGridCast(FieldData->GetGrid());
    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();
    int nBoundFace = grid->GetNBoundFace();
    int NumPts = grid->GetNTotalNode();

    RDouble **qFaceCenter = NewPointer2<RDouble>(nVisualVariables, nBoundFace);
    for (int m = 0; m < nVisualVariables; ++ m)
    {
        int nodeCount = 0;
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            RDouble value = 0.0;
            for (int iNode = 0; iNode < node_number_of_each_face[iFace]; ++ iNode)
            {
                int nodeIndex = face2node[nodeCount + iNode];
                value += qn[m][nodeIndex];
            }

            qFaceCenter[m][iFace] = value / node_number_of_each_face[iFace];
            nodeCount += node_number_of_each_face[iFace];
        }
    }

    int TecioMission = WriteBoundary;
    for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
    {
        if (boundaryGrid[iBoundary]->GetZoneID() != zoneIndex)
        {
            continue;
        }

        int bctype    = boundaryType[iBoundary];
        string bcName = boundaryName[iBoundary];

        PHWrite(cdata, TecioMission);
        PHWrite(cdata, bctype);
        cdata->WriteString(bcName);

        int nTotalNode = boundaryGrid[iBoundary]->GetNTotalNode();
        int nTotalCell = boundaryGrid[iBoundary]->GetNTotalCell();
        PHWrite(cdata, nTotalNode);
        PHWrite(cdata, nTotalCell);

        RDouble *xBoundary = boundaryGrid[iBoundary]->GetX();
        RDouble *yBoundary = boundaryGrid[iBoundary]->GetY();
        RDouble *zBoundary = boundaryGrid[iBoundary]->GetZ();
        PHWrite(cdata, xBoundary, nTotalNode);
        PHWrite(cdata, yBoundary, nTotalNode);
        PHWrite(cdata, zBoundary, nTotalNode);

        for (int m = 0; m < nVisualVariables; ++ m)
        {
            RDouble *qBoundary = new RDouble[nTotalNode];
            for (int iNode = 0; iNode < nTotalNode; ++ iNode)
            {
                int nodeIndex = nodeMap[iBoundary][0][iNode];

                if (nodeIndex < NumPts)
                {
                    qBoundary[iNode] = qn[m][nodeIndex];
                }
                else
                {
                    int faceIndex = nodeIndex - NumPts;
                    qBoundary[iNode] = qFaceCenter[m][faceIndex];
                }
            }

            PHWrite(cdata, qBoundary, nTotalNode);
            delete [] qBoundary;
        }

        int *cell2node = UnstructGridCast(boundaryGrid[iBoundary])->GetCell2Node();
        for (int iNode = 0; iNode < nTotalCell * 4; ++ iNode)
        {
            int nodeIndex = cell2node[iNode] + 1;
            PHWrite(cdata, nodeIndex);
        }
    }

    delete [] qn;    qn = nullptr;
    DelPointer2(qFaceCenter);
}

void Post_WriteVisualFile::StoreUnsFieldVTKData(int zoneIndex, DataContainer *cdata)
{
    if (GetDim() == PHSPACE::TWO_D)
    {
        StoreUnsFieldVTKData2D(zoneIndex, cdata);
    }
    else
    {
        StoreUnsFieldVTKData3D(zoneIndex, cdata);
    }
}

void Post_WriteVisualFile::StoreUnsFieldVTKData2D(int zoneIndex, DataContainer *cData)
{
    int nVisualVariables = flowFieldDataOnNode[zoneIndex]->GetNumberofVisualVariables();
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nVisualVariables = nVisualVariables + 2 * numberOfSpecies;
    }
    RDouble **qn = new RDouble *[nVisualVariables];
    flowFieldDataOnNode[zoneIndex]->GetAllVisualNodeVarPtr(qn);

    Grid *grid = flowFieldDataOnNode[zoneIndex]->GetGrid();
    UnstructGrid *unsGrid = UnstructGridCast(grid);

    int tecioMission = WriteBlock;
    PHWrite(cData, tecioMission);

    int nTotalNode = unsGrid->GetNTotalNode();
    int nTotalCell = unsGrid->GetNTotalCell();
    int nBoundFace = unsGrid->GetNBoundFace();
    int nTotalFace = unsGrid->GetNTotalFace();

    int *face2Node       = unsGrid->GetFace2Node();
    int *leftCellOfFace  = unsGrid->GetLeftCellOfFace();
    int *rightCellOfFace = unsGrid->GetRightCellOfFace();

    RDouble *x = unsGrid->GetX();
    RDouble *y = unsGrid->GetY();
    RDouble *z = unsGrid->GetZ();

    HyList<int> * cell2Face = new HyList <int> [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2Face[iCell].SetAverNnode(4);
    }
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        cell2Face[le].insert(iFace);
    }
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        cell2Face[le].insert(iFace);
        cell2Face[re].insert(iFace);
    }

    //! Find out cellToNode.
    HyList<int> *cell2Node = new HyList <int> [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2Node[iCell].SetAverNnode(4);
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int face0 = cell2Face[iCell][0];
        int ip1   = face2Node[2 * face0];
        int ip2   = face2Node[2 * face0 + 1];
        cell2Node[iCell].insert(ip1);
        cell2Node[iCell].insert(ip2);

        int targetNode = ip2;
        int checkNode  = ip1;

        int nFaceOfCell = cell2Face[iCell].size();
        int isClockWise = false;
        while (! isClockWise)
        {
            for (int iFace = 1; iFace < nFaceOfCell; ++ iFace)
            {
                int face  = cell2Face[iCell][iFace];
                int node1 = face2Node[2 * face];
                int node2 = face2Node[2 * face + 1];

                if (node1 == targetNode)
                {
                    if (node2 != checkNode)
                    {
                        cell2Node[iCell].insert(node2);
                        targetNode = node2;
                    }
                    else
                    {
                        isClockWise = true;
                        break;
                    }
                }
                else if (node2 == targetNode)
                {
                    if (node1 != checkNode)
                    {
                        cell2Node[iCell].insert(node1);
                        targetNode = node1;
                    }
                    else
                    {
                        isClockWise = true;
                        break;
                    }
                }
            }    //! Iterate over all faces of cell.
        }        //! Judge whether the nodes of cell have searched closed.
    }            //! iCell.

    int nTotalCellVTK = 0;
    vector<RDouble> xVTK;
    vector<RDouble> yVTK;
    vector<RDouble> zVTK;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        RDouble xNode = x[iNode];
        RDouble yNode = y[iNode];
        RDouble zNode = z[iNode];
        xVTK.push_back(xNode);
        yVTK.push_back(yNode);
        zVTK.push_back(zNode);
    }

    UnstructGrid *gridVTK = UnstructGridCast(grid);
    RDouble      *xcc     = gridVTK->GetCellCenterX();
    RDouble      *ycc     = gridVTK->GetCellCenterY();
    RDouble      *zcc     = gridVTK->GetCellCenterZ();

    vector<vector<RDouble>> qnVTK;
    qnVTK.resize(nVisualVariables);
    for (int iVisualVariable = 0; iVisualVariable < nVisualVariables; ++ iVisualVariable)
    {
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            qnVTK[iVisualVariable].push_back(qn[iVisualVariable][iNode]);
        }
    }

    vector<int> cell2NodeVTK;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (cell2Node[iCell].size() == 3)
        {
            for (int m = 0; m < 3; ++ m)
            {
                cell2NodeVTK.push_back((cell2Node[iCell].GetData(m) + 1));
            }
            cell2NodeVTK.push_back((cell2Node[iCell].GetData(2) + 1));
            ++ nTotalCellVTK;
        }
        else if (cell2Node[iCell].size() == 4)
        {
            for (int m = 0; m < 4; ++ m)
            {
                cell2NodeVTK.push_back((cell2Node[iCell].GetData(m) + 1));
            }
            ++ nTotalCellVTK;
        }
        else
        {
            int nDividedCell = cell2Node[iCell].size();
            int countNodeID = 0;
            xVTK.push_back(xcc[iCell]);
            yVTK.push_back(ycc[iCell]);
            zVTK.push_back(zcc[iCell]);
            ++ nTotalNode;
            int cellCenterNodeID = nTotalNode - 1;
            RDouble qnVTKcentre;
            for (int iVisualVariable = 0; iVisualVariable < nVisualVariables; ++ iVisualVariable)
            {
                RDouble qnVTKcount = 0;
                for (int iNode = 0; iNode < nDividedCell; ++ iNode)
                {
                    int nodeID = cell2Node[iCell].GetData(iNode);
                    qnVTKcount = qnVTKcount + qnVTK[iVisualVariable][nodeID];
                }
                qnVTKcentre = qnVTKcount / nDividedCell;
                qnVTK[iVisualVariable].push_back(qnVTKcentre);
            }

            for (int iDividedCell = 0; iDividedCell < nDividedCell; ++ iDividedCell)
            {
                for (int iNode = 0; iNode < 2; ++ iNode)
                {
                    cell2NodeVTK.push_back((cellCenterNodeID + 1));
                }
                cell2NodeVTK.push_back((cell2Node[iCell].GetData(iDividedCell) + 1));
                if ((iDividedCell + 1) < nDividedCell)
                {
                    cell2NodeVTK.push_back((cell2Node[iCell].GetData(iDividedCell + 1) + 1));
                }
                else
                {
                    cell2NodeVTK.push_back((cell2Node[iCell].GetData(0) + 1));
                }
                ++ nTotalCellVTK;
            }
        }
    }

    PHWrite(cData, nTotalNode);
    PHWrite(cData, nTotalCellVTK);

    PHWrite(cData, &xVTK[0], nTotalNode);
    PHWrite(cData, &yVTK[0], nTotalNode);
    PHWrite(cData, &zVTK[0], nTotalNode);

    for (int iVisualVariable = 0; iVisualVariable < nVisualVariables; ++ iVisualVariable)
    {
        PHWrite(cData, &qnVTK[iVisualVariable][0], nTotalNode);
    }

    PHWrite(cData, &cell2NodeVTK[0], nTotalCellVTK * 4);

    delete [] qn;           qn        = nullptr;
    delete [] cell2Node;    cell2Node = nullptr;
    delete [] cell2Face;    cell2Face = nullptr;
}

void Post_WriteVisualFile::StoreUnsFieldVTKData3D(int zoneIndex, DataContainer *cdata)
{
    int nVisualVariables = flowFieldDataOnNode[zoneIndex]->GetNumberofVisualVariables();
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nVisualVariables = nVisualVariables + 2 * numberOfSpecies;
    }
    RDouble **qn = new RDouble *[nVisualVariables];
    flowFieldDataOnNode[zoneIndex]->GetAllVisualNodeVarPtr(qn);

    int TecioMission = WriteBlock;
    PHWrite(cdata, TecioMission);

    UnstructGrid *grid = UnstructGridCast(flowFieldDataOnNode[zoneIndex]->GetGrid());
    int nTotalNode = grid->GetNTotalNode();
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();
    int *nodeNumberOfEachCell = grid->GetNodeNumberOfEachCell();
    int *cell2Node = grid->GetCell2Node();

    PHWrite(cdata, nTotalNode);
    PHWrite(cdata, nTotalCell);

    PHWrite(cdata, x, nTotalNode);
    PHWrite(cdata, y, nTotalNode);
    PHWrite(cdata, z, nTotalNode);

    for (int m = 0; m < nVisualVariables; ++ m)
    {
        RDouble *qq = qn[m];
        PHWrite(cdata, qq, nTotalNode);
    }

    PHWrite(cdata, nTotalFace);

    int count = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        count += nodeNumberOfEachCell[iCell];
    }
    PHWrite(cdata, nodeNumberOfEachCell, nTotalCell);
    PHWrite(cdata, cell2Node           , count);

    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();
    int **cell2Face = grid->GetCell2Face();

    PHWrite(cdata, faceNumberOfEachCell, nTotalCell);
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        PHWrite(cdata, cell2Face[iCell], faceNumberOfEachCell[iCell]);
    }

    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
    int *face2Node = grid->GetFace2Node();

    count = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        count += nodeNumberOfEachFace[iFace];
    }
    PHWrite(cdata, nodeNumberOfEachFace, nTotalFace);
    PHWrite(cdata, face2Node, count);

    delete [] qn;
}

void Post_WriteVisualFile::ComputeNodeData()
{
    int nVisualVariables = GlobalDataBase::GetIntParaFromDB("nVisualVariables");
    int visualVariables[100];
    GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int myid = GetCurrentProcessorID();
    for (int iZone =  0; iZone < nZones; ++ iZone)
    {
        Post_Visual *postVisualization = 0;
        int iZoneProcess = GetZoneProcessorID(iZone);
        if (myid == iZoneProcess)
        {
            Grid *grid = GetGrid(iZone);
            postVisualization = new Post_Visual(nVisualVariables, visualVariables, grid, flowType);

            int Nsolver = GlobalSolvers::GetNumberOfSolvers(iZone);
            for (int iSolver = 0; iSolver < Nsolver; ++ iSolver)
            {
                GlobalSolvers::GetSolver(iZone, iSolver)->ComputePostVisualVariables(postVisualization);
            }
        }

        flowFieldDataOnNode[iZone] = postVisualization;
    }
}

bool WantVisualField(Grid *grid)
{
    int plotFieldType = 0;
    GlobalDataBase::GetData("plotFieldType", &plotFieldType, PHINT, 1);

    if (plotFieldType == TEC_SPACE::BlockVisual)
    {
        int dimension = GetDim();

        RDouble *pmin = grid->GetMinBox();
        RDouble *pmax = grid->GetMaxBox();

        RDouble lowerPlotFieldBox[3], upperPlotFieldBox[3];
        GlobalDataBase::GetData("lowerPlotFieldBox", &lowerPlotFieldBox, PHDOUBLE, 3);
        GlobalDataBase::GetData("upperPlotFieldBox", &upperPlotFieldBox, PHDOUBLE, 3);

        if ((pmin[0] > upperPlotFieldBox[0]) || (pmax[0] < lowerPlotFieldBox[0])) return false;

        if ((pmin[1] > upperPlotFieldBox[1]) || (pmax[1] < lowerPlotFieldBox[1])) return false;

        if (dimension == THREE_D)
        {
            if ((pmin[2] > upperPlotFieldBox[2]) || (pmax[2] < lowerPlotFieldBox[2])) return false;
        }
    }
    return grid->GetDim() == TWO_D || (grid->GetDim() == THREE_D && plotFieldType != TEC_SPACE::BoundaryVisual);
}

void WriteVisualFile(int flowTypeIn)
{
    int dumpFlowOnOriginalGrid = 0;
    GlobalDataBase::GetData("dumpFlowOnOriginalGrid", &dumpFlowOnOriginalGrid, PHINT, 1);
    if (dumpFlowOnOriginalGrid)
    {
        Post_WriteTecplotByOriginalGrid *writeTecplotByOriginalGrid = GetWriteTecplotByOriginalGrid();
        writeTecplotByOriginalGrid->Run();
    }
    else
    {
    if (!generateVisualFile)
    {
        int visualfileType = GlobalDataBase::GetIntParaFromDB("visualfileType");
        if (visualfileType == Tecplot_Binary || visualfileType == Tecplot_ASCII)
        {
            generateVisualFile = new Post_WriteTecplot();
        }
        else if (visualfileType == Ensight_Binary || visualfileType == Ensight_ASCII)
        {
            generateVisualFile = new Post_WriteEnsight();
        }
        else if (visualfileType == Paraview)
        {
            generateVisualFile = new Post_WriteParaview();
        }
        else
        {
            TK_Exit::UnexpectedVarValue("visualfileType", visualfileType);
        }
    }
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    if (1 == isAle)
    {
        generateVisualFile->BoundaryGridNeedConstruct();
    }
    generateVisualFile->Run(flowTypeIn);
}

}

}