#include "Post_WriteTecplot.h"
#include "PHIO.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
#include "IO_FileName.h"
#include "Post_Visual.h"
#include "Glb_Dimension.h"
#include "HyList.h"
#include "Solver.h"
#include "Post_WriteVisualFile.h"

namespace PHSPACE
{

Post_WriteTecplotByOriginalGrid* writeTecplotByOriginalGrid = new Post_WriteTecplotByOriginalGrid();

Post_WriteTecplot::Post_WriteTecplot()
{
    this->visualDataList.resize(0);
}

Post_WriteTecplot::~Post_WriteTecplot()
{

}

void Post_WriteTecplot::Initialize()
{
    string visualfile = "tecflow.dat";
    GlobalDataBase::GetData("visualfile", &visualfile, PHSTRING, 1);

    string fname, fext;
    GetNameExt(visualfile, fname, fext, ".");

    VTKFileName = visualfile;
    visualFileName = visualfile;
    VTKvisual = false;

    if (flowType == AverageFlow)
    {
        visualFileName = AddSymbolToFileName(visualfile,"_Average");
    }

    if (flowType == AverageReynoldsStress)
    {
        visualFileName = AddSymbolToFileName(visualfile,"_ReynoldsStress");
    }

    string sentinelfilename = "results/sentinel";
    fstream sentinelfile;
    sentinelfile.open(sentinelfilename.c_str(), ios::in);

    if (!sentinelfile && !visualfile.empty() && (flowType == 0))
    {
        VTKFileName = visualfile + ".bak";
        VTKvisual = true;
    }

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady)
    {
        int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
        visualFileName = AddSymbolToFileName(visualFileName, '_', outnstep);
    }
}

void Post_WriteTecplot::StoreVisualizationData()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        DataContainer *cdata = new DataContainer();
        DataContainer *cBlkData = new DataContainer();    //! Seperate surface and block data into two different cdata elements.
        cdata->MoveToBegin();
        cBlkData->MoveToBegin();

        int send_proc = GetZoneProcessorID(iZone);
        int recv_proc = GetZoneFileID(iZone);

        if (myid == send_proc)
        {
            StoreVisualizationData(iZone, cdata, cBlkData);
        }

        PH_Trade(cdata, send_proc, recv_proc, iZone);
        PH_Trade(cBlkData, send_proc, recv_proc, iZone);

        if (myid == recv_proc)
        {
            visualDataList.push_back(cdata);
            visualDataList.push_back(cBlkData);
        }
        else
        {
            delete cdata;       cdata = nullptr;
            delete cBlkData;    cBlkData = nullptr;
        }
    }
}

void Post_WriteTecplot::StoreVisualizationData(int zoneIndex, DataContainer *cdata, DataContainer *cBlkData)
{
    Grid *grid = flowFieldDataOnNode[zoneIndex]->GetGrid();

    int gridID = grid->GetZoneID();
    int dimension = grid->GetDim();
    int type = grid->Type();

    PHWrite(cdata, gridID);
    PHWrite(cdata, dimension);
    PHWrite(cdata, type);
    PHWrite(cBlkData, gridID);
    PHWrite(cBlkData, dimension);
    PHWrite(cBlkData, type);

    if (grid->Type() == UNSTRUCTGRID)
    {
        StoreUnsBoundaryVisualData(zoneIndex, cdata);

        if (WantVisualField(grid))
        {
            StoreUnsFieldVisualData(zoneIndex, cBlkData);
        }
    }
    else
    {
        StoreStrBoundaryVisualData(zoneIndex, cdata);

        if (WantVisualField(grid))
        {
            StoreStrFieldVisualData(zoneIndex, cBlkData);
        }
    }

    int TecioMission = Writecomplete;
    PHWrite(cdata, TecioMission);
    PHWrite(cBlkData, TecioMission);
}

void Post_WriteTecplot::StoreStrBoundaryVisualData(int zoneIndex, DataContainer *cdata)
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

        int ni = StructGridCast(boundaryGrid[iBoundary])->GetNI();
        int nj = StructGridCast(boundaryGrid[iBoundary])->GetNJ();
        int nk = StructGridCast(boundaryGrid[iBoundary])->GetNK();

        PHWrite(cdata, ni);
        PHWrite(cdata, nj);
        PHWrite(cdata, nk);

        int nTotalNode = boundaryGrid[iBoundary]->GetNTotalNode();
        PHWrite(cdata, nTotalNode);

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
    }

    delete [] qn;
}

void Post_WriteTecplot::StoreStrFieldVisualData(int zoneIndex, DataContainer *cdata)
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

    Grid *grid = flowFieldDataOnNode[zoneIndex]->GetGrid();
    StructGrid *strGrid = StructGridCast(grid);

    int TecioMission = WriteBlock;
    PHWrite(cdata, TecioMission);

    int ni = strGrid->GetNI();
    int nj = strGrid->GetNJ();
    int nk = strGrid->GetNK();

    PHWrite(cdata, ni);
    PHWrite(cdata, nj);
    PHWrite(cdata, nk);

    int nTotalNode = strGrid->GetNTotalNode();
    PHWrite(cdata, nTotalNode);

    RDouble *x = new RDouble[nTotalNode];
    RDouble *y = new RDouble[nTotalNode];
    RDouble *z = new RDouble[nTotalNode];

    RDouble3D &xx = *strGrid->GetStructX();
    RDouble3D &yy = *strGrid->GetStructY();
    RDouble3D &zz = *strGrid->GetStructZ();

    int count = 0;
    for (int k = 1; k <= nk; ++ k)
    {
        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                x[count] = xx(i, j, k);
                y[count] = yy(i, j, k);
                z[count] = zz(i, j, k);
                count ++;
            }
        }
    }

    PHWrite(cdata, x, nTotalNode);
    PHWrite(cdata, y, nTotalNode);
    PHWrite(cdata, z, nTotalNode);

    delete [] x;
    delete [] y;
    delete [] z;

    count = 0;
    for (int m = 0; m < nVisualVariables; ++ m)
    {
        RDouble *qq = new RDouble[nTotalNode];
        for (int k = 1; k <= nk; ++ k)
        {
            for (int j = 1; j <= nj; ++ j)
            {
                for (int i = 1; i <= ni; ++ i)
                {
                    qq[count] = (*qn[m])(i, j, k, 0);
                    count ++;
                }
            }
        }
        PHWrite(cdata, qq, nTotalNode);
        delete [] qq;
        count = 0;
    }

    delete [] qn;
}

void Post_WriteTecplot::StoreUnsBoundaryVisualData(int zoneIndex, DataContainer *cdata)
{
    Post_Visual *FieldData = flowFieldDataOnNode[zoneIndex];
    int nvarplot         = FieldData->GetNumberofVisualVariables();
    int nVisualVariables = nvarplot;
    int nchem            = GlobalDataBase::GetIntParaFromDB("nchem");
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
        int nTotalFace = boundaryGrid[iBoundary]->GetNTotalFace();
        int nTotalCell = boundaryGrid[iBoundary]->GetNTotalCell();
        PHWrite(cdata, nTotalNode);
        PHWrite(cdata, nTotalCell);
        PHWrite(cdata, nTotalFace);

        int TotalNumFaceNodes_Rect = nTotalCell * 4;
        PHWrite(cdata, &TotalNumFaceNodes_Rect, 1);

        int *node_number_of_each_cell = UnstructGridCast(boundaryGrid[iBoundary])->GetNodeNumberOfEachCell();
        PHWrite(cdata, node_number_of_each_cell, nTotalCell);

        PHWrite(cdata, nTotalNode);
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
            delete [] qBoundary;    qBoundary = nullptr;
        }

        int *cell2node = UnstructGridCast(boundaryGrid[iBoundary])->GetCell2Node();
        PHWrite(cdata, cell2node, TotalNumFaceNodes_Rect);
    }

    delete [] qn;    qn = nullptr;
    DelPointer2(qFaceCenter);
}

void Post_WriteTecplot::StoreUnsFieldVisualData(int zoneIndex, DataContainer *cdata)
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

    UnstructGrid *grid = UnstructGridCast(flowFieldDataOnNode[zoneIndex]->GetGrid());
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfFaces = grid->GetNTotalFace();

    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer  = grid->GetFace2Node();

    int *leftCellIndexContainer  = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    int totalNumFaceNodes = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        totalNumFaceNodes += faceNodeNumberContainer[iFace];
    }

    int TecioMission = WriteBlock;
    PHWrite(cdata, TecioMission);
    PHWrite(cdata, numberOfNodes);
    PHWrite(cdata, numberOfCells);
    PHWrite(cdata, numberOfFaces);
    PHWrite(cdata, totalNumFaceNodes);
    PHWrite(cdata, faceNodeNumberContainer, numberOfFaces);

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    PHWrite(cdata, numberOfNodes);
    PHWrite(cdata, x, numberOfNodes);
    PHWrite(cdata, y, numberOfNodes);
    PHWrite(cdata, z, numberOfNodes);

    for (int iVarible = 0; iVarible < nVisualVariables; ++ iVarible)
    {
        RDouble *qq = qn[iVarible];
        PHWrite(cdata, qq, numberOfNodes);
    }

    PHWrite(cdata, faceNodeIndexContainer, totalNumFaceNodes);

    PHWrite(cdata, leftCellIndexContainer , numberOfFaces);
    PHWrite(cdata, rightCellIndexContainer, numberOfFaces);

    delete [] qn;    qn = nullptr;
}

void Post_WriteTecplot::WriteFile()
{
    int intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
    int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (((outIterStep + 1) % intervalStepPlot == 0) && (outIterStep % intervalStepPlot != 0))
    {
        return;
    }

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if (VTKvisual == true && myid == server)
    {
        bool CharacteristicBoundary = false;
        DumpToVTK DumpToVTK(VTKDataList, VTKFileName, CharacteristicBoundary, this->flowType);
        DumpToVTK.Run();

        WriteSentinelFile();
    }

    if (visualDataList.size() > 0)
    {
        int visualfileType = GlobalDataBase::GetIntParaFromDB("visualfileType");
        if (visualfileType == Tecplot_Binary)
        {
#ifdef USE_TecplotLib
            DumpToTecio dumptotecio(visualDataList, visualFileName, this->flowType);
            dumptotecio.Run();
#else
            ostringstream WarningInfo;
            WarningInfo << "Warning : visualfileType = 0 is not right, dump flow data into ASCII file! \n";
            PrintToWindow(WarningInfo);

            DumpToTecplotASCII dumpToTecplotASCII(visualDataList, visualFileName, this->flowType);
            dumpToTecplotASCII.Run();
#endif
        }
        else if (visualfileType == Tecplot_ASCII)
        {
            DumpToTecplotASCII dumpToTecplotASCII(visualDataList, visualFileName, this->flowType);
            dumpToTecplotASCII.Run();
        }
        else
        {
            TK_Exit::UnexpectedVarValue("visualfileType", visualfileType);
        }
    }
}

void Post_WriteTecplot::ClearFieldData()
{
    for (int iData = 0; iData < flowFieldDataOnNode.size(); ++ iData)
    {
        if (flowFieldDataOnNode[iData])
        {
            FreePointer(flowFieldDataOnNode[iData]);
        }
    }

    for (unsigned int iZone = 0; iZone < VTKDataList.size(); ++ iZone)
    {
        delete VTKDataList[iZone];
    }
    VTKDataList.resize(0);

    for (unsigned int iZone = 0; iZone < visualDataList.size(); ++ iZone)
    {
        delete visualDataList[iZone];
    }
    visualDataList.resize(0);
}

Post_WriteTecplotByOriginalGrid::Post_WriteTecplotByOriginalGrid()
{
    nOriginalGrids              = 0;
    OriginalGrid                = 0;
    originalGridProc            = 0;
    originalGridIndexOfEachGrid = 0;
    originalGridData.resize(0);
}

Post_WriteTecplotByOriginalGrid::~Post_WriteTecplotByOriginalGrid()
{
    delete [] originalGridProc;               originalGridProc            = nullptr;
    delete [] originalGridIndexOfEachGrid;    originalGridIndexOfEachGrid = nullptr;
    for (int iZone = 0; iZone < nOriginalGrids; ++ iZone)
    {
        delete originalGridData[iZone];
        delete OriginalGrid[iZone];
    }
    delete [] OriginalGrid;    OriginalGrid = nullptr;
}

void Post_WriteTecplotByOriginalGrid::Run()
{
    string visualFile = "tecflow.dat";
    GlobalDataBase::GetData("visualfile", &visualFile, PHSTRING, 1);
    if (visualFile == "")
    {
        return;
    }

    ActionKey *actKey = new ActionKey();
    actKey->filename    = visualFile + ".bak";
    actKey->tecfilename = visualFile;
    actKey->VTKvisual   = false;

    string sentinelFilename = "results/sentinel";
    fstream sentinelFile;
    sentinelFile.open(sentinelFilename.c_str(), ios::in);
    if (!sentinelFile)
    {
        actKey->VTKvisual = true;
    }

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady)
    {
        int outNStep        = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
        actKey->tecfilename = PHSPACE::AddSymbolToFileName(actKey->tecfilename, '_', outNStep);
    }

    ServerWrite(actKey);
    PostWrite(actKey);
    delete actKey;    actKey = nullptr;
}

void Post_WriteTecplotByOriginalGrid::Initialize()
{
    int myid             = PHMPI::GetCurrentProcessorID();
    int nVisualVariables = GlobalDataBase::GetIntParaFromDB("nVisualVariables");
    int visualVariables[100];
    GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);
    originalGridData.resize(nOriginalGrids);

    for (int iZone = 0; iZone < nOriginalGrids; ++ iZone)
    {
        Post_Visual *postVisualization = 0;
        if (myid == originalGridProc[iZone])
        {
            postVisualization = new Post_Visual(nVisualVariables, visualVariables, OriginalGrid[iZone]);
            set<int> visualVariables = postVisualization->GetVisualVariables();

            for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
            {
                int    variableType = *varIter;
                string varName      = postVisualization->GetVariableName(variableType);
                int    gridType     = OriginalGrid[iZone]->Type();
                if (gridType == UNSTRUCTGRID)
                {
                    int nTotalNode = OriginalGrid[iZone]->GetNTotalNode();
                    RDouble *qn = new RDouble [nTotalNode];
                    postVisualization->UpdateDataPtr(varName, qn);
                }
                else
                {
                    StructGrid *gridStr = StructGridCast(OriginalGrid[iZone]);
                    int ni              = gridStr->GetNI();
                    int nj              = gridStr->GetNJ();
                    int nk              = gridStr->GetNK();
                    Range I(-1, ni + 1);
                    Range J(-1, nj + 1);
                    Range K(-1, nk + 1);
                    if (nk == 1) K.setRange(1, 1);
                    Range M(0, 0);
                    RDouble4D *qn = new RDouble4D(I, J, K, M, fortranArray);
                    postVisualization->UpdateDataPtr(varName, qn);
                }
            }
        }
        originalGridData[iZone] = postVisualization;
    }
}

void Post_WriteTecplotByOriginalGrid::ServerWrite(ActionKey *actKey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int myid   = GetCurrentProcessorID();
    int recvProc;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        DataContainer *cData = actKey->GetData();
        cData->MoveToBegin();
        int originalGridIndex = originalGridIndexOfEachGrid[iZone];
        int sendProc          = GetZoneProcessorID(iZone);
        recvProc = GetOriginalGridProc(originalGridIndex);
        int tag = GetSendRecvTag(actKey, iZone);

        if (myid == sendProc)
        {
            PHWrite(cData, &originalGridIndex, 1);
            ComputerNodeDataOnPartitionGrid(iZone, actKey);
        }
        PH_Trade(actKey, sendProc, recvProc, tag);
        if (myid == recvProc)
        {
            CollectNodeData(actKey);
        }
    }
    StoreVisualizationVar(actKey);
}

void Post_WriteTecplotByOriginalGrid::StoreVisualizationVar(ActionKey *actKey)
{
    int myid = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < nOriginalGrids; ++ iZone)
    {
        int           sendProc = GetOriginalGridProc(iZone);
        int           recvProc = PHMPI::GetServerProcessorID();
        int           tag      = GetSendRecvTag(actKey, iZone);
        DataContainer *cData   = actKey->GetData();
        cData->MoveToBegin();
        DataContainer *tecData = actKey->GetTecData();
        tecData->MoveToBegin();

        if (myid == sendProc)
        {
            Post_Visual *postVisualization = originalGridData[iZone];
            Grid        *grid              = postVisualization->GetGrid();
            if (grid->Type() == UNSTRUCTGRID)
            {
                int nVarPlot         = postVisualization->GetNumberofVisualVariables();
                int nVisualVariables = nVarPlot;
                int nChem            = GlobalDataBase::GetIntParaFromDB("nchem");
                if (nChem)
                {
                    int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
                    nVisualVariables    = nVisualVariables + 2 * numberOfSpecies;
                }
                RDouble **qn = new RDouble * [nVisualVariables];
                postVisualization->GetAllVisualNodeVarPtr(qn);

                if (true == actKey->VTKvisual)
                {
                    int GridID = grid->GetGridID()->GetIndex();
                    PHWrite(cData, GridID);
                    int dimension = grid->GetDim();
                    PHWrite(cData, dimension);
                    int type = grid->Type();
                    PHWrite(cData, &type, 1);
                    BoundaryVTKVisualization(grid, cData, qn, nVisualVariables);
                    if (WantVisualField(grid))
                    {
                        FieldVisualization(grid, cData, qn, nVisualVariables);
                    }
                    int TecioMission = TEC_SPACE::Writecomplete;
                    PHWrite(cData, &TecioMission, 1);
                }

                if (actKey->tecfilename != "")
                {
                    int GridID = grid->GetGridID()->GetIndex();
                    PHWrite(tecData, GridID);
                    int dimension = grid->GetDim();
                    PHWrite(tecData, dimension);
                    int type = grid->Type();
                    PHWrite(tecData, &type, 1);
                    BoundaryVisualization(grid, tecData, qn, nVisualVariables);
                    if (WantVisualField(grid))
                    {
                        SaveDataForTecio(grid, tecData, qn, nVisualVariables);
                    }
                    int TecioMission = TEC_SPACE::Writecomplete;
                    PHWrite(tecData, &TecioMission, 1);
                }
                delete [] qn;    qn = nullptr;
            }
            else if (STRUCTGRID == grid->Type())
            {
                int nVarPlot         = postVisualization->GetNumberofVisualVariables();
                int nVisualVariables = nVarPlot;
                int nChem            = GlobalDataBase::GetIntParaFromDB("nchem");
                if (nChem)
                {
                    int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
                    nVisualVariables    = nVisualVariables + numberOfSpecies;
                }
                RDouble4D **qn = new RDouble4D * [nVisualVariables];
                postVisualization->GetAllVisualNodeVarPtr(qn);

                if (true == actKey->VTKvisual)
                {
                    int GridID = grid->GetGridID()->GetIndex();
                    PHWrite(cData, GridID);
                    int dimension = grid->GetDim();
                    PHWrite(cData, dimension);
                    int type = grid->Type();
                    PHWrite(cData, &type, 1);
                    BoundaryVTKVisualization(grid, cData, qn, nVisualVariables);
                    if (WantVisualField(grid))
                    {
                        FieldVisualization(grid, cData, qn, nVisualVariables);
                    }
                    int TecioMission = TEC_SPACE::Writecomplete;
                    PHWrite(cData, &TecioMission, 1);
                }

                if (actKey->tecfilename != "")
                {
                    int GridID = grid->GetGridID()->GetIndex();
                    PHWrite(tecData, GridID);
                    int dimension = grid->GetDim();
                    PHWrite(tecData, dimension);
                    int type = grid->Type();
                    PHWrite(tecData, &type, 1);
                    BoundaryVisualization(grid, tecData, qn, nVisualVariables);
                    if (WantVisualField(grid))
                    {
                        SaveDataForTecio(grid, tecData, qn, nVisualVariables);
                    }
                    int TecioMission = TEC_SPACE::Writecomplete;
                    PHWrite(tecData, &TecioMission, 1);
                }
                delete [] qn;    qn = nullptr;
            }
            else
            {
                TK_Exit::ExceptionExit("Error!!! Task_GenerateVisualFile::Visualization grid type not right!");
            }
        }
        PHMPI::PH_Trade(actKey, sendProc, recvProc, tag);
        PHMPI::PH_Trade(tecData, sendProc, recvProc, tag);
        if (myid == recvProc)
        {
            vtkDataList.push_back(cData);
            DataContainer *cData0 = new DataContainer();
            actKey->SetData(cData0);
            dataList.push_back(tecData);
            DataContainer *cData1 = new DataContainer();
            actKey->SetTecData(cData1);
        }
    }
}

void Post_WriteTecplotByOriginalGrid::BoundaryVTKVisualization(Grid *gridIn, DataContainer *cData, RDouble4D **qn, int nVarPlot)
{
    StructGrid *grid = StructGridCast(gridIn);
    int s_st[3] = { 0 }, s_ed[3] = { 0 };
    int         nk                 = grid->GetNK();
    RDouble3D   &xx                = *grid->GetStructX();
    RDouble3D   &yy                = *grid->GetStructY();
    RDouble3D   &zz                = *grid->GetStructZ();
    StructBCSet *compositeBCRegion = grid->GetStructBCSet();
    int         nBCRegion          = compositeBCRegion->GetnBCRegion();
    int TecioMission;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *BCRegion = compositeBCRegion->GetBCRegion(iBCRegion);
        int      bcType    = BCRegion->GetBCType();
        if (bcType == PHENGLEI::INTERFACE || bcType < 0) continue;
        TecioMission = TEC_SPACE::WriteBoundary;
        PHWrite(cData, &TecioMission, 1);
        PHWrite(cData, &bcType, 1);
        string bcName = BCRegion->GetBCName();
        cData->WriteString(bcName);
        int *s_lr3d = BCRegion->GetFaceDirectionIndex();
        int nDim    = GetDim();

        for (int m = 0; m < nDim; ++ m)
        {
            s_st[m] = BCRegion->GetStartPoint(m);
            s_ed[m] = BCRegion->GetEndPoint(m);
            if (BCRegion->GetFaceDirection() == m)
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

        Range I(ist, ied);
        Range J(jst, jed);
        Range K(kst, ked);

        if (nDim != THREE_D)
        {
            K.setRange(1, 1);
            kst = 1;
            ked = 1;
        }

        Int3D newIndex(I, J, K, fortranArray);
        int nSurf = BCRegion->GetFaceDirection() + 1;
        int count = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    newIndex(i, j, k) = count;
                    ++ count;
                }
            }
        }

        int nTotalNode = count;
        int nTotalCell = 1;
        for (int m = 0; m < 3; ++ m)
        {
            nTotalCell *= (BCRegion->GetEndPoint(m) - BCRegion->GetStartPoint(m) + 1);
        }

        PHWrite(cData, &nTotalNode, 1);
        PHWrite(cData, &nTotalCell, 1);

        RDouble *xBoundary  = new RDouble [nTotalNode];
        RDouble *yBoundary  = new RDouble [nTotalNode];
        RDouble *zBoundary  = new RDouble [nTotalNode];
        RDouble **qBoundary = new RDouble * [nVarPlot];
        for (int m = 0; m < nVarPlot; ++ m)
        {
            qBoundary[m] = new RDouble [nTotalNode];
        }

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
                    for (int m = 0; m < nVarPlot; ++ m)
                    {
                        qBoundary[m][iNode] = (*qn[m])(i, j, k, 0);
                    }
                    ++ iNode;
                }
            }
        }

        PHWrite(cData, xBoundary, nTotalNode);
        PHWrite(cData, yBoundary, nTotalNode);
        PHWrite(cData, zBoundary, nTotalNode);

        delete [] xBoundary;    xBoundary = nullptr;
        delete [] yBoundary;    yBoundary = nullptr;
        delete [] zBoundary;    zBoundary = nullptr;
        for (int m = 0; m < nVarPlot; ++ m)
        {
            PHWrite(cData, qBoundary[m], nTotalNode);
            delete [] qBoundary[m];    qBoundary[m] = nullptr;
        }
        delete [] qBoundary;    qBoundary = nullptr;

        int il1, jl1, kl1;
        GetNsurfIndex(nSurf, il1, jl1, kl1);
        il1 = 1 - il1;
        jl1 = 1 - jl1;
        kl1 = 1 - kl1;
        if (1 == nk) kl1 = 0;

        int *cell2Node = new int [nTotalCell * 4];
        int nCount = 0;
        for (int k = kst; k <= ked - kl1; ++ k)
        {
            for (int j = jst; j <= jed - jl1; ++ j)
            {
                for (int i = ist; i <= ied - il1; ++ i)
                {
                    int p1, p2, p3, p4;
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    if (1 == nSurf)
                    {
                        p1 = newIndex(i, j, k) + 1;
                        p2 = newIndex(i, jl, k) + 1;
                        p3 = newIndex(i, jl, kl) + 1;
                        p4 = newIndex(i, j, kl) + 1;
                    }
                    else if (2 == nSurf)
                    {
                        p1 = newIndex(i, j, k) + 1;
                        p2 = newIndex(il, j, k) + 1;
                        p3 = newIndex(il, j, kl) + 1;
                        p4 = newIndex(i, j, kl) + 1;
                    }
                    else
                    {
                        p1 = newIndex(i, j, k) + 1;
                        p2 = newIndex(il, j, k) + 1;
                        p3 = newIndex(il, jl, k) + 1;
                        p4 = newIndex(i, jl, k) + 1;
                    }

                    cell2Node[nCount ++] = p1;
                    cell2Node[nCount ++] = p2;
                    cell2Node[nCount ++] = p3;
                    cell2Node[nCount ++] = p4;
                }
            }
        }
        PHWrite(cData, cell2Node, nTotalCell * 4);
        delete [] cell2Node;    cell2Node = nullptr;
    }
}

void Post_WriteTecplotByOriginalGrid::FieldVisualization(Grid *gridIn, DataContainer *cData, RDouble4D **qn, int nVarPlot)
{
    StructGrid *grid = StructGridCast(gridIn);
    int        ni    = grid->GetNI();
    int        nj    = grid->GetNJ();
    int        nk    = grid->GetNK();
    RDouble3D  &xx   = *grid->GetStructX();
    RDouble3D  &yy   = *grid->GetStructY();
    RDouble3D  &zz   = *grid->GetStructZ();
    int        nl    = 5;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);
    int nChem = 0;
    GlobalDataBase::GetData("nchem", &nChem, PHINT, 1);
    RDouble refGama = 1.4;
    GlobalDataBase::GetData("refGama", &refGama, PHDOUBLE, 1);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalNode = grid->GetNTotalNode();

    if (1 == nk)
    {
        int TecioMission = TEC_SPACE::WriteBlock;
        PHWrite(cData, &TecioMission, 1);
        PHWrite(cData, &nTotalNode, 1);
        PHWrite(cData, &nTotalCell, 1);

        RDouble *xBoundary  = new RDouble [nTotalNode];
        RDouble *yBoundary  = new RDouble [nTotalNode];
        RDouble *zBoundary  = new RDouble [nTotalNode];
        RDouble **qBoundary = new RDouble * [nVarPlot];
        for (int m = 0; m < nVarPlot; ++ m)
        {
            qBoundary[m] = new RDouble [nTotalNode];
        }

        int ist, ied, jst, jed, kst, ked;
        grid->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);
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
                    for (int m = 0; m < nVarPlot; ++ m)
                    {
                        qBoundary[m][iNode] = (*qn[m])(i, j, k, 0);
                    }
                    iNode ++;
                }
            }
        }

        PHWrite(cData, xBoundary, nTotalNode);
        PHWrite(cData, yBoundary, nTotalNode);
        PHWrite(cData, zBoundary, nTotalNode);

        delete [] xBoundary;    xBoundary = nullptr;
        delete [] yBoundary;    yBoundary = nullptr;
        delete [] zBoundary;    zBoundary = nullptr;
        for (int m = 0; m < nVarPlot; ++ m)
        {
            PHWrite(cData, qBoundary[m], nTotalNode);
            delete [] qBoundary[m];    qBoundary[m] = nullptr;
        }
        delete [] qBoundary;    qBoundary = nullptr;

        int *cell2Node = new int [nTotalCell * 4];
        int nCount = 0;
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = 1; i <= ni - 1; ++ i)
            {
                cell2Node[nCount ++] = i + (j - 1) * ni;
                cell2Node[nCount ++] = i + 1 + (j - 1) * ni;
                cell2Node[nCount ++] = i + 1 + (j)*ni;
                cell2Node[nCount ++] = i + (j)*ni;
            }
        }
        PHWrite(cData, cell2Node, nTotalCell * 4);
        delete [] cell2Node;    cell2Node = nullptr;
    }
}

void Post_WriteTecplotByOriginalGrid::BoundaryVisualization(Grid *gridIn, DataContainer *cData, RDouble4D **qn, int nVarPlot)
{
    StructGrid *grid = StructGridCast(gridIn);
    int s_st[3] = { 0 }, s_ed[3] = { 0 };
    RDouble3D   &xx                = *grid->GetStructX();
    RDouble3D   &yy                = *grid->GetStructY();
    RDouble3D   &zz                = *grid->GetStructZ();
    StructBCSet *compositeBCRegion = grid->GetStructBCSet();
    int         nBCRegion          = compositeBCRegion->GetnBCRegion();
    int TecioMission;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcRegion = compositeBCRegion->GetBCRegion(iBCRegion);
        int      bcType    = bcRegion->GetBCType();
        if (bcType == PHENGLEI::INTERFACE || bcType < 0) continue;
        TecioMission = TEC_SPACE::WriteBoundary;
        PHWrite(cData, &TecioMission, 1);
        PHWrite(cData, &bcType, 1);
        string bcName = bcRegion->GetBCName();
        cData->WriteString(bcName);
        int *s_lr3d = bcRegion->GetFaceDirectionIndex();
        int nDim    = GetDim();

        for (int m = 0; m < nDim; ++ m)
        {
            s_st[m] = bcRegion->GetStartPoint(m);
            s_ed[m] = bcRegion->GetEndPoint(m);

            if (bcRegion->GetFaceDirection() == m)
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

        Range I(ist, ied);
        Range J(jst, jed);
        Range K(kst, ked);

        if (nDim != THREE_D)
        {
            K.setRange(1, 1);
            kst = 1;
            ked = 1;
        }

        Int3D newIndex(I, J, K, fortranArray);
        int Count = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    newIndex(i, j, k) = Count;
                    ++ Count;
                }
            }
        }

        int ii = ied - ist + 1;
        int jj = jed - jst + 1;
        int kk = ked - kst + 1;

        PHWrite(cData, &ii, 1);
        PHWrite(cData, &jj, 1);
        PHWrite(cData, &kk, 1);
        int nTotalNode = Count;
        PHWrite(cData, &nTotalNode, 1);

        RDouble *x = new RDouble [nTotalNode];
        RDouble *y = new RDouble [nTotalNode];
        RDouble *z = new RDouble [nTotalNode];

        Count = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    x[Count] = xx(i, j, k);
                    y[Count] = yy(i, j, k);
                    z[Count] = zz(i, j, k);
                    Count ++;
                }
            }
        }

        PHWrite(cData, x, nTotalNode);
        PHWrite(cData, y, nTotalNode);
        PHWrite(cData, z, nTotalNode);

        delete [] x;    x = nullptr;
        delete [] y;    y = nullptr;
        delete [] z;    z = nullptr;

        Count = 0;
        for (int n = 0; n < nVarPlot; ++ n)
        {
            RDouble *qq = new RDouble [nTotalNode];
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        qq[Count] = (*qn[n])(i, j, k, 0);
                        Count ++;
                    }
                }
            }
            PHWrite(cData, qq, nTotalNode);
            delete [] qq;    qq = nullptr;
            Count = 0;
        }
    }
}

void Post_WriteTecplotByOriginalGrid::SaveDataForTecio(Grid *gridIn, DataContainer *cData, RDouble4D **qn, int nl)
{
    StructGrid *grid        = StructGridCast(gridIn);
    int        TecioMission = TEC_SPACE::WriteBlock;
    PHWrite(cData, &TecioMission, 1);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    PHWrite(cData, &ni, 1);
    PHWrite(cData, &nj, 1);
    PHWrite(cData, &nk, 1);
    int nTotalNode = grid->GetNTotalNode();
    PHWrite(cData, &nTotalNode, 1);

    RDouble *x = new RDouble [nTotalNode];
    RDouble *y = new RDouble [nTotalNode];
    RDouble *z = new RDouble [nTotalNode];

    RDouble3D &xx = *grid->GetStructX();
    RDouble3D &yy = *grid->GetStructY();
    RDouble3D &zz = *grid->GetStructZ();

    int Count = 0;
    for (int k = 1; k <= nk; ++ k)
    {
        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                x[Count] = xx(i, j, k);
                y[Count] = yy(i, j, k);
                z[Count] = zz(i, j, k);
                Count ++;
            }
        }
    }

    PHWrite(cData, x, nTotalNode);
    PHWrite(cData, y, nTotalNode);
    PHWrite(cData, z, nTotalNode);

    delete [] x;    x = nullptr;
    delete [] y;    y = nullptr;
    delete [] z;    z = nullptr;

    Count = 0;
    for (int m = 0; m < nl; ++ m)
    {
        RDouble *qq = new RDouble [nTotalNode];
        for (int k = 1; k <= nk; ++ k)
        {
            for (int j = 1; j <= nj; ++ j)
            {
                for (int i = 1; i <= ni; ++ i)
                {
                    qq[Count] = (*qn[m])(i, j, k, 0);
                    Count ++;
                }
            }
        }
        PHWrite(cData, qq, nTotalNode);
        delete [] qq;    qq = nullptr;
        Count = 0;
    }
}

void Post_WriteTecplotByOriginalGrid::BoundaryVTKVisualization(Grid *gridIn, DataContainer *cData, RDouble **qn, int nVarPlot)
{
    UnstructGrid *grid                 = UnstructGridCast(gridIn);
    int          nTotalNode            = grid->GetNTotalNode();
    int          nBoundFace            = grid->GetNBoundFace();
    int          *nodeNumberofEachFace = grid->GetNodeNumberOfEachFace();
    int          *face2Node            = grid->GetFace2Node();
    RDouble      *x                    = grid->GetX();
    RDouble      *y                    = grid->GetY();
    RDouble      *z                    = grid->GetZ();
    RDouble      *xFC                  = grid->GetFaceCenterX();
    RDouble      *yFC                  = grid->GetFaceCenterY();
    RDouble      *zFC                  = grid->GetFaceCenterZ();
    RDouble      **qFaceCenter         = NewPointer2<RDouble>(nVarPlot, nBoundFace);

    using namespace PHENGLEI;
    vector<vector<int> > face2NodeList(nBoundFace);
    vector<int> linkMap;
    set<pair<int, string> > bcNameMap;
    UnstructBCSet *unstructBCSet      = grid->GetUnstructBCSet();
    int           *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    int           nodeCount           = 0;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        string     bcName    = bcRegion->GetBCName();
        int        bcType    = bcRegion->GetBCType();
        bcNameMap.insert(pair<int, string>(bcType, bcName));

        for (int m = 0; m < nVarPlot; ++ m)
        {
            RDouble value = 0.0;
            for (int iNode = 0; iNode < nodeNumberofEachFace[iFace]; ++ iNode)
            {
                int nodeIndex = face2Node[nodeCount + iNode];
                value += qn[m][nodeIndex];
            }
            qFaceCenter[m][iFace] = value / nodeNumberofEachFace[iFace];
        }
        nodeCount += nodeNumberofEachFace[iFace];
    }

    int TecioMission;
    set<pair<int, string> >::iterator iter;
    for (iter = bcNameMap.begin(); iter != bcNameMap.end(); ++ iter)
    {
        if ((*iter).first == PHENGLEI::INTERFACE || (*iter).first < 0) continue;
        TecioMission = TEC_SPACE::WriteBoundary;
        PHWrite(cData, &TecioMission, 1);
        int BCtype = (*iter).first;
        PHWrite(cData, &BCtype, 1);
        string bcName = (*iter).second;
        cData->WriteString(bcName);
        face2NodeList.resize(0);
        linkMap.resize(0);
        GetFace2NodeList(grid, *iter, linkMap, face2NodeList);
        int NumPts      = static_cast<int>(linkMap.size());
        int NumElements = static_cast<int>(face2NodeList.size());

        PHWrite(cData, &NumPts, 1);
        PHWrite(cData, &NumElements, 1);

        RDouble *xBoundary = new RDouble [NumPts];
        RDouble *yBoundary = new RDouble [NumPts];
        RDouble *zBoundary = new RDouble [NumPts];

        for (int iNode = 0; iNode < NumPts; ++ iNode)
        {
            int nodeIndex = linkMap[iNode];
            if (nodeIndex < nTotalNode)
            {
                xBoundary[iNode] = x[nodeIndex];
                yBoundary[iNode] = y[nodeIndex];
                zBoundary[iNode] = z[nodeIndex];
            }
            else
            {
                int faceIndex    = nodeIndex - nTotalNode;
                xBoundary[iNode] = xFC[faceIndex];
                yBoundary[iNode] = yFC[faceIndex];
                zBoundary[iNode] = zFC[faceIndex];
            }
        }

        PHWrite(cData, xBoundary, NumPts);
        PHWrite(cData, yBoundary, NumPts);
        PHWrite(cData, zBoundary, NumPts);

        delete [] xBoundary;    xBoundary = nullptr;
        delete [] yBoundary;    yBoundary = nullptr;
        delete [] zBoundary;    zBoundary = nullptr;

        for (int m = 0; m < nVarPlot; ++ m)
        {
            RDouble *qBoundary = new RDouble [NumPts];
            for (int iNode = 0; iNode < NumPts; ++ iNode)
            {
                int nodeIndex = linkMap[iNode];
                if (nodeIndex < nTotalNode)
                {
                    qBoundary[iNode] = qn[m][nodeIndex];
                }
                else
                {
                    int faceIndex    = nodeIndex - nTotalNode;
                    qBoundary[iNode] = qFaceCenter[m][faceIndex];
                }
            }
            PHWrite(cData, qBoundary, NumPts);
            delete [] qBoundary;    qBoundary = nullptr;
        }

        int *cell2Node = new int [NumElements * 4];
        int nCount = 0;
        for (int iCell = 0; iCell < NumElements; ++ iCell)
        {
            uint_t np     = face2NodeList[iCell].size();
            int    index0 = face2NodeList[iCell][0];
            for (int ip = 0; ip < np; ++ ip)
            {
                cell2Node[nCount] = face2NodeList[iCell][ip] + 1;
                nCount ++;
            }
            for (uint_t ip = np; ip < 4; ++ ip)
            {
                cell2Node[nCount] = index0 + 1;
                nCount ++;
            }
        }
        PHWrite(cData, cell2Node, NumElements * 4);
        delete [] cell2Node;    cell2Node = nullptr;
    }
    DelPointer2(qFaceCenter);
}

void Post_WriteTecplotByOriginalGrid::FieldVisualization(Grid *gridIn, DataContainer *cData, RDouble **qn, int nVarPlot)
{
    if (PHSPACE::TWO_D == GetDim())
    {
        FieldVisualizationForVTK(gridIn, cData, qn, nVarPlot);
        return;
    }

    int TecioMission = TEC_SPACE::WriteBlock;
    PHWrite(cData, &TecioMission, 1);
    UnstructGrid *grid                 = UnstructGridCast(gridIn);
    int          nTotalNode            = grid->GetNTotalNode();
    int          nTotalCell            = grid->GetNTotalCell();
    int          nTotalFace            = grid->GetNTotalFace();
    RDouble      *x                    = grid->GetX();
    RDouble      *y                    = grid->GetY();
    RDouble      *z                    = grid->GetZ();
    int          *nodeNumberOfEachCell = grid->GetNodeNumberOfEachCell();
    int          *cell2Node            = grid->GetCell2Node();

    PHWrite(cData, &nTotalNode, 1);
    PHWrite(cData, &nTotalCell, 1);
    PHWrite(cData, x, nTotalNode);
    PHWrite(cData, y, nTotalNode);
    PHWrite(cData, z, nTotalNode);

    for (int m = 0; m < nVarPlot; ++ m)
    {
        RDouble *qq = qn[m];
        PHWrite(cData, qq, nTotalNode);
    }
    PHWrite(cData, &nTotalFace, 1);

    int Count = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        Count += nodeNumberOfEachCell[iCell];
    }
    PHWrite(cData, nodeNumberOfEachCell, nTotalCell);
    PHWrite(cData, cell2Node, Count);

    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();
    int **cell2Face           = grid->GetCell2Face();
    PHWrite(cData, faceNumberOfEachCell, nTotalCell);
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        PHWrite(cData, cell2Face[iCell], faceNumberOfEachCell[iCell]);
    }

    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
    int *face2Node            = grid->GetFace2Node();
    Count = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        Count += nodeNumberOfEachFace[iFace];
    }
    PHWrite(cData, nodeNumberOfEachFace, nTotalFace);
    PHWrite(cData, face2Node, Count);
}

void Post_WriteTecplotByOriginalGrid::FieldVisualizationForVTK(Grid *gridIn, DataContainer *cData, RDouble **qn, int nVarPlot)
{
    int TecioMission = TEC_SPACE::WriteBlock;
    PHWrite(cData, &TecioMission, 1);
    UnstructGrid *grid            = UnstructGridCast(gridIn);
    int          nTotalNode       = grid->GetNTotalNode();
    int          nTotalCell       = grid->GetNTotalCell();
    int          nBoundFace       = grid->GetNBoundFace();
    int          nTotalFace       = grid->GetNTotalFace();
    int          *face2Node       = grid->GetFace2Node();
    int          *leftCellofFace  = grid->GetLeftCellOfFace();
    int          *rightCellofFace = grid->GetRightCellOfFace();
    RDouble      *x               = grid->GetX();
    RDouble      *y               = grid->GetY();
    RDouble      *z               = grid->GetZ();

    //! Find out cellToFace.
    HyList<int> *cell2Face = new HyList <int> [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2Face[iCell].SetAverNnode(4);
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        cell2Face[le].insert(iFace);
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];
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
        int targetNode  = ip2;
        int checkNode   = ip1;
        int nFaceOfCell = cell2Face[iCell].size();
        int isClockWise = false;

        while (!isClockWise)
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

    PHWrite(cData, &nTotalNode, 1);
    PHWrite(cData, &nTotalCell, 1);
    PHWrite(cData, x, nTotalNode);
    PHWrite(cData, y, nTotalNode);
    PHWrite(cData, z, nTotalNode);
    for (int m = 0; m < nVarPlot; ++ m)
    {
        RDouble *qq = qn[m];
        PHWrite(cData, qq, nTotalNode);
    }

    int *cell2NodeVtk = new int [nTotalCell * 4];
    int nCount = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (cell2Node[iCell].size() == 3)
        {
            for (int m = 0; m < 3; ++ m)
            {
                cell2NodeVtk[nCount] = cell2Node[iCell].GetData(m) + 1;
                nCount ++;
            }
            cell2NodeVtk[nCount] = cell2Node[iCell].GetData(2) + 1;
            nCount ++;
        }
        else if (cell2Node[iCell].size() == 4)
        {
            for (int m = 0; m < 4; ++ m)
            {
                cell2NodeVtk[nCount] = cell2Node[iCell].GetData(m) + 1;
                nCount ++;
            }
        }
        else
        {
            cout << "Error: this function only support triangle and quad !\n";
            exit(0);
        }
    }
    PHWrite(cData, cell2NodeVtk, nTotalCell * 4);

    delete [] cell2NodeVtk;    cell2NodeVtk = nullptr;
    delete [] cell2Node;       cell2Node    = nullptr;
    delete [] cell2Face;       cell2Face    = nullptr;
}

void Post_WriteTecplotByOriginalGrid::BoundaryVisualization(Grid *gridIn, DataContainer *cData, RDouble **qn, int nVarPlot)
{
    UnstructGrid *grid                 = UnstructGridCast(gridIn);
    int          nTotalNode            = grid->GetNTotalNode();
    int          nBoundFace            = grid->GetNBoundFace();
    int          *nodeNumberofEachFace = grid->GetNodeNumberOfEachFace();
    int          *face2Node            = grid->GetFace2Node();
    RDouble      *x                    = grid->GetX();
    RDouble      *y                    = grid->GetY();
    RDouble      *z                    = grid->GetZ();
    RDouble      *xFC                  = grid->GetFaceCenterX();
    RDouble      *yFC                  = grid->GetFaceCenterY();
    RDouble      *zFC                  = grid->GetFaceCenterZ();
    RDouble      **qFaceCenter         = NewPointer2<RDouble>(nVarPlot, nBoundFace);

    using namespace PHENGLEI;
    vector<vector<int> > face2NodeList(nBoundFace);
    vector<int> linkMap;
    set<pair<int, string> > bcNameMap;
    UnstructBCSet *unstructBCSet      = grid->GetUnstructBCSet();
    int           *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    int           nodeCount           = 0;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        string     bcName    = bcRegion->GetBCName();
        int        bcType    = bcRegion->GetBCType();
        bcNameMap.insert(pair<int, string>(bcType, bcName));

        for (int m = 0; m < nVarPlot; ++ m)
        {
            RDouble value = 0.0;
            for (int iNode = 0; iNode < nodeNumberofEachFace[iFace]; ++ iNode)
            {
                int nodeIndex = face2Node[nodeCount + iNode];
                value         += qn[m][nodeIndex];
            }
            qFaceCenter[m][iFace] = value / nodeNumberofEachFace[iFace];
        }
        nodeCount += nodeNumberofEachFace[iFace];
    }

    int TecioMission;
    set<pair<int, string> >::iterator iter;
    for (iter = bcNameMap.begin(); iter != bcNameMap.end(); ++ iter)
    {
        if ((*iter).first == PHENGLEI::INTERFACE || (*iter).first < 0) continue;
        face2NodeList.resize(0);
        linkMap.resize(0);
        GetFace2NodeList(grid, *iter, linkMap, face2NodeList);
        if (0 == face2NodeList.size())
        {
            continue;
        }

        TecioMission = TEC_SPACE::WriteBoundary;
        PHWrite(cData, &TecioMission, 1);
        int BCtype = (*iter).first;
        PHWrite(cData, &BCtype, 1);
        string bcName = (*iter).second;
        cData->WriteString(bcName);

        int NumPts      = static_cast<int>(linkMap.size());
        int NumElements = static_cast<int>(face2NodeList.size());
        int NumFaces    = NumElements * 4;

        PHWrite(cData, &NumPts, 1);
        PHWrite(cData, &NumElements, 1);
        PHWrite(cData, &NumFaces, 1);

        int TotalNumFaceNodes_Rect = NumElements * 4;
        int *cell2Node             = new int [TotalNumFaceNodes_Rect];
        int *nodeNumberofEachCell  = new int [NumElements];

        int count = 0;
        for (std::size_t i = 0; i < face2NodeList.size(); ++ i)
        {
            nodeNumberofEachCell[i] = 4;
            uint_t np                   = face2NodeList[i].size();
            int    index0               = face2NodeList[i][0];
            for (int ip = 0; ip < np; ++ ip)
            {
                cell2Node[count] = face2NodeList[i][ip];
                count ++;
            }
            for (uint_t ip = np; ip < 4; ++ ip)
            {
                cell2Node[count] = index0;
                count ++;
            }
        }

        PHWrite(cData, &TotalNumFaceNodes_Rect, 1);
        PHWrite(cData, nodeNumberofEachCell, NumElements);
        PHWrite(cData, &NumPts, 1);

        RDouble *xx = new RDouble [NumPts];
        RDouble *yy = new RDouble [NumPts];
        RDouble *zz = new RDouble [NumPts];

        for (std::size_t iNode = 0; iNode < linkMap.size(); ++ iNode)
        {
            int nodeIndex = linkMap[iNode];
            if (nodeIndex < nTotalNode)
            {
                xx[iNode] = x[nodeIndex];
                yy[iNode] = y[nodeIndex];
                zz[iNode] = z[nodeIndex];
            }
            else
            {
                int faceIndex = nodeIndex - nTotalNode;
                xx[iNode] = xFC[faceIndex];
                yy[iNode] = yFC[faceIndex];
                zz[iNode] = zFC[faceIndex];
            }
        }

        PHWrite(cData, xx, NumPts);
        PHWrite(cData, yy, NumPts);
        PHWrite(cData, zz, NumPts);

        for (int m = 0; m < nVarPlot; ++ m)
        {
            RDouble *qBoundary = new RDouble [NumPts];

            for (int iNode = 0; iNode < NumPts; ++ iNode)
            {
                int nodeIndex = linkMap[iNode];

                if (nodeIndex < nTotalNode)
                {
                    qBoundary[iNode] = qn[m][nodeIndex];
                }
                else
                {
                    int faceIndex = nodeIndex - nTotalNode;
                    qBoundary[iNode] = qFaceCenter[m][faceIndex];
                }
            }
            PHWrite(cData, qBoundary, NumPts);
            delete [] qBoundary;    qBoundary = nullptr;
        }
        PHWrite(cData, cell2Node, TotalNumFaceNodes_Rect);

        delete [] cell2Node;               cell2Node            = nullptr;
        delete [] nodeNumberofEachCell;    nodeNumberofEachCell = nullptr;
        delete [] xx;                      xx                   = nullptr;
        delete [] yy;                      yy                   = nullptr;
        delete [] zz;                      zz                   = nullptr;
    }
    DelPointer2(qFaceCenter);
}

void Post_WriteTecplotByOriginalGrid::SaveDataForTecio(Grid *gridIn, DataContainer *cData, RDouble **qn, int nVarPlot)
{
    UnstructGrid *grid               = UnstructGridCast(gridIn);
    int          numberOfNodes       = grid->GetNTotalNode();
    int          numberOfCells       = grid->GetNTotalCell();
    int          numberOfFaces       = grid->GetNTotalFace();
    int          numberOfNodesInZone = 0;
    int          numberOfFacesInZone = 0;
    int          numberOfCellsInZone = 0;

    vector<int> keyNodes;
    vector<int> keyFaces;
    vector<int> nodeList;
    vector<int> cellList;

    keyNodes.resize(numberOfNodes);
    keyFaces.resize(numberOfFaces);
    nodeList.resize(numberOfNodes);
    cellList.resize(numberOfCells);

    SetField(keyNodes, 0, numberOfNodes);
    SetField(keyFaces, 0, numberOfFaces);
    SetField(nodeList, -1, numberOfNodes);
    SetField(cellList, -1, numberOfCells);

    int *keyActiveOfCells        = grid->GetBlankIndex();
    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();
    int **cellFaceIndexContainer = grid->GetCell2Face();
    int *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();
    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer  = grid->GetFace2Node();
    int *leftCellIndexContainer  = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        if (1 == keyActiveOfCells[iCell])
        {
            cellList[iCell] = numberOfCellsInZone;
            for (int iNodes = 0; iNodes < cellNodeNumberContainer[iCell]; ++ iNodes)
            {
                keyNodes[cellNodeIndexContainer[iCell][iNodes]] = 1;
            }
            for (int iFaces = 0; iFaces < cellFaceNumberContainer[iCell]; ++ iFaces)
            {
                keyFaces[cellFaceIndexContainer[iCell][iFaces]] = 1;
            }
            numberOfCellsInZone ++;
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (1 == keyNodes[iNode])
        {
            nodeList[iNode] = numberOfNodesInZone;
            numberOfNodesInZone ++;
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (1 == keyFaces[iFace])
        {
            numberOfFacesInZone ++;
        }
    }

    int totalNumFaceNodes = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (keyFaces[iFace] != 1) continue;
        totalNumFaceNodes += faceNodeNumberContainer[iFace];
    }

    if (numberOfNodesInZone == 0)
    {
        return;
    }

    int TecioMission = TEC_SPACE::WriteBlock;
    PHWrite(cData, TecioMission);
    PHWrite(cData, numberOfNodesInZone);
    PHWrite(cData, numberOfCellsInZone);
    PHWrite(cData, numberOfFacesInZone);
    PHWrite(cData, totalNumFaceNodes);

    int count = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (keyFaces[iFace] != 1) continue;

        int faceNodeNumber = faceNodeNumberContainer[iFace];
        PHWrite(cData, faceNodeNumber);
    }
    PHWrite(cData, numberOfNodesInZone);

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        PHWrite(cData, x[iNode]);
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        PHWrite(cData, y[iNode]);
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        PHWrite(cData, z[iNode]);
    }

    for (int iVarible = 0; iVarible < nVarPlot; ++ iVarible)
    {
        RDouble *qq = qn[iVarible];
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if (keyNodes[iNode] != 1) continue;
            PHWrite(cData, qq[iNode]);
        }
    }

    count = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];
        if (keyFaces[iFace] != 1)
        {
            count += faceNodeNumber;
            continue;
        }
        for (int iNode = 0; iNode < faceNodeNumberContainer[iFace]; ++ iNode)
        {
            int nodeIndex = nodeList[faceNodeIndexContainer[count++]];
            PHWrite(cData, nodeIndex);
        }
    }

    //! Left cell.
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (keyFaces[iFace] != 1) continue;
        int cellIndex = cellList[leftCellIndexContainer[iFace]];
        PHWrite(cData, cellIndex);
    }

    //! Right cell.
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int cellIndex;
        if (0 == keyFaces[iFace]) continue;
        int rightCell = rightCellIndexContainer[iFace];
        if (rightCell >= numberOfCells || rightCell < 0)
        {
            cellIndex = -1;
        }
        else
        {
            cellIndex = cellList[rightCell];
        }
        PHWrite(cData, cellIndex);
    }
}

void Post_WriteTecplotByOriginalGrid::PostWrite(ActionKey *actKey)
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if (true == actKey->GetVTKVisual() && vtkDataList.size() > 0)
    {
        bool CharacteristicBoundary = false;
        DumpToVTK DumpToVTK(vtkDataList, actKey->GetFileName(), CharacteristicBoundary, 0);
        DumpToVTK.Run();
        WriteSentinelFile();
    }

    if (dataList.size() > 0)
    {
        int visualFileType = GlobalDataBase::GetIntParaFromDB("visualfileType");
        if (0 == visualFileType)
        {
#ifdef USE_TecplotLib
            DumpToTecio dumptotecio(dataList, actKey->tecfilename);
            dumptotecio.Run();
#else
            ostringstream WarningInfo;
            WarningInfo << "Warning : visualfileType = 0 is not right, dump flow data into ASCII file! \n";
            PrintToWindow(WarningInfo);
            DumpToTecplotASCII dumpToTecplotASCII(dataList, actKey->tecfilename);
            dumpToTecplotASCII.Run();
#endif
        }
        else if (1 == visualFileType)
        {
            DumpToTecplotASCII dumpToTecplotASCII(dataList, actKey->tecfilename);
            dumpToTecplotASCII.Run();
        }
        else
        {
            TK_Exit::UnexpectedVarValue("visualfileType = ", visualFileType);
        }
    }

    for (unsigned int iZone = 0; iZone < vtkDataList.size(); ++ iZone)
    {
        delete vtkDataList[iZone];
    }
    vtkDataList.resize(0);

    for (unsigned int iZone = 0; iZone < dataList.size(); ++ iZone)
    {
        delete dataList[iZone];
    }
    dataList.resize(0);
}

void Post_WriteTecplotByOriginalGrid::ComputerNodeDataOnPartitionGrid(int iZone, ActionKey *actKey)
{
    Grid *grid            = GetGrid(iZone, actKey->level);
    int  nVisualVariables = GlobalDataBase::GetIntParaFromDB("nVisualVariables");
    int visualVariables[100];
    GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);
    Post_Visual *postVisualization = new Post_Visual(nVisualVariables, visualVariables, grid);
    int Nsolver = GlobalSolvers::GetNumberOfSolvers(iZone);
    for (int iSolver = 0; iSolver < Nsolver; ++ iSolver)
    {
        GlobalSolvers::GetSolver(iZone, iSolver)->ComputePostVisualVariables(postVisualization);
    }

    int nVarPlot     = postVisualization->GetNumberofVisualVariables();
    nVisualVariables = nVarPlot;
    int nChem        = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nChem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nVisualVariables = nVisualVariables + numberOfSpecies;
    }

    int gridType = grid->Type();
    if (gridType == UNSTRUCTGRID)
    {
        RDouble **qn = new RDouble * [nVisualVariables];
        postVisualization->GetAllVisualNodeVarPtr(qn);
        UpdataNodeData(actKey, grid, qn, nVisualVariables);
        delete [] qn;    qn = nullptr;
    }
    else
    {
        RDouble4D **qn = new RDouble4D * [nVisualVariables];
        postVisualization->GetAllVisualNodeVarPtr(qn);
        UpdataNodeData(actKey, grid, qn, nVisualVariables);
        delete [] qn;    qn = nullptr;
    }
    FreePointer(postVisualization);
}

void Post_WriteTecplotByOriginalGrid::CollectNodeData(ActionKey *actKey)
{
    DataContainer *cData = actKey->GetData();
    cData->MoveToBegin();
    int originalGridIndex = -1;
    PHRead(cData, &originalGridIndex, 1);
    Post_Visual *postVisualization = originalGridData[originalGridIndex];
    Grid        *grid              = OriginalGrid[originalGridIndex];
    int         gridType           = grid->Type();
    if (UNSTRUCTGRID == gridType)
    {
        CollectUnstructedData(actKey, postVisualization);
    }
    else
    {
        CollectStructedData(actKey, postVisualization);
    }
}

void Post_WriteTecplotByOriginalGrid::CollectUnstructedData(ActionKey *actKey, Post_Visual *postVisualization)
{
    UnstructGrid  *grid  = UnstructGridCast(postVisualization->GetGrid());
    DataContainer *cData = actKey->GetData();
    int nTotalNodeOfParGrid;
    PHRead(cData, nTotalNodeOfParGrid);
    int *nodeIndex = new int [nTotalNodeOfParGrid];
    PHRead(cData, nodeIndex, nTotalNodeOfParGrid);

    int nVarPlot         = postVisualization->GetNumberofVisualVariables();
    int nVisualVariables = nVarPlot;
    int nChem            = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nChem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nVisualVariables    = nVisualVariables + numberOfSpecies;
    }

    RDouble **qn = new RDouble * [nVisualVariables];
    postVisualization->GetAllVisualNodeVarPtr(qn);
    for (int iVar = 0; iVar < nVisualVariables; ++ iVar)
    {
        RDouble *qq = new RDouble [nTotalNodeOfParGrid];
        PHRead(cData, qq, nTotalNodeOfParGrid);
        for (int iNode = 0; iNode < nTotalNodeOfParGrid; ++ iNode)
        {
            int currentNodeIndex       = nodeIndex[iNode];
            qn[iVar][currentNodeIndex] = qq[iNode];
        }
        delete [] qq;    qq = nullptr;
    }

    delete [] qn;           qn        = nullptr;
    delete [] nodeIndex;    nodeIndex = nullptr;
}

void Post_WriteTecplotByOriginalGrid::CollectStructedData(ActionKey *actKey, Post_Visual *postVisualization)
{
    StructGrid    *grid  = StructGridCast(postVisualization->GetGrid());
    DataContainer *cData = actKey->GetData();
    int ordinaryDimStartIndex[THREE_D];
    int ordinaryDimEndIndex[THREE_D];
    PHRead(cData, ordinaryDimStartIndex, THREE_D);
    PHRead(cData, ordinaryDimEndIndex, THREE_D);

    int nVarPlot         = postVisualization->GetNumberofVisualVariables();
    int nVisualVariables = nVarPlot;
    int nChem            = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nChem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nVisualVariables    = nVisualVariables + numberOfSpecies;
    }

    RDouble4D **qn = new RDouble4D * [nVisualVariables];
    postVisualization->GetAllVisualNodeVarPtr(qn);
    for (int iVar = 0; iVar < nVisualVariables; ++ iVar)
    {
        for (int kIndex = ordinaryDimStartIndex[2]; kIndex <= ordinaryDimEndIndex[2]; ++ kIndex)
        {
            for (int jIndex = ordinaryDimStartIndex[1]; jIndex <= ordinaryDimEndIndex[1]; ++ jIndex)
            {
                for (int iIndex = ordinaryDimStartIndex[0]; iIndex <= ordinaryDimEndIndex[0]; ++ iIndex)
                {
                    RDouble qVar;
                    PHRead(cData, &qVar, 1);
                    (*qn[iVar])(iIndex, jIndex, kIndex, 0) = qVar;
                }
            }
        }
    }
    delete [] qn;    qn = nullptr;
}

void Post_WriteTecplotByOriginalGrid::UpdataNodeData(ActionKey *actKey, Grid *grid, RDouble **qn, int nVisualVariables)
{
    DataContainer *cData     = actKey->GetData();
    UnstructGrid  *unsGrid   = UnstructGridCast(grid);
    int           nTotalNode = unsGrid->GetNTotalNode();
    PHWrite(cData, nTotalNode);
    int *ordinaryNodeIndex = unsGrid->GetOrdinaryNodeIndex();
    PHWrite(cData, ordinaryNodeIndex, nTotalNode);
    for (int iVar = 0; iVar < nVisualVariables; ++ iVar)
    {
        PHWrite(cData, qn[iVar], nTotalNode);
    }
}

void Post_WriteTecplotByOriginalGrid::UpdataNodeData(ActionKey *actKey, Grid *grid, RDouble4D **qn, int nVisualVariables)
{
    DataContainer *cData                  = actKey->GetData();
    StructGrid    *strGrid                = StructGridCast(grid);
    int           ni                      = strGrid->GetNI();
    int           nj                      = strGrid->GetNJ();
    int           nk                      = strGrid->GetNK();
    int           * ordinaryDimStartIndex = strGrid->GetOrdinaryDimStartIndex();
    int           * ordinaryDimEndIndex   = strGrid->GetOrdinaryDimEndIndex();

    PHWrite(cData, ordinaryDimStartIndex, THREE_D);
    PHWrite(cData, ordinaryDimEndIndex, THREE_D);

    for (int iVar = 0; iVar < nVisualVariables; ++ iVar)
    {
        for (int kIndex = 1; kIndex <= nk; ++ kIndex)
        {
            for (int jIndex = 1; jIndex <= nj; ++ jIndex)
            {
                for (int iIndex = 1; iIndex <= ni; ++ iIndex)
                {
                    RDouble qVar = (*qn[iVar])(iIndex, jIndex, kIndex, 0);
                    PHWrite(cData, &qVar, 1);
                }
            }
        }
    }
}

void Post_WriteTecplotByOriginalGrid::SetOriginalGridProc(int nOriginalGridsIn, int *originalGridProcOut)
{
    int numberOfProcessor = PHMPI::GetNumberOfProcessor();
    this->nOriginalGrids   = nOriginalGridsIn;
    this->originalGridProc = new int[nOriginalGridsIn];
    for (int iZone = 0; iZone < nOriginalGridsIn; ++ iZone)
    {
        originalGridProcOut[iZone] = iZone % numberOfProcessor;
        originalGridProc[iZone]    = iZone % numberOfProcessor;
    }

    int nZones                  = PHMPI::GetNumberofGlobalZones();
    originalGridIndexOfEachGrid = new int[nZones];
    SetField(originalGridIndexOfEachGrid, 0, nZones);
    int *tempIndex = new int[nZones];
    SetField(tempIndex, 0, nZones);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = GetGrid(iZone);
        if (!grid)
        {
            continue;
        }
        int originalGridIndex = grid->GetOrdinaryGridIndex();
        tempIndex[iZone]      = originalGridIndex;
    }
    PH_AllReduce(tempIndex, originalGridIndexOfEachGrid, nZones, PH_SUM);
    delete [] tempIndex;    tempIndex = nullptr;
}

int Post_WriteTecplotByOriginalGrid::GetOriginalGridProc(int iZone)
{
    return originalGridProc[iZone];
}

Post_WriteTecplotByOriginalGrid *GetWriteTecplotByOriginalGrid()
{
    return writeTecplotByOriginalGrid;
}


}