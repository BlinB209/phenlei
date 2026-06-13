#include "UnstructuredOversetConfig.h"
#include "Pre_WalldistCompute.h"
#include "IO_FileName.h"
#include "Glb_Dimension.h"
#include "PHIO.h"
#include "IO_HDF5File.h"
#include "Region.h"
#include "Task.h"
#include "Task_ServerUpdateInterface.h"
#include "OversetKDTree.h"
#include "AleManager.h"
using namespace std;

namespace PHSPACE
{

OversetConfigFactoryMulti *oversetConfigFactory = nullptr;

Region *region;

void RunOversetGridConfig()
{
    if (region == nullptr)
    {
        region = new Region;
    }

    region->InitGridForOversetConfig();

    RunUnstructuredOversetGridConfig();
}

void RunUnstructuredOversetGridConfig()
{
    ReConfigUnstructOversetGrid();

    FreeUnstructOversetGridConfig();
}

void InitializeOversetGridConfig()
{
    if (oversetConfigFactory == nullptr)
    {
        oversetConfigFactory = new OversetConfigFactoryMulti();
        oversetConfigFactory->Initialize();
    }
}

void InitializeExplicitMethod()
{
 
}

void ReConfigUnstructOversetGrid()
{
    TimeSpan testTime;

    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");

    if (isOversetSlip)
    {
        ConfigOversetSlip();

        return;
    }

    int oversetGridConfigMethod = IMPLICIT_CONFIG;
    if (GlobalDataBase::IsExist("oversetGridConfigMethod", PHINT, 1))
    {
        GlobalDataBase::GetData("oversetGridConfigMethod", &oversetGridConfigMethod, PHINT, 1);
    }
    else
    {
        GlobalDataBase::UpdateData("oversetGridConfigMethod", &oversetGridConfigMethod, PHINT, 1);
    }
    if (oversetGridConfigMethod == EXPLICIT_CONFIG)
    {
        ExplicitHoleCut();
    }
    else
    {
        ImplicitHoleCut();
    }

    PrintToWindow("Run ReConfigUnstructOversetGrid sucessfully ...\n");
    testTime.ShowTimeSpanToLogFile("OversetGrid Configing");
}

void ConfigOversetSlip()
{
    InitializeOversetGridConfig();

    oversetConfigFactory->OversetSlipConfigMulti();
    int outTecplotOverset = GlobalDataBase::GetIntParaFromDB("outTecplotOverset");

    if (outTecplotOverset == 1)
    {
        oversetConfigFactory->OutPutOversetFiles();
    }
}

void ExplicitHoleCut()
{
    InitializeExplicitMethod();

    oversetConfigFactory->OversetExConfigMulti();

}

void ImplicitHoleCut()
{
    InitializeOversetGridConfig();

    oversetConfigFactory->OversetImConfigMulti();

    int outTecplotOverset = GlobalDataBase::GetIntParaFromDB("outTecplotOverset");

    if (outTecplotOverset == 1)
    {
        oversetConfigFactory->OutPutOversetFiles();
    }
}

void OutUnstructOversetGridConfig()
{
    oversetConfigFactory->OutPutOversetFiles();
}

void FreeUnstructOversetGridConfig()
{
    if (oversetConfigFactory != nullptr)
    {
        FreePointer(oversetConfigFactory);
    }
}

OversetConfigFactoryMulti *GetOversetConfigFactoryMulti()
{
    return oversetConfigFactory;
}

OversetConfigFactoryMulti::OversetConfigFactoryMulti()
{
    numberOfInterfaceNodesInBlocks = nullptr;
    keyOfInterfaceNodesInBlocks = nullptr;
    hasSymmetryBCInBlocks = nullptr;

    numberOfBlocks = 0;
    numberOfZones = 0;
    hasCalMinMax = false;

    oversetConfigMultiZone = nullptr;
    zonesSearchNodes = nullptr;
    zonesSearchLines = nullptr;
    zonesInterpolateCells = nullptr;

    testTime = 0;
    geometricDimension = GetDim();

    pMinOfBlocks = nullptr;
    pMaxOfBlocks = nullptr;
    pMinOfBlocksByLayer = nullptr;
    pMaxOfBlocksByLayer = nullptr;
}

OversetConfigFactoryMulti::~OversetConfigFactoryMulti()
{
    DelPointer(numberOfInterfaceNodesInBlocks);
    DelPointer2(keyOfInterfaceNodesInBlocks);

    DelPointer2(pMinOfBlocks);
    DelPointer2(pMaxOfBlocks);

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            if (oversetConfigMultiZone[iZone]) delete oversetConfigMultiZone[iZone];
        }
        delete zonesInterpolateCells[iZone];
        zonesInterpolateCells[iZone] = nullptr;
    }

    DelPointer(oversetConfigMultiZone);
    DelPointer(zonesSearchNodes);
    DelPointer(zonesSearchLines);
    DelPointer(zonesInterpolateCells);

    for (unsigned iBlock = 0; iBlock < innerGrid.size(); ++iBlock)
    {
        if (innerGrid[iBlock]) delete innerGrid[iBlock];
    }
    innerGrid.resize(0);
}

void OversetConfigFactoryMulti::Initialize()
{
    if (GetIfOutPutOversetVisualization() == 1)
    {
        PHMakeDir("OversetOutput");
    }
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    numberOfBlocks = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfGridGroups");
    pMinOfBlocks = NewPointer2< RDouble >(numberOfBlocks, geometricDimension);
    pMaxOfBlocks = NewPointer2< RDouble >(numberOfBlocks, geometricDimension);
    numberOfZones = PHMPI::GetNumberofGlobalZones();
    oversetConfigMultiZone = new OversetConfig * [numberOfZones];
    zonesSearchNodes = new ZoneSearchNodes * [numberOfZones];
    zonesInterpolateCells = new InterCells * [numberOfZones];

    for (int iZone = 0; iZone < this->numberOfZones; ++iZone)
    {
        oversetConfigMultiZone[iZone] = nullptr;
        zonesSearchNodes[iZone] = nullptr;
        zonesInterpolateCells[iZone] = nullptr;

        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessor == processorIndex)
        {
            OversetConfig *oversetGrid = new OversetConfig();

            Grid *grid = GetGrid(iZone, 0);
            UnstructGrid *gridUnstr = PHSPACE::UnstructGridCast(grid);
            int iBlock = grid->GetIBlock();

            oversetGrid->Initialize(gridUnstr, this->numberOfZones, iZone, iBlock);
            oversetGrid->BuildBoundaryNodeList();
            this->oversetConfigMultiZone[iZone] = oversetGrid;
        }
    }

    ComputeNodeWalldist();           //! Minimum wall distance of nodes.

    SetNodesMinDistanceInThisBlocksByWalldist();

    ComputeNodeAndCellTopology();    //! Calculate the topological relationship of point to point and volume to point.

    BuildLocalToGlobalInterfaceNodesMaps();

    ReadInInnerGrid();               //! Read in the internal grid of the object to accurately estimate whether it is in the object.

    SetBufferRegionOnOversetBoundary();    //! Set buffer on the boundary of overset.

    OutPutOversetVisualizations("0_1_Initialize");
}


void OversetConfigFactoryMulti::BuildLocalToGlobalInterfaceNodesMaps()
{
    int *numberOfNodesInBlocksLocal = new int[numberOfBlocks];
    int *numberOfNodesInBlocks = new int[numberOfBlocks];
    SetField(numberOfNodesInBlocksLocal, 0, numberOfBlocks);
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            UnstructGrid *grid = oversetConfigMultiZone[iZone]->GetGrid();
            numberOfNodesInBlocksLocal[grid->GetIBlock()] += grid->GetNTotalNode();
        }
    }
    PH_AllReduce(numberOfNodesInBlocksLocal, numberOfNodesInBlocks, numberOfBlocks, MPI_SUM);
    DelPointer(numberOfNodesInBlocksLocal);

    numberOfInterfaceNodesInBlocks = new int[numberOfBlocks];
    SetField(numberOfInterfaceNodesInBlocks, 0, numberOfBlocks);
    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        BuildLocalToGlobalInterfaceNodesMap(iBlock, numberOfNodesInBlocks);
    }

    DelPointer(numberOfNodesInBlocks);

    keyOfInterfaceNodesInBlocks = NewPointer2< int >(numberOfBlocks, numberOfInterfaceNodesInBlocks);
}

void OversetConfigFactoryMulti::BuildLocalToGlobalInterfaceNodesMap(int iBlockIn, int *numberOfNodesInBlocksIn)
{
    int &numberOfNodesInThisBlock = numberOfNodesInBlocksIn[iBlockIn];

    int *keyOfNodesInBlocks = new int[numberOfNodesInThisBlock];
    int *keyOfNodesInBlocksGlobal = new int[numberOfNodesInThisBlock];
    SetField(keyOfNodesInBlocks, 0, numberOfNodesInThisBlock);
    SetField(keyOfNodesInBlocksGlobal, 0, numberOfNodesInThisBlock);

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            UnstructGrid *grid = oversetConfigMultiZone[iZone]->GetGrid();
            int blockIndex = grid->GetIBlock();
            if (blockIndex == iBlockIn)
            {
                vector< int > &localToGlobalNodesMap = grid->GetLocalToGlobalNodesMap();
                vector< int > &interfaceBoundaryNodeList = oversetConfigMultiZone[iZone]->GetInterfaceBoundaryNodeList();
                int             numberOfinterfaceBoundaryNodes = static_cast<int>(interfaceBoundaryNodeList.size());
                for (int iNode = 0; iNode < numberOfinterfaceBoundaryNodes; ++ iNode)
                {
                    keyOfNodesInBlocks[localToGlobalNodesMap[interfaceBoundaryNodeList[iNode]]] = 1;
                }
            }
        }
    }

    PH_AllReduce(keyOfNodesInBlocks, keyOfNodesInBlocksGlobal, numberOfNodesInThisBlock, MPI_MAX);

    for (int iNode = 0; iNode < numberOfNodesInThisBlock; ++ iNode)
    {
        if (keyOfNodesInBlocksGlobal[iNode])
        {
            keyOfNodesInBlocks[iNode] = (numberOfInterfaceNodesInBlocks[iBlockIn]) ++;
        }
    }

    DelPointer(keyOfNodesInBlocksGlobal);

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            UnstructGrid *grid = oversetConfigMultiZone[iZone]->GetGrid();
            int blockIndex = grid->GetIBlock();
            if (blockIndex == iBlockIn)
            {
                vector< int > &localToGlobalNodesMap = grid->GetLocalToGlobalNodesMap();
                int *localToGlobalInterfaceNodesMap = oversetConfigMultiZone[iZone]->GetLocalToGlobalInterfaceNodesMap();
                vector< int > &interfaceBoundaryNodeList = oversetConfigMultiZone[iZone]->GetInterfaceBoundaryNodeList();
                int             numberOfinterfaceBoundaryNodes = static_cast<int>(interfaceBoundaryNodeList.size());

                for (int iNode = 0; iNode < numberOfinterfaceBoundaryNodes; ++ iNode)
                {
                    localToGlobalInterfaceNodesMap[iNode] = keyOfNodesInBlocks[localToGlobalNodesMap[interfaceBoundaryNodeList[iNode]]];
                }
            }
        }
    }

    DelPointer(keyOfNodesInBlocks);
}

void OversetConfigFactoryMulti::BuildNodesMinDistanceInThisBlocks(vector< RDouble > *&wallNodeCoordinateXIn,
    vector< RDouble > *&wallNodeCoordinateYIn,
    vector< RDouble > *&wallNodeCoordinateZIn)
{
    vector< RDouble > *wallNodeCoordinateXInBlocks = new vector< RDouble >[numberOfBlocks];
    vector< RDouble > *wallNodeCoordinateYInBlocks = new vector< RDouble >[numberOfBlocks];
    vector< RDouble > *wallNodeCoordinateZInBlocks = new vector< RDouble >[numberOfBlocks];

    int blockIndexOfThisZone = 0;
    int numberOfWallNodesInThisZones = 0;
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            blockIndexOfThisZone = oversetConfigMultiZone[iZone]->GetGrid()->GetIBlock();
            numberOfWallNodesInThisZones = static_cast<int>(wallNodeCoordinateXIn[iZone].size());
        }
        PHMPI::PH_Bcast(&numberOfWallNodesInThisZones, sizeof(int), processorIndex);

        if (numberOfWallNodesInThisZones)
        {
            if (currentProcessor != processorIndex)
            {
                wallNodeCoordinateXIn[iZone].resize(numberOfWallNodesInThisZones);
                wallNodeCoordinateYIn[iZone].resize(numberOfWallNodesInThisZones);
                wallNodeCoordinateZIn[iZone].resize(numberOfWallNodesInThisZones);
            }

            PHMPI::PH_Bcast(& wallNodeCoordinateXIn[iZone][0], numberOfWallNodesInThisZones * sizeof(RDouble), processorIndex);
            PHMPI::PH_Bcast(& wallNodeCoordinateYIn[iZone][0], numberOfWallNodesInThisZones * sizeof(RDouble), processorIndex);
            PHMPI::PH_Bcast(& wallNodeCoordinateZIn[iZone][0], numberOfWallNodesInThisZones * sizeof(RDouble), processorIndex);

            PHMPI::PH_Bcast(& blockIndexOfThisZone, sizeof(int), processorIndex);
            for (int iNode = 0; iNode < numberOfWallNodesInThisZones; ++ iNode)
            {
                wallNodeCoordinateXInBlocks[blockIndexOfThisZone].push_back(wallNodeCoordinateXIn[iZone][iNode]);
                wallNodeCoordinateYInBlocks[blockIndexOfThisZone].push_back(wallNodeCoordinateYIn[iZone][iNode]);
                wallNodeCoordinateZInBlocks[blockIndexOfThisZone].push_back(wallNodeCoordinateZIn[iZone][iNode]);
            }
        }
    }

    DelPointer(wallNodeCoordinateXIn);
    DelPointer(wallNodeCoordinateYIn);
    DelPointer(wallNodeCoordinateZIn);

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildNodesMinDistanceInThisBlock(wallNodeCoordinateXInBlocks,
                wallNodeCoordinateYInBlocks,
                wallNodeCoordinateZInBlocks);
        }
    }

    DelPointer(wallNodeCoordinateXInBlocks);
    DelPointer(wallNodeCoordinateYInBlocks);
    DelPointer(wallNodeCoordinateZInBlocks);

}
void OversetConfigFactoryMulti::SetNodesMinDistanceInThisBlocksByWalldist()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNodesMinDistanceInThisBlocksByWalldist();
        }
    }
}

void OversetConfigFactoryMulti::SetNodesMinDistanceInOtherBlocks()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNodesMinDistanceInOtherBlocks();
        }
    }
}

bool OversetConfigFactoryMulti::ReadInInnerGridH5(const string &fileName, UnstructGrid *gridIn)
{
    bool HDF5FileRead = false;

    fstream file;
    OpenFile(file, fileName, ios_base::in | ios_base::binary);
    int nZone;
    int sizeInt = sizeof(int);
    file.read((char *)&nZone, sizeInt);
    CloseFile(file);

    if (nZone <= static_cast <int> ( 1e8 ))
    {
        return HDF5FileRead;
    }

    hid_t h5File, grpGrid, grpData;
    string grpName;
    h5File = OpenHDF5File(fileName);

    ostringstream oss;
    oss << "Grid-" << "0";
    grpName = oss.str();
    grpGrid = OpenGroup(h5File, grpName);

    int *iDimensions = ReadIntOneRow(grpGrid, "iDimensions");
    gridIn->SetNTotalNode(iDimensions[0]);
    gridIn->SetNTotalFace(iDimensions[1]);
    gridIn->SetNTotalCell(iDimensions[2]);
    DelPointer(iDimensions);

    //! Read Grid Coordinates.
    grpData = OpenGroup(grpGrid, "GridCoordinates");

    RDouble *CoordinateX = ReadDoubleOneRow(grpData, "CoordinateX");
    RDouble *CoordinateY = ReadDoubleOneRow(grpData, "CoordinateY");
    RDouble *CoordinateZ = ReadDoubleOneRow(grpData, "CoordinateZ");
    gridIn->SetX(CoordinateX);
    gridIn->SetY(CoordinateY);
    gridIn->SetZ(CoordinateZ);

    gridIn->RotateAboutAxis();

    gridIn->ComputeMinMaxBox();

    H5Gclose(grpData);

    //! Read Face Topology.
    grpData = OpenGroup(grpGrid, "FaceTopology");

    int *nodeNumberOfEachFace = ReadIntOneRow(grpData, "nodeNumberOfEachFace");
    gridIn->SetNodeNumberOfEachFace(nodeNumberOfEachFace);

    int *face2Node = ReadIntOneRow(grpData, "face2Node");
    gridIn->SetFace2Node(face2Node);

    int *leftCellOfFace = ReadIntOneRow(grpData, "leftCellOfFace");
    int *rightCellOfFace = ReadIntOneRow(grpData, "rightCellOfFace");
    gridIn->SetLeftCellOfFace(leftCellOfFace);
    gridIn->SetRightCellOfFace(rightCellOfFace);

    int nTotalFace = gridIn->GetNTotalFace();
    int *nodePosi = new int[nTotalFace + 1];
    int nodeSum = 0;
    nodePosi[0] = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        nodeSum += nodeNumberOfEachFace[iFace];
        nodePosi[iFace + 1] = nodeSum;
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if (leftCellOfFace[iFace] < 0)
        {
            std::reverse(face2Node + nodePosi[iFace], face2Node + nodePosi[iFace + 1]);
            SWAP(leftCellOfFace[iFace], rightCellOfFace[iFace]);
        }
    }
    DelPointer(nodePosi);

    int nBoundFace = 0;
    ReadData(grpData, &nBoundFace, "nBoundFace");
    gridIn->SetNBoundFace(nBoundFace);

    int *BcType = ReadIntOneRow(grpData, "BcType");
    UnstructBCSet **bcr = new UnstructBCSet * [nBoundFace];
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcr[iFace] = new UnstructBCSet();
        bcr[iFace]->SetKey(BcType[iFace]);
    }
    gridIn->SetBCRecord(bcr);
    DelPointer(BcType);

    gridIn->SpecifyRightCellofBC();

    H5Gclose(grpData);

    //! Read Interface Info.
    grpData = OpenGroup(grpGrid, "InterfaceInfo");

    InterfaceInfo *interfaceInfo = 0;
    gridIn->SetInterfaceInfo(interfaceInfo);

    int nIFace = 0;
    ReadData(grpData, &nIFace, "nIFace");

    if (nIFace > 0)
    {
        interfaceInfo = new InterfaceInfo(nIFace, gridIn);
        gridIn->SetInterfaceInfo(interfaceInfo);

        int *interFace2ZoneID = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID = interfaceInfo->GetInterFace2InterFaceID();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

        int *interFace2ZoneIDTemp = ReadIntOneRow(grpData, "interFace2ZoneID");
        int *interFace2InterFaceIDTemp = ReadIntOneRow(grpData, "interFace2InterFaceID");
        int *interFace2BoundaryFaceTemp = ReadIntOneRow(grpData, "interFace2BoundaryFace");

        for (int iFace = 0; iFace < nIFace; ++ iFace)
        {
            interFace2ZoneID[iFace] = interFace2ZoneIDTemp[iFace];
            interFace2InterFaceID[iFace] = interFace2InterFaceIDTemp[iFace];
            interFace2BoundaryFace[iFace] = interFace2BoundaryFaceTemp[iFace];
        }

        DelPointer(interFace2ZoneIDTemp);
        DelPointer(interFace2InterFaceIDTemp);
        DelPointer(interFace2BoundaryFaceTemp);
    }

    H5Gclose(grpData);

    //! Read Cell Topology.
    grpData = OpenGroup(grpGrid, "CellTopology");

    int *cell2Node = ReadIntOneRow(grpData, "cell2Node");
    int *nodeNumberOfEachCell = ReadIntOneRow(grpData, "nodeNumberOfEachCell");

    gridIn->SetNodeNumberOfEachCell(nodeNumberOfEachCell);
    gridIn->SetCell2Node(cell2Node);

    H5Gclose(grpData);

    //! Read Unstr BCName.
    grpData = OpenGroup(grpGrid, "BCName");

    char *BCName = ReadCharData(grpData, "BCName");

    int nBFaces = gridIn->GetNBoundFace();
    string *BCNameList = new string[nBFaces];
    uint_long count = 0;
    for (int iFace = 0; iFace < nBFaces; ++ iFace)
    {
        while (BCName[count] != '\0')
        {
            BCNameList[iFace].push_back(BCName[count]);
            ++count;
        }
        ++count;
    }

    for (int iFace = 0; iFace < nBFaces; ++ iFace)
    {
        bcr[iFace]->SetBCName(BCNameList[iFace]);
    }

    DelPointer(BCNameList);
    DelPointer(BCName);

    H5Gclose(grpData);

    H5Gclose(grpGrid);
    H5Fclose(h5File);

    HDF5FileRead = true;
    return HDF5FileRead;
}

void OversetConfigFactoryMulti::ReadInInnerGrid()
{
    if (GetReadInAuxiliaryInnerGrid() != 1)
    {
        return;
    }

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int serverProcessor = PHMPI::GetServerSepMode();

    innerGrid.resize(numberOfBlocks);

    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        ostringstream oss;

        oss << "auxiliaryInnerGrid" << iBlock;

        int fileID = PHMPI::GetFileIndexofCurrentProcessor();
        string fileName = AddSymbolToFileName(GlobalDataBase::GetStrParaFromDB(oss.str()), "_", fileID);

        if (!FileExist(fileName))
        {
            innerGrid[iBlock] = nullptr;
            continue;
        }

        UnstructGrid *grid = new UnstructGrid();
        grid->InitGrid(new GridID(0, 0), 0, GetDim(), UNSTRUCTGRID);
        innerGrid[iBlock] = grid;
        innerGrid[iBlock]->SetIBlock(iBlock);

        bool fileIsH5File = false;
        fileIsH5File = ReadInInnerGridH5(fileName, grid);
        if (fileIsH5File)
        {
            continue;
        }

        fstream file;

        DataContainer *dataContainer = new DataContainer;

        if (currentProcessor == serverProcessor)
        {
            cout << "Read in auxiliary inner grid" << fileName << endl;

            OpenFile(file, fileName, ios_base::in | ios_base::binary);

            int numberOfZones;

            vector < int > zoneProcessorIndexContainer;
            vector < int > zoneIndexContainer;
            vector < int > zoneTypeContainer;

            PHRead(file, numberOfZones);

            zoneProcessorIndexContainer.resize(numberOfZones);
            zoneIndexContainer.resize(numberOfZones);
            zoneTypeContainer.resize(numberOfZones);

            PHRead(file, zoneProcessorIndexContainer, numberOfZones);
            PHRead(file, zoneIndexContainer, numberOfZones);
            PHRead(file, zoneTypeContainer, numberOfZones);

            dataContainer->ReadFile(file);

            CloseFile(file);
        }

        int tag = 0;

        PHMPI::PH_BcastByServerSepMode(dataContainer, tag);

        innerGrid[iBlock]->Decode(dataContainer, 0);

        FreePointer(dataContainer);
    }
}

void OversetConfigFactoryMulti::ComputeNodeWalldist()
{
    int walldistMethod = 3;

    int level = 0;

    Pre_WalldistCompute *wallDistCompute = new Pre_WalldistCompute(walldistMethod, level);

    wallDistCompute->ComputeNodeDistance();

    FreePointer(wallDistCompute);
}

void OversetConfigFactoryMulti::ComputeNodeWalldistInOtherBlock()
{
    int walldistMethod = 3;

    int level = 0;

    Pre_WalldistCompute *wallDistCompute = new Pre_WalldistCompute(walldistMethod, level);

    wallDistCompute->ComputeNodeDistanceInOtherBlock();

    FreePointer(wallDistCompute);
}

void OversetConfigFactoryMulti::ComputeNodeAndCellTopology()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->ComputeNodeAndCellTopology();
        }
    }
}

//! Get the overall globalInterfaceNodes according to the localToGlobalNodeList information of each zone.
void OversetConfigFactoryMulti::GetGlobalInterfaceNodes()
{
    if (numberOfZones == numberOfBlocks)
    {
        return;
    }

    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        GetGlobalInterfaceNodes(iBlock);
    }
}

void OversetConfigFactoryMulti::GetGlobalInterfaceNodes(int iBlock)
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    //! step1: Distinguish globalInterNodes.
    vector< int > keyNodesGlobal;
    vector< int > keyNodesGlobalAllProcessor;

    int numberOfNodesInThisBlock = numberOfGlobalNodes[iBlock];

    keyNodesGlobal.resize(numberOfNodesInThisBlock);
    keyNodesGlobalAllProcessor.resize(numberOfNodesInThisBlock);

    SetField(keyNodesGlobal, 0, numberOfNodesInThisBlock);
    SetField(keyNodesGlobalAllProcessor, 0, numberOfNodesInThisBlock);

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int sendProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == sendProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->MarkGlobalInterfaceNodes(iBlock, keyNodesGlobal);
        }
    }

    PH_AllreduceSepModeForVector(keyNodesGlobal, keyNodesGlobalAllProcessor, numberOfNodesInThisBlock, MPI_INT, MPI_MAX);

    //! step2: keyNodesGlobal becomes the number which points to the new sequence.
    int numberOfGlobalInterfaceNodes = 0;
    for (int iNode = 0; iNode < numberOfNodesInThisBlock; ++ iNode)
    {
        if (keyNodesGlobalAllProcessor[iNode] == 1)
        {
            keyNodesGlobal[iNode] = numberOfGlobalInterfaceNodes;
            numberOfGlobalInterfaceNodes++;
        }
    }

    //! step3: Get localToGlobalInterfaceNodes.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int sendProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == sendProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->CalLocalToGlobalInterfaceNodes(iBlock, keyNodesGlobal);
        }
    }

    keyActiveOfGlobalInterNodes[iBlock].resize(numberOfGlobalInterfaceNodes);
    SetField(&keyActiveOfGlobalInterNodes[iBlock][0], 0, numberOfGlobalInterfaceNodes);
}
void OversetConfigFactoryMulti::OversetSlipConfigMulti()
{
    SetInitialValue();    //! Set initial value.

    FindInterCellsByBoundary();

    CreatGlobalInterCells();

    BuildDefectCellRelationships();

    ResetOversetInformation();    //! Rebuild interpolation information.

    FindActiveCellsByBoundary();

    ResetCellCentersByBoundary();

}
void OversetConfigFactoryMulti::OversetExConfigMulti()
{
}

void OversetConfigFactoryMulti::OversetImConfigMulti()
{
    SetInitialValue();    //! Set initial value.

    NodesOperator();

    CellsOperator();
}

void OversetConfigFactoryMulti::SpecifyCellsInBodyByAuxiliaryInnerGrid()
{
    //TimeSpan *timeSpan = new TimeSpan();
    vector< GridKDTree * > innerAuxiliaryGridKDTrees;

    innerAuxiliaryGridKDTrees.resize(numberOfBlocks);

    //! step1: Interal grid establishment ADT-Tree.
    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        //! Set CurrentActionKey to 0, avoid unnecessary error in ComputMetrices.
        PHSPACE::SetCurrentActionKey(0);

        if (innerGrid[iBlock])
        {
            innerGrid[iBlock]->ComputeMetrics();
            innerAuxiliaryGridKDTrees[iBlock] = new GridKDTree(innerGrid[iBlock]);
        }
        else
        {
            innerAuxiliaryGridKDTrees[iBlock] = nullptr;
        }
    }

    BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTrees(innerAuxiliaryGridKDTrees);

    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        if (innerAuxiliaryGridKDTrees[iBlock] != 0)
        {
            delete innerAuxiliaryGridKDTrees[iBlock];
        }
    }

    OutPutOversetVisualizations("1_3_BuildInBlockCellsAndNodes");
}

void OversetConfigFactoryMulti::BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTrees(vector< GridKDTree * > &innerAuxiliaryGridKDTreesIn)
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTree(innerAuxiliaryGridKDTreesIn);
        }
    }
}

//! Expand a row of active points, points in the innerBody and buffer are set to inactive.
void OversetConfigFactoryMulti::InterFaceBoundaryOptimize()
{
    //TimeSpan *timeSpan = new TimeSpan();
    int keyEnlargeOfInterBoundary = GetKeyEnlargeOfInterBoundary() + 1;

    for (int iter = 0; iter < keyEnlargeOfInterBoundary; ++ iter)
    {
        CharacterNodeToCell();
        CharacterCellToNode();
        SetNegativeForInnerNodes();
    }

    //! Points in buffer are set to inactive.
    SetNegativeForBufferRegion();

    OutPutOversetVisualizations("1_5_InterFaceBoundaryOptimize");
}

void OversetConfigFactoryMulti::SetNegativeForInnerNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNegativeForInnerNodes();
        }
    }
}

void OversetConfigFactoryMulti::SetNegativeForInnerAndInitialInterNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNegativeForInnerAndInitialInterNodes();
        }
    }
}

void OversetConfigFactoryMulti::NodesOperator()
{
    BuildNodesMinDistanceInOtherBlocks();

    if (GetReadInAuxiliaryInnerGrid() == 1)
    {
        SpecifyCellsInBodyByAuxiliaryInnerGrid();
    }
    else
    {
        SpecifySpecicalCellsInteractWithWallBoundary();
    }

    SpecifyInitialInterNodesByDistance();

    FillTheActiveRegion();          //! Physical boundary advancement.

    InterFaceBoundaryOptimize();    //! Hole boundary optimization (expand active areas).
}

void OversetConfigFactoryMulti::ReComputeWallDistance()
{
    int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    int key = 0;

    for (int iBody = 0; iBody < numberOfMovingBodies; ++iBody)
    {
        ostringstream oss;

        oss << "morphing_" << iBody;

        int morphingManner = GlobalDataBase::GetIntParaFromDB(oss.str());

        if (morphingManner != 0)
        {
            key = 1;
            break;
        }
    }

    //! No need to solve the wall distance again if the overall objects don't deform.
    if (key == 0)
    {
        return;
    }

    ComputeNodeWalldist();
}

void OversetConfigFactoryMulti::BuildNodesMinDistanceInOtherBlocks()
{

    BuildKDTreesBySearchRegions();

    BuildNodesMinDistanceInOtherBlocksByKDTrees();

    OutPutOversetVisualizations("1_1_BuildNodesMinDistanceInOtherBlocks");
}

void OversetConfigFactoryMulti::BuildNodesMinDistanceInOtherBlocks2()
{

    BuildKDTreesBySearchRegions();

    BuildZonesSearchNodes();

    ComputeNodeWalldistInOtherBlock();

    SetNodesMinDistanceInOtherBlocks();

    OutPutOversetVisualizations("1_1_BuildNodesMinDistanceInOtherBlocks");
}

void OversetConfigFactoryMulti::BuildNodesMinDistanceInOtherBlocksByKDTrees()
{

    BuildZonesSearchNodes();

    SearchZonesSearchNodesInKDTrees();

    AllReduceNodesMinDistanceInOtherBlock();

    SetNodesMinDistanceInOtherBlockByZonesNodes();
}

void OversetConfigFactoryMulti::SearchZonesSearchNodesInKDTrees()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SearchZonesSearchNodesInKDTree();
        }
    }
}

void OversetConfigFactoryMulti::AllReduceNodesMinDistanceInOtherBlock()
{
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        ZoneSearchNodes *&zoneSearchNodes = zonesSearchNodes[iZone];
        int numberOfNodes = zoneSearchNodes->GetNumberOfNodes();
        RDouble *nodesMinDistanceToOtherBlock = zoneSearchNodes->GetNodesMinDistanceToOtherBlock();

        if (numberOfNodes)
        {
            RDouble *nodesMinDistanceInOtherBlockGlobal = new RDouble[numberOfNodes];

            PH_AllReduce(nodesMinDistanceToOtherBlock, nodesMinDistanceInOtherBlockGlobal, numberOfNodes, PH_MIN);

            SetField(nodesMinDistanceToOtherBlock, nodesMinDistanceInOtherBlockGlobal, numberOfNodes);

            DelPointer(nodesMinDistanceInOtherBlockGlobal);
        }
    }
}

void OversetConfigFactoryMulti::SetNodesMinDistanceInOtherBlockByZonesNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        ZoneSearchNodes *&zoneSearchNodes = zonesSearchNodes[iZone];
        if (currentProcessorIndex == processorIndex)
        {
            RDouble *nodesMinDistanceToOtherBlock = zoneSearchNodes->GetNodesMinDistanceToOtherBlock();
            oversetConfigMultiZone[iZone]->CopyMinDist(nodesMinDistanceToOtherBlock);
        }
        FreePointer(zonesSearchNodes[iZone]);
    }
}

void OversetConfigFactoryMulti::CellsOperator()
{
    FindActiveCells();    //! Identify active cells.

    FindInterCells();     //! Identify interpolation type.

    BuildInterpolateRelationships();

    CheckInterpolateRelationships();

    CollectEffectiveDonorCells();

    BuildDefectCellRelationships();

    ResetOversetInformation();    //! Rebuild interpolation information.
}

void OversetConfigFactoryMulti::SetInitialValue()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetInitialValue();
        }
        if (nullptr != zonesInterpolateCells[iZone])
        {
            FreePointer(zonesInterpolateCells[iZone]);
        }
    }
}

void OversetConfigFactoryMulti::CalMinBoxForSearch()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildMinMaxBoxOfZoneBySearchRegion();
        }
    }
}

void OversetConfigFactoryMulti::SetBufferRegionOnOversetBoundary()
{
    //! Query physical boundary points and overset boundary points.
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetInitialValueForBufferRigion();
        }
    }

    CharacterNodeToCell();    //! A row of points and cells are left in the buffer.
    //! Increase the interpolation area to improve the reconstruction accuracy of interpolation cells, two row of points and cells will be left in the buffer (There's no necessary actually).
    int keyEnlargeOfActiveNodes = GetKeyEnlargeOfNodes();
    for (int iter = 0; iter < keyEnlargeOfActiveNodes; ++ iter)
    {
        CharacterCellToNode();

        CharacterNodeToCell();
    }

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetKeyBufferByKeyActive();
        }
    }
}

//! Identify searchRegion, starting from search point and expand twice.
//! Identify search point by judging whether it's in other blocks.
void OversetConfigFactoryMulti::BuildKDTreesBySearchRegions()
{
    BuildMinMaxBoxOfBlocks();

    SetActivePropertyForSearchRegions();

    //! Why the value is 2 (Ensure that searchRegion covers all blocks except this block).

    for (int i = 0; i < 2; ++ i)
    {
        CharacterNodeToCell();

        CharacterCellToNode();
    }

    SetSearchRegionByActivePropertys();

    BuildMinMaxBoxOfZoneBySearchRegions();

    BuildKDTreeByMinMaxBoxOfZones();

}

void OversetConfigFactoryMulti::SetActivePropertyForSearchRegions()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetActivePropertyForSearchRegion();
        }
    }
}

void OversetConfigFactoryMulti::SetSearchRegionByActivePropertys()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetSearchRegionByActiveProperty();
        }
    }
}

void OversetConfigFactoryMulti::BuildMinMaxBoxOfZoneBySearchRegions()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildMinMaxBoxOfZoneBySearchRegion();
        }
    }
}

void OversetConfigFactoryMulti::BuildKDTreeByMinMaxBoxOfZones()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildKDTreeByMinMaxBoxOfZone();
        }
    }
}

void OversetConfigFactoryMulti::BuildKDTreeByActiveCells() 
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildKDTreeByActiveCell();
        }
    }
}

void OversetConfigFactoryMulti::CollectEffectiveDonorCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->CollectEffectiveDonorCell();
        }
    }
}

void OversetConfigFactoryMulti::BuildMinMaxBoxOfBlocks()
{
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    RDouble **minBoxOfBlocksLocal = NewPointer2< RDouble >(numberOfBlocks, geometricDimension);
    RDouble **maxBoxOfBlocksLocal = NewPointer2< RDouble >(numberOfBlocks, geometricDimension);

    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        if (isAle)
        {
            int RBDMethod = PHSPACE::GetIntegerParameterFromDataBase(iBlock, "RBDMethod");
            if (hasCalMinMax && RBDMethod == 0)
            {
                continue;
            }
        }

        SetField(pMinOfBlocks[iBlock], LARGE, geometricDimension);
        SetField(pMaxOfBlocks[iBlock], -LARGE, geometricDimension);
        SetField(minBoxOfBlocksLocal[iBlock], LARGE, geometricDimension);
        SetField(maxBoxOfBlocksLocal[iBlock], -LARGE, geometricDimension);
    }

    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            int blockIndex = oversetConfigMultiZone[iZone]->GetGrid()->GetIBlock();
            if (isAle)
            {
                int RBDMethod = PHSPACE::GetIntegerParameterFromDataBase(blockIndex, "RBDMethod");
                if (hasCalMinMax && RBDMethod == 0)
                {
                    continue;
                }
            }
            oversetConfigMultiZone[iZone]->BuildMinMaxBoxOfZone();

            RDouble *minBoxOfZone = oversetConfigMultiZone[iZone]->GetPMin();
            RDouble *maxBoxOfZone = oversetConfigMultiZone[iZone]->GetPMax();

            for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
            {
                RDouble &minInThisDimension = minBoxOfBlocksLocal[blockIndex][iDimension];
                RDouble &maxInThisDimension = maxBoxOfBlocksLocal[blockIndex][iDimension];
                minInThisDimension = MIN(minInThisDimension, minBoxOfZone[iDimension]);
                maxInThisDimension = MAX(maxInThisDimension, maxBoxOfZone[iDimension]);
            }
        }
    }

    RDouble eps = GetToleranceForOversetBox();
    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        if (isAle)
        {
            int RBDMethod = PHSPACE::GetIntegerParameterFromDataBase(iBlock, "RBDMethod");
            if (hasCalMinMax && RBDMethod == 0)
            {
                continue;
            }
        }

        PH_AllReduce(minBoxOfBlocksLocal[iBlock], pMinOfBlocks[iBlock], geometricDimension, PH_MIN);
        PH_AllReduce(maxBoxOfBlocksLocal[iBlock], pMaxOfBlocks[iBlock], geometricDimension, PH_MAX);

        for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
        {
            pMinOfBlocks[iBlock][iDimension] -= eps;
            pMaxOfBlocks[iBlock][iDimension] += eps;
        }
    }
    hasCalMinMax = true;
    DelPointer2(minBoxOfBlocksLocal);
    DelPointer2(maxBoxOfBlocksLocal);
}

void OversetConfigFactoryMulti::BuildZonesSearchNodes()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        int *keySearchOfNodes = nullptr;
        int   numberOfNodes = 0;

        int   blockIndex = 0;
        int   numberOfSearchNodes = 0;

        if (currentProcessor == processorIndex)
        {
            keySearchOfNodes = oversetConfigMultiZone[iZone]->GetKeySearchOfNodes();
            numberOfNodes = oversetConfigMultiZone[iZone]->GetGrid()->GetNTotalNode();

            blockIndex = oversetConfigMultiZone[iZone]->GetGrid()->GetIBlock();

            for (int iNode = 0; iNode < numberOfNodes; ++iNode)
            {
                if (keySearchOfNodes[iNode])
                {
                    ++ numberOfSearchNodes;
                }
            }
        }

        PHMPI::PH_Bcast(&numberOfSearchNodes, sizeof(int), processorIndex, 1);

        PHMPI::PH_Bcast(&blockIndex, sizeof(int), processorIndex, 1);

        if (zonesSearchNodes[iZone])
        {
            FreePointer(zonesSearchNodes[iZone]);
        }
        zonesSearchNodes[iZone] = new ZoneSearchNodes(numberOfSearchNodes, blockIndex);

        RDouble *minBoxOfZoneSearchNodes = zonesSearchNodes[iZone]->GetPMin();
        RDouble *maxBoxOfZoneSearchNodes = zonesSearchNodes[iZone]->GetPMax();

        if (currentProcessor == processorIndex)
        {
            RDouble *minBoxOfZone = oversetConfigMultiZone[iZone]->GetPMin();
            RDouble *maxBoxOfZone = oversetConfigMultiZone[iZone]->GetPMax();

            SetField(minBoxOfZoneSearchNodes, minBoxOfZone, geometricDimension);
            SetField(maxBoxOfZoneSearchNodes, maxBoxOfZone, geometricDimension);
        }
        PHMPI::PH_Bcast(minBoxOfZoneSearchNodes, geometricDimension * sizeof(RDouble), processorIndex);
        PHMPI::PH_Bcast(maxBoxOfZoneSearchNodes, geometricDimension * sizeof(RDouble), processorIndex);

        if (numberOfSearchNodes)
        {
            RDouble *searchNodeCoordinateX = zonesSearchNodes[iZone]->GetX();
            RDouble *searchNodeCoordinateY = zonesSearchNodes[iZone]->GetY();
            RDouble *searchNodeCoordinateZ = zonesSearchNodes[iZone]->GetZ();

            if (currentProcessor == processorIndex)
            {
                RDouble *nodeCoordinateX = oversetConfigMultiZone[iZone]->GetGrid()->GetX();
                RDouble *nodeCoordinateY = oversetConfigMultiZone[iZone]->GetGrid()->GetY();
                RDouble *nodeCoordinateZ = oversetConfigMultiZone[iZone]->GetGrid()->GetZ();

                numberOfSearchNodes = 0;

                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    if (keySearchOfNodes[iNode])
                    {
                        searchNodeCoordinateX[numberOfSearchNodes] = nodeCoordinateX[iNode];
                        searchNodeCoordinateY[numberOfSearchNodes] = nodeCoordinateY[iNode];
                        searchNodeCoordinateZ[numberOfSearchNodes] = nodeCoordinateZ[iNode];

                        ++numberOfSearchNodes;
                    }
                }
            }
            PHMPI::PH_Bcast(searchNodeCoordinateX, numberOfSearchNodes * sizeof(RDouble), processorIndex);
            PHMPI::PH_Bcast(searchNodeCoordinateY, numberOfSearchNodes * sizeof(RDouble), processorIndex);
            PHMPI::PH_Bcast(searchNodeCoordinateZ, numberOfSearchNodes * sizeof(RDouble), processorIndex);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void OversetConfigFactoryMulti::SearchZonesSearchNodesInLineBoxes()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SearchZonesSearchNodesInLineBox();
        }
    }
}

void OversetConfigFactoryMulti::SearchZonesSearchNodesInLineBoxesByLayer()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SearchZonesSearchNodesInLineBoxByLayer();
        }
    }
}

void OversetConfigFactoryMulti::SpecifyInitialInterNodesByDistance()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SpecifyInitialInterNodesByDistance();
        }
    }

    OutPutOversetVisualizations("1_2_BuildKeyBasicActiveNodesByDistances");
}

void OversetConfigFactoryMulti::SpecifySpecicalCellsInteractWithWallBoundary()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SpecifySpecicalCellsInteractWithWallBoundary();
        }
    }
    OutPutOversetVisualizations("1_3_BuildInBlockCellsAndNodes");
}

void OversetConfigFactoryMulti::FillTheActiveRegion()
{
    int numberOfNewNodesPresentCycle = 0;
    int sum = 0;

    SetNodeCharacterOnPhysicalBoundary();
    do
    {
        numberOfNewNodesPresentCycle = CharacterNodeToValidNodeByCell();

        //! mpi_reduce.
        PH_AllReduce(&numberOfNewNodesPresentCycle, &sum, 1, PH_SUM);
        numberOfNewNodesPresentCycle = sum;

    } while (numberOfNewNodesPresentCycle > 0);

    FreeKeyBasicActiveNodes();

    OutPutOversetVisualizations("1_4_BuildActiveRegion");
}

void OversetConfigFactoryMulti::SetNodeCharacterOnPhysicalBoundary()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNodeCharacterOnPhysicalBoundary();
        }
    }
}

void OversetConfigFactoryMulti::SetNegativeForBufferRegion()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNegativeForBufferRegion();
        }
    }
}

void OversetConfigFactoryMulti::BuildInterpolateRelationships()
{
    CreatGlobalInterCells();

    FindInterpolateRelationships();

    AllReduceZonesInterpolateCells();
}


void OversetConfigFactoryMulti::BuildDefectCellRelationships()
{
    AllReduceZonesDefectCells();

    FindDefectRelationships();

    AllReduceZonesInterpolateCells();
}

void OversetConfigFactoryMulti::FindInterpolateRelationships()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->FindInterpolateRelationship();
        }
    }
}

void OversetConfigFactoryMulti::FindDefectRelationships()
{
    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            if (isOversetSlip)
            {
                oversetConfigMultiZone[iZone]->FindDefectRelationshipByFaces();
            }
            else
            {
                oversetConfigMultiZone[iZone]->FindDefectRelationshipByExpend();
            }
        }
    }
}

void OversetConfig::FindInterpolateRelationship()
{
    RDouble *nodeCoordinate = new RDouble[THREE_D];
    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    int blockIndex = grid->GetIBlock();
    int zoneIndex = grid->GetZoneID();
    InterCells **zonesInterpolateCells = GetOversetConfigFactoryMulti()->GetZonesInterpolateCells();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *&zoneInterpolateCells = zonesInterpolateCells[iZone];

        if (blockIndex != zoneInterpolateCells->GetBlockIndex())
        {
            int numberOfInterpolateCells = zoneInterpolateCells->GetNumberOfInterCells();

            RDouble *interpolateCellCenterX = zoneInterpolateCells->GetX();
            RDouble *interpolateCellCenterY = zoneInterpolateCells->GetY();
            RDouble *interpolateCellCenterZ = zoneInterpolateCells->GetZ();

            int *donorZoneIndex = zoneInterpolateCells->GetDonorZoneIndex();
            int *donorCellIndex = zoneInterpolateCells->GetDonorCellIndex();

            int *donorLevel = zoneInterpolateCells->GetKey();
            RDouble *cellMinDistance = zoneInterpolateCells->GetMinDist();

            for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
            {
                nodeCoordinate[0] = interpolateCellCenterX[iCell];
                nodeCoordinate[1] = interpolateCellCenterY[iCell];
                nodeCoordinate[2] = interpolateCellCenterZ[iCell];

                int donorCellIndexOfThisCell = gridKDTree->SearchNodeInKDTree(nodeCoordinate);

                if (donorCellIndexOfThisCell != -1)
                {
                    int donorLevelNow = -3;

                    int &donorLevelOfThisCell = donorLevel[iCell];
                    RDouble &cellMinDistanceOfThisCell = cellMinDistance[iCell];

                    switch (keyActiveOfCells[donorCellIndexOfThisCell])
                    {
                        case -1:
                            donorLevelNow = -1;
                            break;
                        case   0:
                            donorLevelNow = -2;
                            break;
                        case   1:
                            donorLevelNow = 1;
                            break;
                    }

                    if (donorLevelNow > donorLevelOfThisCell)
                    {
                        donorLevelOfThisCell = donorLevelNow;
                        cellMinDistanceOfThisCell = ComputeDistanceInDonorZone(nodeCoordinate, donorCellIndexOfThisCell);
                        donorZoneIndex[iCell] = zoneIndex;
                        donorCellIndex[iCell] = donorCellIndexOfThisCell;
                    }
                    else if (donorLevelNow == donorLevelOfThisCell)
                    {
                        RDouble cellMinDistanceNow = ComputeDistanceInDonorZone(nodeCoordinate, donorCellIndexOfThisCell);
                        if (cellMinDistanceNow < cellMinDistanceOfThisCell)
                        {
                            cellMinDistanceOfThisCell = cellMinDistanceNow;
                            donorZoneIndex[iCell] = zoneIndex;
                            donorCellIndex[iCell] = donorCellIndexOfThisCell;
                        }
                    }
                }
                else
                {

                }
            }
        }
    }

    FreePointer(gridKDTree);
    DelPointer(nodeCoordinate);
}

void OversetConfig::FindDefectRelationship()
{
    RDouble *nodeCoordinate = new RDouble[THREE_D];
    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    bool ifWallDist = true;
    if (0 == donorCellMethod)
    {
        ifWallDist = false;
    }
    int blockIndex = grid->GetIBlock();
    int zoneIndex = grid->GetZoneID();
    OversetConfigFactoryMulti *overs = GetOversetConfigFactoryMulti();
    InterCells **zonesInterpolateCells = overs->GetZonesInterpolateCells();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *&zoneInterpolateCells = zonesInterpolateCells[iZone];

        if (blockIndex != zoneInterpolateCells->GetBlockIndex())
        {
            int numberOfInterpolateCells = zoneInterpolateCells->GetNumberOfInterCells();
            vector< int > &interpolateCellList = (overs->GetOversetConfig())[iZone]->GetInterpolateCellList();
            vector< int > &defectCellList = (overs->GetOversetConfig())[iZone]->GetDefectCellList();

            RDouble *interpolateCellCenterX = zoneInterpolateCells->GetX();
            RDouble *interpolateCellCenterY = zoneInterpolateCells->GetY();
            RDouble *interpolateCellCenterZ = zoneInterpolateCells->GetZ();

            int *donorZoneIndex = zoneInterpolateCells->GetDonorZoneIndex();
            int *donorCellIndex = zoneInterpolateCells->GetDonorCellIndex();
            int *isDefectCell   = zoneInterpolateCells->IsDefectCell();

            int *donorLevel = zoneInterpolateCells->GetKey();
            RDouble *cellMinDistance = zoneInterpolateCells->GetMinDist();

            for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
            {
                if (isDefectCell[iCell])
                {
                    nodeCoordinate[0] = interpolateCellCenterX[iCell];
                    nodeCoordinate[1] = interpolateCellCenterY[iCell];
                    nodeCoordinate[2] = interpolateCellCenterZ[iCell];
                    if (!cellKDTree) continue;
                    KDRes *minResults = KDNearest(cellKDTree, nodeCoordinate);
                    int   donorCellIndexOfThisCell = minResults->itemData();

                    if (donorCellIndexOfThisCell != -1)
                    {
                        int donorLevelNow = -3;

                        int &donorLevelOfThisCell = donorLevel[iCell];
                        RDouble &cellMinDistanceOfThisCell = cellMinDistance[iCell];

                        switch (keyActiveOfCells[donorCellIndexOfThisCell])
                        {
                            case -1:
                                donorLevelNow = -1;    //!donorCells is INTERPOLATION
                                break;
                            case   0:
                                donorLevelNow = -2;    //!donorCells is INACTIVE
                                break;
                            case   1:
                                donorLevelNow = 1;     //!donorCells is ACTIVE
                                break;
                        }

                        if (donorLevelNow > donorLevelOfThisCell)
                        {
                            donorLevelOfThisCell = donorLevelNow;
                            cellMinDistanceOfThisCell = ComputeDistanceInDonorZone(nodeCoordinate, donorCellIndexOfThisCell, ifWallDist);
                            donorZoneIndex[iCell] = zoneIndex;
                            donorCellIndex[iCell] = donorCellIndexOfThisCell;
                        }
                        else if (donorLevelNow == donorLevelOfThisCell)
                        {
                            RDouble cellMinDistanceNow = ComputeDistanceInDonorZone(nodeCoordinate, donorCellIndexOfThisCell, ifWallDist);
                            if (cellMinDistanceNow < cellMinDistanceOfThisCell)
                            {
                                cellMinDistanceOfThisCell = cellMinDistanceNow;
                                donorZoneIndex[iCell] = zoneIndex;
                                donorCellIndex[iCell] = donorCellIndexOfThisCell;
                            }
                        }
                    }

                    FreeKDRes(minResults);
                }
            }
        }
    }

    DelPointer(nodeCoordinate);
}

void OversetConfig::FindDefectRelationshipByFaces()
{
    RDouble *nodeCoordinate = new RDouble[THREE_D];
    int *keyActiveOfCells = this->GetKeyActiveOfCells();
    RDouble eps = GlobalDataBase::GetDoubleParaFromDB("toleranceForOversetBox");
   
    int blockIndex      = grid->GetIBlock();
    int zoneIndex       = grid->GetZoneID();
    int *leftCellofFace = grid->GetLeftCellOfFace();
   
    OversetConfigFactoryMulti *overs = GetOversetConfigFactoryMulti();
    InterCells **zonesInterpolateCells = overs->GetZonesInterpolateCells();

    set< int > &interpolateExpendList = GetDonorFaceListList();
    int numberOfExpendListCells = static_cast<int>( ( GetDonorFaceListList() ).size() );

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *&zoneInterpolateCells = zonesInterpolateCells[iZone];

        if (blockIndex != zoneInterpolateCells->GetBlockIndex())
        {
            int numberOfInterpolateCells = zoneInterpolateCells->GetNumberOfInterCells();
            vector< int > &interpolateCellList = (overs->GetOversetConfig())[iZone]->GetInterpolateCellList();
            vector< int > &defectCellList = (overs->GetOversetConfig())[iZone]->GetDefectCellList();

            RDouble *interpolateCellCenterX = zoneInterpolateCells->GetX();
            RDouble *interpolateCellCenterY = zoneInterpolateCells->GetY();
            RDouble *interpolateCellCenterZ = zoneInterpolateCells->GetZ();

            int *donorZoneIndex = zoneInterpolateCells->GetDonorZoneIndex();
            int *donorCellIndex = zoneInterpolateCells->GetDonorCellIndex();
            int *isDefectCell   = zoneInterpolateCells->IsDefectCell();

            int *donorLevel = zoneInterpolateCells->GetKey();
            RDouble *cellMinDistance = zoneInterpolateCells->GetMinDist();

            for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
            {
                if (isDefectCell[iCell])
                {
                    nodeCoordinate[0] = interpolateCellCenterX[iCell];
                    nodeCoordinate[1] = interpolateCellCenterY[iCell];
                    nodeCoordinate[2] = interpolateCellCenterZ[iCell];

                    RDouble &cellMinDistanceOfThisCell = cellMinDistance[iCell];
                    for (set<int>::iterator iter = interpolateExpendList.begin(); iter != interpolateExpendList.end(); ++iter)
                    {
                        int donorFaceIndexOfThisZone = *iter;
                        int le = leftCellofFace[donorFaceIndexOfThisZone];
                        RDouble cellMinDistanceNow = ComputeDistanceInDonorZoneForFace(nodeCoordinate, donorFaceIndexOfThisZone);
                        if (cellMinDistanceNow < cellMinDistanceOfThisCell)
                        {
                            cellMinDistanceOfThisCell = cellMinDistanceNow;
                            donorZoneIndex[iCell] = zoneIndex;
                            donorCellIndex[iCell] = le;
                        }
                    }
                }
            }
        }
    }

    DelPointer(nodeCoordinate);
}

void OversetConfig::FindDefectRelationshipByExpend()
{
    RDouble *nodeCoordinate = new RDouble[THREE_D];
    int *keyActiveOfCells = this->GetKeyActiveOfCells();
    RDouble eps = GlobalDataBase::GetDoubleParaFromDB("toleranceForOversetBox");
    bool ifWallDist = true;
    if (0 == donorCellMethod)
    {
        ifWallDist = false;
    }
    int blockIndex      = grid->GetIBlock();
    int zoneIndex       = grid->GetZoneID();
    int *leftCellofFace = grid->GetLeftCellOfFace();

    OversetConfigFactoryMulti *overs = GetOversetConfigFactoryMulti();
    InterCells **zonesInterpolateCells = overs->GetZonesInterpolateCells();

    set< int > &interpolateExpendList = GetInterCellExpendList();
    int numberOfExpendListCells = static_cast<int>( ( GetInterCellExpendList() ).size() );

    if (0 == numberOfExpendListCells)
    {
        return;
    }
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *&zoneInterpolateCells = zonesInterpolateCells[iZone];

        if (blockIndex != zoneInterpolateCells->GetBlockIndex())
        {
            int numberOfInterpolateCells = zoneInterpolateCells->GetNumberOfInterCells();
            vector< int > &interpolateCellList = (overs->GetOversetConfig())[iZone]->GetInterpolateCellList();
            vector< int > &defectCellList = (overs->GetOversetConfig())[iZone]->GetDefectCellList();

            RDouble *interpolateCellCenterX = zoneInterpolateCells->GetX();
            RDouble *interpolateCellCenterY = zoneInterpolateCells->GetY();
            RDouble *interpolateCellCenterZ = zoneInterpolateCells->GetZ();

            int *donorZoneIndex = zoneInterpolateCells->GetDonorZoneIndex();
            int *donorCellIndex = zoneInterpolateCells->GetDonorCellIndex();
            int *isDefectCell   = zoneInterpolateCells->IsDefectCell();

            int *donorLevel = zoneInterpolateCells->GetKey();
            RDouble *cellMinDistance = zoneInterpolateCells->GetMinDist();

            for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
            {
                if (isDefectCell[iCell])
                {
                    nodeCoordinate[0] = interpolateCellCenterX[iCell];
                    nodeCoordinate[1] = interpolateCellCenterY[iCell];
                    nodeCoordinate[2] = interpolateCellCenterZ[iCell];

                    RDouble &cellMinDistanceOfThisCell = cellMinDistance[iCell];
                    for (set<int>::iterator iter = interpolateExpendList.begin(); iter != interpolateExpendList.end(); ++iter)
                    {
                        int donorCellIndexOfThisZone = *iter;
                        RDouble cellMinDistanceNow = ComputeDistanceInDonorZone(nodeCoordinate, donorCellIndexOfThisZone,ifWallDist);
                        if (cellMinDistanceNow < cellMinDistanceOfThisCell)
                        {
                            cellMinDistanceOfThisCell = cellMinDistanceNow;
                            donorZoneIndex[iCell] = zoneIndex;
                            donorCellIndex[iCell] = donorCellIndexOfThisZone;
                        }
                    }
                    donorLevel[iCell] = ACTIVE;
                }
            }
        }
    }

    DelPointer(nodeCoordinate);

}

void OversetConfigFactoryMulti::CreatGlobalInterCells()
{
    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");
    if (isOversetSlip)
    {
        int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
        for (int iZone = 0; iZone < numberOfZones; ++iZone)
        {
            int   blockIndex = 0;
            int   numberOfInterpolateCells = 0;

            int zoneProcessorIndexOfThisZone = PHMPI::GetZoneProcessorID(iZone);
            if (currentProcessorIndex == zoneProcessorIndexOfThisZone)
            {
                blockIndex = oversetConfigMultiZone[iZone]->GetGrid()->GetIBlock();
                    numberOfInterpolateCells = static_cast<int>( ( oversetConfigMultiZone[iZone]->GetInterpolateFaceList() ).size() );
            }

            PHMPI::PH_Bcast(&blockIndex, sizeof(int), zoneProcessorIndexOfThisZone);
            PHMPI::PH_Bcast(&numberOfInterpolateCells, sizeof(int), zoneProcessorIndexOfThisZone);

            if (zonesInterpolateCells[iZone])
            {
                FreePointer(zonesInterpolateCells[iZone]);
            }

            zonesInterpolateCells[iZone] = new InterCells(numberOfInterpolateCells, blockIndex);

            if (numberOfInterpolateCells)
            {
                RDouble *interpolateCellCenterX = zonesInterpolateCells[iZone]->GetX();
                RDouble *interpolateCellCenterY = zonesInterpolateCells[iZone]->GetY();
                RDouble *interpolateCellCenterZ = zonesInterpolateCells[iZone]->GetZ();

                if (currentProcessorIndex == zoneProcessorIndexOfThisZone)
                {
                    vector< int > &interpolateCellList = oversetConfigMultiZone[iZone]->GetInterpolateFaceList();

                    RDouble *cellCenterX = oversetConfigMultiZone[iZone]->GetGrid()->GetFaceCenterX();
                    RDouble *cellCenterY = oversetConfigMultiZone[iZone]->GetGrid()->GetFaceCenterY();
                    RDouble *cellCenterZ = oversetConfigMultiZone[iZone]->GetGrid()->GetFaceCenterZ();

                    for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
                    {
                        int &faceIndex = interpolateCellList[iCell];

                        interpolateCellCenterX[iCell] = cellCenterX[faceIndex];
                        interpolateCellCenterY[iCell] = cellCenterY[faceIndex];
                        interpolateCellCenterZ[iCell] = cellCenterZ[faceIndex];
                    }
                }

                PHMPI::PH_Bcast(interpolateCellCenterX, numberOfInterpolateCells * sizeof(RDouble), zoneProcessorIndexOfThisZone);
                PHMPI::PH_Bcast(interpolateCellCenterY, numberOfInterpolateCells * sizeof(RDouble), zoneProcessorIndexOfThisZone);
                PHMPI::PH_Bcast(interpolateCellCenterZ, numberOfInterpolateCells * sizeof(RDouble), zoneProcessorIndexOfThisZone);

                int *isDefectCell =  zonesInterpolateCells[iZone]->IsDefectCell();
                SetField(isDefectCell, 1, numberOfInterpolateCells);
            }
        }
    }
    else
    {
        int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
        for (int iZone = 0; iZone < numberOfZones; ++iZone)
        {
            int   blockIndex = 0;
            int   numberOfInterpolateCells = 0;

            int zoneProcessorIndexOfThisZone = PHMPI::GetZoneProcessorID(iZone);
            if (currentProcessorIndex == zoneProcessorIndexOfThisZone)
            {
                blockIndex = oversetConfigMultiZone[iZone]->GetGrid()->GetIBlock();
            numberOfInterpolateCells = static_cast<int>( ( oversetConfigMultiZone[iZone]->GetInterpolateCellList() ).size() );
            }

            PHMPI::PH_Bcast(&blockIndex, sizeof(int), zoneProcessorIndexOfThisZone);
            PHMPI::PH_Bcast(&numberOfInterpolateCells, sizeof(int), zoneProcessorIndexOfThisZone);

            if (zonesInterpolateCells[iZone])
            {
                FreePointer(zonesInterpolateCells[iZone]);
            }

            zonesInterpolateCells[iZone] = new InterCells(numberOfInterpolateCells, blockIndex);

            if (numberOfInterpolateCells)
            {
                RDouble *interpolateCellCenterX = zonesInterpolateCells[iZone]->GetX();
                RDouble *interpolateCellCenterY = zonesInterpolateCells[iZone]->GetY();
                RDouble *interpolateCellCenterZ = zonesInterpolateCells[iZone]->GetZ();

                if (currentProcessorIndex == zoneProcessorIndexOfThisZone)
                {
                    vector< int > &interpolateCellList = oversetConfigMultiZone[iZone]->GetInterpolateCellList();

                    RDouble *cellCenterX = oversetConfigMultiZone[iZone]->GetGrid()->GetCellCenterX();
                    RDouble *cellCenterY = oversetConfigMultiZone[iZone]->GetGrid()->GetCellCenterY();
                    RDouble *cellCenterZ = oversetConfigMultiZone[iZone]->GetGrid()->GetCellCenterZ();

                    for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
                    {
                        int &cellIndex = interpolateCellList[iCell];

                        interpolateCellCenterX[iCell] = cellCenterX[cellIndex];
                        interpolateCellCenterY[iCell] = cellCenterY[cellIndex];
                        interpolateCellCenterZ[iCell] = cellCenterZ[cellIndex];
                    }
                }

                PHMPI::PH_Bcast(interpolateCellCenterX, numberOfInterpolateCells * sizeof(RDouble), zoneProcessorIndexOfThisZone);
                PHMPI::PH_Bcast(interpolateCellCenterY, numberOfInterpolateCells * sizeof(RDouble), zoneProcessorIndexOfThisZone);
                PHMPI::PH_Bcast(interpolateCellCenterZ, numberOfInterpolateCells * sizeof(RDouble), zoneProcessorIndexOfThisZone);
            }
        }
    }
}

void OversetConfigFactoryMulti::AllReduceZonesInterpolateCells()
{
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *zoneInterpolateCells = zonesInterpolateCells[iZone];

        int numberOfInterpolateCells = zoneInterpolateCells->GetNumberOfInterCells();

        if (numberOfInterpolateCells)
        {
            int *donorZoneIndex = zoneInterpolateCells->GetDonorZoneIndex();
            int *donorCellIndex = zoneInterpolateCells->GetDonorCellIndex();
            int *donorLevel = zoneInterpolateCells->GetKey();
            RDouble *cellMinDistance = zoneInterpolateCells->GetMinDist();

            int *donorZoneIndexGlobal = new int[numberOfInterpolateCells];
            int *donorCellIndexGlobal = new int[numberOfInterpolateCells];
            int *donorLevelGlobal = new int[numberOfInterpolateCells];
            RDouble *cellMinDistanceGlobal = new RDouble[numberOfInterpolateCells];

            PH_AllReduce(donorLevel, donorLevelGlobal, numberOfInterpolateCells, PH_MAX);
            for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
            {
                if (donorLevelGlobal[iCell] > donorLevel[iCell])
                {
                    cellMinDistance[iCell] = LARGE;
                    donorZoneIndex[iCell] = -1;
                }
            }

            PH_AllReduce(cellMinDistance, cellMinDistanceGlobal, numberOfInterpolateCells, PH_MIN);
            for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
            {
                if (cellMinDistanceGlobal[iCell] < cellMinDistance[iCell])
                {
                    donorZoneIndex[iCell] = -1;
                }
            }

            PH_AllReduce(donorZoneIndex, donorZoneIndexGlobal, numberOfInterpolateCells, PH_MAX);
            for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
            {
                if (donorZoneIndexGlobal[iCell] != donorZoneIndex[iCell])
                {
                    donorCellIndex[iCell] = -1;
                }
            }

            PH_AllReduce(donorCellIndex, donorCellIndexGlobal, numberOfInterpolateCells, PH_MAX);

            SetField(donorLevel, donorLevelGlobal, numberOfInterpolateCells);
            SetField(cellMinDistance, cellMinDistanceGlobal, numberOfInterpolateCells);
            SetField(donorZoneIndex, donorZoneIndexGlobal, numberOfInterpolateCells);
            SetField(donorCellIndex, donorCellIndexGlobal, numberOfInterpolateCells);
;
            DelPointer(donorLevelGlobal);
            DelPointer(cellMinDistanceGlobal);
            DelPointer(donorZoneIndexGlobal);
            DelPointer(donorCellIndexGlobal);
        }
    }
}

void OversetConfigFactoryMulti::AllReduceZonesDefectCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *zoneInterpolateCells = zonesInterpolateCells[iZone];

        int numberOfInterpolateCells = zoneInterpolateCells->GetNumberOfInterCells();

        if (numberOfInterpolateCells)
        {
            int *isDefectCell =  zoneInterpolateCells->IsDefectCell();
            int *isDefectCellGlobal =  new int[numberOfInterpolateCells];

            PH_AllReduce(isDefectCell, isDefectCellGlobal, numberOfInterpolateCells, PH_MAX);

            SetField(isDefectCell, isDefectCellGlobal, numberOfInterpolateCells);

            DelPointer(isDefectCellGlobal);
        }
    }
}

void OversetConfigFactoryMulti::FindInterCellsByBoundary()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (processorIndex == currentProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->DefineInterCellsByBoundary();
            oversetConfigMultiZone[iZone]->FindInterpolateExpend();
        }
    }
}

void OversetConfigFactoryMulti::FindActiveCellsByBoundary()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (processorIndex == currentProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->DefineActCellsByBoundary();
        }
    }
}

void OversetConfigFactoryMulti::ResetCellCentersByBoundary()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (processorIndex == currentProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->ResetCellCentersByInterCell();
        }
    }
}
//! Identify active cells.
void OversetConfigFactoryMulti::FindActiveCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    //! Interpolation cells are determined.
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (processorIndex == currentProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->FindActiveCells();
        }
    }

    OutPutOversetVisualizations("2_1_FindActiveCells");
}

//! Identify interpolation cells.
void OversetConfigFactoryMulti::FindInterCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    //! Interpolation cells are determined.

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (processorIndex == currentProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->FindInterCells();
        }
    }
    OutPutOversetVisualizations("2_2_FindInterpolateCells");

}

void OversetConfigFactoryMulti::CheckInterpolateRelationships()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        int numberOfInterpolateCellsInThisBlockLocal = 0;

        int numberOfNegativeCellsInThisBlock_1Local = 0;
        int numberOfNegativeCellsInThisBlock_2Local = 0;
        int numberOfNegativeCellsInThisBlock_3Local = 0;

        for (int iZone = 0; iZone < numberOfZones; ++iZone)
        {
            int zoneProcessorIndex = PHMPI::GetZoneProcessorID(iZone);
            if (currentProcessorIndex == zoneProcessorIndex)
            {
                oversetConfigMultiZone[iZone]->CheckInterpolateRelationship(iBlock, zonesInterpolateCells[iZone],
                    numberOfInterpolateCellsInThisBlockLocal,
                    numberOfNegativeCellsInThisBlock_1Local,
                    numberOfNegativeCellsInThisBlock_2Local,
                    numberOfNegativeCellsInThisBlock_3Local);
            }
        }
    }
}

void OversetConfigFactoryMulti::ResetOversetInformation()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    //! Initialize oversetInformation on each zone.
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->InitializeOversetInformation();
        }
    }

    //! Get oversetInterCells and oversetDonorCells at MainProcess, where oversetDonorCells should swap with neighbor zone.
    GetInterAndDonorCellManagersNew();

    //! PostProcessOversetInformation, get storageIndexContainer (Interpolation or contribution cell for a neighbor zone, its serial number in the whole interpolation or contribution cell).
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->PostProcessOversetInformation();
        }
    }
}

void OversetConfigFactoryMulti::GetInterAndDonorCellManagers()
{

}

void OversetConfigFactoryMulti::GetInterAndDonorCellManagersNew()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int sendProcessorIndex = PHMPI::GetZoneProcessorID(iZone);

        //! First get neighborZoneIndexContainer and broadcast to all process.
        PHVectorInt1D neighborZoneIndexContainer;

        GetNeighborZoneIndexContainerForTargetZone(iZone, neighborZoneIndexContainer);

        int numberOfNeighbors = static_cast<int>( neighborZoneIndexContainer.size() );

        //! Loop query all neighbors one by one.
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++iNeighbor)
        {
            DataContainer *dataContainer = new DataContainer;

            int donorZoneIndex = neighborZoneIndexContainer[iNeighbor];
            int receiveProcessorIndex = PHMPI::GetZoneProcessorID(donorZoneIndex);

            int numberOfDonorCell;
            PHVectorInt1D     donorCellIndexContainer;
            PHVectorRDouble1D interCellCenterX;
            PHVectorRDouble1D interCellCenterY;
            PHVectorRDouble1D interCellCenterZ;
            PHVectorRDouble1D interCellVolume;
            if (currentProcessorIndex != receiveProcessorIndex && currentProcessorIndex != sendProcessorIndex)
            {
                continue;
            }

            if (currentProcessorIndex == sendProcessorIndex)
            {
                //! Get interCells ready.
                OversetInformationProxy *oversetInformationProxy = oversetConfigMultiZone[iZone]->GetOversetInformationProxy();
                OversetSwapManager *oversetSwapManager = oversetInformationProxy->GetOversetSwapManager();

                OversetInterplateCell *oversetInterplateCell = new OversetInterplateCell();
                oversetInterplateCell->SetZoneIndex(iZone);
                oversetInterplateCell->SetDonorZoneIndex(donorZoneIndex);
                PHVectorInt1D &cellIndexContainer = oversetInterplateCell->GetCellIndexContainer();
                oversetInformationProxy->GetCellIndexContainer(donorZoneIndex, cellIndexContainer);

                oversetSwapManager->AddOversetInterplateCell(oversetInterplateCell);

                //! Get necessary information ready for donorCells.
                oversetInformationProxy->GetDonorCellIndexContainer(donorZoneIndex, donorCellIndexContainer, interCellCenterX, interCellCenterY, interCellCenterZ, interCellVolume);

                numberOfDonorCell = static_cast<int>( donorCellIndexContainer.size() );

                PHSPACE::PHWrite(dataContainer, numberOfDonorCell);

                PHSPACE::PHWrite(dataContainer, donorCellIndexContainer, numberOfDonorCell);
                PHSPACE::PHWrite(dataContainer, interCellCenterX, numberOfDonorCell);
                PHSPACE::PHWrite(dataContainer, interCellCenterY, numberOfDonorCell);
                PHSPACE::PHWrite(dataContainer, interCellCenterZ, numberOfDonorCell);
                PHSPACE::PHWrite(dataContainer, interCellVolume, numberOfDonorCell);
            }

            PHMPI::PH_Trade(dataContainer, sendProcessorIndex, receiveProcessorIndex);

            if (currentProcessorIndex == receiveProcessorIndex)
            {
                dataContainer->MoveToBegin();

                PHSPACE::PHRead(dataContainer, numberOfDonorCell);

                donorCellIndexContainer.resize(numberOfDonorCell);
                interCellCenterX.resize(numberOfDonorCell);
                interCellCenterY.resize(numberOfDonorCell);
                interCellCenterZ.resize(numberOfDonorCell);
                interCellVolume.resize(numberOfDonorCell);

                PHSPACE::PHRead(dataContainer, donorCellIndexContainer, numberOfDonorCell);
                PHSPACE::PHRead(dataContainer, interCellCenterX, numberOfDonorCell);
                PHSPACE::PHRead(dataContainer, interCellCenterY, numberOfDonorCell);
                PHSPACE::PHRead(dataContainer, interCellCenterZ, numberOfDonorCell);
                PHSPACE::PHRead(dataContainer, interCellVolume, numberOfDonorCell);

                OversetInformationProxy *oversetInformationProxy = oversetConfigMultiZone[donorZoneIndex]->GetOversetInformationProxy();
                OversetSwapManager *oversetSwapManager = oversetInformationProxy->GetOversetSwapManager();

                OversetDonorCell *oversetDonorCell = new OversetDonorCell();
                oversetDonorCell->SetZoneIndex(donorZoneIndex);
                oversetDonorCell->SetInterplateZoneIndex(iZone);
                oversetDonorCell->SetDonorCellIndexContainer(donorCellIndexContainer);
                oversetDonorCell->SetInterCellCenterX(interCellCenterX);
                oversetDonorCell->SetInterCellCenterY(interCellCenterY);
                oversetDonorCell->SetInterCellCenterZ(interCellCenterZ);
                oversetDonorCell->SetInterCellVolume(interCellVolume);

                oversetSwapManager->AddOversetDonorCell(oversetDonorCell);
            }

            FreePointer(dataContainer);
        }
    }
}

void OversetConfigFactoryMulti::GetNeighborZoneIndexContainerForTargetZone(int iZone, PHVectorInt1D &neighborZoneList)
{
    InterCells *zoneInterpolateCells = zonesInterpolateCells[iZone];
    int numberOfInterCells = zoneInterpolateCells->GetNumberOfInterCells();

    int *donorZone = zoneInterpolateCells->GetDonorZoneIndex();

    set< int > neighborZoneIndexSet;

    for (int iCell = 0; iCell < numberOfInterCells; ++iCell)
    {
        neighborZoneIndexSet.insert(donorZone[iCell]);
    }

    for (set< int >::iterator iter = neighborZoneIndexSet.begin(); iter != neighborZoneIndexSet.end(); ++iter)
    {
        neighborZoneList.push_back(*iter);
    }
}

void OversetConfigFactoryMulti::OutPutOversetVisualizations(const string &name)
{
    if (GetIfOutPutOversetVisualization() == 0)
    {
        return;
    }

    string fileDir = "OversetOutput/" + name;
    PHMakeDir(fileDir);
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->OutPutOversetVisualization(name);
        }
    }
}

void OversetConfigFactoryMulti::FreeKeyBasicActiveNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->FreeKeyBasicActiveNode();
        }
    }
}

void OversetConfigFactoryMulti::OutPutOversetFiles()
{
    OutPutBoundaryCells();

    PrintCellInformation();

    OutPutTecplot();

    OutPutOversetInformation();
}

void OversetConfigFactoryMulti::OutPutOversetInformation()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int serverProcessor = PHMPI::GetServerSepMode();

    string oversetGridFileName = GlobalDataBase::GetStrParaFromDB("oversetGridFileName");
    fstream file;
    if (currentProcessor == serverProcessor)
    {
        PrintToWindow("Output Overset Configuration Information: ", oversetGridFileName, "\n");
        OpenFile(file, oversetGridFileName, ios_base::out | ios_base::binary);
    }

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        DataContainer *dataContainer = new DataContainer;

        int sendProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessor == sendProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->GetGrid()->Encode(dataContainer, 1);
        }

        PHMPI::PH_Trade(dataContainer, sendProcessorIndex, serverProcessor);

        if (currentProcessor == serverProcessor)
        {
            dataContainer->MoveToBegin();
            dataContainer->WriteFile(file);
        }
        delete dataContainer;    dataContainer = nullptr;
    }

    if (currentProcessor == serverProcessor)
    {
        CloseFile(file);
    }
}

void OversetConfigFactoryMulti::OutPutTecplot()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->ShowGrid();
        }
    }
}

void OversetConfigFactoryMulti::OutPutBoundaryCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->ShowBoundaryCells();
        }
    }
}

void OversetConfigFactoryMulti::PrintCellInformation()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    int serverProcessorIndex = PHMPI::GetServerSepMode();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int numberOfActiveCells = 0;
        int numberOfNegativeCells = 0;
        int numberOfInterCells = 0;

        int sendProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        DataContainer *dataContainer = new DataContainer;
        int tag = 0;

        if (currentProcessorIndex == sendProcessorIndex)
        {
            int *keyActiveOfCells = oversetConfigMultiZone[iZone]->GetKeyActiveOfCells();

            int numberOfCells = ( oversetConfigMultiZone[iZone]->grid )->GetNTotalCell();

            for (int iCell = 0; iCell < numberOfCells; ++iCell)
            {
                if (keyActiveOfCells[iCell] == ACTIVE)
                {
                    numberOfActiveCells++;
                }
                else if (keyActiveOfCells[iCell] == INACTIVE)
                {
                    numberOfNegativeCells++;
                }
                else if (keyActiveOfCells[iCell] == INTERPOLATION)
                {
                    numberOfInterCells++;
                }
            }

            if (currentProcessorIndex != serverProcessorIndex)
            {
                PHWrite(dataContainer, numberOfActiveCells);
                PHWrite(dataContainer, numberOfNegativeCells);
                PHWrite(dataContainer, numberOfInterCells);
                send(dataContainer, serverProcessorIndex, tag);
            }
        }

        if (currentProcessorIndex == serverProcessorIndex)
        {
            if (currentProcessorIndex != sendProcessorIndex)
            {
                receive(dataContainer, sendProcessorIndex, tag);
                PHRead(dataContainer, numberOfActiveCells);
                PHRead(dataContainer, numberOfNegativeCells);
                PHRead(dataContainer, numberOfInterCells);
            }

            cout << "iZone = " << iZone << endl;
            cout << "       numberOfActiveCells   = " << numberOfActiveCells << endl;
            cout << "       numberOfNegativeCells = " << numberOfNegativeCells << endl;
            cout << "       numberOfInterCells    = " << numberOfInterCells << endl;
        }

        FreePointer(dataContainer);
    }
}

void OversetConfigFactoryMulti::SetCellLayersByActiveNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetCellLayersByActiveNodes();
        }
    }
}

void OversetConfigFactoryMulti::CharacterNodeToCell()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->CharacterNodeToCell();
        }
    }
}

//! For physical boundary advance, advanced to the interpolation boundary and stop, advanced to wallCell and stop too.
int OversetConfigFactoryMulti::CharacterNodeToValidNodeByCell()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    int numberOfNewNodesAllZones = 0;

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            int numberOfNewNodes = oversetConfigMultiZone[iZone]->CharacterNodeToValidNodeByCell();
            numberOfNewNodesAllZones += numberOfNewNodes;
        }
    }

    //! Transmission of point attributes through interface communication. InnerAndInitialInterNodes needs to be constrained.
    CommunicateNodeCharacterOnInterface();

    SetNegativeForInnerAndInitialInterNodes();

    return numberOfNewNodesAllZones;
}

void OversetConfigFactoryMulti::CharacterCellToNode()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->CharacterCellToNode();
        }
    }

    //! Transmission of point attributes through interface communication.
    CommunicateNodeCharacterOnInterface();
}

void OversetConfigFactoryMulti::CommunicateCellIBlank()
{
    using namespace PHMPI;
    ActionKey *actkey = new ActionKey();
    actkey->action = COMMCELLIBLANK;
    actkey->kind = GRID_BASED;
    actkey->level = 0;

    if (IsNeedNonBlockingCommunication(actkey))
    {
        MainTaskNonBlocking(actkey);
    }
    else
    {
        MainTaskBlocking(actkey);
    }

    FreePointer(actkey);
}

//! Through the communication of virtual cell information in interface, realize the indirect communication of point attributes.
void OversetConfigFactoryMulti::CommunicateNodeCharacterOnInterface()
{
    if (numberOfBlocks == numberOfZones)
    {
        return;
    }

    int currentProcessor = PHMPI::GetCurrentProcessorID();

    //! step1: Communication cell properties.
    CommunicateCellIBlank();

    //! step2: Judge the property of points through Interface.
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->CharacterInterfaceToNode();
        }
    }
}

void OversetConfigFactoryMulti::CommunicateNodePropertyByKeyOfInterfaceNodesInBlocks()
{
    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        SetField(keyOfInterfaceNodesInBlocks[iBlock], 0, numberOfInterfaceNodesInBlocks[iBlock]);
    }

    //! step1: Every zone update keyActiveOfGlobalNodes (Summary keyActiveOfNodes).
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->UpdateInterfaceNodeCharacter();
        }
    }

    //! step2: mpiAllreduce.
    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        int &numberOfInterfaceNodesInThisBlock = numberOfInterfaceNodesInBlocks[iBlock];
        int *keyOfInterfaceNodesInThisBlockGlobal = new int[numberOfInterfaceNodesInThisBlock];

        PH_AllReduce(keyOfInterfaceNodesInBlocks[iBlock], keyOfInterfaceNodesInThisBlockGlobal, numberOfInterfaceNodesInThisBlock, MPI_MAX);
        SetField(keyOfInterfaceNodesInBlocks[iBlock], keyOfInterfaceNodesInThisBlockGlobal, numberOfInterfaceNodesInThisBlock);

        delete [] keyOfInterfaceNodesInThisBlockGlobal;
        keyOfInterfaceNodesInThisBlockGlobal = NULL;
    }

    //! step3: Assign the value of after mpiAllreduce to keyActiveOfNodes.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->DownInterfaceNodeCharacter();
        }
    }
}

void OversetConfigFactoryMulti::OutInterCellInformation()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    string fileName = AddSymbolToFileName("interCells.dat", currentProcessorIndex);
    fstream file;
    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *zoneInterpolateCells = zonesInterpolateCells[iZone];
        int numberOfInterCells = zoneInterpolateCells->GetNumberOfInterCells();

        RDouble *minDist = zoneInterpolateCells->GetMinDist();
        int *donorZone = zoneInterpolateCells->GetDonorZoneIndex();
        int *donorCell = zoneInterpolateCells->GetDonorCellIndex();

        file << "iZone = " << iZone << endl;
        for (int iCell = 0; iCell < numberOfInterCells; ++iCell)
        {
            file << "iCell = " << iCell << ":" << endl;
            file << "        " << minDist[iCell] << endl;
            file << "        " << donorZone[iCell] << endl;
            file << "        " << donorCell[iCell] << endl;
        }
    }
    PHSPACE::CloseFile(file);
}

void OversetConfigFactoryMulti::OutKeyActive()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->OutKeyActive();
        }
    }
}

void OversetConfigFactoryMulti::OutGlobalInterCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    string  fileName = AddSymbolToFileName("globalInterCell.dat", currentProcessorIndex);
    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *zoneInterpolateCells = zonesInterpolateCells[iZone];
        int blockIndex = zoneInterpolateCells->GetBlockIndex();
        int numberOfInterCells = zoneInterpolateCells->GetNumberOfInterCells();

        RDouble *x = zoneInterpolateCells->GetX();
        RDouble *y = zoneInterpolateCells->GetX();
        RDouble *z = zoneInterpolateCells->GetX();

        file << "blockIndex         = " << blockIndex << endl;
        file << "numberOfInterCells = " << numberOfInterCells << endl;

        for (int iCell = 0; iCell < numberOfInterCells; ++iCell)
        {
            file << x[iCell] << "  " << y[iCell] << "  " << z[iCell] << endl;
        }
    }

    PHSPACE::CloseFile(file);
}

void OversetConfigFactoryMulti::OutGlobalNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    string fileName = AddSymbolToFileName("globalNodes.dat", currentProcessorIndex);
    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        ZoneSearchNodes *zoneSearchNodes = zonesSearchNodes[iZone];
        int blockIndex = zoneSearchNodes->GetBlockIndex();
        int numberOfNodes = zoneSearchNodes->GetNumberOfNodes();

        RDouble *x = zoneSearchNodes->GetX();
        RDouble *y = zoneSearchNodes->GetX();
        RDouble *z = zoneSearchNodes->GetX();
        RDouble *pMin = zoneSearchNodes->GetPMin();
        RDouble *pMax = zoneSearchNodes->GetPMax();

        file << "blockIndex    = " << blockIndex << endl;
        file << "numberOfNodes = " << numberOfNodes << endl;

        file << pMin[0] << "  " << pMin[1] << "  " << pMin[2] << endl;

        file << pMax[0] << "  " << pMax[1] << "  " << pMax[2] << endl;

        for (int iNode = 0; iNode < numberOfNodes; ++iNode)
        {
            file << x[iNode] << "  " << y[iNode] << "  " << z[iNode] << endl;
        }
    }

    PHSPACE::CloseFile(file);
}

ZoneSearchLines::ZoneSearchLines(int numberOfSearchNodesIn, int numberOfSearchLinesIn,
    int blockIndexIn):
    numberOfNodes(numberOfSearchNodesIn),
    numberOfLines(numberOfSearchLinesIn),
    blockIndex(blockIndexIn)
{
    int geometricDimension = GetDim();
    pMin = new RDouble[geometricDimension];
    pMax = new RDouble[geometricDimension];

    if (numberOfNodes)
    {
        x = new RDouble[numberOfNodes];
        y = new RDouble[numberOfNodes];
        z = new RDouble[numberOfNodes];

        line2Node = new int[numberOfNodes];

        lineIndex = new int[numberOfLines];

        nodesActive =  new int[numberOfNodes];

        SetField(nodesActive, VARIABLEACTIVE, numberOfNodes);

        linesActive =  new int[numberOfLines];

        SetField(linesActive, VARIABLEACTIVE, numberOfLines);
    }
}
ZoneSearchLines ::~ZoneSearchLines()
{
    if (numberOfNodes)
    {
        DelPointer(x);
        DelPointer(y);
        DelPointer(z);

        DelPointer(line2Node);
        DelPointer(lineIndex);

        DelPointer(nodesActive);
        DelPointer(linesActive);

    }

    DelPointer(pMin);
    DelPointer(pMax);
}

ZoneSearchNodes::ZoneSearchNodes(int numberOfSearchNodesIn,
    int blockIndexIn) :
    numberOfNodes(numberOfSearchNodesIn),
    blockIndex(blockIndexIn)
{
    int geometricDimension = GetDim();
    pMin = new RDouble[geometricDimension];
    pMax = new RDouble[geometricDimension];

    if (numberOfNodes)
    {
        x = new RDouble[numberOfNodes];
        y = new RDouble[numberOfNodes];
        z = new RDouble[numberOfNodes];

        nodesMinDistanceToOtherBlock = new RDouble[numberOfNodes];

        SetField(nodesMinDistanceToOtherBlock, LARGE, numberOfNodes);
    }
}

ZoneSearchNodes ::~ZoneSearchNodes()
{
    if (numberOfNodes)
    {
        DelPointer(x);
        DelPointer(y);
        DelPointer(z);
        DelPointer(nodesMinDistanceToOtherBlock);
    }

    DelPointer(pMin);
    DelPointer(pMax);
}

InterCells::InterCells(int numberOfInterpolateCellsIn,
    int blockIndexIn) :
    numberOfInterCells(numberOfInterpolateCellsIn),
    blockIndex(blockIndexIn)
{
    if (numberOfInterCells)
    {
        donorZoneIndex = new int[numberOfInterCells];
        donorCellIndex = new int[numberOfInterCells];
        isDefectCell   = new int[numberOfInterCells];

        donorLevel = new int[numberOfInterCells];
        cellMinDistance = new RDouble[numberOfInterCells];

        x = new RDouble[numberOfInterCells];
        y = new RDouble[numberOfInterCells];
        z = new RDouble[numberOfInterCells];

        SetField(donorZoneIndex, -1, numberOfInterCells);
        SetField(donorCellIndex, -1, numberOfInterCells);
        SetField(isDefectCell,    0, numberOfInterCells);

        SetField(donorLevel, -3, numberOfInterCells);
        SetField(cellMinDistance, LARGE, numberOfInterCells);
    }
}

InterCells::~InterCells()
{
    if (numberOfInterCells)
    {
        DelPointer(donorZoneIndex);
        DelPointer(donorCellIndex);
        DelPointer(isDefectCell);

        if (nullptr != donorLevel)
        {
            DelPointer(donorLevel);
            DelPointer(cellMinDistance);

            DelPointer(x);
            DelPointer(y);
            DelPointer(z);
        }
    }
}

void InterCells::FreeUselessMemery()
{
    if (numberOfInterCells)
    {
        DelPointer(donorLevel);
        DelPointer(cellMinDistance);

        DelPointer(x);
        DelPointer(y);
        DelPointer(z);
    }
}

OversetConfig::OversetConfig()
{
    gridKDTree = nullptr;
    cellKDTree = nullptr;

    minBoxOfZone = nullptr;
    maxBoxOfZone = nullptr;
    minBoxOfZoneByLayer = nullptr;
    maxBoxOfZoneByLayer = nullptr;
    minBoxOfLine = nullptr;
    maxBoxOfLine = nullptr;
    minBoxNodeID = nullptr;
    maxBoxNodeID = nullptr;
    hasCalMinMaxForZone = false;

    keyActiveOfNodes = nullptr;
    keyActiveOfLines = nullptr;
    keyActiveOfFaces = nullptr;

    keySearchOfNodes = nullptr;
    keySearchOfCells = nullptr;
    keySearchOfLines = nullptr;
    keySearchOfNodesByLayer = nullptr;
    keySearchOfCellsByLayer = nullptr;

    keyBufferOfNodes = nullptr;
    keyBufferOfCells = nullptr;
    keySearchOfNodes = nullptr;
    keySearchOfCells = nullptr;
    keyInBlockNodes = nullptr;
    keyBasicActiveNodes = nullptr;

    nodesMinDistanceInThisBlock = nullptr;
    nodesMinDistanceInOtherBlock = nullptr;

    localToGlobalInterfaceNodesMap = nullptr;

    zoneIndex = 0;
    donorCellMethod = 0;
    nLayerOfCellsOfCutting = 0;
}

void OversetConfig::BuildBoundaryNodeList(vector< RDouble > &wallNodeCoordinateXIn,
    vector< RDouble > &wallNodeCoordinateYIn,
    vector< RDouble > &wallNodeCoordinateZIn)
{
    int numberOfBoundaryElements = grid->GetNBoundFace();
    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int **faceNodeIndex = grid->GetFace2NodeArray();

    UnstructBCSet **unstructBCSet = grid->GetBCRecord();

    RDouble *nodeCoordinateX = grid->GetX();
    RDouble *nodeCoordinateY = grid->GetY();
    RDouble *nodeCoordinateZ = grid->GetZ();

    set< int > interfaceBoundaryNodeSet;
    set< int > wallBoundaryNodeSet;
    set< int > oversetBoundaryNodeSet;
    set< int > symmetryBoundaryNodeSet;
    set< int > physicalBoundaryNodeSet;

    for (int iFace = 0; iFace < numberOfBoundaryElements; ++iFace)
    {
        int boundaryConditionOfThisFace = unstructBCSet[iFace]->GetKey();
        int &numberOfNodesInThisFace = faceNodeNumber[iFace];

        switch (boundaryConditionOfThisFace)
        {
            case PHENGLEI::INTERFACE:
                for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                {
                    interfaceBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            case PHENGLEI::SOLID_SURFACE:
                for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                {
                    wallBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            case PHENGLEI::OVERSET:
                for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                {
                    oversetBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            case PHENGLEI::SYMMETRY:
                for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                {
                    symmetryBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            default:
                for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                {
                    physicalBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
        }
    }

    for (set< int >::iterator iter = interfaceBoundaryNodeSet.begin(); iter != interfaceBoundaryNodeSet.end(); ++iter)
    {
        interfaceBoundaryNodeList.push_back(*iter);
    }
    for (set< int >::iterator iter = oversetBoundaryNodeSet.begin(); iter != oversetBoundaryNodeSet.end(); ++iter)
    {
        innerBoundaryNodeList.push_back(*iter);
    }
    for (set< int >::iterator iter = physicalBoundaryNodeSet.begin(); iter != physicalBoundaryNodeSet.end(); ++iter)
    {
        physicalBoundaryNodeList.push_back(*iter);
    }

    for (set< int >::iterator iter = wallBoundaryNodeSet.begin(); iter != wallBoundaryNodeSet.end(); ++iter)
    {
        physicalBoundaryNodeList.push_back(*iter);
        wallNodeCoordinateXIn.push_back(nodeCoordinateX[*iter]);
        wallNodeCoordinateYIn.push_back(nodeCoordinateY[*iter]);
        wallNodeCoordinateZIn.push_back(nodeCoordinateZ[*iter]);
    }

    if (GetSymetryOrNot())
    {
        for (set< int >::iterator iter = symmetryBoundaryNodeSet.begin(); iter != symmetryBoundaryNodeSet.end(); ++iter)
        {
            physicalBoundaryNodeList.push_back(*iter);
        }
    }
}

void OversetConfig::BuildBoundaryNodeList()
{
    int numberOfBoundaryElements = grid->GetNBoundFace();
    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int **faceNodeIndex = grid->GetFace2NodeArray();

    UnstructBCSet **bcr = grid->GetBCRecord();

    set< int > interfaceBoundaryNodeSet;
    set< int > wallBoundaryNodeSet;
    set< int > oversetBoundaryNodeSet;
    set< int > symmetryBoundaryNodeSet;
    set< int > physicalBoundaryNodeSet;

    for (int iFace = 0; iFace < numberOfBoundaryElements; ++iFace)
    {
        int boundaryConditionOfThisFace = bcr[iFace]->GetKey();
        int &numberOfNodesInThisFace = faceNodeNumber[iFace];

        switch (boundaryConditionOfThisFace)
        {
            case PHENGLEI::INTERFACE:
                for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                {
                    interfaceBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            case PHENGLEI::SOLID_SURFACE:
                for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                {
                    wallBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            case PHENGLEI::OVERSET:
                for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                {
                    oversetBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            case PHENGLEI::SYMMETRY:
                for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                {
                    symmetryBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            default:
                for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                {
                    physicalBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
        }
    }

    for (set< int >::iterator iter = interfaceBoundaryNodeSet.begin(); iter != interfaceBoundaryNodeSet.end(); ++iter)
    {
        interfaceBoundaryNodeList.push_back(*iter);
    }
    for (set< int >::iterator iter = oversetBoundaryNodeSet.begin(); iter != oversetBoundaryNodeSet.end(); ++iter)
    {
        innerBoundaryNodeList.push_back(*iter);
    }
    for (set< int >::iterator iter = physicalBoundaryNodeSet.begin(); iter != physicalBoundaryNodeSet.end(); ++iter)
    {
        physicalBoundaryNodeList.push_back(*iter);
    }

    for (set< int >::iterator iter = wallBoundaryNodeSet.begin(); iter != wallBoundaryNodeSet.end(); ++iter)
    {
        physicalBoundaryNodeList.push_back(*iter);
    }

    if (GetSymetryOrNot())
    {
        for (set< int >::iterator iter = symmetryBoundaryNodeSet.begin(); iter != symmetryBoundaryNodeSet.end(); ++iter)
        {
            physicalBoundaryNodeList.push_back(*iter);
        }
    }
    if (nullptr == localToGlobalInterfaceNodesMap)
    {
        localToGlobalInterfaceNodesMap = new int[interfaceBoundaryNodeList.size()];
    }

}

void OversetConfig::BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTree(vector< GridKDTree * > &innerAuxiliaryGridKDTreesIn)
{
    RDouble *nodeCoordinate = new RDouble[THREE_D];

    int blockIndex = grid->GetIBlock();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfBlocks = GetOversetConfigFactoryMulti()->GetNumberOfBlocks();

    RDouble *cellCenterX = grid->GetCellCenterX();
    RDouble *cellCenterY = grid->GetCellCenterY();
    RDouble *cellCenterZ = grid->GetCellCenterZ();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    SetField(keyInBlockNodes, 0, numberOfNodes);

    int donorCellIndex;
    int symetryOrNot = GetSymetryOrNot();
    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        donorCellIndex = -2;

        if (!symetryOrNot || cellCenterZ[iCell] > 1.0e-5)
        {
            donorCellIndex = -1;

            nodeCoordinate[0] = cellCenterX[iCell];
            nodeCoordinate[1] = cellCenterY[iCell];
            nodeCoordinate[2] = cellCenterZ[iCell];

            for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
            {
                if (blockIndex != iBlock)
                {
                    if (innerAuxiliaryGridKDTreesIn[iBlock] != NULL)
                    {
                        donorCellIndex = innerAuxiliaryGridKDTreesIn[iBlock]->SearchNodeInKDTree(nodeCoordinate);
                    }
                }
            }
        }

        if (donorCellIndex != -1)
        {
            for (int iNode = 0; iNode < cellNodeNumberContainer[iCell]; ++iNode)
            {
                keyInBlockNodes[cellNodeIndexContainer[iCell][iNode]] = 1;
            }
        }
    }

    DelPointer(nodeCoordinate);
}

OversetConfig::~OversetConfig()
{
    DelPointer(keyActiveOfNodes);
    DelPointer(keyActiveOfLines);
    DelPointer(keyActiveOfFaces);
    DelPointer(keyBufferOfNodes);
    DelPointer(keyBufferOfCells);
    DelPointer(keySearchOfNodes);
    DelPointer(keySearchOfLines);

    FreeKeyInBlockNodes();

    if (nodesMinDistanceInThisBlock != nullptr)
    {
        DelPointer(nodesMinDistanceInThisBlock);
    }

    if (localToGlobalInterfaceNodesMap != nullptr)
    {
        DelPointer(localToGlobalInterfaceNodesMap);
    }

    DelPointer(nLayerOfCellsFromWall);
    DelPointer(keySearchOfCells);
    DelPointer(minBoxOfZoneByLayer);
    DelPointer(maxBoxOfZoneByLayer);
    DelPointer(minBoxOfZone);
    DelPointer(maxBoxOfZone);
    DelPointer(minBoxNodeID);
    DelPointer(maxBoxNodeID);
}

void OversetConfig::FreeKeyBasicActiveNode()
{
    DelPointer(keyBasicActiveNodes);
}

void OversetConfig::FreeKeyInBlockNodes()
{
    DelPointer(keyInBlockNodes);
}

void OversetConfig::FreeKDTree()
{
    if (cellKDTree)
    {
        cellKDTree->Free();
        FreePointer(cellKDTree);
    }
}

void OversetConfig::Initialize(UnstructGrid *gridIn,
    int nZone, int iZone, int iBlock)
{
    geometricDimension = GetDim();
    numberOfZones = nZone;
    zoneIndex = iZone;
    blockIndex = iBlock;
    grid = gridIn;
    oversetInformationProxy = nullptr;

    donorCellMethod = 0;
    if (GlobalDataBase::IsExist("donorCellMethod", PHINT, 1))
    {
        GlobalDataBase::GetData("donorCellMethod", &donorCellMethod, PHINT, 1);
    }
    nLayerOfCellsOfCutting = PHSPACE::GlobalDataBase::GetIntParaFromDB("nLayerOfCellsOfCutting");

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    keyActiveOfNodes = new int[numberOfNodes];
    SetField(keyActiveOfNodes, 0, numberOfNodes);

    nodesMinDistanceInThisBlock = new RDouble[numberOfNodes];
    SetField(nodesMinDistanceInThisBlock, LARGE, numberOfNodes);

    keyInBlockNodes = new int[numberOfNodes];

    nLayerOfCellsFromWall = new int[numberOfCells];
    SetField(nLayerOfCellsFromWall, MAXCELLLAYER, numberOfCells);

    minBoxOfZone = new RDouble[geometricDimension];
    maxBoxOfZone = new RDouble[geometricDimension];

    minBoxNodeID = new int[geometricDimension];
    maxBoxNodeID = new int[geometricDimension];
}

void OversetConfig::SetInitialValue()
{
    vector<int>().swap(interpolateCellList);
    vector<int>().swap(interpolateFaceList);
    vector<int>().swap(defectCellList);

    interpolateExpendList.clear();
    donorFaceList.clear();
}

//! Calculate the minimum box according to the cell.
void OversetConfig::BuildMinMaxBoxOfZoneBySearchRegion()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    vector< RDouble * > nodeCoordinate;
    RDouble *x = grid->GetX();
    nodeCoordinate.push_back(x);
    RDouble *y = grid->GetY();
    nodeCoordinate.push_back(y);
    RDouble *z = grid->GetZ();
    nodeCoordinate.push_back(z);

    SetField(minBoxOfZone, LARGE, geometricDimension);
    SetField(maxBoxOfZone, -LARGE, geometricDimension);

    RDouble eps = GetToleranceForOversetBox();
    for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
    {
        RDouble &minInThisDimension = minBoxOfZone[iDimension];
        RDouble &maxInThisDimension = maxBoxOfZone[iDimension];
        for (int iCell = 0; iCell < numberOfCells; ++iCell)
        {
            if (keySearchOfCells[iCell])
            {
                int cellNodeNumber = cellNodeNumberContainer[iCell];

                for (int iNode = 0; iNode < cellNodeNumber; ++iNode)
                {
                    int &nodeIndex = cellNodeIndexContainer[iCell][iNode];

                    RDouble &nodeCoordinateInThisDimension = nodeCoordinate[iDimension][nodeIndex];
                    minInThisDimension = MIN(minInThisDimension, nodeCoordinateInThisDimension);
                    maxInThisDimension = MAX(maxInThisDimension, nodeCoordinateInThisDimension);
                }
            }
        }
        minInThisDimension -= eps;
        maxInThisDimension += eps;
    }
}

void OversetConfig::SetInitialValueForBufferRigion()
{
    //! Because communicate_cell_iblank is the direct communication keyActiveOfCells.
    //! keyBuffer is determined by keyActiveOfCells array.
    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    SetField(keyActiveOfNodes, 0, numberOfNodes);
    SetField(keyActiveOfCells, 0, numberOfCells);

    int numberOfInnerBoundaryNodes = static_cast<int>(innerBoundaryNodeList.size());
    for (int iNode = 0; iNode < numberOfInnerBoundaryNodes; ++iNode)
    {
        keyActiveOfNodes[innerBoundaryNodeList[iNode]] = ACTIVE;
    }
}

void OversetConfig::BuildNodesMinDistanceInThisBlock(vector< RDouble > *wallNodeCoordinateXInBlocksIn,
    vector< RDouble > *wallNodeCoordinateYInBlocksIn,
    vector< RDouble > *wallNodeCoordinateZInBlocksIn)
{
    int blockIndex = grid->GetIBlock();
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfWallNodesInBlock = static_cast<int>(wallNodeCoordinateXInBlocksIn[blockIndex].size());

    nodesMinDistanceInThisBlock = new RDouble[numberOfNodes];

    if (!numberOfWallNodesInBlock)
    {
        double wallDistanceOfBackground = GlobalDataBase::GetDoubleParaFromDB("walldistMainZone");
        SetField(nodesMinDistanceInThisBlock, wallDistanceOfBackground, numberOfNodes);
        return;
    }

    RDouble *nodeCoordinateX = grid->GetX();
    RDouble *nodeCoordinateY = grid->GetY();
    RDouble *nodeCoordinateZ = grid->GetZ();

    SetField(nodesMinDistanceInThisBlock, LARGE, numberOfNodes);
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        for (int iWallNode = 0; iWallNode < numberOfWallNodesInBlock; ++iWallNode)
        {
            RDouble dx = nodeCoordinateX[iNode] - wallNodeCoordinateXInBlocksIn[blockIndex][iWallNode];
            RDouble dy = nodeCoordinateY[iNode] - wallNodeCoordinateYInBlocksIn[blockIndex][iWallNode];
            RDouble dz = nodeCoordinateZ[iNode] - wallNodeCoordinateZInBlocksIn[blockIndex][iWallNode];
            RDouble tempDistance = DISTANCE(dx, dy, dz);

            nodesMinDistanceInThisBlock[iNode] = MIN(tempDistance, nodesMinDistanceInThisBlock[iNode]);
        }
    }
}

void OversetConfig::SetNodesMinDistanceInThisBlocksByWalldist()
{
    int numberOfNodes = grid->GetNTotalNode();
    vector<RDouble> &nodeWalldist = grid->GetNodeWallDist();

    SetField(nodesMinDistanceInThisBlock, &nodeWalldist[0], numberOfNodes);
}

void OversetConfig::SetNodesMinDistanceInOtherBlocks()
{
    int numberOfNodes = grid->GetNTotalNode();
    vector<RDouble> &nodeWalldistInOtherBlocks = grid->GetNodeWallDistInOtherBlock();
    nodesMinDistanceInOtherBlock = new RDouble[numberOfNodes];
    SetField(nodesMinDistanceInOtherBlock, &nodeWalldistInOtherBlocks[0], numberOfNodes);
}

void OversetConfig::SetKeyBufferByKeyActive()
{
    //! Because communicate_cell_iblank is the direct communication keyActiveOfCells.
    //! keyBuffer is determined by keyActiveOfCells array.
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    keyBufferOfNodes = new int[numberOfNodes];
    keyBufferOfCells = new int[numberOfCells];

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    SetField(keyBufferOfCells, keyActiveOfCells, numberOfCells);
    SetField(keyBufferOfNodes, keyActiveOfNodes, numberOfNodes);
}

void OversetConfig::ComputeNodeAndCellTopology()
{
    grid->GetCell2NodeArray();
    grid->GetCell2Face();
}

void OversetConfig::CopyMinDist(RDouble *nodesMinDistanceInOtherBlockIn)
{
    int numberOfNodes = grid->GetNTotalNode();
    nodesMinDistanceInOtherBlock = new RDouble[numberOfNodes];
    SetField(nodesMinDistanceInOtherBlock, LARGE, numberOfNodes);

    int iSearchNode = 0;

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keySearchOfNodes[iNode])
        {
            nodesMinDistanceInOtherBlock[iNode] = nodesMinDistanceInOtherBlockIn[iSearchNode++];
        }
    }

    DelPointer(keySearchOfNodes);
}

void OversetConfig::CopyNodesActive(int* nodesActiveIn, int* nodes2NodeIn, int numberOfLinesIn)
{
    int numberOfNodes = grid->GetNTotalNode();
    keyActiveOfNodes = new int[numberOfNodes];
    SetField(keyActiveOfNodes, VARIABLEACTIVE, numberOfNodes);

    int countIndex = 0;
    int countActive = 0;
    for (int iSearchLine = 0; iSearchLine < numberOfLinesIn; iSearchLine++)
    {
        int nodeIndex0 = nodes2NodeIn[countIndex++];
        int nodeIndex1 = nodes2NodeIn[countIndex++];

        keyActiveOfNodes[nodeIndex0] = nodesActiveIn[countActive++];
        keyActiveOfNodes[nodeIndex1] = nodesActiveIn[countActive++]; 
    }

    DelPointer(keyActiveOfNodes);
}

void OversetConfig::BuildMinMaxBoxOfZone()
{
    vector< RDouble * > nodeCoordinate;
    RDouble *x = grid->GetX();
    nodeCoordinate.push_back(x);
    RDouble *y = grid->GetY();
    nodeCoordinate.push_back(y);
    RDouble *z = grid->GetZ();
    nodeCoordinate.push_back(z);

    if (hasCalMinMaxForZone == false)
    {
    SetField(minBoxOfZone, LARGE, geometricDimension);
    SetField(maxBoxOfZone, -LARGE, geometricDimension);

        SetField(minBoxNodeID, 0, geometricDimension);
        SetField(maxBoxNodeID, 0, geometricDimension);
    }
    else
    {
        for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
        {
            minBoxOfZone[iDimension] = nodeCoordinate[iDimension][minBoxNodeID[iDimension]];
            maxBoxOfZone[iDimension] = nodeCoordinate[iDimension][maxBoxNodeID[iDimension]];
        }
        return;
    }

    int numberOfNodes = grid->GetNTotalNode();
    for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
    {
        RDouble &minInThisDimension = minBoxOfZone[iDimension];
        RDouble &maxInThisDimension = maxBoxOfZone[iDimension];
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            RDouble &nodeCoordinateInThisDimension = nodeCoordinate[iDimension][iNode];
            if (minInThisDimension >= nodeCoordinateInThisDimension)
            {
                 minBoxNodeID[iDimension] = iNode;
    }
            if (maxInThisDimension <= nodeCoordinateInThisDimension)
            {
                 maxBoxNodeID[iDimension] = iNode;
            }
            minInThisDimension = MIN(minInThisDimension, nodeCoordinateInThisDimension);
            maxInThisDimension = MAX(maxInThisDimension, nodeCoordinateInThisDimension);
        }
    }
    hasCalMinMaxForZone = true;
}

void OversetConfig::BuildKDTreeByMinMaxBoxOfZone()
{
    gridKDTree = new GridKDTree(grid, keySearchOfCells, minBoxOfZone, maxBoxOfZone);
}

void OversetConfig::BuildKDTreeByActiveCell()
{
    int nDim = grid->GetDim();

    int numberOfCells = grid->GetNTotalCell();
    int nBoundaryFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int *keyActiveOfCells = grid->GetBlankIndex();

    bool hasActiveCells = false, hasInterCells = false;
    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        if (keyActiveOfCells[iCell] == ACTIVE)
        {
            hasActiveCells = true;
            break;
        }
    }

    RDouble *nodeCoordinate = new RDouble[nDim];
    RDouble *cellCenterX = grid->GetCellCenterX();
    RDouble *cellCenterY = grid->GetCellCenterY();
    RDouble *cellCenterZ = grid->GetCellCenterZ();

    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");
    if (isOversetSlip)
    {
        return;
    }
    else 
    {
        if (hasActiveCells && !cellKDTree)
        {
            //!There is a memory leak in the cellKDTree, which needs to be addressed
            cellKDTree = CreatKDTree(nDim); 
        }

        for (int iCell = 0; iCell < numberOfCells; ++iCell)
        {
            if (keyActiveOfCells[iCell] == ACTIVE)
            {
                nodeCoordinate[0] = cellCenterX[iCell];
                nodeCoordinate[1] = cellCenterY[iCell];
                if (nDim == 3)nodeCoordinate[2] = cellCenterZ[iCell];
        
                KDInsert(cellKDTree, nodeCoordinate, iCell);
            }
        }
    }

    DelPointer(nodeCoordinate);
}

void OversetConfig::SearchZonesSearchNodesInLineBox()
{
    int blockIndex = grid->GetIBlock();
    int numberOfCells = grid->GetNTotalCell();
    int **cellFaceIndexContainer  = grid->GetCell2Face();
    int * cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    minBoxOfLine = new RDouble[geometricDimension];
    maxBoxOfLine = new RDouble[geometricDimension];

    RDouble *minBoxOfTemp = new RDouble[geometricDimension];
    RDouble *maxBoxOfTemp = new RDouble[geometricDimension];

    int numberOfNodes = grid->GetNTotalNode();
    int *keyActiveOfCells = GetKeyActiveOfCells();

    SetField(keyActiveOfNodes, 1, numberOfNodes);
    SetField(keyActiveOfCells, 0, numberOfCells);

    vector< RDouble * > nodeCoordinate;
    RDouble *x = grid->GetX();
    nodeCoordinate.push_back(x);
    RDouble *y = grid->GetY();
    nodeCoordinate.push_back(y);
    RDouble *z = grid->GetZ();
    nodeCoordinate.push_back(z);

    Vector3D s1P0, s1P1, s2P0, s2P1;
    RDouble eps = GetToleranceForOversetBox();
    //RDouble eps = 0.0;
    int **faceNodeIndex = grid->GetFace2NodeArray();
    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();

    ZoneSearchNodes **zonesSearchNodes = GetOversetConfigFactoryMulti()->GetZonesSearchNodes();
    OversetConfig **oversetConfigMultiZone = GetOversetConfigFactoryMulti()->GetOversetConfig();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        ZoneSearchNodes *&zoneSearchNodes = zonesSearchNodes[iZone];
        OversetConfig *&oversetConfig = oversetConfigMultiZone[iZone];

        UnstructGrid *otherGrid = oversetConfig->GetGrid(); 
        int numberOfBoundary = otherGrid->GetNBoundFace();
        UnstructBCSet **unstructBCSet = otherGrid->GetBCRecord();
        int *NodeNumber = otherGrid->GetNodeNumberOfEachFace();
        int **NodeIndex = otherGrid->GetFace2NodeArray();
        int *nLayerOfCellsFromWall = oversetConfig->GetLayerOfCells(); 
        int nTotalCells = otherGrid->GetNTotalCell();
        int **cellToFace = otherGrid->GetCell2Face();
        int *faceNumberOfEachCell = otherGrid->GetFaceNumberOfEachCell();

        int nodeID = 0;
        if (blockIndex != zoneSearchNodes->GetBlockIndex())
        {
            vector< RDouble * > otherNodeCoordinate;
            RDouble *x = otherGrid->GetX();
            otherNodeCoordinate.push_back(x);
            RDouble *y = otherGrid->GetY();
            otherNodeCoordinate.push_back(y);
            RDouble *z = otherGrid->GetZ();
            otherNodeCoordinate.push_back(z);

            for (int iCell = 0; iCell < numberOfCells; ++iCell)
            {
                if (keySearchOfCellsByLayer[iCell])
                {
                    for(int iFace = 0; iFace < cellFaceNumberContainer[ iCell ]; ++ iFace)
                    {
                        SetField(minBoxOfLine, LARGE, geometricDimension);
                        SetField(maxBoxOfLine, - LARGE, geometricDimension);

                        int faceIndex = cellFaceIndexContainer[iCell][iFace];
                        int nNode = faceNodeNumber[faceIndex];
                        for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
                        {
                            RDouble &minInThisDimension = minBoxOfLine[iDimension];
                            RDouble &maxInThisDimension = maxBoxOfLine[iDimension];
                            for (int iNode = 0; iNode < nNode; ++iNode)
                            {
                                nodeID = faceNodeIndex[faceIndex][iNode];
                                RDouble &nodeCoordinateInThisDimension = nodeCoordinate[iDimension][nodeID];
                                minInThisDimension = MIN(minInThisDimension, nodeCoordinateInThisDimension);
                                maxInThisDimension = MAX(maxInThisDimension, nodeCoordinateInThisDimension);
                            }
                            minInThisDimension -= eps;
                            maxInThisDimension += eps;
                        }

                        s1P0.x = nodeCoordinate[0][faceNodeIndex[faceIndex][0]];
                        s1P0.y = nodeCoordinate[1][faceNodeIndex[faceIndex][0]];
                        s1P1.x = nodeCoordinate[0][faceNodeIndex[faceIndex][1]];
                        s1P1.y = nodeCoordinate[1][faceNodeIndex[faceIndex][1]];

                        if (nLayerOfCellsOfCutting == 0)
                        {
                            for (int jFace = 0; jFace < numberOfBoundary; ++jFace)
                            {
                                int bcType = unstructBCSet[jFace]->GetKey();
                                if (bcType == PHENGLEI::SOLID_SURFACE)
                                {
                                    SetField(minBoxOfTemp, LARGE, geometricDimension);
                                    SetField(maxBoxOfTemp, - LARGE, geometricDimension);
                                    int &numberOfNodesInThisFace = NodeNumber[jFace];
                                    for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
                                    {
                                        RDouble &minInThisDimension = minBoxOfTemp[iDimension];
                                        RDouble &maxInThisDimension = maxBoxOfTemp[iDimension];
                                        for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                                        {
                                            RDouble &nodeCoordinateInThisDimension = otherNodeCoordinate[iDimension][NodeIndex[jFace][iNode]];
                                            minInThisDimension = MIN(minInThisDimension, nodeCoordinateInThisDimension);
                                            maxInThisDimension = MAX(maxInThisDimension, nodeCoordinateInThisDimension);
                                        }
                                        minInThisDimension -= eps;
                                        maxInThisDimension += eps;
                                    }
                                    if (IfWallFacesIntersectionWithThisZone(minBoxOfTemp,maxBoxOfTemp))
                                    {
                                        s2P0.x = otherNodeCoordinate[0][NodeIndex[jFace][0]];
                                        s2P0.y = otherNodeCoordinate[1][NodeIndex[jFace][0]];
                                        s2P1.x = otherNodeCoordinate[0][NodeIndex[jFace][1]];
                                        s2P1.y = otherNodeCoordinate[1][NodeIndex[jFace][1]];

                                        if (IfWallFaceIntersectionWithThisZone(s1P0 ,s1P1, s2P0, s2P1))
                                        {
                                            keyActiveOfNodes[faceNodeIndex[faceIndex][0]] = 0;
                                            keyActiveOfNodes[faceNodeIndex[faceIndex][1]] = 0;
                                            intersectLinesWithWall.insert(faceIndex);
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            for (int iCellInOtherGrid = 0; iCellInOtherGrid < nTotalCells; iCellInOtherGrid++)
                            {
                                if (nLayerOfCellsFromWall[iCellInOtherGrid] == nLayerOfCellsOfCutting)
                                {
                                    for (int iFace = 0; iFace < faceNumberOfEachCell[iCellInOtherGrid]; ++iFace)
                                    {
                                        SetField(minBoxOfTemp, LARGE, geometricDimension);
                                        SetField(maxBoxOfTemp, - LARGE, geometricDimension);

                                        int jFace = cellToFace[iCellInOtherGrid][iFace];
                                        int &numberOfNodesInThisFace = NodeNumber[jFace];

                                        for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
                                        {
                                            RDouble& minInThisDimension = minBoxOfTemp[iDimension];
                                            RDouble& maxInThisDimension = maxBoxOfTemp[iDimension];
                                            for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
                                            {
                                                RDouble& nodeCoordinateInThisDimension = otherNodeCoordinate[iDimension][NodeIndex[jFace][iNode]];
                                                minInThisDimension = MIN(minInThisDimension, nodeCoordinateInThisDimension);
                                                maxInThisDimension = MAX(maxInThisDimension, nodeCoordinateInThisDimension);
                                            }
                                            minInThisDimension -= eps;
                                            maxInThisDimension += eps;
                                        }
                                        if (IfWallFacesIntersectionWithThisZone(minBoxOfTemp,maxBoxOfTemp))
                                        {
                                            s2P0.x = otherNodeCoordinate[0][NodeIndex[jFace][0]];
                                            s2P0.y = otherNodeCoordinate[1][NodeIndex[jFace][0]];
                                            s2P1.x = otherNodeCoordinate[0][NodeIndex[jFace][1]];
                                            s2P1.y = otherNodeCoordinate[1][NodeIndex[jFace][1]];

                                            if (IfWallFaceIntersectionWithThisZone(s1P0, s1P1, s2P0, s2P1))
                                            {
                                                keyActiveOfNodes[faceNodeIndex[faceIndex][0]] = 0;
                                                keyActiveOfNodes[faceNodeIndex[faceIndex][1]] = 0;
                                                intersectLinesWithWall.insert(faceIndex);
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }  
        }
    }

    DelPointer(minBoxOfLine);
    DelPointer(maxBoxOfLine);

    DelPointer(minBoxOfTemp);
    DelPointer(maxBoxOfTemp);
}

void OversetConfig::SearchZonesSearchNodesInLineBoxByLayer()
{

}

void OversetConfig::CollectEffectiveDonorCell()
{
    int nBoundaryFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int *keyActiveOfCells = grid->GetBlankIndex();

    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");
    if (isOversetSlip)
    {
        return;
    }
    else 
    {

        for (int iFace = 0; iFace < nBoundaryFace; iFace++)
        {
            int le = leftCellofFace[iFace];
            int re = rightCellofFace[iFace];
            if (keyActiveOfCells[re] == INTERPOLATION && keyActiveOfCells[le] == ACTIVE)
            {
                interpolateExpendList.insert(le);
            }
        }

        for (int iFace = nBoundaryFace; iFace < nTotalFace; iFace++)
        {
            int le = leftCellofFace[iFace];
            int re = rightCellofFace[iFace];
            if (keyActiveOfCells[re] == INTERPOLATION && keyActiveOfCells[le] == ACTIVE)
            {
                interpolateExpendList.insert(le);
            }

            if (keyActiveOfCells[le] == INTERPOLATION && keyActiveOfCells[re] == ACTIVE)
            {
                interpolateExpendList.insert(re);
            }
        }
    }
}

RDouble OversetConfig::ComputeDistanceInDonorZone(RDouble *nodeCoordinateIn, int donorCell, bool wallDist)
{
    if (donorCell < 0)
    {
        return LARGE;
    }
    if (!wallDist)
    {
        RDouble *coorX = grid->GetCellCenterX();
        RDouble *coorY = grid->GetCellCenterY();
        RDouble *coorZ = grid->GetCellCenterZ();

        RDouble dx = nodeCoordinateIn[0] - coorX[donorCell];
        RDouble dy = nodeCoordinateIn[1] - coorY[donorCell];
        RDouble dz = nodeCoordinateIn[2] - coorZ[donorCell];
        RDouble dist = DISTANCE(dx, dy, dz);

        return dist;
    }
    else
    {
        int **cellNodeIndexContainer = grid->GetCell2NodeArray();
        int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

        RDouble *coordinateX = grid->GetX();
        RDouble *coordinateY = grid->GetY();
        RDouble *coordinateZ = grid->GetZ();

        int &numberOfNodeInCell = cellNodeNumberContainer[donorCell];

        RDouble weightDist = 0.0;
        RDouble weight = 0.0;

        for (int iNode = 0; iNode < numberOfNodeInCell; ++iNode)
        {
            int &nodeIndex = cellNodeIndexContainer[donorCell][iNode];

            RDouble dx = nodeCoordinateIn[0] - coordinateX[nodeIndex];
            RDouble dy = nodeCoordinateIn[1] - coordinateY[nodeIndex];
            RDouble dz = nodeCoordinateIn[2] - coordinateZ[nodeIndex];
            RDouble dist = DISTANCE(dx, dy, dz) + TINY;

            weight += 1.0 / dist;
            weightDist += 1.0 / dist * nodesMinDistanceInThisBlock[nodeIndex];
        }

        weightDist /= ( weight + TINY );

        return weightDist;
    }
}

RDouble OversetConfig::ComputeDistanceInDonorZoneForFace(RDouble* nodeCoordinateIn, int donorFace)
{
    RDouble *coorX = grid->GetFaceCenterX();
    RDouble *coorY = grid->GetFaceCenterY();
    RDouble *coorZ = grid->GetFaceCenterZ();

    RDouble dx = nodeCoordinateIn[0] - coorX[donorFace];
    RDouble dy = nodeCoordinateIn[1] - coorY[donorFace];
    RDouble dz = nodeCoordinateIn[2] - coorZ[donorFace];
    RDouble dist = DISTANCE(dx, dy, dz);

    return dist;

}
//! Interpolation boundary is determined according to the local distance and neighbor distance.
//! Overset boundary can be ignored.
void OversetConfig::SpecifyInitialInterNodesByDistance()
{
    int numberOfNodes = grid->GetNTotalNode();

    keyBasicActiveNodes = new int[numberOfNodes];
    SetField(keyBasicActiveNodes, 0, numberOfNodes);

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (nodesMinDistanceInOtherBlock[iNode] > nodesMinDistanceInThisBlock[iNode])
        {
            keyBasicActiveNodes[iNode] = 1;
        }
    }

    DelPointer(nodesMinDistanceInOtherBlock);
}

//! Any two points of an cell pass through the physical plane, the cell is a boundary cell.
//! This method isn't applicable in the presence of background grid.
void OversetConfig::SpecifySpecicalCellsInteractWithWallBoundary()
{
    int numberOfCells = grid->GetNTotalCell();
    int numberOfNodes = grid->GetNTotalNode();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int *innerCellOrNot = new int[numberOfCells];

    SetField(innerCellOrNot, 0, numberOfCells);
    SetField(keyInBlockNodes, 0, numberOfNodes);

    RDouble *xNode = grid->GetX();
    RDouble *yNode = grid->GetY();
    RDouble *zNode = grid->GetZ();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        int node1 = cellNodeIndexContainer[iCell][0];

        int numberOfNodeInCell = cellNodeNumberContainer[iCell];

        for (int iNode = 1; iNode < numberOfNodeInCell; ++iNode)
        {
            int node2 = cellNodeIndexContainer[iCell][iNode];

            RDouble dist1 = nodesMinDistanceInOtherBlock[node1];
            RDouble dist2 = nodesMinDistanceInOtherBlock[node2];

            RDouble dx = xNode[node1] - xNode[node2];
            RDouble dy = yNode[node1] - yNode[node2];
            RDouble dz = zNode[node1] - zNode[node2];

            RDouble length = DISTANCE(dx, dy, dz);

            if (dist1 < length || dist2 < length)
            {
                if (dist2 > 0.5 * LARGE || dist1 > 0.5 * LARGE)
                {
                    innerCellOrNot[iCell] = 1;
                    break;
                }
            }
        }
    }

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        if (innerCellOrNot[iCell] == 1)
        {
            for (int iNode = 0; iNode < cellNodeNumberContainer[iCell]; ++iNode)
            {
                int nodeIndex = cellNodeIndexContainer[iCell][iNode];
                keyInBlockNodes[nodeIndex] = 1;
            }
        }
    }

    DelPointer(innerCellOrNot);
}

void OversetConfig::SetNegativeForBufferRegion()
{
    int numberOfNodes = grid->GetNTotalNode();
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyBufferOfNodes[iNode])
        {
            keyActiveOfNodes[iNode] = INACTIVE;
        }
    }
}

void OversetConfig::SetNegativeForInnerNodes()
{
    int numberOfNodes = grid->GetNTotalNode();

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyInBlockNodes[iNode])
        {
            keyActiveOfNodes[iNode] = INACTIVE;
        }
    }
}

void OversetConfig::SetNegativeForInnerAndInitialInterNodes()
{
    int numberOfNodes = grid->GetNTotalNode();

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (!keyBasicActiveNodes[iNode] || keyInBlockNodes[iNode])
        {
            keyActiveOfNodes[iNode] = INACTIVE;
        }
    }
}

void OversetConfig::SetNodeCharacterOnPhysicalBoundary()
{
    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    SetField(keyActiveOfNodes, 0, numberOfNodes);
    SetField(keyActiveOfCells, 0, numberOfCells);

    int numberOfPhysicalBoundaryNodes = static_cast<int>(physicalBoundaryNodeList.size());
    for (int iNode = 0; iNode < numberOfPhysicalBoundaryNodes; ++iNode)
    {
        keyActiveOfNodes[physicalBoundaryNodeList[iNode]] = ACTIVE;
    }
}

void OversetConfig::FillTheRemainingNodes()
{
    int numberOfNodes = grid->GetNTotalNode();

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyActiveOfNodes[iNode] == 0)
        {
            keyActiveOfNodes[iNode] = INTERPOLATION;
        }
    }
}

void OversetConfig::SetActivePropertyForSearchRegion()
{
    //! Because communicate_cell_iblank is the direct communication keyActiveOfCells.
    //! keyBuffer is determined by keyActiveOfCells array.
    OversetConfigFactoryMulti *overs = GetOversetConfigFactoryMulti();
    RDouble **pMinOfBlocks = overs->GetMinBoxOfBlocks();
    RDouble **pMaxOfBlocks = overs->GetMaxBoxOfBlocks();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int *keyActiveOfCells = GetKeyActiveOfCells();

    SetField(keyActiveOfNodes, 0, numberOfNodes);
    SetField(keyActiveOfCells, 0, numberOfCells);

    vector< RDouble * > nodeCoordinate;
    RDouble *x = grid->GetX();
    nodeCoordinate.push_back(x);
    RDouble *y = grid->GetY();
    nodeCoordinate.push_back(y);
    RDouble *z = grid->GetZ();
    nodeCoordinate.push_back(z);

    int blockIndex = grid->GetIBlock();
    int numberOfBlocks = overs->GetNumberOfBlocks();
    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        if (iBlock != blockIndex)
        {
            for (int iNode = 0; iNode < numberOfNodes; ++iNode)
            {
                if (IfNodeIsInBox(nodeCoordinate, iNode, pMinOfBlocks[iBlock], pMaxOfBlocks[iBlock], geometricDimension))
                {
                    keyActiveOfNodes[iNode] = ACTIVE;
                }
            }
        }
    }
}

void OversetConfig::SetSearchRegionByActiveProperty()
{
    //! keySearch is determined by keyActive array.
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    keySearchOfNodes = new int[numberOfNodes];
    keySearchOfCells = new int[numberOfCells];

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    SetField(keySearchOfNodes, keyActiveOfNodes, numberOfNodes);
    SetField(keySearchOfCells, keyActiveOfCells, numberOfCells);
}

void OversetConfig::DefineInterCellsByBoundary()
{
    int numberOfCells = grid->GetNTotalCell();
    int nBoundaryFace = grid->GetNBoundFace();
    int nTotal = numberOfCells + nBoundaryFace;;

    UnstructBCSet **unstructBCSet = grid->GetBCRecord();
    int *rightCellofFace = grid->GetRightCellOfFace();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    SetField(keyActiveOfCells, ACTIVE, nTotal);

    for (int iFace = 0; iFace < nBoundaryFace; iFace++)
    {
        int bcType = unstructBCSet[iFace]->GetKey();
        if (bcType == PHENGLEI::OVERSET)
        {
            int re = rightCellofFace[iFace];
            keyActiveOfCells[re] = INTERPOLATION;
            interpolateCellList.push_back(re);
            interpolateFaceList.push_back(iFace);
        }
    }
}

void OversetConfig::DefineActCellsByBoundary()
{
    int numberOfCells = grid->GetNTotalCell();
    int nBoundaryFace = grid->GetNBoundFace();

    UnstructBCSet **unstructBCSet = grid->GetBCRecord();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    for (int iFace = 0; iFace < nBoundaryFace; iFace++)
    {
        int bcType = unstructBCSet[iFace]->GetKey();
        if (bcType == PHENGLEI::OVERSET)
        {
            int re = rightCellofFace[iFace];
            keyActiveOfCells[re] = ACTIVE;
        }
    }
}

void OversetConfig::ResetCellCentersByInterCell()
{
    OversetDataProxy *oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForSend();
    int numberOfNeighbors = oversetDataProxy->GetNumberOfNeighbors();
    if (numberOfNeighbors == 0)
    {
        return;
    }

    PHVector1D < OversetCell * > *oversetCellContainer = oversetInformationProxy->GetOversetCellContainer();

    RDouble *xcc, *ycc, *zcc, *vol;
    xcc = grid->GetCellCenterX();
    ycc = grid->GetCellCenterY();
    zcc = grid->GetCellCenterZ();
    vol = grid->GetCellVolume();

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int *leftCellofFace = grid->GetLeftCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
    {
        int  numberOfCells        = oversetDataProxy->GetNumberOfCells(iNeighbor);
        int *cellIndexContainer   = oversetDataProxy->GetCellIndexContainer(iNeighbor);
        RDouble *interCellCenterX = oversetDataProxy->GetInterCellCenterX(iNeighbor);
        RDouble *interCellCenterY = oversetDataProxy->GetInterCellCenterY(iNeighbor);
        RDouble *interCellCenterZ = oversetDataProxy->GetInterCellCenterZ(iNeighbor);
        RDouble *interCellVolume  = oversetDataProxy->GetInterCellVolume(iNeighbor);
        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            int cellIndex = cellIndexContainer[iCell];

            for (int j = 0; j < faceNumberOfEachCell[cellIndex]; ++j)
            {
                int face = cell2Face[cellIndex][j];
                if (nBoundFace <= face)
                {
                    continue;
                }

                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[face]);
                int bcType = bcRegion->GetBCType();
                if (bcType != PHENGLEI::OVERSET)
                {
                    continue;
                }
                int le = leftCellofFace[face];
                int re = rightCellofFace[face];

                xcc[re] = interCellCenterX[iCell];
                ycc[re] = interCellCenterY[iCell];
                zcc[re] = interCellCenterZ[iCell];
                vol[re] = interCellVolume[iCell];
            }
        }
    }
}

void OversetConfig::FindInterpolateExpend()
{
    int numberOfCells = grid->GetNTotalCell();
    int nBoundaryFace = grid->GetNBoundFace();
    int *leftCellofFace = grid->GetLeftCellOfFace();

    UnstructBCSet **unstructBCSet = grid->GetBCRecord();
    for (int iFace = 0; iFace < nBoundaryFace; iFace++)
    {
        int bcType = unstructBCSet[iFace]->GetKey();
        if (bcType == PHENGLEI::OVERSET)
        {
            int le = leftCellofFace[iFace];
            donorFaceList.insert(iFace);
        }
    }
}

//! It is an active cell if the all cell nodes are 1.
void OversetConfig::FindActiveCells()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndex = grid->GetCell2NodeArray();
    int *cellNodeNumber = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    SetField(keyActiveOfCells, ACTIVE, numberOfCells);

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        int &numberOfNodeInCell = cellNodeNumber[iCell];
        for (int iNode = 0; iNode < numberOfNodeInCell; ++iNode)
        {
            if (!keyActiveOfNodes[cellNodeIndex[iCell][iNode]])
            {
                keyActiveOfCells[iCell] = INACTIVE;
                break;
            }
        }
    }
}

void OversetConfig::SearchZonesSearchNodesInKDTree()
{
    RDouble *nodeCoordinate = new RDouble[THREE_D];
    int blockIndex = grid->GetIBlock();

    ZoneSearchNodes **zonesSearchNodes = GetOversetConfigFactoryMulti()->GetZonesSearchNodes();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        ZoneSearchNodes *&zoneSearchNodes = zonesSearchNodes[iZone];
        RDouble *pMin = zoneSearchNodes->GetPMin();
        RDouble *pMax = zoneSearchNodes->GetPMax();

        if (blockIndex != zoneSearchNodes->GetBlockIndex())
        {
            if (IfZoneSearchNodesIntersectionWithThisZone(pMin, pMax))
            {
                int numberOfSearchNodes = zoneSearchNodes->GetNumberOfNodes();

                RDouble *searchNodeCoordinateX = zoneSearchNodes->GetX();
                RDouble *searchNodeCoordinateY = zoneSearchNodes->GetY();
                RDouble *searchNodeCoordinateZ = zoneSearchNodes->GetZ();

                RDouble *nodesMinDistanceToOtherBlock = zoneSearchNodes->GetNodesMinDistanceToOtherBlock();

                for (int iNode = 0; iNode < numberOfSearchNodes; ++iNode)
                {
                    RDouble &nodesMinDistanceInOtherBlockOfThisNode = nodesMinDistanceToOtherBlock[iNode];

                    nodeCoordinate[0] = searchNodeCoordinateX[iNode];
                    nodeCoordinate[1] = searchNodeCoordinateY[iNode];
                    nodeCoordinate[2] = searchNodeCoordinateZ[iNode];

                    nodesMinDistanceInOtherBlockOfThisNode = MIN(nodesMinDistanceInOtherBlockOfThisNode,
                        ComputeDistanceInDonorZone(nodeCoordinate, FindDonorCell(nodeCoordinate)));
                }
            }
        }
    }

    DelPointer(nodeCoordinate);
}

int OversetConfig::FindDonorCell(RDouble *nodeCoordinateIn)
{
    int donorCellIndex = gridKDTree->SearchNodeInKDTree(nodeCoordinateIn);
    if (donorCellIndex != -1)
    {
        if (keyBufferOfCells[donorCellIndex]) donorCellIndex = -1;
    }

    return donorCellIndex;
}

//! It is an interpolation cell if the cell node has both 1 and 0.
void OversetConfig::FindInterCells()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        if (keyActiveOfCells[iCell] == INACTIVE)
        {
            int numberOfNodeInCell = cellNodeNumberContainer[iCell];
            for (int iNode = 0; iNode < numberOfNodeInCell; ++iNode)
            {
                int nodeIndex = cellNodeIndexContainer[iCell][iNode];
                if (keyActiveOfNodes[nodeIndex])
                {
                    keyActiveOfCells[iCell] = INTERPOLATION;
                    interpolateCellList.push_back(iCell);
                    break;
                }
            }
        }
    }
}

void OversetConfig::InitializeOversetInformation()
{
    //! Movement of the mesh, because of the ressembly of the mesh, it is necessary to delete oversetInformationProxy and regenerate it.
    if (oversetInformationProxy != nullptr)
    {
        FreePointer(oversetInformationProxy);
    }

    oversetInformationProxy = new OversetInformationProxy();
    oversetInformationProxy->SetGrid(grid);
    grid->SetOversetInformationProxy(oversetInformationProxy);

    PHVector1D < OversetCell * > *oversetCellContainer = new PHVector1D < OversetCell * >();

    InterCells **zonesInterpolateCells = GetOversetConfigFactoryMulti()->GetZonesInterpolateCells();
    InterCells *zoneInterpolateCells = zonesInterpolateCells[zoneIndex];

    int numberOfInterCells = zoneInterpolateCells->GetNumberOfInterCells();

    int *donorCellIndex = zoneInterpolateCells->GetDonorCellIndex();
    int *donorZoneIndex = zoneInterpolateCells->GetDonorZoneIndex();

    for (int iCell = 0; iCell < numberOfInterCells; ++iCell)
    {
        int interCell = interpolateCellList[iCell];
        int donorCell = donorCellIndex[iCell];
        int donorZone = donorZoneIndex[iCell];

        OversetCell *oversetCell = new OversetCell();
        oversetCell->SetCellIndex(interCell);
        oversetCell->SetZoneIndex(zoneIndex);
        oversetCell->SetDonorCellIndex(donorCell);
        oversetCell->SetDonorZoneIndex(donorZone);

        oversetCellContainer->push_back(oversetCell);
    }

    oversetInformationProxy->InitializeReceivingInformation(oversetCellContainer);
}

void OversetConfig::PostProcessOversetInformation()
{
    oversetInformationProxy->PostProcessOversetInformation();
}

void OversetConfig::MarkGlobalInterfaceNodes(int iBlock, vector< int > &keyNodesGlobal)
{
    if (blockIndex != iBlock) return;

    int numberOfInterNodes = static_cast<int>(interfaceBoundaryNodeList.size());
    for (int iNode = 0; iNode < numberOfInterNodes; ++iNode)
    {
        int localNodeIndex = interfaceBoundaryNodeList[iNode];
        int globalNodeIndex = localToGlobalNodesIndex[localNodeIndex];
        keyNodesGlobal[globalNodeIndex] = 1;
    }
}

void OversetConfig::CalLocalToGlobalInterfaceNodes(int iBlock, vector< int > &keyNodesGlobal)
{
    if (blockIndex != iBlock) return;

    int numberOfInterNodes = static_cast<int>(interfaceBoundaryNodeList.size());
    localToGlobalInterNodesIndex.resize(numberOfInterNodes);

    for (int iNode = 0; iNode < numberOfInterNodes; ++iNode)
    {
        int localNodeIndex = interfaceBoundaryNodeList[iNode];
        int globalNodeIndex = localToGlobalNodesIndex[localNodeIndex];
        int globelInterNodeIndex = keyNodesGlobal[globalNodeIndex];

        //! Connection from local interNodes to global sequence.
        localToGlobalInterNodesIndex[iNode] = globelInterNodeIndex;
    }
}

void OversetConfig::SetCellLayersByActiveNodes()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndex = grid->GetCell2NodeArray();
    int *cellNodeNumber = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        int &numberOfNodeInCell = cellNodeNumber[iCell];
        for (int iNode = 0; iNode < numberOfNodeInCell; ++iNode)
        {
            if (keyActiveOfNodes[cellNodeIndex[iCell][iNode]])
            {
                keyActiveOfCells[iCell] = 1;
                break;
            }
        }
    }
}

//! If a cell has a node of 1, its attributes is 1, otherwise it is 0.
void OversetConfig::CharacterNodeToCell()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndex = grid->GetCell2NodeArray();
    int *cellNodeNumber = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        int &numberOfNodeInCell = cellNodeNumber[iCell];
        for (int iNode = 0; iNode < numberOfNodeInCell; ++iNode)
        {
            if (keyActiveOfNodes[cellNodeIndex[iCell][iNode]])
            {
                keyActiveOfCells[iCell] = ACTIVE;
                break;
            }
        }
    }
}

//! Property of the point on interface is judged by the property of the virtual cell on interface.
void OversetConfig::CharacterInterfaceToNode()
{
    int numberOfBoundaryFaces = grid->GetNBoundFace();

    UnstructBCSet **unstructBCSet = grid->GetBCRecord();

    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndex = grid->GetFace2Node();

    int *leftCellIndex = grid->GetLeftCellOfFace();
    int *rightCellIndex = grid->GetRightCellOfFace();

    int *keyActiveOfCells = grid->GetBlankIndex();

    int iCount = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++iFace)
    {
        int boundaryType = unstructBCSet[iFace]->GetKey();

        int numberOfNodesInFace = faceNodeNumber[iFace];

        int leftCell = leftCellIndex[iFace];
        int rightCell = rightCellIndex[iFace];

        if (boundaryType == PHENGLEI::INTERFACE)
        {
            if (keyActiveOfCells[leftCell] != ACTIVE && keyActiveOfCells[rightCell] != ACTIVE)
            {
                iCount += numberOfNodesInFace;
                continue;
            }

            for (int iNode = 0; iNode < numberOfNodesInFace; ++iNode)
            {
                int nodeIndex = faceNodeIndex[iCount];
                keyActiveOfNodes[nodeIndex] = ACTIVE;
                iCount++;
            }
        }
        else
        {
            iCount += numberOfNodesInFace;
        }
    }
}

void OversetConfig::CharacterCellToNode()
{
    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        if (keyActiveOfCells[iCell])
        {
            int &numberOfNodeInCell = cellNodeNumberContainer[iCell];
            for (int iNode = 0; iNode < numberOfNodeInCell; ++iNode)
            {
                int &nodeIndex = cellNodeIndexContainer[iCell][iNode];
                keyActiveOfNodes[nodeIndex] = ACTIVE;
            }
        }
    }
}

//! Property of the point on interface is judged by the property of the virtual cell on interface.
int OversetConfig::CharacterNodeToValidNodeByCell()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    int numberOfNewNodes = 0;

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        int cellNodeNumber = cellNodeNumberContainer[iCell];

        for (int iNode = 0; iNode < cellNodeNumber; ++iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];

            if (keyActiveOfNodes[nodeIndex])
            {
                keyActiveOfCells[iCell] = ACTIVE;
                break;
            }
        }

        if (keyActiveOfCells[iCell])
        {
            for (int iNode = 0; iNode < cellNodeNumber; ++iNode)
            {
                int &nodeIndex = cellNodeIndexContainer[iCell][iNode];
                if (keyBasicActiveNodes[nodeIndex] && !keyInBlockNodes[nodeIndex] && !keyActiveOfNodes[nodeIndex])
                {
                    keyActiveOfNodes[nodeIndex] = ACTIVE;
                    numberOfNewNodes++;
                }
            }
        }
    }
    return numberOfNewNodes;
}

void OversetConfig::CheckInterpolateRelationship(int iBlockIn,
    InterCells *zoneInterpolateCellsIn,
    int &numberOfInterpolateCellsInThisBlockIn,
    int &numberOfNegativeCellsInThisBlock_1In,
    int &numberOfNegativeCellsInThisBlock_2In,
    int &numberOfNegativeCellsInThisBlock_3In)
{
    if (iBlockIn == zoneInterpolateCellsIn->GetBlockIndex())
    {
        int zoneIndex = grid->GetZoneID();

        int numberOfInterpolateCells = zoneInterpolateCellsIn->GetNumberOfInterCells();

        int *donorLevel = zoneInterpolateCellsIn->GetKey();

        int *donorZoneIndex = zoneInterpolateCellsIn->GetDonorZoneIndex();
        int *donorCellIndex = zoneInterpolateCellsIn->GetDonorCellIndex();
        int *isDefectCell = zoneInterpolateCellsIn->IsDefectCell();

        RDouble *interpolateCellCenterX = zoneInterpolateCellsIn->GetX();
        RDouble *interpolateCellCenterY = zoneInterpolateCellsIn->GetY();
        RDouble *interpolateCellCenterZ = zoneInterpolateCellsIn->GetZ();

        for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
        {
            ++numberOfInterpolateCellsInThisBlockIn;
            int &donorLevelOfThisCell = donorLevel[iCell];
            switch (donorLevelOfThisCell)
            {
                case -3:
                    ++numberOfNegativeCellsInThisBlock_3In;
                    isDefectCell[iCell] = 1;
                    break;
                case -2:
                    ++numberOfNegativeCellsInThisBlock_2In;
                    isDefectCell[iCell] = 1;
                    break;
                case -1:
                    ++numberOfNegativeCellsInThisBlock_1In;
                    isDefectCell[iCell] = 1;
                    break;
            }
        }
    }
}

void OversetConfig::DownInterfaceNodeCharacter()
{
    int blockIndex = grid->GetIBlock();
    int **keyOfInterfaceNodesInBlocks = GetOversetConfigFactoryMulti()->GetKeyOfInterfaceNodesInBlocks();

    int numberOfInterfaceBoundaryNodes = static_cast<int>(interfaceBoundaryNodeList.size());
    for (int iNode = 0; iNode < numberOfInterfaceBoundaryNodes; ++ iNode)
    {
        keyActiveOfNodes[interfaceBoundaryNodeList[iNode]] =
        keyOfInterfaceNodesInBlocks[blockIndex][localToGlobalInterfaceNodesMap[iNode]];
    }
}

void OversetConfig::UpdateInterfaceNodeCharacter()
{
    int blockIndex = grid->GetIBlock();
    int **keyOfInterfaceNodesInBlocks = GetOversetConfigFactoryMulti()->GetKeyOfInterfaceNodesInBlocks();

    int numberOfInterfaceBoundaryNodes = static_cast<int>(interfaceBoundaryNodeList.size());
    for (int iNode = 0; iNode < numberOfInterfaceBoundaryNodes; ++ iNode)
    {
        int &valueGlobal = keyOfInterfaceNodesInBlocks[blockIndex][localToGlobalInterfaceNodesMap[iNode]];
        valueGlobal = MAX(valueGlobal, keyActiveOfNodes[interfaceBoundaryNodeList[iNode]]);
    }
}

void OversetConfig::ShowBoundaryCells()
{
    string fileName = AddSymbolToFileName("./results/boundaryCells.dat", zoneIndex);

    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfCells = grid->GetNTotalCell();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer = grid->GetFace2Node();
    int *leftCellIndexContainer = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    int numberOfNodesInZone = 0;
    int numberOfFacesInZone = 0;
    int numberOfCellsInZone = 0;

    vector< int > keyNodes;
    vector< int > keyFaces;
    vector< int > nodeList;
    vector< int > cellList;

    keyNodes.resize(numberOfNodes);
    keyFaces.resize(numberOfFaces);
    nodeList.resize(numberOfNodes);
    cellList.resize(numberOfCells);

    SetField(keyNodes, 0, numberOfNodes);
    SetField(keyFaces, 0, numberOfFaces);
    SetField(nodeList, -1, numberOfNodes);
    SetField(cellList, 0, numberOfCells);

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int **cellFaceIndexContainer = grid->GetCell2Face();

    int *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        bool innerCellOrNot = false;
        for (int iNode = 0; iNode < cellNodeNumberContainer[iCell]; ++iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];
            if (keyInBlockNodes[nodeIndex])
            {
                innerCellOrNot = true;
                continue;
            }
        }
        if (innerCellOrNot)
        {
            numberOfCellsInZone++;
            cellList[iCell] = numberOfCellsInZone;

            for (int iNodes = 0; iNodes < cellNodeNumberContainer[iCell]; ++iNodes)
            {
                keyNodes[cellNodeIndexContainer[iCell][iNodes]] = 1;
            }

            for (int iFaces = 0; iFaces < cellFaceNumberContainer[iCell]; ++iFaces)
            {
                keyFaces[cellFaceIndexContainer[iCell][iFaces]] = 1;
            }
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] == 1)
        {
            numberOfNodesInZone++;
            nodeList[iNode] = numberOfNodesInZone;
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 1)
        {
            numberOfFacesInZone++;
        }
    }

    if (numberOfCellsInZone == 0)
    {
        return;
    }

    file << "title     = \"overset grid \"\n";
    file << "variables = \"x\" \"y\" \"z\"\n";

    file << "ZONE T=\"Boundary Cells\" \n";

    if (PHSPACE::GetDim() == THREE_D)
    {
        file << "ZoneType = FEPolyhedron\n";
    }
    else
    {
        file << "ZoneType = FEPolygon\n";
    }

    int totalNumFaceNodes = 0;

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] != 1) continue;
        totalNumFaceNodes += faceNodeNumberContainer[iFace];
    }

    file << "Nodes    = " << numberOfNodesInZone << "\n";
    file << "Faces    = " << numberOfFacesInZone << "\n";
    file << "Elements = " << numberOfCellsInZone << "\n";
    file << "TotalNumFaceNodes = " << totalNumFaceNodes << "\n";
    file << "NumConnectedBoundaryFaces = 0\n";
    file << "TotalNumBoundaryConnections = 0\n";

    int numberOfWordsInEachLine = 5;
    int iPrint = 0;

    //! nodes coordinate x
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << x[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! nodes coordinate y
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << y[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! nodes coordinate z
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << z[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! face nodes number
    iPrint = 0;
    if (PHSPACE::GetDim() == THREE_D)
    {
        for (int iFace = 0; iFace < numberOfFaces; ++iFace)
        {
            if (keyFaces[iFace] != 1) continue;
            file << faceNodeNumberContainer[iFace] << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << "\n";
    }

    //! face nodes list
    int count = 0;
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];

        if (keyFaces[iFace] == 0)
        {
            count += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumberContainer[iFace]; ++iNode)
        {
            int nodeIndex = nodeList[faceNodeIndexContainer[count++]];
            file << nodeIndex << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! left cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 0) continue;

        int cellIndex = cellList[leftCellIndexContainer[iFace]];
        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! right cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int cellIndex;
        if (keyFaces[iFace] == 0) continue;
        int rightCell = rightCellIndexContainer[iFace];

        if (rightCell >= numberOfCells || rightCell < 0)
        {
            cellIndex = 0;
        }
        else
        {
            cellIndex = cellList[rightCellIndexContainer[iFace]];
        }

        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    PHSPACE::CloseFile(file);
}

void OversetConfig::ShowDifferentLayersCell(int iLayer)
{
    string fileName = AddSymbolToFileName("./results/layerCells.dat", iLayer);

    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfCells = grid->GetNTotalCell();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer = grid->GetFace2Node();
    int *leftCellIndexContainer = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    int numberOfNodesInZone = 0;
    int numberOfFacesInZone = 0;
    int numberOfCellsInZone = 0;

    vector< int > keyNodes;
    vector< int > keyFaces;
    vector< int > nodeList;
    vector< int > cellList;

    keyNodes.resize(numberOfNodes);
    keyFaces.resize(numberOfFaces);
    nodeList.resize(numberOfNodes);
    cellList.resize(numberOfCells);

    SetField(keyNodes, 0, numberOfNodes);
    SetField(keyFaces, 0, numberOfFaces);
    SetField(nodeList, -1, numberOfNodes);
    SetField(cellList, 0, numberOfCells);

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int **cellFaceIndexContainer = grid->GetCell2Face();

    int *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        bool innerCellOrNot = false;

        if (hasBeenSearched(iCell))
        {
            innerCellOrNot = true;
        }
        if (innerCellOrNot)
        {
            numberOfCellsInZone++;
            cellList[iCell] = numberOfCellsInZone;

            for (int iNodes = 0; iNodes < cellNodeNumberContainer[iCell]; ++iNodes)
            {
                keyNodes[cellNodeIndexContainer[iCell][iNodes]] = 1;
            }

            for (int iFaces = 0; iFaces < cellFaceNumberContainer[iCell]; ++iFaces)
            {
                keyFaces[cellFaceIndexContainer[iCell][iFaces]] = 1;
            }
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] == 1)
        {
            numberOfNodesInZone++;
            nodeList[iNode] = numberOfNodesInZone;
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 1)
        {
            numberOfFacesInZone++;
        }
    }

    if (numberOfCellsInZone == 0)
    {
        return;
    }

    file << "title     = \"overset grid \"\n";
    file << "variables = \"x\" \"y\" \"z\"\n";
    file << "ZONE T=\" Cells' Layer "<< iLayer <<"\"\n";
    if (PHSPACE::GetDim() == THREE_D)
    {
        file << "ZoneType = FEPolyhedron\n";
    }
    else
    {
        file << "ZoneType = FEPolygon\n";
    }

    int totalNumFaceNodes = 0;

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] != 1) continue;
        totalNumFaceNodes += faceNodeNumberContainer[iFace];
    }

    file << "Nodes    = " << numberOfNodesInZone << "\n";
    file << "Faces    = " << numberOfFacesInZone << "\n";
    file << "Elements = " << numberOfCellsInZone << "\n";
    file << "TotalNumFaceNodes = " << totalNumFaceNodes << "\n";
    file << "NumConnectedBoundaryFaces = 0\n";
    file << "TotalNumBoundaryConnections = 0\n";

    int numberOfWordsInEachLine = 5;
    int iPrint = 0;

    //! nodes coordinate x
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << x[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! nodes coordinate y
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << y[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! nodes coordinate z
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << z[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! face nodes number
    iPrint = 0;
    if (PHSPACE::GetDim() == THREE_D)
    {
        for (int iFace = 0; iFace < numberOfFaces; ++iFace)
        {
            if (keyFaces[iFace] != 1) continue;
            file << faceNodeNumberContainer[iFace] << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << "\n";
    }

    //! face nodes list
    int count = 0;
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];

        if (keyFaces[iFace] == 0)
        {
            count += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumberContainer[iFace]; ++iNode)
        {
            int nodeIndex = nodeList[faceNodeIndexContainer[count++]];
            file << nodeIndex << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! left cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 0) continue;

        int cellIndex = cellList[leftCellIndexContainer[iFace]];
        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! right cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int cellIndex;
        if (keyFaces[iFace] == 0) continue;
        int rightCell = rightCellIndexContainer[iFace];

        if (rightCell >= numberOfCells || rightCell < 0)
        {
            cellIndex = 0;
        }
        else
        {
            cellIndex = cellList[rightCellIndexContainer[iFace]];
        }

        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    PHSPACE::CloseFile(file);
}

void OversetConfig::ShowIntersectCells()
{
    string fileName = AddSymbolToFileName("./results/boundaryCells.dat", zoneIndex);

    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfCells = grid->GetNTotalCell();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer = grid->GetFace2Node();
    int *leftCellIndexContainer = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    int numberOfNodesInZone = 0;
    int numberOfFacesInZone = 0;
    int numberOfCellsInZone = 0;

    vector< int > keyNodes;
    vector< int > keyFaces;
    vector< int > nodeList;
    vector< int > cellList;

    keyNodes.resize(numberOfNodes);
    keyFaces.resize(numberOfFaces);
    nodeList.resize(numberOfNodes);
    cellList.resize(numberOfCells);

    SetField(keyNodes, 0, numberOfNodes);
    SetField(keyFaces, 0, numberOfFaces);
    SetField(nodeList, -1, numberOfNodes);
    SetField(cellList, 0, numberOfCells);

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int **cellFaceIndexContainer = grid->GetCell2Face();

    int *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        bool innerCellOrNot = false;
        for (int iNode = 0; iNode < cellNodeNumberContainer[iCell]; ++iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];
            if (keyActiveOfNodes[nodeIndex] == 0)
            {
                innerCellOrNot = true;
                continue;
            }
        }
        if (innerCellOrNot)
        {
            numberOfCellsInZone++;
            cellList[iCell] = numberOfCellsInZone;

            for (int iNodes = 0; iNodes < cellNodeNumberContainer[iCell]; ++iNodes)
            {
                keyNodes[cellNodeIndexContainer[iCell][iNodes]] = 1;
            }

            for (int iFaces = 0; iFaces < cellFaceNumberContainer[iCell]; ++iFaces)
            {
                keyFaces[cellFaceIndexContainer[iCell][iFaces]] = 1;
            }
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] == 1)
        {
            numberOfNodesInZone++;
            nodeList[iNode] = numberOfNodesInZone;
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 1)
        {
            numberOfFacesInZone++;
        }
    }

    if (numberOfCellsInZone == 0)
    {
        return;
    }

    file << "title     = \"overset grid \"\n";
    file << "variables = \"x\" \"y\" \"z\"\n";

    file << "ZONE T=\"Boundary Cells\" \n";

    if (PHSPACE::GetDim() == THREE_D)
    {
        file << "ZoneType = FEPolyhedron\n";
    }
    else
    {
        file << "ZoneType = FEPolygon\n";
    }

    int totalNumFaceNodes = 0;

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] != 1) continue;
        totalNumFaceNodes += faceNodeNumberContainer[iFace];
    }

    file << "Nodes    = " << numberOfNodesInZone << "\n";
    file << "Faces    = " << numberOfFacesInZone << "\n";
    file << "Elements = " << numberOfCellsInZone << "\n";
    file << "TotalNumFaceNodes = " << totalNumFaceNodes << "\n";
    file << "NumConnectedBoundaryFaces = 0\n";
    file << "TotalNumBoundaryConnections = 0\n";

    int numberOfWordsInEachLine = 5;
    int iPrint = 0;

    //! nodes coordinate x
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << x[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! nodes coordinate y
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << y[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! nodes coordinate z
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << z[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! face nodes number
    iPrint = 0;
    if (PHSPACE::GetDim() == THREE_D)
    {
        for (int iFace = 0; iFace < numberOfFaces; ++iFace)
        {
            if (keyFaces[iFace] != 1) continue;
            file << faceNodeNumberContainer[iFace] << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << "\n";
    }

    //! face nodes list
    int count = 0;
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];

        if (keyFaces[iFace] == 0)
        {
            count += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumberContainer[iFace]; ++iNode)
        {
            int nodeIndex = nodeList[faceNodeIndexContainer[count++]];
            file << nodeIndex << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! left cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 0) continue;

        int cellIndex = cellList[leftCellIndexContainer[iFace]];
        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! right cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int cellIndex;
        if (keyFaces[iFace] == 0) continue;
        int rightCell = rightCellIndexContainer[iFace];

        if (rightCell >= numberOfCells || rightCell < 0)
        {
            cellIndex = 0;
        }
        else
        {
            cellIndex = cellList[rightCellIndexContainer[iFace]];
        }

        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    PHSPACE::CloseFile(file);
}

void OversetConfig::ShowIntersectLines(fstream& flowFile)
{
    int nTotalNode = grid->GetNTotalNode();
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfBoundary = grid->GetNBoundFace();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    UnstructBCSet **unstructBCSet = grid->GetBCRecord();

    int **cellFaceIndexContainer = grid->GetCell2Face();
    int *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();
    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int **faceNodeIndex = grid->GetFace2NodeArray();

    vector<vector<int> > face2nodelist(numberOfBoundary);
    vector<int> linkmap;
    pair<int, string> iter;
    for (int jFace = 0; jFace < numberOfBoundary; ++jFace)
    {
        int bcType = unstructBCSet[jFace]->GetKey();
        string bcName = unstructBCSet[jFace]->GetBCName();
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            (iter).first = bcType;
            (iter).second = bcName;
        }
    }

    GetFace2NodeListForLine(grid, iter, linkmap, face2nodelist, intersectLinesWithWall);
    int NumPts = static_cast<int>(linkmap.size());
    int numElements = static_cast<int>(face2nodelist.size());

    RDouble *xBoundary = new RDouble [NumPts];
    RDouble *yBoundary = new RDouble [NumPts];
    RDouble *zBoundary = new RDouble [NumPts];
    RDouble *valueBoundary = new RDouble [NumPts];
    for (int iNode = 0; iNode < NumPts; ++ iNode)
    {
        int nodeIndex = linkmap[iNode];

        if (nodeIndex < nTotalNode)
        {
            xBoundary[iNode] = x[nodeIndex];
            yBoundary[iNode] = y[nodeIndex];
            zBoundary[iNode] = z[nodeIndex];
            valueBoundary[iNode] = keyActiveOfNodes[nodeIndex];
        }
    }

    int *cell2node = new int [numElements * 4];
    int *node_number_of_each_cell = new int [numElements];

    int count = 0;
    for (std::size_t i = 0; i < face2nodelist.size(); ++ i)
    {
        node_number_of_each_cell[i] = 4;
        uint_t np = face2nodelist[i].size();
        int index0 = face2nodelist[i][0];
        for (int ip = 0; ip < np; ++ ip)
        {
            cell2node[count] = face2nodelist[i][ip];
            count ++;
        }

        for (uint_t ip = np; ip < 4; ++ ip)
        {
            cell2node[count] = index0;
            count ++;
        }
    }
    if (NumPts == 0)
    {
        DelPointer(xBoundary);
        DelPointer(yBoundary);
        DelPointer(zBoundary);
        DelPointer(valueBoundary);
        DelPointer(cell2node);
        DelPointer(node_number_of_each_cell);

        return;
    }
    flowFile << "N = " << NumPts << "\n";
    flowFile << "E = " << numElements << "\n";
    flowFile << "f = FEPOINT" << "\n";
    flowFile << "ET = quadrilateral" << "\n";

    int wordWidth = 20;
    for (int iNode = 0; iNode < NumPts; ++ iNode)
    {
        flowFile << setw(wordWidth) << xBoundary[iNode]
            << setw(wordWidth) << yBoundary[iNode]
            << setw(wordWidth) << zBoundary[iNode]
            << setw(wordWidth) << valueBoundary[iNode];
        flowFile << "\n";
    }

    count = 0;
    for (int iCell = 0; iCell < numElements; ++ iCell)
    {
        for (int iNode = 0; iNode < 4; ++ iNode)
        {
            flowFile << cell2node[count] + 1 << "  ";
            count ++;
        }
        flowFile << "\n";
    }
    flowFile << endl;

    DelPointer(xBoundary);
    DelPointer(yBoundary);
    DelPointer(zBoundary);
    DelPointer(valueBoundary);
    DelPointer(cell2node);
    DelPointer(node_number_of_each_cell);
}

void OversetConfig::OutPutOversetVisualization(const string &name)
{
    ostringstream oss;
    oss << "OversetOutput/" << name << "/" << grid->GetZoneID() << ".dat";

    string fileName = oss.str();

    fstream file;

    OpenFile(file, fileName, ios_base::out | ios_base::trunc);

    OutPutOversetVisualizationHeader(file);
    OutPutOversetVisualizationCell(file);
    OutPutOversetVisualizationPoint(file);
    OutPutOversetVisualizationTail(file);

    CloseFile(file);
}

void OversetConfig::OutPutOversetVisualizationHeader(fstream &file)
{
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfCells = grid->GetNTotalCell();

    RDouble *nodeCoordinateX = grid->GetX();
    RDouble *nodeCoordinateY = grid->GetY();
    RDouble *nodeCoordinateZ = grid->GetZ();

    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();

    file << "Title = \"OverSetVisualization\"" << endl << endl;

    file << "Variables = \"X\" \"Y\" \"Z\"" << endl;
    file << "\"keyActiveOfCells\" \"keyBufferOfCells\" \"keySearchOfCells\" \"interpolateCellList\"" << endl;
    file << "\"keyActiveOfNodes\" \"keyBufferOfNodes\" \"keySearchOfNodes\" \"keyInBlockNodes\" " << endl;
    file << "\"keyBasicActiveNodes\"  \"nodesMinDistanceInOtherBlock\" " << endl;
    file << endl;

    file << "Zone T    = \"Zone_" << grid->GetZoneID() << "\"" << endl << endl;
    file << "VARLOCATION = ( [1-3,8-13]=NODAL, [4-7]=CELLCENTERED )" << endl << endl;
    if (THREE_D == geometricDimension)
    {
        file << "ZoneType = FEPolyhedron" << endl << endl;
    }
    else
    {
        file << "ZoneType = FEPolygon" << endl << endl;
    }

    int totalNumberOfFaceNodes = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        totalNumberOfFaceNodes += faceNodeNumber[iFace];
    }
    file << "Nodes                       = " << numberOfNodes << endl;
    file << "Faces                       = " << numberOfFaces << endl;
    file << "Elements                    = " << numberOfCells << endl;
    file << "TotalNumFaceNodes           = " << totalNumberOfFaceNodes << endl;
    file << "NumConnectedBoundaryFaces   = 0" << endl;
    file << "TotalNumBoundaryConnections = 0" << endl;
    file << endl;

    OutPutValueName(file, "nodeCoordinateX");
    OutPutPointerValue(file, numberOfNodes, nodeCoordinateX);

    OutPutValueName(file, "nodeCoordinateY");
    OutPutPointerValue(file, numberOfNodes, nodeCoordinateY);

    OutPutValueName(file, "nodeCoordinateZ");
    OutPutPointerValue(file, numberOfNodes, nodeCoordinateZ);
}

void OversetConfig::OutPutOversetVisualizationCell(fstream &file)
{
    int numberOfCells = grid->GetNTotalCell();

    OutPutValueName(file, "keyActiveOfCells");
    OutPutPointerValue(file, numberOfCells, GetKeyActiveOfCells());

    OutPutValueName(file, "keyBufferOfCells");
    OutPutPointerValue(file, numberOfCells, keyBufferOfCells);

    OutPutValueName(file, "keySearchOfCells");
    OutPutPointerValue(file, numberOfCells, keySearchOfCells);

    OutPutValueName(file, "interpolateCellList");
    OutPutVectorValue(file, numberOfCells, interpolateCellList);
}

void OversetConfig::OutPutOversetVisualizationPoint(fstream &file)
{
    int numberOfNodes = grid->GetNTotalNode();

    OutPutValueName(file, "keyActiveOfNodes");
    OutPutPointerValue(file, numberOfNodes, keyActiveOfNodes);

    OutPutValueName(file, "keyBufferOfNodes");
    OutPutPointerValue(file, numberOfNodes, keyBufferOfNodes);

    OutPutValueName(file, "keySearchOfNodes");
    OutPutPointerValue(file, numberOfNodes, keySearchOfNodes);

    OutPutValueName(file, "keyInBlockNodes");
    OutPutPointerValue(file, numberOfNodes, keyInBlockNodes);

    OutPutValueName(file, "keyBasicActiveNodes");
    OutPutPointerValue(file, numberOfNodes, keyBasicActiveNodes);

    OutPutValueName(file, "nodesMinDistanceInOtherBlock");
    OutPutPointerValue(file, numberOfNodes, nodesMinDistanceInOtherBlock);

}

void OversetConfig::OutPutOversetVisualizationTail(fstream &file)
{
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfCells = grid->GetNTotalCell();

    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int **faceNodeIndex = grid->GetFace2NodeArray();
    int *leftCellIndex  = grid->GetLeftCellOfFace();
    int *rightCellIndex = grid->GetRightCellOfFace();

    if (THREE_D == geometricDimension)
    {
        OutPutValueName(file, "faceNodeNumber");
        OutPutPointerValue(file, numberOfFaces, faceNodeNumber);
    }

    int numberOfWordsInEachLine = 5;
    int iPrint = 0;

    OutPutValueName(file, "faceNodeIndex");
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        for (int iNode = 0; iNode < faceNodeNumber[iFace]; ++iNode)
        {
            int nodeIndex = faceNodeIndex[iFace][iNode] + 1;
            file << nodeIndex << " ";
            ++iPrint;
            if (iPrint % numberOfWordsInEachLine == 0) file << endl;
        }
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << endl;
    file << endl;

    OutPutValueName(file, "leftCellIndex");
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int cellIndex = leftCellIndex[iFace] + 1;
        if (cellIndex > numberOfCells || cellIndex < 0) cellIndex = 0;
        file << cellIndex << " ";
        ++iPrint;
        if (iPrint % numberOfWordsInEachLine == 0) file << endl;
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << endl;
    file << endl;

    OutPutValueName(file, "rightCellIndex");
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int cellIndex = rightCellIndex[iFace] + 1;
        if (cellIndex > numberOfCells || cellIndex < 0) cellIndex = 0;
        file << cellIndex << " ";
        ++iPrint;
        if (iPrint % numberOfWordsInEachLine == 0) file << endl;
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << endl;
    file << endl;
}

void OversetConfig::ShowGrid()
{
    string fileName = AddSymbolToFileName("./results/oversetGrid.dat", zoneIndex);

    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    PrintGridByOverSetConditionCell(file, ACTIVE);
    PrintGridByOverSetConditionCell(file, INACTIVE);
    PrintGridByOverSetConditionCell(file, INTERPOLATION);

    PHSPACE::CloseFile(file);
}

void OversetConfig::PrintGridByOverSetConditionCell(fstream &file, int key)
{
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfCells = grid->GetNTotalCell();
    int zoneIndex = grid->GetZoneID();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer = grid->GetFace2Node();
    int *leftCellIndexContainer = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    int numberOfNodesInZone = 0;
    int numberOfFacesInZone = 0;
    int numberOfCellsInZone = 0;

    vector< int > keyNodes;
    vector< int > keyFaces;
    vector< int > nodeList;
    vector< int > cellList;

    keyNodes.resize(numberOfNodes);
    keyFaces.resize(numberOfFaces);
    nodeList.resize(numberOfNodes);
    cellList.resize(numberOfCells);

    SetField(keyNodes, 0, numberOfNodes);
    SetField(keyFaces, 0, numberOfFaces);
    SetField(nodeList, -1, numberOfNodes);
    SetField(cellList, 0, numberOfCells);

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int **cellFaceIndexContainer = grid->GetCell2Face();
    int *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        if (keyActiveOfCells[iCell] == key)
        {
            numberOfCellsInZone++;
            cellList[iCell] = numberOfCellsInZone;

            for (int iNodes = 0; iNodes < cellNodeNumberContainer[iCell]; ++iNodes)
            {
                keyNodes[cellNodeIndexContainer[iCell][iNodes]] = 1;
            }

            for (int iFaces = 0; iFaces < cellFaceNumberContainer[iCell]; ++iFaces)
            {
                keyFaces[cellFaceIndexContainer[iCell][iFaces]] = 1;
            }
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] == 1)
        {
            numberOfNodesInZone++;
            nodeList[iNode] = numberOfNodesInZone;
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 1)
        {
            numberOfFacesInZone++;
        }
    }

    if (numberOfCellsInZone == 0)
    {
        return;
    }

    file << "title     = \"overset grid \"\n";
    file << "variables = \"x\" \"y\" \"z\"\n";

    if (key == INACTIVE)
    {
        file << "ZONE T=\" Zone " << zoneIndex << " Inactive Cells \"\n";
    }
    else if (key == ACTIVE)
    {
        file << "ZONE T=\" Zone " << zoneIndex << " Active Cells \"\n";
    }
    else if (key == INTERPOLATION)
    {
        file << "ZONE T=\" Zone " << zoneIndex << " Interpolation Cells \"\n";
    }
    else if (key == DONORCELL)
    {
        file << "ZONE T=\" Zone " << zoneIndex << " Donor Cells \"\n";
    }

    if (PHSPACE::GetDim() == THREE_D)
    {
        file << "ZoneType = FEPolyhedron\n";
    }
    else
    {
        file << "ZoneType = FEPolygon\n";
    }

    int totalNumFaceNodes = 0;

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] != 1) continue;
        totalNumFaceNodes += faceNodeNumberContainer[iFace];
    }

    file << "Nodes    = " << numberOfNodesInZone << "\n";
    file << "Faces    = " << numberOfFacesInZone << "\n";
    file << "Elements = " << numberOfCellsInZone << "\n";
    file << "TotalNumFaceNodes = " << totalNumFaceNodes << "\n";
    file << "NumConnectedBoundaryFaces = 0\n";
    file << "TotalNumBoundaryConnections = 0\n";

    int numberOfWordsInEachLine = 5;
    int iPrint = 0;

    //! nodes coordinate x
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << x[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! nodes coordinate y
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << y[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! nodes coordinate z
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << z[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! face nodes number
    iPrint = 0;
    if (PHSPACE::GetDim() == THREE_D)
    {
        for (int iFace = 0; iFace < numberOfFaces; ++iFace)
        {
            if (keyFaces[iFace] != 1) continue;
            file << faceNodeNumberContainer[iFace] << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << "\n";
    }

    //! face nodes list
    int count = 0;
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];

        if (keyFaces[iFace] == 0)
        {
            count += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumberContainer[iFace]; ++iNode)
        {
            int nodeIndex = nodeList[faceNodeIndexContainer[count++]];
            file << nodeIndex << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! left cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 0) continue;

        int cellIndex = cellList[leftCellIndexContainer[iFace]];
        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //! right cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int cellIndex;
        if (keyFaces[iFace] == 0) continue;
        int rightCell = rightCellIndexContainer[iFace];

        if (rightCell >= numberOfCells || rightCell < 0)
        {
            cellIndex = 0;
        }
        else
        {
            cellIndex = cellList[rightCellIndexContainer[iFace]];
        }

        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";
}

void OversetConfig::PrintGridByOverSetConditionPoint(fstream &file, int key)
{
    int numberOfNodes = grid->GetNTotalNode();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int numberOfSpecifyNodes = 0;

    vector< int > nodeList;

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyActiveOfNodes[iNode] == key)
        {
            numberOfSpecifyNodes++;
            nodeList.push_back(iNode);
        }
    }

    if (numberOfSpecifyNodes != 0)
    {
        file << "Title=\"THE SPECIFY NODES\" \n";
        file << "Variables=\"x\" \"y\" \"z\" \n";

        if (key == 1)
        {
            file << "ZONE T=\"Active Nodes\" \n";
        }
        else if (key == -1)
        {
            file << "ZONE T=\"Negative Nodes\" \n";
        }

        file << "N=" << numberOfSpecifyNodes << " E=" << numberOfSpecifyNodes << " F=FEPOINT, ET=quadrilateral\n";

        for (int iSpecifyNodes = 0; iSpecifyNodes < numberOfSpecifyNodes; ++iSpecifyNodes)
        {
            int activeNodeIndex = nodeList[iSpecifyNodes];

            file << x[activeNodeIndex] << " ";
            file << y[activeNodeIndex] << " ";
            file << z[activeNodeIndex] << " ";
            file << "\n";
        }

        for (int iSpecifyNodes = 0; iSpecifyNodes < numberOfSpecifyNodes; ++iSpecifyNodes)
        {
            file << iSpecifyNodes + 1 << " ";
            file << iSpecifyNodes + 1 << " ";
            file << iSpecifyNodes + 1 << " ";
            file << iSpecifyNodes + 1 << " ";
            file << "\n";
        }
    }
}

void OversetConfig::OutKeyBuffer()
{
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    string fileName = AddSymbolToFileName("keyBuffer.dat", zoneIndex);
    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        file << keyBufferOfCells[iCell] << endl;
    }

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        file << keyBufferOfNodes[iNode] << endl;
    }

    CloseFile(file);
}

void OversetConfig::OutKeyActive()
{
    int numberOfNodes = grid->GetNTotalNode();
    string fileName = AddSymbolToFileName("keyActive.dat", zoneIndex);
    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    file << "keyActiveOfNodes" << endl;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        file << keyActiveOfNodes[iNode] << endl;
    }

    CloseFile(file);
}

void OversetConfig::OutputLocalNodeWalldistance()
{
    int numberOfNodes = grid->GetNTotalNode();

    string fileName = AddSymbolToFileName("localNodeWalldistance.dat", zoneIndex);
    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        file << nodeWalldistDistance[iNode] << endl;
    }

    PHSPACE::CloseFile(file);
}

int *OversetConfig::GetKeyActiveOfCells()
{
    int *keyActiveOfCells = this->grid->GetBlankIndex();
    return keyActiveOfCells;
}

RDouble GetToleranceForOversetBox()
{
    return GlobalDataBase::GetDoubleParaFromDB("toleranceForOversetBox");
}

int GetKeyEnlargeOfInterBoundary()
{
    return GlobalDataBase::GetIntParaFromDB("keyEnlargeOfInterBoundary");
}

int GetNCommunicateCellIBlank()
{
    return GlobalDataBase::GetIntParaFromDB("nCommunicateCellIBlank");
}

int GetKeyEnlargeOfNodes()
{
    return GlobalDataBase::GetIntParaFromDB("keyEnlargeOfBufferNodes");
}

int GetKeyEnlargeOfSearchNodes()
{
    return GlobalDataBase::GetIntParaFromDB("keyEnlargeOfSearchNodes");
}

int GetReadInAuxiliaryInnerGrid()
{
    return GlobalDataBase::GetIntParaFromDB("readInAuxiliaryInnerGrid");
}

int GetReadInAuxiliaryOuterGrid()
{
    return GlobalDataBase::GetIntParaFromDB("readInAuxiliaryOuterGrid");
}

int GetReadInSklFileOrNot()
{
    return GlobalDataBase::GetIntParaFromDB("readInSklFileOrNot");
}

int GetSymetryOrNot()
{
    return GlobalDataBase::GetIntParaFromDB("symetryOrNot");
}

int GetIfOutPutOversetVisualization()
{
    return GlobalDataBase::GetIntParaFromDB("outPutOversetVisualization");
}

template< typename T >
void PH_AllreduceSepModeForVector(vector< T > &sendingBuffer, vector< T > &receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op)
{
    T *send = &sendingBuffer[0];
    T *recv = &receivingBuffer[0];

    PH_AllreduceSepMode(send, recv, numberOfElements, mpiDataType, op);
}

}
