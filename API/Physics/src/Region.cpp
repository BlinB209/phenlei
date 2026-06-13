#include "Geometry.h"
#include "IO_FileName.h"
#include "TK_Parse.h"
#include "TK_Time.h"
#include "Pre_WalldistCompute.h"
#include "Geo_Sample.h"
#include "Task_ServerUpdateInterface.h"
#include "Post_WriteVisualFile.h"
#include "TK_Warning.h"
#include "OversetInformation.h"
#include "UnstructuredOversetConfig.h"
#include "Pre_HDF5File.h"
#include "Glb_Dimension.h"
#include "Task.h"
#include "Solver.h"
#include "Zone.h"
#include "Geo_UnstructBC.h"
#include "Pre_GridConversion.h"
#include "Geo_StructGrid.h"
#include "Pre_StructToUnstruct.h"
#include "PHIO.h"
#include <omp.h>
#include <Geo_StructMerge.h>
#include <Geo_UnstructMerge.h>
#include <Post_WriteTecplot.h>
#include <GridType.h>


namespace PHSPACE
{
LIB_EXPORT Region::Region()
{
    zoneContainer = new vector <Zone *>;
    GlobalZones = new vector <Zone *>;
    numberOfZoneInLocalProcess = 0;
}

LIB_EXPORT Region::~Region()
{
    for (std::size_t iZone = 0; iZone < zoneContainer->size(); ++ iZone)
    {
        delete (*zoneContainer)[iZone];    (*zoneContainer)[iZone] = NULL;
        (*GlobalZones)[iZone] = NULL;
    }

    delete zoneContainer;    zoneContainer = NULL;
    delete GlobalZones;    GlobalZones = NULL;

    DeleteLabelContainer();
}

LIB_EXPORT void Region::ReadGrid()
{
    int numberOfGridGroups = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfGridGroups");

    PHVectorString1D gridGroupNameList;

    gridGroupNameList.push_back(PHSPACE::GlobalDataBase::GetStrParaFromDB("gridfile"));

    for (int iGridFile = 1; iGridFile < numberOfGridGroups; ++ iGridFile)
    {
        ostringstream oss;

        oss << "gridfile" << iGridFile;

        gridGroupNameList.push_back(PHSPACE::GlobalDataBase::GetStrParaFromDB(oss.str()));
    }

    bool HDF5FileRead = ReadHDF5File(this, gridGroupNameList);

    if (!HDF5FileRead)
    {
        PrintToWindow("Grid type is old fts type !\n");
        ReadOrdinaryGrid(gridGroupNameList);

        int isOldGrid = 1;
        GlobalDataBase::UpdateData("isOldGrid", &isOldGrid, PHINT, 1);
    }

    //! Read parameters from file of "boundary_condition.hypara". Notice that the five reference variables for the component aerodynamic coefficients integration should be updated with "reynoldsReferenceLengthDimensional".
    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();

    int isOverset = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    if (isOverset)
    {
        AddUnStructInfoToStructGrid();
    }

    //! iBlank (for unstructured grid) or cellTypeContainer (for structured grid) should be allocated and initialized
    //! whenever the overset grid is considered or not!
    AllocateOversetGrid();

    GridValidityCheck();
}

LIB_EXPORT void Region::InitGridForOversetConfig()
{
    ReadGrid();

    UpdateAllZonesPara();

    InitGeometry();

    SwapGeometry();
}

void ResetAllZonesOversetGrid()
{
    int isOverset = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    if (!isOverset)
    {
        return;
    }

    PHSPACE::ReConfigUnstructOversetGrid();
}

LIB_EXPORT void Region::ReadOrdinaryGrid(string fileName)
{
    using namespace PHMPI;

    PHVectorString1D gridGroupNameList;

    gridGroupNameList.push_back(fileName);

    if (! PHMPI::IsParallelRun())
    {
        ReadOrdinaryGridCSMode(gridGroupNameList);
    }
    else
    {
        ReadOrdinaryGridP2PMode(gridGroupNameList);
    }
}

LIB_EXPORT void Region::ReadOrdinaryGrid(string fileName, int fileMode)
{
    TK_Warning::FunctionDeprecation("Region::ReadOrdinaryGrid(string fileName, int fileMode)", "Region::ReadOrdinaryGrid(string fileName)");
    ReadOrdinaryGrid(fileName);
}

LIB_EXPORT void Region::ReadOrdinaryGrid(PHVectorString1D &gridGroupNameList)
{
    using namespace PHMPI;

    if (! PHMPI::IsParallelRun())
    {
        ReadOrdinaryGridCSMode(gridGroupNameList);
    }
    else
    {
        ReadOrdinaryGridP2PMode(gridGroupNameList);
    }
}

LIB_EXPORT int Region::GetProcessGlobalZoneIndexToLocalZoneIndex(int globalZoneIndex)
{
    int localZoneIndex = -1;

    int nLocalZones = PHMPI::GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        if (GetProcessLocalZoneIndexToGlobalZoneIndex(iZone) == globalZoneIndex)
        {
            localZoneIndex = iZone;
            break;
        }
    }

    if (localZoneIndex == -1)
    {
        TK_Exit::PrintDebugInfoExit("Error: Could not find out the local zone index !");
    }

    return localZoneIndex;
}

void Region::AddZone(Zone *zone)
{
    zoneContainer->push_back(zone);
    GlobalZones->push_back(zone);
    if (zone)
    {
        int globalIndex = static_cast<int>(zoneContainer->size()) - 1;
        PHMPI::InsertGlobalZones(globalIndex);
        ++ numberOfZoneInLocalProcess;
        PHMPI::SetNumberofLocalZones(numberOfZoneInLocalProcess);
    }
}

LIB_EXPORT void Region::BuildAndDumpMixedGrid()
{
    string outGridFile = GlobalDataBase::GetStrParaFromDB("out_gfile");
    int nZones = PHMPI::GetNumberofGlobalZones();

    Grid **gridArray = new Grid * [nZones];
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        gridArray[iZone] = GetGrid(iZone);
    }

    //! Construct connections.
    ConstructGlobalInterfaceLink(gridArray, nZones);

    //! Dump out grid.
    DumpGrid(outGridFile, gridArray, nZones);
    WriteFaceBoundaryName(gridArray, nZones, outGridFile);

    delete [] gridArray;
}

void Region::CreateGridFileNameList(vector< string > &filenamelist)
{
    if (! PHMPI::IsParallelRun())
    {
        CreateGridFileNameListCSMode(filenamelist);
    }
    else
    {
        CreateGridFileNameListP2PMode(filenamelist);
    }
}

LIB_EXPORT void Region::WriteFaceBoundaryName(Grid **grids, int nZones, const string &gridFileName)
{
    set< pair<int, string> > bcNamePairSet;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = grids[iZone];
        if (grid->Type() == UNSTRUCTGRID)
        {
            UnstructGrid *unstrGrid = UnstructGridCast(grid);

            UnstructBCSet *unstructBCSet = unstrGrid->GetUnstructBCSet();
            int nBCRegion = unstructBCSet->GetnBCRegion();

            for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
            {
                UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
                string bcName = bcRegion->GetBCName();
                int bcType = bcRegion->GetBCType();

                pair<int, string> bcNamePair(bcType, bcName);
                bcNamePairSet.insert(bcNamePair);
            }
        }
        else if (grid->Type() == STRUCTGRID)
        {
            StructGrid *strGrid = StructGridCast(grid);
            StructBCSet *structBCSet = strGrid->GetStructBCSet();
            int  nBCRegion = structBCSet->GetnBCRegion();

            for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion ++)
            {
                StructBC *bcRegion = structBCSet->GetBCRegion(iBCRegion);

                int bcType = bcRegion->GetBCType();
                string bcName = bcRegion->GetBCName();
                if (bcName == "Connection") continue;

                pair<int, string> bcNamePair(bcType, bcName);
                bcNamePairSet.insert(bcNamePair);
            }
        }
    }

    DumpBoundaryInfo(bcNamePairSet);
}

void Region::CreateGridFileNameListCSMode(vector< string > &filenamelist)
{
    using namespace PHMPI;
    string gridfile = "grid.grd";
    GlobalDataBase::GetData("gridfile", &gridfile, PHSTRING, 1);

    string gridfile_proc = gridfile;
    GetFileNameofThisProcessor(gridfile_proc);

    filenamelist.push_back(gridfile_proc);
}

void Region::CreateGridFileNameListP2PMode(vector< string > &filenamelist)
{
    using namespace PHMPI;
    string gridFileName;
    GlobalDataBase::GetData("gridfile", &gridFileName, PHSTRING, 1);

    string word;
    FindNextWord(gridFileName, word, ".");
    int numberOfParts = PHSPACE::ExtractNumberOfPartFromFileName(gridFileName);

    if (GetNumberOfProcessor() == 1)
    {
        for (int iPart = 0; iPart < numberOfParts; ++ iPart)
        {
            string gridfileProc = gridFileName;
            gridfileProc = PHSPACE::AddSymbolToFileName(gridfileProc, "_", iPart);
            filenamelist.push_back(gridfileProc);
        }
    }
    else if (GetNumberOfProcessor() == numberOfParts)
    {
        string gridfileProc = gridFileName;
        gridfileProc = PHSPACE::AddSymbolToFileName(gridfileProc, "_", GetCurrentProcessorID());
        filenamelist.push_back(gridfileProc);
    }
    else
    {
        TK_Exit::PrintDebugInfoExit("Error: have not considered the case when nProcessor less than nParts !"); 
    }
}

void Region::InitZoneLayoutInformation(vector< string > &filenamelist)
{
    zoneStart = 0;
    for (std::size_t ifile = 0; ifile < filenamelist.size(); ++ ifile)
    {
        GridGroup *grid_group = new GridGroup(zoneStart);
        grid_group->InitZoneLayoutInformation(filenamelist[ifile]);
        zoneStart += grid_group->GetNZones();
        delete grid_group;

        if (PHMPI::IsParallelRun())
        {
            // For self mode, only read the first file to get number of zones.
            break;
        }
    }
}

LIB_EXPORT void Region::PreprocessLayoutInformation(vector< string > &filenamelist)
{
    GetZoneLayoutInformation(filenamelist);
    AllocateZoneLayoutInformation();
    InitZoneLayoutInformation(filenamelist);
}

LIB_EXPORT void Region::ReadOrdinaryGridP2PMode(PHVectorString1D &gridGroupNameList)
{
    // Read Grid from file to grid group.
    IO_GridReader *gridReader = new IO_GridReader(gridGroupNameList, PHSPACE::GetDim());

    gridReader->ReadGridsByP2PMode();
    
    // Classify the grid system type: unstructured, structured and mixed.\n
    // Based on the block type.
    //ClassifyGridSystem(); // The grid type has been classified by GridReader.
    
    // Create Zones from grid group which is included in the Grid Reader.
    //CreateZonesFromGridReader(gridReader);
    CreateZonesFromGridGroups(gridReader);

    if (gridGroupNameList.size() == 1)
    {
        gridReader->ReadInterpointInfoByP2PMode();

        gridReader->ReadCellToNodeByP2PMode();

        gridReader->ReadWallDistByP2PMode();

        gridReader->ReadFaceBCByP2PMode();
    }

    delete gridReader;

}

LIB_EXPORT void Region::ReadOrdinaryGridCSMode(PHVectorString1D &gridGroupNameList)
{
    IO_GridReader *gridReader = new IO_GridReader(gridGroupNameList, PHSPACE::GetDim());

    gridReader->ReadGridsByCSMode();

    //ClassifyGridSystem();  // The grid type has been classified by GridReader.
    
    CreateZonesFromGridGroups(gridReader);

    gridReader->ReadCellToNodeByCSMode();

    gridReader->ReadInterpointInfoByCSMode();

    gridReader->ReadWallDistByCSMode();

    gridReader->ReadFaceBCByCSMode();

    if (GetDim() == THREE_D)
    {
        gridReader->ReadFaceBcDirByCSMode();
    }

    delete gridReader;
}

void Region::ReadUnstructuredOversetGrid()
{
    string fileName;
    GlobalDataBase::GetData("oversetGridFileName", &fileName, PHSTRING, 1);

    fstream file;

    int currentProcessor = PHMPI::GetCurrentProcessorID();

    int serverProcesspr  = PHMPI::GetServerProcessorID();

    if (currentProcessor == serverProcesspr)
    {
        OpenFile(file, fileName, ios_base::in|ios_base::binary);
    }

    int nZones = PHMPI::GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        DataContainer *dataContainer = new DataContainer();

        int zoneProcessor = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessor == serverProcesspr)
        {
            ReadFile(file, dataContainer);
        }

        PHMPI::PH_Trade(dataContainer, serverProcesspr, zoneProcessor, iZone);

        if (currentProcessor == zoneProcessor)
        {
            UnstructGridCast(GetGrid(iZone, 0))->Decode(dataContainer, 1);
        }

        delete dataContainer;
    }

    if (currentProcessor == serverProcesspr)
    {
        CloseFile(file);
    }
}

void Region::ReadStructuredOversetGrid()
{
    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();

    PHMPI::CreateZoneStartPointContainer(nZones);

    PHMPI::CreateZoneStartCenterContainer(nZones);

    string linkFileName = PHSPACE::GlobalDataBase::GetStrParaFromDB("linkFileName");

    ReadGridTopologyData(linkFileName);    //! Only server mode is considered here.

    return;
}

void Region::ReadGridTopologyData(const string &fileName)
{
    using namespace PHMPI;
    int *zoneStartPointContainer  = PHMPI::GetZoneStartPointContainer();
    int *zoneStartCenterContainer = PHMPI::GetZoneStartCenterContainer();

    int serverProcessorIndex = GetServerProcessorID();
    fstream file;
    int globalZoneNumber, globalPointNumber, globalCenterNumber;
    int **linkMap = new int * [3];

    ParallelOpenFile(file, fileName, ios_base::in|ios_base::binary);
    PH_ReadBcastInteger(file, &globalZoneNumber, 1, serverProcessorIndex);
    PH_ReadBcastInteger(file, zoneStartPointContainer, globalZoneNumber, serverProcessorIndex);

    PH_ReadBcastInteger(file, &globalPointNumber, 1, serverProcessorIndex);

    for (int iDirection = 0; iDirection < 3; ++iDirection)
    {
        linkMap[iDirection] = new int [globalPointNumber];
        
        PH_ReadBcastInteger(file, linkMap[iDirection], globalPointNumber, serverProcessorIndex);
    }

    GenerateGlobalBinaryTree(zoneStartPointContainer, globalZoneNumber, globalPointNumber);

    PH_ReadBcastInteger(file, zoneStartCenterContainer, globalZoneNumber, serverProcessorIndex);
    PH_ReadBcastInteger(file, &globalCenterNumber, 1, serverProcessorIndex);

    int *characterPointContainer = new int [globalCenterNumber];
    RDouble *cellCenterXContainer = new RDouble [globalCenterNumber];
    RDouble *cellCenterYContainer = new RDouble [globalCenterNumber];
    RDouble *cellCenterZContainer = new RDouble [globalCenterNumber];

    PH_ReadBcastInteger(file, characterPointContainer, globalCenterNumber, serverProcessorIndex);
    PH_ReadBcastDouble(file, cellCenterXContainer, globalCenterNumber, serverProcessorIndex);
    PH_ReadBcastDouble(file, cellCenterYContainer, globalCenterNumber, serverProcessorIndex);
    PH_ReadBcastDouble(file, cellCenterZContainer, globalCenterNumber, serverProcessorIndex);

    ParallelCloseFile(file);

    for (int iZone = 0; iZone < globalZoneNumber; ++ iZone)
    {
        Grid *gridIn = PHSPACE::GetGrid(iZone);

        if (!gridIn) continue;

        StructGrid *grid = PHSPACE::StructGridCast(gridIn);
        grid->SetZoneStartPointLabel(zoneStartPointContainer[iZone]);

        grid->CreateLinkPointStructure();

        for (int iDirection = 0; iDirection < 3; ++ iDirection)
        {
            int *iLabel = &linkMap[iDirection][zoneStartPointContainer[iZone]];
            grid->SetLinkPointLabel(iLabel, iDirection);
        }

        grid->SetZoneStartCenterLabel(zoneStartCenterContainer[iZone]);

        grid->CreateCoreParameterContainer();

        int *iLabel = &characterPointContainer[zoneStartCenterContainer[iZone]];
        grid->SetHingedPointContainer(iLabel);

        RDouble *xLabel = &cellCenterXContainer[zoneStartCenterContainer[iZone]];
        grid->SetCoreCoordinateContainer(xLabel, 0);

        RDouble *yLabel = &cellCenterYContainer[zoneStartCenterContainer[iZone]];
        grid->SetCoreCoordinateContainer(yLabel, 1);

        RDouble *zLabel = &cellCenterZContainer[zoneStartCenterContainer[iZone]];
        grid->SetCoreCoordinateContainer(zLabel, 2);
    }

    for (int iDirection = 0; iDirection < 3; ++ iDirection)
    {
        delete [] linkMap[iDirection];
    }

    delete [] linkMap;
    delete [] characterPointContainer;
    delete [] cellCenterXContainer;
    delete [] cellCenterYContainer;
    delete [] cellCenterZContainer;

    linkMap = NULL;
    characterPointContainer = NULL;
    cellCenterXContainer = NULL;
    cellCenterYContainer = NULL;
    cellCenterZContainer = NULL;
}

void Region::DeleteLabelContainer()
{
    using namespace PHMPI;

    DeleteGlobalBinaryTree();

    DeleteZoneStartPointContainer();

    DeleteZoneStartCenterContainer();

    return;
}

LIB_EXPORT void Region::ReadOversetGrid()
{
    ReadUnstructuredOversetGrid();
}

LIB_EXPORT void Region::SetOversetGrid()
{
    using namespace PHMPI;

    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);
        Grid *grid = GetGrid(zoneID, 0);
        if (!grid) continue;
        if (grid->Type() == UNSTRUCTGRID) continue;
        StructGridCast(grid)->CopyOversetFromUnstructGrid();
    }
}

void Region::GetZoneLayoutInformation(vector< string > &gridfilelist)
{
    numberOfTotalZones = 0;
    for (std::size_t ifile = 0; ifile < gridfilelist.size(); ++ ifile)
    {
        GetZoneLayoutInformation(gridfilelist[ifile]);

        if (PHMPI::IsParallelRun())
        {
            //! For self mode, only read the first file to get number of zones.
            break;
        }
    }

    PrintToWindow("Number Of Total Zones = ", numberOfTotalZones, "\n");
}

void Region::GetZoneLayoutInformation(const string &filename)
{
    fstream file;
    //ParallelOpenFile(file, filename, ios_base::in|ios_base::binary);
    OpenSerialFile(file, filename, ios_base::in|ios_base::binary);

    ReadNumberOfZones(file);

    ParallelCloseFile(file);
}

void Region::ReadNumberOfZones(fstream &file)
{
    using namespace PHMPI;

    int file_proc;
    GetFileProc(file_proc);

    int nzones = 0;
    PH_Read_Bcast(file, &nzones, sizeof(int), file_proc);

    numberOfTotalZones += nzones;
}

void Region::AllocateZoneLayoutInformation()
{
    using namespace PHMPI;

    SetNumberOfGlobalZones(numberOfTotalZones);

    CreateZoneProcessorIDDump(numberOfTotalZones);
    CreateZoneProcessorID(numberOfTotalZones);
    CreateZoneGridID     (numberOfTotalZones);
    CreateZoneGridType   (numberOfTotalZones);
    CreateZoneFileID     (numberOfTotalZones);
}

void Region::ReadOversetGrid(const string &filename)
{
    using namespace PHMPI;

    fstream file;
    ParallelOpenFile(file, filename, ios_base::in|ios_base::binary);

    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ReadSingleOversetGrid(file, iZone);
    }

    ParallelCloseFile(file);
}

void Region::ReadGrid(fstream &file)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ReadSingleGrid(file, iZone);
    }
}

void Region::AllocateOversetGrid()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);
        Grid *grid = GetGrid(zoneID, 0);
        if (!grid) continue;
        grid->AllocateOversetGrid();
    }
}

void Region::AddUnStructInfoToStructGrid()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    Grid **gridArray = new Grid * [nLocalZones];
    for (int iZone = 0; iZone < nLocalZones; iZone ++)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);
        gridArray[iZone] = GetGrid(zoneID, 0);
    }

    Pre_StructToUnstruct *StructToUnstruct = new Pre_StructToUnstruct(gridArray, nLocalZones);
    StructToUnstruct->Run();
    Grid **unstructGridArray = StructToUnstruct->GetUnstructGrids();
    for (int iZone = 0; iZone < nLocalZones; iZone ++)
    {
        Grid *grid = gridArray[iZone];
        if (grid->Type() != STRUCTGRID) continue;
        StructGrid *strGrid = StructGridCast(grid);
        UnstructGrid *unstrGrid = UnstructGridCast(unstructGridArray[iZone]);
        strGrid->SetUnstrGrid(unstrGrid);
        strGrid->GetUnstrGrid()->ComputeMetrics();
    }
}

void Region::ReadSingleOversetGrid(fstream &file, int iZone)
{
    using namespace PHMPI;

    int send_proc = 0;
    int recv_proc = 0;

    GetSendRecvProc(iZone, send_proc, recv_proc);

    int myid = GetCurrentProcessorID();

    //! If the process neither sends nor receives, it will return.
    if (myid != send_proc && myid != recv_proc) return;

    DataContainer *cdata = new DataContainer();

    //! Read grid file ,and send to the corresponding process.
    ReadAbstractData(file, cdata, send_proc, recv_proc);

    if (myid == recv_proc)
    {
        Grid *grid = GetGrid(iZone, 0);
        grid->Decode(cdata, 1);
    }

    delete cdata;    cdata = nullptr;
}

void Region::CreateZonesFromGridReader(IO_GridReader *gridReader)
{
    numberOfTotalZones = gridReader->GetNumberofGlobalZones();

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfTotalZones; ++ iZone)
    {
        int send_proc = 0;
        int recv_proc = 0;
        GetSendRecvProc(iZone, send_proc, recv_proc);

        if (myid != send_proc && myid != recv_proc)
        {
            AddZone(0);
            continue;
        }

        if (myid == recv_proc)
        {
            Zone *zone = new Zone(iZone);
            Grid *grid = gridReader->GetGrid(iZone);
            int nLocalZones = PHMPI::GetNumberofLocalZones();
            grid->GetGridID()->SetLocalIndex(nLocalZones);

            zone->AddGridAndCopyZonePara(grid);
            zone->AddGridToGlobal(iZone);
            AddZone(zone);
        }
        else 
        {
            //myid == sendproc myid != recv_proc
            AddZone(0);
        }
    }
}

void Region::ExcangeProcessorIDForGridProcessor()
{
    using namespace PHMPI;

    int *block_procs      = GetZoneProcessorID();
    int *block_procs_grid = GetZoneProcessorIDForGrid();

    if (PHMPI::CurrentProcessorIsGridProcessor())
    {
        int nZones = GetNumberofGlobalZones();

        SetField(block_procs, block_procs_grid, nZones);
    }
}

void Region::CreateZonesFromGridGroups(IO_GridReader *gridReader)
{
    //ExcangeProcessorIDForGridProcessor();

    using namespace PHMPI;

    int myid = GetCurrentProcessorID();

    vector< GridGroup * > &gridGroups = gridReader->GetGridGroups();

    uint_t numberOfGridGroups = gridGroups.size();

    int zoneStart = 0;

    for (int iGroup = 0; iGroup < numberOfGridGroups; ++ iGroup)
    {
        GridGroup *gridGroup = gridGroups[iGroup];

        int nZones = gridGroup->GetNZones();

        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int globalZoneIndex = zoneStart + iZone;

            int recv_proc      = GetZoneProcessorID(globalZoneIndex);
            int recv_proc_grid = GetZoneProcessorIDForGrid(globalZoneIndex);

            if (myid == recv_proc || myid == recv_proc_grid)
            {
                Zone *zone = new Zone(globalZoneIndex);

                Grid *grid = gridGroup->GetGrid(iZone);

                int nLocalZones = PHMPI::GetNumberofLocalZones();

                grid->GetGridID()->SetLocalIndex(nLocalZones);

                zone->AddGridAndCopyZonePara(grid);

                zone->AddGridToGlobal(globalZoneIndex);

                AddZone(zone);
            }
            else
            {
                AddZone(0);
                continue;
            }
        }
        zoneStart += nZones;
    }
}

void Region::CreateZoneAndGrid(int iZone)
{
    using namespace PHMPI;

    int send_proc = 0;
    int recv_proc = 0;

    GetSendRecvProc(iZone, send_proc, recv_proc);

    int myid = GetCurrentProcessorID();

    //! If the process neither sends nor receives, it will return.
    if (myid != send_proc && myid != recv_proc)
    {
        //! There is no problem here,but the design is ReadGrid to add zone,so it is responsible here.
        AddZone(0);
        return;
    }

    if (myid == recv_proc)
    {
        Zone *zone = new Zone(iZone);
        Grid *grid = CreateGridGeneral(GetZoneGridType(iZone), new GridID(iZone, PHMPI::GetNumberofLocalZones()), 0, GetDim());

        zone->AddGridAndCopyZonePara(grid);
        zone->AddGridToGlobal(iZone);
        AddZone(zone);
    }
    else //myid == sendproc myid != recv_proc
    {
        AddZone(0);
    }
}

void Region::InterpretDataToGrid(DataContainer *cdata, int iZone)
{
    using namespace PHMPI;

    int send_proc = 0;
    int recv_proc = 0;

    GetSendRecvProc(iZone, send_proc, recv_proc);

    int myid = GetCurrentProcessorID();

    if (myid != recv_proc) return;

    //! Decode the read-in data to grid data.
    Grid *grid = GetGrid(iZone, 0);
    
    grid->Decode(cdata, 0);
}

LIB_EXPORT void Region::ReviseInterZoneIndexForOverSetGrids(int iZone)
{
    using namespace PHMPI;

    int sendProcessorIndex = 0;
    int receiveProcessorIndex = 0;

    GetSendRecvProc(iZone, sendProcessorIndex, receiveProcessorIndex);

    int currentProcessorIndex = GetCurrentProcessorID();

    if (currentProcessorIndex != receiveProcessorIndex) return;

    //! Decode the read-in data to grid data.
    Grid *grid = PHSPACE::GetGrid(iZone, 0);

    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();

    if (!interfaceInfo)
    {
        return;
    }

    int nmberOfInterfaces = interfaceInfo->GetNIFace();

    int *NeighborZoneIndexContainer = interfaceInfo->GetInterFace2ZoneID();

    for (int iInterface = 0; iInterface < nmberOfInterfaces; ++ iInterface)
    {
        NeighborZoneIndexContainer[iInterface] += zoneStart;
    }
}

void Region::RecoverGridFromFile(fstream &file, int iZone)
{
    using namespace PHMPI;

    int send_proc = 0;
    int recv_proc = 0;

    GetSendRecvProc(iZone, send_proc, recv_proc);

    DataContainer *cdata = new DataContainer();

    //! Read grid file, and send to the corresponding process.
    //! Read in nlen and the remaining file content.

    ReadAbstractData(file, cdata, send_proc, recv_proc);

    InterpretDataToGrid(cdata, iZone);

    FreePointer(cdata);
}

LIB_EXPORT void Region::ReadSingleGrid(fstream &file, int iZone)
{
    CreateZoneAndGrid(iZone);
    RecoverGridFromFile(file, iZone);
    ReviseInterZoneIndexForOverSetGrids(iZone);
}

LIB_EXPORT void Region::ClassifyGridSystem()
{
    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();

    set<int> grid_type;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        grid_type.insert(GetZoneGridType(iZone));
    }

    int sys_gridtype = GetSystemGridType(grid_type);
    GlobalDataBase::UpdateData("sys_gridtype", &sys_gridtype, PHINT, 1);

    vector<string> sys_gridtype_name;
    sys_gridtype_name.push_back("UNSTRUCTGRID");
    sys_gridtype_name.push_back("STRUCTGRID");
    sys_gridtype_name.push_back("MIXGRID");

    PrintToWindow("Grid Type: ", sys_gridtype_name[sys_gridtype], "\n");
}

void Region::InitVariableWallTemperature()
{
    string wallTemperatureFile = "";
    GlobalDataBase::GetData("wallTemperaturefile", &wallTemperatureFile, PHSTRING, 1);

    int metabolicWallTemperature = 0;
    GlobalDataBase::UpdateData("metabolicWallTemperature", &metabolicWallTemperature, PHINT, 1);

    if (wallTemperatureFile == "")
    {
        return;
    }

    if (!FileExist(wallTemperatureFile))
    {
        return;
    }

    metabolicWallTemperature = 1;
    GlobalDataBase::UpdateData("metabolicWallTemperature", &metabolicWallTemperature, PHINT, 1);

    PrintToWindow("\n  Temperature file exist, read wall temperature file!\n\n");

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int myid = GetCurrentProcessorID();

    for (int iZone =  0; iZone < nZones; ++ iZone)
    {
        int zoneProc = GetZoneProcessorID(iZone);
        if (myid == zoneProc)
        {
            Grid *grid = GetGrid(iZone);
            grid->InitVariableWallTemperature();
        }
    }

    RDouble wallTemperature = -1;
    GlobalDataBase::GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
    if (wallTemperature < 0.0)
    {
        wallTemperature = 288.0;
        GlobalDataBase::UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
    }
}

void Region::InitOriginalGrid()
{
    int dumpFlowOnOriginalGrid = 0;
    GlobalDataBase::GetData("dumpFlowOnOriginalGrid", &dumpFlowOnOriginalGrid, PHINT, 1);
    if (!dumpFlowOnOriginalGrid)
    {
        return;
    }

    Post_WriteTecplotByOriginalGrid *writeTecplotByOriginalGrid = GetWriteTecplotByOriginalGrid();
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int nDim           = GetDim();
    int nTotalZone     = PHMPI::GetNumberofGlobalZones();
    int *ordinaryGridIndexOfEachZone = new int [nTotalZone];
    for (int iZone = 0; iZone < nTotalZone; ++ iZone)
    {
        ordinaryGridIndexOfEachZone[iZone] = -1;
    }
    for (int iZone = 0; iZone < nTotalZone; ++ iZone)
    {
        if (!GetGrid(iZone, 0)) continue;
        Grid *grid             = GetGrid(iZone, 0);
        int  ordinaryGridIndex = grid->GetOrdinaryGridIndex();
        ordinaryGridIndexOfEachZone[iZone] = ordinaryGridIndex;
    }
    int *ordinaryGridIndexOfEachZoneGlobal = new int [nTotalZone];
    PH_AllReduce(ordinaryGridIndexOfEachZone, ordinaryGridIndexOfEachZoneGlobal, nTotalZone, MPI_MAX);
    DelPointer(ordinaryGridIndexOfEachZone);

    int nTotalblock = 0;
    for (int iZone = 0; iZone < nTotalZone; ++ iZone)
    {
        nTotalblock = MAX(nTotalblock, ordinaryGridIndexOfEachZoneGlobal[iZone]);
    }
    ++ nTotalblock;

    vector < vector < Grid* > > originalgridofBlock;
    originalgridofBlock.resize(nTotalblock);
    int    myid     = PHMPI::GetCurrentProcessorID();
    string gridFile = GlobalDataBase::GetStrParaFromDB("gridfile");

    string gridfileProc = AddSymbolToFileName(gridFile, "_", myid);
    bool   fileExist    = FileExist(gridfileProc);
    int nFileProcessors = 0;
    if (!fileExist)
    {
        if (myid < nTotalblock && UNSTRUCTGRID == systemGridType)
        {
            //! The server processor must be file processor.
            TK_Exit::FileOpenErrorExit(gridfileProc);
        }
    }
    else
    {
        nFileProcessors = 1;
    }
    int globalFileNumber = 0;
    PH_AllReduce(&nFileProcessors, &globalFileNumber, 1, MPI_SUM);
    PHMPI::SetNumberofFiles(globalFileNumber);
    WriteLogFile("Global grid files number: ", globalFileNumber);

    int *fileIDOfEachZone = new int [nTotalZone];
    for (int iZone = 0; iZone < nTotalZone; ++ iZone)
    {
        fileIDOfEachZone[iZone] = 0;
    }
    if (UNSTRUCTGRID == systemGridType)
    {
        for (int iFile = 0; iFile < globalFileNumber; ++iFile)
        {
            IO_ReadHDF5Grid *reader = new IO_ReadHDF5Grid(AddSymbolToFileName(gridFile, "_", iFile));
            reader->OpenFile();
            int nTotalBlocksInFile = 0;
            reader->GetNumberOfZones(nTotalBlocksInFile);
            int *blockIdx = new int[nTotalBlocksInFile];
            reader->GetZoneIndex(blockIdx);
            for (int iZone = 0; iZone < nTotalBlocksInFile; ++iZone)
            {
                fileIDOfEachZone[blockIdx[iZone]] = iFile;
            }

            DelPointer(blockIdx);
            reader->CloseFile();
            FreePointer(reader);
        }
    }

    int *originalGridProc = new int [nTotalblock];
    writeTecplotByOriginalGrid->SetOriginalGridProc(nTotalblock, originalGridProc);

    for (int iZone = 0; iZone < nTotalZone; ++ iZone)
    {
        int fileID = fileIDOfEachZone[iZone];
        IO_ReadHDF5Grid *gridReader = new IO_ReadHDF5Grid(AddSymbolToFileName(gridFile, "_", fileID));
        gridReader->OpenFile();
        int indexofBlock = ordinaryGridIndexOfEachZoneGlobal[iZone];
        if (myid == originalGridProc[indexofBlock])
        {
            Grid *singleGrid = 0;
            GridID *index = new GridID(iZone);
            if (UNSTRUCTGRID == systemGridType)
            {
                singleGrid = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
            }
            else if (STRUCTGRID == systemGridType)
            {
                singleGrid = CreateGridGeneral(STRUCTGRID, index, 0, GetDim());
            }
            gridReader->GetGridInfo(singleGrid);
            originalgridofBlock[indexofBlock].push_back(singleGrid);
        }
        gridReader->CloseFile();
        FreePointer(gridReader);
    }
    DelPointer(fileIDOfEachZone);

    if (UNSTRUCTGRID == systemGridType)
    {
        OriginalUnstructGridMerge originalUnstructGridMerge;
        originalUnstructGridMerge.SetOriginalGridOfBlock(originalgridofBlock, nTotalblock, nTotalZone, nDim);
        originalUnstructGridMerge.SetOrdinaryGridIndexOfEachZoneGlobal(ordinaryGridIndexOfEachZoneGlobal);
        originalUnstructGridMerge.Initial();
        originalUnstructGridMerge.Run();
        int nBlocks = originalUnstructGridMerge.GetnTotalblock();
        UnstructGrid **mergegrid = new UnstructGrid * [nBlocks];
        originalUnstructGridMerge.GetMergeGrid(mergegrid);
        Grid **grids = new Grid * [nBlocks];
        for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
        {
            Grid *singleGrid = 0;
            if (mergegrid[iBlock])
            {
                singleGrid = mergegrid[iBlock];
            }
            grids[iBlock] = singleGrid;
        }
        writeTecplotByOriginalGrid->SetOriginalGrid(grids);
        writeTecplotByOriginalGrid->Initialize();
        DelPointer(mergegrid);
    }
    else if (STRUCTGRID == systemGridType)
    {
        OriginalStructGridMerge originalStructGridMerge;
        originalStructGridMerge.Initial(originalgridofBlock, nTotalblock, nTotalZone, nDim);
        originalStructGridMerge.Run();
        int nBlocks = originalStructGridMerge.GetnTotalblock();
        StructGrid **mergegrid = new StructGrid * [nBlocks];
        originalStructGridMerge.GetMergeGrid(mergegrid);

        Grid **grids = new Grid * [nBlocks];
        for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
        {
            Grid *singleGrid = 0;
            if (mergegrid[iBlock])
            {
                singleGrid = mergegrid[iBlock];
            }
            grids[iBlock] = singleGrid;
        }
        writeTecplotByOriginalGrid->SetOriginalGrid(grids);
        writeTecplotByOriginalGrid->Initialize();
        DelPointer(mergegrid);
    }

    DelPointer(ordinaryGridIndexOfEachZoneGlobal);
}

LIB_EXPORT void Region::DataMonitor()
{
    if (!GlobalDataBase::IsExist("ifSetDataMonitor", PHINT, 1))
    {
        return;
    }

    int ifSetDataMonitor = 0;
    GlobalDataBase::GetData("ifSetDataMonitor", &ifSetDataMonitor, PHINT, 1);
    if (!ifSetDataMonitor) return;

    SampleLocationInfo *globalSampleLocation = new SampleLocationInfo();
    globalSampleLocation->Run();
}

LIB_EXPORT void Region::SetupMultiBlockRelationship()
{
    MultiBlocksInterface();

    MultiBlocksInterpoint();

    //! Changed by Guo Yongheng 20160825.
    //MultiBlocksOverInterface();
}

void Region::MultiBlocksInterface()
{
    InitZoneNeighborsInformation();
    SwapZoneNeighborsInformation();
}

void Region::MultiBlocksInterpoint()
{
    if (zoneConnectivityForPoint)
    {
        InitZoneNeighborsForPointInformation();
        SwapZoneNeighborsForPointInformation();
    }
}

void Region::MultiBlocksOverInterface()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    zone_over_link.resize(nZones);

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone) continue;

        OversetInfoProxy *oversetInfoProxy = GetOversetInfoProxy(iZone, 0);
        if (oversetInfoProxy == 0) continue;

        vector<int> &neighborZoneIndexArray = oversetInfoProxy->GetNeighborZoneIndexArray();

        //! zone_over_link[iZone][j] means the neighbor zones of iZone, zone numbers in ascending order.
        zone_over_link[iZone] = neighborZoneIndexArray;
    }

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = GetZoneProcessorID(iZone);

        //! Get the number of neighbor zone of zone_over_link[iZone].
        int numberOfNeighbor = static_cast<int>(zone_over_link[iZone].size());

        //! If it is parallel,exchange numberOfNeighbor, otherwise do nothing.
        PH_Bcast(&numberOfNeighbor, sizeof(int), proc, iZone);

        //! Set the zone_over_link[iZone] size to numberOfNeighbor, and do nothing if serial.
        zone_over_link[iZone].resize(numberOfNeighbor);

        //! When numberOfNeighbor = 0, zone_over_link[iZone][0] will cross the border.
        if (numberOfNeighbor <= 0) continue;

        //! If it is parallel,exchange zone_over_link[iZone], otherwise do nothing.
        PH_Bcast(&zone_over_link[iZone][0], numberOfNeighbor*sizeof(int), proc, iZone);
    }
}

void Region::InitZoneNeighborsInformation()
{
    //! These codes find all the neighbor zones of each zone, if this block has interface boundary.
    //! Here the neighbors also include this zone itself.
    //! Finally,the zone numbers of all neighbors are stored in zoneid according to ascending order.
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    zoneConnectivity = new ZoneConnectivity();
    zoneConnectivity->CreateAllZones(nZones);

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone) continue;

        InterfaceInfo *interfaceInfo = GetInterfaceInfo(iZone, 0);

        if (!interfaceInfo) continue;
        interfaceInfo->InitNeighbors();
    }

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone) continue;

        InterfaceInfo *interfaceInfo = GetInterfaceInfo(iZone, 0);
        if (!interfaceInfo) continue;

        int numberOfNeighbor = interfaceInfo->GetNumberOfNeighbor();

        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);

        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            //! zoneNeighbor stores the neighbor zones of iZone,zone numbers in ascending order.
            zoneNeighbor->AddZoneIndexOfNeighbor(interfaceInfo->GetZoneIndexOfNeighbor(iNeighbor));
        }
    }
}

void Region::InitZoneNeighborsForPointInformation()
{
    //! These codes find all the neighbor zones of each zone,it is not just interface information here,but also concurrent information.
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    //zoneConnectivityForPoint = new ZoneConnectivityForPoint();
    zoneConnectivityForPoint->CreateAllZones(nZones);

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }
        InterpointInformation *interpointInformation = GetInterpointInfo(iZone, 0);

        if (interpointInformation == 0)
        {
            continue;
        }
            interpointInformation->InitNeighbors();
    }

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (! zone)
        {
            continue;
        }
        InterpointInformation *interpointInformation = GetInterpointInfo(iZone, 0);
        if (interpointInformation == 0)
        {
            continue;
        }
        int numberOfNeighbor = interpointInformation->GetNumberOfNeighbor();

        ZoneNeighbor *zoneNeighborForPoint = zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);

        for (int ineighbor = 0; ineighbor < numberOfNeighbor; ++ ineighbor)
        {
            //! zoneNeighbor stores the neighbor zones of iZone,zone numbers in ascending order.
            zoneNeighborForPoint->AddZoneIndexOfNeighbor(interpointInformation->GetZoneIndexOfNeighbor(ineighbor));
        }
    }
}

void Region::GetProcessorZones(int iZone, ZoneNeighbor *zoneNeighbor, vector< int > *neighborProcessors, vector< set< int > > *processorZones)
{
    using namespace PHMPI;
    uint_t numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();

    //! Insert the neighbor processors into set,
    //! including itself if the neighbor zone belong to the same processor.
    // int currentProcessor = GetZoneProcessorID(iZone);
    set< int > neighborProcessorsList;
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(iNeighbor);
        int neighborProcessor = GetZoneProcessorID(neighborZoneIndex);

        neighborProcessorsList.insert(neighborProcessor);
    }

    int numberOfNeighborProcessors = static_cast<int>(neighborProcessorsList.size());
    processorZones->resize(numberOfNeighborProcessors);
    neighborProcessors->resize(numberOfNeighborProcessors);

    int count = 0;
    for (set< int >::iterator iter = neighborProcessorsList.begin(); iter != neighborProcessorsList.end(); ++ iter)
    {
        (*neighborProcessors)[count ++] = *iter;
    }

    //! Construct the zones belong to each neighbor processor.
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(iNeighbor);
        int neighborProcessor = GetZoneProcessorID(neighborZoneIndex);

        int iNeighborProc = -1;
        for (int iProc = 0; iProc < numberOfNeighborProcessors; ++ iProc)
        {
            if ((* neighborProcessors)[iProc] == neighborProcessor)
            {
                iNeighborProc = iProc;
                break;
            }
        }
        
        (*processorZones)[iNeighborProc].insert(neighborZoneIndex);
    }
}

void Region::BroadcastZoneNeighborInformationCollectBcast()
{
    WriteLogFile("Start broadcasting neighbor information to each zones by collect-bcast method ... ");

    //*********************************************************//
    // Step1: collection from each other processors by server.
    //*********************************************************//
    ServerCollectionNeighborInfor();

    //*********************************************************//
    // Step2: broadcasting to each other processors by server.
    //*********************************************************//
    ServerBcastNeighborInfor();

    //*********************************************************//
    // Step3: delete the neighbor zone information that have 
    //        nothing to do with the current zones.
    //*********************************************************//
    CleanNeighborInfor();

    WriteLogFile("End broadcasting neighbor information to each zones ... ");
}

void Region::ServerCollectionNeighborInfor()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int serverTmp = GetServerProcessorID();
    int myid   = GetCurrentProcessorID();
    int number_of_processor = GetNumberOfProcessor();

    //! Grid processors do not send information to server.
    if (PHMPI::CurrentProcessorIsGridProcessor())
    {
        return;
    }

    DataContainer *cdata = new DataContainer();
    cdata->MoveToBegin();

    //! Write the number of zones in this processor.
    int nZonesInThisProc = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = GetZoneProcessorID(iZone);
        if (myid == proc)
        {
            ++ nZonesInThisProc;
        }
    }
    PHWrite(cdata, nZonesInThisProc);

    //! Compress all of neighbor information into data-container.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = GetZoneProcessorID(iZone);

        if (myid == proc)
        {
            PHWrite(cdata, &iZone, 1);

            ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
            int numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();

            PHWrite(cdata, &numberOfNeighbor, 1);

            if (numberOfNeighbor <= 0) continue;

            int *neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor();
            PHWrite(cdata, neighborZoneIndex, numberOfNeighbor);
        }
    }

    //! Send the data-container to server.
    if (myid != serverTmp)
    {
        int tag = myid;
        send(cdata, serverTmp, tag);
    }
    else
    {
        for (int iProc = 0; iProc < number_of_processor; ++ iProc)
        {
            int tag = iProc;

            if (myid != iProc)
            {
                receive(cdata, iProc, tag);
            }

            int nZonesInProc;
            cdata->MoveToBegin();
            PHRead(cdata, &nZonesInProc, 1);
            for (int iZone = 0; iZone < nZonesInProc; ++ iZone)
            {
                int zoneIndex;
                PHRead(cdata, &zoneIndex, 1);

                int numberOfNeighbor;
                PHRead(cdata, &numberOfNeighbor, 1);

                ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(zoneIndex);
                zoneNeighbor->SetNumberOfNeighbor(numberOfNeighbor);

                if (numberOfNeighbor <= 0) continue;

                int *neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor();
                PHRead(cdata, neighborZoneIndex, numberOfNeighbor);
            }

            delete cdata;    cdata = nullptr;
            cdata = new DataContainer();
        }
    }
    delete cdata;    cdata = nullptr;
}

void Region::ServerBcastNeighborInfor()
{
    using namespace PHMPI;
    int number_of_processor = GetNumberOfProcessor();
    if (number_of_processor <= 1) return;

    int tag = 0;
    int nZones = GetNumberofGlobalZones();
    int serverTmp = GetServerProcessorID();
    int myid   = GetCurrentProcessorID();

    DataContainer *cdata = new DataContainer();

    if (myid == serverTmp)
    {
        cdata->MoveToBegin();

        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
            int numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();

            PHWrite(cdata, &numberOfNeighbor, 1);

            if (numberOfNeighbor <= 0) continue;

            int *neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor();
            PHWrite(cdata, neighborZoneIndex, numberOfNeighbor);
        }
    }

    PH_BcastByServerToClinet(cdata, 2, tag);

    if (myid != serverTmp)
    {
        cdata->MoveToBegin();

        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int numberOfNeighbor;
            PHRead(cdata, &numberOfNeighbor, 1);

            ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
            zoneNeighbor->SetNumberOfNeighbor(numberOfNeighbor);

            if (numberOfNeighbor <= 0) continue;

            int *neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor();
            PHRead(cdata, neighborZoneIndex, numberOfNeighbor);
        }
    }

    delete cdata;    cdata = nullptr;
}

void Region::CleanNeighborInfor()
{
    using namespace PHMPI;
    int number_of_processor = GetNumberOfProcessor();
    if (number_of_processor <= 1) return;

    int nZones = GetNumberofGlobalZones();
    int myid   = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = GetZoneProcessorIDSepMode(iZone);
        if (proc == myid) continue;

        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
        uint_t numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();

        bool isNeighbor = false;
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int zoneIndexOfNeighbor = zoneNeighbor->GetZoneIndexOfNeighbor(iNeighbor);
            int neighborProc = GetZoneProcessorIDSepMode(zoneIndexOfNeighbor);
            if (neighborProc == myid)
            {
                isNeighbor = true;
                break;
            }
        }
        if (!isNeighbor)
        {
            zoneNeighbor->SetNumberOfNeighbor(0);
        }
    }
}

void Region::BroadcastZoneNeighborInformation()
{
    //! This is old code.
    bool isStructuredGrid = false;
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if (GetZone(iZone))
        {
            if (GetZone(iZone)->GetGeometry()->GetGrid()->Type() == STRUCTGRID)
            {
                isStructuredGrid = true;
                break;
            }
        }
    }

    int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (isOverset && sysGridType != STRUCTGRID)
    {
        BroadcastZoneNeighborInformationCollectBcast();
    }
    else
    {
        SwapZoneConnectivityBetweenNeighborZones();
    }
}

void Region::BroadcastZoneNeighborForPointInformation()
{
    int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (isOverset && sysGridType != STRUCTGRID)
    {
    }
    else
    {
        SwapZoneConnectivityForPointBetweenNeighborZones();
    }
}
void Region::BroadcastZoneNeighborInformationToEachProcessor()
{
    WriteLogFile("Start broadcasting neighbor information to each zones by old method... ");

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = GetZoneProcessorID(iZone);

        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();

        //! If it is parallel,exchange numberOfNeighbor, otherwise do nothing.
        PH_Bcast(&numberOfNeighbor, sizeof(int), proc, iZone);

        //! Set the zoneNeighbor data  size to numberOfNeighbor, and do nothing if serial.
        zoneNeighbor->SetNumberOfNeighbor(numberOfNeighbor);

        //! When numberOfNeighbor = 0,it will cross the border.
        if (numberOfNeighbor <= 0) continue;

        //! If it is parallel,exchange zoneNeighbor, otherwise do nothing.
        int *neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor();
        PH_Bcast(neighborZoneIndex, numberOfNeighbor*sizeof(int), proc, iZone);
    }

    WriteLogFile("End broadcasting neighbor information to each zones ... ");
}

void Region::BroadcastZoneNeighborInformationToNeighborZones()
{
    WriteLogFile("Start broadcasting neighbor information to neighbor zones ... ");

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int sizeOfInt = sizeof(int);
    int tag = 0;
    for (int zoneIndex = 0; zoneIndex < nZones; ++ zoneIndex)
    {
        if (!GetZone(zoneIndex)) continue;
        
        int currentZoneIndex = zoneIndex;

        ZoneNeighbor *zoneNeighborOfCurrentZone = zoneConnectivity->GetZoneNeighbor(zoneIndex);
        uint_t numberOfNeighborZones = zoneNeighborOfCurrentZone->GetNumberOfNeighbor();
        int *neighborZonesOfCurrentZone = zoneNeighborOfCurrentZone->GetZoneIndexOfNeighbor();

        //! Find the mid-zone.
        int midN, iZone;
        for (iZone = 0; iZone < numberOfNeighborZones; ++ iZone)
        {
            if (zoneNeighborOfCurrentZone->GetZoneIndexOfNeighbor(iZone) >= currentZoneIndex)
            {
                break;
            }
        }
        midN = iZone;
        
        //! Part1: receiving messages from neighbors with zone number less than myZone-1.
        for (int iNeighbor = 0; iNeighbor < midN; ++ iNeighbor)
        {
            int neighborZoneIndex  = zoneNeighborOfCurrentZone->GetZoneIndexOfNeighbor(iNeighbor);
            int sendProcessorIndex = GetZoneProcessorID(neighborZoneIndex);

            ZoneNeighbor *zoneNeighborOfNeighborZone = zoneConnectivity->GetZoneNeighbor(neighborZoneIndex);
            int size;
            PH_Receive(&size, 1 * sizeOfInt, sendProcessorIndex, tag);
            zoneNeighborOfNeighborZone->SetNumberOfNeighbor(size);

            int *neighborZonesOfNeighborZone = zoneNeighborOfNeighborZone->GetZoneIndexOfNeighbor();
            PH_Receive(neighborZonesOfNeighborZone, size * sizeOfInt, sendProcessorIndex, tag);
        }

        //! Part2: sending messages to neighbors with zone number big than myZone-1.
        if (midN < numberOfNeighborZones)
        {
            for (int iNeighbor = midN; iNeighbor < numberOfNeighborZones; ++ iNeighbor)
            {
                int neighborZoneIndex     = zoneNeighborOfCurrentZone->GetZoneIndexOfNeighbor(iNeighbor);
                int receiveProcessorIndex = GetZoneProcessorID(neighborZoneIndex);

                int size = zoneNeighborOfCurrentZone->GetNumberOfNeighbor();
                PH_Send(&size, 1 * sizeOfInt, receiveProcessorIndex, tag);
                PH_Send(neighborZonesOfCurrentZone, size * sizeOfInt, receiveProcessorIndex, tag);
            }
        }

        //! Part3: receiving messages from neighbors with zone number big than myZone-1.
        for (int iNeighbor = midN; iNeighbor < numberOfNeighborZones; ++ iNeighbor)
        {
            int neighborZoneIndex  = zoneNeighborOfCurrentZone->GetZoneIndexOfNeighbor(iNeighbor);
            int sendProcessorIndex = GetZoneProcessorID(neighborZoneIndex);

            ZoneNeighbor *zoneNeighborOfNeighborZone = zoneConnectivity->GetZoneNeighbor(neighborZoneIndex);
            int size;
            PH_Receive(&size, 1 * sizeOfInt, sendProcessorIndex, tag);
            zoneNeighborOfNeighborZone->SetNumberOfNeighbor(size);

            int *neighborZonesOfNeighborZone = zoneNeighborOfNeighborZone->GetZoneIndexOfNeighbor();
            PH_Receive(neighborZonesOfNeighborZone, size * sizeOfInt, sendProcessorIndex, tag);
        }

        //! Part4: sending messages to neighbors with zone number less than myZone-1.
        if (midN > 0)
        {
            for (int iNeighbor = 0; iNeighbor < midN; ++ iNeighbor)
            {
                int neighborZoneIndex     = zoneNeighborOfCurrentZone->GetZoneIndexOfNeighbor(iNeighbor);
                int receiveProcessorIndex = GetZoneProcessorID(neighborZoneIndex);

                int size = zoneNeighborOfCurrentZone->GetNumberOfNeighbor();
                PH_Send(&size, 1 * sizeOfInt, receiveProcessorIndex, tag);
                PH_Send(neighborZonesOfCurrentZone, size * sizeOfInt, receiveProcessorIndex, tag);
            }
        }
    }

    WriteLogFile("End broadcasting neighbor information to neighbor zones ... ");
}

void Region::SwapZoneConnectivityBetweenNeighborZones()
{
    WriteLogFile("Start broadcasting neighbor information to neighbor zones using new method... ");
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int currentProcessor = GetCurrentProcessorID();

    int sizeOfInt = sizeof(int);
    int tag = 0;
    for (int zoneIndex = 0; zoneIndex < nZones;  zoneIndex ++)
    {
        //! "This" refers to zoneIndex, "current" refers to information on current processor.
        ZoneNeighbor *zoneNeighborOfThisZone = zoneConnectivity->GetZoneNeighbor(zoneIndex);

        //! The proc id of this zone.
        int send_proc = GetZoneProcessorID(zoneIndex);

        //! This zone is on current processor, it may send message to others.
        if (currentProcessor == send_proc)
        {
            int numberOfNeighborZones = zoneNeighborOfThisZone->GetNumberOfNeighbor();
            //! if this zone do not have neighbor, do nothing.
            if (numberOfNeighborZones <= 0) continue;
            int *neighborZonesOfThisZone = zoneNeighborOfThisZone->GetZoneIndexOfNeighbor();
            set<int> processorHasBeenSend;
            for (int iNeighbor = 0; iNeighbor < numberOfNeighborZones; iNeighbor ++)
            {
                int zoneIndexOfNeighbor = zoneNeighborOfThisZone->GetZoneIndexOfNeighbor(iNeighbor);
                int recv_proc = GetZoneProcessorID(zoneIndexOfNeighbor);
                if (send_proc == recv_proc)
                {
                    //! This neighborZone is on the same processor of this zone.
                    continue;
                }
                else
                {
                    if (processorHasBeenSend.find(recv_proc) == processorHasBeenSend.end())
                    {
                        //! One target processor only send once.
                        processorHasBeenSend.insert(recv_proc);
                        int size = numberOfNeighborZones;
                        //WriteLogFile("Send #zone to #proc ", zoneIndex, recv_proc);
                        PH_Send(&size, 1 * sizeOfInt, recv_proc, tag);
                        PH_Send(neighborZonesOfThisZone, size * sizeOfInt, recv_proc, tag);
                    }
                }
            }
        }

        //! Current processor do not contain this zone, but it may contains this zone's neighbor zone.
        else
        {
            int nLocalZones = PHMPI::GetNumberofLocalZones();

            //! Loop all zones on current proc, check if has neighbor zone of zoneIndex.
            //! If has, receive zoneIndex's zoneNeighbor data.
            //! But only need receive once.
            for (int localZoneIndex = 0; localZoneIndex < nLocalZones; localZoneIndex ++)
            {
                int currentGlobalZone = GetLocalZoneIDToGlobalZoneID(localZoneIndex);
                ZoneNeighbor *zoneNeighborOfCurrentZone = zoneConnectivity->GetZoneNeighbor(currentGlobalZone);
                if (zoneNeighborOfCurrentZone->hasNeighborZone(zoneIndex))
                {
                    //WriteLogFile("Receive #Zone from #proc ", zoneIndex, send_proc, true);
                    int size=0;
                    PH_Receive(&size, 1 * sizeOfInt, send_proc, tag);
                    zoneNeighborOfThisZone->SetNumberOfNeighbor(size);
                    int *neighborZonesOfNeighborZone = zoneNeighborOfThisZone->GetZoneIndexOfNeighbor();
                    //! No allocate memory? through SetNumberOfNeighbor's resize.
                    PH_Receive(neighborZonesOfNeighborZone, size * sizeOfInt, send_proc, tag);
                    break;
                }
            }
        }
    }
    WriteLogFile("End broadcasting neighbor information to neighbor zones by new method... ");
}

void Region::SwapZoneConnectivityForPointBetweenNeighborZones()
{
    WriteLogFile("Start broadcasting neighbor information for point ");
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int currentProcessor = GetCurrentProcessorID();

    int sizeOfInt = sizeof(int);
    int tag = 0;
    for (int zoneIndex = 0; zoneIndex < nZones;  zoneIndex ++)
    {
        //! "This" refers to zoneIndex, "current" refers to information on current processor.
        ZoneNeighbor *zoneNeighborOfThisZone = zoneConnectivityForPoint->GetZoneNeighborForPoint(zoneIndex);
            
        //! The proc id of this zone.
        int send_proc = GetZoneProcessorID(zoneIndex);

        //! This zone is on current processor, it may send message to others.
        if (currentProcessor == send_proc)
        {
            int numberOfNeighborZones = zoneNeighborOfThisZone->GetNumberOfNeighbor();
            //! If this zone do not have neighbor, do nothing.
            if (numberOfNeighborZones <= 0) continue;
            int *neighborZonesOfThisZone = zoneNeighborOfThisZone->GetZoneIndexOfNeighbor();
            set<int> processorHasBeenSend; 
            for (int iNeighbor = 0; iNeighbor < numberOfNeighborZones; iNeighbor ++)
            {
                int zoneIndexOfNeighbor = zoneNeighborOfThisZone->GetZoneIndexOfNeighbor(iNeighbor);
                int recv_proc = GetZoneProcessorID(zoneIndexOfNeighbor);
                if (send_proc == recv_proc)
                {
                    //! This neighborZone is on the same processor of this zone.
                    continue;
                }
                else
                {
                    if (processorHasBeenSend.find(recv_proc) == processorHasBeenSend.end())
                    {
                        //! One target processor only send once.
                        processorHasBeenSend.insert(recv_proc);
                        int size =  numberOfNeighborZones;
                        //WriteLogFile("Send #zone to #proc ", zoneIndex, recv_proc);
                        PH_Send(&size, 1 * sizeOfInt, recv_proc, tag);
                        PH_Send(neighborZonesOfThisZone, size * sizeOfInt, recv_proc, tag);
                    }
                }
            }
        }

        //! Current processor do not contain this zone, but it may contains this zone's neighbor zone.
        else
        {
            int nLocalZones = PHMPI::GetNumberofLocalZones();

            //! Loop all zones on current proc, check if has neighbor zone of zoneIndex.
            //! If has,  receive zoneIndex's zoneNeighbor data.
            //! But only need receive once.
            for (int localZoneIndex = 0; localZoneIndex < nLocalZones; localZoneIndex ++)
            {
                int currentGlobalZone = GetLocalZoneIDToGlobalZoneID(localZoneIndex);
                ZoneNeighbor *zoneNeighborOfCurrentZone = zoneConnectivityForPoint->GetZoneNeighborForPoint(currentGlobalZone);
                if (zoneNeighborOfCurrentZone->hasNeighborZone(zoneIndex))
                {
                    //WriteLogFile("Receive #Zone from #proc ", zoneIndex, send_proc, true);
                    int size = 0;
                    PH_Receive(&size, 1 * sizeOfInt, send_proc, tag);
                    zoneNeighborOfThisZone->SetNumberOfNeighbor(size);
                    int *neighborZonesOfNeighborZone = zoneNeighborOfThisZone->GetZoneIndexOfNeighbor();
                    //! No allocate memory? through SetNumberOfNeighbor's resize.
                    PH_Receive(neighborZonesOfNeighborZone, size * sizeOfInt, send_proc, tag);
                    break;
                }
            }
        }
    }
    WriteLogFile("End broadcasting neighbor information for point ");
}

void Region::SwapZoneNeighborsInformation()
{
    WriteLogFile("Start Swapping Zone Neighbors Information!");

    TimeSpan timeSpan;
    BroadcastZoneNeighborInformation();
    timeSpan.ShowTimeSpanToLogFile("Neighbor zone information BCast "); 

    timeSpan.ResetTime();
    SwapNeighborsSendContent();
    //SwapNeighborsSendContentNonblocking();
    timeSpan.ShowTimeSpanToLogFile("Swap Neighbor Send Content ");

    WriteLogFile("End Swapping Zone Neighbors Information!");
}

void Region::SwapZoneNeighborsForPointInformation()
{
    BroadcastZoneNeighborForPointInformation();
    SwapNeighborsSendContentForPoint();
}

void Region::SwapNeighborsSendContentNonblocking()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector< PH_Request > requestContainer;
    vector< vector< DataContainer * > > receivedDataBuffer;
    vector< vector< DataContainer * > > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    //! Step 0: Compressing data firstly.
    WriteLogFile("        Compressing data ...");
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if (currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate and compress the buffers for sending.
            DataContainer *sendBuffer = 0;
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone][iNeighbor] = sendBuffer;

                //! Compress the send information into the actkey.
                //int localZoneIndex = this->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
                //compressData(zoneGridSolver, sendBuffer, neighborZone);
                InterfaceInfo *interfaceInfo = GetInterfaceInfo(iZone, 0);
                int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(iNeighbor);
                int *faceIndexForSend = interfaceInfo->ComputeFaceIndexForSend(iNeighbor);
                sendBuffer->MoveToBegin();
                PHWrite(sendBuffer, nIFaceOfNeighbor);
                PHWrite(sendBuffer, faceIndexForSend, nIFaceOfNeighbor);
            }
        }
    }

    //! Step 1: Communication.
    WriteLogFile("        Communication ...");
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        //! Communicating.
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor    = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[iZone].push_back(receivedBuffer);
            }
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //if (sendProcessor == receiveProcessor) continue;

            //! Communication.
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->Size();

                    //! Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                //! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                int neighborOrer = -1;
                ZoneNeighbor *globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                int numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (int neighborID = 0; neighborID < numberOfNeighborTemp; ++ neighborID)
                {
                    int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(neighborID);
                    if (neighborZoneTemp == iZone)
                    {
                        neighborOrer = neighborID;
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    //! The data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->Size();

                    //! Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    SWAP(receivedDataBuffer[iZone].back(), sendDataBuffer[iZone][iNeighbor]);
                }
            }
        }
    }

    //! Step2: MPI waiting for Non-blocking communications.
    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    //! Step3: Translate the data container.
    WriteLogFile("        Translating data ...");
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                receiveData->MoveToBegin();
                int nIFaceOfNeighbor;
                PHRead(receiveData, nIFaceOfNeighbor);
                int *faceIndexForSend = new int [nIFaceOfNeighbor];
                PHRead(receiveData, faceIndexForSend, nIFaceOfNeighbor);

                InterfaceInfo *tiinfo = GetInterfaceInfo(neighborZone, 0);
                tiinfo->FillFaceIndexForSend(iZone, faceIndexForSend);

                delete [] faceIndexForSend;

                ++ count;
            }
        }
    }

    WriteLogFile("        Free buffer data ...");
    //! Step4: Free the buffers.
    for (std::size_t iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < receivedDataBuffer[ iDim ].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

void Region::SwapNeighborsSendContent()
{
    WriteLogFile("Start SwapNeighborsSendContent ...");
    using namespace PHMPI;

    int myid = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();

    int level = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = GetZoneProcessorIDSepMode(iZone);

        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
        uint_t numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();

        for (int ineighbor = 0; ineighbor < numberOfNeighbor; ++ ineighbor)
        {
            int zoneIndexOfNeighbor = zoneNeighbor->GetZoneIndexOfNeighbor(ineighbor);
            int recv_proc = GetZoneProcessorIDSepMode(zoneIndexOfNeighbor);
            int nIFaceOfNeighbor = 0;
            int *faceIndexForSend = 0;

            if (myid == send_proc)
            {
                InterfaceInfo *interfaceInfo = GetInterfaceInfo(iZone, level);
                nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
                faceIndexForSend = interfaceInfo->ComputeFaceIndexForSend(ineighbor);
            }

            PH_Trade(&nIFaceOfNeighbor, sizeof(int), send_proc, recv_proc, iZone + level * nZones);

            if (myid == recv_proc && send_proc != recv_proc) faceIndexForSend = new int [nIFaceOfNeighbor]();
            PH_Trade(faceIndexForSend, nIFaceOfNeighbor * sizeof(int), send_proc, recv_proc, iZone + level * nZones);

            if (myid == recv_proc)
            {
                InterfaceInfo *tiinfo = GetInterfaceInfo(zoneIndexOfNeighbor, level);
                tiinfo->FillFaceIndexForSend(iZone, faceIndexForSend);
            }
            if (faceIndexForSend != nullptr)
            {
                delete [] faceIndexForSend;    faceIndexForSend = nullptr;
            }
        }
    }

    WriteLogFile("End SwapNeighborsSendContent ...");
}

void Region::SwapNeighborsSendContentForPoint()
{
    using namespace PHMPI;

    int myID = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();

    int level = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessor = GetZoneProcessorIDSepMode(iZone);

        ZoneNeighbor *zoneNeighborForPoint = zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);
        uint_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();

        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int zoneIndexOfNeighbor = zoneNeighborForPoint->GetZoneIndexOfNeighbor(iNeighbor);
            int receiveProcessor = GetZoneProcessorIDSepMode(zoneIndexOfNeighbor);
            int numberOfInterpointsForNeighbor = 0;
            int *pointIndexForSend = 0;

            if (myID == sendProcessor)
            {
                InterpointInformation *interpointInformation = GetInterpointInfo(iZone, level);
                numberOfInterpointsForNeighbor = interpointInformation->GetNumberOfInterpointsForNeighbor(iNeighbor);
                pointIndexForSend = interpointInformation->ComputePointIndexForSend(iNeighbor);
            }

            PH_Trade(&numberOfInterpointsForNeighbor, sizeof(int), sendProcessor, receiveProcessor, iZone + level * nZones);

            if (myID == receiveProcessor && sendProcessor != receiveProcessor)
            {
                pointIndexForSend = new int [numberOfInterpointsForNeighbor];
            }
            PH_Trade(pointIndexForSend, numberOfInterpointsForNeighbor*sizeof(int), sendProcessor, receiveProcessor, iZone + level * nZones);

            if (myID == receiveProcessor)
            {
                InterpointInformation *interpointInformation = GetInterpointInfo(zoneIndexOfNeighbor, level);
                interpointInformation->FillPointIndexForSend(iZone, pointIndexForSend);
            }

            delete [] pointIndexForSend;    pointIndexForSend = NULL;
        }
    }
}

LIB_EXPORT void Region::MultiGridSetup()
{
    CheckMeshMultigridByFTS();

    CoarseGrids();

    ShowCoarseGrids();
}

LIB_EXPORT void Region::InitGrids()
{
    SetNGrids();

    ComputeMetrics();

    InitCellBoundaryType();

#ifdef AI_RandomForestRegressor
    DumpCellCenterMetrics();
#endif
}

#ifdef AI_RandomForestRegressor
void Region::DumpCellCenterMetrics()
{
    ofstream dataFile;
    dataFile.open("./grid/CellCenter.txt", ofstream::app);

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = GetGrid(iZone);
        int nTotalCell = grid->GetNTotalCell();
        UnstructGrid *unsGrid = UnstructGridCast(grid);
        RDouble *xcc = unsGrid->GetCellCenterX();
        RDouble *ycc = unsGrid->GetCellCenterY();
        RDouble *zcc = unsGrid->GetCellCenterZ();

        vector<int> *cell2cell = unsGrid->GetCell2Cell();

        dataFile << iZone << " " << nTotalCell << endl;
        using namespace std;
        for(int i = 0; i<nTotalCell; ++i)
        {
            dataFile << xcc[i] << " " << ycc[i] << " " <<zcc[i] << "    ";

            for(int j = 0; j<cell2cell[i].size(); ++j)
            {
                dataFile << cell2cell[i][j] << " ";
            }
            dataFile << endl;
        }
    }
    dataFile.close();
}
#endif

void Region::DumpWallFaceCenter()
{
    int dumpWallFaceCenter = GlobalDataBase::GetIntParaFromDB("dumpWallFaceCenter");
    if (!dumpWallFaceCenter)
    {
        return;
    }

    using namespace PHMPI;
    ActionKey *actkeyDumpFaceCenter = new ActionKey();
    string faceCenterFile = "./grid/wallFaceCenter.dat";

    actkeyDumpFaceCenter->filename = faceCenterFile;
    ios_base::openmode openmode = ios_base::out|ios_base::trunc;
    actkeyDumpFaceCenter->openmode = openmode;

    int currentProcessorID = GetCurrentProcessorID();
    int serverTmp = GetServerSepMode();
 
    if (currentProcessorID == serverTmp)
    {
        ParallelOpenFile(actkeyDumpFaceCenter);
    }

    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        DataContainer *cdata = actkeyDumpFaceCenter->GetData();
        delete cdata;    cdata = nullptr;

        DataContainer *cdataNew = new DataContainer();
        actkeyDumpFaceCenter->SetData(cdataNew);

        int sendProcess = GetZoneProcessorIDSepMode(iZone);
        int recvProcess = GetServerSepMode();

        int myProcessorID = GetCurrentProcessorID();
        int sendRecvTag = GetSendRecvTag(actkeyDumpFaceCenter, iZone);

        if (myProcessorID == sendProcess)
        {
            //! If zone i belongs to the current process.
            Grid *grid = PHSPACE::GetGrid(iZone, 0);
            grid->DumpWallFaceCenter(actkeyDumpFaceCenter);
        }

        PH_Trade(actkeyDumpFaceCenter, sendProcess, recvProcess, sendRecvTag);

        if (myProcessorID == recvProcess)
        {
            //! If the current process is the server process.
            WriteASCIIFile(actkeyDumpFaceCenter);
        }
    }

    if (currentProcessorID == serverTmp)
    {
        ParallelCloseFile(actkeyDumpFaceCenter);
    }

    delete actkeyDumpFaceCenter;    actkeyDumpFaceCenter = nullptr;
}

void Region::SetNGrids()
{
    int nMGLevel = 1;
    GlobalDataBase::GetData("nMGLevel", &nMGLevel, PHINT, 1);

    nGrids = nMGLevel;
}

void Region::CheckMeshMultigridByFTS()
{
    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid   = PHMPI::GetCurrentProcessorID();

    int km = 0, NN = 0;
    int km_grid = 1024;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = PHMPI::GetZoneProcessorID(iZone);

        if (myid == procID)
        {
            Grid *gridofThisProcessor = GetGrid(iZone, 0);

            int gridtype = gridofThisProcessor->Type();

            if (gridtype == PHSPACE::STRUCTGRID)
            {
                StructGrid *grid = StructGridCast(gridofThisProcessor);

                StructBCSet *structBCSet = grid->GetStructBCSet();

                int nBCRegion = structBCSet->GetnBCRegion();
                for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion ++)
                {
                    km = 1;
                    NN = 2;

                    StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
                    int *s_st = structBC->GetStartPoint();
                    int *s_ed = structBC->GetEndPoint();

                    while (true)
                    {
                        int CandoMultigrid = 0;
                        for (int i = 0; i < 3; ++ i)
                        {
                            if (((abs(s_ed[i] - s_st[i]) + 1) % NN) == 0)
                            {
                                ++ CandoMultigrid;
                            }
                        }

                        if ((GetDim() == THREE_D && CandoMultigrid < 2) || (GetDim() == TWO_D && CandoMultigrid < 1))
                        {
                            break;
                        }

                        NN = NN * 2;
                        km = km + 1;
                    }
                    km_grid = min(km_grid, km);
                }
            }
        }
    }

    PH_CompareMaxMin(km_grid, 2);

    int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");

    if (nMGLevel > km_grid && km_grid < 1024)
    {
        GlobalDataBase::UpdateData("nMGLevel", &km_grid, PHINT, 1);
        this->nGrids = km_grid;
    }

    if (myid == server && km_grid < 1024)
    {
        PrintToWindow("\nThe most valid multi-grid level is: ", km_grid, "\n\n");
        WriteLogFile("\nThe most valid multi-grid level is: ", km_grid, "\n");
    }
}

void Region::CoarseGrids()
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone) continue;

        zone->CoarseGrids();
    }
}

void Region::ShowCoarseGrids()
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone); 
        if (!zone) continue;
        zone->ShowCoarseGrids();
    }
}

void Region::InitCellBoundaryType(int level)
{
    if (level != 0)
    {
        return;
    }

    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = GetZoneProcessorIDSepMode(iZone);
        int myid = PHMPI::GetCurrentProcessorID();
        if (myid == send_proc)
        {
            Grid *grid = GetGrid(iZone);
            grid->ComputeCellBoundaryType();
        }
    }
}

void Region::ComputeMetrics(int level)
{
    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();
    ActionKey *actkey = new ActionKey();
    actkey->action = COMPUTEMETRICS;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        //int send_proc = GetZoneProcessorIDSepMode(iZone);
        int send_proc = GetZoneProcessorIDSepMode(iZone);
        int recv_proc = GetServerSepMode();

        int tag = GetSendRecvTag(GRID_BASED, 0, iZone);
        //int proc = GetZoneProcessorID(iZone);
        int myid = PHMPI::GetCurrentProcessorID();
        if (myid == send_proc)
        {
            //! If the ProcessorID of i-th zone equal to the running ProcessorID, compute grid metrics.
            Grid *grid = GetGrid(iZone, level);
            grid->ComputeMetrics(actkey);
        }

        //! Collect all info out by ComputeMetrics function to cdata.
        PH_Trade(actkey, send_proc, recv_proc, tag);

        if (myid == recv_proc)
        {
            //! If ProcessorID is serve ProcessorID, output the info in cdata.
            WriteScreen(actkey);
        }
    }
    delete actkey;    actkey = nullptr;
}

void Region::GlobalGridStatic(int level)
{
    //! Check if need to compute the global grid size.
    int sys_gridtype;
    GlobalDataBase::GetData("sys_gridtype", &sys_gridtype, PHINT, 1);

    if (sys_gridtype == STRUCTGRID) 
    {
        ComputeGlobalCellSize(level);
    }
    else if (sys_gridtype == UNSTRUCTGRID)
    {
        ComputeGlobalGridSize(level);

#ifdef USE_GMRESSOLVER
        int originaltscheme = GlobalDataBase::GetIntParaFromDB("OriginalTscheme");
        if (originaltscheme == GMRES)
        {
            GetCurrentShiftCellSizeToGlobalSize(level);  // GMRESParallel
        }
#endif

        int compressible = GlobalDataBase::GetIntParaFromDB("compressible");
        if (compressible != COMPRESSIBLE)
        {
            buildGlobalCellsArrary(level);
        }
    }
    else if (sys_gridtype == MIXGRID)
    {
        ComputeGlobalGridCellSize(level);
    }
}

void Region::ComputeGlobalCellSize(int level)
{
    int minNTCell = 10000000;
    int maxNTCell = 0;
    RDouble TotalNodes = 0.0;
    RDouble TotalFaces = 0.0;
    RDouble TotalCells = 0.0;
    int TotalWallFace = 0;

    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }
        StructGrid *grid = StructGridCast(GetGrid(iZone, level));

        TotalNodes += grid->GetNTotalNode();
        TotalFaces += grid->GetNTotalFace();
        TotalCells += grid->GetNTotalCell();

        TotalWallFace += grid->GetNumberOfWallCell();

        minNTCell = MIN(minNTCell, TotalCells);
        maxNTCell = MAX(maxNTCell, TotalCells);
    }
    
    RDouble globalTotalNodes = TotalNodes * 1.0;
    RDouble globalTotalFaces = TotalFaces * 1.0;
    RDouble globalTotalCells = TotalCells * 1.0;
    int globalTotalWallFace = TotalWallFace;

    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_AllReduce(&TotalNodes, &globalTotalNodes, 1, PH_SUM);
        PH_AllReduce(&TotalFaces, &globalTotalFaces, 1, PH_SUM);
        PH_AllReduce(&TotalCells, &globalTotalCells, 1, PH_SUM);
        PH_AllReduce(&TotalWallFace, &globalTotalWallFace, 1, PH_SUM);
    }

    GlobalDataBase::UpdateData("GlobalTotalCells", &globalTotalCells, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("GlobalTotalWallFace", &globalTotalWallFace, PHINT, 1);

    PrintToWindow("The Global Information of StructGrid ", "\n");
    PrintToWindow("  GlobalTotalNodes: ", globalTotalNodes, "\n");
    PrintToWindow("  GlobalTotalFaces: ", globalTotalFaces, "\n");
    PrintToWindow("  GlobalTotalCells: ", globalTotalCells, "\n");

    RDouble minNTCellCompare = minNTCell * 1.0;
    RDouble maxNTCellCompare = maxNTCell * 1.0;

    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_CompareMaxMin(minNTCellCompare, 2);
        PH_CompareMaxMin(maxNTCellCompare, 1);
        WriteLogFile("Min && Max nTCell of all zones: ", minNTCellCompare, maxNTCellCompare);
    }
}

void Region::ComputeGlobalGridSize(int level)
{
    RDouble TotalNodes = 0.0;
    RDouble TotalFaces = 0.0;
    RDouble TotalCells = 0.0;
    RDouble TotalVolume = 0.0;
    RDouble averageVolume = 0.0;
    int TotalWallFace = 0;

    int nZones = PHMPI::GetNumberofGlobalZones();
    int nMGLevel = 1;
    GlobalDataBase::GetData("nMGLevel", &nMGLevel, PHINT, 1);
    for (int iLevel = 0; iLevel < nMGLevel; ++ iLevel)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            if (STRUCTGRID == PHMPI::GetZoneGridType(iZone))
            {
                continue;
            }
            Zone *zone = GetZone(iZone);
            if (!zone) continue;
            UnstructGrid *grid = UnstructGridCast(GetGrid(iZone, iLevel));
            RDouble *vol = grid->GetCellVolume();

            int nTotalNode = grid->GetNTotalNode();
            int nTotalFace = grid->GetNTotalFace();
            int nTotalCell = grid->GetNTotalCell();
            TotalNodes += nTotalNode * 1.0;
            TotalFaces += nTotalFace * 1.0;
            TotalCells += nTotalCell * 1.0;

            TotalWallFace += grid->GetNumberOfWallCell();

            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                TotalVolume += vol[iCell];
            }
        }

        RDouble globalTotalNodes = TotalNodes * 1.0;
        RDouble globalTotalFaces = TotalFaces * 1.0;
        RDouble globalTotalCells = TotalCells * 1.0;
        RDouble globalTotalVolume = TotalVolume;
        int globalTotalWallFace = TotalWallFace;

        if (PHMPI::GetNumberOfProcessor() > 1)
        {
            PH_AllReduce(&TotalNodes, &globalTotalNodes, 1, PH_SUM);
            PH_AllReduce(&TotalFaces, &globalTotalFaces, 1, PH_SUM);
            PH_AllReduce(&TotalCells, &globalTotalCells, 1, PH_SUM);
            PH_AllReduce(&TotalVolume, &globalTotalVolume, 1, PH_SUM);
            PH_AllReduce(&TotalWallFace, &globalTotalWallFace, 1, PH_SUM);
        }

        GlobalDataBase::UpdateData("GlobalTotalCells", &globalTotalCells, PHDOUBLE, 1);
        GlobalDataBase::UpdateData("UnstrGlobalTotalCells", &globalTotalCells, PHDOUBLE, 1);
        GlobalDataBase::UpdateData("UnstrGlobalTotalVolume", &globalTotalVolume, PHDOUBLE, 1);
        GlobalDataBase::UpdateData("GlobalTotalWallFace", &globalTotalWallFace, PHINT, 1);

        PrintToWindow("The Global Information of UnstructGrid ", "\n");
        PrintToWindow("  GlobalTotalNodes: ", globalTotalNodes, "\n");
        PrintToWindow("  GlobalTotalFaces: ", globalTotalFaces, "\n");
        PrintToWindow("  GlobalTotalCells: ", globalTotalCells, "\n");

        averageVolume = globalTotalVolume / globalTotalCells;
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            if (STRUCTGRID == PHMPI::GetZoneGridType(iZone))
            {
                continue;
            }
            Zone *zone = GetZone(iZone);
            if (!zone) continue;
            UnstructGrid *grid = UnstructGridCast(GetGrid(iZone, iLevel));
            grid->SetAverageVolume(averageVolume);
        }
    }
}
#ifdef USE_GMRESSOLVER
//! GMRESParallel
void Region::GetCurrentShiftCellSizeToGlobalSize(int level)
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    int nMGLevel = 1;
    GlobalDataBase::GetData("nMGLevel", &nMGLevel, PHINT, 1);

    RDouble *TotalCells = new RDouble[nZones];
    int *ShiftCellIndex = new int[nZones];

    for (int iLevel = 0; iLevel < nMGLevel; iLevel++)
    {
        for (int iZone = 0; iZone < nZones; iZone++)
        {
            TotalCells[iZone] = -1;
        }

        for (int iZone = 0; iZone < nZones; iZone++)
        {
            if (STRUCTGRID == PHMPI::GetZoneGridType(iZone))
            {
                continue;
            }
            Zone *zone = GetZone(iZone);
            if (!zone) continue;
            UnstructGrid *grid = UnstructGridCast(GetGrid(iZone, iLevel));

            int nTotalCell = grid->GetNTotalCell();
            TotalCells[iZone] += nTotalCell * 1.0;
        }

        if (PHMPI::GetNumberOfProcessor() > 1)
        {
            for (int iZone = 0; iZone < nZones; iZone++)
            {
                RDouble targetTotalCell = -1;
                PH_AllReduce(&TotalCells[iZone], &targetTotalCell, 1, PH_MAX);
                TotalCells[iZone] = targetTotalCell;
            }
        }

       //! determine whether the information is assigned
       for (int iZone = 0; iZone < nZones; iZone++)
       {
            if (STRUCTGRID == PHMPI::GetZoneGridType(iZone))
            {
                continue;
            }
            Zone *zone = GetZone(iZone);
            if (!zone) continue;

            if(TotalCells[iZone] < 0.0)
                TK_Exit::ExceptionExit("Global mesh information loading incorrectly");
       }

       //! get the shift cell index
       for (int iZone = 0; iZone < nZones; iZone++)
       {
            if (STRUCTGRID == PHMPI::GetZoneGridType(iZone))
            {
                continue;
            }

            int globalCellIDShift = 0;
            for (int zone = 0; zone < iZone; zone++)
            {
                globalCellIDShift += (int)TotalCells[zone];
            }
            globalCellIDShift += iZone;
            ShiftCellIndex[iZone] = globalCellIDShift;
            // printf("The number of cells in zone %d is %d\n", iZone, grid->GetNTotalCell());
            // printf("The shift cell index to global in zone %d is %d\n", iZone, globalCellIDShift);
       }
    }

    PH_Barrier();

    GlobalDataBase::UpdateDataPtr("shiftCellIndexToGlobal", ShiftCellIndex);

    delete[] TotalCells; TotalCells = NULL;

    //! get the local cell index
    for (int iLevel = 0; iLevel < nMGLevel; iLevel++)
    {
       for (int iZone = 0; iZone < nZones; iZone++)
       {
            if (STRUCTGRID == PHMPI::GetZoneGridType(iZone))
            {
                continue;
            }
            Zone *zone = GetZone(iZone);
            if (!zone)
            {
                continue;
            }

            UnstructGrid *grid = UnstructGridCast(GetGrid(iZone, iLevel));
            InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
            if(!interfaceInfo) 
            {
                continue;
            }

            int nNeighbor = interfaceInfo->GetNumberOfNeighbor();
            int *leftCellOfFace = grid->GetLeftCellOfFace();
            int *rightCellOfFace = grid->GetRightCellOfFace();
            for (int ineighbor = 0; ineighbor < nNeighbor; ineighbor++)
            {
                InterfacePatch *interFacePatch = interfaceInfo->GetInterfacePatch(ineighbor);

                int *faceID4Recv = interfaceInfo->GetFaceIndexForRecv(ineighbor);
                int zoneID = interfaceInfo->GetZoneIndexOfNeighbor(ineighbor);
                int nFace = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);

                int *cellIndexOfLocal = interFacePatch->GetCellIndexOfLocal();
                int *cellIndexOfLocalGhost = interFacePatch->GetCellIndexOfLocalGhost();
                int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

                for (int face = 0; face < nFace; face++)
                {
                    int realface = interFace2BoundaryFace[faceID4Recv[face]];
                    cellIndexOfLocal[face] = leftCellOfFace[realface];
                    cellIndexOfLocalGhost[face] = rightCellOfFace[realface];
                }
            }
       }
    }

    int* shiftCellID = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("shiftCellIndexToGlobal"));

    //! get the neighbor cell index
    for (int iLevel = 0; iLevel < nMGLevel; iLevel++)
    {
        int count = 0;
        for (int iZone = 0; iZone < nZones; iZone++)
        {
            if(STRUCTGRID ==PHMPI::GetZoneGridType(iZone))
            {
                continue;
            }
            Zone *zone = GetZone(iZone);
            if(!zone) continue;

            UnstructGrid *grid = UnstructGridCast(GetGrid(iZone, iLevel));
            InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
            if(!interfaceInfo) 
            {
                continue;
            }
            int nNeighbor = interfaceInfo->GetNumberOfNeighbor();
            for (int ineighbor = 0; ineighbor < nNeighbor; ineighbor++)
            {
                InterfacePatch *interFacePatch = interfaceInfo->GetInterfacePatch(ineighbor);

                int zoneID = interfaceInfo->GetZoneIndexOfNeighbor(ineighbor);
                int targetRank = PHMPI::GetZoneProcessorID(zoneID);

                if(targetRank != PHMPI::GetCurrentProcessorID())
                {
                    count++;
                }
            }
        }
        MPI_Request *request = new MPI_Request[count * 2];
        MPI_Status *status = new MPI_Status[count * 2];

        int icount = 0;
        for (int iZone = 0; iZone < nZones; iZone++)
        {
            if (STRUCTGRID == PHMPI::GetZoneGridType(iZone))
            {
                continue;
            }
            Zone *zone = GetZone(iZone);
            if (!zone)
            {
                continue;
            }

            UnstructGrid *grid = UnstructGridCast(GetGrid(iZone, iLevel));
            InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
            if(!interfaceInfo) 
            {
                continue;
            }

            int nIFace = interfaceInfo->GetNIFace();
            int *neighbor = interfaceInfo->GetInterFace2ZoneID();
            int nNeighbor = interfaceInfo->GetNumberOfNeighbor();
            int* leftCellOfFace = grid->GetLeftCellOfFace();
            int *rightCellOfFace = grid->GetRightCellOfFace();
            for (int ineighbor = 0; ineighbor < nNeighbor; ineighbor++)
            {
                InterfacePatch *interFacePatch = interfaceInfo->GetInterfacePatch(ineighbor);

                int zoneID = interfaceInfo->GetZoneIndexOfNeighbor(ineighbor);
                int nFace = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);

                int *cellIndexOfLocal = interFacePatch->GetCellIndexOfLocal();
                int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

                // data communication
                int targetRank = PHMPI::GetZoneProcessorID(zoneID);
                if (targetRank == PHMPI::GetCurrentProcessorID())
                {
                    UnstructGrid *gridTarget = UnstructGridCast(GetGrid(zoneID, iLevel));
                    InterfaceInfo *interfaceInfoTarget = gridTarget->GetInterfaceInfo();
                    int matchedNeighbor = interfaceInfoTarget->FindIthNeighbor(iZone);
                    InterfacePatch *interFacePatchTarget = interfaceInfoTarget->GetInterfacePatch(matchedNeighbor);
                    int matchednFace = interFacePatchTarget->GetNumberOfFace();
                    if(matchednFace != nFace)
                    {
                        printf("!idea: %d, matched: %d\n", nFace, matchednFace);
                        TK_Exit::ExceptionExit("Unmatched connectivy informatiion");
                    }
                    int *cellIndexTarget = interFacePatchTarget->GetCellIndexOfLocal();
                    interFacePatch->SetCellIndexOfNeighbor(cellIndexTarget);
                }
                else // not in one processor
                {
                    int *cellIndexOfNeighbor = new int[nFace];
                    MPI_Isend(cellIndexOfLocal, nFace, MPI_INT, targetRank, iZone, MPI_COMM_WORLD, &request[icount]);
                    MPI_Irecv(cellIndexOfNeighbor, nFace, MPI_INT, targetRank, zoneID, MPI_COMM_WORLD, &request[icount + 1]);
                    interFacePatch->SetCellIndexOfNeighbor(cellIndexOfNeighbor);
                    icount += 2;
                }
            }
        }

        MPI_Waitall(count * 2, request, status);
        delete[] request;
        request = NULL;
        delete[] status;
        status = NULL;

        for (int iZone = 0; iZone < nZones; iZone++)
        {
            if (STRUCTGRID == PHMPI::GetZoneGridType(iZone))
            {
                continue;
            }
            Zone *zone = GetZone(iZone);
            if (!zone)
            {
                continue;
            }

            UnstructGrid *grid = UnstructGridCast(GetGrid(iZone, iLevel));
            grid->SetGlobalCellIndexShift(shiftCellID[iZone]);
            InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
            if(!interfaceInfo) 
            {
                continue;
            }

            int nIFace = interfaceInfo->GetNIFace();
            int nNeighbor = interfaceInfo->GetNumberOfNeighbor();
            int* rightCellOfFace = grid->GetLeftCellOfFace();
            int *globalNeighborCellIndex = new int[nIFace];
            int *localInterfaceCellIndex = new int[nIFace];
            int *localInterfacePhysicalCellIndex = new int[nIFace];
            int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
            int iFace = 0;

            for (int ineighbor = 0; ineighbor < nNeighbor; ineighbor++)
            {
                InterfacePatch *interFacePatch = interfaceInfo->GetInterfacePatch(ineighbor);

                int zoneID = interfaceInfo->GetZoneIndexOfNeighbor(ineighbor);
                int nFace = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);

                int *cellIndexOfLocal = interFacePatch->GetCellIndexOfLocal();
                int *cellIndexOfNeighbor = interFacePatch->GetCellIndexOfNeighbor();
                int *cellIndexOfLocalGhost = interFacePatch->GetCellIndexOfLocalGhost();

                for (int n = 0; n < nFace; n++)
                {
                    globalNeighborCellIndex[iFace] = cellIndexOfNeighbor[n] +shiftCellID[zoneID];
                    localInterfaceCellIndex[iFace] = cellIndexOfLocalGhost[n];
                    localInterfacePhysicalCellIndex[iFace] = cellIndexOfLocal[n];
                    // int realFace = interFace2BoundaryFace[iFace];
                    // localInterfaceCellIndex[iFace] = leftCellOfFace[realFace];
                    // printf("!local zone %d, neighbor zone %d, face ID %d, target cell %d, local cell %d, shiftcell %d\n", iZone, zoneID, realFace, cellIndexOfNeighbor[n], localInterfaceCellIndex[iFace], globalNeighborCellIndex[iFace]);
                    iFace++;
                }
            }

            // UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
            // int nBCRegion = unstructBCSet->GetnBCRegion();
            // iFace = 0;
            // for (int ibc = 0; ibc < nBCRegion; ibc++)
            // {
            //     UnstructBC *bcRegion = unstructBCSet->GetBCRegion(ibc);
            //     int bcType = bcRegion->GetBCType();
            //     if(IsInterface(bcType))
            //     {
            //         vector<int> * face = bcRegion->GetFaceIndex();
            //         for (vector<int>::iterator iter = face->begin(); iter < face->end(); iter++)
            //         {
            //             localInterfaceCellIndex[iFace] = rightCellOfFace[*iter];
            //             // printf("!local zone %d, face ID %d, local cell %d, shiftcell %d\n", iZone,  *iter, localInterfaceCellIndex[iFace], globalNeighborCellIndex[iFace]);
            //             iFace++;
            //         }
            //     }
            // }

            interfaceInfo->SetGlobalNeighborCellIndex(globalNeighborCellIndex);
            interfaceInfo->SetLocalInterfaceCellIndex(localInterfaceCellIndex);
            interfaceInfo->SetLocalInterfacePhysicalCellIndex(localInterfacePhysicalCellIndex);
            // if(iZone == 0)
            // {
            //     for (int iter = 0; iter < nIFace; iter++)
            //     {
            //         printf("!left cell %d matched cell %d\n", localInterfaceCellIndex[iter] + 1, globalNeighborCellIndex[iter] + 1);
            //     }
            // }
        }
    }

    PH_Barrier();

    // TK_Exit::ExceptionExit("This is the debug information\n");
}
#endif
void Region::buildGlobalCellsArrary(int level)
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    int* cellNumbersOfEachZone = new int[nZones];
    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        cellNumbersOfEachZone[iZone] = 0;
    }

    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        if (STRUCTGRID == PHMPI::GetZoneGridType(iZone))
        {
            continue;
        }
        Zone* zone = GetZone(iZone);
        if (!zone) continue;
        UnstructGrid* grid = UnstructGridCast(GetGrid(iZone, 0));
        int nTotalCell = grid->GetNTotalCell();
        cellNumbersOfEachZone[iZone] = nTotalCell;
    }
    int* cellNumbersOfEachZoneGlobal = new int[nZones];
    PH_AllReduce(cellNumbersOfEachZone, cellNumbersOfEachZoneGlobal, nZones, MPI_SUM);

    int * globalCellNoOfGlobalZones = new int[nZones + 1];
    globalCellNoOfGlobalZones[0] = 0;
    int cellCount = 0;
    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        cellCount += cellNumbersOfEachZoneGlobal[iZone];
        globalCellNoOfGlobalZones[iZone + 1] = cellCount;
    }

    GlobalDataBase::UpdateDataPtr("globalCellNoOfGlobalZones", globalCellNoOfGlobalZones);

    delete [] cellNumbersOfEachZone;    cellNumbersOfEachZone = nullptr;
    delete [] cellNumbersOfEachZoneGlobal;    cellNumbersOfEachZoneGlobal = nullptr;
}

void Region::ComputeGlobalGridCellSize(int level)
{
    RDouble TotalNodes_unstr = 0.0;
    RDouble TotalFaces_unstr = 0.0;
    RDouble TotalCells_unstr = 0.0;

    RDouble TotalNodes_str = 0.0;
    RDouble TotalFaces_str = 0.0;
    RDouble TotalCells_str = 0.0;

    RDouble TotalVolume_unstr = 0.0;
    RDouble averageVolume_unstr = 0.0;

    int minNTCell_str = 10000000;
    int maxNTCell_str = 0;

    int nZones = PHMPI::GetNumberofGlobalZones();
    
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone) continue;
        
        if (STRUCTGRID == PHMPI::GetZoneGridType(iZone))
        {
            StructGrid *grid = StructGridCast(GetGrid(iZone, level));

            TotalNodes_str += grid->GetNTotalNode();
            TotalFaces_str += grid->GetNTotalFace();
            TotalCells_str += grid->GetNTotalCell();

            minNTCell_str = MIN(minNTCell_str, TotalCells_str);
            maxNTCell_str = MAX(maxNTCell_str, TotalCells_str);
        }
        else
        {
            UnstructGrid *grid = UnstructGridCast(GetGrid(iZone, level));

            TotalNodes_unstr += grid->GetNTotalNode();
            TotalFaces_unstr += grid->GetNTotalFace();
            TotalCells_unstr += grid->GetNTotalCell();

            RDouble *vol = grid->GetCellVolume();
            for (int iCell = 0; iCell < (grid->GetNTotalCell()); ++ iCell)
            {
                TotalVolume_unstr += vol[iCell];
            }
        }
    }
    
    RDouble globalTotalNodes_unstr = TotalNodes_unstr * 1.0;
    RDouble globalTotalFaces_unstr = TotalFaces_unstr * 1.0;
    RDouble globalTotalCells_unstr = TotalCells_unstr * 1.0;
    RDouble globalTotalVolume_unstr = TotalVolume_unstr;
    
    RDouble globalTotalNodes_str = TotalNodes_str * 1.0;
    RDouble globalTotalFaces_str = TotalFaces_str * 1.0;
    RDouble globalTotalCells_str = TotalCells_str * 1.0;

    RDouble minNTCellCompare_str = minNTCell_str * 1.0;
    RDouble maxNTCellCompare_str = maxNTCell_str * 1.0;
    
    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_AllReduce(&TotalNodes_unstr,  &globalTotalNodes_unstr, 1, PH_SUM);
        PH_AllReduce(&TotalFaces_unstr,  &globalTotalFaces_unstr, 1, PH_SUM);
        PH_AllReduce(&TotalCells_unstr,  &globalTotalCells_unstr, 1, PH_SUM);
        PH_AllReduce(&TotalVolume_unstr, &globalTotalVolume_unstr, 1, PH_SUM);
        
        PH_AllReduce(&TotalNodes_str, &globalTotalNodes_str, 1, PH_SUM);
        PH_AllReduce(&TotalFaces_str, &globalTotalFaces_str, 1, PH_SUM);
        PH_AllReduce(&TotalCells_str, &globalTotalCells_str, 1, PH_SUM);

        PH_CompareMaxMin(minNTCellCompare_str, 2);
        PH_CompareMaxMin(maxNTCellCompare_str, 1);
        WriteLogFile("Min && Max nTCell of all zones: ", minNTCellCompare_str, maxNTCellCompare_str);
    }

    RDouble globalTotalNodes = globalTotalNodes_unstr + globalTotalNodes_str;
    RDouble globalTotalFaces = globalTotalFaces_unstr + globalTotalFaces_str;
    RDouble globalTotalCells = globalTotalCells_unstr + globalTotalCells_str;

    GlobalDataBase::UpdateData("GlobalTotalCells", &globalTotalCells, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("UnstrGlobalTotalCells", &globalTotalCells_unstr, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("UnstrGlobalTotalVolume", &globalTotalVolume_unstr, PHDOUBLE, 1);
    
    PrintToWindow("The Global Information of MixGrid ", "\n");
    PrintToWindow("  GlobalTotalNodes: ", globalTotalNodes, "\n");
    PrintToWindow("  GlobalTotalFaces: ", globalTotalFaces, "\n");
    PrintToWindow("  GlobalTotalCells: ", globalTotalCells, "\n");
    
    PrintToWindow("\n");
    
    PrintToWindow("The Global Information of UnstrGrid ", "\n");
    PrintToWindow("  GlobalTotalNodes_unstr: ", globalTotalNodes_unstr, "\n");
    PrintToWindow("  GlobalTotalFaces_unstr: ", globalTotalFaces_unstr, "\n");
    PrintToWindow("  GlobalTotalCells_unstr: ", globalTotalCells_unstr, "\n");
    
    PrintToWindow("\n");
    
    PrintToWindow("The Global Information of StructGrid ", "\n");
    PrintToWindow("  GlobalTotalNodes_str: ", globalTotalNodes_str, "\n");
    PrintToWindow("  GlobalTotalFaces_str: ", globalTotalFaces_str, "\n");
    PrintToWindow("  GlobalTotalCells_str: ", globalTotalCells_str, "\n");
    
    averageVolume_unstr = globalTotalVolume_unstr / globalTotalCells_unstr;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if (STRUCTGRID == PHMPI::GetZoneGridType(iZone))
        {
            continue;
        }
        Zone *zone = GetZone(iZone);
        if (!zone) continue;
        UnstructGrid *grid = UnstructGridCast(GetGrid(iZone, level));
        grid->SetAverageVolume(averageVolume_unstr);
    }
}

void Region::InitMovingGrids(int level)
{
    int nsimutask = GlobalDataBase::GetIntParaFromDB("nsimutask");

    if (nsimutask == CAL_WALL_DIST) return;

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

    if (!isUnsteady)
    {
        return;
    }

    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone) continue;

        Grid *grid = GetGrid(iZone, level);
        grid->InitMovingGrids();
    }
}

LIB_EXPORT void Region::UpdateAllZonesPara()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone) continue;
        zone->UpdateAllData();
    }
}

LIB_EXPORT void Region::SwapGeometry()
{
    CommunicateCellCenterData();

    ComputeMetrics_BC();

    ComputeWeight();

    TestReconstruction();

    WriteLogFile("Swap geommetry sucessfully ...");
}

void Region::CommunicateCellCenterData()
{
    using namespace PHMPI;
    for (int level = 0; level < nGrids; ++ level)
    {
        ActionKey *actkey = new ActionKey();
        actkey->action = COMMCELLCENTERDATA;
        actkey->kind = GRID_BASED;
        actkey->level = level;

        //Task_ServerUpdateInterface *task = new Task_ServerUpdateInterface();
        if (IsNeedNonBlockingCommunication(actkey))
        {
            //task->MainTaskNonBlocking(actkey);
            MainTaskNonBlocking(actkey);
        }
        else
        {
            //task->MainTaskBlocking(actkey);
            MainTaskBlocking(actkey);
        }

        //delete task;
        delete actkey;    actkey = nullptr;
    }
}

void Region::CommunicateCellIBlank()
{
    using namespace PHMPI;
    //for (int level = 0; level < nGrids; ++ level)
    {
        ActionKey *actkey = new ActionKey();
        actkey->action = COMMCELLIBLANK;
        actkey->kind = GRID_BASED;
        actkey->level = 0;

        //Task_ServerUpdateInterface *task = new Task_ServerUpdateInterface();
        if (IsNeedNonBlockingCommunication(actkey))
        {
            //task->MainTaskNonBlocking(actkey);
            MainTaskNonBlocking(actkey);
        }
        else
        {
            //task->MainTaskBlocking(actkey);
            MainTaskBlocking(actkey);
        }

        //delete task;
        delete actkey;    actkey = nullptr;
    }
}

void Region::ComputeMetrics_BC()
{
    int isFVMOrFDM = PHSPACE::GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FD_METHOD)
    {
        ComputeMetricsandJacobinWithBCforStructHighOrder();
        return;
    }

    bool Rwithghost = true;

    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    int nMGLevel = PHSPACE::GlobalDataBase::GetIntParaFromDB("nMGLevel");
    int compressible = GlobalDataBase::GetIntParaFromDB("compressible");
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int ifLowSpeedPrecon = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");
    int wallFunctionType = GlobalDataBase::GetIntParaFromDB("wallFunctionType");

    if (tscheme != LU_SGS) Rwithghost = false;
    if (systemGridType != 1) Rwithghost = false; 
    if (isOverset) Rwithghost = false;
    if (nMGLevel != 1) Rwithghost = false;
    if (compressible != COMPRESSIBLE) Rwithghost = false;
    if (isUnsteady) Rwithghost = false;
    if (nChemical != 0) Rwithghost = false;
    if (ifLowSpeedPrecon != 0) Rwithghost = false;
    if (wallFunctionType != WALLFUNCTION::NONE) Rwithghost = false;

    if (Rwithghost == false) return;

    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = GetZoneProcessorIDSepMode(iZone);
        int myid = PHMPI::GetCurrentProcessorID();
        if (myid == send_proc)
        {
            StructGrid *grid = StructGridCast(GetGrid(iZone, 0));
            grid->ComputeMetrics_BC();
        }
    }
}

void Region::ComputeMetricsandJacobinWithBCforStructHighOrder()
{
    using namespace PHMPI;

    int mission_CommXYZghost        = 1;
    int mission_CommXYZghostofghost = 2;
    int mission_CommMetrics         = 3;
    int mission_CommJacobian        = 4;

    string highordersolvername = GlobalDataBase::GetStrParaFromDB("str_highorder_solver");
    if (highordersolvername.substr(0, 4) == "WCNS")
    {
        CommGridInfoOnlyforStructHighOrder(mission_CommXYZghost);

        CommGridInfoOnlyforStructHighOrder(mission_CommXYZghostofghost);

        UpdateDiffOrderforStructHighOrder();

        UpdateMetricsWithBCforStructHighOrder();

        CommGridInfoOnlyforStructHighOrder(mission_CommMetrics);

        UpdateJacobianWithBCforStructHighOrder();

        CommGridInfoOnlyforStructHighOrder(mission_CommJacobian);
    }
    else if (highordersolvername.substr(0, 4) == "HDCS")
    {
        CommGridInfoOnlyforStructHighOrder(mission_CommXYZghost);
        //! xcc,ycc,zcc,vol
    }
}

void Region::UpdateDiffOrderforStructHighOrder()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if (STRUCTGRID != GetZoneGridType(iZone))
        {
            continue;
        }
        int send_proc = GetZoneProcessorIDSepMode(iZone);

        int myid = PHMPI::GetCurrentProcessorID();
        if (myid == send_proc)
        {
            StructGrid *grid = StructGridCast(GetGrid(iZone, 0));
            grid->InitialDiffOrderforStructHighOrder();

            ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
            int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

            for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
            {
                int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
                
                if (STRUCTGRID != GetZoneGridType(neighborZone))
                {
                    grid->CorrectDiffOrderNearUnstructBlockforStructHighOrder(neighborZone);
                }
            }
        }
    }    
}

void Region::UpdateMetricsWithBCforStructHighOrder()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if (STRUCTGRID != GetZoneGridType(iZone))
        {
            continue;
        }
        int send_proc = GetZoneProcessorIDSepMode(iZone); 
        int myid = PHMPI::GetCurrentProcessorID();
        if (myid == send_proc)
        {
            StructGrid *grid = StructGridCast(GetGrid(iZone, 0));
            grid->ComputeMetricsWithBCforStructHighOrder();
        }
    }
}

void Region::UpdateJacobianWithBCforStructHighOrder()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if (STRUCTGRID != GetZoneGridType(iZone))
        {
            continue;
        }
        int send_proc = GetZoneProcessorIDSepMode(iZone); 
        int myid = PHMPI::GetCurrentProcessorID();
        if (myid == send_proc)
        {
            StructGrid *grid = StructGridCast(GetGrid(iZone, 0));
            grid->ComputeJacobianforStructHighOrder();
        }
    }
}

void Region::CommGridInfoOnlyforStructHighOrder(int mission)
{
    if (GetDim() == TWO_D && mission == 2)
    {
        return;
    }

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector <PH_Request > requestContainer;
    vector <vector <DataContainer *> > receivedDataBuffer;
    vector <vector <DataContainer *> > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    //! Step 0: Compressing data firstly.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int gridtypeofiZone = PHMPI::GetZoneGridType(iZone);

        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if (currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
            
            int gridtypeofneighborZone = PHMPI::GetZoneGridType(neighborZone);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate and compress the buffers for sending.
            DataContainer *sendBuffer = 0;
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone][iNeighbor] = sendBuffer;

                // Compress the send information into the actkey.
                Grid *grid = GetGrid(iZone, 0);

                InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
                if (!interfaceInformation)
                {
                }
                else
                {
                    if (gridtypeofiZone == UNSTRUCTGRID || gridtypeofneighborZone == UNSTRUCTGRID)
                    {
                        RDouble uploadData = 0.0;
                        sendBuffer->MoveToBegin();
                        PHWrite(sendBuffer, uploadData);
                    }
                    else
                    {
                        StructGrid *strgrid = StructGridCast(grid);
                        strgrid->UploadGridInfoOnlyforStructHighOrder(sendBuffer, neighborZone, mission);
                    }
                }
            }
        }
    }

    //! Step 1: Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        //! Communicating.
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
            
            int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[iZone].push_back(receivedBuffer);
            }
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //! Communication.
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->Size();

                    //! Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                //! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                int neighborOrer = -1;
                ZoneNeighbor *globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                int numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (int neighborID = 0; neighborID < numberOfNeighborTemp; ++ neighborID)
                {
                    int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(neighborID);
                    if (neighborZoneTemp == iZone)
                    {
                        neighborOrer = neighborID;
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    //! The data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->Size();

                    //! Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    SWAP(receivedDataBuffer[iZone].back(), sendDataBuffer[iZone][iNeighbor]);
                }
            }
        }
    }

    //! Step2: MPI waiting for Non-blocking communications.
    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    //! Step3: Translate the data container.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int gridtypeofiZone = PHMPI::GetZoneGridType(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
            
            int gridtypeofneighborZone = PHMPI::GetZoneGridType(neighborZone);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                Grid *grid = GetGrid(neighborZone, 0);

                InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
                if (!interfaceInformation)
                {
                    //! can not use return here, because ++count in the after
                }
                else
                {
                    if (gridtypeofiZone == UNSTRUCTGRID || gridtypeofneighborZone == UNSTRUCTGRID)
                    {
                        RDouble downdData = 0.0;
                        receiveData->MoveToBegin();
                        PHRead(receiveData, downdData);
                    }
                    else
                    {
                        StructGrid *strgrid = StructGridCast(grid);
                        strgrid->DownloadGridInfoOnlyforStructHighOrder(receiveData, iZone, mission);
                    }
                }
                ++ count;
            }
        }
    }

    //! Step4: Free the buffers.
    for (int iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (int iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

void Region::TestReconstruction()
{
    using namespace PHMPI;

    ActionKey *actkey = new ActionKey();
    actkey->action = TEST_RECONSTRUCTION;
    actkey->level  = 0;

    int nZones = GetNumberofGlobalZones();

    int myid = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = GetZoneProcessorIDSepMode(iZone);

        if (myid == proc)
        {
            Grid *grid = GetGrid(iZone, actkey->level);
            grid->SkewnessSummary(actkey);
        }
    }
    delete actkey;    actkey = nullptr;
}

void Region::ComputeWeight()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone) continue;
        zone->ComputeWeight();
    }
}

LIB_EXPORT void Region::PostSimulation()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone) continue;
        int Nsolver = GlobalSolvers::GetNumberOfSolvers(iZone);
        for (int iSolver = 0; iSolver < Nsolver; ++ iSolver)
        {
            GlobalSolvers::GetSolver(iZone, iSolver)->Post();
        }
    }

    GlobalSolve();
#ifdef USE_DEMOSOLVER
    return;
#endif
    //! Temporarily, put it here. In the future, add to Post().
    WriteVisualFile();
    //! end
}

LIB_EXPORT void Region::InitGeometry()
{
    SetupMultiBlockRelationship();

    InitGrids();

    MultiGridSetup();

    GlobalGridStatic();

    DumpWallFaceCenter();

    //! Map the OLD bcType to a NEW bcType.
    ChangeBoundaryCondition();

    WriteLogFile("Change boundary condition ...");

    DataMonitor();

    InitVariableWallTemperature();

    InitOriginalGrid();

    WriteLogFile("Monitoring probes ...");
}

//! If need, change the default boundary condition to other type.
LIB_EXPORT void Region::ChangeBoundaryCondition()
{
    int nMGLevel = 1;
    GlobalDataBase::GetData("nMGLevel", &nMGLevel, PHINT, 1);

    string fileName = "bin/BCMap.hypara";
    ifstream infile(fileName.c_str(), ios::in);
    if (!infile) return;

    set< int > changedBC;

    Ustar2Fantasy bcMap(fileName);

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iLevel = 0; iLevel < nMGLevel; ++ iLevel)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            Zone *zone = GetZone(iZone);
            if (!zone) continue;
            Grid *grid = GetGrid(iZone, iLevel);
            int gridtype = grid->Type();
            if (gridtype == STRUCTGRID)
            {
                continue;
            }
            UnstructGrid *unstructGrid = UnstructGridCast(grid);

            UnstructBCSet *unstructBCSet = unstructGrid->GetUnstructBCSet();
            int nBCRegion = unstructBCSet->GetnBCRegion();
            for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
            {
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);

                int originalBC = bcRegion->GetBCType();
                int newBC = bcMap.TranslateBC(originalBC);
                if (newBC != originalBC)
                {
                    changedBC.insert(originalBC);
                }
                bcRegion->SetBCType(newBC);

            }
        }
    }

    if (changedBC.size())
    {
        cout << "The Number of Changed Boundary Condition Type is: " << changedBC.size() << endl;
        int count = 1;
        for (set< int >::iterator iter = changedBC.begin(); iter != changedBC.end(); ++ iter)
        {
            int originalBC = (*iter);
            int newBC = bcMap.TranslateBC(originalBC);
            cout << "Changed BC " << count ++ << ": " << originalBC << " --> " << newBC << endl;
        }
    }
}

void Region::DumpProcessorGridNumber(void)
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    int number_of_processor = GetNumberOfProcessor();
    int nZones = GetNumberofGlobalZones();
    int *nZonesInProc = GetNZonesInProcessor();

    if (! nZonesInProc)
    {
        nZonesInProc = new int [number_of_processor];
    }

    int *numberOfCellInProc = new int [number_of_processor];
    for (int iProc = 0; iProc < number_of_processor; ++ iProc)
    {
        nZonesInProc[iProc] = 0;
        numberOfCellInProc[iProc] = 0;
    }

    int *block_proc = GetZoneProcessorID();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = this->GetZone(iZone);
        if (!zone) continue;
        int nTotalCell = zone->GetGeometry()->GetGrid(0) -> GetNTotalCell();
        int processor = block_proc[iZone];
        numberOfCellInProc[processor] += nTotalCell;
        ++ nZonesInProc[processor];
    }
    
    for (int iProc = 0; iProc < number_of_processor; ++ iProc)
    {
        if (iProc == myid)
        {
            ostringstream oss;
            oss << "In Processor " << iProc << ", Number of Zones = " << nZonesInProc[iProc] <<
                ", Number of Cells = " << numberOfCellInProc[iProc] << endl;
            WriteLogFile(oss);
        }
    }

    delete [] numberOfCellInProc;
    delete [] nZonesInProc;
}

void Region::ReSetBoundaryConditionByGlobalBC()
{
    vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if (globalBCList == 0) return;

    int nMGLevel = 1;
    if (GetTaskCode() == SOLVE_FIELD)
    {
        nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");
    }
    int nZones = PHMPI::GetNumberofGlobalZones();

    for (int iLevel = 0; iLevel < nMGLevel; ++ iLevel)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            Zone *zone = GetZone(iZone);
            if (!zone) continue;

            Grid *grid = GetGrid(iZone, iLevel);

            grid->ReSetBoundaryConditionByGlobalBC();
        }
    }
}

void Region::GridValidityCheck()
{
    using namespace PHMPI;

    RDouble minBox[3] = {LARGE, LARGE, LARGE};
    RDouble maxBox[3] = {-LARGE, -LARGE, -LARGE};

    int nLocalZones = PHMPI::GetNumberofLocalZones();
    int myid = PHMPI::GetCurrentProcessorID();

    //! Check OpenMP threads information firstly.
#ifdef USE_OMP
    int systemGridType = GlobalDataBase::GetIntParaFromDB("sys_gridtype");

    if (systemGridType== MIXGRID || systemGridType== STRUCTGRID)
    {
        TK_Exit::ExceptionExit("The USE_OMP is only for UNSTRUCTGRID!\n");
    }
#ifdef WIN32
    omp_set_num_threads(4);
#endif
#pragma omp parallel for

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        if (myid == 0)
        {
            ostringstream oss;
            oss << "OMP Info: zone ID = " << zoneID << ", thread ID = " << omp_get_thread_num()
                << " of total thread num = " << omp_get_num_threads() << ", MPI ID = "
                << myid << " of total processor num = " << omp_get_num_procs() << endl;
            PrintToWindow(oss);
        }
    }
#endif

    //! Checking and computing global computational fields bounds box.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        int level = 0;
        Grid *grid = GetGrid(zoneID, level);
        grid->GridValidityCheck();

        RDouble *gridMinBox = grid->GetMinBox();
        RDouble *gridMaxBox = grid->GetMaxBox();
        minBox[0] = MIN(minBox[0], gridMinBox[0]);
        minBox[1] = MIN(minBox[1], gridMinBox[1]);
        minBox[2] = MIN(minBox[2], gridMinBox[2]);

        maxBox[0] = MAX(maxBox[0], gridMaxBox[0]);
        maxBox[1] = MAX(maxBox[1], gridMaxBox[1]);
        maxBox[2] = MAX(maxBox[2], gridMaxBox[2]);
    }

    RDouble globalMinBox[3], globalMaxBox[3];
    PHSPACE::SetField(globalMinBox,  PHSPACE::LARGE, 3);
    PHSPACE::SetField(globalMaxBox, -PHSPACE::LARGE, 3);
    if (GetTaskCode() == CAL_WALL_DIST)
    {
        //! For wall distance computing, each processor need know the box.
        PH_AllReduce(minBox, globalMinBox, 3, MPI_MIN);
        PH_AllReduce(maxBox, globalMaxBox, 3, MPI_MAX);
    }
    else
    {
       if (IsReadWalldist() == true)
        {
            //! For general CFD computing, only server processor need know the box.
            PH_Reduce(minBox, globalMinBox, 3, MPI_MIN);
            PH_Reduce(maxBox, globalMaxBox, 3, MPI_MAX);
        }
        else
        {
            //! For wall distance computing, each processor need know the box.
            PH_AllReduce(minBox, globalMinBox, 3, MPI_MIN);
            PH_AllReduce(maxBox, globalMaxBox, 3, MPI_MAX);
        }
    }

    SetGlobalMinMaxBox(globalMinBox, globalMaxBox);

    if (PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
    {
        cout << "Server Processor ID: " << PHMPI::GetServerProcessorID() << endl;
        ostringstream oss;
        oss << setiosflags(ios::right);
        oss << setprecision(6);
        oss << setiosflags(ios::scientific);
        oss << setiosflags(ios::showpoint);
        int wordwidth = 8;
        oss << "Minimum and maximum boundary box of computational field:" << endl;
        oss << "    xMin = " << setw(wordwidth) << globalMinBox[0] << ",  xMax = " << setw(wordwidth) << globalMaxBox[0] << "\n";
        oss << "    yMin = " << setw(wordwidth) << globalMinBox[1] << ",  yMax = " << setw(wordwidth) << globalMaxBox[1] << "\n";
        oss << "    zMin = " << setw(wordwidth) << globalMinBox[2] << ",  zMax = " << setw(wordwidth) << globalMaxBox[2] << "\n";

        PrintToWindow(oss);
    }
}

}
