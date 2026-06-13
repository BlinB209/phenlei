#include "IO_GridReader.h"
#include "Math_BasisFunction.h"
#include "Geo_UnstructBC.h"
#include "IO_FileName.h"
#include "IO_HDF5File.h"
#include "TK_Time.h"
#include "PHIO.h"
#include "GridType.h"

namespace PHSPACE
{
using namespace PHMPI;

LIB_EXPORT IO_GridReader::IO_GridReader(PHVectorString1D & gridGroupNameListIn, int dimension)
{
    globalFileNumber   = 0;
    numberOfTotalZones = 0;
    presentGridGroup   = 0;
    this->gridGroupNameList = gridGroupNameListIn;
    gridGroup = 0;
    this->dimension = dimension;

    myProcID = GetCurrentProcessorID();
    serverProcID = 0;
    if (myProcID == serverProcID)
    {
        this->fileID = 0;
    }
    else
    {
        this->fileID = -1;
    }
    nZonesList = new int [ PHMPI::GetNumberOfProcessor() + PHMPI::GetNumberOfGridProcessor() ];

    cellnodeFileExist = false;
    interpointInfoFileExist = false;
    faceBCFileExist   = false;
    walldistFileExist = false;
}

LIB_EXPORT IO_GridReader::~IO_GridReader()
{
    //delete gridGroup;
    delete [] nZonesList;
}

void IO_GridReader::InitGridFileNameListAndFileProcessor()
{
    //! Global file number is determined by gridGroupNameList[0].
    //! Each grid group should have the same number of files.
    globalFileNumber = CheckNumberofFiles();

    uint_t numberOfGridGroups = gridGroupNameList.size();

    int fileProcessorInterval = PHMPI::GetIntervalofProcessorFiles();

    bool isFileProcessor = (myProcID % fileProcessorInterval == 0 && myProcID < globalFileNumber * fileProcessorInterval);
    
    if (isFileProcessor)
    {
        //! This processor is the file server.
        this->fileID = PHMPI::GetFileIndexofCurrentProcessor();

        //! Set file processor index.
        for (int iGridGroup = 0; iGridGroup < numberOfGridGroups; ++ iGridGroup)
        {
            string gridfileProc = PHSPACE::AddSymbolToFileName(gridGroupNameList[ iGridGroup ], "_", this->fileID);

            gridGroupNameList[ iGridGroup ] = gridfileProc;
        }
    }
    else
    {
        //! This processor is NOT file server.
        this->fileID = -1;
    }

    WriteLogFile("File ID of current processor: ", fileID);

    //! Check whether the cellToNodeFile or the faceBCFile exits.
    cellnodeFileExist = true;

    interpointInfoFileExist = true;

    faceBCFileExist   = true;

    walldistFileExist = true;

    if (myProcID == serverProcID)
    {
        for (int iGridGroup = 0; iGridGroup < numberOfGridGroups; ++ iGridGroup)
        {
            for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
            {
                string cellToNodeFileName = ChangeExtensionOfFileName(gridGroupNameList[ iGridGroup ], "c2n");

                if (! FileExist(cellToNodeFileName))
                {
                    cellnodeFileExist = false;
                    break;
                }
            }

            for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
            {
                string interpointInfoFileName = ChangeExtensionOfFileName(gridGroupNameList[iGridGroup], "interpoint");
                
                if (! FileExist(interpointInfoFileName))
                {
                    interpointInfoFileExist = false;
                    break;
                }
            }

            for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
            {
                string faceBCFileName = ChangeExtensionOfFileName(gridGroupNameList[ iGridGroup ], "bc");

                if (! FileExist(faceBCFileName))
                {
                    faceBCFileExist = false;
                    break;
                }
            }

            for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
            {
                string walldistFileName = ChangeExtensionOfFileName(gridGroupNameList[ iGridGroup ], "wdt");

                if (! FileExist(walldistFileName))
                {
                    walldistFileExist = false;
                    break;
                }
            }
        }
    }

    PH_BcastToClinet(& cellnodeFileExist, sizeof(bool), serverProcID, 2);

    PH_BcastToClinet(& interpointInfoFileExist, sizeof(bool), serverProcID, 2);

    PH_BcastToClinet(& faceBCFileExist,   sizeof(bool), serverProcID, 2);

    PH_BcastToClinet(& walldistFileExist, sizeof(bool), serverProcID, 2);

}

void IO_GridReader::ReadNumberOfTotalZonesAndInitialGlobalValue()
{
    numberOfTotalZones = 0;

    int numberOfProcessor = GetNumberOfProcessor();

    for (unsigned int iGridGroup = 0; iGridGroup < gridGroupNameList.size(); ++ iGridGroup)
    {
        fstream file;

        int nZones = 0;

        Open(file, gridGroupNameList[ iGridGroup ], ios_base::in|ios_base::binary);

        ReadNumberOfZones(file, nZones);

        numberOfTotalZones += nZones;

        Close(file);
    }

    //! Allocate global zone layout information.
    SetNumberOfGlobalZones   (numberOfTotalZones);
    CreateZoneProcessorIDDump(numberOfTotalZones);
    CreateZoneProcessorID    (numberOfTotalZones);
    CreateZoneGridID         (numberOfTotalZones);
    CreateZoneGridType       (numberOfTotalZones);
    CreateZoneFileID         (numberOfTotalZones);

    if (numberOfTotalZones < numberOfProcessor)
    {
        ostringstream oss;
        oss << "   ---- numberOfProcessor = " << numberOfProcessor << ", numberOfTotalZones = " << numberOfTotalZones << " ----\n"
             << "Abnormally Exit Program!!!\n";
        TK_Exit::ExceptionExit(oss.str());
    }
}

int IO_GridReader::CheckNumberofFiles()
{
    //! P2P model.
    string gridFileName = gridGroupNameList.front();

    //! Set file processor index.
    string gridfileProc = gridFileName;
    gridfileProc = PHSPACE::AddSymbolToFileName(gridfileProc, "_", myProcID);
    ifstream infile(gridfileProc.c_str(), ios_base::in);

    int nFileProcessors = 0;
    if (!infile)
    {
        if (myProcID == serverProcID)
        {
            //! The server processor must be file processor.
            TK_Exit::FileOpenErrorExit(gridfileProc);
        }
    }
    else
    {
        nFileProcessors = 1;

        PHSPACE::CloseFile(infile);
    }

    PH_AllReduce(&nFileProcessors, &globalFileNumber, 1, MPI_SUM);
    
    PHMPI::SetNumberofFiles(globalFileNumber);
    WriteLogFile("Global grid files number: ", globalFileNumber);
    
    return globalFileNumber;
}

bool IO_GridReader::CheckGridConvertionMixGrid()
{
    int nsimutask = GlobalDataBase::GetIntParaFromDB("nsimutask");
    if (nsimutask != 1)
    {
        return false;
    }

    int gridobj = GlobalDataBase::GetIntParaFromDB("gridobj");
    if (gridobj != 1)
    {
        return false;
    }

    int gridtype = GlobalDataBase::GetIntParaFromDB("gridtype");
    if (gridtype != 2)
    {
        return false;
    }

    return true;
}

LIB_EXPORT void IO_GridReader::ReadGridsByP2PMode()
{
    WriteLogFile("Start to read grid ... ");

    //! Init grid groups and files.
    InitGridFileNameListAndFileProcessor();

    //! Read number of total zones and initial global value.
    ReadNumberOfTotalZonesAndInitialGlobalValue();

    //! Build gridGroups. Read in layout and distribute.
    int startZoneIndex = 0;

    for (unsigned int iGridGroup = 0; iGridGroup < gridGroupNameList.size(); ++ iGridGroup)
    {
        ostringstream oss;
        oss << "startZoneIndex = " << startZoneIndex << endl;
        WriteLogFile(oss.str());

        presentGridGroup = iGridGroup;

        gridGroup = new GridGroup(startZoneIndex);

        gridGroups.push_back(gridGroup);

        ReadGridGroup(gridGroupNameList[ iGridGroup ], startZoneIndex);

        startZoneIndex += gridGroup->GetNZones();
    }

    //ReadGrid(gridGroupNameList.front(), startZoneIndex);

    WriteLogFile("End reading grid ... ");
}

//! Read in multi files by CS Mode.
LIB_EXPORT void IO_GridReader::ReadGridsByCSMode()
{
    WriteLogFile("Read in grid(s) by CSMode ... ");

    bool gridconvertionmixgrid = CheckGridConvertionMixGrid();

    int numberOfOriginalProcessor = GetNumberOfProcessor();

    uint_t numberOfGridFiles = gridGroupNameList.size();

    numberOfTotalZones = 0;

    if (myProcID == serverProcID)
    {
        this->fileID = 0;
    }
    else
    {
        this->fileID = -1;
    }

    this->numberOfZones.resize(0);

    for (int iFile = 0; iFile < numberOfGridFiles; ++ iFile)
    {
        fstream file;

        int nZones = 0;

        //Open(file, gridGroupNameList[ iFile ], ios_base::in|ios_base::binary);
        OpenSerialFile(file, gridGroupNameList[ iFile ], ios_base::in|ios_base::binary);

        Read_Bcast(file, &nZones, sizeof(int));

        this->numberOfZones.push_back(nZones);

        numberOfTotalZones += nZones;

        Close(file);
    }

    if (numberOfTotalZones < numberOfOriginalProcessor)
    {
        ostringstream oss;
        cout << "   ---- numberOfProcessor = " << numberOfOriginalProcessor << ", numberOfTotalZones = " << numberOfTotalZones << " ----\n"
             << "Abnormally Exit Program!!!\n";

        TK_Exit::ExceptionExit(oss);
    }

    //! Allocate global zone layout information.
    SetNumberOfGlobalZones   (numberOfTotalZones);
    CreateZoneProcessorIDDump(numberOfTotalZones);
    CreateZoneProcessorID    (numberOfTotalZones);
    CreateZoneGridID         (numberOfTotalZones);
    CreateZoneGridType       (numberOfTotalZones);
    CreateZoneFileID         (numberOfTotalZones);

    if (gridconvertionmixgrid)
    {
        CreatetZoneProcessorID_INP(numberOfTotalZones);
    }

    int startZoneIndex = 0;

    for (int iFile = 0; iFile < numberOfGridFiles; ++ iFile)
    {
        presentGridGroup = iFile;

        fstream file;

        OpenSerialFile(file, gridGroupNameList[ iFile ], ios_base::in|ios_base::binary);

        int nZones = 0;

        Read_Bcast(file, &nZones, sizeof(int));

        gridGroup = new GridGroup(startZoneIndex);

        gridGroups.push_back(gridGroup);

        gridGroup->SetNZones(nZones);

        int * block_proc      = new int [nZones];
        int * block_proc_grid = new int [nZones];
        int * block_proc_dump = new int [nZones];
        int * block_idx       = new int [nZones];
        int * block_type      = new int [nZones];
        int * block_file      = new int [nZones];
        int * block_proc_tmpformix = new int [nZones];

        Read_Bcast(file, block_proc, nZones * sizeof(int));
        Read_Bcast(file, block_idx , nZones * sizeof(int));
        Read_Bcast(file, block_type, nZones * sizeof(int));

        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            block_proc_tmpformix[iZone] = block_proc[iZone];
            block_proc_dump[iZone] = block_proc[iZone];
            block_file[iZone]      = serverProcID;
        }

        gridGroup->SetBlockProc    (block_proc     );
        gridGroup->SetBlockProcGrid(block_proc_grid);
        gridGroup->SetBlockProcDump(block_proc_dump);
        gridGroup->SetBlockIndex   (block_idx      );
        gridGroup->SetBlockType    (block_type     );
        gridGroup->SetBockFileProc (block_file     );

        //! Attribute the grid zone to processor
        //for (int iZone = 0; iZone < nZones; ++ iZone)
        //{
        //    block_proc[ iZone ] = (startZoneIndex + iZone) % numberOfOriginalProcessor;
        //}

        //if (numberOfGridProcessor != 0)
        //{
        //    for (int iZone = 0; iZone < nZones; ++ iZone)
        //    {
        //        block_proc_grid[ iZone ] = (startZoneIndex + iZone) % numberOfGridProcessor + numberOfOriginalProcessor;
        //    }
        //}
        //else
        //{
        //    for (int iZone = 0; iZone < nZones; ++ iZone)
        //    {
        //        block_proc_grid[ iZone ] = -1;
        //    }
        //}

        RedistributeZonesIfNeed(startZoneIndex);

        gridGroup->SetGlobalZoneLayout();

        if (gridconvertionmixgrid)
        {
            int *global_block_proc_inp = GetZoneProcessorID_INP();
            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                int multi_izone = startZoneIndex + iZone;
                global_block_proc_inp[multi_izone] = block_proc_tmpformix[iZone];
            }
        }

        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            ReadOneZoneGrid(file, iZone);
        }

        startZoneIndex += nZones;

        Close(file);

        delete [] block_proc_tmpformix;
    }

    ClassifyGridSystem();

    WriteLogFile("End reading grid ... ");
}

void IO_GridReader::ReadGridGroup(const string &gridFileName, int startZoneIndex)
{
    int nZones = 0;

    fstream file;

    Open(file, gridFileName, ios_base::in|ios_base::binary);

    TimeSpan partTime;

    ostringstream oss;

    ReadNumberOfZones(file, nZones);

    oss << "ReadNumberOfZones:" << partTime.GetTimeSpan() << endl; partTime.ResetTime();

    ReadZoneLayout(file, nZones);
    oss << "ReadZoneLayout:" << partTime.GetTimeSpan() << endl; partTime.ResetTime();

    gridGroup->SetNZones(nZones);

    RedistributeZonesIfNeed(startZoneIndex);

    gridGroup->SetGlobalZoneLayout();

    ClassifyGridSystem();

    //! All of the global zones need to be looped over.
    //! The zone not belong to this processor would be added by 'null'.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ReadOneZoneGrid(file, iZone);
    }

    oss << "ReadZoneGrid:" << partTime.GetTimeSpan() << endl; partTime.ResetTime();

    WriteLogFile(oss);

    Close(file);
}

void IO_GridReader::ReadNumberOfZones(fstream & file, int & nZones)
{
    int nZonesLocal = 0;

    Read(file, &nZonesLocal, sizeof(int));

    int nTProcessor = GetNumberOfProcessor() + GetNumberOfGridProcessor();
    
    //! old method.
    //FillGobalData(nZonesLocal, nZonesList);

    //* New method.
    int * localSizeTemp = new int [ nTProcessor ];

    SetField(localSizeTemp, 1, nTProcessor);

    int ** nZonesListTemp = new int * [ nTProcessor ];

    PHMPI::FillGobalDataCollectBcast(&nZonesLocal, localSizeTemp, nZonesListTemp);
    
    for (int iProc = 0; iProc < nTProcessor; ++ iProc)
    {
        nZonesList[ iProc ] = nZonesListTemp[ iProc ][ 0 ];
    }

    DelPointer2(nZonesListTemp);
    DelPointer(localSizeTemp);

    nZones = 0;
    for (int iProc = 0; iProc < nTProcessor; ++ iProc)
    {
        nZones += nZonesList[ iProc ];
    }

    if (nZones == 0)
    {
        TK_Exit::ExceptionExit("Error: number of zones is 0");
    }
}

void IO_GridReader::ReadZoneLayout(fstream & file, int nZones)
{
    int nTProcessors = PHMPI::GetNumberOfProcessor() + PHMPI::GetNumberOfGridProcessor() ;
    int nZonesLocal = nZonesList[myProcID];

    int * block_proc_local = new int [ nZonesLocal ];
    int * block_proc_dump_local = new int [nZonesLocal];
    int * block_idx_local  = new int [ nZonesLocal ];
    int * block_type_local = new int [ nZonesLocal ];

    for (int iZone = 0; iZone < nZonesLocal; ++ iZone)
    {
        block_proc_local[iZone] = -1;
        block_idx_local[iZone]  = -1;
        block_type_local[iZone] = -1;
    }

    //! Read the local layout information in this file.
    Read(file, block_proc_local, nZonesLocal*sizeof(int));
    Read(file, block_idx_local , nZonesLocal*sizeof(int));
    Read(file, block_type_local, nZonesLocal*sizeof(int));

    //! Compress the three information into ONE array,\n
    //! to save the communication time.
    int * tempArray = new int[3 * nZonesLocal];
    int * tempSize  = new int[nTProcessors];
    int count = 0;

    for (int iZone = 0; iZone < nZonesLocal; ++ iZone)
    {
        tempArray[count++] = block_proc_local[iZone];
        tempArray[count++] = block_idx_local[iZone];
        tempArray[count++] = block_type_local[iZone];
    }
    
    int ** globalData = new int * [nTProcessors];
    for (int iProc = 0; iProc < nTProcessors; ++ iProc)
    {
        globalData[iProc] = 0;
        tempSize[iProc] = nZonesList[iProc] * 3;
    }

    //! Bcast the zone layout information to each other processor.\n
    //! Then the layout information of each processor obtained.
    //PHMPI::FillGobalData(tempArray, tempSize, globalData);
    PHSPACE::PHMPI::FillGobalDataCollectBcast(tempArray, tempSize, globalData);

    //! Compute the global layout information.
    int * block_proc = new int [ nZones ];
    int * block_idx  = new int [ nZones ];
    int * block_type = new int [ nZones ];
    int * block_fileProc = new int [nZones];
    int * block_proc_dump = new int [nZones];

    int count1 = 0;
    for (int iProcessor = 0; iProcessor < nTProcessors; ++ iProcessor)
    {
        int count2 = 0;
        for (int jData = 0; jData < nZonesList[iProcessor]; ++ jData)
        {
            block_proc[count1] = globalData[iProcessor][count2++];
            block_idx[count1]  = globalData[iProcessor][count2++];
            block_type[count1] = globalData[iProcessor][count2++];

            block_fileProc[block_idx[count1]] = iProcessor;
            ++ count1;
        }
    }

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        block_proc_dump[iZone] = block_proc[iZone];
    }

    //! Free memory.
    for (int iProc = 0; iProc < nTProcessors; ++ iProc)
    {
        int *temp = globalData[iProc];
        if (temp) delete [] temp;
    }
    delete [] globalData;
    delete [] tempArray;
    delete [] tempSize;
    delete [] block_proc_local;
    delete [] block_proc_dump_local;
    delete [] block_idx_local;
    delete [] block_type_local;

    if (!gridGroup->IsZoneLayoutInitialized()) 
    {
        gridGroup->SetBlockProc(block_proc);
        gridGroup->SetBlockProcDump(block_proc_dump);
        gridGroup->SetBlockIndex(block_idx);
        gridGroup->SetBlockType(block_type);
        gridGroup->SetBockFileProc(block_fileProc);
    }
    else
    {
        delete [] block_proc;
        delete [] block_proc_dump;
        delete [] block_idx;
        delete [] block_type;
        delete [] block_fileProc;
        return;
    }
}

void IO_GridReader::ClassifyGridSystem()
{
    using namespace PHMPI;

    set<int> grid_type;
    for (int iZone = 0; iZone < numberOfTotalZones; ++ iZone)
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

void IO_GridReader::RedistributeZonesIfNeed(int startZoneIndex)
{
    int nZones = gridGroup->GetNZones();
    int number_of_processor = PHSPACE::GetNumberOfProcessor();
    int * block_type      = gridGroup->GetBlockType();
    int * block_proc      = gridGroup->GetBlockProc();
    int * block_proc_grid = gridGroup->GetBlockProcGrid();
    int * block_proc_dump = gridGroup->GetBlockProcDump();

    int zoneDistributionMethod;
    int maxproc = block_proc[0];
    
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        maxproc = MAX(maxproc, block_proc[iZone]);
    }
    maxproc += 1;
    if (number_of_processor == maxproc)
    {
        zoneDistributionMethod = DETERMIN_BY_PARTITION;
    }else
    {
        zoneDistributionMethod = REDISTRIBUTION;
    }

    PHMPI::SetZoneDistributionMethod(zoneDistributionMethod);

    ostringstream oss;
    oss << "Max processor ID : " << maxproc << " of total zone " << nZones;
    if (zoneDistributionMethod == DETERMIN_BY_PARTITION)
    {
        oss << ", Zones' distribution is specified by partition ...";
    }
    else
    {
        oss << ", Zones' distribution is specified randomly ...";
    }
    oss << endl;
    WriteLogFile(oss);

    if (zoneDistributionMethod == REDISTRIBUTION)
    {
        if (nZones % number_of_processor == 0)
        {
            int nZonesPerProcessor = nZones / number_of_processor;
            int processorID = 0;
            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                block_proc[iZone] = processorID;
                if ((iZone + 1) % nZonesPerProcessor == 0)
                {
                    processorID ++;
                }
            }
        }
        else
        {
            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                block_proc[iZone] = (startZoneIndex + iZone) % number_of_processor;
            }
        }
    }
    else
    {
        if (number_of_processor == 1)
        {
            //! Bell 20131124 add.
            if (IsConvertGridToMixGrid())
            {
                //! For structured + unstructured = mix grid,move the specified process of  unstructured grid partition backward.
                if (block_type[0] == PHSPACE::STRUCTGRID)
                {
                    maxproc = block_proc[0];
                    for (int iZone = 0; iZone < nZones; ++ iZone)
                    {
                        if (block_proc[iZone] > maxproc)
                        {
                            maxproc = block_proc[iZone];
                        }
                    }
                    maxproc += 1;
                    SetNumberOfProcStructUsed(maxproc);
                }
                else
                {
                    for (int iZone = 0; iZone < nZones; ++ iZone)
                    {
                        block_proc_dump[iZone] += GetNumberOfProcStructUsed();
                    }
                }
            }
            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                block_proc[iZone] = 0;
            }
        }
    }

    if (!block_proc_grid)
    {
        block_proc_grid = new int[ nZones ];
    }

    gridGroup->SetBlockProcGrid(block_proc_grid);

    int numberOfGridProcessor = GetNumberOfGridProcessor();

    if (numberOfGridProcessor != 0)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            block_proc_grid[ iZone ] = (startZoneIndex + iZone) % numberOfGridProcessor + number_of_processor;
        }
    }
    else
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            block_proc_grid[ iZone ] = -1;
        }
    }
}

void IO_GridReader::ReadOneZoneGrid(fstream & file, int & iZone)
{
    int send_proc      = 0;
    int recv_proc      = 0;
    int recv_proc_grid = 0;

    GetSendRecvProcID(iZone, send_proc, recv_proc, recv_proc_grid);

    DataContainer * gridData = new DataContainer();

    if (myProcID == send_proc)
    {
        ReadFile(file, gridData);
    }

    PH_Trade(gridData, send_proc, recv_proc);

    if (recv_proc_grid != -1)
    {
        PH_Trade(gridData, send_proc, recv_proc_grid);
    }

    //! Decompress the grid data.
    DecompressDataToGrid(gridData, iZone);

    delete gridData;
}

void IO_GridReader::DecompressDataToGrid(DataContainer *cdata, int iZone)
{
    int send_proc      = 0;
    int recv_proc      = 0;
    int recv_proc_grid = 0;

    GetSendRecvProcID(iZone, send_proc, recv_proc, recv_proc_grid);

    if (!(myProcID == recv_proc || myProcID == recv_proc_grid))
    {
        gridGroup->AddGrid(0);
        return;
    }

    //! Create new grid.
    int * block_type = gridGroup->GetBlockType();
    int gridType = block_type[iZone];

    //! The id of the zone needs to be added with zoneStart(cxh 20180515).
    int zoneStart = gridGroup->GetZoneStart();
    GridID * index = new GridID(iZone + zoneStart);

    Grid *grid = 0;
    if (gridType == PHSPACE::UNSTRUCTGRID)
    {
        grid = new UnstructGrid();
    }
    else
    {
        grid = new StructGrid();
    }

    grid->InitGrid(index, 0, dimension, gridType);

    grid->SetIBlock(presentGridGroup);

    //! Decompress.
    grid->Decode(cdata, 0);

    //! Add zoneStart to grid.
    grid->ReviseNeighborZoneIndex(zoneStart);

    //! Add to grid group.
    gridGroup->AddGrid(grid);
}

void IO_GridReader::GetSendRecvProcID(int iZone, int &send_proc, int &recv_proc)
{
    int *block_proc = PHMPI::GetZoneProcessorID();
    if (! IsParallelRun())
    {
        send_proc = this->serverProcID;
        recv_proc = this->serverProcID;
    }
    else
    {
        int *block_fileProc = PHMPI::GetZoneFileID();
        send_proc = block_fileProc[iZone];
        recv_proc = block_proc[iZone];
    }
}

void IO_GridReader::GetSendRecvProcID(int localZoneID, int &send_proc, int &recv_proc, int &recv_proc_grid)
{
    int * block_proc      = PHMPI::GetZoneProcessorID();
    int * block_proc_grid = PHMPI::GetZoneProcessorIDForGrid();

    recv_proc      = block_proc[localZoneID];
    recv_proc_grid = block_proc_grid[localZoneID];

    if (! IsParallelRun())
    {
        send_proc      = this->serverProcID;
    }
    else
    {
        int * block_fileProc = PHMPI::GetZoneFileID();
        send_proc = block_fileProc[localZoneID];
    }
}

void IO_GridReader::Read(fstream &file, void *data, int size)
{
    if (this->fileID != -1)
    {
        PH_Read(file,data,size);
    }
}

void IO_GridReader::Read_Bcast(fstream &file, void *data, int size)
{
    Read(file, data, size);
    PH_BcastToClinet(data, size, serverProcID, 2);
}

void IO_GridReader::Bcast(void *data, int size, int proc, int tag)
{
    if (myProcID == proc)
    {
        int number_of_processor = GetNumberOfProcessor();
        for (int i = 0; i < number_of_processor; ++ i)
        {
            if (i == proc) continue;
            PHSPACE::PH_Send(data,size,i,tag);
        }
    }
    else
    {
        PHSPACE::PH_Receive(data,size,proc,tag);
    }
}

void IO_GridReader::ReadCompressedData(fstream &file, DataContainer *cdata, int send_proc, int recv_proc, int tag)
{
    using namespace PHMPI;

    if (myProcID == send_proc)
    {
        ReadFile(file, cdata);
    }

    PH_Trade(cdata, send_proc, recv_proc, tag);
}

void IO_GridReader::Open(fstream &file, const string &filename, const ios_base::openmode &openmode)
{
    using namespace PHMPI;

    if (fileID == -1) return;

    OpenFile(file, filename, openmode);
}

void IO_GridReader::OpenSerialFile(fstream &file, const string &filename, const ios_base::openmode &openmode)
{
    string actualFileName = AddSymbolToFileName(filename, "_", 0);
    this->Open(file, actualFileName, openmode);
}

void IO_GridReader::Close(fstream &file)
{
    if (myProcID != fileID) return;

    file.close();
    file.clear();
}

LIB_EXPORT void IO_GridReader::ReadWallDistByCSMode()
{
    bool allFilesExitOrNot = true;

    for (unsigned int iFile = 0; iFile < gridGroupNameList.size(); ++ iFile)
    {
        string gridFileName = gridGroupNameList[ iFile ];
        //string actualFileName = AddSymbolToFileName(gridFileName, "_", 0);
        string cellToNodeFileName = ChangeExtensionOfFileName(gridFileName, "wdt");

        if (! SerialFileExist(cellToNodeFileName))
        {
            allFilesExitOrNot = false;
            break;
        }
    }

    if (! allFilesExitOrNot) 
    {
        return;
    }

    fstream file;
    string gridFileName =ChangeExtensionOfFileName(gridGroupNameList[0], "wdt");
    OpenSerialFile(file, gridFileName, ios_base::in | ios_base::binary);
    int numberOfcells;
    int sizeInt = sizeof(int);
    file.read((char *)&numberOfcells, sizeInt);
    CloseFile(file);

    if (numberOfcells >= static_cast <int> (1e8) || numberOfcells < 0)
    {
        return;
    }

    int zoneStart = 0;

    bool gridIsHdf = false;
    if (this->gridGroups.size() == 0)
    {
        gridIsHdf = true;
        if (this->numberOfZones.size() == 0)
        {
        this->numberOfZones.push_back(PHMPI::GetNumberofGlobalZones());
    }
    }

    for (unsigned int iFile = 0; iFile < gridGroupNameList.size(); ++ iFile)
    {
        if (!gridIsHdf)
        {
            this->gridGroup = this->gridGroups[iFile];
        }

        string gridFileNamei = gridGroupNameList[ iFile ];

        string cellToNodeFileName = ChangeExtensionOfFileName(gridFileNamei, "wdt");

        fstream file1;

        OpenSerialFile(file1, cellToNodeFileName, ios_base::in|ios_base::binary);

        int nZones = this->numberOfZones[iFile];

        for (int iZone = zoneStart; iZone < nZones + zoneStart; ++ iZone)
        {
            int send_proc = 0;
            int recv_proc = 0;

            GetSendRecvProcID(iZone, send_proc, recv_proc);

            DataContainer * WallDistData = new DataContainer();

            ReadCompressedData(file1, WallDistData, send_proc, recv_proc);

            DecompressDataToWallDist(WallDistData, iZone);

            delete WallDistData;
        }

        zoneStart += nZones;

        Close(file1);
    }
}

LIB_EXPORT void IO_GridReader::ReadCellToNodeByCSMode()
{
    bool allFilesExitOrNot = true;

    for (unsigned int iFile = 0; iFile < gridGroupNameList.size(); ++ iFile)
    {
        string gridFileName = gridGroupNameList[ iFile ];
        //string actualFileName = AddSymbolToFileName(gridFileName, "_", 0);
        string cellToNodeFileName = ChangeExtensionOfFileName(gridFileName, "c2n");

        if (! SerialFileExist(cellToNodeFileName))
        {
            allFilesExitOrNot = false;
            break;
        }
    }

    if (! allFilesExitOrNot) 
    {
        return;
    }

    int zoneStart = 0;

    for (unsigned int iFile = 0; iFile < gridGroupNameList.size(); ++ iFile)
    {
        this->gridGroup = this->gridGroups[iFile];

        string gridFileName = gridGroupNameList[ iFile ];

        string cellToNodeFileName = ChangeExtensionOfFileName(gridFileName, "c2n");

        fstream file;

        OpenSerialFile(file, cellToNodeFileName, ios_base::in|ios_base::binary);

        int nZones = 0;

        Read_Bcast(file, & nZones, sizeof(int));

        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int send_proc = 0;
            int recv_proc = 0;
            int recv_proc_grid = 0;

            GetSendRecvProcID(iZone, send_proc, recv_proc, recv_proc_grid);

            DataContainer * c2nData = new DataContainer();

            //! Read the compressed grid data.
            ReadCompressedData(file, c2nData, send_proc, recv_proc);

            //! Decompress the grid data.
            DecompressDataToCellNode(c2nData, iZone);

            delete c2nData;
        }

        zoneStart += nZones;

        Close(file);
    }
}

LIB_EXPORT void IO_GridReader::ReadInterpointInfoByCSMode()
{
    bool allFilesExitOrNot = true;

    for (unsigned int iFile = 0; iFile < gridGroupNameList.size(); ++ iFile)
    {
        string gridFileName = gridGroupNameList[iFile];

        string interpointInfoFileName = ChangeExtensionOfFileName(gridFileName, "interpoint");

        if (! SerialFileExist(interpointInfoFileName))
        {
            allFilesExitOrNot = false;
            break;
        }
    }

    if (! allFilesExitOrNot)
    {
        return;
    }

    int zoneStart = 0;

    bool gridIsHdf = false;
    if (numberOfZones.size() == 0)
    {
        gridIsHdf = true;
    }

    for (unsigned int iFile = 0; iFile < gridGroupNameList.size(); ++ iFile)
    {
        if (!gridIsHdf)
        {
            this->gridGroup = this->gridGroups[iFile];
        }

        string gridFileName = gridGroupNameList[iFile];

        string interpointInfoFileName = ChangeExtensionOfFileName(gridFileName, "interpoint");

        fstream file;

        OpenSerialFile(file, interpointInfoFileName, ios_base::in|ios_base::binary);

        int nZones = 0;

        Read_Bcast(file, & nZones, sizeof(int));

        if (gridIsHdf)
        {
            this->numberOfZones.push_back(nZones);
        }

        for (int iZone = 0; iZone < nZones + zoneStart; ++ iZone)
        {
            int send_proc = 0;
            int recv_proc = 0;
            int recv_proc_grid = 0;

            GetSendRecvProcID(iZone, send_proc, recv_proc, recv_proc_grid);

            DataContainer *interpointData = new DataContainer();

            //! Read the compressed grid data.
            ReadCompressedData(file, interpointData, send_proc, recv_proc);

            //! Decompress the grid data.
            DecompressInterpointInfo(interpointData, iZone);

            delete interpointData;
        }

        zoneStart += nZones;

        Close(file);
    }

    zoneConnectivityForPoint = new ZoneConnectivityForPoint();
}

LIB_EXPORT void IO_GridReader::ReadWallDistByP2PMode()
{
    if (!walldistFileExist) return;

    string gridFileName, WallDistFileName;

    if (fileID != -1)
    {
        gridFileName = gridGroupNameList.front();
        WallDistFileName = ChangeExtensionOfFileName(gridFileName, "wdt");
    }

    //! Check if the walldist file exist.
    if (fileID != -1 && !FileExist(WallDistFileName)) return;

    int wdtFileIsH5 = 0;
    if (fileID != -1)
    {
        hid_t fileId;
        fileId = OpenHDF5File(WallDistFileName);

        if (!(fileId < 0))
        {
            H5Fclose(fileId);
            wdtFileIsH5 = 1;
        }
    }

    PH_CompareMaxMin(wdtFileIsH5, 1);
    if (wdtFileIsH5)
    {
        return;
    }

    fstream file;
    Open(file, WallDistFileName, ios_base::in|ios_base::binary);

    //Number of zones in each .fts file, 1 is default.
    int nZones = PHMPI::GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = 0;
        int recv_proc = 0;
        GetSendRecvProcID(iZone, send_proc, recv_proc);
        
        DataContainer * walldistData = new DataContainer();

        ReadCompressedData(file, walldistData, send_proc, recv_proc);
        DecompressDataToWallDist(walldistData, iZone);
        
        delete walldistData;
    }

    Close(file);
}

LIB_EXPORT void IO_GridReader::ReadCellToNodeByP2PMode()
{
    if (!cellnodeFileExist) return;
    
    string gridFileName, cellToNodeFileName;
    if (fileID != -1)
    {
        gridFileName = gridGroupNameList.front();
        cellToNodeFileName = ChangeExtensionOfFileName(gridFileName, "c2n");
    }
    
    //! Check if the c2n file exist.
    if (fileID != -1 && !FileExist(cellToNodeFileName)) return;
    
    fstream file;
    Open(file, cellToNodeFileName, ios_base::in|ios_base::binary);

    int nZones = 0;
    Read_Bcast(file, &nZones, sizeof(int));

    //if (nZones > 1) return;
    if (nZones == 0)
    {
        nZones = this->numberOfTotalZones;
    }

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = 0;
        int recv_proc = 0;
        GetSendRecvProcID(iZone, send_proc, recv_proc);
        
        DataContainer * c2nData = new DataContainer();

        //! Read the compressed grid data.
        ReadCompressedData(file, c2nData, send_proc, recv_proc);

        //! Decompress the grid data.
        DecompressDataToCellNode(c2nData, iZone);

        delete c2nData;
    }

    Close(file);
}

LIB_EXPORT void IO_GridReader::ReadInterpointInfoByP2PMode()
{
    if (!interpointInfoFileExist)
    {
        return;
    }
    string gridFileName, interpointInfoFileName;
    if (fileID != -1)
    {
        gridFileName = gridGroupNameList.front();
        interpointInfoFileName = ChangeExtensionOfFileName(gridFileName, "interpoint");
    }
    
    //! Check if the c2n file exist.
    if (fileID != -1 && ! FileExist(interpointInfoFileName))
    {
        return;
    }
    
    fstream file;
    Open(file, interpointInfoFileName, ios_base::in|ios_base::binary);

    int nZones = 0;
    Read_Bcast(file, &nZones, sizeof(int));

    //if (nZones > 1) return;
    if (nZones == 0)
    {
        nZones = this->numberOfTotalZones;
    }
    
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = 0;
        int recv_proc = 0;
        GetSendRecvProcID(iZone, send_proc, recv_proc);
        
        DataContainer *interpointInformation = new DataContainer();

        //! Read the compressed grid data.
        ReadCompressedData(file, interpointInformation, send_proc, recv_proc);

        //! Decompress the grid data.
        DecompressInterpointInfo(interpointInformation, iZone);

        delete interpointInformation;
    }
    Close(file);

    zoneConnectivityForPoint = new ZoneConnectivityForPoint();
}


LIB_EXPORT void IO_GridReader::ReadFaceBCByCSMode()
{
    bool allFilesExitOrNot = true;

    for (unsigned int iFile = 0; iFile < gridGroupNameList.size(); ++ iFile)
    {
        string gridFileName = gridGroupNameList[ iFile ];
        //string actualFileName = AddSymbolToFileName(gridFileName, "_", 0);
        string faceBCFileName = ChangeExtensionOfFileName(gridFileName, "bc");

        if (!SerialFileExist(faceBCFileName))
        {
            allFilesExitOrNot = false;
            break;
        }
    }

    if (!allFilesExitOrNot)
    {
        return;
    }

    int zoneStart = 0;

    for (unsigned int iFile = 0; iFile < gridGroupNameList.size(); ++ iFile)
    {
        this->gridGroup = this->gridGroups[iFile];

        string gridFileName = gridGroupNameList[ iFile ];

        string faceBCFileName = ChangeExtensionOfFileName(gridFileName, "bc");

        fstream file;

        OpenSerialFile(file, faceBCFileName, ios_base::in|ios_base::binary);

        int nZones = 0;

        Read_Bcast(file, & nZones, sizeof(int));

        for (int iZone = zoneStart; iZone < nZones + zoneStart; ++ iZone)
        {
            int send_proc = 0;
            int recv_proc = 0;
            int recv_proc_grid = 0;

            GetSendRecvProcID(iZone, send_proc, recv_proc, recv_proc_grid);

            DataContainer * faceBCData = new DataContainer();

            //! Read the compressed grid data.
            ReadCompressedData(file, faceBCData, send_proc, recv_proc);

            //! Decompress the grid data.
            DecompressDataToFaceBC(faceBCData, iZone);

            delete faceBCData;
        }

        zoneStart += nZones;

        Close(file);
    }
}

LIB_EXPORT void IO_GridReader::ReadFaceBcDirByCSMode()
{
    bool allFilesExitOrNot = true;

    for (unsigned int iFile = 0; iFile < gridGroupNameList.size(); ++ iFile)
    {
        string gridFileName = gridGroupNameList[iFile];
        //string actualFileName = AddSymbolToFileName(gridFileName, "_", 0);
        string faceBCFileName = ChangeExtensionOfFileName(gridFileName, "bcdir");

        if (!SerialFileExist(faceBCFileName))
        {
            allFilesExitOrNot = false;
            break;
        }
    }

    if (!allFilesExitOrNot) 
    {
        return;
    }

    int zoneStart = 0;

    for (unsigned int iFile = 0; iFile < gridGroupNameList.size(); ++ iFile)
    {
        string gridFileName = gridGroupNameList[iFile];
        //string actualFileName = AddSymbolToFileName(gridFileName, "_", 0);
        string faceBCFileName = ChangeExtensionOfFileName(gridFileName, "bcdir");

        fstream file;

        OpenSerialFile(file, faceBCFileName, ios_base::in|ios_base::binary);

        int nZones = 0;

        Read_Bcast(file, &nZones, sizeof(int));

        for (int iZone = 0; iZone < nZones + zoneStart; ++ iZone)
        {
            int send_proc = 0;
            int recv_proc = 0;
            int recv_proc_grid = 0;

            GetSendRecvProcID(iZone, send_proc, recv_proc, recv_proc_grid);

            DataContainer * faceBcDirData = new DataContainer();

            //! Read the compressed grid data.
            ReadCompressedData(file, faceBcDirData, send_proc, recv_proc);

            //! Decompress the grid data.
            DecompressDataToFaceBcDir(faceBcDirData, iZone);

            delete faceBcDirData;
        }

        zoneStart += nZones;

        Close(file);
    }
}

LIB_EXPORT void IO_GridReader::ReadFaceBCByP2PMode()
{
    if (!faceBCFileExist) return;

    string gridFileName, faceBCFileName;
    if (fileID != -1)
    {
        gridFileName = gridGroupNameList.front();
        faceBCFileName = ChangeExtensionOfFileName(gridFileName, "bc");
    }

    //! Check if the face BC file exist.
    if (fileID != -1 && !FileExist(faceBCFileName)) return;

    fstream file;
    Open(file, faceBCFileName, ios_base::in|ios_base::binary);

    int nZones = 0;
    Read_Bcast(file, &nZones, sizeof(int));
    //if (nZones > 1) return;
    if (nZones == 0)
    {
        nZones = this->numberOfTotalZones;
    }

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = 0;
        int recv_proc = 0;
        GetSendRecvProcID(iZone, send_proc, recv_proc);

        DataContainer * faceBCData = new DataContainer();

        //! Read the compressed grid data.
        ReadCompressedData(file, faceBCData, send_proc, recv_proc);

        //! Decompress the grid data.
        DecompressDataToFaceBC(faceBCData, iZone);

        delete faceBCData;
    }

    Close(file);
}

LIB_EXPORT void IO_GridReader::ReadFaceBcDirByP2PMode()
{
    if (!faceBCFileExist) return;

    string gridFileName, faceBCFileName;
    if (fileID != -1)
    {
        gridFileName = gridGroupNameList.front();
        faceBCFileName = ChangeExtensionOfFileName(gridFileName, "bcdir");
    }

    //! Check if the face BC file exist.
    if (fileID != -1 && !FileExist(faceBCFileName)) return;

    fstream file;
    Open(file, faceBCFileName, ios_base::in|ios_base::binary);

    int nZones = 0;
    Read_Bcast(file, &nZones, sizeof(int));
    //if (nZones > 1) return;
    if (nZones == 0)
    {
        nZones = this->numberOfTotalZones;
    }

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = 0;
        int recv_proc = 0;
        GetSendRecvProcID(iZone, send_proc, recv_proc);

        DataContainer * faceBcDirData = new DataContainer();

        //! Read the compressed grid data.
        ReadCompressedData(file, faceBcDirData, send_proc, recv_proc);

        //! Decompress the grid data.
        DecompressDataToFaceBcDir(faceBcDirData, iZone);

        delete faceBcDirData;
    }

    Close(file);
}

void IO_GridReader::DecompressDataToWallDist(DataContainer *cdata, int iZone)
{
    int send_proc = 0;
    int recv_proc = 0;

    GetSendRecvProcID(iZone, send_proc, recv_proc);

    if (myProcID != recv_proc)
    {
        return;
    }

    int zoneStart = 0;
    if (gridGroups.size() != 0)
    {
        zoneStart = this->gridGroup->GetZoneStart();
    }

    Grid *grid = PHSPACE::GetGrid(iZone + zoneStart, 0);

    int * block_type = PHMPI::GetZoneGridType();
    int gridType = block_type[iZone];

    if (gridType == PHSPACE::UNSTRUCTGRID)
    {
        UnstructGrid * unstructuredGrid = UnstructGridCast(grid);
        unstructuredGrid->InitWallDist();
        unstructuredGrid->DecodeWallDist(cdata);
    }
    else
    {
        StructGrid * structuredgrid = StructGridCast(grid);
        structuredgrid->InitWallDist();
        structuredgrid->DecodeWallDist(cdata);
    }
}

void IO_GridReader::DecompressDataToCellNode(DataContainer *cdata, int iZone)
{
    int send_proc = 0;
    int recv_proc = 0;

    GetSendRecvProcID(iZone, send_proc, recv_proc);

    if (myProcID != recv_proc)
    {
        return;
    }

    //! Decode the read-in data to cell-node.
    int zoneStart = this->gridGroup->GetZoneStart();
    Grid *grid = PHSPACE::GetGrid(iZone + zoneStart, 0);
    UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

    unstructuredGrid->DecodeCellNode(cdata);
}

void IO_GridReader::DecompressInterpointInfo(DataContainer *cdata, int iZone)
{
    int send_proc = 0;
    int recv_proc = 0;

    GetSendRecvProcID(iZone, send_proc, recv_proc);

    if (myProcID != recv_proc)
    {        
        return;
    }

    //! Decode the read-in data to cell-node.
    int zoneStart = 0;
    if (gridGroups.size() != 0)
    {
        zoneStart = this->gridGroup->GetZoneStart();
    }
    Grid *grid = PHSPACE::GetGrid(iZone + zoneStart, 0);
    UnstructGrid * unstructuredGrid = UnstructGridCast(grid);    

    unstructuredGrid->DecodeInterpointInfo(cdata);
}

void IO_GridReader::DecompressDataToFaceBC(DataContainer *cdata, int iZone)
{
    int send_proc = 0;
    int recv_proc = 0;

    GetSendRecvProcID(iZone, send_proc, recv_proc);

    if (myProcID != recv_proc)
    {
        return;
    }

    //! Decode the data into face BC information.
    //int zoneStart = this->gridGroup->GetZoneStart();
    //Grid *grid = PHSPACE::GetGrid(iZone + zoneStart, 0);

    Grid *grid = PHSPACE::GetGrid(iZone, 0);
    if (grid->Type() == UNSTRUCTGRID)
    {
        UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

        unstructuredGrid->DecodeBCFace(cdata);
    }
    else
    {
        StructGrid *structuredGrid = StructGridCast(grid);

        structuredGrid->DecodeBCFace(cdata);
    }
}

void IO_GridReader::DecompressDataToFaceBcDir(DataContainer *cdata, int iZone)
{
    int send_proc = 0;
    int recv_proc = 0;

    GetSendRecvProcID(iZone, send_proc, recv_proc);

    if (myProcID != recv_proc)
    {
        return;
    }

    //! Decode the data into face BC information.
    Grid *grid = PHSPACE::GetGrid(iZone, 0);
    if (grid->Type() == UNSTRUCTGRID)
    {
        return;
    }
    else
    {
        StructGrid *structuredGrid = StructGridCast(grid);

        structuredGrid->DecodeBcDirFace(cdata);
    }
}

}