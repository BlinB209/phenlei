#ifdef PH_PARALLEL
    #include "mpi.h"
#endif

#include <stack>
#include "PHMpi.h"
#include "Constants.h"
#include "Math_BasisFunction.h"
#include "GlobalDataBase.h"

using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{

namespace PHMPI
{
//! Kernel Global Variables.
//! To control the whole PHengLEI software architecture,
//! especially for parallel computational.
//! The zones (grid blocks are including) information is defined here.
//! -# myID               :  Index of the current processor.
//! -# numberOfProcessor  :  The total number of processors used by this procedure.
//! -# server             :  Index of the server processor, default 0.
//! -# nZones             :  The total number of global zones that used in parallel approach.
//!                          It is equal to the sum of zones in all processors.
//! Zones lay out information:
//! -# block_proc :  The processor index of each zone(grid).
//! -# block_type :  The grid type of each zone(grid), structured or unstructured.
//! -# block_idx  :  The zone(block) index of each zone(grid).
//! -# block_file :  The processor index that manage data in each zone(grid).
int csmode     = 0;    // Initialization()

int server     = 0;    // Initialization()
int myID       = 0;    // Initialization()
int gridServer = 0;
int nZones     = 0;     //! The number of global zones.

//! Some special zones information for HYBRID SOLVER.
int *nZones_in_proc = 0;    // Initialization()

VVInt zone_over_link;

int numberOfProcessor          = 1;    // Initialization()
int numberOfGridProcessor      = 0;
int numberOfZonesInThisProcess = 0;
int numberofFiles              = 0;    // Initialization()

int *block_proc = NULL;               // Initialization()
int *block_proc_grid = NULL;
int *block_type = NULL;               // Initialization()
int *block_idx  = NULL;               // Initialization()
int *block_file = NULL;
int *block_proc_inp = NULL;           // Initialization()

int *block_proc_dump = NULL;

int *block_group = NULL;

//! Is current processor is file processor:
//! -#  -1: has not been judged.
//! -#  0 : no.
//! -#  1 : yes.
int isFileServer = -1;

//! zoneDistributionMethod: zone block distribution method, default 0.
//! -# 0 : each zone is assigned to the one that defined in grid partition procedure.
//! -# 1 : random assigned for each zone or by some else ways.
int zoneDistributionMethod = 0;

vector<int> *localZoneIDToGlobalZoneID = NULL;

//! Processor tree.
PHProcessorTree *processorTree = NULL;
PHProcessorTree *processorTreeForAllProcessor  = NULL;
PHProcessorTree *processorTreeForGridProcessor = NULL;

int numberOfProcessorStructgridUsed = 0;

int *zoneStartPointContainer = NULL;
int *zoneStartCenterContainer = NULL;

DataStruct_BinaryTree *binaryTree = NULL;

void Initialization()
{
    myID       = 0;
    server     = 0;
    nZones_in_proc = 0;
    block_proc = 0;
    block_type = 0;
    block_idx  = 0;
    block_proc_inp = 0;

    //! Default model.
    csmode     = 0;

#ifdef PH_PARALLEL
    int argc = 0;
    char ***argv = 0;

    PH_Init(argc, argv);
    PH_Rank(myID);

    int numberOfProcessorTmp;
    PH_Size(numberOfProcessorTmp);
    SetNumberOfProcessor(numberOfProcessorTmp);
    PHMPI::BuildProcessorTree();
    string proc_name;
    PH_GetProcessorName(proc_name);

    ostringstream oss;
    oss << "Process " << myID << " numberOfProcessor = " << numberOfProcessorTmp << " name = " << proc_name << endl;

    if (GetCurrentProcessorID() == server)
    {
        std::cout << oss.str();
    }

    oss.clear();
    oss.str("");

    PHMPI::SetNumberofFiles(1);
#endif
}

void SetServerProcessorID(int serverId)
{
    server = serverId;
}

void SetGridServerProcessorID(int serverId)
{
    gridServer = serverId;
}

LIB_EXPORT void SetNumberOfGlobalZones(int nblock)
{
    nZones = nblock;
}

int GetParallelStrategy()
{
    return PHSPACE::GlobalDataBase::GetIntParaFromDB("parallelStrategy");
}

void SetParallelStrategy(int parallelStrategy)
{
    PHSPACE::GlobalDataBase::UpdateData("parallelStrategy", &parallelStrategy, PHINT, 1);
}

void SetNumberOfProcessor(int nproc)
{
    numberOfProcessor = nproc;
}

void SetNumberOfGridProcessor(int nproc)
{
    numberOfGridProcessor = nproc;
}

LIB_EXPORT void SetNumberofLocalZones(int nblock)
{
    numberOfZonesInThisProcess = nblock;
}

void SetNumberofFiles(int fileNumber)
{
    numberofFiles = fileNumber;
}

LIB_EXPORT void SetZoneProcessorID(int iZone, int id)
{
    block_proc[iZone] = id;
}

void SetZoneProcessID_INP(int *data)
{
    block_proc_inp = data;
}

void SetIsFileServer(int isServer)
{
    isFileServer = isServer;
}

void SetZoneDistributionMethod(int method)
{
    zoneDistributionMethod = method;
}

void InsertGlobalZones(int globalZoneID)
{
    if (!localZoneIDToGlobalZoneID)
    {
        localZoneIDToGlobalZoneID = new vector <int>;
    }

    localZoneIDToGlobalZoneID->push_back(globalZoneID);
}

void BuildProcessorTree()
{
    int numberOfOriginalProcessor = GetNumberOfProcessor();
    RDouble start = 0.0;
    RDouble end = static_cast<RDouble>(numberOfOriginalProcessor);

    int dim = 1;
    processorTree = new PHProcessorTree(dim, &start, &end);

    for (int iProc = 0; iProc < numberOfOriginalProcessor; ++ iProc)
    {
        RDouble coord = static_cast<RDouble>(iProc);
        PHProcessorNode *node = new PHProcessorNode(dim, &coord, iProc);
        processorTree->AddNode(node);
    }

    int numberOfGridProcessorTmp  = GetNumberOfGridProcessor();
    int numberOfTotalProcessor = numberOfGridProcessorTmp + numberOfOriginalProcessor;

    if (numberOfGridProcessorTmp == 0)
    {
        processorTreeForAllProcessor  = processorTree;
        processorTreeForGridProcessor = processorTree;
        return;
    }

    //! processorTreeForAllProcessor.
    start = 0;
    end = numberOfTotalProcessor;
    processorTreeForAllProcessor = new PHProcessorTree(dim, &start, &end);

    for (int iProc = 0; iProc < numberOfTotalProcessor; ++ iProc)
    {
        RDouble coord = static_cast<RDouble>(iProc);
        PHProcessorNode *node = new PHProcessorNode(dim, &coord, iProc);
        processorTreeForAllProcessor->AddNode(node);
    }

    //! processorTreeForGridProcessor.
    start = 0;
    end = numberOfGridProcessorTmp;
    processorTreeForGridProcessor = new PHProcessorTree(dim, &start, &end);

    for (int iProc = 0; iProc < numberOfGridProcessorTmp; ++ iProc)
    {
        RDouble coord = static_cast<RDouble>(iProc);
        int item = iProc + numberOfOriginalProcessor;

        PHProcessorNode *node = new PHProcessorNode(dim, &coord, item);
        processorTreeForGridProcessor->AddNode(node);
    }
}

void SetNumberOfProcStructUsed(int number)
{
    numberOfProcessorStructgridUsed = number;
}

int GetServerProcessorID()
{
    return server;
}

LIB_EXPORT int GetCurrentProcessorID()
{
    return myID;
}

int GetGridServerProcessorID()
{
    return gridServer;
}

LIB_EXPORT int GetNumberofGlobalZones()
{
    return nZones;
}

int * GetNZonesInProcessor()
{
    return nZones_in_proc;
}

LIB_EXPORT int GetNumberOfProcessor()
{
    return numberOfProcessor;
}

LIB_EXPORT int GetNumberOfGridProcessor()
{
    return numberOfGridProcessor;
}

LIB_EXPORT int GetNumberofLocalZones()
{
    return numberOfZonesInThisProcess;
}

LIB_EXPORT int GetNumberofFiles()
{
    return numberofFiles;
}

LIB_EXPORT int * GetZoneProcessorID()
{
    return block_proc;
}

LIB_EXPORT int GetZoneProcessorID(int iZone)
{
    return block_proc[iZone];
}

LIB_EXPORT int * GetZoneProcessorIDForGrid()
{
    return block_proc_grid;
}

LIB_EXPORT int GetZoneProcessorIDForGrid(int iZone)
{
    return block_proc_grid[iZone];
}

LIB_EXPORT int * GetZoneGridType()
{
    return block_type;
}

LIB_EXPORT int GetZoneGridType(int iZone)
{
    return block_type[iZone];
}

LIB_EXPORT int * GetZoneGridID()
{
    return block_idx;
}

LIB_EXPORT int * GetZoneFileID()
{
    return block_file;
}

LIB_EXPORT int GetZoneFileID(int iZone)
{
    return block_file[iZone];
}

LIB_EXPORT int * GetZoneProcessorID_INP()
{
    return block_proc_inp;
}

LIB_EXPORT int * GetZoneProcessorIDDump()
{
    return block_proc_dump;
}

LIB_EXPORT int * GetZoneGroupID()
{
    return block_group;
}

LIB_EXPORT int GetZoneDistributionMethod()
{
    return zoneDistributionMethod;
}

LIB_EXPORT int GetLocalZoneIDToGlobalZoneID(int localZoneID)
{
    if (!localZoneIDToGlobalZoneID)
    {
        ostringstream oss;
        oss << "  Error: this situation has not been considered, for " << "local zone ID" << " = " << localZoneID << endl;

        PHMPI::FreeBlockData();
        cout << oss.str() << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    return (*localZoneIDToGlobalZoneID)[localZoneID];
}

LIB_EXPORT int GetGlobalZoneIDToLocalZoneID(int globalZoneIndex)
{
    int localZoneIndex = -1;

    int numberOfZones = PHMPI::GetNumberofLocalZones();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        if (GetLocalZoneIDToGlobalZoneID(iZone) == globalZoneIndex)
        {
            localZoneIndex = iZone;
            break;
        }
    }

    if (localZoneIndex == -1)
    {
        cout << "Error: Could not find out the local zone index !" << endl;
    }

    return localZoneIndex;
}


int GetNumberOfProcStructUsed()
{
    return numberOfProcessorStructgridUsed;
}

void CreateZoneProcessorID(int nZonesTmp)
{
    if (block_proc)
    {
        delete [] block_proc;
        //cout << "  Warning: block_proc array has been created, it would be deleted here !\n";
    }
    block_proc = new int [nZonesTmp];

    if (block_proc_grid)
    {
        delete [] block_proc_grid;
        //cout << "  Warning: block_proc array has been created, it would be deleted here !\n";
    }
    block_proc_grid = new int [nZonesTmp];
}

void CreatetZoneProcessorID_INP(int nZonesTmp)
{
    if (block_proc_inp)
    {
        delete [] block_proc_inp;
        //cout << "  Warning: block_proc_inp array has been created, it would be deleted here !\n";
    }

    block_proc_inp = new int [nZonesTmp];
}

void CreateZoneGridID(int nZonesTmp)
{
    if (block_idx)
    {
        delete [] block_idx;
        //cout << "  Warning: block_idx array has been created, it would be deleted here !\n";
    }
    block_idx = new int [nZonesTmp];
}

void CreateZoneGridType(int nZonesTmp)
{
    if (block_type)
    {
        delete [] block_type;
        //cout << "  Warning: block_type array has been created, it would be deleted here !\n";
    }

    block_type = new int [nZonesTmp];
}

void CreateZoneFileID(int nZonesTmp)
{
    if (block_file)
    {
        delete [] block_file;
        //cout << "  Warning: block_file array has been created, it would be deleted here !\n";
    }

    block_file = new int [nZonesTmp];
}

void CreateZoneProcessorIDDump(int nZonesTmp)
{
    if (block_proc_dump)
    {
        delete [] block_proc_dump;
        //cout << "Warning: CreateZoneProcessorIDDumphas been created, it would be deleted here !\n";
    }
    block_proc_dump = new int [nZonesTmp];
}

void CreateZoneGroupID(int nZonesTmp)
{
    if (block_group)
    {
        delete [] block_group;
        //cout << "Warning: CreateZoneGroupID has been created, it would be deleted here !\n";
    }
    block_group = new int [nZonesTmp];
}

void FreeBlockData()
{
    if (block_proc)
    {
        delete [] block_proc;
    }
    block_proc = NULL;

    if (block_type)
    {
        delete [] block_type;
    }
    block_type = NULL;

    if (block_idx)
    {
        delete [] block_idx;
    }
    block_idx = NULL;

    if (block_file)
    {
        delete [] block_file;
    }
    block_file = NULL;

    if (block_proc_inp)
    {
        delete [] block_proc_inp;
    }
    block_proc_inp = NULL;

    if (nZones_in_proc)
    {
        delete [] nZones_in_proc;
    }
    nZones_in_proc = NULL;

    if (block_proc_dump)
    {
        delete [] block_proc_dump;
    }
    block_proc_dump = NULL;

    if (localZoneIDToGlobalZoneID)
    {
        delete localZoneIDToGlobalZoneID;
    }
    localZoneIDToGlobalZoneID = NULL;
}

void CreateZoneStartPointContainer(int nZonesTmp)
{
    using namespace PHMPI;

    zoneStartPointContainer = new int [nZonesTmp];
}

int *GetZoneStartPointContainer()
{
    return zoneStartPointContainer;
}

void DeleteZoneStartPointContainer()
{
    delete [] zoneStartPointContainer;
    zoneStartPointContainer = NULL;
}

void CreateZoneStartCenterContainer(int nZonesTmp)
{
    using namespace PHMPI;

    zoneStartCenterContainer = new int [nZonesTmp];
}

int * GetZoneStartCenterContainer()
{
    return zoneStartCenterContainer;
}

void DeleteZoneStartCenterContainer()
{
    delete [] zoneStartCenterContainer;
    zoneStartCenterContainer = NULL;
}

void GenerateGlobalBinaryTree(int *zoneStartPointContainerTmp, int nZonesTmp, int numberOfPoints)
{
    binaryTree = new DataStruct_BinaryTree(zoneStartPointContainerTmp, nZonesTmp, numberOfPoints);
    return;
}

int GetZoneIndexAccordingToGlobalPointLabel(int globalPointLabel)
{
    return binaryTree->ComputeZoneLabel(globalPointLabel);
}

void DeleteGlobalBinaryTree()
{
    if (binaryTree)
    {
        delete binaryTree;    binaryTree = NULL; 
    }
    return;
}

LIB_EXPORT bool IsParallelRun()
{
    return numberOfProcessor > 1;
}

LIB_EXPORT bool CurrentProcessorIsGridProcessor()
{
    int numberOfTotalProcessor = numberOfProcessor + numberOfGridProcessor;

    if (myID >= numberOfProcessor && myID < numberOfTotalProcessor)
    {
        return true;
    }
    else
    {
        return false;
    }
}

LIB_EXPORT int GetIntervalofProcessorFiles()
{
    return PHMPI::GetNumberOfProcessor() / PHMPI::GetNumberofFiles();
}

LIB_EXPORT int GetFileIndexofCurrentProcessor()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    int fileProcessorInterval = PHMPI::GetIntervalofProcessorFiles();
    int fileID = currentProcessorIndex / fileProcessorInterval;
    return fileID;
}

LIB_EXPORT int GetZoneProcessorIDSepMode(int iZone)
{
    if (CurrentProcessorIsGridProcessor())
    {
        return block_proc_grid[iZone];
    }
    else
    {
        return block_proc[iZone];
    }
}

LIB_EXPORT void GetSendRecvProc(int iZone, int &sendProcessor, int &recvProcessor)
{
    sendProcessor = block_file[iZone];
    recvProcessor = block_proc[iZone];
}

int GetServerSepMode()
{
    if (CurrentProcessorIsGridProcessor())
    {
        return gridServer;
    }
    else
    {
        return server;
    }
}

LIB_EXPORT int * GetZoneProcessorIDSepMode()
{
    if (CurrentProcessorIsGridProcessor())
    {
        return block_proc_grid;
    }
    else
    {
        return block_proc;
    }
}

LIB_EXPORT void GetFileProc(int &file_proc)
{
    if (!IsParallelRun())
    {
        file_proc = server;
        if (myID == PHMPI::GetServerProcessorID())
        {
            SetIsFileServer(1);
        }
    }
    else
    {
        //file_proc = myID;
        if (!block_file)
        {
            //! The file_proc has not been initialized.
            file_proc = myID;
            SetIsFileServer(1);
        }
        else
        {
            //! The file_proc has been initialized.
            file_proc = -1;
            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                if (block_file[iZone] == myID)
                {
                    file_proc = myID;
                    break;
                }
            }

            if (file_proc == -1)
            {
                SetIsFileServer(0);
            }
            else
            {
                SetIsFileServer(1);
            }
        }
    }
}

LIB_EXPORT int IsCurrentProcessorFileServer()
{
    if (isFileServer >= 0)
    {
        return isFileServer;
    }
    else
    {
        int fileProcessorIndex;
        GetFileProc(fileProcessorIndex);

        if (GetCurrentProcessorID() == fileProcessorIndex)
            return 1;
        else
            return 0;
    }
}

void PH_Read(fstream &file, void *data, int size, int proc)
{
    if (myID == proc)
    {
        PH_Read(file, data, size);
    }
}

void PH_Read(fstream &file, void *data, int size)
{
    file.read(reinterpret_cast<char *>(data), size);
}

void PH_ReadBcastInteger(fstream &file, void *data, int size, int proc)
{
    PH_Read_Bcast(file, data, size * sizeof(int), proc);
}

void PH_ReadBcastDouble(fstream &file, void *data, int size, int proc)
{
    PH_Read_Bcast(file, data, size * sizeof(double), proc);
}

void PH_Read_Bcast(fstream &file, void *data, int size, int proc)
{
    PH_Read(file, data, size, proc);

    PH_Bcast(data, size, proc);
}

void PH_BcastToClinet(void *data, int size, int proc, int key, int tag)
{
    int numberOfOriginalProcessor = GetNumberOfProcessor();
    int numberOfGridProcessorTmp     = GetNumberOfGridProcessor();

    int numberOfProcessorTmp = numberOfOriginalProcessor;
    int startProcessor = 0;

    //! 0: to original processors; 1: to grid processors; 2: to all processors.
    if (key == 0)
    {
        numberOfProcessorTmp = numberOfOriginalProcessor;
        startProcessor = 0;
    }
    else if (key == 1)
    {
        numberOfProcessorTmp = numberOfOriginalProcessor + numberOfGridProcessorTmp;
        startProcessor = numberOfOriginalProcessor;
    }
    else if (key == 2)
    {
        numberOfProcessorTmp = numberOfOriginalProcessor + numberOfGridProcessorTmp;
        startProcessor = 0;
    }

    bool keyReceive = true;

    if (key == 0 && myID >= numberOfOriginalProcessor)
    {
        keyReceive = false;
    }

    if (key == 1 && myID < numberOfOriginalProcessor)
    {
        keyReceive = false;
    }

    if (myID == proc)
    {
        for (int i = startProcessor; i < numberOfProcessorTmp; ++ i)
        {
            if (i == proc) continue;
            PHSPACE::PH_Send(data, size, i, tag);
        }
    }
    else
    {
        if (keyReceive)
        {
            PHSPACE::PH_Receive(data, size, proc, tag);
        }
    }
}

void PH_BcastSepMode(void *data, int size, int proc, int tag)
{
    int startProcessor = 0;
    int endProcessor = GetNumberOfProcessor();

    if (CurrentProcessorIsGridProcessor())
    {
        startProcessor = GetNumberOfProcessor();
        endProcessor = startProcessor + GetNumberOfGridProcessor();
    }

    if (myID == proc)
    {
        for (int i = startProcessor; i < endProcessor; ++ i)
        {
            if (i == proc) continue;
            PHSPACE::PH_Send(data, size, i, tag);
        }
    }
    else
    {
        PHSPACE::PH_Receive(data, size, proc, tag);
    }
}

void PH_Bcast(void *data, int size, int proc, int tag)
{
    if (myID == proc)
    {
        int numberOfProcessorTmp = GetNumberOfProcessor();
        for (int i = 0; i < numberOfProcessorTmp; ++ i)
        {
            if (i == proc)
            {
                continue;
            }
            PHSPACE::PH_Send(data, size, i, tag);
        }
    }
    else
    {
        PHSPACE::PH_Receive(data, size, proc, tag);
    }
}

void PH_Bcast(DataContainer *cdata, int proc, int tag)
{
    if (myID == proc)
    {
        int numberOfProcessorTmp = GetNumberOfProcessor();
        for (int i = 0; i < numberOfProcessorTmp; ++ i)
        {
            if (i == proc)
            {
                continue;
            }
            send(cdata, i, tag);
        }
    }
    else
    {
        receive(cdata, proc, tag);
    }
}

void PH_BcastByServer(DataContainer *cdata, int tag)
{
    queue <PHProcessorNode *> nodeQueue;
    PHProcessorNode *start = processorTree->GetRoot();

    if (!start)
    {
        return;
    }

    if (start->GetData() != GetServerProcessorID())
    {
        ostringstream oss;
        oss << "\n";
        oss << "++++++++++++++++++  StopProgram Information  +++++++++++++++++++++++++++\n";
        oss <<  "Error: the server should be processor 0!" << "\n";
        oss << " The stop filename is : " << __FILE__ << "\n";
        oss << " at line " << __LINE__ << "\n";
        oss << " Compiled On line " << __DATE__ << " at " << __TIME__ << "\n";
        oss << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        PHMPI::FreeBlockData();
        cout << oss.str() << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    nodeQueue.push(start);
    while (!nodeQueue.empty())
    {
        start = nodeQueue.front();
        nodeQueue.pop();

        int sendProcessor = start->GetData();

        if (start->left)
        {
            nodeQueue.push(start->left);
            int recvProcessor = start->left->GetData();

            if (sendProcessor != recvProcessor)
            {
                if (myID == sendProcessor)
                {
                    send(cdata, recvProcessor, tag);
                }
                if (myID == recvProcessor)
                {
                    receive(cdata, sendProcessor, tag);
                }
            }
        }

        if (start->right)
        {
            nodeQueue.push(start->right);

            int recvProcessor = start->right->GetData();

            if (sendProcessor != recvProcessor)
            {
                if (myID == sendProcessor)
                {
                    send(cdata, recvProcessor, tag);
                }
                if (myID == recvProcessor)
                {
                    receive(cdata, sendProcessor, tag);
                }
            }
        }
    }
}

void PH_GatherByServer(DataContainer *cdata, int tag)
{
    stack <PHProcessorNode *> nodeStack;
    PHProcessorNode *start = processorTree->GetRoot();
    PHProcessorNode *cur = start, *pre = NULL;
    if (!start)
    {
        return;
    }

    if (start->GetData() != GetServerProcessorID())
    {
        ostringstream oss;
        oss << "\n";
        oss << "++++++++++++++++++  StopProgram Information  +++++++++++++++++++++++++++\n";
        oss <<  "Error: the server should be processor 0!" << "\n";
        oss << " The stop filename is : " << __FILE__ << "\n";
        oss << " at line " << __LINE__ << "\n";
        oss << " Compiled On line " << __DATE__ << " at " << __TIME__ << "\n";
        oss << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        PHMPI::FreeBlockData();
        cout << oss.str() << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    nodeStack.push(start);

    while (!nodeStack.empty())
    {
        cur = nodeStack.top();

        int recvProcessor = cur->GetData();

        if ((cur->left == NULL && cur->right == NULL) || ((pre == cur->left || pre == cur->right) && pre != NULL))
        {
            nodeStack.pop();
            pre = cur;
        } 
        else
        {
            if (cur->right)
            {
                nodeStack.push(cur->right);
                int sendProcessor = cur->right->GetData();

                if (sendProcessor != recvProcessor)
                {
                    if (myID == sendProcessor)
                    {
                        send(cdata, recvProcessor, tag);
                    }
                    if (myID == recvProcessor)
                    {
                        receive(cdata, sendProcessor, tag);
                    }
                }
            }
            if (cur->left)
            {
                nodeStack.push(cur->left);
                int sendProcessor = cur->right->GetData();

                if (sendProcessor != recvProcessor)
                {
                    if (myID == sendProcessor)
                    {
                        send(cdata, recvProcessor, tag);
                    }
                    if (myID == recvProcessor)
                    {
                        receive(cdata, sendProcessor, tag);
                    }
                }
            }
        }
    }
}

void PH_BcastByServerSepMode(DataContainer *cdata, int tag)
{
    queue <PHProcessorNode *> nodeQueue;

    PHProcessorNode *start = processorTree->GetRoot();

    if (CurrentProcessorIsGridProcessor())
    {
        start = processorTreeForGridProcessor->GetRoot();
    }

    if (!start)
    {
        return;
    }

    if (start->GetData() != GetServerSepMode())
    {
        ostringstream oss;
        oss << "\n";
        oss << "++++++++++++++++++  StopProgram Information  +++++++++++++++++++++++++++\n";
        oss <<  "Error: the server should be processor 0!" << "\n";
        oss << " The stop filename is : " << __FILE__ << "\n";
        oss << " at line " << __LINE__ << "\n";
        oss << " Compiled On line " << __DATE__ << " at " << __TIME__ << "\n";
        oss << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        PHMPI::FreeBlockData();
        cout << oss.str() << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    nodeQueue.push(start);
    while (!nodeQueue.empty())
    {
        start = nodeQueue.front();
        nodeQueue.pop();

        int sendProcessor = start->GetData();

        if (start->left)
        {
            nodeQueue.push(start->left);
            int recvProcessor = start->left->GetData();

            if (sendProcessor != recvProcessor)
            {
                if (myID == sendProcessor)
                {
                    send(cdata, recvProcessor, tag);
                }
                if (myID == recvProcessor)
                {
                    receive(cdata, sendProcessor, tag);
                }
            }
        }

        if (start->right)
        {
            nodeQueue.push(start->right);

            int recvProcessor = start->right->GetData();

            if (sendProcessor != recvProcessor)
            {
                if (myID == sendProcessor)
                {
                    send(cdata, recvProcessor, tag);
                }
                if (myID == recvProcessor)
                {
                    receive(cdata, sendProcessor, tag);
                }
            }
        }
    }
}

void PH_BcastByServerToClinet(DataContainer *cdata, int key, int tag)
{
    queue <PHProcessorNode *> nodeQueue;

    PHProcessorNode *start = processorTree->GetRoot();

    if (key == 1)
    {
        start = processorTreeForGridProcessor->GetRoot();
    }
    else if (key == 2)
    {
        start = processorTreeForAllProcessor->GetRoot();
    }

    if (!start)
    {
        return;
    }

    if (start->GetData() != GetServerProcessorID())
    {
        ostringstream oss;
        oss << "\n";
        oss << "++++++++++++++++++  StopProgram Information  +++++++++++++++++++++++++++\n";
        oss <<  "Error: the server should be processor 0!" << "\n";
        oss << " The stop filename is : " << __FILE__ << "\n";
        oss << " at line " << __LINE__ << "\n";
        oss << " Compiled On line " << __DATE__ << " at " << __TIME__ << "\n";
        oss << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        PHMPI::FreeBlockData();
        cout << oss.str() << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    nodeQueue.push(start);
    while (!nodeQueue.empty())
    {
        start = nodeQueue.front();
        nodeQueue.pop();

        int sendProcessor = start->GetData();

        if (start->left)
        {
            nodeQueue.push(start->left);
            int recvProcessor = start->left->GetData();

            if (sendProcessor != recvProcessor)
            {
                if (myID == sendProcessor)
                {
                    send(cdata, recvProcessor, tag);
                }
                if (myID == recvProcessor)
                {
                    receive(cdata, sendProcessor, tag);
                }
            }
        }

        if (start->right)
        {
            nodeQueue.push(start->right);

            int recvProcessor = start->right->GetData();

            if (sendProcessor != recvProcessor)
            {
                if (myID == sendProcessor)
                {
                    send(cdata, recvProcessor, tag);
                }
                if (myID == recvProcessor)
                {
                    receive(cdata, sendProcessor, tag);
                }
            }
        }
    }
}

void PH_BcastNonBlocking(DataContainer *&cdata, int bcastProcessor, int tag)
{
    using namespace PHMPI;

    vector<DataContainer *> sendDataBuffer;
    vector<DataContainer *> receivedDataBuffer;
    vector<PH_Request> requestContainer;

    int sendProcessor = bcastProcessor;
    int myIDTmp = GetCurrentProcessorID();
    int numberOfProcessorTmp = GetNumberOfProcessor();

    CharVecSizeType dataLength;
    if (myIDTmp == bcastProcessor)
    {
        dataLength = cdata->Size();
    }
    PH_BcastNonBlocking(&dataLength, 1, bcastProcessor, 0);

    if (myIDTmp == bcastProcessor)
    {
        sendDataBuffer.reserve(numberOfProcessorTmp);
        requestContainer.reserve(numberOfProcessorTmp);

        //! Step 1: Server sending.
        for (int iProc = 0; iProc < numberOfProcessorTmp; ++ iProc)
        {
            if (iProc == bcastProcessor)
            {
                continue;
            }

            int receiveProcessor = iProc;

            //! Allocate the buffers for sending and receiving.
            DataContainer *sendBuffer = new DataContainer (cdata);
            sendDataBuffer.push_back(sendBuffer);

            //! Send the data to others.
            send(sendBuffer, receiveProcessor, dataLength, requestContainer, tag);
        }
    }
    else
    {
        DataContainer *receivedBuffer = new DataContainer();
        receivedDataBuffer.push_back(receivedBuffer);

        //! Receive the data from server.
        receive(receivedBuffer, sendProcessor, dataLength, requestContainer, tag);
    }

    //! Step2: MPI waiting for Non-blocking communications.
    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    //! Step3: Translating.
    if (myIDTmp != bcastProcessor)
    {
        if (cdata) delete cdata;    cdata = nullptr;

        cdata = receivedDataBuffer[0];
    }

    //! Step4: Free the buffers.
    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        delete sendDataBuffer[iDim];
    }

    requestContainer.clear();
    receivedDataBuffer.clear();
    sendDataBuffer.clear();
}

void PH_BcastNlenData(int &nlen, char *&data, int proc, int tag)
{
    PH_Bcast(&nlen, sizeof(int), proc, tag);

    if (myID != proc)
    {
        data = new char [nlen];
    }

    PH_Bcast(data, nlen, proc, tag);
}

//! proc process obtains data and sends them to all processes, other processes receive data and then deal with them.
void PH_BcastComposite(DATA_COMPRESS dt_c, DATA_DECOMPRESS dt_d, int proc, int tag)
{
    DataContainer *cdata = new DataContainer();

    if (myID == proc)
    {
        //! Compress the data to cdata.
        dt_c(cdata);
    }

    //! Send the cdata to other processes which need it.
    PH_Bcast(cdata, proc, tag);

    if (myID != proc)
    {
        //! Decompress data of the cdata to obtain the required information.
        dt_d(cdata);
    }

    delete cdata;    cdata = nullptr;
}

void PH_BcastString(string &cs, int processor, int tag)
{
    if (myID == processor)
    {
        int numberOfProcessorTmp = GetNumberOfProcessor();
        for (int i = 0; i < numberOfProcessorTmp; ++ i)
        {
            if (i == processor)
            {
                continue;
            }
            PH_SendString(cs, processor, tag);
        }
    }
    else
    {
        PH_ReceiveString(cs, processor, tag);
    }
}

void PH_CollectString(string &cs, int processor, int tag)
{
    if (myID != processor)
    {
        //! Processes other than proc are responsible for sending information to proc process.
        PH_SendString(cs, processor, tag);
        //WriteLogFile(cs);
    }
    else
    {
        //! myID == proc, all information is received and handled by proc process.
        int numberOfProcessorTmp = GetNumberOfProcessor();
        for (int i = 0; i < numberOfProcessorTmp; ++ i)
        {
            if (i != processor)
            {
                PH_ReceiveString(cs, i, tag);
            }
            if (GetCurrentProcessorID() == server)
            {
                cout << cs;
            }
        }
    }
}

void ServerRead(fstream &file, void *data, int size)
{
    PH_Read(file, data, size, server);
    PH_Bcast(data, size, server);
}

void PH_Trade(void *data, int size, int sendProcessor, int recvProcessor, int tag)
{
    //! Do nothing if the sendProcessID equal to the recieveProcessID.
    if (sendProcessor == recvProcessor) return;

    if (myID == sendProcessor)
    {
        //! Send the data to the receive process.
        PHSPACE::PH_Send(data, size, recvProcessor, tag);
    }
    else if (myID == recvProcessor)
    {
        //! Receive the data from the send process.
        PHSPACE::PH_Receive(data, size, sendProcessor, tag);
    }
}

void PH_Trade(DataContainer *cdata, int sendProcessor, int recvProcessor, int tag)
{
    if (sendProcessor == recvProcessor) return;

    if (myID == sendProcessor)
    {
        send(cdata, recvProcessor, tag);
    }
    else if (myID == recvProcessor)
    {
        receive(cdata, sendProcessor, tag);
    }
}

void PH_Trade(DataContainer *cdata, int sendProcessor, int recvProcessor, streamsize nlen, int tag)
{
    if (sendProcessor == recvProcessor) return;

    if (myID == sendProcessor)
    {
        send(cdata, recvProcessor, static_cast<CharVecSizeType>(nlen), tag);
    }
    if (myID == recvProcessor)
    {
        receive(cdata, sendProcessor, static_cast<CharVecSizeType>(nlen), tag);
    }
}

void PH_Trade(ActionKey *actkey, int sendProcessor, int recvProcessor, int tag)
{
    DataContainer *cdata = actkey->GetData();

    PH_Trade(cdata, sendProcessor, recvProcessor, tag);
}

void PH_Trade(ActionKey *actkey, int sendProcessor, int recvProcessor, streamsize nlen, int tag)
{
    DataContainer *cdata = actkey->GetData();
    PH_Trade(cdata, sendProcessor, recvProcessor, nlen, tag);
}

}

bool IsNeedNonBlockingCommunication(const ActionKey *actkey)
{
    using namespace PHMPI;

    if (actkey->action == UPDATE_INTERFACE_DATA)
    {
        return true;
    }
    //if (actkey->action == COMMCELLIBLANK)
    //{
    //    return true;
    //}
    if (actkey->action == UPDATE_INTERPOINT_DATA)
    {
        return true;
    }
    return false;
}

void StreamToActionKey(ActionKey *actkey, ostringstream &oss)
{
    if (actkey)
    {
        DataContainer *cdata = actkey->GetData();
        cdata->MoveToBegin();
        cdata->Write(&oss);
    }
}

void DataCharToString(DataContainer *cdata, string &str)
{
    cdata->ToString(str);
}

int GetSendRecvTag(ActionKey *actkey, int iZone)
{
    using namespace PHMPI;
    int nZonesTmp = GetNumberofGlobalZones();

    int tag = 0;

    if (nZonesTmp <= 0)
    {
        //! It is easy to make mistakes here.
        //! nZones is not assigned here.
        PHMPI::FreeBlockData();
        cout << "FATAL ERROR in GetSendRecvTag\n" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (actkey->kind == SOLVER_BASED)
    {
        tag = iZone + actkey->solverID * nZonesTmp;
    }
    else
    {
        tag = iZone + actkey->level * nZonesTmp;
    }

    return tag;
}

int GetSendRecvTag(int tagKind, int index, int iZone)
{
    using namespace PHMPI;
    int nZonesTmp = GetNumberofGlobalZones();

    int tag = 0;

    if (nZonesTmp <= 0)
    {
        //! It is easy to make mistakes here.
        //! nZones is not assigned here.
        PHMPI::FreeBlockData();
        cout << "FATAL ERROR in GetSendRecvTag\n" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    tag = iZone + index * nZonesTmp;

    return tag;
}

void PH_AllreduceSepMode(int *sendingBuffer, int *receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op)
{
    using namespace PHMPI;

#ifdef PH_PARALLEL
    int numberOfGridProcessorTmp = GetNumberOfGridProcessor();

    //ostringstream oss;
    //oss << "sendingBuffer = " << sendingBuffer[ 0 ] << endl;
    //WriteLogFile(oss.str());

    if (numberOfGridProcessorTmp == 0)
    {
        MPI_Allreduce(sendingBuffer, receivingBuffer, numberOfElements, mpiDataType, op, MPI_COMM_WORLD);
        return;
    }

    int serverTmp = PHMPI::GetServerSepMode();

    int numberOfOriginalProcessor = GetNumberOfProcessor();

    int startProcessor = 0;
    int endProcessor = numberOfOriginalProcessor;

    if (CurrentProcessorIsGridProcessor())
    {
        startProcessor = numberOfOriginalProcessor;
        endProcessor = numberOfGridProcessorTmp + numberOfOriginalProcessor;
    }

    int myIDTmp = GetCurrentProcessorID();

    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        receivingBuffer[iElement] = sendingBuffer[iElement];
    }

    for (int iProc = startProcessor; iProc < endProcessor; ++ iProc)
    {
        if (iProc == serverTmp)
        {
            continue;
        }

        int tag = iProc;

        if (myIDTmp == iProc)
        {
            PH_SmartSend(sendingBuffer, numberOfElements, serverTmp, tag);
        }

        if (myIDTmp == serverTmp)
        {
            int *tmpBuffer = new int[numberOfElements];

            PH_SmartRecv(tmpBuffer, numberOfElements, iProc, tag);

            if (op == MPI_MAX)
            {
                for (int iElement = 0; iElement < numberOfElements; ++ iElement)
                {
                    receivingBuffer[iElement] = MAX(receivingBuffer[iElement], tmpBuffer[iElement]);
                }
            }
            else if (op == MPI_MIN)
            {
                for (int iElement = 0; iElement < numberOfElements; ++ iElement)
                {
                    receivingBuffer[iElement] = MIN(receivingBuffer[iElement], tmpBuffer[iElement]);
                }
            }
            else if (op == MPI_SUM)
            {
                for (int iElement = 0; iElement < numberOfElements; ++ iElement)
                {
                    receivingBuffer[iElement] += tmpBuffer[iElement];
                }
            }

            delete [] tmpBuffer;
        }
    }

    //! Bcast the final results to each processors.

    DataContainer *cdata;

    cdata = new DataContainer;

    if (myIDTmp == serverTmp)
    {
        cdata->Write(receivingBuffer, numberOfElements * sizeof(int));
    }

    PH_BcastByServerSepMode(cdata);

    cdata->MoveToBegin();

    cdata->Read(receivingBuffer, numberOfElements * sizeof(int));

    delete cdata;        cdata = nullptr;

    //ostringstream oss1;
    //oss1 << "After all reduce, sendingBuffer = " << receivingBuffer[0] << endl;
    //WriteLogFile(oss1.str());

#endif
}

void PH_AllreduceSepMode(float *sendingBuffer, float *receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op)
{
    using namespace PHMPI;

#ifdef PH_PARALLEL
    int numberOfGridProcessorTmp = GetNumberOfGridProcessor();

    //ostringstream oss;
    //oss << "sendingBuffer = " << sendingBuffer[0] << endl;
    //WriteLogFile(oss.str());

    if (numberOfGridProcessorTmp == 0)
    {
        MPI_Allreduce(sendingBuffer, receivingBuffer, numberOfElements, mpiDataType, op, MPI_COMM_WORLD);
        return;
    }

    int serverTmp = PHMPI::GetServerSepMode();

    int numberOfOriginalProcessor = GetNumberOfProcessor();

    int startProcessor = 0;
    int endProcessor = numberOfOriginalProcessor;

    if (CurrentProcessorIsGridProcessor())
    {
        startProcessor = numberOfOriginalProcessor;
        endProcessor = numberOfGridProcessorTmp + numberOfOriginalProcessor;
    }

    int myIDTmp = GetCurrentProcessorID();

    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        receivingBuffer[iElement] = sendingBuffer[iElement];
    }

    for (int iProc = startProcessor; iProc < endProcessor; ++ iProc)
    {
        if (iProc == serverTmp)
        {
            continue;
        }

        int tag = iProc;

        if (myIDTmp == iProc)
        {
            PH_SmartSend(sendingBuffer, numberOfElements, serverTmp, tag);
        }

        if (myIDTmp == serverTmp)
        {
            float *tmpBuffer = new float [numberOfElements];

            PH_SmartRecv(tmpBuffer, numberOfElements, iProc, tag);

            if (op == MPI_MAX)
            {
                for (int iElement = 0; iElement < numberOfElements; ++ iElement)
                {
                    receivingBuffer[iElement] = MAX(receivingBuffer[iElement], tmpBuffer[iElement]);
                }
            }
            else if (op == MPI_MIN)
            {
                for (int iElement = 0; iElement < numberOfElements; ++ iElement)
                {
                    receivingBuffer[iElement] = MIN(receivingBuffer[iElement], tmpBuffer[iElement]);
                }
            }
            else if (op == MPI_SUM)
            {
                for (int iElement = 0; iElement < numberOfElements; ++ iElement)
                {
                    receivingBuffer[iElement] += tmpBuffer[iElement];
                }
            }

            delete [] tmpBuffer;
        }
    }

    //! Bcast the final results to each processors.

    DataContainer *cdata;

    cdata = new DataContainer;

    if (myIDTmp == serverTmp)
    {
        cdata->Write(receivingBuffer, numberOfElements * sizeof(int));
    }

    PH_BcastByServerSepMode(cdata);

    cdata->MoveToBegin();

    cdata->Read(receivingBuffer, numberOfElements * sizeof(int));

    delete cdata;    cdata = nullptr;

    //ostringstream oss1;
    //oss1 << "After all reduce, sendingBuffer = " << receivingBuffer[0] << endl;
    //WriteLogFile(oss1.str());

#endif
}

void PH_AllreduceSepMode(double *sendingBuffer, double *receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op)
{
    using namespace PHMPI;

#ifdef PH_PARALLEL
    int numberOfGridProcessorTmp = GetNumberOfGridProcessor();

    //ostringstream oss;
    //oss << "sendingBuffer = " << sendingBuffer[0] << endl;
    //WriteLogFile(oss.str());

    if (numberOfGridProcessorTmp == 0)
    {
        MPI_Allreduce(sendingBuffer, receivingBuffer, numberOfElements, mpiDataType, op, MPI_COMM_WORLD);
        return;
    }

    int serverTmp = PHMPI::GetServerSepMode();

    int numberOfOriginalProcessor = GetNumberOfProcessor();

    int startProcessor = 0;
    int endProcessor = numberOfOriginalProcessor;

    if (CurrentProcessorIsGridProcessor())
    {
        startProcessor = numberOfOriginalProcessor;
        endProcessor = numberOfGridProcessorTmp + numberOfOriginalProcessor;
    }

    int myIDTmp = GetCurrentProcessorID();

    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        receivingBuffer[iElement] = sendingBuffer[iElement];
    }

    for (int iProc = startProcessor; iProc < endProcessor; ++ iProc)
    {
        if (iProc == serverTmp)
        {
            continue;
        }

        int tag = iProc;

        if (myIDTmp == iProc)
        {
            PH_SmartSend(sendingBuffer, numberOfElements, serverTmp, tag);
        }

        if (myIDTmp == serverTmp)
        {
            double *tmpBuffer = new double[numberOfElements];

            PH_SmartRecv(tmpBuffer, numberOfElements, iProc, tag);

            if (op == MPI_MAX)
            {
                for (int iElement = 0; iElement < numberOfElements; ++ iElement)
                {
                    receivingBuffer[iElement] = MAX(receivingBuffer[iElement], tmpBuffer[iElement]);
                }
            }
            else if (op == MPI_MIN)
            {
                for (int iElement = 0; iElement < numberOfElements; ++ iElement)
                {
                    receivingBuffer[iElement] = MIN(receivingBuffer[iElement], tmpBuffer[iElement]);
                }
            }
            else if (op == MPI_SUM)
            {
                for (int iElement = 0; iElement < numberOfElements; ++ iElement)
                {
                    receivingBuffer[iElement] += tmpBuffer[iElement];
                }
            }

            delete [] tmpBuffer;
        }
    }

    //! Bcast the final results to each processors.

    DataContainer *cdata;

    cdata = new DataContainer;

    if (myIDTmp == serverTmp)
    {
        cdata->Write(receivingBuffer, numberOfElements * sizeof(double));
    }

    PH_BcastByServerSepMode(cdata);

    cdata->MoveToBegin();

    cdata->Read(receivingBuffer, numberOfElements * sizeof(double));

    delete cdata;    cdata = nullptr;

    //ostringstream oss1;
    //oss1 << "After all reduce, sendingBuffer = " << receivingBuffer[0] << endl;
    //WriteLogFile(oss1.str());

#endif
}

int PH_Reduce(int *sendbuf, int *recvbuf, int count, PH_Op op)
{
    int errcode = 0;
    int root = PHMPI::GetServerProcessorID();

#ifdef PH_PARALLEL
    errcode = MPI_Reduce(sendbuf, recvbuf, count, MPI_INT, op, root, MPI_COMM_WORLD);
#endif

    return errcode;
}

int PH_Reduce(float *sendbuf, float *recvbuf, int count, PH_Op op)
{
    int errcode = 0;
    int root = PHMPI::GetServerProcessorID();

#ifdef PH_PARALLEL
    errcode = MPI_Reduce(sendbuf, recvbuf, count, MPI_FLOAT, op, root, MPI_COMM_WORLD);
#endif

    return errcode;
}

int PH_Reduce(double *sendbuf, double *recvbuf, int count, PH_Op op)
{
    int errcode = 0;
    int root = PHMPI::GetServerProcessorID();

    if (PHMPI::GetNumberOfProcessor() == 1)
    {
        for (int iData = 0; iData < count; ++ iData)
        {
            recvbuf[iData] = sendbuf[iData];
        }
        return 1;
    }

#ifdef PH_PARALLEL
    errcode = MPI_Reduce(sendbuf, recvbuf, count, MPI_DOUBLE, op, root, MPI_COMM_WORLD);
#endif

    return errcode;
}

//! flag:
//! MAX = 1;
//! MIN = 2;
void PH_CompareMaxMin(int &local_data, int flag)
{
    const int MAX = 1;
    const int MIN = 2;

    if (PHMPI::GetNumberOfProcessor() == 1)
    {
        return;
    }

#ifdef PH_PARALLEL
    int global_data = 0;
    if (flag == MAX)
    {
        MPI_Allreduce(&local_data, &global_data, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    else if (flag == MIN)
    {
        MPI_Allreduce(&local_data, &global_data, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    }
    local_data = global_data;
#endif
}

//! flag:
//! MAX = 1;
//! MIN = 2;
void PH_CompareMaxMin(RDouble &local_data, int flag)
{
    const int MAX = 1;
    const int MIN = 2;

#ifdef PH_PARALLEL
    RDouble global_data = 0.0;
    if (flag == MAX)
    {
        MPI_Allreduce(&local_data, &global_data, 1, PH_MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        //PH_AllreduceSepMode(&local_data, &global_data, 1, PH_MPI_DOUBLE, MPI_MAX);
    }
    else if (flag == MIN)
    {
        MPI_Allreduce(&local_data, &global_data, 1, PH_MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        //PH_AllreduceSepMode(&local_data, &global_data, 1, PH_MPI_DOUBLE, MPI_MIN);
    }
    local_data = global_data;
#endif
}

int PH_BarrierSepMode()
{
    int test = 1;
    int testAll = 1;

    PH_AllreduceSepMode(&test, &testAll, 1, MPI_INT, MPI_SUM);

    return testAll;
}

}