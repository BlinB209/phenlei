#include <iostream>
#include "PHMpi.h"
#include "Precision.h"
#include "Constants.h"
#include "AleManager.h"
#include "AleForceManager.h"
#include "Solver.h"
#include "Geometry.h"
#include "SixDofManager.h"
#include "DeformingSolverManager.h"
#include "MovingMesh.h"
#include "Zone.h"
#include "GridType.h"
#include "Region.h"
#include "OversetInformation.h"
#include "GlobalDataBase.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "IO_FileReader.h"
#include "TK_Exit.h"
#include "TK_Log.h"
#include "Task_ServerUpdateInterface.h"
#include "Geo_UnstructGrid.h"
#include "Glb_Dimension.h"
#include "Pre_WalldistCompute.h"

namespace PHSPACE
{

AleManager * aleManager = NULL;

void CreateAleManager()
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (!isUnsteady)
    {
        return;
    }

    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    if (!isAle)
    {
        return;
    }

    aleManager = new AleManager;

    aleManager->Initialize();
}

AleManager * GetAleManager()
{
    return aleManager;
}

void FreeAleManager()
{
    int IsUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    if (IsUnsteady  && isAle)
    {
        if (aleManager)
        {
            delete aleManager;
            aleManager = NULL;
        }
    }
}

void SolveSingleInnerIterationStepForAleSolvers()
{
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    if (!isAle)
    {
        return;
    }

    //! it does not need to solve Ale for quasi steady questions.
    int aleStartStrategy = GlobalDataBase::GetIntParaFromDB("aleStartStrategy");
    if (aleStartStrategy == -1)
    {
        return;
    }

    //! backup the old meshes information (n time) first in the inner iterations.
    int innStep = GlobalDataBase::GetIntParaFromDB("innstep");
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    if (innStep == 1)
    {
        BackupAllZonesOldGrid();

        RDouble aleSimulationTime = GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");
        RDouble aleTimeStep = GetAleTimeStep();

        aleSimulationTime += aleTimeStep;
        GlobalDataBase::UpdateData("aleSimulationTime", &aleSimulationTime, PHDOUBLE, 1);
    }

    //! if it is the first step of the inner step or implicit methods, update the mesh information.
    if (innStep == 1 || IsImplicitAleMethod())
    {
        ComputeAerodynamicForceForALE();

        aleManager->Run();

        //ComputeMetrics();
        aleManager->ComputeMetrics();

        aleManager->CommunicateCellCenterData();

        if (viscousType > ALGEBRAIC && numberOfMovingBodies > 1)
        {
            ComputeWalldist();
        }

        //! computing the faceNormalVelocityContainer. we do not compute the movement velocity for quasi steady questions.
        int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
        if (methodOfDualTime != 0)
        {
            aleManager->ComputeGridFaceVelocity();
        }
    }

    if (innStep == 1)
    {
        //! if it is the overset mesh, should dig the hole.
        int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
        if (isOverset)
        {
            //! clear up the store memory before assembling.
            DeAllocateOversetStorage();

            ResetAllZonesOversetGrid();

            InitializeOversetInterfaceTopology();

            //! the memory of the variables need to re open, when the overset mesh changes.
            AllocateOversetStorage();
        }
    }
}

void PostprocessAfterInnerIterationStepForAleSolvers()
{
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    if (!isAle)
    {
        return;
    }

    aleManager->UpdateAleForce();

    aleManager->Post();

    aleManager->DumpRestart();
}

AleManager::AleManager()
{
    aleForceManager      = NULL;
    sixDofSolverManager  = NULL;
    deformingManager     = NULL;
    movingMeshManager    = NULL;
    currentAleFile       = NULL;
}

AleManager::~AleManager()
{
    if (aleForceManager    ) delete aleForceManager;
    if (sixDofSolverManager) delete sixDofSolverManager;
    if (deformingManager   ) delete deformingManager;
    if (movingMeshManager  ) delete movingMeshManager;
}

void AleManager::Initialize()
{
    SychonizeAleTime();

    Allocate();

    int aleStartStrategy = GlobalDataBase::GetIntParaFromDB("aleStartStrategy");
    if (aleStartStrategy == 1)
    {
        Restart();
    }

    aleStartStrategy = GlobalDataBase::GetIntParaFromDB("aleStartStrategy");
    if (aleStartStrategy == -1)
    {
        aleManager->Run();

        aleManager->ComputeMetrics();

        aleManager->CommunicateCellCenterData();

        int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
        if (isOverset)
        {
            ResetAllZonesOversetGrid();

            InitializeOversetInterfaceTopology();
            //! the memory of the variables need to re open, when the overset mesh changes.
            AllocateOversetStorage();
        }
    }

    Dump();
}

void AleManager::ComputeMetrics()
{
    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();
    ActionKey *actkey = new ActionKey();
    actkey->action = COMPUTEMETRICS;

    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        int send_proc = GetZoneProcessorIDSepMode(iZone);
        int recv_proc = GetServerSepMode();

        int tag = GetSendRecvTag(GRID_BASED, 0, iZone);
        int myid = PHMPI::GetCurrentProcessorID();
        if (myid == send_proc)
        {
            //! If the ProcessorID of i-th zone equal to the running ProcessorID, compute grid metrics.
            Grid *grid = GetGrid(iZone, 0);
            grid->ComputeMetrics(actkey);
            if (grid->Type() == PHSPACE::STRUCTGRID)
            {
                UnstructGrid *gridUnstr = UnstructGridCast(grid);
                if (gridUnstr)
                {
                    gridUnstr->ComputeMetrics(actkey);
                }
            }
        }

        //! collect all info out by ComputeMetrics function to cdata.
        PH_Trade(actkey, send_proc, recv_proc, tag);

        if (myid == recv_proc)
        {
            //! If ProcessorID ProcessorID is serve ProcessorID, output the info in cdata.
            WriteScreen(actkey);
        }
    }
    delete actkey;    actkey = nullptr;
}

void AleManager::CommunicateCellCenterData()
{
    using namespace PHMPI;
    int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");
    for (int level = 0; level < nMGLevel; ++level)
    {
        ActionKey *actkey = new ActionKey();
        actkey->action = COMMCELLCENTERDATA;
        actkey->kind = GRID_BASED;
        actkey->level = level;

        if (IsNeedNonBlockingCommunication(actkey))
        {
            MainTaskNonBlocking(actkey);
        }
        else
        {
            MainTaskBlocking(actkey);
        }

        delete actkey;    actkey = nullptr;
    }
}

void AleManager::ComputeGridFaceVelocity()
{
    int numberOfZones = PHMPI::GetNumberofGlobalZones();

    int gridLevelIndex = 0;
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        Grid *grid = PHSPACE::GetGrid(iZone, gridLevelIndex);
        if (!grid) continue;

        grid->ComputeGridFaceVelocity();
    }
}

void AleManager::Allocate()
{
    InitializeAerodynamicForce();

    aleForceManager     = new AleForceManager       ;
    deformingManager    = new DeformingSolverManager;
    sixDofSolverManager = new SixDofSolverManager   ;
    movingMeshManager   = new MovingMeshManager     ;

    RDouble aleSimulationTime = 0.0;
    GlobalDataBase::UpdateData("aleSimulationTime", &aleSimulationTime, PHDOUBLE, 1);

    aleForceManager    ->Initialize();
    deformingManager   ->Initialize();
    sixDofSolverManager->Initialize();
    movingMeshManager  ->Initialize();
}

void AleManager::Restart()
{
    fstream aleFile;
    string aleFileName = "./results/aleRestartFile.dat";

    ifstream infile(aleFileName.c_str(), ios::in);
    if (!infile)
    {
        int aleStartStrategy = 0;
        GlobalDataBase::UpdateData("aleStartStrategy", &aleStartStrategy, PHINT, 1);
        return;
    }

    ios_base::openmode openMode = ios_base::in;
    PHSPACE::OpenFile(aleFile, aleFileName, openMode);
    currentAleFile = & aleFile;
    SetIsReadAleFile(true);

    RDouble aleSimulationTime = 0.0;
    //PHRead(aleFile, aleSimulationTime);
    aleFile >> aleSimulationTime;
    GlobalDataBase::UpdateData("aleSimulationTime", &aleSimulationTime, PHDOUBLE, 1);

    //! the mesh  need to transform to n time from 0 time, when restart the computing.
    aleForceManager    ->Restart();
    deformingManager   ->Restart();
    sixDofSolverManager->Restart();
    movingMeshManager  ->Run();
    sixDofSolverManager->ResetMassCenter();
    movingMeshManager  ->Restart();

    //! the mesh moves, it needs to recompute the mesh information.
    aleManager->ComputeMetrics();

    aleManager->CommunicateCellCenterData();

    PHSPACE::CloseFile(aleFile);
}

void AleManager::InitializeDimensionalReferenceTime()
{
    RDouble referenceVelocity = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("referenceVelocity");
    RDouble referenceLength   = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("referenceLength");

    RDouble dimensionalReferenceTime = referenceLength / referenceVelocity;
    PHSPACE::GlobalDataBase::UpdateData("dimensionalReferenceTime", &dimensionalReferenceTime, PHDOUBLE, 1);
}

void AleManager::SychonizeAleTime()
{
    RDouble aleSimulationTime = 0.0;
    PHSPACE::GlobalDataBase::UpdateData("aleSimulationTime", &aleSimulationTime, PHDOUBLE, 1);

    int aleStartStrategy = GlobalDataBase::GetIntParaFromDB("aleStartStrategy");
    if (aleStartStrategy == 0)
    {
        int outIterStep = 0;
        PHSPACE::GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);
        int innerIterStep = 0;
        PHSPACE::GlobalDataBase::UpdateData("innstep", &innerIterStep, PHINT, 1);
    }

    InitializeDimensionalReferenceTime();

    RDouble dimensionalReferenceTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalReferenceTime");

    RDouble physicalTime = aleSimulationTime / dimensionalReferenceTime;

    GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);
}

void AleManager::UpdateAleForce()
{
    //! update the aerodynamic force after the iterations to prevent get the wrong energy consumption when restarting.
    aleForceManager->Run();
}

void AleManager::Run()
{
    deformingManager    ->Run();
    sixDofSolverManager ->Run();
    movingMeshManager   ->Run();
}

void AleManager::Post()
{
    aleForceManager     ->Post();
    deformingManager    ->Post();
    sixDofSolverManager ->Post();
    movingMeshManager   ->Post();
}

void AleManager::Dump()
{
    deformingManager      ->Dump();
    sixDofSolverManager   ->Dump();
}

void AleManager::DumpRestart()
{
    using namespace PHMPI;
    int currentProcessor = GetCurrentProcessorID();
    int serverProcessor = GetServerProcessorID();
    if (currentProcessor != serverProcessor)
    {
        return;
    }

    fstream aleFile;

    int outerStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int numberOfFieldSaveSteps = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");
    int maxSimuStep = GlobalDataBase::GetIntParaFromDB("maxSimuStep");

    if (outerStep % numberOfFieldSaveSteps == 0 || outerStep == maxSimuStep)
    {
        string aleFileName = "./results/aleRestartFile.dat";
        ios_base::openmode openMode = ios_base::out;
        PHSPACE::OpenFile(aleFile, aleFileName, openMode);

        RDouble aleTime = GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");

        aleFile << aleTime << "\n";
        //PHWrite(aleFile, aleTime);

        aleForceManager->DumpRestartAleForceOfBodyContainer(aleFile);
        deformingManager->DumpRestart();
        sixDofSolverManager->DumpRestart(aleFile);

        PHSPACE::CloseFile(aleFile);
    }
}

void BackupAllZonesOldGrid()
{
    int numberOfZones = PHMPI::GetNumberofGlobalZones();

    int gridLevelIndex = 0;
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Grid *grid = GetGrid(iZone, gridLevelIndex);
        if (! grid) continue;
        grid->BackUpOldGrid();
    }
}

void DeAllocateOversetStorage()
{
    int nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");
    for (int iSolver = 0; iSolver < nSolver; ++iSolver)
    {
        int numberOfZones = PHMPI::GetNumberofGlobalZones();
        int gridLevelIndex = 0;
        for (int iZone = 0; iZone < numberOfZones; ++iZone)
        {
            Grid *grid = GetGrid(iZone, gridLevelIndex);
            if (!grid) continue;
            PHSolver *solver = GlobalSolvers::GetSolver(iZone, iSolver);
            solver->DeAllocateOversetInterfaceVar(grid);
        }
    }
}

void AllocateOversetStorage()
{
    int nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");
    for (int iSolver = 0; iSolver < nSolver; ++iSolver)
    {
        int numberOfZones = PHMPI::GetNumberofGlobalZones();
        int gridLevelIndex = 0;
        for (int iZone = 0; iZone < numberOfZones; ++iZone)
        {
            Grid *grid = GetGrid(iZone, gridLevelIndex);
            if (!grid) continue;
            PHSolver *solver = GlobalSolvers::GetSolver(iZone, iSolver);
            solver->RegisterOversetField();
        }
    }
}

bool isReadAleFile = false;

bool IsReadAleFile()
{
    return isReadAleFile;
}

void SetIsReadAleFile(bool isReadAleFileIn)
{
    PHSPACE::isReadAleFile = isReadAleFileIn;
}

fstream & GetCurrentAleFile()
{
    return aleManager->GetCurrentAleFile();
}

bool IsImplicitAleMethod()
{
    int keyMethod = PHSPACE::GlobalDataBase::GetIntParaFromDB("methodForKineticEquation");

    if (keyMethod == KINETIC_SCHEME::ADMAS_BASHFORTH_FIRST_ORDER)   return false;
    if (keyMethod == KINETIC_SCHEME::ADMAS_BASHFORTH_SECOND_ORDER)  return false;

    return true;
}

bool IsExplicitAleMethod()
{
    return ! IsImplicitAleMethod();
}

int GetIntegerParameterFromDataBase(int bodyIndex, const string& parameterName)
{
    ostringstream oss;
    oss << parameterName << "_" << bodyIndex;
    return PHSPACE::GlobalDataBase::GetIntParaFromDB(oss.str());
}

int GetIntegerParameterFromDataBaseIfExist(int bodyIndex, const string& parameterName, int defaultValue)
{
    int value = defaultValue;
    ostringstream oss;
    oss << parameterName << "_" << bodyIndex;
    if (GlobalDataBase::IsExist(oss.str(), PHINT, 1))
    {
        value = PHSPACE::GlobalDataBase::GetIntParaFromDB(oss.str());
    }
    return value;
}

RDouble GetRDoubleParameterFromDataBase(int bodyIndex, const string& parameterName)
{
    ostringstream oss;
    oss << parameterName << "_" << bodyIndex;
    return PHSPACE::GlobalDataBase::GetDoubleParaFromDB(oss.str());
}

RDouble GetRDoubleParameterFromDataBaseIfExist(int bodyIndex, const string &parameterName, RDouble defaultValue)
{
    RDouble value = defaultValue;
    ostringstream oss;
    oss << parameterName << "_" << bodyIndex;
    if (GlobalDataBase::IsExist(oss.str(), PHDOUBLE, 1))
    {
        value = PHSPACE::GlobalDataBase::GetDoubleParaFromDB(oss.str());
    }
    return value;
}

string GetStringParameterFromDataBase(int bodyIndex, const string &parameterName)
{
    ostringstream oss;
    oss << parameterName << "_" << bodyIndex;
    return PHSPACE::GlobalDataBase::GetStrParaFromDB(oss.str());
}

int * GetIntegerArrayFromDataBase(const string &nameOfIntegerArray, int numberOfElements)
{
    int *integerArray = new int[numberOfElements];

    GlobalDataBase::GetData(nameOfIntegerArray, integerArray, PHINT, numberOfElements);

    return integerArray;
}

void GetRDoubleVectorFromDataBase(vector< RDouble >& parameter, int bodyIndex, const string& parameterName, int numberOfElements)
{
    ostringstream oss;
    oss << parameterName << "_" << bodyIndex;
    PHSPACE::GlobalDataBase::GetRDoubleVectorFromDB(parameter, oss.str(), numberOfElements);
}

RDouble GetAleTimeStep()
{
    //!dimensional real physical time step.
    RDouble aleTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStepDimensional");

    return aleTimeStep;
}

}