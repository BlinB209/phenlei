#include "CFDSolver.h"
#include "Glb_Dimension.h"
#include "IO_FileName.h"
#include "TK_Parse.h"
#include "Constants.h"
#include "Residual.h"
#include "Force.h"
#include "AleForceManager.h"
#include "Post_WriteVisualFile.h"
#include "Task_ServerUpdateInterface.h"
#include "Task_ServerUpdateInterpoint.h"
#include "OversetInformation.h"
#include "Param_CFDSolver.h"
#include "FieldProxy.h"
#include "PHIO.h"
#include "TK_Exit.h"
#include "TK_Log.h"
#include <cstdio>
#include <stdio.h>
#include "TK_Time.h"

#ifdef USE_TecplotLib
#include "TECXXX.h"
#endif

#ifdef USE_CUDA
#include "TemporaryOperations.h"
#include "BasicDeviceVariables.h"
#include "GPUMPICommunication.h"
using namespace GPUMemory;
using namespace GPUControlVariables;
using namespace TemporaryOperations;
#endif

#include "Gas.h"

using namespace std;

namespace PHSPACE
{
bool cfd_debug = false;

CFDSolver::CFDSolver()
{
    postVisualization = 0;
    dqProxy           = 0;
    controlParameters = 0;
    postVisualWall = 0;
}

CFDSolver::~CFDSolver()
{
    if (postVisualization != NULL)
    {
        FreePointer(postVisualization);
    }
    if (postVisualWall != NULL)
    {
        FreePointer(postVisualWall);
    }
}

void CFDSolver::InitMemory()
{
    InitControlParameters();

    ReadParameter();

    AllocateGlobalVariables();
}

void CFDSolver::ReleaseMemory()
{
    DeAllocateGlobalVariables();

    FreeControlParameters();
}

void CFDSolver::InitFlow()
{
    int nstart = 0;

    //! This method would be wrong when turbulent computational from laminar flow field.
    //! But this situation breaks the law, so it is not considered.
    int readFlowFile = INITINFLOW;
    if (JudgeIfRestart())
    {
        readFlowFile = RESTARTFLIE;

        nstart = 1;
        GlobalDataBase::UpdateData("nstart", &nstart, PHINT, 1);
    }

    if (JudgeIfInterpolate())
    {
        readFlowFile = INTERPOLATEFLIE;

        nstart = 1;
        GlobalDataBase::UpdateData("nstart", &nstart, PHINT, 1);
    }
 

    if (!readFlowFile)
    {
        //int outerInterStep = 0;
        //GlobalDataBase::UpdateData("outnstep", &outerInterStep, PHINT, 1);

        nstart = 0;
        GlobalDataBase::UpdateData("nstart", &nstart, PHINT, 1);

        InitFlowAsRestart();
    }
    else if(readFlowFile == RESTARTFLIE)
    {
        //! To read data in restart file for continual simulation.
        ActionKey *actkeyReadRestartData = new ActionKey();
        FillActionKey(actkeyReadRestartData, READ_RESTART, 0);
        InitFlowAsReadingRestart(actkeyReadRestartData);
        FreePointer(actkeyReadRestartData);
    }
    else if(readFlowFile == INTERPOLATEFLIE)
    {
        //! To read data in restart file for continual simulation.
        ActionKey *actkeyReadInterpolateData = new ActionKey();
        FillActionKey(actkeyReadInterpolateData, READ_INTERPOLATE, 0);
        InitFlowByReadingInterpolate(actkeyReadInterpolateData);
        FreePointer(actkeyReadInterpolateData);
    }

    RDouble globalMinTimeStep = -LARGE;
    GlobalDataBase::UpdateData("globalMinTimeStep", &globalMinTimeStep, PHDOUBLE, 1);

    //! Compute gama, temperature, viscosity, and set boundary condition once.
    InitDependentVariables();
}

void CFDSolver::InitProtectedFlow()
{
    int nstart = 0;

    //! This method would be wrong when turbulent computational from laminar flow field.
    //! But this situation breaks the law, so it is not considered.
    bool readFlowFile = false;
    if (JudgeIfProtectedRestart())
    {
        readFlowFile = true;

        nstart = 1;
        GlobalDataBase::UpdateData("nstart", &nstart, PHINT, 1);
    }

    if (!readFlowFile)
    {
        //int outerInterStep = 0;
        //GlobalDataBase::UpdateData("outnstep", &outerInterStep, PHINT, 1);

        nstart = 0;
        GlobalDataBase::UpdateData("nstart", &nstart, PHINT, 1);

        InitFlowAsRestart();
    }
    else
    {
        //! To read data in restart file for continual simulation.
        ActionKey *actkeyReadRestartData = new ActionKey();
        FillActionKey(actkeyReadRestartData, READ_RESTART, 0);
        InitFlowAsReadingProtectedRestart(actkeyReadRestartData);
        FreePointer(actkeyReadRestartData);
    }

    RDouble globalMinTimeStep = -LARGE;
    GlobalDataBase::UpdateData("globalMinTimeStep", &globalMinTimeStep, PHDOUBLE, 1);

    //! Compute gama, temperature, viscosity, and set boundary condition once.
    InitDependentVariables();
}

void CFDSolver::InitCoarseGridsFlow ()
{

}

void CFDSolver::InitFlowAsReadingRestart()
{

}

bool CFDSolver::JudgeIfReadAverage()
{
    return false;
}

bool CFDSolver::JudgeIfRestart()
{
    return false;
}

bool CFDSolver::JudgeIfInterpolate()
{
    return false;
}

bool CFDSolver::JudgeIfProtectedRestart()
{
    return false;
}

void CFDSolver::InitControlParameters()
{

}

void CFDSolver::InitDependentVariables()
{

}

void CFDSolver::FreeControlParameters()
{
    if (controlParameters)
    {
        FreePointer(controlParameters);
    }
}

void CFDSolver::ReadParameter()
{

}

void CFDSolver::InitFlowAsRestart()
{

}

RDouble CFDSolver::UnsteadyConvergence(Grid *grid)
{ 
    return 0.0;
}

void CFDSolver::UpdateUnsteadyFlow(Grid *grid)
{

}

void CFDSolver::UpdateUnsteadyVariable(Grid *grid)
{

}

void CFDSolver::UpdateDt()
{

}

LIB_EXPORT Param_CFDSolver * CFDSolver::GetControlParameters() const
{
    return controlParameters;
}

void CFDSolver::ComputePreconditionCoefficient(Grid *grid)
{

}

void CFDSolver::ComputePreconditionCoefficientUnsteady(Grid *grid)
{

}

void CFDSolver::ConservativePreconMatrix(Grid *grid, int sign)
{

}
void CFDSolver::InviscidFlux(Grid *grid)
{

}

void CFDSolver::ViscousFlux(Grid *grid)
{

}

void CFDSolver::ComputeGradient(Grid *grid)
{

}

void CFDSolver::ComputeLimiter(Grid *grid)
{

}

void CFDSolver::Crossing(Grid *grid)
{

}

void CFDSolver::Blending(Grid *grid)
{

}

void CFDSolver::ComputeGradientCellCenter(Grid *grid)
{

}

void CFDSolver::ComputeGradientCellCenterForLES(Grid *grid)
{

}

void CFDSolver::SourceFlux(Grid *grid)
{

}

void CFDSolver::TimeStep(Grid *grid)
{

}

#ifdef USE_GMRESSOLVER
void CFDSolver::GMRESSolver(Grid *gridIn, FieldProxy *dqProxy)
{

}  
#endif

void CFDSolver::ComputeMinTimeStep(Grid *grid, RDouble &minDt, RDouble &maxDt)
{

}

void CFDSolver::ReduceMaxTimeStep(Grid *grid, RDouble globalMinDt)
{

}

void CFDSolver::PostSolve(Grid *grid, int stage, int level)
{

}

void CFDSolver::DumpResultFile(Grid *grid, int level)
{

}

void CFDSolver::CheckResult(Grid *grid, int level)
{

}

void CFDSolver::CheckResiduals(Grid *grid, int level)
{

}

void CFDSolver::DumpCFLNumber(Grid *grid, int level)
{

}

void CFDSolver::CheckSurfaceHeatingChange(Grid *grid, int level)
{

}

void CFDSolver::DumpLeakageInformation(Grid *grid)
{

}

void CFDSolver::ZeroResiduals(Grid *grid)
{

}

void CFDSolver::UpdateResiduals(Grid *grid) 
{

}

void CFDSolver::RestrictAllQ(Grid *fineGrid, Grid *coarseGrid)
{

}

void CFDSolver::RestrictDefect(Grid *fineGrid, Grid *coarseGrid)
{

}

void CFDSolver::PutCorrectionBack(Grid *grid, FieldProxy *qProxy)
{

}

void CFDSolver::CorrectFineGrid(Grid *fineGrid, Grid *coarseGrid)
{

}

void CFDSolver::InterpolatFineGrid(Grid *fineGrid, Grid *coarseGrid)
{

}

void CFDSolver::StoreRhsByResidual(Grid *grid, FieldProxy *rhsProxy) 
{

}

void CFDSolver::InitResidual(Grid *grid, FieldProxy *rhsProxy)
{

}

void CFDSolver::RecoverResidual(Grid *grid, FieldProxy *rhsProxy)
{

}

void CFDSolver::LoadQ(Grid *grid, FieldProxy *qProxy)
{

}

FieldProxy *CFDSolver::CreateFieldProxy(Grid *grid) 
{
    return 0;
}


void CFDSolver::Solve()
{

}

void CFDSolver::DumpRestart(ActionKey *actkey)
{

}

void CFDSolver::DumpRestartH5(ActionKey *actkey)
{

}

void CFDSolver::ReadRestart(ActionKey *actkey)
{

}

void CFDSolver::ReadRestartH5(ActionKey *actkey)
{

}

void CFDSolver::InterpolateFromRestartH5(ActionKey* actkey)
{

}

void CFDSolver::AllocateGlobalVariables()
{
    Grid *grid = GetGrid();
    while (grid)
    {
        AllocateGlobalVar(grid);
        grid = grid->GetCoarseGrid();
    }

    RegisterCFDSolverInterfaceField();
    RegisterCFDSolverInterpointField();

    Param_CFDSolver *parameters = GetControlParameters();
    int isOverset = parameters->GetIsOverLapping();

    if (isOverset)
    {
        RegisterOversetField();
    }
}

void CFDSolver::DeAllocateGlobalVariables()
{
    Grid *grid = GetGrid(0);
    int isOverset = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    while (grid)
    {
        DeAllocateGlobalVar(grid);
        if (isOverset)
        {
            DeAllocateOversetInterfaceVar(grid);
        }        
        grid = grid->GetCoarseGrid();
    }

    ReleaseInterfaceField();
}

void CFDSolver::AllocateGlobalVar(Grid *grid)
{
    //! When multi-grid, only need allocate once for postVisualization.
    if (postVisualization == 0)
    {
        int nVisualVariables = GlobalDataBase::GetIntParaFromDB("nVisualVariables");
        int visualVariables[100] = {0};
        GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);
        postVisualization = new Post_Visual(nVisualVariables, visualVariables);
        SetPostVisualization(postVisualization);
    }

    if (postVisualWall == 0)
    {
        int nVisualWallVariables = GlobalDataBase::GetIntParaFromDB("nVisualWallVariables");
        int visualWallVariables[30] = {0};
        GlobalDataBase::GetData("visualWallVariables", visualWallVariables, PHINT, nVisualWallVariables);

        Post_VisualWall *postVis = new Post_VisualWall(nVisualWallVariables, visualWallVariables);
        this->SetPostVisualizationWall(postVis);
    }
}

void CFDSolver::DeAllocateGlobalVar(Grid *grid)
{
    if (postVisualization != 0)
    {
        FreePointer(postVisualization);
    }

    if (postVisualWall != 0)
    {
        FreePointer(postVisualWall);
    }
}

void CFDSolver::RegisterCFDSolverInterfaceField()
{

}

void CFDSolver::RegisterCFDSolverInterpointField()
{

}

void CFDSolver::RegisterOversetField()
{

}

void CFDSolver::RegisterInterfaceField(const string &name, const int type, const int dimesion)
{
    int solverID = this->GetIndex();
    Grid *grid = GetGrid(0);
    while (grid)
    {
        if (grid->GetInterfaceInfo())
        {
            grid->RegisterInterfaceField(name, type, dimesion, solverID);
        }
        grid = grid->GetCoarseGrid();
    }
}

void CFDSolver::ReleaseInterfaceField()
{
    Grid *grid = GetGrid(0);
    while (grid)
    {
        if (grid->GetInterfaceInfo())
        {
            grid->ReleaseInterfaceField();
        }
        grid = grid->GetCoarseGrid();
    }
}

void CFDSolver::ReleaseInterfaceField(string varName)
{
    Grid *grid = GetGrid(0);
    while (grid)
    {
        if (grid->GetInterfaceInfo())
        {
            grid->ReleaseInterfaceField(varName);
        }
        grid = grid->GetCoarseGrid();
    }
}

void CFDSolver::RegisterInterpointField(const string & name, const int type, const int dimesion)
{
    int solverID = this->GetIndex();
    Grid *grid = GetGrid(0);
    while (grid)
    {
        if (grid->GetInterpointInfo())
        {
            grid->RegisterInterpointField(name, type, dimesion, solverID);
        }
        grid = grid->GetCoarseGrid();
    }
}

LIB_EXPORT void CFDSolver::RegisterOversetField(const string &name, const int type, const int dimesion)
{
    Grid *grid = GetGrid(0);
    while (grid)
    {
        UnstructGridCast(grid)->RegisterOversetField(name, type, dimesion);
        grid = grid->GetCoarseGrid();
    }
}

void CFDSolver::DeAllocateOversetInterfaceVar(Grid *grid)
{
    UnstructGrid *gridUnstr = UnstructGridCast(grid);
    if (gridUnstr == nullptr)
    {
        return;
    }

    OversetInformationProxy *oversetInformationProxy = gridUnstr->GetOversetInformationProxy();
    if (oversetInformationProxy)
    {
        DeAllocateOversetInterfaceVar(oversetInformationProxy);
    }
    }

void CFDSolver::DeAllocateOversetInterfaceVar(OversetInformationProxy *oversetInformationProxy)
{
    DeAllocateOversetInterfaceVar(oversetInformationProxy->GetOversetDataProxyForReceive()->GetDataStorage());
    DeAllocateOversetInterfaceVar(oversetInformationProxy->GetOversetDataProxyForSend()->GetDataStorage());
}

void CFDSolver::DeAllocateOversetInterfaceVar(Data_ParamFieldSuite *dataStore)
{
}

void CFDSolver::GetInterfaceData(ActionKey *actkey)
{
    int level = actkey->level;
    Grid *grid = GetGrid(level);
    InterfaceInfo *iterfaceInfo = grid->GetInterfaceInfo();
    if (iterfaceInfo == 0)
    {
        return;
    }

    InterfaceDataProxy *interfaceDataProxy = new InterfaceDataProxy();

    for (int indexGhostCellLayer = GetNumberOfGhostCellLayers() - 1; indexGhostCellLayer >= 0; -- indexGhostCellLayer)
    {
        PrepareInterfaceData(grid, iterfaceInfo->GetSendDataStorage(indexGhostCellLayer), interfaceDataProxy);
    }

    GetInterfaceData(grid->GetFinestGrid()->GetInterfaceInfo(),actkey,interfaceDataProxy);

    FreePointer(interfaceDataProxy);
}

void CFDSolver::GetInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    Grid *grid = GetGrid(level);
    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation == 0)
    {
        return;
    }
    InterpointDataProxy *interpointDataProxy = new InterpointDataProxy();

    for (int i = GetNumberOfGhostCellLayers() - 1; i >= 0; -- i)
    {
        PrepareInterpointData(grid, interpointInformation->GetSendDataStorage(i), interpointDataProxy);
    }

    //PHSPACE::GetInterpointData(interfaceInfo,actkey,interpointDataProxy);
    GetInterpointData(grid->GetFinestGrid()->GetInterpointInfo(),actkey,interpointDataProxy);

    FreePointer(interpointDataProxy);
}

void CFDSolver::TranslateInterfaceData(ActionKey *actkey)
{
    int level = actkey->level;
    Grid *grid = GetGrid(level);
    InterfaceInfo *iterfaceInfo = grid->GetInterfaceInfo();
    if (iterfaceInfo == 0)
    {
        return;
    }

    InterfaceDataProxy *interfaceDataProxy = new InterfaceDataProxy();

    for (int indexGhostCellLayer = GetNumberOfGhostCellLayers() - 1; indexGhostCellLayer >= 0; -- indexGhostCellLayer)
    {
        PrepareInterfaceData(grid, iterfaceInfo->GetRecvDataStorage(indexGhostCellLayer), interfaceDataProxy);
    }

    TranslateInterfaceData(grid->GetFinestGrid()->GetInterfaceInfo(),actkey,interfaceDataProxy);
    FreePointer(interfaceDataProxy);
}

void CFDSolver::TranslateInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    Grid *grid = GetGrid(level);
    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation == 0)
    {
        return;
    }
    InterpointDataProxy *interpointDataProxy = new InterpointDataProxy();

    for (int i = GetNumberOfGhostCellLayers() - 1; i >= 0; -- i)
    {
        PrepareInterpointData(grid, interpointInformation->GetReceiveDataStorage(i), interpointDataProxy);
    }

    TranslateInterpointData(grid->GetFinestGrid()->GetInterpointInfo(),actkey,interpointDataProxy);
    FreePointer(interpointDataProxy);
}

streamsize CFDSolver::TranslateInterfaceDataLength(ActionKey *actkey)
{
    int level = actkey->level;
    Grid *grid = GetGrid(level);
    InterfaceInfo *iterfaceInfo = grid->GetInterfaceInfo();
    if (iterfaceInfo == 0)
    {
        return 0;
    }
    InterfaceDataProxy *interfaceDataProxy = new InterfaceDataProxy();

    for (int indexGhostCellLayer = GetNumberOfGhostCellLayers() - 1; indexGhostCellLayer >= 0; -- indexGhostCellLayer)
    {
        PrepareInterfaceData(grid, iterfaceInfo->GetRecvDataStorage(indexGhostCellLayer), interfaceDataProxy);
    }

    streamsize nlen = TranslateInterfaceDataLength(grid->GetFinestGrid()->GetInterfaceInfo(),actkey,interfaceDataProxy);

    FreePointer(interfaceDataProxy);

    return nlen;
}

streamsize CFDSolver::TranslateInterpointDataLength(ActionKey *actkey)
{
    int level = actkey->level;
    Grid *grid = GetGrid(level);
    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation == 0)
    {
        return 0;
    }
    InterpointDataProxy *interpointDataProxy = new InterpointDataProxy();

    for (int indexGhostCellLayer = GetNumberOfGhostCellLayers() - 1; indexGhostCellLayer >= 0; -- indexGhostCellLayer)
    {
        PrepareInterpointData(grid, interpointInformation->GetReceiveDataStorage(indexGhostCellLayer), interpointDataProxy);
    }

    streamsize nlen = TranslateInterpointDataLength(grid->GetFinestGrid()->GetInterpointInfo(),actkey,interpointDataProxy);

    FreePointer(interpointDataProxy);

    return nlen;
}

void CFDSolver::GetOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    Grid *grid = GetGrid(level);

    OversetInformationProxy *oversetInformationProxy = UnstructGridCast(grid)->GetOversetInformationProxy();

    InterfaceDataProxy *interfaceDataProxy = new InterfaceDataProxy();

    OversetDataProxy *oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForSend();

    Data_ParamFieldSuite *dataStorageForSend = oversetDataProxy->GetDataStorage();

    PrepareOversetInterfaceData(dataStorageForSend, interfaceDataProxy);

    GetOversetInterfaceData(oversetInformationProxy, actkey, interfaceDataProxy);

    FreePointer(interfaceDataProxy);
}

void CFDSolver::TranslateOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    Grid *grid = GetGrid(level);

    OversetInformationProxy *oversetInformationProxy = UnstructGridCast(grid)->GetOversetInformationProxy();

    InterfaceDataProxy *interfaceDataProxy = new InterfaceDataProxy();

    OversetDataProxy *oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForReceive();

    Data_ParamFieldSuite *dataStorageForReveive = oversetDataProxy->GetDataStorage();

    PrepareOversetInterfaceData(dataStorageForReveive, interfaceDataProxy);

    TranslateOversetInterfaceData(oversetInformationProxy, actkey, interfaceDataProxy);

    FreePointer(interfaceDataProxy);
}

void CFDSolver::Boundary(Grid *grid)
{

}

void CFDSolver::CalculateBoundaryData()
{

}

void CFDSolver::ModifyResiduals(Grid *grid)
{

}

void CFDSolver::GetResidual(ActionKey *actkey)
{

}

void CFDSolver::GetResidual(ActionKey *actkey, RDouble &localMaxRes)
{

}

void CFDSolver::GetResidual(Grid *gridIn, RDouble &localMaxRes)
{

}

void CFDSolver::ObtainCFLNumber(Grid *gridIn, RDouble &globalMinCFL, RDouble &globalMaxCFL)
{

}

void CFDSolver::GetSurfaceHeatingChange(ActionKey *actkey, RDouble &localMaxHeatChange)
{

}

void CFDSolver::PrepareInterfaceData(Grid *grid, Data_ParamFieldSuite *dataStore, InterfaceDataProxy *interfaceDataProxy)
{
    InterfaceFields *interfaceFields = grid->GetInterfaceFields();
    if (!interfaceFields)
    {
        return;
    }

    vector <RDouble **> &vectorData = interfaceDataProxy->GetVectorData();
    vector <int> &vectorDimension = interfaceDataProxy->GetVectorDimension();
    vector<string> &vectorName = interfaceDataProxy->GetVectorName();

    for (int iData = 0; iData < interfaceFields->Size(); ++ iData)
    {
        int solverID = interfaceFields->GetSolverIndex(iData);
        if (solverID != this->GetIndex())
        {
            //! This data is not belong to this solver.
            continue;
        }

        int dataDimension = interfaceFields->GetDim(iData);
        string &dataName = interfaceFields->GetName(iData);
        RDouble **data = reinterpret_cast <RDouble **> (dataStore->GetDataPtr(dataName));

        vectorDimension.push_back(dataDimension);
        vectorData.push_back(data);
        vectorName.push_back(dataName);
    }
}

void CFDSolver::PrepareInterpointData(Grid *grid, Data_ParamFieldSuite *datastore, InterpointDataProxy *interpointDataProxy)
{
    InterpointFields *interpointFields = grid->GetInterpointFields();
    if (!interpointFields)
    {
        return;
    }
    vector<RDouble **> &vectorData = interpointDataProxy->GetVectorData();
    vector<int>        &vectorDimension = interpointDataProxy->GetVectorDimension();

    for (int iData = 0; iData < interpointFields->Size(); ++ iData)
    {
        int solverID = interpointFields->GetSolverIndex(iData);
        if (solverID != this->GetIndex())
        {
            //! This data is not belong to this solver.
            continue;
        }

        int dataDimension = interpointFields->GetDimension(iData);
        string &dataName = interpointFields->GetName(iData);
        RDouble **data = reinterpret_cast<RDouble **>(datastore->GetDataPtr(dataName));

        vectorDimension.push_back(dataDimension);
        vectorData.push_back(data);
    }
}

void CFDSolver::PrepareOversetInterfaceData(Data_ParamFieldSuite *dataStore, InterfaceDataProxy *interfaceDataProxy)
{

}

void CFDSolver::ActionReflect(ActionTag *acttag)
{

}

void CFDSolver::VisualizationAverageFlow(ActionKey *actkey)
{

}

void CFDSolver::CommunicationInterfaceData()
{
    Param_CFDSolver *parameters = GetControlParameters();
    int isOverset = parameters->GetIsOverLapping();

    //! Firstly the overset information should be exchanged.
    //! The structured overset grid should be taken into consideration further!!!
    if (isOverset)
    {
        UploadOversetData();
        UpdateOversetData();
        DownloadOversetData();
    }

    UploadInterfaceData();
    UpdateInterfaceData();
    DownloadInterfaceData();
}

void CFDSolver::CommunicationOversetSlipData()
{
    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");
    int isCommGradientforSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isCommGradientforSlip");
    if (isOversetSlip && isCommGradientforSlip)
    {
        UploadOversetData();
        UpdateOversetData();
        DownloadOversetData();
    }
}
void CFDSolver::CommunicationInterfaceData(int iLevel)
{
    Param_CFDSolver *parameters = GetControlParameters();
    int isOverLapping = parameters->GetIsOverLapping();

    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");

    //! Firstly the overset information should be exchanged.
    //! The structured overset grid should be taken into consideration further!!!
    if (isOverLapping && systemGridType != STRUCTGRID)
    {
        TK_Exit::ExceptionExit("Consider each level for over-set grid.");
        //UploadOversetData();
        //UpdateOversetData();
        //DownloadOversetData();
    }

    UploadInterfaceData(iLevel);
    UpdateInterfaceData(iLevel);
    DownloadInterfaceData(iLevel);
}

void CFDSolver::CommunicationInterpointData()
{
    UploadInterpointData();
    UpdateInterpointData();
    DownloadInterpointData();

    CommunicationInterpointWeight();
}

void UpdateInterfaceData()
{
    using namespace PHMPI;

    for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
    {
        ActionKey *actkey = new ActionKey();
        actkey->action = UPDATE_INTERFACE_DATA;
        actkey->kind = SOLVER_BASED;
        actkey->level = iLevel;
        actkey->solverID = 0;

        //Task_ServerUpdateInterface *task = new Task_ServerUpdateInterface();
        if (IsNeedNonBlockingCommunication(actkey))
        {
            MainTaskNonBlocking(actkey);
        }
        else
        {
            MainTaskBlocking(actkey);
        }
        //delete task;

        FreePointer(actkey);
    }
}

void CFDSolver::UploadInterfaceValue(ActionKey *actkey)
{

}

void CFDSolver::DownloadInterfaceValue(ActionKey *actkey)
{

}

void CFDSolver::UploadInterfaceData(ActionKey *actkey)
{

}

void CFDSolver::DownloadInterfaceData(ActionKey *actkey)
{

}

void CFDSolver::UploadInterpointData(ActionKey *actkey)
{

}

void CFDSolver::DownloadInterpointData(ActionKey *actkey)
{

}

void CFDSolver::CommunicationInterpointWeight(ActionKey *actkey)
{

}

void CFDSolver::DownloadInterpointWeight(ActionKey *actkey)
{

}

void CFDSolver::UploadOversetData(ActionKey *actkey)
{

}

void CFDSolver::DownloadOversetData(ActionKey *actkey)
{

}

void CFDSolver::FillField(Grid *grid, FieldProxy *fieldProxy, RDouble value)
{

}

void CFDSolver::FillField(Grid *grid, FieldProxy *field1Proxy, FieldProxy *field2Proxy)
{

}

FieldProxy * CFDSolver::GetFieldProxy(Grid *grid, const string &fieldName) 
{
    return 0;
}

void CFDSolver::LoadResiduals(Grid *grid, FieldProxy *rhsProxy)
{

}

void CFDSolver::RungeKuttaResidual(Grid *grid, FieldProxy *dqProxy, RDouble coef)
{

};

void CFDSolver::LHSMatrixs(Grid *grid)
{

}

void CFDSolver::LineLUSGSMatrixs(Grid *grid)
{

}

FieldProxy *CFDSolver::GetResidualProxy(Grid *grid) 
{
    return 0; 
}

void CFDSolver::UpdateFlowField(Grid *grid, FieldProxy *qProxy, FieldProxy *dqProxy)
{

}

void CFDSolver::CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn)
{

}

void CFDSolver::DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn)
{

}

void CFDSolver::CompressAnInterfaceVar(DataContainer *&dataContainer, RDouble **fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn)
{

}

void CFDSolver::DecompressAnInterfaceVar(DataContainer *&dataContainer, RDouble **fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn)
{

}

void CFDSolver::DumpInterfaceInfoBell() 
{

}

void CFDSolver::SolveOnGridSteady(Grid *grid)
{
    //! Zero the residuals of the finest grid.
    ZeroResiduals(grid);

    //! Solve for all grids starting from the finest grid.
    SolveOnGrid(grid);
}

void CFDSolver::SolveOnGridUnsteady(Grid *grid)
{
    RDouble tol_sub_iter = 0.01;
    GlobalDataBase::GetData("tol_sub_iter", &tol_sub_iter, PHDOUBLE, 1);

    int max_sub_iter = 1;
    GlobalDataBase::GetData("max_sub_iter", &max_sub_iter, PHINT, 1);

    int neqn = GetNumberOfEquations(); 

    FieldProxy *qnProxy = CreateFieldProxy(grid);    //! It represents the value of q at time step n. 
    FieldProxy *q1Proxy = CreateFieldProxy(grid);    //! It represents the q value of p-th iteration at time step n+1. 
    FieldProxy *q2Proxy = GetFieldProxy(grid, CastName("q"));    //! It represents the q value of p+1-th iteration at time step n+1.
    //! Notice: q_2 changes with iteration,it represents q(n+1)(p+1).

    LoadQ(grid, qnProxy);

    int iner_iter = 0;
    RDouble norm0, norm1, cvg;
    while (true)
    {
        iner_iter = iner_iter + 1;
        LoadQ(grid, q1Proxy);
        ZeroResiduals(grid);
        SolveOnGrid(grid);
        
        if (iner_iter == 1)
        {
            cvg = one;
        }
        else
        {
            grid->CVGNorm(qnProxy, q2Proxy, neqn, norm0);
            grid->CVGNorm(q1Proxy, q2Proxy, neqn, norm1);
            cvg = norm1 / (norm0 + SMALL);
        }

        cout << "iner_iter = " << iner_iter << " cvg = " << cvg << "\n";

        if (cvg <= tol_sub_iter)
        {
            break;
        }
        if (iner_iter >= max_sub_iter)
        {
            break;
        }
    }

    FreePointer(qnProxy);
    FreePointer(q1Proxy);
    FreePointer(q2Proxy);
}

//! Execute the multi-grid flow solver on grid level.
void CFDSolver::SolveOnGrid(Grid *grid)
{

}

void CFDSolver::ComputeResidual(Grid *grid, FieldProxy *rhsProxy)
{
    LoadResiduals(grid, rhsProxy);

    UpdateResiduals(grid);
}

void CFDSolver::RungeKutta(Grid *grid, FieldProxy *rhsProxy)
{
    //! Explicit multi-stage time marching method.
    int RKStage = GlobalDataBase::GetIntParaFromDB("RKStage");

    RDouble *lamda = new RDouble [RKStage];
    GlobalDataBase::GetData("lamda", lamda, PHDOUBLE, RKStage);

    //! Using get or create temporarily, since the UpdateFlowField function
    //! has not been consistent for NS/Turb solver.
    FieldProxy *qProxy = 0;
    if (this->GetIndex() == 0)
    {
        //! NS solver.
        qProxy = CreateFieldProxy(grid);
    }else
    {
        //! Turb solver.
        qProxy = GetFieldProxy(grid, CastName("q"));
    }

    FieldProxy *q0ProxyTemp = CreateFieldProxy(grid);

    //! The real residual.
    FieldProxy *resProxy = GetResidualProxy(grid);

    //! The right term of RK, named as RK-residual, which is lamda*dt*res/volume.
    FieldProxy *resRK = CreateFieldProxy(grid);

#ifdef USE_CUDA
    using namespace GPUMemory;
    using namespace GPUControlVariables;
    using namespace TemporaryOperations;
    SetGPUNSFillField(d_FillField_ns, "q_proxyToq_proxy_temp");
    SetGPUTurbFillField(d_FillField_turb, "q_proxyToq_proxy_temp");
#endif

    LoadQ(grid, qProxy);

    //! Back up the Q0.
    FillField(grid, q0ProxyTemp, qProxy); 

    for (int iStage = 0; iStage < RKStage; ++ iStage)
    {
        LoadResiduals(grid, rhsProxy);

        UpdateResiduals(grid);

    #ifdef USE_CUDA
        using namespace GPUMemory;
        using namespace GPUControlVariables;
        using namespace TemporaryOperations;
        SetGPUNSFillField(d_FillField_ns, "empty");
        SetGPUTurbFillField(d_FillField_turb, "empty");
    #endif 

        FillField(grid, resRK, resProxy);

        //! Compute the DQ: (Dt * Res / vol).
        RungeKuttaResidual(grid, resRK, lamda[iStage]);

    #ifdef USE_CUDA
        SetGPUNSFillField(d_FillField_ns, "q_proxy_tempToq_proxy");
        SetGPUTurbFillField(d_FillField_turb, "q_proxy_tempToq_proxy");
    #endif

        //! Give the back-up Q0 back to Q, in order to compute new Q.
        FillField(grid, qProxy, q0ProxyTemp);

        //! Q = Q0 - DQ.
        //! Warning: the targets are qProxy and grid-q in NS and turbulence solver, respectively.
        UpdateFlowField(grid, qProxy, resRK);
    }

    DelPointer(lamda);
    FreePointer(qProxy);
    FreePointer(q0ProxyTemp);
    FreePointer(resProxy);
    FreePointer(resRK);
}

void CFDSolver::SolveLUSGSForward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{

}

void CFDSolver::SolveLUSGSBackward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{

}

void CFDSolver::SolveBLUSGSForward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{

}

void CFDSolver::SolveBLUSGSBackward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{

}

void CFDSolver::SolveMatrixLUSGSForward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{

}

void CFDSolver::SolveMatrixLUSGSBackward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{

}

void CFDSolver::LineLUSGS(Grid *grid, int level)
{

}

void CFDSolver::SolveLineLUSGSForward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{

}

void CFDSolver::SolveLineLUSGSBackward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{

}

void CFDSolver::DetermineCFLNumber(Grid *grid, FieldProxy *dqProxy, FieldProxy *dqProxy1)
{

}

//! Put correction on the coarse grid.
void CFDSolver::PutCorrection(Grid *grid, FieldProxy *qProxy)
{
    int neqn = GetNumberOfEquations();

    FieldProxy *old_q_proxy = GetFieldProxy(grid, CastName("q"));

    PutCorrection(grid, qProxy, old_q_proxy, neqn);

    FreePointer(old_q_proxy);
}

void CFDSolver::Action(ActionKey *actkey)
{
    switch (actkey->action)
    {
        case UPDATE_INTERFACE_DATA:
            GetInterfaceData(actkey);
            break;
        case UPDATE_INTERPOINT_DATA:
            GetInterpointData(actkey);
            break;
        case UPDATE_OVERSET_DATA:
            GetOversetData(actkey);
            break;
        case VISUALIZATION_AVERAGE_FLOW:
            VisualizationAverageFlow(actkey);
            break;
        default:
            break;
    } 
}

void CFDSolver::TranslateAction(ActionKey *actkey)
{
    switch (actkey->action)
    {
        case UPDATE_INTERFACE_DATA:
            TranslateInterfaceData(actkey);
            break;
        case UPDATE_OVERSET_DATA:
            TranslateOversetData(actkey);
            break;
        case UPDATE_INTERPOINT_DATA:
            TranslateInterpointData(actkey);
            break;
        default:
            break;
    }
}

streamsize CFDSolver::TranslateActionLength(ActionKey *actkey)
{
    if (actkey->action == UPDATE_INTERFACE_DATA)
    {
        return TranslateInterfaceDataLength(actkey);
    }
    else if (actkey->action == UPDATE_INTERPOINT_DATA)
    {
        return TranslateInterpointDataLength(actkey);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("action", actkey->action);
        return 0;
    }
}

bool CFDSolver::IsNeedStatistics()
{
    Param_CFDSolver *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        if (ifStaticsFlowField > 0)
        {
            int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
            int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");

            if (outnstep >= startStatisticStep)
            {
                return true;
            }
        }
    }
    return false;
}

bool CFDSolver::IsNeedReynoldsStressStatistics()
{
    Param_CFDSolver *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        int ifStaticsReynoldsStress = parameters->GetIfStaticsReynoldsStress();
        if (ifStaticsReynoldsStress > 0)
        {
            int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
            int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");

            if (outnstep >= startStatisticStep)
            {
                return true;
            }
        }
    }
    return false;
}

void CFDSolver::GetInterfaceData(InterfaceInfo *iterfaceInfo, ActionKey *actkey, InterfaceDataProxy *interfaceDataProxy)
{
    if (iterfaceInfo == 0)
    {
        return;
    }

    vector <RDouble **> &vectorData = interfaceDataProxy->GetVectorData();
    vector <int> &vectorDimension = interfaceDataProxy->GetVectorDimension();
    vector<string> &vectorName = interfaceDataProxy->GetVectorName();

    int ineighbor = iterfaceInfo->FindIthNeighbor(actkey->ipos);

    int dim_sum = 0;

    for (std::size_t iDimension = 0; iDimension < vectorDimension.size(); ++ iDimension)
    {
        dim_sum += vectorDimension[iDimension];
    }

    int nIFaceOfNeighbor = iterfaceInfo->GetNIFaceOfNeighbor(ineighbor);
    int *faceIndexForSend = iterfaceInfo->GetFaceIndexForSend(ineighbor);

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();
    for (std::size_t iDim = 0; iDim < vectorDimension.size(); ++ iDim)
    {
        string dataName = vectorName[iDim];

        int iFace;
        int nm = vectorDimension[iDim];
        RDouble **q = vectorData[iDim];

        //! Bell 20130129 mod from Hexin.
        int nlen = nm * nIFaceOfNeighbor;

        RDouble *qtmp = new RDouble [nlen];

        int icount = 0;
        for (int jFaceOfNeighbor = 0; jFaceOfNeighbor < nIFaceOfNeighbor; ++ jFaceOfNeighbor)
        {
            iFace = faceIndexForSend[jFaceOfNeighbor];
            for (int m = 0; m < nm; ++ m)
            {
                qtmp[icount ++] = q[m][iFace];
            }
        }
        cdata->Write(qtmp, nlen * sizeof(RDouble));
        DelPointer(qtmp);
    }
}

void CFDSolver::GetInterpointData(InterpointInformation *interpointInformation, ActionKey *actkey, InterpointDataProxy *interpointDataProxy)
{
    if (interpointInformation == 0)
    {
        return;
    }
    vector<RDouble **>&vectorData = interpointDataProxy->GetVectorData();
    vector<int> vectorDimension = interpointDataProxy->GetVectorDimension();

    int iNeighbor = interpointInformation->FindIthNeighbor(actkey->ipos);

    int sumDimension = 0;

    for (std::size_t i = 0; i < vectorDimension.size(); ++ i)
    {
        sumDimension += vectorDimension[i];
    }

    int numberOfInterpointsForNeighbor = interpointInformation->GetNumberOfInterpointsForNeighbor(iNeighbor);
    int *pointIndexForSend = interpointInformation->GetPointIndexForSend(iNeighbor);

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();
    for (std::size_t idim = 0; idim < vectorDimension.size(); ++ idim)
    {
        int ipoint;
        int nm = vectorDimension[idim];
        RDouble **q = vectorData[idim];

        int nlen = nm * numberOfInterpointsForNeighbor;

        RDouble *qtmp = new RDouble [nlen];

        int icount = 0;
        for (int j = 0; j < numberOfInterpointsForNeighbor; ++ j)
        {
            ipoint = pointIndexForSend[j];
            for (int m = 0; m < nm; ++ m)
            {
                qtmp[icount ++] = q[m][ipoint];
            }
        }
        cdata->Write(qtmp, nlen * sizeof(RDouble));
        DelPointer(qtmp);
    }
}

void CFDSolver::TranslateInterfaceData(InterfaceInfo *iterfaceInfo, ActionKey *actkey, InterfaceDataProxy *interfaceDataProxy)
{
    if (iterfaceInfo == 0)
    {
        return;
    }

    vector <RDouble **> &vectorData = interfaceDataProxy->GetVectorData();
    vector <int> &vectorDimension = interfaceDataProxy->GetVectorDimension();

    //! According to the design,actkey->ipos stores zone i,the current zone is the j-th neighbor of zone i.
    //! Here we need to find out that which neighbor of the current zone is zone i,and the number is ineighbor.

    int iNeighbor = iterfaceInfo->FindIthNeighbor(actkey->ipos);
    int nIFaceOfNeighbor = iterfaceInfo->GetNIFaceOfNeighbor(iNeighbor);
    int *faceIndexForRecv = iterfaceInfo->GetFaceIndexForRecv(iNeighbor);

    //! Receive according to your needs.
    DataContainer *cData = actkey->GetData();
    cData->MoveToBegin();
    for (std::size_t iDim = 0; iDim < vectorDimension.size(); ++ iDim)
    {
        int iFace;
        int nm = vectorDimension[iDim];
        RDouble **q = vectorData[iDim];

        //! Bell 20130129 mod from Hexin.
        int nlen = nm * nIFaceOfNeighbor;

        RDouble *qTemp = new RDouble[nlen];

        cData->Read(qTemp, nlen * sizeof(RDouble));

        int icount = 0;
        for (int jFaceOfNeighbor = 0; jFaceOfNeighbor < nIFaceOfNeighbor; ++ jFaceOfNeighbor)
        {
            iFace = faceIndexForRecv[jFaceOfNeighbor];

            for (int m = 0; m < nm; ++ m)
            {
                q[m][iFace] = qTemp[icount ++];
            }
        }
        DelPointer(qTemp);
    }
}

void CFDSolver::TranslateInterpointData(InterpointInformation *interpointInformation, ActionKey *actkey, InterpointDataProxy * interpointDataProxy)
{
    if (interpointInformation == 0)
    {
        return;
    }
    vector<RDouble **> &vectorData = interpointDataProxy->GetVectorData();
    vector<int> &vectorDimension = interpointDataProxy->GetVectorDimension();

    //! According to the design,actkey->ipos stores zone i,the current zone is the j-th neighbor of zone i.
    //! Here we need to find out that which neighbor of the current zone is zone i,and the number is ineighbor.

    int iNeighbor = interpointInformation->FindIthNeighbor(actkey->ipos);
    int numberOfInterpointsForNeighbor = interpointInformation->GetNumberOfInterpointsForNeighbor(iNeighbor);
    int *pointIndexForReceive = interpointInformation->GetPointIndexForReceive(iNeighbor);

    //! Receive according to your needs.
    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();
    for (std::size_t iDim = 0; iDim < vectorDimension.size(); ++ iDim)
    {
        int iPoint;
        int nm = vectorDimension[iDim];
        RDouble **q = vectorData[iDim];

        int nlen = nm * numberOfInterpointsForNeighbor;

        RDouble *qtmp = new RDouble[nlen];

        cdata->Read(qtmp, nlen * sizeof(RDouble));

        int iCount = 0;
        for (int j = 0; j < numberOfInterpointsForNeighbor; ++ j)
        {
            iPoint = pointIndexForReceive[j];

            for (int m = 0; m < nm; ++ m)
            {
                q[m][iPoint] = qtmp[iCount ++];
            }
        }
        DelPointer(qtmp);
    }
}

streamsize CFDSolver::TranslateInterfaceDataLength(InterfaceInfo *iterfaceInfo, ActionKey *actkey, InterfaceDataProxy *interfaceDataProxy)
{
    streamsize nlen = 0;
    if (iterfaceInfo == 0)
    {
        return nlen;
    }
    vector <int> &vectorDimension = interfaceDataProxy->GetVectorDimension();

    int iNeighbor = iterfaceInfo->FindIthNeighbor(actkey->ipos);
    int nIFaceOfNeighbor = iterfaceInfo->GetNIFaceOfNeighbor(iNeighbor);

    for (std::size_t iDim = 0; iDim < vectorDimension.size(); ++ iDim)
    {
        int nm = vectorDimension[iDim];

        int nlenOfThisDim = nm * nIFaceOfNeighbor;
        nlen += nlenOfThisDim * sizeof(RDouble);
    }

    return nlen;
}

streamsize CFDSolver::TranslateInterpointDataLength(InterpointInformation *interpointInformation, ActionKey *actkey, InterpointDataProxy *interpointDataProxy)
{
    streamsize nlen = 0;
    if (interpointInformation == 0)
    {
        return nlen;
    }
    vector <int> &vectorDimension = interpointDataProxy->GetVectorDimension();

    int iNeighbor = interpointInformation->FindIthNeighbor(actkey->ipos);
    int nIpointOfNeighbor = interpointInformation->GetNumberOfInterpointsForNeighbor(iNeighbor);

    for (std::size_t iDim = 0; iDim < vectorDimension.size(); ++ iDim)
    {
        int nm = vectorDimension[iDim];

        int nlenOfThisDim = nm * nIpointOfNeighbor;
        nlen += nlenOfThisDim * sizeof(RDouble);
    }

    return nlen;
}

void CFDSolver::GetOversetInterfaceData(OversetInformationProxy *oversetInformationProxy, ActionKey *actkey, InterfaceDataProxy *interfaceDataProxy)
{
    vector <RDouble **> &vectorData = interfaceDataProxy->GetVectorData();
    vector <int> &vectorDimension = interfaceDataProxy->GetVectorDimension();

    OversetDataProxy *oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForSend();

    int neighborIndex = oversetDataProxy->FindNeighborIndex(actkey->ipos);

    int dim_sum = 0;

    for (std::size_t indexVector = 0; indexVector < vectorDimension.size(); ++ indexVector)
    {
        dim_sum += vectorDimension[indexVector];
    }

    int numberOfCellForSend = oversetDataProxy->GetNumberOfCells(neighborIndex);
    int *cellIndexMapForSend = oversetDataProxy->GetCellStorageIndexMapping(neighborIndex);
    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    for (std::size_t iDim = 0; iDim < vectorDimension.size(); ++ iDim)
    {
        int nm = vectorDimension[iDim];
        RDouble **q = vectorData[iDim];

        int nlen = nm * numberOfCellForSend;
        RDouble *qtmp = new RDouble [nlen];
        int icount = 0;

        for (int iCellLocal = 0; iCellLocal < numberOfCellForSend; ++ iCellLocal)
        {
            int iCell = cellIndexMapForSend[iCellLocal];
            for (int m = 0; m < nm; ++ m)
            {
                qtmp[icount++] = q[m][iCell];
            }
        }
        cdata->Write(qtmp, nlen * sizeof(RDouble));
        DelPointer(qtmp);
    }
}

void CFDSolver::TranslateOversetInterfaceData(OversetInformationProxy *oversetInformationProxy, ActionKey *actkey, InterfaceDataProxy *interfaceDataProxy)
{
    vector < RDouble **> &vectorData = interfaceDataProxy->GetVectorData();
    vector <int> &vectorDimension = interfaceDataProxy->GetVectorDimension();

    //! According to the design,actkey->ipos stores zone i,the current zone is the j-th neighbor of zone i.
    //! Here we need to find out that which neighbor of the current zone is zone i,and the number is ineighbor.

    OversetDataProxy *oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForReceive();
    int iNeighbor = oversetDataProxy->FindNeighborIndex(actkey->ipos);

    int numberOfCells = oversetDataProxy->GetNumberOfCells(iNeighbor);
    int *cellStorageIndexMapping = oversetDataProxy->GetCellStorageIndexMapping(iNeighbor);

    //! Receive according to your needs.
    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    for (std::size_t iDim = 0; iDim < vectorDimension.size(); ++ iDim)
    {
        int nm = vectorDimension[iDim];
        RDouble **q = vectorData[iDim];

        for (int iCellLocal = 0; iCellLocal < numberOfCells; ++ iCellLocal)
        {
            int iCell = cellStorageIndexMapping[iCellLocal];

            for (int m = 0; m < nm; ++ m)
            {
                cdata->Read(&q[m][iCell], sizeof(RDouble));
            }
        }
    }
}

void CFDSolver::PutCorrection(Grid *grid, FieldProxy *fieldProxy, FieldProxy *oldFieldProxy, int totalEqn)
{
    if (grid->Type() == STRUCTGRID)
    {
        RDouble4D &field = fieldProxy->GetField_STR();
        RDouble4D &oldfield = oldFieldProxy->GetField_STR();

        oldfield -= field;
    }
    else
    {
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal = nTotalCell + nBoundFace;

        RDouble **field = fieldProxy->GetField_UNS();
        RDouble **oldfield = oldFieldProxy->GetField_UNS();

        for (int indexEquation = 0; indexEquation < totalEqn; ++ indexEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                oldfield[indexEquation][iCell] -= field[indexEquation][iCell];
            }
        }
    }
}

void CFDSolver::ComputeDualTimeCoefficient(int methodOfDualTime, RDouble *dualTimeCoefficient)
{
    RDouble theta = 0.0, xi = 0.0, phi = 0.0;
    const int EULEREXPLICIT           = 0; //! One order.
    const int BACKWARDEULER           = 1; //! One order, A-stable.
    const int ONESTEPTRAPEZOIDAL      = 2; //! Two order, A-stable, also called Crank-Nicolson.
    const int BACKWARDDIFFERENTIATION = 3; //! Two order, A-stable.
    const int ADAMSTYPE               = 4; //! Two order, A-stable.
    const int LEESTYPE                = 5; //! Two order, A-stable.
    const int TWOSTEPTRAPEZOIDAL      = 6; //! Two order, A-stable.
    const int LEAPFROG                = 7; //! Two order.
    const int ADAMSBASHFORTH          = 8; //! Two order, forward difference.
    const int THIRDORDERIMPLICIT      = 9; //! Three order.
    const int ADAMSMOULTON            = 10;//! Three order.
    const int MILNE                   = 11;//! Four order.

    switch (methodOfDualTime)
    {
        case EULEREXPLICIT:
            theta =  0.0;
            xi    =  0.0;
            phi   =  0.0;
            break;
        case BACKWARDEULER:
            theta =  1.0;
            xi    =  0.0;
            phi   =  0.0;
            break;
        case ONESTEPTRAPEZOIDAL:
            theta =  0.5;
            xi    =  0.0;
            phi   =  0.0;
            break;
        case BACKWARDDIFFERENTIATION:
            theta =  1.0;
            xi    =  0.5;
            phi   =  0.0;
            break;
        case ADAMSTYPE:
            theta =  0.75;
            xi    =  0.0;
            phi   = -0.25;
            break;
        case LEESTYPE:
            theta =  1.0/3.0;
            xi    = -1.0/2.0;
            phi   = -1.0/3.0;
            break;
        case TWOSTEPTRAPEZOIDAL:
            theta =  1.0/2.0;
            xi    = -1.0/2.0;
            phi   = -1.0/2.0;
            break;
        case LEAPFROG:
            theta =  0.0;
            xi    = -0.5;
            phi   =  0.0;
            break;
        case ADAMSBASHFORTH:
            theta =  0.0;
            xi    =  0.0;
            phi   =  0.5;
            break;
        case THIRDORDERIMPLICIT:
            theta =  1.0/3.0;
            xi    = -1.0/6.0;
            phi   =  0.0;
            break;
        case ADAMSMOULTON:
            theta =  5.0/12.0;
            xi    =  0.0;
            phi   =  1.0/12.0;
            break;
        case MILNE:
            theta =  1.0/6.0;
            xi    = -1.0/2.0;
            phi   = -1.0/6.0;
            break;
        default:
            TK_Exit::UnexpectedVarValue("methodOfDualTime", methodOfDualTime);
            break;
    }

    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    //! Computation of dualtime coefficients, including three coefficients for Residual of R(p), R(n) and R(n-1),
    //! and three coefficients for conservative variables of qcsv(p), qcsv(n), and qcsv(n-1),
    //! the fourth coefficient is also used for the real    time part of the diagonal data for the matrix of the left term when LUSGS is running,
    //! the last   coefficient is      used for the virtual time part of the diagonal data for the matrix of the left term when LUSGS is running.
    if (methodOfDualTime == EULEREXPLICIT || methodOfDualTime == LEAPFROG || methodOfDualTime == ADAMSBASHFORTH)   //! theta = 0, explicit method.
    {
        dualTimeCoefficient[0] = 0.0;
        dualTimeCoefficient[1] = 1.0 + phi;
        dualTimeCoefficient[2] =     - phi;
        dualTimeCoefficient[3] = - (1.0 +       xi) / physicalTimeStep;
        dualTimeCoefficient[4] =   (1.0 + 2.0 * xi) / physicalTimeStep; //! dualTimeCoefficient[ 4 ] = - dualTimeCoefficient[ 3 ] - dualTimeCoefficient[ 5 ];
        dualTimeCoefficient[5] =               - xi   / physicalTimeStep;
        dualTimeCoefficient[6] = 1.0;
    }
    else //! theta /= 0, implicit method
    {
        dualTimeCoefficient[0] = 1.0;
        dualTimeCoefficient[1] = (1.0 - theta + phi) / theta;
        dualTimeCoefficient[2] =               - phi   / theta;
        dualTimeCoefficient[3] = - (1.0 +       xi)  / theta / physicalTimeStep;
        dualTimeCoefficient[4] =   (1.0 + 2.0 * xi)  / theta / physicalTimeStep; //! dualTimeCoefficient[ 4 ] = - dualTimeCoefficient[ 3 ] - dualTimeCoefficient[ 5 ];
        dualTimeCoefficient[5] =               - xi    / theta / physicalTimeStep;
        dualTimeCoefficient[6] = 1.0                   / theta;
    }
}

void CFDSolver::DumpResidual(ActionKey *actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0) 
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    actkey->filename = GetResidualFileName();
    if (actkey->filename == "")
    {
        return;
    }

    ios_base::openmode openmode = ios_base::out|ios_base::app;
    actkey->openmode = openmode;

    InitMaxResidual();

    //! First: Initialize Residual, set residual on Each Zone equal zero.
    PHSPACE::InitResidual();

    //! Second: Computing residual on Each Zone.
    int nZones = GetNumberofGlobalZones();
    actkey->GetData()->MoveToBegin();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int zoneProcessorID = GetZoneProcessorID(iZone);
        int currentProcessorID = GetCurrentProcessorID();
        if (currentProcessorID == zoneProcessorID)
        {
            //! If zone i belongs to the current process.
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->GetResidual(actkey);
        }
    }

    //! Third: Collecting residual on Each Zone to serve.
    PHSPACE::CollectResidual(actkey);

    //! Forth: Output residual on serve.
    int currentProcessorID = GetCurrentProcessorID();
    if (currentProcessorID == server)
    {
        PHSPACE::PostDumpResidual(actkey);
    }
}

void CFDSolver::CheckResidual(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int nCurrentProcessorID = GetCurrentProcessorID();
    RDouble localMaxResidual = 0.0, localZoneResidual = 0.0;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int zoneProcessorID = GetZoneProcessorID(iZone);
        if (zoneProcessorID == nCurrentProcessorID)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->GetResidual(actkey, localZoneResidual);
            if (localZoneResidual > localMaxResidual)
            {
                localMaxResidual = localZoneResidual;
            }
        }
    }

    RDouble globalMaxResidual = localMaxResidual;
    int nProcessorNumber = GetNumberOfProcessor();
    if (nProcessorNumber > 1)
    {
        //! Synchronization of communication.
        MPI_Barrier(MPI_COMM_WORLD);

        //! Communication.
        MPI_Allreduce(&localMaxResidual, &globalMaxResidual, 1, PH_MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    GlobalDataBase::UpdateData("monitorResidual", &globalMaxResidual, PHDOUBLE, 1);
}

void CFDSolver::ObtainCFLNumber(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int nCurrentProcessorID = GetCurrentProcessorID();
    RDouble localMaxCFL = 0.0, localZoneCFLMax = 0.0, localMinCFL = LARGE, localZoneCFLMin = 0.0;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int zoneProcessorID = GetZoneProcessorID(iZone);
        if (zoneProcessorID == nCurrentProcessorID)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->ObtainCFLNumber(GetGrid(actkey->level), localZoneCFLMin, localZoneCFLMax);
            if (localZoneCFLMin < localMinCFL)
            {
                localMinCFL = localZoneCFLMin;
            }
            if (localZoneCFLMax > localMaxCFL)
            {
                localMaxCFL = localZoneCFLMax;
            }
        }
    }

    RDouble globalMinCFL = localMinCFL;
    RDouble globalMaxCFL = localMaxCFL;
    int nProcessorNumber = GetNumberOfProcessor();
    if (nProcessorNumber > 1)
    {
        //! Synchronization of communication.
        MPI_Barrier(MPI_COMM_WORLD);

        //! Communication.
        MPI_Allreduce(&localMinCFL, &globalMinCFL, 1, PH_MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        MPI_Allreduce(&localMaxCFL, &globalMaxCFL, 1, PH_MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    GlobalDataBase::UpdateData("globalMinCFL", &globalMinCFL, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("globalMaxCFL", &globalMaxCFL, PHDOUBLE, 1);
}

void CFDSolver::ComputeSurfaceHeatingChange(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int nCurrentProcessorID = GetCurrentProcessorID();
    RDouble localMaxHeatChange = 0.0, localZoneHeatChange = 0.0;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int zoneProcessorID = GetZoneProcessorID(iZone);
        if (zoneProcessorID == nCurrentProcessorID)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->GetSurfaceHeatingChange(actkey, localZoneHeatChange);
            if (localZoneHeatChange > localMaxHeatChange)
            {
                localMaxHeatChange = localZoneHeatChange;
            }
        }
    }

    RDouble globalMaxHeatChange = localMaxHeatChange;
    int nProcessorNumber = GetNumberOfProcessor();
    if (nProcessorNumber > 1)
    {
        //! Synchronization of communication.
        MPI_Barrier(MPI_COMM_WORLD);

        //! Communication.
        MPI_Allreduce(&localMaxHeatChange, &globalMaxHeatChange, 1, PH_MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    GlobalDataBase::UpdateData("monitorHeatingChange", &globalMaxHeatChange, PHDOUBLE, 1);
}

//! Get the file name for Residual dumping, NS-solver by default.
const string CFDSolver::GetResidualFileName()
{
    string resSaveFile = "res.dat";
    GlobalDataBase::GetData("resSaveFile", &resSaveFile, PHSTRING, 1);
    return resSaveFile;
}

void CFDSolver::DumpRestartData(ActionKey *actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0) 
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    actkey->filename = GetRestartFileName();
    if (actkey->filename == "")
    {
        return;
    }

    if (PHMPI::IsParallelRun())
    {
        actkey->filename = PHSPACE::AddSymbolToFileName(actkey->filename, "_", 0);
    }

    hid_t file;
    file = CreateHDF5File(actkey->filename);
    actkey->filepos = file;

    CreateH5RestartFile(actkey);

    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessID = GetZoneProcessorID(iZone);
        if (currentProcessorID == sendProcessID)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->DumpRestartH5(actkey);
        }
    }

    H5Fclose(file);
    actkey->filepos = 0;
}

void CFDSolver::DumpProtectData(ActionKey *actkey, int &RestartFile)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0) 
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    string protectionFile0 = "./results/flow0.dat";
    GlobalDataBase::GetData("protectionFile0", &protectionFile0, PHSTRING, 1);

    string protectionFile1 = "./results/flow1.dat";
    GlobalDataBase::GetData("protectionFile1", &protectionFile1, PHSTRING, 1);

    if (PHMPI::IsParallelRun())
    {
        protectionFile0 = PHSPACE::AddSymbolToFileName(protectionFile0, "_", 0);
        protectionFile1 = PHSPACE::AddSymbolToFileName(protectionFile1, "_", 0);
    }

    ifstream infile_0(protectionFile0.c_str(), ios::in);
    if (infile_0)
    {
        RestartFile = 1;
    }

    ifstream infile_1(protectionFile1.c_str(), ios::in);
    if (infile_1)
    {
        RestartFile = 0;
    }

    if (RestartFile == 0)
    {
        actkey->filename = protectionFile0;
    }
    else
    {
        actkey->filename = protectionFile1;
    }

    if (actkey->filename == "")
    {
        return;
    }

    hid_t file;
    file = CreateHDF5File(actkey->filename);
    actkey->filepos = file;

    CreateH5RestartFile(actkey);

    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessID = GetZoneProcessorID(iZone);
        if (currentProcessorID == sendProcessID)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->DumpRestartH5(actkey);
        }
    }

    H5Fclose(file);
    actkey->filepos = 0;
}

void CFDSolver::DumpTurbProtectData(ActionKey *actkey, int &RestartFile)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0) 
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    string protectionTurbFile0 = "./results/turb0.dat";
    GlobalDataBase::GetData("protectionTurbFile0", &protectionTurbFile0, PHSTRING, 1);

    string protectionTurbFile1 = "./results/turb1.dat";
    GlobalDataBase::GetData("protectionTurbFile1", &protectionTurbFile1, PHSTRING, 1);

    if (PHMPI::IsParallelRun())
    {
        protectionTurbFile0 = PHSPACE::AddSymbolToFileName(protectionTurbFile0, "_", 0);
        protectionTurbFile1 = PHSPACE::AddSymbolToFileName(protectionTurbFile1, "_", 0);
    }

    ifstream infile_0(protectionTurbFile0.c_str(), ios::in);
    if (infile_0)
    {
        RestartFile = 1;
    }

    ifstream infile_1(protectionTurbFile1.c_str(), ios::in);
    if (infile_1)
    {
        RestartFile = 0;
    }

    if (RestartFile == 0)
    {
        actkey->filename = protectionTurbFile0;
    }
    else
    {
        actkey->filename = protectionTurbFile1;
    }

    if (actkey->filename == "")
    {
        return;
    }

    hid_t file;
    file = CreateHDF5File(actkey->filename);
    actkey->filepos = file;

    CreateH5RestartFile(actkey);

    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessID = GetZoneProcessorID(iZone);
        if (currentProcessorID == sendProcessID)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->DumpRestartH5(actkey);
        }
    }

    H5Fclose(file);
    actkey->filepos = 0;
}

void CFDSolver::DumpTransitionProtectData(ActionKey *actkey, int &RestartFile)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0) 
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    string protectionTransitionFile0 = "./results/transition0.dat";
    GlobalDataBase::GetData("protectionTransitionFile0", &protectionTransitionFile0, PHSTRING, 1);

    string protectionTransitionFile1 = "./results/transition1.dat";
    GlobalDataBase::GetData("protectionTransitionFile1", &protectionTransitionFile1, PHSTRING, 1);

    if (PHMPI::IsParallelRun())
    {
        protectionTransitionFile0 = PHSPACE::AddSymbolToFileName(protectionTransitionFile0, "_", 0);
        protectionTransitionFile1 = PHSPACE::AddSymbolToFileName(protectionTransitionFile1, "_", 0);
    }

    ifstream infile_0(protectionTransitionFile0.c_str(), ios::in);
    if (infile_0)
    {
        RestartFile = 1;
    }

    ifstream infile_1(protectionTransitionFile1.c_str(), ios::in);
    if (infile_1)
    {
        RestartFile = 0;
    }

    if (RestartFile == 0)
    {
        actkey->filename = protectionTransitionFile0;
    }
    else
    {
        actkey->filename = protectionTransitionFile1;
    }

    if (actkey->filename == "")
    {
        return;
    }

    hid_t file;
    file = CreateHDF5File(actkey->filename);
    actkey->filepos = file;

    CreateH5RestartFile(actkey);

    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessID = GetZoneProcessorID(iZone);
        if (currentProcessorID == sendProcessID)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->DumpRestartH5(actkey);
        }
    }

    H5Fclose(file);
    actkey->filepos = 0;
}

void CFDSolver::CreateH5RestartFile(ActionKey *actkey)
{
    Param_CFDSolver *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    int isWennScheme = parameters->GetWennSchemeFlag();

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ostringstream oss;
        oss << "Group" << iZone;

        string grpName;

        grpName = oss.str();
        CreateGroup(actkey->filepos, grpName);

        int length = 0;
        int nWallBC = 0, nCurMaxCellNum = 0;
        int sendProcessID = GetZoneProcessorID(iZone);
        if (currentProcessorID == sendProcessID)
        {
            int gridTypeofLocalZone = GetZoneGridType(iZone);
            if (gridTypeofLocalZone == STRUCTGRID)
            {
                StructGrid *grid = StructGridCast(PHSPACE::GetGrid(iZone, actkey->level));
                int ni = grid->GetNI();
                int nj = grid->GetNJ();
                int nk = grid->GetNK();

                if (GetDim() == THREE_D)
                {
                    length = (ni + 3) * (nj + 3) * (nk + 3);
                }
                else
                {
                    length = (ni + 3) * (nj + 3);
                }

                if (isWennScheme == 1)
                {
                    if (GetDim() == THREE_D)
                    {
                        length = (ni + 5) * (nj + 5) * (nk + 5);
                    }
                    else
                    {
                        length = (ni + 5) * (nj + 5);
                    }
                }
                //! Count the cells on wall.
                grid->GetWallBCInfo(nWallBC, nCurMaxCellNum);
            }
            else
            {
                Grid *grid = PHSPACE::GetGrid(iZone, actkey->level);
                int nBoundFace = grid->GetNBoundFace();
                int nTotalCell = grid->GetNTotalCell();

                length = nBoundFace + nTotalCell;
            }

        }

        int getherLength;
        MPI_Allreduce(&length, &getherLength, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        int nTotalWallCell, nCurTotalWallCell;
        nCurTotalWallCell = nWallBC * nCurMaxCellNum;
        MPI_Allreduce(&nCurTotalWallCell, &nTotalWallCell, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        hid_t  grpData;
        grpData = OpenGroup(actkey->filepos, grpName);

        if (actkey->solver == NS_EQUATION)
        {
            int nEquation = GetNumberOfEquations();
            CreateEmptyData(grpData, "q", nEquation, getherLength, PHDOUBLE);

            int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
            int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
            int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
            int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");

            //! Create the memories of temperatures.
            //if (nChemical > 0)
            //{
                CreateEmptyData(grpData, "t", nTemperatureModel, getherLength, PHDOUBLE);
            //}

            if (nTotalWallCell > 0)    //! The current grid includes the wall boundary conditions.
            {
                //RDouble wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
                CreateEmptyData(grpData, "surfaceTemperature", 1, nTotalWallCell, PHDOUBLE);

                if (nChemical > 0 && numberOfSpecies > 0)
                {
                    CreateEmptyData(grpData, "surfaceMassFraction", numberOfSpecies, nTotalWallCell, PHDOUBLE);
                }

                if (nSlipBCModel > 0)
                {
                    CreateEmptyData(grpData, "surfaceSlipVariables", numberOfSpecies + 6, nTotalWallCell, PHDOUBLE);
                }
            }

            if (isUnsteady == 1)
            {
                CreateEmptyData(grpData, "q_unsteady_n1", nEquation, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "q_unsteady_n2", nEquation, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "res", nEquation, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "res_unsteady_n1", nEquation, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "res_unsteady_n2", nEquation, getherLength, PHDOUBLE);

                CreateEmptyData(grpData, "qAverage", nEquation, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "tauAverage", 6, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "q2Average", 6, getherLength, PHDOUBLE);
            }
            CreateEmptyData(grpData, "surfaceTemperatureFlag", 1, 1, PHDOUBLE);
            CreateEmptyData(grpData, "catalyticCoefFlag", 1, 1, PHDOUBLE);

            CreateEmptyData(grpData, "nIsChemicalFlow", 1, 1, PHINT);
            CreateEmptyData(grpData, "numberOfSpecies", 1, 1, PHINT);
            CreateEmptyData(grpData, "nTemperatureModel", 1, 1, PHINT);
            CreateEmptyData(grpData, "speciesNameList", 16, 1, PHINT);
        }
        else if (actkey->solver == TURBULENCE)
        {
            int n_turb   = GlobalDataBase::GetIntParaFromDB("n_turb");
            CreateEmptyData(grpData, "q_turb", n_turb, getherLength, PHDOUBLE);
            CreateEmptyData(grpData, "visturb", 1    , getherLength, PHDOUBLE);

            int transitionType   = GlobalDataBase::GetIntParaFromDB("transitionType");
            if (transitionType == IREGAMA)
            {
                CreateEmptyData(grpData, "SpSdRatio", 1, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "gamaeff", 1    , getherLength, PHDOUBLE);
            }

            if (isUnsteady == 1)
            {
                CreateEmptyData(grpData, "q_turb_unsteady_n1", n_turb, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "q_turb_unsteady_n2", n_turb, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "res_turb", n_turb, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "res_turb_unsteady_n1", n_turb, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "res_turb_unsteady_n2", n_turb, getherLength, PHDOUBLE);

                CreateEmptyData(grpData, "qAverageTurb", n_turb, getherLength, PHDOUBLE);
            }
        }
        else if (actkey->solver == TRANSITION)
        {
            int n_transition   = GlobalDataBase::GetIntParaFromDB("n_transition");
            CreateEmptyData(grpData, "q_transition", n_transition, getherLength, PHDOUBLE);

            if (isUnsteady == 1)
            {
                CreateEmptyData(grpData, "q_transition_unsteady_n1", n_transition, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "q_transition_unsteady_n2", n_transition, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "res_transition", n_transition, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "res_transition_unsteady_n1", n_transition, getherLength, PHDOUBLE);
                CreateEmptyData(grpData, "res_transition_unsteady_n2", n_transition, getherLength, PHDOUBLE);

                CreateEmptyData(grpData, "qAverageTransition", n_transition, getherLength, PHDOUBLE);
            }
        }

        CreateEmptyData(grpData, "nTotalCell", 1, 1, PHINT);

        H5Gclose(grpData);
    }

    CreateEmptyData(actkey->filepos, "Version", 1, 1, PHINT);
    CreateEmptyData(actkey->filepos, "outnstep", 1, 1, PHINT);
    if (isUnsteady == 1 && actkey->solver == NS_EQUATION)
    {
        CreateEmptyData(actkey->filepos, "physicalTime", 1, 1, PHDOUBLE);
    }

    CreateEmptyData(actkey->filepos, "nStatisticalStep", 1, 1, PHINT);
}

void CFDSolver::InitFlowAsReadingRestart(ActionKey *actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    actkey->filename = GetRestartFileName();
    if (actkey->filename == "")
    {
        return;
    }

    if (PHMPI::IsParallelRun())
    {
        actkey->filename = PHSPACE::AddSymbolToFileName(actkey->filename, "_", 0);
    }

    hid_t file;
    file = OpenHDF5File(actkey->filename);
    actkey->filepos = file;

    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int recvProcess = GetZoneProcessorID(iZone);

        if (currentProcessorID == recvProcess)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->ReadRestartH5(actkey);
        }
    }

    H5Fclose(file);
    actkey->filepos = 0;
}

void CFDSolver::InitFlowByReadingInterpolate(ActionKey* actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;
  
    actkey->gridfilename = GetRestartGridFileName();

    actkey->filename = GetRestartFlowFileName();

    if (actkey->gridfilename == "" || actkey->filename == "")
    {
        return;
    }

    hid_t gridfile;
    gridfile = OpenHDF5File(actkey->gridfilename);
    actkey->gridfilepos = gridfile;

    hid_t file;
    file = OpenHDF5File(actkey->filename);
    actkey->filepos = file;

    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int recvProcess = GetZoneProcessorID(iZone);

        if (currentProcessorID == recvProcess)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->InterpolateFromRestartH5(actkey);
        }
    }

    H5Fclose(gridfile);
    actkey->gridfilepos = 0;

    H5Fclose(file);
    actkey->filepos = 0;

}
void CFDSolver::InitFlowAsReadingProtectedRestart(ActionKey *actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;
    int RestartFile = 0;
    string protectionFile0 = "./results/flow0.dat";
    GlobalDataBase::GetData("protectionFile0", &protectionFile0, PHSTRING, 1);

    string protectionFile1 = "./results/flow1.dat";
    GlobalDataBase::GetData("protectionFile1", &protectionFile1, PHSTRING, 1);

    if (PHMPI::IsParallelRun())
    {
        protectionFile0 = PHSPACE::AddSymbolToFileName(protectionFile0, "_", 0);
        protectionFile1 = PHSPACE::AddSymbolToFileName(protectionFile1, "_", 0);
    }

    ifstream infile_0(protectionFile0.c_str(), ios::in);
    if (infile_0)
    {
        RestartFile = 0;
    }

    ifstream infile_1(protectionFile1.c_str(), ios::in);
    if (infile_1)
    {
        RestartFile = 1;
    }

    if (infile_0 && infile_1)     //! If both flow0 and flow1 exist,determine which file to use by the last modified time.
    {
        RDouble TimeFile0, TimeFile1;
        TIME_SPACE::GetFileModifiedTime(protectionFile0, TimeFile0);
        TIME_SPACE::GetFileModifiedTime(protectionFile1, TimeFile1);
        if (TimeFile0 < TimeFile1)
        {
            RestartFile = 0;
        }
        else if (TimeFile0 > TimeFile1)
        {
            RestartFile = 1;
        }
    }

    if (RestartFile == 0)
    {
        actkey->filename = protectionFile0;
    }
    else
    {
        actkey->filename = protectionFile1;
    }

    if (actkey->filename == "")
    {
        return;
    }

    hid_t file;
    file = OpenHDF5File(actkey->filename);
    actkey->filepos = file;

    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int recvProcess = GetZoneProcessorID(iZone);

        if (currentProcessorID == recvProcess)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->ReadRestartH5(actkey);
        }
    }

    H5Fclose(file);
    actkey->filepos = 0;
}

//! Get the file name for Restart info dumping, NS-solver by default.
const string CFDSolver::GetRestartFileName()
{
    string restartNSFile = "flow.dat";
    GlobalDataBase::GetData("restartNSFile", &restartNSFile, PHSTRING, 1);
    return restartNSFile;
}
//! Get the grid file name for Restart info by interpolated.
const string CFDSolver::GetRestartGridFileName()
{
    string restartGridFile = "";
    GlobalDataBase::GetData("restartGridFile", &restartGridFile, PHSTRING, 1);
    return restartGridFile;
}

//! Get the flow file name for Restart info by interpolated.
const string CFDSolver::GetRestartFlowFileName()
{
    string restartNSVarFile = "";
    GlobalDataBase::GetData("restartNSVarFile", &restartNSVarFile, PHSTRING, 1);
    return restartNSVarFile;
}

void CFDSolver::DumpAirForceCoef(ActionKey *actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0) 
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    string aircoeffile = GetAirCoefFileName();
    if (aircoeffile == "")
    {
        return;
    }

    actkey->filename = aircoeffile;
    ios_base::openmode openmode = ios_base::out | ios_base::app;
    actkey->openmode = openmode;

    //! First: Initialize AirForceCoef: allocate array for AirForceCoef variables.
    PHSPACE::InitForce();

    //! Second: Computing AirForceCoef on Each Zone.
    int nZones = GetNumberofGlobalZones();
    actkey->GetData()->MoveToBegin();
    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        int ZoneProcessorID = GetZoneProcessorID(iZone);
        int currentProcessorID = GetCurrentProcessorID();

        if (currentProcessorID == ZoneProcessorID)
        {
            //! If zone i belongs to the current process.
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->AirForceCoef(actkey);
        }
    }

    //! Third: Collecting AirForceCoef on Each Zone to serve.
    PHSPACE::CollectionForce(actkey);

    //! Forth: Output AirForceCoef to aircoeffile on serve.
    int currentProcessorID = GetCurrentProcessorID();
    if (currentProcessorID == server)
    {
        PHSPACE::DumpForce(actkey);
    }
}

//! Get the file name for AirCoef dumping, NS-solver by default.
const string CFDSolver::GetAirCoefFileName()
{
    string airCoefFile = "aircoef.dat";
    GlobalDataBase::GetData("aircoeffile", &airCoefFile, PHSTRING, 1);
    return airCoefFile;
}

void CFDSolver::DumpCpDistri(ActionKey *actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0) 
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    string wall_aircoefile   = "wall_aircoef.dat";
    GlobalDataBase::GetData("wall_aircoefile", &wall_aircoefile, PHSTRING, 1);
    if (wall_aircoefile == "")
    {
        return;
    }

    actkey->filename = wall_aircoefile;

    vector <DataContainer *> datalist;
    datalist.resize(0);

    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcess = GetZoneProcessorIDSepMode(iZone);
        int recvProcess = GetServerSepMode();

        int currentProcessorID = GetCurrentProcessorID();
        int sendRecvTag = GetSendRecvTag(actkey, iZone);

        if (currentProcessorID == sendProcess)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->CpDistriCoef(actkey);
            int TecioMission = 0;
            PHWrite(actkey->GetData(), &TecioMission, 1);
        }

        PH_Trade(actkey, sendProcess, recvProcess, sendRecvTag);

        if (currentProcessorID == recvProcess)
        {
            datalist.push_back(actkey->GetData());
            DataContainer *cdata0 = new DataContainer();
            actkey->SetData(cdata0);
        }
    }

    int currentProcessorID = GetCurrentProcessorID();
    int serverTmp = GetServerSepMode();

    if (currentProcessorID == serverTmp)
    {
        int visualfileType = GlobalDataBase::GetIntParaFromDB("visualfileType");
        if (visualfileType == 0)
        {
            #ifdef USE_TecplotLib
                WriteWallAircoef(actkey, datalist);
            #else
                ostringstream WarningInfo;
                WarningInfo << "Warning : visualfileType = 0 is not right, dump wall aircoef data into ASCII file! \n";
                PrintToWindow(WarningInfo);

                WriteWallAircoefASCII(actkey, datalist);
            #endif
        }
        //else if (visualfileType == 1)
        else
        {
            WriteWallAircoefASCII(actkey, datalist);
        }
        /*else
        {
            TK_Exit::UnexpectedVarValue("visualfileType = ", visualfileType);
        }*/
    }

    for (unsigned int iZone = 0 ; iZone < datalist.size(); ++ iZone)
    {
        FreePointer(datalist[iZone]);
    }
}

void CFDSolver::DumpSurfaceInfo(ActionKey *actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;
    string wall_varfile = "surface_tec.dat";
    GlobalDataBase::GetData("wall_varfile", &wall_varfile, PHSTRING, 1);

    if (wall_varfile == "")
    {
        return;
    }

    actkey->filename = wall_varfile;
    ios_base::openmode openmode = ios_base::out|ios_base::trunc;
    actkey->openmode = openmode;

    int currentProcessorID = GetCurrentProcessorID();
    int serverTmp = GetServerSepMode();
 
    if (currentProcessorID == serverTmp)
    {
        ParallelOpenFile(actkey);
    }

    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcess = GetZoneProcessorIDSepMode(iZone);
        int recvProcess = GetServerSepMode();

        int myProcessorID = GetCurrentProcessorID();
        int sendRecvTag = GetSendRecvTag(actkey, iZone);

        if (myProcessorID == sendProcess)
        {
            //! If zone i belongs to the current process.
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->ExportSurfaceVariables(actkey);
        }

        PH_Trade(actkey, sendProcess, recvProcess, sendRecvTag);

        if (myProcessorID == recvProcess)
        {
            //! If the current process is the server process.
            WriteASCIIFile(actkey);
        }
    }

    if (currentProcessorID == serverTmp)
    {
        ParallelCloseFile(actkey);
    }
}

void CFDSolver::WriteWallAircoefASCII(ActionKey *actkey, vector <DataContainer *> datalist)
{
    int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");

    ostringstream ossTitle;
    ossTitle << "title=\"Wall Aircoef Data of PHengLEI\"" << endl;
    ossTitle << "variables=\"x\", \"y\", \"z\"";

    int nVariableNumber = postVisualWall->GetVisualVariablesNumber();
    int dataNumber = nVariableNumber;
    string varName;
    int varType;
    for (int iVar = 0; iVar < nVariableNumber; ++ iVar)
    {
        varType = postVisualWall->GetVisualVariablesType(iVar);
        if (varType == VISUAL_WALL_NS)
        {
            dataNumber += numberOfSpecies - 1;
        }
        else
        {
            varName = postVisualWall->GetVariableName(varType);
            ossTitle << ", \"" << varName.c_str() << "\"";
        }
    }

    //! Export the mass fractions of species.
    if (nVariableNumber != dataNumber)
    {
        using namespace GAS_SPACE;
        string *speciesName = gas->GetNameOfSpecies();
        for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
        {
            ossTitle << ", \"" << speciesName[iSpecies].c_str() << "\"";
        }
    }

    string flowFileName = actkey->filename;
    fstream wallAircoefFile;
    OpenFile(wallAircoefFile, flowFileName, ios_base::out|ios_base::trunc);
    wallAircoefFile << setiosflags(ios::left);
    wallAircoefFile << setiosflags(ios::scientific);
    wallAircoefFile << setprecision(10);
    wallAircoefFile << ossTitle.str() << endl;

    uint_t nzones = datalist.size();
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        int gridTypeofLocalZone = PHMPI::GetZoneGridType(iZone);

        DataContainer *cdata = datalist[iZone];
        cdata->MoveToBegin();

        int IMxOrNumPts, JMxOrNumElements, KMxOrNumFaces;

        int oriGridIndex = -1;
        PHRead(cdata, &oriGridIndex, 1);

        int GridID;
        PHRead(cdata, &GridID, 1);

        int TecioMission;
        PHRead(cdata, &TecioMission, 1);

        while (true)
        {
            if (TecioMission == 0) break;

            string bcName;
            cdata->ReadString(bcName);

            PHRead(cdata, &IMxOrNumPts, 1);
            PHRead(cdata, &JMxOrNumElements, 1);
            PHRead(cdata, &KMxOrNumFaces, 1);

            int coor_count = 0;
            PHRead(cdata, coor_count);

            RDouble *x = new RDouble[coor_count];
            RDouble *y = new RDouble[coor_count];
            RDouble *z = new RDouble[coor_count];
            PHRead(cdata, x, coor_count);
            PHRead(cdata, y, coor_count);
            PHRead(cdata, z, coor_count);

            vector <int> data_count;
            data_count.resize(0);
            vector <RDouble *> vector_data;
            vector_data.resize(0);
            for (int iData = 0; iData < dataNumber; ++ iData)
            {
                int idata_count = 0;
                PHRead(cdata, idata_count);
                data_count.push_back(idata_count);

                int iValueLocation = 0;
                PHRead(cdata, iValueLocation);

                RDouble *qq = new RDouble [idata_count];
                PHRead(cdata, qq, idata_count);

                vector_data.push_back(qq);
            }

            ostringstream ossZoneTitle;

            if (oriGridIndex >= 0)
            {
                ossZoneTitle << "\"oriZone_" << oriGridIndex << " ";
                ossZoneTitle << "Zone_" << GridID << "_" << bcName << "\"";
            }
            else
            {
                ossZoneTitle << "\"Zone_" << GridID << "_" << bcName << "\"";
            }

            wallAircoefFile << "zone T = " << ossZoneTitle.str() << "\n";

            if (gridTypeofLocalZone == STRUCTGRID)
            {
                wallAircoefFile << "I = " << IMxOrNumPts << "\n";
                wallAircoefFile << "J = " << JMxOrNumElements << "\n";
                wallAircoefFile << "K = " << KMxOrNumFaces << "\n";
                wallAircoefFile << "f = BLOCK" << "\n";

                int wordWidth = 20;
                int nWordOfLine = 5;
                for (int iNode = 0; iNode < coor_count; ++ iNode)
                {
                    wallAircoefFile << setw(wordWidth) << x[iNode];
                    if ((iNode + 1) % nWordOfLine == 0) 
                    {
                        wallAircoefFile << "\n";
                    }
                }
                if (coor_count % nWordOfLine != 0) wallAircoefFile << "\n";

                for (int iNode = 0; iNode < coor_count; ++ iNode)
                {
                    wallAircoefFile << setw(wordWidth) << y[iNode];
                    if ((iNode + 1) % nWordOfLine == 0) 
                    {
                        wallAircoefFile << "\n";
                    }
                }
                if (coor_count % nWordOfLine != 0) wallAircoefFile << "\n";

                for (int iNode = 0; iNode < coor_count; ++ iNode)
                {
                    wallAircoefFile << setw(wordWidth) << z[iNode];
                    if ((iNode + 1) % nWordOfLine == 0) 
                    {
                        wallAircoefFile << "\n";
                    }
                }
                if (coor_count % nWordOfLine != 0) wallAircoefFile << "\n";

                for (int iData = 0; iData < dataNumber; ++ iData)
                {
                    for (int iNode = 0; iNode < coor_count; ++ iNode)
                    {
                        wallAircoefFile << setw(wordWidth) << vector_data[iData][iNode];
                        if ((iNode + 1) % nWordOfLine == 0) 
                        {
                            wallAircoefFile << "\n";
                        }
                    }
                    if (coor_count % nWordOfLine != 0) wallAircoefFile << "\n";

                    DelPointer(vector_data[iData]);
                }
            }
            else
            {
                int TotalNumFaceNodes_Rect = 0;
                PHRead(cdata, &TotalNumFaceNodes_Rect, 1);

                int *face2node = new int[TotalNumFaceNodes_Rect];
                PHRead(cdata, face2node, TotalNumFaceNodes_Rect);

                wallAircoefFile << "N = " << IMxOrNumPts << "\n";
                wallAircoefFile << "E = " << JMxOrNumElements << "\n";
                wallAircoefFile << "f = FEPOINT" << "\n";
                wallAircoefFile << "ET = quadrilateral" << "\n";

                int wordWidth = 20;
                for (int iNode = 0; iNode < IMxOrNumPts; ++ iNode)
                {
                    wallAircoefFile << setw(wordWidth) << x[iNode]
                                    << setw(wordWidth) << y[iNode]
                                    << setw(wordWidth) << z[iNode];
                    for (int iData = 0; iData < dataNumber; ++ iData)
                    {
                        wallAircoefFile << setw(wordWidth) << vector_data[iData][iNode];
                    }
                    wallAircoefFile << "\n";
                }

                for (int iData = 0; iData < dataNumber; ++ iData)
                {
                    DelPointer(vector_data[iData]);
                }

                int count = 0;
                for (int iCell = 0; iCell < JMxOrNumElements; ++ iCell)
                {
                    for (int iNode = 0; iNode < 4; ++ iNode)
                    {
                        wallAircoefFile << face2node[count] << "  ";
                        count ++;
                    }
                    wallAircoefFile << "\n";
                }

                DelPointer(face2node);
            }

            DelPointer(x);
            DelPointer(y);
            DelPointer(z);

            PHRead(cdata, &TecioMission, 1);
        }
    }

    CloseFile(wallAircoefFile);
}

void CFDSolver::DumpHeatFlux(ActionKey *actkey)
{
    //int zoneIDofThisSolver = GetZoneLocalID();
    //if (zoneIDofThisSolver != 0)
    //{
    //    //! Very Important: only the first zone on a processor need to implement the task.
    //    //! other wise, the tasks would be implemented several times.
    //    return;
    //}

    using namespace PHMPI;
    string wall_heatfluxfile = "wall_heatflux.dat";
    GlobalDataBase::GetData("wall_heatfluxfile", &wall_heatfluxfile, PHSTRING, 1);
    if (wall_heatfluxfile == "")
    {
        return;
    }

    actkey->filename = wall_heatfluxfile;

    ios_base::openmode openmode = ios_base::out|ios_base::app;
    actkey->openmode = openmode;

    int nZones = GetNumberofGlobalZones();
    actkey->GetData()->MoveToBegin();

    RDouble *HeatFlux = new RDouble[7];
    HeatFlux[0] = 0.0, HeatFlux[1] = 0.0, HeatFlux[2] = 0.0, HeatFlux[3] = 0.0, HeatFlux[4] = 0.0, HeatFlux[5] = 0.0, HeatFlux[6] = 0.0;
    // [0]:maxHeatFlux, [1]:averageHeatFlux, [2]:totalHeatFlux, [3]:total_Tw, [4]:average_Tw, [5]:total_Pw, [6]:average_Pw
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int zoneProcessorID = GetZoneProcessorID(iZone);
        int currentProcessorID = GetCurrentProcessorID();
        if (currentProcessorID == zoneProcessorID)
        {
            //! If zone i belongs to the current process.
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->ComputeHeatFlux(actkey, HeatFlux);
        }
    }
    RDouble localMaxHeatFlux = HeatFlux[0];
    RDouble localAverageHeatFlux = HeatFlux[2];
    RDouble localAverageTw = HeatFlux[3];
    RDouble localAveragePw = HeatFlux[5];

    RDouble globalMaxHeatFlux = HeatFlux[0];
    RDouble globalAverageHeatFlux = HeatFlux[2];
    RDouble globalAverageTw = HeatFlux[3];
    RDouble globalAveragePw = HeatFlux[5];
    int nProcessorNumber = GetNumberOfProcessor();
    if (nProcessorNumber > 1)
    {
        //synchronization of communication.
        MPI_Barrier(MPI_COMM_WORLD);

        //Communication.
        MPI_Allreduce(&localMaxHeatFlux, &globalMaxHeatFlux, 1, PH_MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&localAverageHeatFlux, &globalAverageHeatFlux, 1, PH_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localAverageTw, &globalAverageTw, 1, PH_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localAveragePw, &globalAveragePw, 1, PH_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    int GlobalTotalWallFace = GlobalDataBase::GetIntParaFromDB("GlobalTotalWallFace");
    HeatFlux[0] = globalMaxHeatFlux;
    HeatFlux[1] = globalAverageHeatFlux / GlobalTotalWallFace;
    HeatFlux[4] = globalAverageTw / GlobalTotalWallFace;
    HeatFlux[6] = globalAveragePw / GlobalTotalWallFace;

    int currentProcessorID = GetCurrentProcessorID();
    int serverTmp = GetServerSepMode();
    if (currentProcessorID == serverTmp)
    {
        WriteHeatFluxASCII(actkey, HeatFlux);    //! Write data to file.
    }

    DelPointer(HeatFlux);
}

void CFDSolver::WriteHeatFluxASCII(ActionKey *actkey, RDouble *HeatFlux)
{
    if (!actkey->IsNeedOpenFile())
    {
        return;
    }

    ParallelOpenFile(actkey);

    std::ostringstream oss;
    fstream &file = *(actkey->file);

    //! To obtain the reference values.
    Param_CFDSolver *parameters = GetControlParameters();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    RDouble referenceTemperature = parameters->GetRefDimensionalTemperature();
    RDouble refVelocity = refDimensionalSonicSpeed * refMachNumber;
    RDouble refSquareVelocity = refVelocity * refVelocity;
    RDouble refEnergy = refDimensionalDensity * refVelocity * refSquareVelocity;
    RDouble refHeatFluxDimension = refEnergy * 0.001;
    HeatFlux[0] = HeatFlux[0] * refHeatFluxDimension;
    HeatFlux[1] = HeatFlux[1] * refHeatFluxDimension;
    HeatFlux[4] = HeatFlux[4] * referenceTemperature;
    HeatFlux[6] = HeatFlux[6] * refDimensionalDensity * refSquareVelocity;

    if (IfFileEmpty(file))
    {
        oss << "Title=\"The Convergence of wall parameters\"" << endl;
        oss << "Variables=" << endl;
        oss << "iter" << endl;
        oss << "Max_Qw(kW/m^2)" << endl;
        oss << "Average_Qw(kW/m^2)" << endl;
        oss << "Average_Tw(K)" << endl;
        oss << "Average_Pw(Pa)" << endl;
    }

    int iterationStep;
    GlobalDataBase::GetData("outnstep", &iterationStep, PHINT, 1);
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    oss << iterationStep << "    ";
    oss << HeatFlux[0] << "    ";
    oss << HeatFlux[1] << "    ";
    oss << HeatFlux[4] << "    ";
    oss << HeatFlux[6] << "    ";
    oss << "\n";

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkey);
}

#ifdef USE_TecplotLib
void CFDSolver::WriteWallAircoef(ActionKey *actkey, vector <DataContainer *> datalist)
{
    ostringstream oss;
    oss << "x y z";

    int systemgridtype = GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");

    int dataNumber;
    dataNumber = 6;
    oss << " " << "cp";
    oss << " " << "cf";
    oss << " " << "yplus";
    oss << " " << "Q_NonDim";
    oss << " " << "Q_Dim(kW/m2)";
    oss << " " << "pw(Pa)";

    if (systemgridtype == STRUCTGRID && nChemical == 1)
    {
        dataNumber = 9;
        oss << " " << "Qs(kW/m2)";
        oss << " " << "Qv(kW/m2)";
        oss << " " << "Qe(kW/m2)";
    }

    string Variable = oss.str();
    const char *Variables = Variable.c_str();

    const char *pltname = actkey->filename.c_str();

    INTEGER4 Debug = 0;
    INTEGER4 VIsDouble = 1;
    INTEGER4 FileType = 0;
    INTEGER4 I;

    I = TECINI112((char*)"PHengLEI Grid Videotex",
                  (char*)Variables,
                  (char*)pltname,
                  (char*)".",
                  &FileType,
                  &Debug,
                  &VIsDouble);

    bool FileIsFull = false;
    uint_t nzones = datalist.size();
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        int gridTypeofLocalZone = PHMPI::GetZoneGridType(iZone);

        DataContainer *cdata = datalist[iZone];
        cdata->MoveToBegin();

        INTEGER4 IMxOrNumPts, JMxOrNumElements, KMxOrNumFaces;
        INTEGER4 ZoneType;

        if (gridTypeofLocalZone == STRUCTGRID)
        {
            ZoneType = TEC_SPACE::ORDERED;
        }
        else
        {
            ZoneType = TEC_SPACE::FEQUADRILATERAL;
        }

        int oriGridIndex = -1;
        PHRead(cdata, &oriGridIndex, 1);

        int GridID;
        PHRead(cdata, &GridID, 1);

        int TecioMission;
        PHRead(cdata, &TecioMission, 1);

        while (true)
        {
            if (TecioMission == 0) break;

            string bcName;
            cdata->ReadString(bcName);

            PHRead(cdata, &IMxOrNumPts, 1);
            PHRead(cdata, &JMxOrNumElements, 1);
            PHRead(cdata, &KMxOrNumFaces, 1);

            int coor_count = 0;
            PHRead(cdata, coor_count);

            RDouble *x = new RDouble[coor_count];
            RDouble *y = new RDouble[coor_count];
            RDouble *z = new RDouble[coor_count];
            PHRead(cdata, x, coor_count);
            PHRead(cdata, y, coor_count);
            PHRead(cdata, z, coor_count);

            INTEGER4 *ValueLocation = new INTEGER4[dataNumber + 3];
            ValueLocation[0] = 1;
            ValueLocation[1] = 1;
            ValueLocation[2] = 1;

            vector < int > data_count;
            data_count.resize(0);
            vector < RDouble *> vector_data;
            vector_data.resize(0);
            for (int iData = 0; iData < dataNumber; ++ iData)
            {
                int idata_count = 0;
                PHRead(cdata, idata_count);
                data_count.push_back(idata_count);

                int iValueLocation = 0;
                PHRead(cdata, iValueLocation);
                ValueLocation[3 + iData] = iValueLocation;

                RDouble *qq = new RDouble [idata_count];
                PHRead(cdata, qq, idata_count);

                vector_data.push_back(qq);
            }

            INTEGER4 ICellMax                   = 0;
            INTEGER4 JCellMax                   = 0;
            INTEGER4 KCellMax                   = 0;
            double   SolutionTime               = 360.0;
            INTEGER4 StrandID                   = 0;
            INTEGER4 ParentZone                 = 0;
            INTEGER4 IsBlock                    = 1;
            INTEGER4 NumFaceConnections         = 0;
            INTEGER4 FaceNeighborMode           = 0;

            INTEGER4 TotalNumFaceNodes_Rect     = 0;

            if (gridTypeofLocalZone == UNSTRUCTGRID)
            {
                PHRead(cdata, &TotalNumFaceNodes_Rect, 1);
            }

            INTEGER4 NumConnBndryFaces_Rect  = 0;
            INTEGER4 TotalNumBndryConns_Rect = 0;
            INTEGER4 SharConn                = 0;

            ostringstream oss1;

            if (oriGridIndex >= 0)
            {
                oss1 << "oriZone_" << oriGridIndex << " ";
            }
            oss1 << "Zone_" << GridID << "_" << bcName;
            string Zonetitle = oss1.str();
            const char * ZoneTitle = Zonetitle.c_str();

            I = TECZNE112((char*)ZoneTitle,
                          &ZoneType,
                          &IMxOrNumPts,
                          &JMxOrNumElements,
                          &KMxOrNumFaces,
                          &ICellMax,
                          &JCellMax,
                          &KCellMax,
                          &SolutionTime,
                          &StrandID,
                          &ParentZone,
                          &IsBlock,
                          &NumFaceConnections,
                          &FaceNeighborMode,
                          &TotalNumFaceNodes_Rect,
                          &NumConnBndryFaces_Rect,
                          &TotalNumBndryConns_Rect,
                          NULL,
                          ValueLocation,
                          NULL,
                          &SharConn);

            INTEGER4 IsDouble = 1;
            I = TECDAT112(&coor_count, x, &IsDouble);
            I = TECDAT112(&coor_count, y, &IsDouble);
            I = TECDAT112(&coor_count, z, &IsDouble);

            for (int iData = 0; iData < dataNumber; ++ iData)
            {
                I = TECDAT112(&data_count[iData], vector_data[iData], &IsDouble);

                DelPointer(vector_data[iData]);
            }

            if (gridTypeofLocalZone == UNSTRUCTGRID)
            {
                int *face2node = new INTEGER4[TotalNumFaceNodes_Rect];
                PHRead(cdata, face2node, TotalNumFaceNodes_Rect);

                I = TECNOD112(face2node);

                DelPointer(face2node);
            }

            if (FileIsFull == false)
            {
                FileIsFull = true;
            }

            DelPointer(x);
            DelPointer(y);
            DelPointer(z);

            PHRead(cdata, &TecioMission, 1);
            DelPointer(ValueLocation);
        }
    }

    I = TECEND112();
    
    if (FileIsFull == false)
    {
        remove(pltname);
    }

    return;
}
#endif

void CFDSolver::UploadInterfaceData()
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
    {
        ActionKey *actkey = new ActionKey();
        FillActionKey(actkey, UPLOAD_INTERFACE_DATA, iLevel);

        int nZones = GetNumberofGlobalZones();

        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int ZoneProcessorID = GetZoneProcessorIDSepMode(iZone);
            int CurrentProcessorID = GetCurrentProcessorID();

            if (CurrentProcessorID == ZoneProcessorID)
            {
                CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
                solver->UploadInterfaceData(actkey);
            }
        }
        FreePointer(actkey);
    }
}

void CFDSolver::UpdateInterfaceData()
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
    {
        ActionKey *actkey = new ActionKey();
        FillActionKey(actkey, UPDATE_INTERFACE_DATA, iLevel);

        //Task_ServerUpdateInterface *task = new Task_ServerUpdateInterface();
        if (IsNeedNonBlockingCommunication(actkey))
        {
            MainTaskNonBlocking(actkey);
        }
        else
        {
            MainTaskBlocking(actkey);
        }
        //delete task;
        FreePointer(actkey);
    }
}

void CFDSolver::DownloadInterfaceData()
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
    {
        ActionKey *actkey = new ActionKey();
        FillActionKey(actkey, DOWNLOAD_INTERFACE_DATA, iLevel);

        int nZones = GetNumberofGlobalZones();
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int zoneProcessorID = GetZoneProcessorIDSepMode(iZone);
            int currentProcessorID = GetCurrentProcessorID();

            if (currentProcessorID == zoneProcessorID)
            {
                CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
                solver->DownloadInterfaceData(actkey);
            }
        }

        FreePointer(actkey);
    }
}

void CFDSolver::UploadInterfaceData(int iLevel)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;
    ActionKey *actkey = new ActionKey();
    FillActionKey(actkey, UPLOAD_INTERFACE_DATA, iLevel);

    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int ZoneProcessorID = GetZoneProcessorIDSepMode(iZone);
        int CurrentProcessorID = GetCurrentProcessorID();

        if (CurrentProcessorID == ZoneProcessorID)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->UploadInterfaceData(actkey);
        }
    }

    FreePointer(actkey);
}

void CFDSolver::UpdateInterfaceData(int iLevel)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;
    ActionKey *actkey = new ActionKey();
    FillActionKey(actkey, UPDATE_INTERFACE_DATA, iLevel);

    //Task_ServerUpdateInterface *task = new Task_ServerUpdateInterface();
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

void CFDSolver::DownloadInterfaceData(int iLevel)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;
    ActionKey *actkey = new ActionKey();
    FillActionKey(actkey, DOWNLOAD_INTERFACE_DATA, iLevel);

    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int zoneProcessorID = GetZoneProcessorIDSepMode(iZone);
        int currentProcessorID = GetCurrentProcessorID();

        if (currentProcessorID == zoneProcessorID)
        {
            CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->DownloadInterfaceData(actkey);
        }
    }

    FreePointer(actkey);
}

void CFDSolver::UploadInterpointData()
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
    {
        ActionKey *actkey = new ActionKey();
        FillActionKey(actkey, UPLOAD_INTERPOINT_DATA, iLevel);

        int nZones = GetNumberofGlobalZones();
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int ZoneProcessorID = GetZoneProcessorIDSepMode(iZone);
            int CurrentProcessorID = GetCurrentProcessorID();

            if (CurrentProcessorID == ZoneProcessorID)
            {
                CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
                solver->UploadInterpointData(actkey);
            }
        }
        FreePointer(actkey);
    }
}

void CFDSolver::UpdateInterpointData()
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
    {
        ActionKey *actkey = new ActionKey();
        FillActionKey(actkey, UPDATE_INTERPOINT_DATA, iLevel);

        TaskServerUpdateInterpoint *task = new TaskServerUpdateInterpoint();
        if (IsNeedNonBlockingCommunication(actkey))
        {
            task->MainTaskNonBlocking(actkey);
        }
        else
        {
            task->MainTaskBlocking(actkey);
        }

        FreePointer(actkey);
        FreePointer(task);
    }
}

void CFDSolver::DownloadInterpointData()
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
    {
        ActionKey *actkey = new ActionKey();
        FillActionKey(actkey, DOWNLOAD_INTERPOINT_DATA, iLevel);

        int nZones = GetNumberofGlobalZones();
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int zoneProcessorID = GetZoneProcessorIDSepMode(iZone);
            int currentProcessorID = GetCurrentProcessorID();

            if (currentProcessorID == zoneProcessorID)
            {
                CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
                solver->DownloadInterpointData(actkey);
            }
        }

        FreePointer(actkey);
    }
}

void CFDSolver::CommunicationInterpointWeight()
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
    {
        ActionKey *actkey = new ActionKey();
        FillActionKey(actkey, DOWNLOAD_INTERPOINT_DATA, iLevel);

        CommunicationInterpointWeight(actkey);
        DownloadInterpointWeight(actkey);

        FreePointer(actkey);
    }
}

void CFDSolver::UploadOversetData()
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
    {
        ActionKey *actkey = new ActionKey();
        FillActionKey(actkey, UPLOAD_OVERSET_DATA, iLevel);

        int nZones = GetNumberofGlobalZones();
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int zoneProcessorID = GetZoneProcessorIDSepMode(iZone);
            int currentProcessorID = GetCurrentProcessorID();

            if (currentProcessorID == zoneProcessorID)
            {
                CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
                solver->UploadOversetData(actkey);
            }
        }

        FreePointer(actkey);
    }
}

void CFDSolver::UpdateOversetData()
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
    {
        ActionKey *actkey = new ActionKey();
        FillActionKey(actkey, UPDATE_OVERSET_DATA, iLevel);

        int nZones = GetNumberofGlobalZones();
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            OversetTopologyManager *oversetTopologyManager = PHSPACE::GetOversetTopologyManager();
            uint_t numberOfOversetNeighbors = oversetTopologyManager->GetNumberOfOversetNeighbors(iZone);

            for (int iNeighbor = 0; iNeighbor < numberOfOversetNeighbors; ++ iNeighbor)
            {
                int oversetNeighborZoneIndex = oversetTopologyManager->GetNeighborOversetZoneIndex(iZone, iNeighbor);
                PH_Interface(actkey, iZone, oversetNeighborZoneIndex);
            }
        }
        FreePointer(actkey);
    }
}

void CFDSolver::DownloadOversetData()
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    for (int level = 0; level < GetNumberOfMultiGrid(); ++ level)
    {
        ActionKey *actkey = new ActionKey();
        FillActionKey(actkey, DOWNLOAD_OVERSET_DATA, level);

        int nZones = GetNumberofGlobalZones();
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int zoneProcessorID = GetZoneProcessorIDSepMode(iZone);
            int currentProcessorID = GetCurrentProcessorID();

            if (currentProcessorID == zoneProcessorID)
            {
                CFDSolver *solver = static_cast <CFDSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
                solver->DownloadOversetData(actkey);
            }
        }
        FreePointer(actkey);
    }
}

void CFDSolver::TecOutAverageFlow(ActionKey *actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    WriteVisualFile(AverageFlow);
}

void CFDSolver::TecOutAverageReynoldsStress(ActionKey *actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    WriteVisualFile(AverageReynoldsStress);
}

//! Get the file name for Restart info dumping, NS-solver by default.
const string CFDSolver::GetFlowTecFileName()
{
    string visualfile = "flow.dat";
    GlobalDataBase::GetData("visualfile", &visualfile, PHSTRING, 1);
    return visualfile;
}

#ifdef USE_GMRESSOLVER
PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void *dummy)
{
       Vec x;

       /*
          Build the solution vector
       */
       PetscCall(KSPBuildSolution(ksp, NULL, &x));

       /*
          Write the solution vector and residual norm to stdout.
           - PetscPrintf() handles output for multiprocessor jobs
             by printing from only one processor in the communicator.
           - The parallel viewer PETSC_VIEWER_STDOUT_WORLD handles
             data from multiple processors so that the output
             is not jumbled.
       */
       //PetscCall(PetscPrintf(PETSC_COMM_WORLD, "iteration %" PetscInt_FMT " solution vector:\n", n));
       //PetscCall(VecView(x, PETSC_VIEWER_STDOUT_WORLD));
       PetscCall(PetscPrintf(PETSC_COMM_WORLD, "iteration %" PetscInt_FMT " KSP Residual norm %14.12e \n", n, (double)rnorm));
       return 0;
}
#endif

}

