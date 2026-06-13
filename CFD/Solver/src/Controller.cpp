#include "Controller.h"
#include "CFDSolver.h"
#include "Region.h"
#include "Math_BasisFunction.h"
#include "Task.h"
#include "TK_Log.h"

#ifdef USE_CUDA
#include "BasicDeviceVariables.h"
#ifdef CUDAUNITTEST
#include "OutputDebug.h"
#endif
#ifdef NEWMPI
#include "GPUMPICommunication.h"
#endif 
#endif

using namespace std;

namespace PHSPACE
{
Controller::Controller(int iSolver, Region *region)
{
    this->solverID = iSolver;
    this->region = region;

    qProxy = 0;

    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();

    rhsProxy = new ZoneFieldProxy(iSolver, GetNumberOfMultiGrid(), nZones);

    if (GetNumberOfMultiGrid() > 1)
    {
        qProxy = new ZoneFieldProxy(iSolver, GetNumberOfMultiGrid(), nZones);
    }
}

Controller::~Controller()
{
    if (GetNumberOfMultiGrid() > 1)
    {
        delete qProxy;
    }

    delete rhsProxy;
}

LIB_EXPORT void Controller::InitCoarseGridsFlow(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        Grid *coarseGrid = grid->GetCoarseGrid();
        if (coarseGrid)
        {
            if (iSolver != NS_SOLVER) continue;
            CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
            //solver->InitCGrid(grid, coarseGrid);
            solver->ComputeGamaAndTemperature(coarseGrid);
            solver->Boundary(coarseGrid);
            solver->ComputeGamaAndTemperature(coarseGrid);
            solver->ComputeNodeValue(coarseGrid);
        }
    }
}

LIB_EXPORT void Controller::Initialize(int iSolver)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int nProtectData = 0;
    if (GlobalDataBase::IsExist("nProtectData", PHINT, 1))
    {
        nProtectData = GlobalDataBase::GetIntParaFromDB("nProtectData");
    }

    int level = 0;
    //! First: Initialize memory.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        PHSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->InitMemory();
    }

    //! Second: Initialize flow variables on the finest mesh.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        PHSolver *solver = GetCFDSolver(zoneID, iSolver);

    #ifdef USE_CUDA
        using namespace GPUMemory;
        using namespace GPUGeomVariables;
        using namespace GPUFlowVariables;
        if (iSolver == 0)    //! NS solver
        {
            GPUFlowVarMemAlloc(iSolver);    //! for establish TransferBuffer, iSolver is required
        }
        else    //! Turb solver
        {
            GPUFlowVarTurbMemAlloc(iSolver);
        }
    #endif

        if (nProtectData == 0)
        {
            solver->InitFlow();
        }
        else
        {
            solver->InitProtectedFlow();
        }
    }

    //! Third: Initialize flow variables on coarser grids, if multi-grid is used.
    int nMGLevel = PHSPACE::GlobalDataBase::GetIntParaFromDB("nMGLevel");
    for (int iLevel = 0; iLevel < nMGLevel - 1; ++ iLevel)
    {
        RestrictAllQ(iSolver, iLevel);
        InitCoarseGridsFlow(iSolver, iLevel);
    }

    //! Finally: Interfaces communications.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        PHSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->CommunicationInterfaceData();
        InterpointInformation *interPointInfor = grid->GetInterpointInfo();
        if (interPointInfor)
        {
            solver->CommunicationInterpointData();
        }
    }

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone) continue;

        for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); iLevel ++)
        {
            Grid *grid = GetGrid(zoneID, iLevel);
            CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
            solver->ModifyNodeValue(grid);
            solver->ModifyQTurbNodeValue(grid);
        }
    }

    int compressible = GlobalDataBase::GetIntParaFromDB("compressible");
    if (COMPRESSIBLE == compressible)
    {
        ComputeNodeBCType();
    }

    //! Initialization of boundary and dependent variables for first iteration.
    //! Turn on it when gradient pre-computing for turbulent solver is used.
    //SetBoundaryValueOfField(iSolver, level);
    //UpdateIndependentField(iSolver, level);
}

void Controller::CleanUp(int iSolver)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    //! First: Initialize memory.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        PHSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->ReleaseMemory();
    }
}

void Controller::ComputeNodeBCType()
{
    int sysGridType = GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (sysGridType != UNSTRUCTGRID)
    {
        return;
    }

    using namespace PHMPI;

    int needCommunication = 1;
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int iZoneProcessorID = GetZoneProcessorID(iZone);
        int myid = GetCurrentProcessorID();
        if (myid == iZoneProcessorID)
        {
            //! Compute node bctype in each zone.
            Grid *grid = GetGrid(iZone, 0);
            grid->ComputeNodeBCType();
            if (!(grid->GetInterpointInfo()))
            {
                needCommunication = 0;
            }
        }
    }

    PH_CompareMaxMin(needCommunication, 2);
    if (!needCommunication)
    {
        return;
    }
    CommunicationNodeBCType();
}

void Controller::CommunicationNodeBCType()
{
    using namespace PHMPI;
    using namespace PHENGLEI;

    map<int, int> bcTypeMap;
    bcTypeMap.insert(pair<int,int>(SOLID_SURFACE   , 4));
    bcTypeMap.insert(pair<int,int>(FARFIELD        , 3));
    bcTypeMap.insert(pair<int,int>(SYMMETRY        , 2));
    bcTypeMap.insert(pair<int,int>(INTERFACE       , 1));
    bcTypeMap.insert(pair<int,int>(INFLOW          , 1));
    bcTypeMap.insert(pair<int,int>(OUTFLOW         , 1));
    bcTypeMap.insert(pair<int,int>(PRESSURE_INLET  , 1));
    bcTypeMap.insert(pair<int,int>(PRESSURE_OUTLET , 1));
    bcTypeMap.insert(pair<int,int>(OUTFLOW_CONFINED, 1));

    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int myid = GetCurrentProcessorID();
        int tag = iZone;
        int sendProc = GetZoneProcessorID(iZone);

        ZoneNeighbor *zoneNeighborForPoint = zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);
        std::size_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();

        for (std::size_t iNeighborZone = 0; iNeighborZone < numberOfNeighbor; ++ iNeighborZone)
        {
            int neighborZoneIndex = zoneNeighborForPoint->GetZoneIndexOfNeighbor(iNeighborZone);
            int recv_proc = GetZoneProcessorID(neighborZoneIndex);

            DataContainer *cdata = new DataContainer();

            if (myid == sendProc)
            {
                UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(iZone, 0));
                InterpointInformation *interpointInformation = grid->GetInterpointInfo();
                int *nodeBCType = reinterpret_cast <int *> (grid->GetDataPtr("nodeBCType"));

                cdata->MoveToBegin();

                int iNeighbor = interpointInformation->FindIthNeighbor(neighborZoneIndex);
                int numberOfInterpointsForNeighbor = interpointInformation->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForSend = interpointInformation->GetPointIndexForSend(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                int BCType;
                int *qtmp = new int [numberOfInterpointsForNeighbor];
                for (int ipointLocal = 0; ipointLocal < numberOfInterpointsForNeighbor; ++ ipointLocal)
                {
                    iPoint = pointIndexForSend[ipointLocal];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    BCType = nodeBCType[globalPoint];

                    qtmp[ipointLocal] = BCType;
                }
                cdata->Write(qtmp, numberOfInterpointsForNeighbor * sizeof(int));
                delete [] qtmp;    qtmp = nullptr;
            }

            PH_Trade(cdata, sendProc, recv_proc, tag);

            if (myid == recv_proc)
            {
                UnstructGrid *gridNeighbor = UnstructGridCast(PHSPACE::GetGrid(neighborZoneIndex, 0));
                InterpointInformation *interpointInformationNeighbor = gridNeighbor->GetInterpointInfo();
                int *nodeBCType = reinterpret_cast<int *>(gridNeighbor->GetDataPtr("nodeBCType"));

                cdata->MoveToBegin();

                int iNeighbor = interpointInformationNeighbor->FindIthNeighbor(iZone);
                int numberOfInterpointsForNeighbor = interpointInformationNeighbor->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForReceive = interpointInformationNeighbor->GetPointIndexForReceive(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformationNeighbor->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                int recvBCType;
                for (int ipointLocal = 0; ipointLocal < numberOfInterpointsForNeighbor; ++ ipointLocal)
                {
                    iPoint = pointIndexForReceive[ipointLocal];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    cdata->Read(&recvBCType, sizeof(int));

                    int localBCIndex = bcTypeMap[nodeBCType[globalPoint]];
                    int recvBCIndex  = bcTypeMap[recvBCType];
                    if (recvBCIndex > localBCIndex)
                    {
                        nodeBCType[globalPoint] = recvBCType;
                    }

                }
            }

            delete cdata;    cdata = nullptr;
        }
    }
}

LIB_EXPORT void Controller::SolveSteadyField()
{
    int iSolver = this->GetSolverIndex();
    int compressible = GlobalDataBase::GetIntParaFromDB("compressible");
    if (COMPRESSIBLE == compressible)
    {
        //! Zero the residuals of the finest grid.
        ZeroResiduals(iSolver);

        MultiGrid(iSolver);
    }
    else if (INCOMPRESSIBLE == compressible)
    {
        SolveIncomSteadyField(iSolver);
    }
}

LIB_EXPORT void Controller::SolveIncomSteadyField(int iSolver, int level)
{
    int nLocalZones  = PHMPI::GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForPressureBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        CFDSolver* solver = GetCFDSolver(zoneID, iSolver);

        solver->SolveIncomSteadyField();
    }
}

LIB_EXPORT void Controller::SolveHyrbridDBPBSteadyField(int iSolver, int level)
{
    ZeroResiduals(iSolver);
    MultiGrid(iSolver);
    SolveIncomSteadyField(iSolver);
}

LIB_EXPORT void Controller::SolveUnsteadyField()
{
    int iSolver = this->GetSolverIndex();
    int inerStepInterval = GlobalDataBase::GetIntParaFromDB("innstep");

    while (true)
    {
        inerStepInterval += 1;
        GlobalDataBase::UpdateData("innstep", &inerStepInterval, PHINT, 1);

        SolveSteadyField();

        if (CVGNorm(iSolver))
        {
            UpdateUnsteadyFlow(iSolver);
            break;
        }

        PHSolver *solver = GetCFDSolver(0, iSolver);
        solver->CommunicationInterpointData();
    }
}

void Controller::MultiGrid(int iSolver, int level)
{
    CreateRHSProxy(iSolver, level);
    StoreRhsByResidual(iSolver, level);
    int outerIteration = GlobalDataBase::GetIntParaFromDB("outnstep");
    int MGPreIteration = GlobalDataBase::GetIntParaFromDB("MGPreIteration");
    bool computeOnFinestGrid = (outerIteration % MGPreIteration != 0);

    //! If grid is in the coarsest level.
    if (ISCoarsestGrid(level) || iSolver != NS_SOLVER || computeOnFinestGrid)
    {
        Relaxation(iSolver, level);
    }
    else
    {
        int cLevel = level + 1;

        //! Smooth the solutions before corrections.
        Relaxation(iSolver, level);

        CommunicationInterfaceDataOneLevel(iSolver, level);

        //! Restrict from fine to coarse grid for all q.
        //! Corresponding to nsmb: wsav = restr(w).
        //! Call the following function to update the Q value(cq) on the coarse grid.
        RestrictAllQ(iSolver, level);

        CommunicationInterfaceDataOneLevel(iSolver, cLevel);

        //! Store coarse grid solution in temporary array cqsav.
        //! cqsav represents the variables on the coarse grid.
        CreateQProxy(iSolver, cLevel);

        //! After LoadQ(cgrid, cqsav), cqsav is actually wsav(Corresponding to NSMB).
        //! LoadQ  gets out the value of q(cq) on coarse grid and assigns it to cqsav.
        LoadQ(iSolver, cLevel);

        //! Defect.
        //! The followings corresponding to nsmb d = RL-1(wsav) - restr(RL(w)-f).
        //! in which (*res) = - (*rhs) correspond to - f
        //! Note that the value of res is on fine grid, rhs ought be f��and this should be checked.
        //! Here should be ni and nj,rather than cni and cnj!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        InitResidual(iSolver, level);

        //! After UpdateResiduals,(*res) = RL(w) - f.
        UpdateResiduals(iSolver, level);

        //! After RestrictDefect,(*cres) = - restr(RL(w) - f).
        RestrictDefect(iSolver, level);

        //! After UpdateResiduals,(*cres) = RL-1(wsav) - restr(RL(w) - f).
        UpdateResiduals(iSolver, cLevel);

        //! Solve on coarse grid.
        int MGFasType = GlobalDataBase::GetIntParaFromDB("MGFasType");
        int ifLowSpeedPrecon = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");
        if (ifLowSpeedPrecon == 1)
        {
            MGFasType = 1;    //! Precondition method is only available for V cycle.
        }

        //! After MultiGrid(cgrid), the value of q on coarse grid has been updated.
        //! At this time,the value of q on coarse grid equals to w0 of NSMB.
        //! MGFasType = 1 V cycle
        //! MGFasType = 2 W cycle
        for (int i = 0; i < MGFasType; ++ i)
        {
            MultiGrid(iSolver, cLevel);
        }

        CommunicationInterfaceDataOneLevel(iSolver, cLevel);

        //! Put correction back to the coarse grid.
        //! Actually, PutCorrection has done the operation of (w0 - wsav).
        PutCorrection(iSolver, cLevel);

        //! Communicate coarse grid here!!!

        //! Correct variables in fine grid.
        //! Actually CorrectFineGrid has done the operation of w = w + prol(w0 - wsav).
        CorrectFineGrid(iSolver, level);

        //! Smooth the solutions after corrections.
        int n_post = GlobalDataBase::GetIntParaFromDB("n_post");
        for (int iter = 0; iter < n_post; ++ iter)
        {
            Relaxation(iSolver, level);
        }

        //! The following function actually recover the value of q on coarse grid.
        PutCorrectionBack(iSolver, cLevel);
        //! Delete the temporary variables of QProxy.
        DestroyQProxy(cLevel);
    }
    //! Put RHS back into its places.
    //! Here, we just recover the initial value of res.
    RecoverResidual(iSolver, level);

    //! Delete the temporary variables of rhsProxy.
    DestroyRHSProxy(level);
}

void Controller::UpdateUnsteadyFlow(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        solver->UpdateUnsteadyFlow(grid);
    }
}

void Controller::Relaxation(int iSolver, int level)
{
    using namespace PHSPACE;
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FD_METHOD)
    {
        PreSolveforStructHighOrder(iSolver, level);
    }

    int isUseLocalCFL = 0;
    if (GlobalDataBase::IsExist("isUseLocalCFL", PHINT, 1))
    {
        isUseLocalCFL = GlobalDataBase::GetIntParaFromDB("isUseLocalCFL");
    }

    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    int ifLowSpeedPrecon = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");

    if (tscheme == LU_SGS || level != 0)
    {
        //! LU-SGS implicit time integration method.
        if (!ifLowSpeedPrecon)
        {
            LoadResiduals(iSolver, level);

            UpdateResiduals(iSolver, level);

            SpectrumRadius(iSolver, level);

            TimeStep(iSolver, level);

            if (isUseLocalCFL == 1)
            {
                DualTimeStepsLUSGS(iSolver, level);
            }
            else
            {
                LUSGS(iSolver, level);
            }
        }
        else
        {
            SpectrumRadius(iSolver, level);

            TimeStep(iSolver, level);

            LoadResiduals(iSolver, level);

            UpdateResiduals(iSolver, level);

            if (isUseLocalCFL == 1)
            {
                DualTimeStepsLUSGS(iSolver, level);
            }
            else
            {
                LUSGS(iSolver, level);
            }
        }
    }
    else if (tscheme == GMRES)
    {
        //! Generalized Minimal Residual method added by Sen Zhang, Boqian Wang, et al. 
        LoadResiduals(iSolver, level);

        UpdateResiduals(iSolver, level);

        SpectrumRadius(iSolver, level);

        TimeStep(iSolver, level);

#ifdef USE_GMRESSOLVER
        GMRESSolver(iSolver, level);
#endif

    }
    else if (tscheme == MULTI_STAGE)
    {
        //! Explicit multi-stage time marching method.
        SetBoundaryValueOfField(iSolver, level);

        //! To obtain the viscous coff in spectrum radium term
        ComputeViscousCoeff(iSolver, level);

        SpectrumRadius(iSolver, level);

        TimeStep(iSolver, level);

        RungeKutta(iSolver, level);
    }
    else if (tscheme == Matrix_LU_SGS)
    {
        LoadResiduals(iSolver, level);

        UpdateResiduals(iSolver, level);

        SpectrumRadius(iSolver, level);

        TimeStep(iSolver, level);

        MatrixLUSGS(iSolver, level);
    }
    else
    {
        //! default
        LoadResiduals(iSolver, level);

        UpdateResiduals(iSolver, level);

        SpectrumRadius(iSolver, level);

        TimeStep(iSolver, level);

        if (isUseLocalCFL == 1)
        {
            DualTimeStepsLUSGS(iSolver, level);
        }
        else
        {
            LUSGS(iSolver, level);
        }
    }
}

void Controller::TimeStep(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    if (iSolver != NS_SOLVER)
    {
        //! Only get time step for NS solver.
        return;
    }

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->TimeStep(grid);
    }

    ReduceMaxTimeStep(iSolver, level);
}

RDouble Controller::ComputeMinTimeStep(int iSolver, int level)
{
    using namespace PHMPI;
    int outnstep      = GlobalDataBase::GetIntParaFromDB("outnstep");
    int allReduceStep = GlobalDataBase::GetIntParaFromDB("allReduceStep");
    RDouble globalMinTimeStep = GlobalDataBase::GetDoubleParaFromDB("globalMinTimeStep");

    if (outnstep % allReduceStep != 0 && globalMinTimeStep > 0.0)
    {
        //! Compute global minimum time step each 'allReduceStep' iterations.
        return globalMinTimeStep;
    }

    int nLocalZones = GetNumberofLocalZones();
    RDouble globalMinDt =  LARGE;
    RDouble localMaxDt  = -LARGE;

    //! Important: Can not use OMP directly, due to the MIN/MAX comparison!
    // #ifdef USE_OMP
    // #pragma omp parallel for
    // #endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        RDouble minDt, maxDt;
        solver->ComputeMinTimeStep(grid, minDt, maxDt);

        globalMinDt = MIN(globalMinDt, minDt);
        localMaxDt  = MAX(localMaxDt, maxDt);
    }


    PH_CompareMaxMin(globalMinDt, 2);
    GlobalDataBase::UpdateData("globalMinTimeStep", &globalMinDt, PHDOUBLE, 1);

    return globalMinDt;
}

void Controller::ReduceMaxTimeStep(int iSolver, int level)
{
    using namespace PHMPI;
    int ifLocalTimeStep = GlobalDataBase::GetIntParaFromDB("ifLocalTimeStep");
    if (ifLocalTimeStep != 0)
    {
        //! NON local time step method, return.
        return;
    }
    if (iSolver != NS_SOLVER)
    {
        //! Only reduce time step for NS solver.
        return;
    }

    RDouble globalMinDt = ComputeMinTimeStep(iSolver, level);

    int nLocalZones = GetNumberofLocalZones();
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->ReduceMaxTimeStep(grid, globalMinDt);
    }
}

void Controller::Diagonal(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->Diagonal(grid);
    }
}

void Controller::SpectrumRadius(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    if (iSolver != NS_SOLVER)
    {
        return;
    }

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->SpectrumRadius(grid);
    }
}

void Controller::ZeroResiduals(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->ZeroResiduals(grid);
    }
}

void Controller::MixingPlaneBoundary(int iSolver,int level)
{
    using namespace PHMPI;

    int nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");

    //! MixingPlane Boundary.
    if (referenceFrame == ROTATIONAL_FRAME)
    {
        int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");

        if (nTurboZone > 1)
        {
            //! first compute flux of all MixingPlane.
            for (int iZone = 0; iZone < nTurboZone; ++iZone)
            {
                int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

                bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
                if (!ifCalTheZone)
                {
                    continue;
                }

                Grid *grid = GetGrid(zoneID, level);
                CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

                solver->AverageMixingPlane(grid);
            }

            //! second communicate data.
            for (int iZone = 0; iZone < nTurboZone; ++iZone)
            {
                int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

                bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
                if (!ifCalTheZone)
                {
                    continue;
                }

                Grid *grid = GetGrid(zoneID, level);
                CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

                int upStreamZone, downStreamZone;

                if (zoneID == 0)
                {
                    downStreamZone = zoneID + 1;
                    Grid *NeighborGrid = GetGrid(downStreamZone, level);
                    solver->MixingPlaneDataTransfer(grid, NeighborGrid);
                }
                else if (zoneID == nTurboZone - 1)
                {
                    upStreamZone = zoneID - 1;
                    Grid *NeighborGrid = GetGrid(upStreamZone, level);
                    solver->MixingPlaneDataTransfer(grid, NeighborGrid);
                }
                else
                {
                    upStreamZone = zoneID - 1;
                    downStreamZone = zoneID + 1;
                    Grid *NeighborGrid1 = GetGrid(upStreamZone, level);
                    Grid *NeighborGrid2 = GetGrid(downStreamZone, level);
                    solver->MixingPlaneDataTransfer(grid, NeighborGrid1);
                    solver->MixingPlaneDataTransfer(grid, NeighborGrid2);
                }
            }

            //! third NonReflective
            for (int iZone = 0; iZone < nTurboZone; ++iZone)
            {
                int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

                bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
                if (!ifCalTheZone)
                {
                    continue;
                }

                Grid *grid = GetGrid(zoneID, level);
                CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

                int upStreamZone, downStreamZone;

                if (zoneID == 0)
                {
                    downStreamZone = zoneID + 1;
                    Grid *NeighborGrid = GetGrid(downStreamZone, level);
                    solver->NonReflective(grid, NeighborGrid);
                }
                else if (zoneID == nTurboZone - 1)
                {
                    upStreamZone = zoneID - 1;
                    Grid *NeighborGrid = GetGrid(upStreamZone, level);
                    solver->NonReflective(grid, NeighborGrid);
                }
                else
                {
                    upStreamZone = zoneID - 1;
                    downStreamZone = zoneID + 1;
                    Grid *NeighborGrid1 = GetGrid(upStreamZone, level);
                    Grid *NeighborGrid2 = GetGrid(downStreamZone, level);
                    solver->NonReflective(grid, NeighborGrid1);
                    solver->NonReflective(grid, NeighborGrid2);
                }
            }

            //! Fourth SetMixingPlaneData
            for (int iZone = 0; iZone < nTurboZone; ++iZone)
            {
                int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

                bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
                if (!ifCalTheZone)
                {
                    continue;
                }

                Grid *grid = GetGrid(zoneID, level);
                CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

                solver->SetMixingPlaneData(grid);
            }
        }
    }
}


void Controller::UpdateResiduals(int iSolver, int level)
{
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    string highordersolvername = GlobalDataBase::GetStrParaFromDB("str_highorder_solver");

    if ((systemGridType == MIXGRID && iSolver == NS_SOLVER && isFVMOrFDM == FV_METHOD) || (iSolver == NS_SOLVER && isFVMOrFDM == FD_METHOD && highordersolvername.substr(0, 4) == "HDCS"))
    {
        //! ugly codes, modify later.
        //! for mixgrid, more variables need communication for conservation.
        //! InvscidFlux() and ViscousFlux() must be breaked into a lot of small functions.
        UpdateResidualsOnlyforMixGrid(iSolver, level);
        return;
    }

    using namespace PHMPI;

#ifdef USE_GMRESSOLVER
    // GMRESCoupled
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if( viscousType == ONE_EQU && iSolver == NS_SOLVER )
    {
        SetBoundaryValueOfField(NS_SOLVER, level);
        SetBoundaryValueOfField(TURB_SOLVER, level);

        UpdateIndependentField(NS_SOLVER, level);
        UpdateIndependentField(TURB_SOLVER, level);

        SetBoundaryValueOfField(TURB_SOLVER, level);

    }
    else if ( viscousType == INVISCID || viscousType == LAMINAR || viscousType == TWO_EQU )
    {
        SetBoundaryValueOfField(iSolver, level);

        UpdateIndependentField(iSolver, level);
    }

#else 
    SetBoundaryValueOfField(iSolver, level);

    UpdateIndependentField(iSolver, level);
#endif

    int nLocalZones = GetNumberofLocalZones();

    //!Overset communications for independent variables such sa gradients.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        PHSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->CommunicationOversetSlipData();
    }

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->UpdateResiduals(grid);
    }
}

void Controller::UpdateResidualsOnlyforMixGrid(int iSolver, int level)
{
    using namespace PHMPI;

    SetBoundaryValueOfField(iSolver, level);

    UpdateIndependentField(iSolver, level);

    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->UpdateQlQrOnlyforMixGrid(grid);
    }
    CommunicateQlQrOnlyforMixGrid();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->UpdateInviscidfluxOnlyforMixGrid(grid);

        solver->UpdateViscousfluxOnlyforMixGrid(grid);
    }
    CommunicateFacefluxOnlyforMixGrid();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }
        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
        solver->LoadGlobalfacefluxtoResOnlyforMixGrid(grid);
    }

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->SourceFlux(grid);
    }
}


void Controller::LoadResiduals(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    FieldProxy **rhsProxy = GetRHSProxy(level);

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        CFDSolver *solver = GetCFDSolver(zoneID,iSolver);

        solver->LoadResiduals(grid, rhsProxy[zoneID]);
    }
}

void Controller::PressureFactor(int iSolver, int level)
{
    int nEntropyFixMethod = GlobalDataBase::GetIntParaFromDB("roeEntropyFixMethod");

    if (NS_SOLVER != iSolver || (nEntropyFixMethod != 2 && nEntropyFixMethod != 6))
    {
        return;
    }

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->CalPressureFactor(grid);
    }

    CommunicateSpecificArray(level, "rtem", 0);
}

void Controller::ComputeGradient(int iSolver, int level)
{
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (viscousType == INVISCID && STRUCTGRID == systemGridType)
    {
        return;
    }

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        int gridTypeofLocalZone = GetZoneGridType(zoneID);

        if (gridTypeofLocalZone == STRUCTGRID)
        {
            string sgsmodel;
            GlobalDataBase::GetData("sgsmodel", &sgsmodel, PHSTRING, 1);
            solver->ComputeGradientCellCenter(grid);

            if (MIXGRID == systemGridType && isFVMOrFDM == FV_METHOD)
            {
                solver->ComputeGradientCellCenterOfruvwptOnlyForMixGrid(grid);
            }
        }
        else if (gridTypeofLocalZone == UNSTRUCTGRID)
        {
            solver->ComputeGradient(grid);
            if (MIXGRID == systemGridType && isFVMOrFDM == FV_METHOD)
            {
                solver->ComputelimiterOnlyForMixGrid(grid);
            }
        }
    }

    CommunicateGradientArray(iSolver, level);
}

void Controller::ComputeLimiter(int iSolver, int level)
{
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (UNSTRUCTGRID != systemGridType || NS_SOLVER != iSolver)
    {
        return;
    }

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
        solver->ComputeLimiter(grid);
    }

    int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    CommunicateSpecificArray(level, "limit2D", nLaminar);
}

void Controller::ComputeGamaAndTField(int iSolver, int level)
{
    if (NS_SOLVER != iSolver)
    {
        return;
    }

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->ComputeGamaAndTemperature(grid);
    }
}

void Controller::ComputeViscousCoeff(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if (viscousType == INVISCID)
    {
        return;
    }

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->ComputeViscousCoeff(grid);
    }
}

void Controller::ComputeBlendFunctionSST(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (iSolver == NS_SOLVER || iSolver == TRANSITION_SOLVER || viscousType != TWO_EQU || 0 != level || isFVMOrFDM == FD_METHOD)
    {
        return;
    }

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->Crossing(grid);
        solver->Blending(grid);
    }
    CommunicateSpecificArray(level, "blend", 0);
}

//! Compute the precondition coefficient of low mach flow.
void Controller::ComputePreconditionCoefficient(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    int ifLowSpeedPrecon =  GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (iSolver != NS_SOLVER  || 0 != level || isFVMOrFDM == FD_METHOD || ifLowSpeedPrecon != 1)
    {
        return;
    }
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        if (isUnsteady == 0)
        {
            solver->ComputePreconditionCoefficient(grid);
        }
        else
        {
            solver->ComputePreconditionCoefficientUnsteady(grid);
        }

        //solve M * T-, when sign = -1;
        solver->ConservativePreconMatrix(grid, -1);
    }

    CommunicateSpecificArray(level, "preconCoefficient", 0);

    if (isUnsteady)
    {
        CommunicateSpecificArray(level, "timeCoefficient", 0);
        CommunicateSpecificArray(level, "timeCoefficientInverse", 0);
    }
}

void Controller::SetBoundaryValueOfField(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->Boundary(grid);
    }
    MixingPlaneBoundary(iSolver, level);
}

void Controller::PreSolveforStructHighOrder(int iSolver, int level)
{
    //int mission_Comm(q_FD                                    ) =  5;
    //int mission_Comm(rtem                                    ) =  6;
    //int mission_Comm(rtem, dqdxyznokxi,dqdxyznoeta,dqdxyznocta) =  7;
    //int mission_Comm(rtem, gradUVWTCellCenterXYZ             ) =  8;
    //int mission_Comm(vist, gradqturb                         ) =  9;
    //int mission_Comm(vist, gradqturb, blending               ) = 10;

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");

    if (iSolver == NS_SOLVER)
    {
        GetIndependentVariablesforStructHighOrder(iSolver, level);
        this->region->CommGridInfoOnlyforStructHighOrder(5);

        GetDependentVariablesforStructHighOrder(iSolver, level);
        if (viscousType == INVISCID)
        {
            //! Euler, no gradient communicated.
            this->region->CommGridInfoOnlyforStructHighOrder(6);
        }
        else
        {
            string gradient_method = GlobalDataBase::GetStrParaFromDB("structhighordergradient");
            if (gradient_method.substr(0,10) == "chain_rule")
            {
                this->region->CommGridInfoOnlyforStructHighOrder(7);
            }
            else
            {
                this->region->CommGridInfoOnlyforStructHighOrder(8);
            }
        }
    }
    else if (iSolver == TURB_SOLVER)
    {
        GetDependentVariablesforStructHighOrder(iSolver, level);
        if (viscousType == ONE_EQU)
        {
            this->region->CommGridInfoOnlyforStructHighOrder(9);
        }
        else if (viscousType == TWO_EQU)
        {
            this->region->CommGridInfoOnlyforStructHighOrder(10);
        }
        if (iLES == LES_SOLVER)
        {
            this->region->CommGridInfoOnlyforStructHighOrder(8);

            this->region->CommGridInfoOnlyforStructHighOrder(11);
        }
    }
}

void Controller::GetIndependentVariablesforStructHighOrder(int iSolver, int level)
{    
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->GetIndependentVariablesforStructHighOrder(grid);
    }
}

void Controller::GetDependentVariablesforStructHighOrder(int iSolver, int level)
{    
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->GetDependentVariablesforStructHighOrder(grid);
    }   
}

void Controller::RungeKuttaThirdTVD(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    FieldProxy **q0Proxy = new FieldProxy *[nLocalZones];

    CreatqProxyfor3rdRK(iSolver, level, q0Proxy);

    Stagefor3rdRK(iSolver, level, q0Proxy, 0);

    Stagefor3rdRK(iSolver, level, q0Proxy, 1);

    Stagefor3rdRK(iSolver, level, q0Proxy, 2);
    
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        delete q0Proxy[iZone];
    }
    delete [] q0Proxy;
}

void Controller::CreatqProxyfor3rdRK(int iSolver, int level, FieldProxy **q0Proxy)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
        
        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        q0Proxy[iZone] = solver->CreateFieldProxy(grid);
        RDouble4D &q0 = q0Proxy[iZone]->GetField_STR();
        q0 = q;
    }
}

void Controller::Stagefor3rdRK(int iSolver, int level, FieldProxy **q0Proxy, int nstage)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();    
    
    if (nstage != 0)
    {
        //! already solved before timestep;
        PreSolveforStructHighOrder(iSolver, level);
    }

    LoadResiduals(iSolver, level);
    UpdateResiduals(iSolver, level);      

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);

        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
        FieldProxy *resRKProxy = solver->CreateFieldProxy(grid);

        if (nstage == 0)
        {
            solver->RungeKuttaResidual(grid, resRKProxy, 1.0);
            solver->UpdateFlowField(grid, q0Proxy[iZone], resRKProxy);            
        }
        else if (nstage == 1)
        {
            FieldProxy *qtempProxy = solver->CreateFieldProxy(grid);
            
            RDouble4D &qtemp   = qtempProxy->GetField_STR();
            RDouble4D &q0   = q0Proxy[iZone]->GetField_STR();
            RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

            qtemp = 0.75 * q0 + 0.25 * q;
            solver->RungeKuttaResidual(grid, resRKProxy, 0.25);
            solver->UpdateFlowField(grid, qtempProxy, resRKProxy);

            delete qtempProxy;
        }
        else if (nstage == 2)
        {
            FieldProxy *qtempProxy = solver->CreateFieldProxy(grid);
            
            RDouble4D &qtemp   = qtempProxy->GetField_STR();
            RDouble4D &q0   = q0Proxy[iZone]->GetField_STR();
            RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

            qtemp = (q0 + 2.0 * q) / 3.0;
            solver->RungeKuttaResidual(grid, resRKProxy, 2.0 / 3.0);
            solver->UpdateFlowField(grid, qtempProxy, resRKProxy);

            delete qtempProxy;
        }
        delete resRKProxy;
    }
}

void Controller::UpdateIndependentField(int iSolver, int level)
{
    this->ComputeGamaAndTField(iSolver, level);

    this->PressureFactor(iSolver, level);

    this->ComputeViscousCoeff(iSolver, level);

    this->ComputeGradient(iSolver, level);

    this->ComputeLimiter(iSolver, level);

    this->ComputeBlendFunctionSST(iSolver, level);

    this->ComputePreconditionCoefficient(iSolver, level);
}

void Controller::LUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy **dqProxy, FieldProxy **LUplusDQ)
{
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");

    int nLocalZones = PHMPI::GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        //! Memory allocating.
        dqProxy[iZone] = solver->CreateFieldProxy(grid);
        LUplusDQ[iZone]= solver->CreateFieldProxy(grid);

        int gridtype = grid->Type();
        if (isFVMOrFDM == FD_METHOD && gridtype == STRUCTGRID)
        {
            solver->LUSGSInitializationStructHighOrder(grid, dqProxy[iZone], LUplusDQ[iZone]);
        }
        else
        {
            //! The delta Flux is initialized by zero, which will be used in the forward sweep.
            solver->FillField(grid, LUplusDQ[iZone], zero);

            //! Warning: it is important to research the effect of the initialization in future!
            //! e.g., using rhs to init when nsweep > 1.
            //! Initialization.
            if (nSweep == 1)
            {
                //! Single sweep, init by res.
                FieldProxy *rhsProxy = solver->GetFieldProxy(grid, solver->CastName("res")); 
                solver->FillField(grid, dqProxy[iZone], rhsProxy);
                delete rhsProxy;
            }
            else
            {
                //! Multi-sweep, init by zero.
                solver->FillField(grid, dqProxy[iZone], zero);
            }
        }
    }
}

void Controller::BLUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy **dqProxy, FieldProxy **LUplusDQ)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        //! Memory allocating.
        dqProxy[iZone] = solver->CreateFieldProxy(grid);
        LUplusDQ[iZone]= solver->CreateFieldProxy(grid); 

        //! The delta Flux is initialized by zero, which will be used in the forward sweep.
        solver->FillField(grid, LUplusDQ[iZone], zero);

        //! Warning: it is important to research the effect of the initialization in future!
        //! e.g., using rhs to init when nsweep > 1.
        //! Initialization.
        if (nSweep == 1)
        {
            //! Single sweep, init by res.
            FieldProxy *rhsProxy = solver->GetFieldProxy(grid, solver->CastName("res")); 
            solver->FillField(grid, dqProxy[iZone], rhsProxy);
            delete rhsProxy;
        }
        else
        {
            //! Multi-sweep, init by zero.
            solver->FillField(grid, dqProxy[iZone], zero);
        }
    }
}

void Controller::MatrixLUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy **dqProxy, FieldProxy **LUplusDQ)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        //! Memory allocating.
        dqProxy[iZone] = solver->CreateFieldProxy(grid);
        LUplusDQ[iZone]= solver->CreateFieldProxy(grid);

        //! The delta Flux is initialized by zero, which will be used in the forward sweep.
        solver->FillField(grid, LUplusDQ[iZone], zero);

        //! Warning: it is important to research the effect of the initialization in future!
        //! e.g., using rhs to init when nsweep > 1.
        //! Initialization.
        if (nSweep == 1)
        {
            //! Single sweep, init by res.
            FieldProxy *rhsProxy = solver->GetFieldProxy(grid, solver->CastName("res")); 
            solver->FillField(grid, dqProxy[iZone], rhsProxy);
            delete rhsProxy;
        }
        else
        {
            //! Multi-sweep, init by zero.
            solver->FillField(grid, dqProxy[iZone], zero);
        }
        //*/    
    }
}

void Controller::LUSGS(int iSolver, int level)
{
    Diagonal(iSolver, level);

    int nLocalZones = PHMPI::GetNumberofLocalZones();
    int nLUSGSSweeps = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");

    if (level != 0)
    {
        nLUSGSSweeps = 1;
    }

    RDouble LUSGSTolerance = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");

    FieldProxy **dqProxy = new FieldProxy *[nLocalZones];

    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy **LUplusDQ= new FieldProxy *[nLocalZones];

    LUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy, LUplusDQ);

    RDouble sweepNormalInit = 0.0, sweepEpsilon = 1.1, globleNormal = 0.0;
    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++ iSweep)
    {
        RDouble sweepNormal = 0.0;

        if (systemGridType == MIXGRID || isFVMOrFDM == FD_METHOD)
        {
            if (iSolver == NS_SOLVER)
            {
                CommunicationDQ(iSolver, level, dqProxy);
            }

            ForwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);

            BackwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);
        }
        else
        {
            ForwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);

            CommunicationDQ(iSolver, level, dqProxy);

            BackwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);
        }

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if (iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }
        else
        {
            sweepEpsilon = sweepNormal / sweepNormalInit;
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

    UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; ++ iLocalZone)
    {
        delete dqProxy[iLocalZone];    dqProxy[iLocalZone] = nullptr;
        delete LUplusDQ[iLocalZone];    LUplusDQ[iLocalZone] = nullptr;
    }
    delete [] dqProxy;    dqProxy = nullptr;
    delete [] LUplusDQ;    LUplusDQ = nullptr;
}

void Controller::DualTimeStepsLUSGS(int iSolver, int level)
{
    Diagonal(iSolver, level);

    int nLocalZones = PHMPI::GetNumberofLocalZones();
    int nLUSGSSweeps = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");

    if (level != 0)
    {
        nLUSGSSweeps = 1;
    }

    RDouble LUSGSTolerance = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");

    FieldProxy **dqProxy = new FieldProxy *[nLocalZones];
    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy **LUplusDQ= new FieldProxy *[nLocalZones];

    FieldProxy **dqProxy1 = new FieldProxy *[nLocalZones];
    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy **LUplusDQ1= new FieldProxy *[nLocalZones];

    LUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy, LUplusDQ);
    LUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy1, LUplusDQ1);

    //! The basic iteration step.
    RDouble sweepNormalInit = 0.0, sweepEpsilon = 1.1, globleNormal = 0.0;
    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++ iSweep)
    {
        RDouble sweepNormal = 0.0;

        if (systemGridType == MIXGRID || isFVMOrFDM == FD_METHOD)
        {
            if (iSolver == NS_SOLVER)
            {
                CommunicationDQ(iSolver, level, dqProxy);
            }

            ForwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);

            BackwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);
        }
        else
        {
            ForwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);

            CommunicationDQ(iSolver, level, dqProxy);

            BackwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);
        }

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if (iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }
        else
        {
            sweepEpsilon = sweepNormal / sweepNormalInit;
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

    //! The advancing time iteration step.
    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++ iSweep)
    {
        RDouble sweepNormal = 0.0;

        if (systemGridType == MIXGRID || isFVMOrFDM == FD_METHOD)
        {
            if (iSolver == NS_SOLVER)
            {
                CommunicationDQ(iSolver, level, dqProxy1);
            }

            ForwardSweep(iSolver, level, dqProxy1, LUplusDQ1, sweepNormal, true);

            BackwardSweep(iSolver, level, dqProxy1, LUplusDQ1, sweepNormal, true);
        }
        else
        {
            ForwardSweep(iSolver, level, dqProxy1, LUplusDQ1, sweepNormal, true);

            CommunicationDQ(iSolver, level, dqProxy1);

            BackwardSweep(iSolver, level, dqProxy1, LUplusDQ1, sweepNormal, true);
        }

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if (iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }
        else
        {
            sweepEpsilon = sweepNormal / sweepNormalInit;
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

    //! Determine the CFL number of the next iteration step.
    DetermineCFLNumber(iSolver, level, dqProxy, dqProxy1);

    UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; ++ iLocalZone)
    {
        delete dqProxy[iLocalZone];    dqProxy[iLocalZone] = nullptr;
        delete LUplusDQ[iLocalZone];    LUplusDQ[iLocalZone] = nullptr;
        delete dqProxy1[iLocalZone];    dqProxy1[iLocalZone] = nullptr;
        delete LUplusDQ1[iLocalZone];    LUplusDQ1[iLocalZone] = nullptr;
    }
    delete [] dqProxy;    dqProxy = nullptr;
    delete [] LUplusDQ;    LUplusDQ = nullptr;
    delete [] dqProxy1;    dqProxy1 = nullptr;
    delete [] LUplusDQ1;    LUplusDQ1 = nullptr;
}

void Controller::ForwardSweep(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones]();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        CFDSolver *solver = GetCFDSolver(zoneID,iSolver);

        zoneNormal[iZone] = 0.0;
        solver->SolveLUSGSForward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone], iAdvanceStep);
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete [] zoneNormal;    zoneNormal = nullptr;
}

void Controller::BackwardSweep(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones]();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        CFDSolver *solver = GetCFDSolver(zoneID,iSolver);

        zoneNormal[iZone] = 0.0;
        solver->SolveLUSGSBackward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone], iAdvanceStep);
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete [] zoneNormal;    zoneNormal = nullptr;
}

void Controller::UpdateFlowField(int iSolver, int level, FieldProxy **dqProxy)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }
        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        int gridtype = grid->Type();
        if (iSolver == NS_SOLVER && isFVMOrFDM == FD_METHOD && gridtype == STRUCTGRID)
        {
            FieldProxy *qtemp_Proxy = solver->CreateFieldProxy(grid);
            FieldProxy *qfield      = solver->GetFieldProxy(grid, solver->CastName("q"));
            solver->FillField(grid, qtemp_Proxy, qfield);

            solver->UpdateFlowField(grid, qtemp_Proxy, dqProxy[iZone]);
            delete qtemp_Proxy;
        }
        else
        {
            FieldProxy *qProxy = solver->GetFieldProxy(grid, solver->CastName("q"));
            solver->UpdateFlowField(grid, qProxy, dqProxy[iZone]);
            delete qProxy;
        }
    }
}

void Controller::DetermineCFLNumber(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **dqProxy1)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->DetermineCFLNumber(grid, dqProxy[iZone], dqProxy1[iZone]);
    }
}

void Controller::CommunicateGradientandLimitOnlyForMixGrid()
{
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

                //CFDSolver *solver = GetCFDSolver(iZone, iSolver);
                CFDSolver *solver = GetCFDSolver(iZone, 0);

                //! Compress the send information into the actkey.
                Grid *grid = PHSPACE::GetGrid(iZone, 0);

                solver->CompressGradientandLimitToInterfaceOnlyForMixGrid(sendBuffer, grid, neighborZone);
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

            //if (sendProcessor == receiveProcessor)
            //{
            //    continue;
            //}

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
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                //CFDSolver *solver = GetCFDSolver(neighborZone, iSolver);
                CFDSolver *solver = GetCFDSolver(neighborZone, 0);

                Grid *grid = PHSPACE::GetGrid(neighborZone, 0);

                solver->DecompressGradientandLimitToInterfaceOnlyForMixGrid(receiveData, grid, iZone);

                ++count;
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

void Controller::CommunicateQlQrOnlyforMixGrid()
{
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (systemGridType != MIXGRID)
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

                //CFDSolver *solver = GetCFDSolver(iZone, iSolver);
                CFDSolver *solver = GetCFDSolver(iZone, 0);
                //! Compress the send information into the actkey.
                Grid *grid = PHSPACE::GetGrid(iZone, 0);

                InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
                if (!interfaceInformation)
                {
                }
                else
                {
                    if (gridtypeofiZone == gridtypeofneighborZone)
                    {
                        RDouble uploadData = 0.0;
                        sendBuffer->MoveToBegin();
                        PHWrite(sendBuffer, uploadData);
                    }
                    else
                    {
                        solver->CompressQlQrToInterfaceOnlyForMixGrid(sendBuffer, grid, neighborZone);
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

            //if (sendProcessor == receiveProcessor)
            //{
            //    continue;
            //}

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
                //CFDSolver *solver = GetCFDSolver(neighborZone, iSolver);
                CFDSolver *solver = GetCFDSolver(neighborZone, 0);

                Grid *grid = PHSPACE::GetGrid(neighborZone, 0);

                InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
                if (!interfaceInformation)
                {
                    //! can not use return here, because ++count in the after
                }
                else
                {
                    if (gridtypeofiZone == gridtypeofneighborZone)
                    {
                        RDouble downdData = 0.0;
                        receiveData->MoveToBegin();
                        PHRead(receiveData, downdData);
                    }
                    else
                    {
                        solver->DecompressQlQrToInterfaceOnlyForMixGrid(receiveData, grid, iZone);
                    }
                }

                ++count;
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

void Controller::CommunicateFacefluxOnlyforMixGrid()
{
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (systemGridType != MIXGRID)
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

                //CFDSolver *solver = GetCFDSolver(iZone, iSolver);
                CFDSolver *solver = GetCFDSolver(iZone, 0);
                //! Compress the send information into the actkey.
                Grid *grid = PHSPACE::GetGrid(iZone, 0);

                solver->CompressFacefluxToInterfaceOnlyForMixGrid(sendBuffer, grid, neighborZone);
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
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                //CFDSolver *solver = GetCFDSolver(neighborZone, iSolver);
                CFDSolver *solver = GetCFDSolver(neighborZone, 0);

                Grid *grid = PHSPACE::GetGrid(neighborZone, 0);

                solver->DecompressFacefluxToInterfaceOnlyForMixGrid(receiveData, grid, iZone);

                ++count;
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

void Controller::CommunicationDQ(int iSolver, int level, FieldProxy **dqProxy)
{
    if (iSolver != NS_SOLVER)
    {
        //! Only consider the communication of NS-Solver.
        return;
    }

    int swapDq = GlobalDataBase::GetIntParaFromDB("swapDq");
    if (swapDq == 0)
    {
        return;
    }

    int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(0);
    PHSolver *solver = GetCFDSolver(zoneID, iSolver);
    int nEquation  = solver->GetNumberOfEquations();

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

                int localZoneID = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
                CFDSolver *solver1 = GetCFDSolver(iZone, iSolver);

                //! Compress the send information into the actkey.
                Grid *grid = PHSPACE::GetGrid(iZone, level);
                solver1->CompressDQ(sendBuffer, dqProxy[localZoneID], grid, iZone, neighborZone, nEquation);
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
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                int localZoneID = region->GetProcessGlobalZoneIndexToLocalZoneIndex(neighborZone);
                CFDSolver *solver2 = GetCFDSolver(neighborZone, iSolver);

                Grid *grid = PHSPACE::GetGrid(neighborZone, level);
                solver2->DecompressDQ(receiveData, dqProxy[localZoneID], grid, neighborZone, iZone, nEquation);

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

void Controller::CommunicateGenericArray(int iSolver, int level)
{
    using namespace PHMPI;
    int zoneID = GetLocalZoneIDToGlobalZoneID(0);
    CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

    if (iSolver == NS_SOLVER)
    {
        int nEquation = solver->GetNumberOfEquations();
        CommunicateSpecificArray(level, "q", nEquation);

        int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
        CommunicateSpecificArray(level, "t", nTemperatureModel);
    }
    else if (iSolver == TURB_SOLVER)
    {
        int nEquation = solver->GetNumberOfEquations();
        CommunicateSpecificArray(level, "q_turb", nEquation);
    }
    else if (iSolver == TRANSITION_SOLVER)
    {
        int nEquation = solver->GetNumberOfEquations();
        CommunicateSpecificArray(level, "q_transition", nEquation);
    }
}

void Controller::CommunicateGradientArray(int iSolver, int level)
{
    using namespace PHMPI;
    int zoneID = GetLocalZoneIDToGlobalZoneID(0);
    CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int isFVMOrFDM = PHSPACE::GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    int iLES = PHSPACE::GlobalDataBase::GetIntParaFromDB("iLES");

    //! gradient has been communicated in another place for isFVMOrFDM == FDM.
    if (systemGridType == STRUCTGRID && isFVMOrFDM == FV_METHOD)
    {
        int nZones = PHMPI::GetNumberofGlobalZones();
        int currentProcessor = PHMPI::GetCurrentProcessorID();

        if (iSolver == NS_SOLVER)
        {
            //! NSSolver, GradientCommunication for struct high order has been done at another place, with another function.
            int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
            int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
            int nNSEquation = GlobalDataBase::GetIntParaFromDB("nm");
            int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
            int nSpeciesNumber = nLaminar + nChemical - nNSEquation;

            int nTotalVariable = nTemperatureModel + 3;
            CommunicateSpecificArray(level, "gradUVWTCellCenterX", nTotalVariable);
            CommunicateSpecificArray(level, "gradUVWTCellCenterY", nTotalVariable);
            CommunicateSpecificArray(level, "gradUVWTCellCenterZ", nTotalVariable);

            if (nChemical > 0 && nSpeciesNumber > 0)
            {
                CommunicateSpecificArray(level, "dcdx", nSpeciesNumber);
                CommunicateSpecificArray(level, "dcdy", nSpeciesNumber);
                CommunicateSpecificArray(level, "dcdz", nSpeciesNumber);
            }
            if (periodicType == ROTATIONAL_PERIODICITY)
            {
                for (int iZone = 0; iZone < nZones; ++ iZone)
                {
                    int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
                    ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                    int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
                    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
                    {
                        int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
                        int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);
                        if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                        {
                            continue;
                        }
                        if (currentProcessor == sendProcessor)
                        {
                            CFDSolver *solver1 = GetCFDSolver(iZone, iSolver);
                            Grid *grid = PHSPACE::GetGrid(iZone, level);
                            solver1->RotateVectorFromInterface(grid, iZone, nTotalVariable);
                        }
                    }
                }
            }
        }
        else if (iSolver == TURB_SOLVER && iLES == NOLES_SOLVER)
        {
            int nEquation = solver->GetNumberOfEquations();
            CommunicateSpecificArray(level, "gradTurbulenceCellCenterX", nEquation);
            CommunicateSpecificArray(level, "gradTurbulenceCellCenterY", nEquation);
            CommunicateSpecificArray(level, "gradTurbulenceCellCenterZ", nEquation);
            if (periodicType == ROTATIONAL_PERIODICITY)
            {
                for (int iZone = 0; iZone < nZones; ++ iZone)
                {
                    int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
                    ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                    int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
                    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
                    {
                        int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
                        int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);
                        if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                        {
                            continue;
                        }
                        if (currentProcessor == sendProcessor)
                        {
                            CFDSolver *solver1 = GetCFDSolver(iZone, iSolver);
                            Grid *grid = PHSPACE::GetGrid(iZone, level);
                            solver1->RotateVectorFromInterface(grid, iZone, nEquation);
                        }
                    }
                }
            }
        }
        else if (iSolver == TRANSITION_SOLVER)
        {
            int nEquation = solver->GetNumberOfEquations();
            CommunicateSpecificArray(level, "gradTransitionCellCenterX", nEquation);
            CommunicateSpecificArray(level, "gradTransitionCellCenterY", nEquation);
            CommunicateSpecificArray(level, "gradTransitionCellCenterZ", nEquation);
            if (periodicType == ROTATIONAL_PERIODICITY)
            {
                for (int iZone = 0; iZone < nZones; ++ iZone)
                {
                    int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
                    ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                    int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
                    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
                    {
                        int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
                        int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);
                        if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                        {
                            continue;
                        }
                        if (currentProcessor == sendProcessor)
                        {
                            CFDSolver *solver1 = GetCFDSolver(iZone, iSolver);
                            Grid *grid = PHSPACE::GetGrid(iZone, level);
                            solver1->RotateVectorFromInterface(grid, iZone, nEquation);
                        }
                    }
                }
            }
        }
        else if (iSolver == TURB_SOLVER && iLES == LES_SOLVER)
        {
            //! For LES, Gradient Communication should be done another time after NSSolver ends.
            int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
            int nTotalVariable = nTemperatureModel + 3;

            CommunicateSpecificArray(level, "gradUVWTCellCenterX", nTotalVariable);
            CommunicateSpecificArray(level, "gradUVWTCellCenterY", nTotalVariable);
            CommunicateSpecificArray(level, "gradUVWTCellCenterZ", nTotalVariable);

            if (periodicType == ROTATIONAL_PERIODICITY)
            {
                for (int iZone = 0; iZone < nZones; ++ iZone)
                {
                    int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
                    ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                    int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
                    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
                    {
                        int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
                        int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);
                        if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                        {
                            continue;
                        }
                        if (currentProcessor == sendProcessor)
                        {
                            CFDSolver *solver1 = GetCFDSolver(iZone, iSolver);
                            Grid *grid = PHSPACE::GetGrid(iZone, level);
                            solver1->RotateVectorFromInterface(grid, iZone, nTotalVariable);
                        }
                    }
                }
            }
        }
    }
    else if (systemGridType == MIXGRID && isFVMOrFDM == FV_METHOD)
    {
        if (iSolver == NS_SOLVER)
        {
            CommunicateGradientandLimitOnlyForMixGrid();
        }
        else if (iSolver == TURB_SOLVER)
        {
        }
    }
    if (systemGridType == UNSTRUCTGRID)
    {
        int nZones = PHMPI::GetNumberofGlobalZones();
        int currentProcessor = PHMPI::GetCurrentProcessorID();

        if (iSolver == NS_SOLVER)
        {
            //! NSSolver, GradientCommunication for struct high order has been done at another place, with another function.
            int nNSEquation = GlobalDataBase::GetIntParaFromDB("nm");
            int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
            CommunicateSpecificArray(level, "gradPrimtiveVarX", nNSEquation);
            CommunicateSpecificArray(level, "gradPrimtiveVarY", nNSEquation);
            CommunicateSpecificArray(level, "gradPrimtiveVarZ", nNSEquation);

            CommunicateSpecificArray(level, "gradTemperatureX", nTemperatureModel);
            CommunicateSpecificArray(level, "gradTemperatureY", nTemperatureModel);
            CommunicateSpecificArray(level, "gradTemperatureZ", nTemperatureModel);

            if (periodicType == ROTATIONAL_PERIODICITY)
            {
                for (int iZone = 0; iZone < nZones; ++iZone)
                {
                    int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
                    ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                    int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

                    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++iNeighbor)
                    {
                        int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
                        int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

                        if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                        {
                            continue;
                        }

                        if (currentProcessor == sendProcessor)
                        {
                            CFDSolver *solver1 = GetCFDSolver(iZone, iSolver);
                            Grid *grid = PHSPACE::GetGrid(iZone, level);
                            solver1->RotateVectorFromInterface(grid, iZone, nNSEquation);
                            solver1->RotateVectorFromInterface(grid, iZone, nTemperatureModel);
                        }
                    }
                }
            }

        }
        else if (iSolver == TURB_SOLVER && iLES == NOLES_SOLVER)
        {
            int nEquation = solver->GetNumberOfEquations();
            CommunicateSpecificArray(level, "gradTurbulenceX", nEquation);
            CommunicateSpecificArray(level, "gradTurbulenceY", nEquation);
            CommunicateSpecificArray(level, "gradTurbulenceZ", nEquation);

            if (periodicType == ROTATIONAL_PERIODICITY)
            {
                for (int iZone = 0; iZone < nZones; ++iZone)
                {
                    int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
                    ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                    int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
                    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++iNeighbor)
                    {
                        int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
                        int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);
                        if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                        {
                            continue;
                        }
                        if (currentProcessor == sendProcessor)
                        {
                            CFDSolver *solver1 = GetCFDSolver(iZone, iSolver);
                            Grid *grid = PHSPACE::GetGrid(iZone, level);
                            solver1->RotateVectorFromInterface(grid, iZone, nEquation);
                        }
                    }
                }
            }
        }
        else if (iSolver == TRANSITION_SOLVER)
        {
            int nEquation = solver->GetNumberOfEquations();
            CommunicateSpecificArray(level, "gradTransitionX", nEquation);
            CommunicateSpecificArray(level, "gradTransitionY", nEquation);
            CommunicateSpecificArray(level, "gradTransitionZ", nEquation);

            if (periodicType == ROTATIONAL_PERIODICITY)
            {
                for (int iZone = 0; iZone < nZones; ++iZone)
                {
                    int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
                    ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                    int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
                    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++iNeighbor)
                    {
                        int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
                        int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);
                        if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                        {
                            continue;
                        }
                        if (currentProcessor == sendProcessor)
                        {
                            CFDSolver *solver1 = GetCFDSolver(iZone, iSolver);
                            Grid *grid = PHSPACE::GetGrid(iZone, level);
                            solver1->RotateVectorFromInterface(grid, iZone, nEquation);
                        }
                    }
                }
            }
        }
    }
}

void Controller::CommunicationInterfaceData(int iSolver)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    //! Face Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        solver->CommunicationInterfaceData();
    }
}

void Controller::CommunicationInterfaceDataOneLevel(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    //! Face Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        solver->CommunicationInterfaceData(level);
    }
}

void Controller::BLUSGS(int iSolver, int level)
{
    //! If the solver is Turb Solver, it will take the LUSGS methods.
    if (iSolver != NS_SOLVER)
    {
        LUSGS(iSolver, level);
        return;
    }
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    int nLUSGSSweeps = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");
    RDouble LUSGSTolerance = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");

    FieldProxy **dqProxy = new FieldProxy *[nLocalZones];

    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy **LUplusDQ= new FieldProxy *[nLocalZones];

    BLUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy, LUplusDQ);

    RDouble sweepNormalInit = 0.0, sweepEpsilon = 1.1, globleNormal = 0.0;
    
    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++ iSweep)
    {
        RDouble sweepNormal = 0.0;
        ForwardSweepBLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal,iSweep);

        CommunicationDQ(iSolver, level, dqProxy);

        BackwardSweepBLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal);

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if (iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }else
        {
            sweepEpsilon = sqrt(sweepNormal / sweepNormalInit);
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

    UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; ++ iLocalZone)
    {
        delete dqProxy[iLocalZone];    dqProxy[iLocalZone] = nullptr;
        delete LUplusDQ[iLocalZone];    LUplusDQ[iLocalZone] = nullptr;
    }
    delete [] dqProxy;    dqProxy = nullptr;
    delete [] LUplusDQ;    LUplusDQ = nullptr;
}

void Controller::ForwardSweepBLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal , int &iSweep)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones]();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        if (iSweep == 0)
        {
            solver->LHSMatrixs(grid);
        }

        zoneNormal[iZone] = 0.0;
        solver->SolveBLUSGSForward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone]);//! The Block LU-SGS method.
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete [] zoneNormal;    zoneNormal = NULL;
}

void Controller::BackwardSweepBLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones]();
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }
        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        zoneNormal[iZone] = 0.0;
        solver->SolveBLUSGSBackward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone]);//! The Block LU-SGS method
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete [] zoneNormal;    zoneNormal = NULL;
}

void Controller::MatrixLUSGS(int iSolver, int level)
{
    //! If the solver is Turb Solver, it will take the LUSGS methods.
    if (iSolver != NS_SOLVER)
    {
        LUSGS(iSolver, level);
        return;
    }
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    int nLUSGSSweeps = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");
    RDouble LUSGSTolerance = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");

    FieldProxy **dqProxy = new FieldProxy *[nLocalZones];

    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy **LUplusDQ= new FieldProxy *[nLocalZones];

    MatrixLUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy, LUplusDQ);

    RDouble sweepNormalInit = 1.0, sweepEpsilon = 1.1, globleNormal;
    
    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++ iSweep)
    {
        RDouble sweepNormal = 0.0;
        ForwardSweepMatrixLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal,iSweep);

        CommunicationDQ(iSolver, level, dqProxy);

        BackwardSweepMatrixLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal);

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if (iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }else
        {
            sweepEpsilon = sqrt(sweepNormal / sweepNormalInit);
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

    UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; ++ iLocalZone)
    {
        delete dqProxy[iLocalZone];    dqProxy[iLocalZone] = nullptr;
        delete LUplusDQ[iLocalZone];    LUplusDQ[iLocalZone] = nullptr;
    }
    delete [] dqProxy;    dqProxy = nullptr;
    delete [] LUplusDQ;    LUplusDQ = nullptr;
}

void Controller::ForwardSweepMatrixLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal , int &iSweep)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones];

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int izone = 0; izone < nLocalZones; ++ izone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        if (iSweep == 0)
        {
            solver->LHSMatrixs(grid);
        }

        zoneNormal[izone] = 0.0;
        solver->SolveMatrixLUSGSForward(grid, dqProxy[izone], LUplusDQ[izone], zoneNormal[izone]); //! The Matrix LU-SGS method.
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete [] zoneNormal;
}

void Controller::BackwardSweepMatrixLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones];
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int izone = 0; izone < nLocalZones; ++ izone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }
        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        zoneNormal[izone] = 0.0;
        solver->SolveMatrixLUSGSBackward(grid, dqProxy[izone], LUplusDQ[izone], zoneNormal[izone]);//! The Matrix LU-SGS method
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete [] zoneNormal;
}

void Controller::LineLUSGS(int iSolver, int level)
{
    //! If the solver is Turb Solver, it will take the LUSGS methods.
    if (iSolver != NS_SOLVER)
    {
        LUSGS(iSolver, level);
        return;
    }
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    int nLUSGSSweeps = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");
    RDouble LUSGSTolerance = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");

    FieldProxy** dqProxy = new FieldProxy * [nLocalZones];

    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy** LUplusDQ = new FieldProxy * [nLocalZones];

    LineLUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy, LUplusDQ);

    RDouble sweepNormalInit = 0.0, sweepEpsilon = 1.1, globleNormal = 0.0;

    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++iSweep)
    {
        RDouble sweepNormal = 0.0;
        ForwardSweepLineLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal, iSweep);

        CommunicationDQ(iSolver, level, dqProxy);

        BackwardSweepLineLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal);

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if (iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }
        else
        {
            sweepEpsilon = sweepNormal / sweepNormalInit;
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

    UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; ++iLocalZone)
    {
        delete dqProxy[iLocalZone];    dqProxy[iLocalZone] = nullptr;
        delete LUplusDQ[iLocalZone];    LUplusDQ[iLocalZone] = nullptr;
    }
    delete [] dqProxy;    dqProxy = nullptr;
    delete [] LUplusDQ;    LUplusDQ = nullptr;
}

void Controller::LineLUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy** dqProxy, FieldProxy** LUplusDQ)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid* grid = GetGrid(zoneID, level);
        CFDSolver* solver = GetCFDSolver(zoneID, iSolver);

        //! Memory allocating.
        dqProxy[iZone] = solver->CreateFieldProxy(grid);
        LUplusDQ[iZone] = solver->CreateFieldProxy(grid);

        //! The delta Flux is initialized by zero, which will be used in the forward sweep.
        solver->FillField(grid, LUplusDQ[iZone], zero);

        //! Warning: it is important to research the effect of the initialization in future!
        //! e.g., using rhs to init when nsweep > 1.
        //! Initialization.
        if (nSweep == 1)
        {
            //! Single sweep, init by res.
            FieldProxy* rhsProxy = solver->GetFieldProxy(grid, solver->CastName("res"));
            solver->FillField(grid, dqProxy[iZone], rhsProxy);
            delete rhsProxy;
        }
        else
        {
            //! Multi-sweep, init by zero.
            solver->FillField(grid, dqProxy[iZone], zero);
        }
    }
}

void Controller::ForwardSweepLineLUSGS(int iSolver, int level, FieldProxy** dqProxy, FieldProxy** LUplusDQ, RDouble& sweepNormal, int& iSweep)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble* zoneNormal = new RDouble[nLocalZones]();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid* grid = GetGrid(zoneID, level);
        CFDSolver* solver = GetCFDSolver(zoneID, iSolver);

        if (iSweep == 0)
        {
            solver->LineLUSGSMatrixs(grid);
        }

        zoneNormal[iZone] = 0.0;
        solver->SolveLineLUSGSForward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone]);//! The Block LU-SGS method.
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete [] zoneNormal;    zoneNormal = nullptr;
}

void Controller::BackwardSweepLineLUSGS(int iSolver, int level, FieldProxy** dqProxy, FieldProxy** LUplusDQ, RDouble& sweepNormal)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble* zoneNormal = new RDouble[nLocalZones]();
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }
        Grid* grid = GetGrid(zoneID, level);
        CFDSolver* solver = GetCFDSolver(zoneID, iSolver);

        zoneNormal[iZone] = 0.0;
        solver->SolveLineLUSGSBackward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone]);//The Block LU-SGS method
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete [] zoneNormal;    zoneNormal = nullptr;
}

void Controller::RestrictAllQ(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        Grid *coarseGrid = grid->GetCoarseGrid();
        if (coarseGrid)
        {
            CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

            solver->RestrictAllQ(grid, coarseGrid);
        }
    }
}

void Controller::RestrictDefect(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        Grid * coarseGrid = grid->GetCoarseGrid();
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->RestrictDefect(grid, coarseGrid);
    }
}

void Controller::PutCorrection(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **qProxy = GetQProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        solver->PutCorrection(grid, qProxy[iZone]);
    }
}

void Controller::PutCorrectionBack(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **qProxy = GetQProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        solver->PutCorrectionBack(grid, qProxy[iZone]);
    }
}

void Controller::CorrectFineGrid(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->CorrectFineGrid(grid, grid->GetCoarseGrid());
    }
}

Zone * Controller::GetZone(int iZone)
{
    return region->GetZone(iZone);
}

CFDSolver * Controller::GetCFDSolver(int iZone, int iSolver)
{
    return static_cast<CFDSolver *>(GlobalSolvers::GetSolver(iZone, iSolver));
}

#ifdef USE_GMRESSOLVER
//! For GMRES -- Boqian Wang
void Controller::GMRESSolver(int iSolver, int level)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();

    FieldProxy **dqProxy = new FieldProxy *[nLocalZones];

    for (int iZone = 0; iZone < nLocalZones; iZone++)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if(!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        dqProxy[iZone] = solver->CreateFieldProxy(grid);
        solver->FillField(grid, dqProxy[iZone], zero);  // difference of conservative variables, init by zero
        solver->GMRESSolver(grid, dqProxy[iZone]);
    }

    // GMRESWholeJacobian
    UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; iLocalZone++)
    {
        delete dqProxy[iLocalZone];
    }
    delete []dqProxy;
}

//! GMRES Coupled for NS and Turb
void Controller::GMRESSolver_Coupled(int level)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    if( tscheme != GMRES )
    {
        return;
    }

    FieldProxy **dqProxyNS      = new FieldProxy *[nLocalZones];
    FieldProxy **dqProxyTurb    = new FieldProxy *[nLocalZones];

    for (int iZone = 0; iZone < nLocalZones; iZone++)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if(!zone)
        {
            continue;
        }

        Grid *grid0 = GetGrid(zoneID, level);

        CFDSolver *solverNS     = GetCFDSolver(zoneID, 0);
        CFDSolver *solverTurb   = GetCFDSolver(zoneID, 1);

        dqProxyNS[iZone] = solverNS->CreateFieldProxy(grid0);
        solverNS->FillField(grid0, dqProxyNS[iZone], zero);  // difference of conservative variables, init by zero

        dqProxyTurb[iZone] = solverTurb->CreateFieldProxy(grid0);
        solverTurb->FillField(grid0, dqProxyTurb[iZone], zero);  // difference of turbulence variables, init by zero
        

        UnstructGrid* grid = static_cast<UnstructGrid *>(grid0);
        
        RDouble **dRdq                  = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdq"));
        RDouble **dRdq_turb             = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdq_turb"));
        RDouble **res                   = reinterpret_cast<RDouble **>(grid->GetDataPtr("res"));
        RDouble **res_turb              = reinterpret_cast<RDouble **>(grid->GetDataPtr("res_turb"));
        RDouble** dRdqCoupledTerm       = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdqCoupledTerm"));
        RDouble** dRdqCoupledTerm_turb  = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdqCoupledTerm_turb"));

        RDouble **dDdP      = reinterpret_cast<RDouble**>(grid->GetDataPtr("dDdP"));

        // RDouble **dRdqCoupledTerm = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdqCoupledTerm"));

        int nTotalCells         = grid->GetNTotalCell();
        int nBoundFace          = grid->GetNBoundFace();
        int nEquations          = solverNS->GetNumberOfEquations(); 
        int nTurbulenceEquation = solverTurb->GetNumberOfEquations();
        
        RDouble* vol = grid->GetCellVolume();
        RDouble* dt = reinterpret_cast<RDouble*> (grid->GetDataPtr("dt"));
        vector<int> AI = grid->GetJacobianAI4GMRES();  // GMRESSparse GMRESCSR
        vector<int> AJ = grid->GetJacobianAJ4GMRES();  // GMRESSparse GMRESCSR

        // GMRESCSR
        for (int icell = 0; icell < nTotalCells; icell++)
        {
            vector<int>::iterator result = find(AJ.begin() + AI[icell], AJ.begin() + AI[icell + 1], icell); // find qcell
            int index = distance(AJ.begin(), result);
            
            int indexNS = index*nEquations;
            for (int iequ = 0; iequ < nEquations; iequ++)
            {
                 dRdq[iequ][indexNS + iequ] += 1.0 / dt[icell];
            }

            int indexTurb = index*nTurbulenceEquation;
            for (int iequ = 0; iequ < nTurbulenceEquation; iequ++)
            {
                 dRdq_turb[iequ][indexTurb + iequ] += 1.0 / dt[icell];
            }
        }

        // GMRESJac1st
        vector<int> AI1st = grid->GetJacobianAI1st4GMRES(); 
        vector<int> AJ1st = grid->GetJacobianAJ1st4GMRES();
        RDouble **dRdq1st = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdq1st"));
        int JacOrder    = grid->GetJacobianOrder();

        // GMRESCoupled1st
        RDouble** dRdqCoupledTerm1st        = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdqCoupledTerm1st"));
        RDouble** dRdqCoupledTerm_turb1st   = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdqCoupledTerm_turb1st"));
        RDouble** dRdq_turb1st              = reinterpret_cast<RDouble**> (grid->GetDataPtr("dRdq_turb1st"));       

        if (JacOrder == 2)
        {
            int jmax1st = AI1st[nTotalCells + nBoundFace];
            int rcell1st = 0;
            int rcell2nd = 0;
            int j1st = 0;
            do
            {
                for (j1st = AI1st[rcell1st]; j1st < AI1st[rcell1st + 1]; j1st++)
                {
                    int col1st = AJ1st[j1st];
                    vector<int>::iterator result = find(AJ.begin() + AI[rcell2nd], AJ.begin() + AI[rcell2nd + 1], col1st); // find qcell
                    int index = distance(AJ.begin(), result);

                    int indexNS     = index*nTurbulenceEquation;
                    int indexNS1st  = j1st*nTurbulenceEquation;
                    for (int iequI = 0; iequI < nEquations; iequI++)
                    {
                        dRdqCoupledTerm1st[iequI][indexNS1st] = dRdqCoupledTerm[iequI][indexNS];
                    }
                    dRdq_turb1st[0][indexNS1st] = dRdq_turb[0][indexNS];

                    int indexTurb           = index*nEquations;
                    int indexTurb1st        = j1st*nEquations;
                    for (int iequI = 0; iequI < nEquations; iequI++)
                    {
                        dRdqCoupledTerm_turb1st[0][indexTurb1st+iequI] = dRdqCoupledTerm_turb[0][indexTurb+iequI];
                    }
                }
                rcell1st++;
                rcell2nd++;

            } while (j1st < jmax1st);
            

            for (int icell = 0; icell < nTotalCells; icell++)
            {
                vector<int>::iterator result = find(AJ1st.begin() + AI1st[icell], AJ1st.begin() + AI1st[icell + 1], icell); // find qcell
                int index = distance(AJ1st.begin(), result);
                if (result == AJ1st.begin() + AI1st[icell + 1])
                {
                    printf("\nrcell %d cannot find qcell %d from \n", icell, icell);
                    for (int i = AI[icell]; i < AI[icell + 1]; i++)
                    {
                        printf("%d ", AJ1st[i]);
                    }
                    printf("\n");
                    TK_Exit::ExceptionExit("Sparse matrix index is wrong");
                }

                int indexNS = index*nEquations;
                for (int iequ = 0; iequ < nEquations; iequ++)
                {
                     dRdq1st[iequ][indexNS + iequ] += 1.0 / dt[icell];
                }

                int indexTurb = index*nTurbulenceEquation;
                for (int iequ = 0; iequ < nTurbulenceEquation; iequ++)
                {
                     dRdq_turb1st[iequ][indexTurb + iequ] += 1.0 / dt[icell];
                }
            }
        }
      
        // entry to PETSc for NS equations\n
        Vec x, b;
        Mat A;
        Mat APC; // GMRESJac1st
        KSP ksp;                        // define solver
        PC pc;                          // define matrix
        PetscInt nTotalSize;            // dimension size of matrix and vector
        PetscInt nLevels = 3;           // define number of levels used by preconditioner ILU
        PetscInt MaxStepsRestart = 500;  // define maximum steps required by [restart]
        PetscInt MaxSteps = 500;        // define maximum steps
        PetscReal rTolerance = 1.e-1;   // define tolerance of relative error // -1
        PetscReal dropTolerance = 1.e5; // define tolerance of drop error

        int blockSize = nEquations + nTurbulenceEquation;
       
        int globalTotalCells;
        globalTotalCells = (int)GlobalDataBase::GetDoubleParaFromDB("GlobalTotalCells");
        int nInterfaceCells;
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        int globalCellIndexShift = grid->GetGlobalCellIndexShift();
        if (interfaceInfo != nullptr)
        {
            nInterfaceCells = interfaceInfo->GetNIFace();
        }
        nTotalSize = blockSize * globalTotalCells;

        // initialize PETSc
        PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);

        // create vectors;
        VecCreateMPI(PETSC_COMM_WORLD, nTotalCells * blockSize, nTotalSize, &x);
        PetscObjectSetName((PetscObject)x, "solution");
        VecSetUp(x);
        VecDuplicate(x, &b);

        // create matrix
        // GMRESCSR
        PetscInt nnz[nTotalCells];
        for (int icell = 0; icell < nTotalCells; icell++)
        {
            int count = 0;
            for (int j = AI[icell]; j < AI[icell + 1]; j++)
            {
                if(AJ[j] < nTotalCells)
                {
                    count++;
                }
            }
            nnz[icell] = count;
        }
        MatCreate(PETSC_COMM_WORLD, &A);
        MatSetSizes(A, nTotalCells * blockSize, nTotalCells * blockSize,  nTotalSize, nTotalSize);
        MatSetType(A, MATMPIBAIJ);
        MatSetBlockSize(A, blockSize);
        MatMPIBAIJSetPreallocation(A, blockSize, PETSC_DEFAULT, nnz, 20, PETSC_NULL);
        PetscInt firstRow = globalTotalCells * blockSize;
        PetscInt lastOneMoreRow = (globalTotalCells + nTotalCells) * blockSize;
        MatGetOwnershipRange(A, &firstRow, &lastOneMoreRow);
        MatSetUp(A);
        MatZeroEntries(A);

        MatCreate(PETSC_COMM_WORLD, &APC);
        MatSetSizes(APC, nTotalCells * blockSize, nTotalCells * blockSize, nTotalSize, nTotalSize);
        MatSetType(APC, MATMPIBAIJ);
        MatSetBlockSize(APC, blockSize);

        // GMRESJac1st
        if( JacOrder == 2)
        {
            PetscInt nnz1st[nTotalCells];
            for (int icell = 0; icell < nTotalCells; icell++)
            {
                int count = 0;
                for (int j = AI1st[icell]; j < AI1st[icell + 1]; j++)
                {
                    if(AJ1st[j] < nTotalCells)
                    {
                        count++;
                    }
                }
                nnz1st[icell] = count;
            }
            MatMPIBAIJSetPreallocation(APC, blockSize, PETSC_DEFAULT, nnz1st, 0, PETSC_NULL);
        }
        else
        {
            MatMPIBAIJSetPreallocation(APC, blockSize, PETSC_DEFAULT, nnz, 0, PETSC_NULL);
        }
        MatGetOwnershipRange(APC, &firstRow, &lastOneMoreRow);
        MatSetUp(APC);
        MatZeroEntries(APC);

        // ----------------
        // |   A       |B |
        // |           |  |
        // |           |  |
        // |           |  |
        // ----------------
        // | C         |D |
        // ----------------
        // improved approach to assemble matrix CSR
        PetscInt petsccol[nEquations];
        PetscInt petscrow;
        int jmax = AI[nTotalCells];
        int rcell = 0;
        int j = AI[rcell];
        do
        {
            int rowidx = (rcell + globalCellIndexShift) * blockSize;
            for (j = AI[rcell]; j < AI[rcell + 1]; j++)
            {
                int qcell = AJ[j]; // get qcell (colume index of matrix)
                // start to insert values: limited to dimension of nTotalSize
                if(qcell < nTotalCells)
                {
                    int colidx = (qcell + globalCellIndexShift) * blockSize;

                    for (int irow = 0; irow < nEquations; irow++)
                    {
                        petscrow = rowidx + irow;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }
                        MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][j * nEquations], ADD_VALUES);  // part A
                        
                        PetscInt colidxturb =  colidx + nEquations;
                        MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdqCoupledTerm[irow][j*nTurbulenceEquation], ADD_VALUES);    // part B
                    }

                    petscrow = rowidx + nEquations;
                    
                    for (int n = 0; n < nEquations; n++)
                    {
                        petsccol[n] = colidx + n;
                    }

                    MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdqCoupledTerm_turb[0][j * nEquations], ADD_VALUES);  // part C
                        
                    PetscInt colidxturb =  colidx + nEquations;
                    MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdq_turb[0][j*nTurbulenceEquation], ADD_VALUES);   // part D

                }
            }
            rcell++; // update rcell, go to the next row index of matrix
        } while (j < jmax);

        //! consider the second order template for 2nd order jacobian matrix for domain decomposition strategy
        if (interfaceInfo != nullptr && JacOrder == 2)
        {
            vector<int> *neighborCells = grid->GMRESGetNeighborCells();
            rcell = nTotalCells;
            jmax = AI[nTotalCells + nBoundFace];
            j = AI[rcell];
            do
            {
                int rowidx = interfaceInfo->MatchedGlobalNeighborCellIndex(rcell) * blockSize;
                if (rowidx >= 0) //! rcell is the ghost cell, rowidx >=0 means this cell locates on the interface boundary
                {
                    int leftcell = interfaceInfo->MatchedLocalPhysicalCellIndex(rcell); // return its corresponding physical cell index
                    int colidx = (leftcell + globalCellIndexShift) * blockSize;
                    vector<int>::iterator result = find(AJ.begin() + AI[rcell], AJ.begin() + AI[rcell + 1], leftcell); // find leftcell
                    int index = distance(AJ.begin(), result);

                    for (int irow = 0; irow < nEquations; irow++)
                    {
                        petscrow = rowidx + irow;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }
                        MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][index * nEquations], ADD_VALUES); // part A

                        PetscInt colidxturb = colidx + nEquations;
                        MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdqCoupledTerm[irow][index * nTurbulenceEquation], ADD_VALUES); // part B
                    }

                    petscrow = rowidx + nEquations;
                    for (int  n = 0; n < nEquations; n++)
                    {
                        petsccol[n] = colidx + n;
                    }
                    MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdqCoupledTerm_turb[0][index * nEquations], ADD_VALUES); // part C
                    PetscInt colidxturb = colidx + nEquations;
                    MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdq_turb[0][index * nTurbulenceEquation], ADD_VALUES); // part D

                    // obtain the list of the neighbor physical cell for the left cell
                    for (int id = 0; id < neighborCells[leftcell].size(); id++)
                    {
                        int neighborIndex = neighborCells[leftcell][id];
                        if (neighborIndex >= nTotalCells)
                        {
                            continue;
                        }

                        int colidx = (neighborIndex + globalCellIndexShift) * blockSize;
                        vector<int>::iterator result = find(AJ.begin() + AI[rcell], AJ.begin() + AI[rcell + 1], neighborIndex); // find neighborIndex
                        int index = distance(AJ.begin(), result);

                        for (int irow = 0; irow < nEquations; irow++)
                        {
                            petscrow = rowidx + irow;
                            for (int n = 0; n < nEquations; n++)
                            {
                                petsccol[n] = colidx + n;
                            }
                            MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][index * nEquations], ADD_VALUES); // part A
                            PetscInt colidxturb = colidx + nEquations;
                            MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdqCoupledTerm[irow][index * nTurbulenceEquation], ADD_VALUES); // part B
                        }

                        petscrow = rowidx + nEquations;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }
                        MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdqCoupledTerm_turb[0][index * nEquations], ADD_VALUES); // part C
                        PetscInt colidxturb = colidx + nEquations;
                        MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdq_turb[0][index * nTurbulenceEquation], ADD_VALUES); // part D
                    }

                    //! obtain the R_physical/q_ghost caused by the neighbor cell at interface bc
                    result = find(AJ.begin() + AI[leftcell], AJ.begin() + AI[leftcell + 1], rcell); // find rcell
                    index = distance(AJ.begin(), result);
                    for (int irow = 0; irow < nEquations; irow++)
                    {
                        petscrow = (leftcell + globalCellIndexShift) * blockSize + irow;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = rowidx + n;
                        }
                        MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][index * nEquations], ADD_VALUES); // part A
                        PetscInt colidxturb = rowidx + nEquations;
                        MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdqCoupledTerm[irow][index * nTurbulenceEquation], ADD_VALUES); // part B
                    }

                    petscrow = (leftcell + globalCellIndexShift) * blockSize + nEquations;
                    for (int n = 0; n < nEquations; n++)
                    {
                        petsccol[n] = rowidx + n;
                    }
                    MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdqCoupledTerm_turb[0][index * nEquations], ADD_VALUES); // part C
                    colidxturb = rowidx + nEquations;
                    MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdq_turb[0][index * nTurbulenceEquation], ADD_VALUES); // part D

                    //! obtain the list of the neighbor physical cell for the left cell
                    for (int id = 0; id < neighborCells[leftcell].size(); id++)
                    {
                        int neighborIndex = neighborCells[leftcell][id];
                        if (neighborIndex >= nTotalCells)
                        {
                            continue;
                        }
                        vector<int>::iterator result = find(AJ.begin() + AI[neighborIndex], AJ.begin() + AI[neighborIndex + 1], rcell);
                        int index = distance(AJ.begin(), result);
                        
                        for (int irow = 0; irow < nEquations; irow++)
                        {
                            petscrow = (neighborIndex + globalCellIndexShift) * blockSize + irow;
                            for (int n = 0; n < nEquations; n++)
                            {
                                petsccol[n] = rowidx + n;
                            }
                            MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][index * nEquations], ADD_VALUES); // part A

                            PetscInt colidxturb = rowidx + nEquations;
                            MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdqCoupledTerm[irow][index * nTurbulenceEquation], ADD_VALUES); // part B
                        }

                        petscrow = (neighborIndex + globalCellIndexShift) * blockSize + nEquations;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = rowidx + n;
                        }
                        MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdqCoupledTerm_turb[0][index * nEquations], ADD_VALUES); // part C
                        PetscInt colidxturb = rowidx + nEquations;
                        MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdq_turb[0][index * nTurbulenceEquation], ADD_VALUES); // part D
                    }
                }
                rcell++;
            } while (rcell < nTotalCells + nBoundFace);
        }
        else if (interfaceInfo != nullptr && JacOrder != 2)
        {
            vector<int> *neighborCells = grid->GMRESGetNeighborCells();
            rcell = nTotalCells;
            jmax = AI[nTotalCells + nBoundFace];
            j = AI[rcell];
            do
            {
                int rowidx = interfaceInfo->MatchedGlobalNeighborCellIndex(rcell);
                if (rowidx >= 0) 
                {
                    int leftcell = interfaceInfo->MatchedLocalPhysicalCellIndex(rcell); // return its corresponding physical cell index
                    int colidx = (leftcell + globalCellIndexShift) * blockSize;
                    vector<int>::iterator result = find(AJ.begin() + AI[rcell], AJ.begin() + AI[rcell + 1], leftcell);
                    int index = distance(AJ.begin(), result);

                    rowidx *= blockSize;
                    for (int irow = 0; irow < nEquations; irow++)
                    {
                        petscrow = rowidx + irow;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }
                        MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][index * nEquations], ADD_VALUES); // part A
                        PetscInt colidxturb = colidx + nEquations;
                        MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdqCoupledTerm[irow][index * nTurbulenceEquation], ADD_VALUES); // part B
                    }

                    petscrow = rowidx + nEquations;
                    for (int n = 0; n < nEquations; n++)
                    {
                        petsccol[n] = colidx + n;
                    }
                    MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdqCoupledTerm_turb[0][index * nEquations], ADD_VALUES); // part C
                    PetscInt colidxturb = colidx + nEquations;
                    MatSetValues(A, 1, &petscrow, 1, &colidxturb, &dRdq_turb[0][index * nTurbulenceEquation], ADD_VALUES); // part D
                }
                rcell++;
            } while (rcell < nTotalCells + nBoundFace);
        }

        if (JacOrder == 2)
        {
            int jmax = AI1st[nTotalCells];
            int rcell = 0;
            int j = AI1st[rcell];
            do
            {
                int rowidx = (rcell + globalCellIndexShift) * blockSize;
                for (int j = AI1st[rcell]; j < AI1st[rcell + 1]; j++)
                {
                    int qcell = AJ1st[j]; // get qcell (colume index of matrix)
                    // start to insert values: limited to dimension of nTotalSize
                    if (qcell < nTotalCells)
                    {
                        int colidx = (qcell + globalCellIndexShift) * blockSize;
                        for (int irow = 0; irow < nEquations; irow++)
                        {
                            petscrow = rowidx + irow;
                            for (int n = 0; n < nEquations; n++)
                            {
                                petsccol[n] = colidx + n;
                            }
                            MatSetValues(APC, 1, &petscrow, nEquations, petsccol, &dRdq1st[irow][j * nEquations], INSERT_VALUES); // part A
                            PetscInt colidxturb = colidx + nEquations;
                            MatSetValues(APC, 1, &petscrow, 1, &colidxturb, &dRdqCoupledTerm1st[irow][j * nTurbulenceEquation], INSERT_VALUES); // part B
                        }

                        petscrow = rowidx + nEquations;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }
                        MatSetValues(APC, 1, &petscrow, nEquations, petsccol, &dRdqCoupledTerm_turb1st[0][j * nEquations], INSERT_VALUES); // part C
                        PetscInt colidxturb = colidx + nEquations;
                        MatSetValues(APC, 1, &petscrow, 1, &colidxturb, &dRdq_turb1st[0][j * nTurbulenceEquation], INSERT_VALUES); // part D
                    }
                }
                rcell++;
            } while (rcell < nTotalCells);
        }
        else
        {
            int jmax = AI[nTotalCells];
            int rcell = 0;
            int j = AI[rcell];
            do
            {
                int rowidx = (rcell + globalCellIndexShift) * blockSize;
                for (j = AI[rcell]; j < AI[rcell + 1]; j++)
                {
                    int qcell = AJ[j]; // get qcell (colume index of matrix)
                    // start to insert values: limited to dimension of nTotalSize
                    if(qcell < nTotalCells)
                    {
                        int colidx = (qcell + globalCellIndexShift) * blockSize;

                        for (int irow = 0; irow < nEquations; irow++)
                        {
                            petscrow = rowidx + irow;
                            for (int n = 0; n < nEquations; n++)
                            {
                                petsccol[n] = colidx + n;
                            }
                            MatSetValues(APC, 1, &petscrow, nEquations, petsccol, &dRdq[irow][j * nEquations], INSERT_VALUES);  // part A

                            PetscInt colidxturb =  colidx + nEquations;
                            MatSetValues(APC, 1, &petscrow, 1, &colidxturb, &dRdqCoupledTerm[irow][j*nTurbulenceEquation], INSERT_VALUES);    // part B
                        }

                        petscrow = rowidx + nEquations;
                    
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }

                        MatSetValues(APC, 1, &petscrow, nEquations, petsccol, &dRdqCoupledTerm_turb[0][j * nEquations], INSERT_VALUES);  // part C

                        PetscInt colidxturb =  colidx + nEquations;
                        MatSetValues(APC, 1, &petscrow, 1, &colidxturb, &dRdq_turb[0][j*nTurbulenceEquation], INSERT_VALUES);   // part D
                    
                    }
                }
                rcell++; // update rcell, go to the next row index of matrix
            } while (rcell < nTotalCells);
        }
        MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(APC, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(APC, MAT_FINAL_ASSEMBLY);

        // assemble vector
        VecSet(b, zero);
        for (PetscInt icell = 0; icell < nTotalCells; icell++)
        {
            for (PetscInt iequ = 0; iequ < nEquations; iequ++)
            {
                PetscInt isize = (icell + globalCellIndexShift) * blockSize + iequ;
                VecSetValues(b, 1, &isize, &res[iequ][icell], ADD_VALUES);
            }

            PetscInt isize = (icell + globalCellIndexShift) * blockSize + nEquations;
            VecSetValues(b, 1, &isize, &res_turb[0][icell], ADD_VALUES);
        }

        VecAssemblyBegin(b);
        VecAssemblyEnd(b);

        // define solver
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetType(ksp, KSPGMRES);
        KSPSetOperators(ksp, A, APC);

        KSPSetPCSide(ksp, PC_RIGHT);
        KSPSetUp(ksp);
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCBJACOBI);

        // set sub-ksp and sub-pc
        KSP *subksp;
        PC subpc;
        int localTotalSize = nTotalCells * blockSize;
        PCBJacobiSetLocalBlocks(pc, 1, &localTotalSize);

        int nLocal = 1;
        int firstLocal = grid->GetZoneID();
        PCBJacobiGetSubKSP(pc, &nLocal, &firstLocal, &subksp);
        KSPGetPC(subksp[0], &subpc);
        PCSetType(subpc, PCILU);
        PCFactorSetLevels(subpc, nLevels);
        KSPGMRESSetRestart(ksp, MaxStepsRestart);
       
        KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
        KSPSetTolerances(ksp, rTolerance, 1.e-20, dropTolerance, MaxSteps);

        // solution
        KSPSolve(ksp, b, x);

        //! GMRESResidual
        RDouble **dqNS      = dqProxyNS[iZone]->GetField_UNS();
        RDouble **dqTurb    = dqProxyTurb[iZone]->GetField_UNS();

        // convert x to dqProxy
        for (PetscInt icell = 0; icell < nTotalCells; icell++)
        {
            for (PetscInt iequ = 0; iequ < nEquations; iequ++)
            {
                PetscInt isize = (icell + globalCellIndexShift) * blockSize + iequ;
                PetscReal value;
                VecGetValues(x, 1, &isize, &dqNS[iequ][icell]);
            }

            PetscInt isize = (icell + globalCellIndexShift) * blockSize + nEquations;
            VecGetValues(x, 1, &isize, &dqTurb[0][icell]);
        }

        // free work space
        VecDestroy(&x);
        VecDestroy(&b);
        MatDestroy(&A);
        MatDestroy(&APC);
        KSPDestroy(&ksp);
        PetscFinalize();
    
        
    }
    
    // GMRESWholeJacobian
    UpdateFlowField(0, level, dqProxyNS);
    UpdateFlowField(1, level, dqProxyTurb);


    for (int iLocalZone = 0; iLocalZone < nLocalZones; iLocalZone++)
    {
        delete dqProxyNS[iLocalZone];
        delete dqProxyTurb[iLocalZone];
    }
    delete []dqProxyNS;
    delete []dqProxyTurb;
}

//! For GMRES -- Boqian Wang
void Controller::GMRESSolver_backupversion(int iSolver, int level)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();

    FieldProxy **dqProxy = new FieldProxy *[nLocalZones];

    for (int iZone = 0; iZone < nLocalZones; iZone++)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if(!zone)
        {
            continue;
        }

        Grid *gridin = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        UnstructGrid* grid = static_cast<UnstructGrid *>(gridin);
        //! Memory allocating
        RDouble **dRdq = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdq"));
        RDouble **dDdP = reinterpret_cast<RDouble**>(grid->GetDataPtr("dDdP"));
        
        int nTotalCells = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nEquations = solver->GetNumberOfEquations(); 
        
        RDouble **res = reinterpret_cast<RDouble **>(grid->GetDataPtr("res"));
        dqProxy[iZone] = solver->CreateFieldProxy(grid);
        solver->FillField(grid, dqProxy[iZone], zero);  // difference of conservative variables, init by zero
        UnstructGrid* gridunstr = UnstructGridCast(grid);
        RDouble* vol = gridunstr->GetCellVolume();
        RDouble* dt = reinterpret_cast<RDouble*> (grid->GetDataPtr("dt"));
        vector<int> AI = gridunstr->GetJacobianAI4GMRES();  // GMRESSparse GMRESCSR
        vector<int> AJ = gridunstr->GetJacobianAJ4GMRES();  // GMRESSparse GMRESCSR

        // Output res  GMRESWholeJacobian
        string Outputfilename4 = "./results/res1.dat";
        ofstream file4(Outputfilename4.c_str(), ios::trunc);
            
        for (int i = 0; i < nTotalCells; i++)
        {
            for (int j = 0; j < nEquations; j++)
            {
               file4 << res[j][i] << " ";
            }
            file4 << "\n";
        }
        file4.close();


    
        printf("adding diagoal...\n");

        
        // GMRESCSR
        for (int icell = 0; icell < nTotalCells; icell++)
        {
            vector<int>::iterator result = find(AJ.begin() + AI[icell], AJ.begin() + AI[icell + 1], icell); // find qcell
            int index = distance(AJ.begin(), result);
            index *= nEquations;
            for (int iequ = 0; iequ < nEquations; iequ++)
            {
                dRdq[iequ][index + iequ] += 1.0 / dt[icell];
            }
        }

        printf("entry to PETSc\n");
        
        
        // Output dRdq
        string Outputfilename = "./results/nolim_dRdq.dat";
        ofstream file(Outputfilename.c_str(), ios::trunc);
        int total = nTotalCells + nBoundFace;
        int TotalSize = AI[total]*nEquations;
            
        for (int i = 0; i < nEquations; i++)
        {
            for (int j = 0; j < TotalSize; j++)
            {
               file << dRdq[i][j] << " ";
            }
            file << "\n";
        }
        file.close();

       
        // GMRESCSR
        GMRESForZone(dRdq, AI, AJ, res, dqProxy[iZone], nTotalCells, nEquations);
    }
    
    // GMRESWholeJacobian
    UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; iLocalZone++)
    {
        delete dqProxy[iLocalZone];
    }
    delete []dqProxy;
    
}

/// @brief A GMRES solver with PETSc to solve dq for each zone [Ax = b]
/// @param dRdq       Jacobin matrix (A)
/// @param res         residual (b)
/// @param dqProxy     difference of conservative variables (x)
/// @param nTotalCells the real size of the matrix
/// @param nEquations  the number of equations
void Controller::GMRESForZone(RDouble **dRdq, vector<int>& AI, vector<int>& AJ, RDouble **res, FieldProxy *dqProxy, const int& nTotalcells, const int& nEquations)
{
   
    const RDouble LITTLE = 1.0E-14;
    Vec x, b;
    Mat A;
    KSP ksp;                        // define solver
    PC pc;                          // define matrix
    PetscInt nTotalSize;            // dimension size of matrix and vector
    PetscInt nLevels = 3;           // define number of levels used by preconditioner ILU
    PetscInt MaxStepsRestart = 10;  // define maximum steps required by [restart]
    PetscInt MaxSteps = 100;          // define maximum steps
    PetscReal rTolerance = 1.e-3;  // define tolerance of relative error // -1
    PetscReal dropTolerance = 1.e5; // define tolerance of drop error

    nTotalSize = nEquations * nTotalcells;

    // initialize PETSc
    //PetscFunctionBeginUser;
    PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);

    // create vectors;
    VecCreate(PETSC_COMM_SELF, &x);
    PetscObjectSetName((PetscObject)x, "solution");
    VecSetSizes(x, PETSC_DECIDE, nTotalSize);
    VecSetUp(x);
    VecDuplicate(x, &b);

    // create matrix
    // GMRESCSR
    PetscInt nnz[nTotalcells];
    for (int icell = 0; icell < nTotalcells; icell++)
    {
        int count = 0;
        for (int j = AI[icell]; j < AI[icell + 1]; j++)
        {
            if(AJ[j] < nTotalcells)
            {
                count++;
            }
        }
        //count *= nEquations;
        /* for (int n = 0; n < nEquations; n++)
        {
            nnz[icell * nEquations + n] = count;
        } */
        nnz[icell] = count;
    }
    MatCreate(PETSC_COMM_SELF, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nTotalSize, nTotalSize);
    MatSetType(A, MATBAIJ);
    MatSetBlockSize(A, nEquations);
    MatSeqBAIJSetPreallocation(A, nEquations, PETSC_DEFAULT, nnz);
    MatSetUp(A);

    // improved approach to assemble matrix CSR
    PetscInt petsccol[nEquations];
    PetscInt petscrow;
    int jmax = AI[nTotalcells];
    int rcell = 0;
    int j = AI[rcell];
    do
    {
        int rowidx = rcell * nEquations;
        for (j = AI[rcell]; j < AI[rcell + 1]; j++)
        {
            int qcell = AJ[j]; // get qcell (colume index of matrix)
            // start to insert values: limited to dimension of nTotalSize
            if(qcell < nTotalcells)
            {
                int colidx = qcell * nEquations;

                for (int irow = 0; irow < nEquations; irow++)
                {
                    petscrow = rowidx + irow;
                   
                    for (int n = 0; n < nEquations; n++)
                    {
                        petsccol[n] = colidx + n;
                    }
                    MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][j * nEquations], INSERT_VALUES);
                }
            }
        }
        rcell++; // update rcell, go to the next row index of matrix
    } while (j < jmax);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    printf("Finish assemabling matrix...\n");

    // assemble vector
    VecSet(b, zero);
    for (PetscInt icell = 0; icell < nTotalcells; icell++)
    {
        for (PetscInt iequ = 0; iequ < nEquations; iequ++)
        {
            PetscInt isize = icell * nEquations + iequ;
            VecSetValues(b, 1, &isize, &res[iequ][icell], ADD_VALUES);
        }
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    printf("Finish assemabling vector...\n");

    // define solver
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCILU);
    PCFactorSetLevels(pc, nLevels);
    KSPGMRESSetRestart(ksp, MaxStepsRestart);
    KSPGMRESSetOrthogonalization(ksp, KSPGMRESClassicalGramSchmidtOrthogonalization);
    KSPGMRESSetCGSRefinementType(ksp, KSP_GMRES_CGS_REFINE_ALWAYS);
    KSPSetTolerances(ksp, rTolerance, 1.e-20, dropTolerance, MaxSteps);
    KSPMonitorSet(ksp, MyKSPMonitor, NULL, 0);
    // solution
    KSPSolve(ksp, b, x);

    // convert x to dqProxy
    RDouble **dq = dqProxy->GetField_UNS();
    for (PetscInt icell = 0; icell < nTotalcells; icell++)
    {
        for (PetscInt iequ = 0; iequ < nEquations; iequ++)
        {
            PetscInt isize = icell * nEquations + iequ;
            PetscReal value;
            VecGetValues(x, 1, &isize, &dq[iequ][icell]);
        }
    }

    // free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
    PetscFinalize();
    return;
}
#endif

bool Controller::CVGNorm(int iSolver, int level)
{
    int innstep = 0;
    GlobalDataBase::GetData("innstep", &innstep, PHINT, 1);

    RDouble tol_sub_iter = 0.01;
    GlobalDataBase::GetData("tol_sub_iter", &tol_sub_iter, PHDOUBLE, 1);

    int min_sub_iter = 1;
    GlobalDataBase::GetData("min_sub_iter", &min_sub_iter, PHINT, 1);

    int max_sub_iter = 1;
    GlobalDataBase::GetData("max_sub_iter", &max_sub_iter, PHINT, 1);

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    RDouble cvgmax = zero, cvg;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        cvg = solver->UnsteadyConvergence(grid);
        cvgmax = MAX(cvgmax, cvg);
    }

    //cout << "iner_iter = " << iner_iter << " cvg = " << cvg << "\n";

    GlobalSolve();

    PH_CompareMaxMin(cvgmax, 1);

    if (innstep <  min_sub_iter) return false;
    if (cvgmax  <= tol_sub_iter) return true ;
    if (innstep >= max_sub_iter) return true ;
    return false;
}

RDouble Controller::ComputeSubIterationNorm(int iSolver, int level)
{
    RDouble subIterationNorm    =  zero;
    RDouble maxSubIterationNorm = -LARGE;

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        subIterationNorm = solver->UnsteadyConvergence(grid);

        maxSubIterationNorm = MAX(maxSubIterationNorm, subIterationNorm);
    }

    //PH_CompareMaxMin(averageQuotientMAX, 1);

    return maxSubIterationNorm;
}

int Controller::GetNPostSolve(int iSolver)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        return solver->GetNPostSolve();
    }

    return 1;
}

void Controller::PostSolve(int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int iSolver = this->GetSolverIndex();
    int nPostSolve = GetNPostSolve(iSolver);

    //! Face Communication.
    for (int n = 0; n < nPostSolve; ++ n)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
            if (!ifCalTheZone)
            {
                continue;
            }

            Grid *grid = GetGrid(iZone,level);
            CFDSolver *solver = GetCFDSolver(iZone, iSolver);

            solver->PostSolve(grid, n, level);
        }
    }

    //! The following two functions are used for gradient pre-computing for turbulent solver.
    //! Initialization of boundary and dependent variables for 
    //! the next iteration and other coupled solvers, such as turbulence solver.
    SetBoundaryValueOfField(iSolver, level);

    //UpdateIndependentField(iSolver, level);

    //! Post tasks, such as residual computing.
    int zoneIndex = GetLocalZoneIDToGlobalZoneID(0);
    Grid *grid = GetGrid(zoneIndex, level);
    CFDSolver *solver = static_cast<CFDSolver *>(GlobalSolvers::GetSolver(zoneIndex, iSolver));
    solver->DumpResultFile(grid, level);

     //! Exam the change of key variables such as the pressure, the translational temperature and the vibrational temperature, etc.
    int isAdaptiveSolver = GlobalDataBase::GetIntParaFromDB("isAdaptiveSolver");
    if (isAdaptiveSolver > 0)
    {
        solver->CheckResult(grid, level);
    }

    //To exam the residual of key variable.
    int isUseHybridMethod = 0;
    if (GlobalDataBase::IsExist("isUseHybridMethod", PHINT, 1))
    {
        isUseHybridMethod = GlobalDataBase::GetIntParaFromDB("isUseHybridMethod");
    }
    if (isUseHybridMethod > 0)
    {
        solver->CheckResiduals(grid, level);
    }

    //! Dump the surface heating change.
    int nSurfHeatMonitor = 0;
    if (GlobalDataBase::IsExist("nSurfHeatMonitor", PHINT, 1))
    {
        nSurfHeatMonitor = GlobalDataBase::GetIntParaFromDB("nSurfHeatMonitor");
    }
    if (nSurfHeatMonitor > 0)
    {
        solver->CheckSurfaceHeatingChange(grid, level);
    }

    int nLeakageMonitor = 0;
    if (GlobalDataBase::IsExist("nLeakageMonitor", PHINT, 1))
    {
        nLeakageMonitor = GlobalDataBase::GetIntParaFromDB("nLeakageMonitor");
    }

    if (nLeakageMonitor > 0)
    {
        solver->DumpLeakageInformation(grid);
    }

    //! Dump the CFL number to the file.
    int nDumpCFLNumber = 0;
    if (GlobalDataBase::IsExist("nDumpCFLNumber", PHINT, 1))
    {
        nDumpCFLNumber = GlobalDataBase::GetIntParaFromDB("nDumpCFLNumber");
    }
    if (nDumpCFLNumber > 0)
    {
        solver->DumpCFLNumber(grid, level);
    }

    //! Node communication: part1.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        CFDSolver *solver1 = GetCFDSolver(iZone, iSolver);
        for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
        {
            Grid *grid1 = GetGrid(iZone, iLevel);
            solver1->ComputeNodeValue(grid1);
            solver1->ComputeQTurbNodeValue(grid1);
        }
    }

    //! Node communication: part2.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }
        Grid *grid1 = GetGrid(iZone,level);
        CFDSolver *solver1 = GetCFDSolver(iZone, iSolver);
        InterpointInformation *inpterPointInfor = grid1->GetInterpointInfo();
        if (inpterPointInfor)
        {
            solver1->CommunicationInterpointData();
        }
    }

    //! Node communication: part3.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (! zone)
        {
            continue;
        }

        CFDSolver *solver1 = GetCFDSolver(iZone, iSolver);
        for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
        {
            Grid *grid1 = GetGrid(iZone, iLevel);
            solver1->ModifyNodeValue(grid1);
            solver1->ModifyQTurbNodeValue(grid1);
        }
    }
}

void Controller::FillField(int iSolver, int level, FieldProxy **targetProxy, FieldProxy **sourceProxy)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        CFDSolver *solver = GetCFDSolver(zoneID,iSolver);

        solver->FillField(grid, targetProxy[iZone], sourceProxy[zoneID]);
    }
}

void Controller::RungeKutta(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        solver->RungeKutta(grid, rightHandSideProxy[iZone]);
    }
}

void Controller::InitResidual(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        solver->InitResidual(grid, rightHandSideProxy[iZone]);
    }
}

void Controller::CreateQProxy(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **qProxy = GetQProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        qProxy[iZone] = solver->CreateFieldProxy(grid);
    }
}

void Controller::DestroyQProxy(int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **qProxy = GetQProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, 0);
        if (!ifCalTheZone)
        {
            continue;
        }
        delete qProxy[iZone];
        qProxy[iZone] = NULL;
    }
}

void Controller::CreateRHSProxy(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);

#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        rightHandSideProxy[zoneID] = solver->CreateFieldProxy(grid);
    }
}

void Controller::DestroyRHSProxy(int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);
#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, 0);
        if (!ifCalTheZone)
        {
            continue;
        }

        delete rightHandSideProxy[zoneID];
        rightHandSideProxy[zoneID] = nullptr;
    }
}

void Controller::StoreRhsByResidual(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);

#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->StoreRhsByResidual(grid, rightHandSideProxy[zoneID]);
    }
}

void Controller::RecoverResidual(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);
#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(zoneID, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        CFDSolver *solver = GetCFDSolver(zoneID,iSolver);

        solver->RecoverResidual(grid,rightHandSideProxy[zoneID]);
    }
}

bool Controller::JudgeIfCalTheZoneForDensityBasedMethod(int zoneID, int iSolver)
{
    using namespace PHMPI;
    Zone *zone = GetZone(zoneID);
    bool IfCalTheZone = true;
    if (!zone)
    {
        IfCalTheZone = false;
    }
    else
    {
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
        int solverType = solver->GetKey();
        if (PRESSURE_BASED_SOLVER == solverType)
        {
            IfCalTheZone = false;
        }
    }
    return IfCalTheZone;
}

bool Controller::JudgeIfCalTheZoneForPressureBasedMethod(int zoneID, int iSolver)
{
    using namespace PHMPI;
    Zone* zone = GetZone(zoneID);
    bool  IfCalTheZone = true;
    if (!zone)
    {
        IfCalTheZone = false;
    }
    else
    {
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
        int solverType = solver->GetKey();
        if (DENSITY_BASED_SOLVER == solverType)
        {
            IfCalTheZone = false;
        }
    }
    return IfCalTheZone;
}


FieldProxy **Controller::GetQProxy(int level)
{
    return qProxy->GetFieldProxy(level);
}

FieldProxy **Controller::GetRHSProxy(int level)
{
    return rhsProxy->GetFieldProxy(level);
}

void Controller::LoadQ(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **qProxy = GetQProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        solver->LoadQ(grid, qProxy[iZone]);
    }
}

bool Controller::ISCoarsestGrid(int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, 0);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        if (grid->GetCoarseGrid())
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    return false;
}

LIB_EXPORT void Controller::InterpolatFineGrid(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        bool ifCalTheZone = JudgeIfCalTheZoneForDensityBasedMethod(iZone, iSolver);
        if (!ifCalTheZone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        Grid *coarseGrid = grid->GetCoarseGrid();

        solver->InterpolatFineGrid(grid, coarseGrid);
        solver->CommunicationInterfaceData();
    }
}

LIB_EXPORT void Controller::MultiGridInitFlow(int level)
{
    int iSolver = this->GetSolverIndex();

    ZeroResiduals(iSolver, level);

    //! Create Right Hand Side term.
    CreateRHSProxy(iSolver, level);
    StoreRhsByResidual(iSolver, level);

    //! Relaxation.
    Relaxation(iSolver, level);

    //! Recover the residual.
    DestroyRHSProxy(level);
}

MultiGridFieldProxy::MultiGridFieldProxy(int iSolver, int level, int nZone)
{
    this->iSolver = iSolver;
    this->level   = level;
    this->nzone   = nZone;

    field_proxy = new FieldProxy *[nZone];
    for (int iZone = 0; iZone < nZone; ++ iZone)
    {
        field_proxy[iZone] = 0;
    }
}

MultiGridFieldProxy::~MultiGridFieldProxy()
{
    for (int iZone = 0; iZone < nzone; ++ iZone)
    {
        delete field_proxy[iZone];
    }
    delete [] field_proxy;
}

ZoneFieldProxy::ZoneFieldProxy(int iSolver, uint_t nLevel, int nZone)
{
    this->iSolver = iSolver;
    this->nzone   = nZone;
    this->nlevel  = nLevel;

    field_proxy = new MultiGridFieldProxy *[nLevel];

    for (int level = 0; level < nLevel; ++ level)
    {
        field_proxy[level] = new MultiGridFieldProxy(iSolver, level, nZone);
    }
}

ZoneFieldProxy::~ZoneFieldProxy()
{
    for (int level = 0; level < nlevel; ++ level)
    {
        delete field_proxy[level];
    }
    delete [] field_proxy;
}

FieldProxy **ZoneFieldProxy::GetFieldProxy(int level)
{
    return field_proxy[level]->GetFieldProxy();
}

}
