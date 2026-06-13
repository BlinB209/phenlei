#include "TransitionSolver.h"
#include "Constants.h"
#include "TK_Exit.h"
#include "TK_Log.h"

namespace PHSPACE
{

void Transition_MxDQ(RDouble *mat,RDouble *dq,int numberOfTransitionEquation, RDouble *df)
{
    for (int m = 0; m < numberOfTransitionEquation; ++ m)
    {
        df[m] = mat[m] * dq[m];
    }
}

TransitionSolver::TransitionSolver()
{
}

TransitionSolver::~TransitionSolver()
{
}

string TransitionSolver::CastName(const string &name)
{
    if (name == "res")
    {
        return "res_transition";
    }
    else if (name == "q")
    {
        return "q_transition";
    }

    return name;
}

int TransitionSolver::GetNumberOfEquations()
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    return nTransitionEquation;
}

void TransitionSolver::InitDependentVariables()
{
    Grid *fineGrid = GetGrid(0);

    ComputeQTransitionNodeValue(fineGrid);
}

void TransitionSolver::InitCoarseGridsFlow()
{
}

void TransitionSolver::Action(ActionKey *actkey)
{
    CFDSolver::Action(actkey);
}

void TransitionSolver::FillActionKey(ActionKey *actkey, int action, int level)
{
    actkey->action   = action;
    actkey->solver   = TRANSITION;
    actkey->solverID = this->GetIndex();
    actkey->kind     = SOLVER_BASED;
    actkey->level    = level;
}

void TransitionSolver::Solve()
{
    int outnstep = 0;
    //! Get iteration and zonal number.
    GlobalDataBase::GetData("outnstep", &outnstep, PHINT, 1);

    Grid *grid = GetGrid(0);

    //! Zero the residuals of the finest grid.
    ZeroResiduals(grid);

    //! Solve NS for all grids starting from the finest grid.
    SolveOnGrid(grid);

    CommunicationInterfaceData();

    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation)
    {
        CommunicationInterpointData();
    }
    GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

    int intervalStepFlow = 100;
    GlobalDataBase::GetData("intervalStepFlow", &intervalStepFlow, PHINT, 1);

    if (outnstep % intervalStepFlow == 0)
    {
        //! To dump restart data for continual simulation.
        ActionKey *actkeyDumpRestartData = new ActionKey();
        FillActionKey(actkeyDumpRestartData, DUMP_RESTART, 0);
        DumpRestartData(actkeyDumpRestartData);
        delete actkeyDumpRestartData;
    }

    int intervalStepRes;
    GlobalDataBase::GetData("intervalStepRes", &intervalStepRes, PHINT, 1);
    if (outnstep % intervalStepRes == 0)
    {
        //! To dump residual data.
        ActionKey *actkeyDumpResidual = new ActionKey();
        FillActionKey(actkeyDumpResidual, DUMP_RESIDUAL, 0);
        DumpResidual(actkeyDumpResidual);
        delete actkeyDumpResidual;
    }
}

void TransitionSolver::PostSolve(Grid *grid_in, int stage, int level)
{
    CommunicationInterfaceData();
}

void TransitionSolver::DumpResultFile(Grid *grid, int level)
{
    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int intervalStepFlow = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");

    bool isSubIterationDump = true; //! It is true default in steady simulations.

    Param_TransitionSolver *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        GlobalDataBase::GetData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);
    }

    int nProtectData = 0;
    if (GlobalDataBase::IsExist("nProtectData", PHINT, 1))
    {
        nProtectData = GlobalDataBase::GetIntParaFromDB("nProtectData");
    }

    if (outnstep % intervalStepFlow == 0 && level == 0 && isSubIterationDump)
    {

        if (nProtectData ==0)
        {
            //! To dump restart data for continual simulation.
            ActionKey *actkeyDumpRestartData = new ActionKey();
            FillActionKey(actkeyDumpRestartData, DUMP_RESTART, 0);
            DumpRestartData(actkeyDumpRestartData);
            delete actkeyDumpRestartData;
        }
        else
        {
            int RestartFile = 0;
            ActionKey *actkeyDumpRestartData = new ActionKey();
            FillActionKey(actkeyDumpRestartData, DUMP_RESTART, level);
            DumpTransitionProtectData(actkeyDumpRestartData, RestartFile);
            delete actkeyDumpRestartData;
        
            if (outnstep > 1)
            {
                if (RestartFile == 0)
                {
                    string path = "./results/transition1.dat";
                    remove(path.c_str());

                    if (PHMPI::IsParallelRun())
                    {
                        string path1 = "./results/transition1_0.dat";
                        remove(path1.c_str());
                    }
                }

                else if (RestartFile == 1)
                {
                    string path = "./results/transition0.dat";
                    remove(path.c_str());

                    if (PHMPI::IsParallelRun())
                    {
                        string path2 = "./results/transition0_0.dat";
                        remove(path2.c_str());
                    }
                }
            }
        }
    }

    int intervalStepRes = GlobalDataBase::GetIntParaFromDB("intervalStepRes");
    if (outnstep == 1 || outnstep % intervalStepRes == 0 && (!isUnsteady || !isSubIterationDump))
    {
        //! To dump residual data.
        ActionKey * actkeyDumpResidual = new ActionKey();
        FillActionKey(actkeyDumpResidual, DUMP_RESIDUAL, level);
        DumpResidual(actkeyDumpResidual);
        delete actkeyDumpResidual;
    }
}

void TransitionSolver::Post()
{
    //! To dump restart data for continual simulation.
    ActionKey *actkeyDumpRestartData = new ActionKey();
    FillActionKey(actkeyDumpRestartData, DUMP_RESTART, 0);
    DumpRestartData(actkeyDumpRestartData);
    delete actkeyDumpRestartData;
}

void TransitionSolver::RegisterCFDSolverInterfaceField()
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    
    CFDSolver::RegisterInterfaceField("transition::q", PHDOUBLE, nTransitionEquation);

    CFDSolver::RegisterInterfaceField("transition::dq", PHDOUBLE, nTransitionEquation);

    CFDSolver::RegisterInterfaceField("transition::dqdx", PHDOUBLE, nTransitionEquation);
    CFDSolver::RegisterInterfaceField("transition::dqdy", PHDOUBLE, nTransitionEquation);
    CFDSolver::RegisterInterfaceField("transition::dqdz", PHDOUBLE, nTransitionEquation);
}

void TransitionSolver::RegisterCFDSolverInterpointField()
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    CFDSolver::RegisterInterpointField("qTransitionNode", PHDOUBLE, nTransitionEquation);
}

void TransitionSolver::RegisterOversetField()
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    CFDSolver::RegisterOversetField("transition::q", PHDOUBLE, nTransitionEquation);
}

void TransitionSolver::DeAllocateOversetInterfaceVar(Data_ParamFieldSuite *dataStore)
{
    RDouble **q = reinterpret_cast <RDouble **> (dataStore->GetDataPtr("transition::q"));
    DelPointer2(q);
}

void TransitionSolver::PrepareOversetInterfaceData(Data_ParamFieldSuite *datastore, InterfaceDataProxy *interfaceDataProxy)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **qTransition = reinterpret_cast<RDouble **> (datastore->GetDataPtr("transition::q"));

    vector<RDouble **> &vectorData = interfaceDataProxy->GetVectorData();
    vector<int>   &vectorDimension = interfaceDataProxy->GetVectorDimension();

    vectorDimension.push_back(nTransitionEquation);
    vectorData.push_back(qTransition);
}

//! Drive functions to calculate residuals for different orders and dimensions.
void TransitionSolver::UpdateResiduals(Grid *grid)
{
    //! Set values for flow variables in ghost cells.
    Boundary(grid);

    //! Zero the other variables of the grid.
    InitSpectrum(grid);

    //! Compute the source term.
    SourceFlux(grid);

    //! Compute the spectrum radius.
    Diagonal(grid);

    //! Compute the inviscid flux.
    InviscidFlux(grid);

    //! Compute the viscous flux.
    ViscousFlux(grid);

    DualTimeSource(grid);

    ZeroResidualOfSpecialCells(grid); 

    FreeGradientProxy(grid);
}

void TransitionSolver::SourceFlux(Grid *grid)
{
    SourceFluxTwoEquation(grid);
}

//! Compare the outnstep in the flowfield file between the NS and transition.
void TransitionSolver::CompareOutStepOfFlowfieldFile(int outnstepofNS, int outnstepofTransition) const
{
    if (outnstepofNS == outnstepofTransition) return;

    ostringstream oss;
    oss << "outnstepofNS   = " << outnstepofNS   << endl;
    oss << "outnstepofTransition = " << outnstepofTransition << endl;
    oss << "The NS flowfield file is not matched with the transition." << endl;

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (!isUnsteady)
    {
        oss << "Warning: Cluster writes files error! It needs more iteration step." << endl;
        PrintToWindow(oss.str());
        WriteLogFile (oss.str());

        GlobalDataBase::UpdateData("outnstep", &outnstepofTransition, PHINT, 1);
    }
    else
    {
        TK_Exit::ExceptionExit("Error: Cluster writes files error! Please restart the calculation. \n");
    }
}

//! Get the file name for Residual dumping.
const string TransitionSolver::GetResidualFileName()
{
    string transitionResFile = "transitionRes.dat";
    GlobalDataBase::GetData("transitionResFile", &transitionResFile, PHSTRING, 1);
    return transitionResFile;
}

//! Get the file name for relay calculation.
const string TransitionSolver::GetRestartFileName()
{
    string transitionFile = "transition.dat";
    GlobalDataBase::GetData("transitionFile", &transitionFile, PHSTRING, 1);
    return transitionFile;
}

LIB_EXPORT Param_TransitionSolver *TransitionSolver::GetControlParameters()
{
    return static_cast <Param_TransitionSolver *> (controlParameters);
}

}
