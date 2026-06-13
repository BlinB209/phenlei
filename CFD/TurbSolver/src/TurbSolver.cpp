#include "TurbSolver.h"
#include "Constants.h"
#include "TK_Exit.h"
#include "TK_Log.h"

namespace PHSPACE
{

void Turb_MxDQ(RDouble *mat,RDouble *dq,int numberOfTurbulenceEquation, RDouble *df)
{
    for (int m = 0; m < numberOfTurbulenceEquation; ++ m)
    {
        df[m] = mat[m] * dq[m];
    }
}

TurbSolver::TurbSolver()
{
}

TurbSolver::~TurbSolver()
{
}

string TurbSolver::CastName(const string &name)
{
    if (name == "res")
    {
        return "res_turb";
    }
    else if (name == "q")
    {
        return "q_turb";
    }

    return name;
}

int TurbSolver::GetNumberOfEquations()
{
    Param_TurbSolver *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    return nTurbulenceEquation;
}

void TurbSolver::InitDependentVariables()
{
    Grid *fineGrid = GetGrid(0);
    //! Don't be used for STRUCTGRID because of turb viscosity being stored.
    //if (fineGrid->Type() != STRUCTGRID)
    //{
    //    ComputeViscousCoeff(fineGrid);
    //}

    ComputeQTurbNodeValue(fineGrid);
}

void TurbSolver::InitCoarseGridsFlow()
{
}

void TurbSolver::Action(ActionKey *actkey)
{
    CFDSolver::Action(actkey);
}

void TurbSolver::FillActionKey(ActionKey *actkey, int action, int level)
{
    actkey->action   = action;
    actkey->solver   = TURBULENCE;
    actkey->solverID = this->GetIndex();
    actkey->kind     = SOLVER_BASED;
    actkey->level    = level;
}

void TurbSolver::Solve()
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
        FreePointer(actkeyDumpRestartData);
    }

    int intervalStepRes;
    GlobalDataBase::GetData("intervalStepRes", &intervalStepRes, PHINT, 1);
    if (outnstep % intervalStepRes == 0)
    {
        //! To dump residual data.
        ActionKey *actkeyDumpResidual = new ActionKey();
        FillActionKey(actkeyDumpResidual, DUMP_RESIDUAL, 0);
        DumpResidual(actkeyDumpResidual);
        FreePointer(actkeyDumpResidual);
    }
}

void TurbSolver::PostSolve(Grid *grid_in, int stage, int level)
{
    CommunicationInterfaceData();
}

void TurbSolver::DumpResultFile(Grid *grid, int level)
{
    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int intervalStepFlow = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");

    bool isSubIterationDump = true;    //! It is true default in steady simulations.

    Param_TurbSolver *parameters = GetControlParameters();
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
            FreePointer(actkeyDumpRestartData);
        }
        else
        {
            int RestartFile = 0;
            ActionKey *actkeyDumpRestartData = new ActionKey();
            FillActionKey(actkeyDumpRestartData, DUMP_RESTART, level);
            DumpTurbProtectData(actkeyDumpRestartData, RestartFile);
            FreePointer(actkeyDumpRestartData);

            if (outnstep > 1)
            {
                if (RestartFile == 0)
                {
                    string path = "./results/turb1.dat";
                    remove(path.c_str());

                    if (PHMPI::IsParallelRun())
                    {
                        string path1 = "./results/turb1_0.dat";
                        remove(path1.c_str());
                    }
                }

                else if (RestartFile == 1)
                {
                    string path = "./results/turb0.dat";
                    remove(path.c_str());

                    if (PHMPI::IsParallelRun())
                    {
                        string path2 = "./results/turb0_0.dat";
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
        ActionKey *actkeyDumpResidual = new ActionKey();
        FillActionKey(actkeyDumpResidual, DUMP_RESIDUAL, level);
        DumpResidual(actkeyDumpResidual);
        FreePointer(actkeyDumpResidual);
    }
}

void TurbSolver::Post()
{
    //! To dump restart data for continual simulation.
    ActionKey *actkeyDumpRestartData = new ActionKey();
    FillActionKey(actkeyDumpRestartData, DUMP_RESTART, 0);
    DumpRestartData(actkeyDumpRestartData);
    FreePointer(actkeyDumpRestartData);
}

void TurbSolver::RegisterCFDSolverInterfaceField()
{
    Param_TurbSolver *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    
    CFDSolver::RegisterInterfaceField("turb::q", PHDOUBLE, nTurbulenceEquation);

    CFDSolver::RegisterInterfaceField("turb::dq", PHDOUBLE, nTurbulenceEquation);

    CFDSolver::RegisterInterfaceField("turb::dqdx", PHDOUBLE, nTurbulenceEquation);
    CFDSolver::RegisterInterfaceField("turb::dqdy", PHDOUBLE, nTurbulenceEquation);
    CFDSolver::RegisterInterfaceField("turb::dqdz", PHDOUBLE, nTurbulenceEquation);

    CFDSolver::RegisterInterfaceField("vist", PHDOUBLE, 1);
    CFDSolver::RegisterInterfaceField("blend", PHDOUBLE, 1);
}

void TurbSolver::RegisterCFDSolverInterpointField()
{
    Param_TurbSolver *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    CFDSolver::RegisterInterpointField("qTurbNode", PHDOUBLE, nTurbulenceEquation);
}

void TurbSolver::RegisterOversetField()
{
    Param_TurbSolver *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    CFDSolver::RegisterOversetField("turb::q", PHDOUBLE, nTurbulenceEquation);
    CFDSolver::RegisterOversetField("vist", PHDOUBLE, 1);

    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");
    int isCommGradientforSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isCommGradientforSlip");
    if (isOversetSlip && isCommGradientforSlip)
    {
        CFDSolver::RegisterOversetField("turb::dqdx", PHDOUBLE, nTurbulenceEquation);
        CFDSolver::RegisterOversetField("turb::dqdy", PHDOUBLE, nTurbulenceEquation);
        CFDSolver::RegisterOversetField("turb::dqdz", PHDOUBLE, nTurbulenceEquation);
    }
}

void TurbSolver::DeAllocateOversetInterfaceVar(Data_ParamFieldSuite *dataStore)
{
    RDouble **q = reinterpret_cast <RDouble **> (dataStore->GetDataPtr("turb::q"));
    DelPointer2(q);

    RDouble **vist = reinterpret_cast <RDouble **> (dataStore->GetDataPtr("vist"));
    DelPointer2(vist);

    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");
    int isCommGradientforSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isCommGradientforSlip");
    if (isOversetSlip && isCommGradientforSlip)
    {
        RDouble **gradientTurbulenceX = reinterpret_cast <RDouble **> (dataStore->GetDataPtr("turb::dqdx"));
        RDouble **gradientTurbulenceY = reinterpret_cast <RDouble **> (dataStore->GetDataPtr("turb::dqdy"));
        RDouble **gradientTurbulenceZ = reinterpret_cast <RDouble **> (dataStore->GetDataPtr("turb::dqdz"));

        DelPointer2(gradientTurbulenceX);
        DelPointer2(gradientTurbulenceY);
        DelPointer2(gradientTurbulenceZ);
    }
}

void TurbSolver::PrepareOversetInterfaceData(Data_ParamFieldSuite *datastore, InterfaceDataProxy *interfaceDataProxy)
{
    Param_TurbSolver *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (datastore->GetDataPtr("turb::q"));
    RDouble **viscousTurbulence = reinterpret_cast<RDouble **> (datastore->GetDataPtr("vist"));

    vector<RDouble **> &vectorData = interfaceDataProxy->GetVectorData();
    vector<int>   &vectorDimension = interfaceDataProxy->GetVectorDimension();

    int isSolve = GlobalDataBase::GetIntParaFromDB("isSolve");
    if (!isSolve)
    {
    vectorDimension.push_back(nTurbulenceEquation);
    vectorData.push_back(qTurbulence);

        vectorDimension.push_back(1);
        vectorData.push_back(viscousTurbulence);
    }
    else
    {
        RDouble **gradientTurbulenceX = reinterpret_cast <RDouble**> (datastore->GetDataPtr("turb::dqdx"));
        RDouble **gradientTurbulenceY = reinterpret_cast <RDouble**> (datastore->GetDataPtr("turb::dqdy"));
        RDouble **gradientTurbulenceZ = reinterpret_cast <RDouble**> (datastore->GetDataPtr("turb::dqdz"));
        vectorDimension.push_back(nTurbulenceEquation);
        vectorData.push_back(gradientTurbulenceX);
        vectorDimension.push_back(nTurbulenceEquation);
        vectorData.push_back(gradientTurbulenceY);
        vectorDimension.push_back(nTurbulenceEquation);
        vectorData.push_back(gradientTurbulenceZ);
    }
}

//! Drive functions to calculate residuals for different orders and dimensions.
void TurbSolver::UpdateResiduals(Grid *grid)
{
    //! Set values for flow variables in ghost cells.
    Boundary(grid);

    //! Modify the variable value in the cell center with wall boundary face.
    //ResetWallScalar(grid);

    //! Zero the other variables of the grid.
    InitSpectrum(grid);

    //! Compute the source term.
    SourceFlux(grid);

    //! Compute the spectrum radius.
    //Diagonal(grid);

    //! Compute the inviscid flux.
    InviscidFlux(grid);

    //! Compute the viscous flux.
    ViscousFlux(grid);

    DualTimeSource(grid);

    ZeroResidualOfSpecialCells(grid);    //! Added By Guo Yongheng 20160825.

    ModifyResiduals(grid);

    FreeGradientProxy(grid);
}

void TurbSolver::SourceFlux(Grid *grid)
{
    Param_TurbSolver *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    if (viscousType == ONE_EQU)
    {
        //! One equation turbulence model.
        SourceFluxOneEquation(grid);
    }
    else if (viscousType == TWO_EQU)
    {
        //! Two equation turbulence model.
        SourceFluxTwoEquation(grid);
    }
}

//! Compare the outnstep in the flowfield file between the NS and turbulence.
void TurbSolver::CompareOutStepOfFlowfieldFile(int outnstepofNS, int outnstepofTurb) const
{
    if (outnstepofNS == outnstepofTurb) return;

    ostringstream oss;
    oss << "outnstepofNS   = " << outnstepofNS   << endl;
    oss << "outnstepofTurb = " << outnstepofTurb << endl;
    oss << "The NS flowfield file is not matched with the turbulence." << endl;

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (!isUnsteady)
    {
        oss << "Warning: Cluster writes files error! It needs more iteration step." << endl;
        PrintToWindow(oss.str());
        WriteLogFile (oss.str());

        GlobalDataBase::UpdateData("outnstep", &outnstepofTurb, PHINT, 1);
    }
    else
    {
        TK_Exit::ExceptionExit("Error: Cluster writes files error! Please restart the calculation. \n");
    }
}

//! Get the file name for Residual dumping.
const string TurbSolver::GetResidualFileName()
{
    string turbresFile = "turbres.dat";
    GlobalDataBase::GetData("turbresfile", &turbresFile, PHSTRING, 1);
    return turbresFile;
}

//! Get the file name for relay calculation.
const string TurbSolver::GetRestartFileName()
{
    string turbFile = "turb.dat";
    GlobalDataBase::GetData("turbfile", &turbFile, PHSTRING, 1);
    return turbFile;
}

LIB_EXPORT Param_TurbSolver *TurbSolver::GetControlParameters()
{
    return static_cast <Param_TurbSolver *> (controlParameters);
}

}
