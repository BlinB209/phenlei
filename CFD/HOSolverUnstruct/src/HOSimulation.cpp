#include <cstdlib>
#include "PHMpi.h"
#include "PHHeader.h"
#include "Region.h"
#include "Zone.h"
//#include "AleModel.h"
#include "Glb_Dimension.h"
#include "Constants.h"
#include "TK_Log.h"
#include "Controller.h"
#include "TK_Exit.h"
#include "TK_Time.h"
#include "Geo_SimpleBC.h"
//#include "UnstructuredOversetConfig.h"
#include "OversetInformation.h"
#include "Post_WriteVisualFile.h"
#include "Task.h"
#include "PostProcess.h"
#include "TK_Exit.h"
#include "Math_BasisFunction.h"
#include "HOSimulation.h"
#include "Geometry.h"
#include "GridType.h"
#include "HOSolverUnstruct.h"
#include "TurbSolverUnstr.h"

using namespace std;

namespace HOUnstruct
{
using namespace PHSPACE;

HOSimulation::HOSimulation(int dimension)
{
    using namespace PHMPI;
    //InitFantasyMPI();

    SetDim(dimension);

    region = new PHSPACE::Region();

    isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
}

HOSimulation::~HOSimulation()
{
    using namespace PHMPI;

    delete region;

    TK_Exit::ExitPHengLEI();
}

void HOSimulation::InitializeProtoTable()
{

}

void HOSimulation::Run()
{
    ReadGeometry();

    UpdateAllZonesPara();

    InitGeometry();
    
    SwapGeometry();

    InitOversetGrid();

    if (PHMPI::CurrentProcessorIsGridProcessor())
    {
        return;
    }

    InitializeProtoTable();
     
    //! Define solver names for all zones and create the newly defined solvers.
    //! IMPORTANT for solver developers:
    //! Define your solver name which is same with the solver class name here! 
    DefineAndCreateSolvers();

    //! Initialize each solver on each zone.
    InitializeSolvers();

    InitGlobalBoundaryCondition();

    MultiGridInitFlow();
    
    SolveSimu();

    PostSimu();
}

void HOSimulation::SolveSimu()
{
    WriteLogFile("In SolveSimu, Start iterating ...");

    int outnstep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
    int maxSimuStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxSimuStep");

    while (true)
    {
        UpdateOuterIterationStep(outnstep);

        SolveOneOuterStep();

        if (outnstep >= maxSimuStep && maxSimuStep > 0) break;
    }
}

void HOSimulation::UpdateOuterIterationStep(int &outnstep)
{
    outnstep = outnstep + 1;
    GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

    int newnstep = GlobalDataBase::GetIntParaFromDB("newnstep");
    newnstep ++;
    GlobalDataBase::UpdateData("newnstep", &newnstep, PHINT, 1);

    //int isUnsteady = 0;
    //GlobalDataBase::GetData("iunsteady", &isUnsteady, PHINT, 1);
    if (isUnsteady == HO_STEADY) return;

    int ifStaticsFlowField = GlobalDataBase::GetIntParaFromDB("ifStaticsFlowField");
    if (ifStaticsFlowField == 1)
    {
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
        if (outnstep >= startStatisticStep)
        {
            int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
            ++ nStatisticalStep;
            GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
        }
    }
}

void HOSimulation::SolveOneOuterStep()
{
    //PreSolve();

    MainSolve();

    PostSolve();
}

void HOSimulation::PostSolve()
{
    int intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
    int outnstep  = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (outnstep % intervalStepPlot == 0)
    {
        WriteVisualFile();
    }
}

void HOSimulation::MainSolve()
{
    //int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

    if (isUnsteady == HO_STEADY)
    {
        SolveSteady();
    }
    else
    {
        SolveUnsteady();
    }
}

void HOSimulation::PostSimu()
{
    //this->region->PostSimu();
}

void HOSimulation::ComputeWallDistance()
{
    ReadGeometry();

    UpdateAllZonesPara();

    InitGeometry();
}

void HOSimulation::SwapGeometry()
{
    this->region->SwapGeometry();
}

void HOSimulation::InitGeometry()
{
    this->region->InitGeometry();
}

void HOSimulation::ReadGeometry()
{
    PHSPACE::TimeSpan testTime;
    ReadGrid();
    testTime.ShowTimeSpanToLogFile("Grid Reading");
}

void HOSimulation::InitializeSolvers()
{
    WriteLogFile("Creating and initialization solvers start\n");

    int nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");
    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        Controller * controller = new Controller(iSolver, this->region);

        controller->Initialize(iSolver);

        delete controller;
    }

    WriteLogFile("Creating and initialization solvers over\n");
}

void HOSimulation::DefineAndCreateSolvers()
{
    InitGlobalValuesOfAllZones();

    int maxNumOfSolver = 0;
    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        PHSPACE::Zone * zone = this->region->GetZone(iZone);
        if (!zone) continue;

        //! For each zone, define the solver names on it.
        //! If you are solver developer, define your solve name here,
        //! according your solver logic.

        vector<PHSolver *> solverList;

        PHGeometry *geometry = zone->GetGeometry();

        int gridType = geometry->GetGrid()->Type();
        if (gridType != UNSTRUCTGRID)
        {
            TK_Exit::ExceptionExit("HOSolverUnstruct support unstr grid only!");
        }

        PHSolver *solver = new HOSolverUnstruct();
        solver->SetName("HOSolverUnstruct");
        solver->SetKey(DENSITY_BASED_SOLVER);
        solverList.push_back(solver);

        int viscousType = LAMINAR;
        zone->GetData("viscousType",&viscousType, PHINT, 1);

        if (viscousType > LAMINAR)
        {
            PHSolver *solverT = new TurbSolverUnstr();
            solverT->SetName("TurbSolverUnstr");
            solverT->SetKey(DENSITY_BASED_SOLVER);
            solverList.push_back(solverT);
        }

        for (std::size_t iSolver = 0; iSolver < solverList.size(); ++ iSolver)
        {
            PHSolver *solveri = solverList[iSolver];

            solveri->SetGeometry(geometry);
            solveri->SetIndex(static_cast<int> (iSolver));
            zone->AddSolver(solveri);
        }

        if (maxNumOfSolver < solverList.size())
        {
            maxNumOfSolver = solverList.size();
        }
    }

    GlobalDataBase::UpdateData("maxNumOfSolver", &maxNumOfSolver, PHINT, 1);
}

void HOSimulation::PreSolve()
{
    //FYAleModel::UnsteadyPreSolve();
}

void HOSimulation::ReadGrid()
{
    region->ReadGrid();
}

void HOSimulation::UpdateAllZonesPara()
{
    this->region->UpdateAllZonesPara();
}

void HOSimulation::InitGlobalBoundaryCondition()
{
    this->region->ReSetBoundaryConditionByGlobalBC();
}

void HOSimulation::MultiGridInitFlow()
{
    WriteLogFile("Multi-Grid Init Flow ...");

    int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");
    if (nMGLevel <= 1) return;

    int outnstep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
    if (outnstep > 0)
    {
//      int newnstep = outnstep;
//      GlobalDataBase::UpdateData("newnstep",&newnstep, PHINT, 1);
        return; //! Read from flow file.
    }

    int mgInitStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("mgInitStep");
    if (mgInitStep <= 0) return;

    int mgvist = GlobalDataBase::GetIntParaFromDB("mgvist");
    int nSolver = 1;    //! Consider NS solver only.
    if (mgvist) 
    {
        nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");
    }

    int isInMGIniting = 1;
    GlobalDataBase::UpdateData("isInMGIniting", &isInMGIniting, PHINT, 1);

    for (int level = nMGLevel - 1; level > 0; -- level)
    {
        int levelStep = 0;
        int totalLevelStep = mgInitStep * level;

        ostringstream oss;
        oss << "Start init flow on level " << level << ", " << totalLevelStep << " iterations ...";
        WriteLogFile(oss.str()); oss << endl;
        PrintToWindow(oss.str());

        while (levelStep < totalLevelStep)
        {
            ++ levelStep;

            UpdateOuterIterationStep(outnstep);

            //using namespace PHMPI;
            for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
            {
                Controller *controller = new Controller(iSolver, this->region);
                controller->MultiGridInitFlow(level);
                controller->PostSolve(level);
                delete controller;
            }
        }

        //! Interpolate from coarse grid to fine grid.
        for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
        {
            Controller *controller = new Controller(iSolver, this->region);
            controller->InterpolatFineGrid(iSolver, level - 1);
            delete controller;
        }
    }

    isInMGIniting = 0;
    GlobalDataBase::UpdateData("isInMGIniting", &isInMGIniting, PHINT, 1);

    int newnstep = 0;
    GlobalDataBase::UpdateData("newnstep",&newnstep, PHINT, 1);

    PrintToWindow("End init flow by coarse grid ...\n");
    WriteLogFile("End init flow by coarse grid ...");
}

void HOSimulation::SolveSteady()
{
    int nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");
    //! Get iteration number.
    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int turbInterval = GlobalDataBase::GetIntParaFromDB("turbInterval");

    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        //! iSolver = 1 or 2, is supposed as turbulence computing.
        if ((iSolver != NS_SOLVER) && (outnstep % turbInterval != 0))
        {
            return;
        }

        Controller *controller = new Controller(iSolver, this->region);

        controller->SolveSteadyField();

        controller->PostSolve();

        delete controller;
    }
}

void HOSimulation::SolveSteadyField()
{
    //MultiGrid(0);
}

void HOSimulation::SolveUnsteady()
{
    int     min_sub_iter = GlobalDataBase::GetIntParaFromDB("min_sub_iter");
    int     max_sub_iter = GlobalDataBase::GetIntParaFromDB("max_sub_iter");
    RDouble tol_sub_iter = GlobalDataBase::GetDoubleParaFromDB("tol_sub_iter");

    int innstep = 0;
    int nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");

    bool isSubIterationDump = false;
    GlobalDataBase::UpdateData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);

    bool isSubIterationConvergence = false;
    while(!isSubIterationConvergence)
    {
        innstep += 1;
        GlobalDataBase::UpdateData("innstep", &innstep, PHINT, 1);

        RDouble subIterationNormOneSolver     =  zero;
        RDouble maxSubIterationNormAllSolvers = -LARGE;

        for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
        {
            Controller * controller = new Controller(iSolver, this->region);

            controller->SolveSteadyField();

            //! Here is in the sub-iteration, PostSolve is only for the CommunicationInterfaceData and dumps the "res.dat", which contains index of the sub-iteration residual.
            controller->PostSolve();

            subIterationNormOneSolver     = controller->ComputeSubIterationNorm(iSolver, 0);
            maxSubIterationNormAllSolvers = MAX(maxSubIterationNormAllSolvers, subIterationNormOneSolver);

            delete controller;
        }

        PH_CompareMaxMin(maxSubIterationNormAllSolvers, 1);

        //! Sub-iteration continue or not judgement.
        isSubIterationConvergence = IsSubIterationConvergence(innstep, min_sub_iter, max_sub_iter, tol_sub_iter, maxSubIterationNormAllSolvers);
    }

    //! Here is already (innstep == max_sub_iter || isSubIterationConvergence).
    isSubIterationDump = true;
    GlobalDataBase::UpdateData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);

    //! Update physical time step flowfield after sub-iteration convergence.
    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        Controller * controller = new Controller(iSolver, this->region);

        controller->UpdateUnsteadyFlow(iSolver, 0);

        //! Here is out of the sub-iteration, PostSolve is to dump various files according to the parameters in the "cfd_para.hypara", such as intervalStepFlow, intervalStepPlot, intervalStepForce, intervalStepRes, et al.
        controller->PostSolve();

        delete controller;
    }
}

bool HOSimulation::IsSubIterationConvergence(const int & innstep, const int & min_sub_iter, const int & max_sub_iter, const RDouble & tol_sub_iter, const RDouble & averageQuotientMaxInSolvers)
{
    bool isSubIterationConvergence = false;

    if (innstep < min_sub_iter)
    {
        isSubIterationConvergence = false;
    }
    else if (averageQuotientMaxInSolvers <= tol_sub_iter)
    {
        isSubIterationConvergence = true;
    }
    else if (innstep >= max_sub_iter)
    {
        isSubIterationConvergence = true;
    }
    else
    {
        isSubIterationConvergence = false;
    }

    return isSubIterationConvergence;
}

}
