#include <cmath>
#include "TransitionSolverStruct.h"
#include "OversetGridFactory.h"
#include "Geo_StructBC.h"
#include "GradientOperation.h"
#include "Math_Limiter.h"
#include "TK_Exit.h"
#include "Param_TransitionSolver.h"
#include "Residual.h"
#include "FieldProxy.h"
#include "Glb_Dimension.h"
#include "TK_Log.h"
#include "IO_FileName.h"
#include "Transition.h"
#include <cstdio>
#include <stdio.h>
#include "TK_Time.h"


using namespace std;

namespace PHSPACE
{
TransitionSolverStr::TransitionSolverStr()
{
}

TransitionSolverStr::~TransitionSolverStr()
{
    DeAllocateGlobalVariables();

    FreeControlParameters();
}

void TransitionSolverStr::AllocateGlobalVar(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    if (nTransitionEquation == 0)
    {
        return;
    }

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range IFACE, JFACE, KFACE;
    GetRange(ni, nj, nk, -2, 2, IFACE, JFACE, KFACE);

    Range M(0, nTransitionEquation - 1);

    RDouble4D *qTransition = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *residualTransition = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *spectrumTransition = new RDouble4D(I, J, K, M, fortranArray);

    grid->UpdateDataPtr("q_transition", qTransition);    //! Unknown variable of transition model.
    grid->UpdateDataPtr("res_transition", residualTransition);    //! Residual or right-hand side.
    grid->UpdateDataPtr("spec_transition", spectrumTransition);    //! Spectral radius of transition.

    RDouble4D *gradTransitionX = new RDouble4D(IFACE, JFACE, KFACE, M, fortranArray);
    RDouble4D *gradTransitionY = new RDouble4D(IFACE, JFACE, KFACE, M, fortranArray);
    RDouble4D *gradTransitionZ = new RDouble4D(IFACE, JFACE, KFACE, M, fortranArray);
    grid->UpdateDataPtr("gradTransitionX", gradTransitionX);    //! Gradient of unknown variable at face, for x direction. in x-face.
    grid->UpdateDataPtr("gradTransitionY", gradTransitionY);    //! Gradient of unknown variable at face, for y direction. in y-face.
    grid->UpdateDataPtr("gradTransitionZ", gradTransitionZ);    //! Gradient of unknown variable at face, for z direction. in z-face.

    RDouble4D *gradTransitionCellCenterX = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *gradTransitionCellCenterY = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *gradTransitionCellCenterZ = new RDouble4D(I, J, K, M, fortranArray);
    grid->UpdateDataPtr("gradTransitionCellCenterX", gradTransitionCellCenterX);    //! Gradient of unknown variable at cell-center in x-direction.
    grid->UpdateDataPtr("gradTransitionCellCenterY", gradTransitionCellCenterY);    //! Gradient of unknown variable at cell-center in y-direction.
    grid->UpdateDataPtr("gradTransitionCellCenterZ", gradTransitionCellCenterZ);    //! Gradient of unknown variable at cell-center in z-direction.

    //! Temporary variable for periodic boundary condition.
    RDouble4D *rotTransitiongradValueX = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *rotTransitiongradValueY = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *rotTransitiongradValueZ = new RDouble4D(I, J, K, M, fortranArray);
    grid->UpdateDataPtr("rotTransitiongradValueX", rotTransitiongradValueX);
    grid->UpdateDataPtr("rotTransitiongradValueY", rotTransitiongradValueY);
    grid->UpdateDataPtr("rotTransitiongradValueZ", rotTransitiongradValueZ);

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D *qTransitionUnsteadyn1 = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *qTransitionUnsteadyn2 = new RDouble4D(I, J, K, M, fortranArray);

        RDouble4D *residualTransitionUnsteadyn1 = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *residualTransitionUnsteadyn2 = new RDouble4D(I, J, K, M, fortranArray);

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble4D *residualTransitionUnsteadyTemporary = new RDouble4D(I, J, K, M, fortranArray);

        grid->UpdateDataPtr("q_transition_unsteady_n1", qTransitionUnsteadyn1);
        grid->UpdateDataPtr("q_transition_unsteady_n2", qTransitionUnsteadyn2);
        grid->UpdateDataPtr("res_transition_unsteady_n1", residualTransitionUnsteadyn1);
        grid->UpdateDataPtr("res_transition_unsteady_n2", residualTransitionUnsteadyn2);
        grid->UpdateDataPtr("res_transition_unsteady_tmp", residualTransitionUnsteadyTemporary);

        //! Statistical variables for DES simulation.
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        if (ifStaticsFlowField > 0)
        {
            RDouble4D *qAverageTransition = new RDouble4D(I, J, K, M, fortranArray);
            grid->UpdateDataPtr("qAverage_transition", qAverageTransition);
        }
    }

    int numberOfDimensions = GetDim();
    Range D(1, numberOfDimensions);

    RDouble5D *matrixTransitionLeft, *matrixTransitionRight;

    matrixTransitionLeft = new RDouble5D(I, J, K, M, D, fortranArray);
    matrixTransitionRight = new RDouble5D(I, J, K, M, D, fortranArray);

    grid->UpdateDataPtr("mat_transitionl", matrixTransitionLeft);
    grid->UpdateDataPtr("mat_transitionr", matrixTransitionRight);
}

void TransitionSolverStr::DeAllocateGlobalVar(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    if (nTransitionEquation == 0)
    {
        return;
    }

    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D *qTransition = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    RDouble4D *residualTransition = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));
    RDouble4D *spectrumTransition = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_transition"));

    delete qTransition;
    delete residualTransition;
    delete spectrumTransition;

    RDouble4D *gradTransitionX = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionX"));
    RDouble4D *gradTransitionY = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionY"));
    RDouble4D *gradTransitionZ = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionZ"));

    delete gradTransitionX;
    delete gradTransitionY;
    delete gradTransitionZ;

    RDouble4D *gradTransitionCellCenterX = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionCellCenterX"));
    RDouble4D *gradTransitionCellCenterY = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionCellCenterY"));
    RDouble4D *gradTransitionCellCenterZ = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionCellCenterZ"));

    delete gradTransitionCellCenterX;
    delete gradTransitionCellCenterY;
    delete gradTransitionCellCenterZ;

    RDouble4D *rotTransitiongradValueX = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTransitiongradValueX"));
    RDouble4D *rotTransitiongradValueY = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTransitiongradValueY"));
    RDouble4D *rotTransitiongradValueZ = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTransitiongradValueZ"));
    delete rotTransitiongradValueX;
    delete rotTransitiongradValueY;
    delete rotTransitiongradValueZ;

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D *qTransitionUnsteadyn1 = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n1"));
        RDouble4D *qTransitionUnsteadyn2 = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n2"));
        RDouble4D *residualTransitionUnsteadyn1 = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n1"));
        RDouble4D *residualTransitionUnsteadyn2 = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n2"));

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble4D *residualTransitionUnsteadyTemporary = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_tmp"));

        delete qTransitionUnsteadyn1;
        delete qTransitionUnsteadyn2;
        delete residualTransitionUnsteadyn1;
        delete residualTransitionUnsteadyn2;
        delete residualTransitionUnsteadyTemporary;

        //! Statistical variables for DES simulation.
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        if (ifStaticsFlowField > 0)
        {
            RDouble4D *qAverageTransition = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("qAverage_transition"));
            delete qAverageTransition;
        }
    }

    RDouble5D *matrixTransitionLeft, *matrixTransitionRight;
    matrixTransitionLeft = reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_transitionl"));
    matrixTransitionRight = reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_transitionr"));

    delete matrixTransitionLeft;
    delete matrixTransitionRight;
}

bool TransitionSolverStr::JudgeIfRestart()
{
    string transitionFile = ".\results\transition.dat";
    GlobalDataBase::GetData("transitionFile", &transitionFile, PHSTRING, 1);

    if (PHMPI::IsParallelRun())
    {
        transitionFile = PHSPACE::AddSymbolToFileName(transitionFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(transitionFile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

bool TransitionSolverStr::JudgeIfProtectedRestart()
{
    string protectionTransitionFile0 = "./results/transition0.dat";
    GlobalDataBase::GetData("protectionTransitionFile0", &protectionTransitionFile0, PHSTRING, 1);
    if (PHMPI::IsParallelRun())
    {
        protectionTransitionFile0 = PHSPACE::AddSymbolToFileName(protectionTransitionFile0, "_", 0);
    }

    string protectionTransitionFile1 = "./results/transition1.dat";
    GlobalDataBase::GetData("protectionTransitionFile1", &protectionTransitionFile1, PHSTRING, 1);
    if (PHMPI::IsParallelRun())
    {
        protectionTransitionFile1 = PHSPACE::AddSymbolToFileName(protectionTransitionFile1, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile_0(protectionTransitionFile0.c_str(), ios::in);
    if (infile_0)
    {
        restart_flag = true;
        infile_0.close();
        infile_0.clear();
    }

    ifstream infile_1(protectionTransitionFile1.c_str(), ios::in);
    if (infile_1)
    {
        restart_flag = true;
        infile_1.close();
        infile_1.clear();
    }

    return restart_flag;
}

bool TransitionSolverStr::JudgeIfReadAverage()
{
    string transitionFile = ".\results\transition.dat";
    GlobalDataBase::GetData("transitionFile", &transitionFile, PHSTRING, 1);
    string averageTransitionFile = AddSymbolToFileName(transitionFile, "_", "Average");

    if (PHMPI::IsParallelRun())
    {
        averageTransitionFile = PHSPACE::AddSymbolToFileName(averageTransitionFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(averageTransitionFile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

void TransitionSolverStr::InitFlowAsReadingProtectedRestart(ActionKey *actkey)
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
        RestartFile = 0;
    }

    ifstream infile_1(protectionTransitionFile1.c_str(), ios::in);
    if (infile_1)
    {
        RestartFile = 1;
    }

    if (infile_0 && infile_1)     //If both flow0 and flow1 exist,determine which file to use by the last modified time.
    {
        RDouble TimeFile0, TimeFile1;
        TIME_SPACE::GetFileModifiedTime(protectionTransitionFile0, TimeFile0);
        TIME_SPACE::GetFileModifiedTime(protectionTransitionFile1, TimeFile1);
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

void TransitionSolverStr::DumpRestartH5(ActionKey *actkey)
{
    using namespace PHMPI;
    int currentProcessorID = GetCurrentProcessorID();

    int version = 1;
    int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (currentProcessorID == GetServerProcessorID())
    {
        WriteData(actkey->filepos, &version, "Version");
        WriteData(actkey->filepos, &outIterStep, "outnstep");

        if (IsNeedStatistics())
        {
            int nStatisticalStep = 0;
            GlobalDataBase::GetData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            WriteData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
        }
    }

    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    StructGrid *strGrid = StructGridCast(GetGrid(actkey->level));
    int nTotalCell = strGrid->GetNTotalCell();

    int ist, ied, jst, jed, kst, ked;
    strGrid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    Range M(0, nTransitionEquation - 1);
    int mst = M.first();

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    WriteData(grploc, &nTotalCell, "nTotalCell");
    WriteData(grploc, &qTransition(ist, jst, kst, mst), "q_transition");

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady == 1)
    {
        RDouble4D &qTransitionUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n1"));
        WriteData(grploc, &qTransitionUnsteadyn1(ist, jst, kst, mst), "q_transition_unsteady_n1");

        RDouble4D &qTransitionUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n2"));
        WriteData(grploc, &qTransitionUnsteadyn2(ist, jst, kst, mst), "q_transition_unsteady_n2");

        RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));
        WriteData(grploc, &residualTransition(ist, jst, kst, mst), "res_transition");

        RDouble4D &residualTransitionUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n1"));
        WriteData(grploc, &residualTransitionUnsteadyn1(ist, jst, kst, mst), "res_transition_unsteady_n1");

        RDouble4D &residualTransitionUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n2"));
        WriteData(grploc, &residualTransitionUnsteadyn2(ist, jst, kst, mst), "res_transition_unsteady_n2");

        if (IsNeedStatistics())
        {
            RDouble4D &qAverageTransition = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage_transition"));
            WriteData(grploc, &qAverageTransition(ist, jst, kst, mst), "qAverageTransition");
        }
    }

    H5Gclose(grploc);
}

void TransitionSolverStr::ReadRestartH5(ActionKey *actkey)
{
    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    StructGrid *strGrid = StructGridCast(GetGrid(actkey->level));
    int nTotalCell = strGrid->GetNTotalCell();

    int ist, ied, jst, jed, kst, ked;
    strGrid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    Range M(0, nTransitionEquation - 1);
    int mst = M.first();

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));

    int outIterstepofNS = GlobalDataBase::GetIntParaFromDB("outnstep");
    int outIterStepofTransition = 0;
    ReadData(actkey->filepos, &outIterStepofTransition, "outnstep");

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    int nTotalCellRestart = 0;
    ReadData(grploc, &nTotalCellRestart, "nTotalCell");
    if (nTotalCellRestart != nTotalCell)
    {
        ostringstream erroeInfo;
        erroeInfo << " Error: the cell number in transition.dat is not equal to the cell number in grid file !" << endl;
        TK_Exit::ExceptionExit(erroeInfo.str());
    }

    ReadData(grploc, &qTransition(ist, jst, kst, mst), "q_transition");

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D &qTransitionUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n1"));
        RDouble4D &qTransitionUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n2"));
        RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));
        RDouble4D &residualTransitionUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n1"));
        RDouble4D &residualTransitionUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n2"));
        RDouble4D &qAverageTransition = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage_transition"));

        int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();

        if (ifStartFromSteadyResults)
        {
            if (grid->GetGridID()->GetIndex() == 0)
            {
                PrintToWindow("Restart from steady Transition flow field, reset outer step to be zero!\n");
            }

            outIterStepofTransition = 0;
            GlobalDataBase::UpdateData("outnstep", &outIterStepofTransition, PHINT, 1);

            qTransitionUnsteadyn1 = qTransition;
            qTransitionUnsteadyn2 = qTransition;

            residualTransition = 0.0;
            residualTransitionUnsteadyn1 = 0.0;
            residualTransitionUnsteadyn2 = 0.0;
        }
        else
        {
            ReadData(grploc, &qTransitionUnsteadyn1(ist, jst, kst, mst), "q_transition_unsteady_n1");
            ReadData(grploc, &qTransitionUnsteadyn2(ist, jst, kst, mst), "q_transition_unsteady_n2");
            ReadData(grploc, &residualTransition(ist, jst, kst, mst), "res_transition");
            ReadData(grploc, &residualTransitionUnsteadyn1(ist, jst, kst, mst), "res_transition_unsteady_n1");
            ReadData(grploc, &residualTransitionUnsteadyn2(ist, jst, kst, mst), "res_transition_unsteady_n2");
        }

        int nStatisticalStep = 0;
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
        bool isReadAverageFlow = ifStaticsFlowField && (outIterStepofTransition >= startStatisticStep);
        if (isReadAverageFlow)
        {
            if (ifStartFromSteadyResults)
            {
                qAverageTransition = qTransition;
                GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            }
            else
            {
                ReadData(grploc, &qAverageTransition(ist, jst, kst, mst), "qAverageTransition");

                ReadData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
                GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            }
        }
    }

    H5Gclose(grploc);

    CompareOutStepOfFlowfieldFile(outIterstepofNS, outIterStepofTransition);
}

void TransitionSolverStr::GetResidual(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid *grid = StructGridCast(GetGrid(level));

    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    PHSPACE::ComputeResidualonGrid(grid, actkey, residualTransition, nTransitionEquation);
}

void TransitionSolverStr::ZeroResiduals(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    if (nTransitionEquation == 0)
    {
        return;
    }

    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("res_transition"));
    residualTransition = 0;
}

void TransitionSolverStr::InitSpectrum(Grid *gridIn)
{
    RDouble4D &spectrumTransition = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("spec_transition"));
    spectrumTransition = 0;
}

void TransitionSolverStr::InitFlowAsRestart()
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    if (nTransitionEquation == 0)
    {
        return;
    }

    StructGrid *grid = StructGridCast(GetGrid(0));

    //int outIterStep = 0;
    //GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

    RDouble *freeStreamTransitionVar = parameters->GetFreeStreamTransitionVar();
    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Range M(0, nTransitionEquation - 1);
    int mst = M.first();
    int med = M.last();

    for (int m = mst; m <= med; ++m)
    {
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    qTransition(i, j, k, m) = freeStreamTransitionVar[m];
                }
            }
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D &qTransitionUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n1"));
        RDouble4D &qTransitionUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n2"));

        qTransitionUnsteadyn1 = qTransition;
        qTransitionUnsteadyn2 = qTransitionUnsteadyn1;

        RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));
        RDouble4D &residualTransitionUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n1"));
        RDouble4D &residualTransitionUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n2"));

        residualTransition = 0.0;
        residualTransitionUnsteadyn1 = 0.0;
        residualTransitionUnsteadyn2 = 0.0;

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble4D &residualTransitionUnsteadyTemporary = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_tmp"));

        residualTransitionUnsteadyTemporary = 0.0;
    }
}

void TransitionSolverStr::FillField(Grid *gridIn, FieldProxy *fieldProxy, RDouble value)
{
    RDouble4D &field = fieldProxy->GetField_STR();

    field = value;
}

void TransitionSolverStr::FillField(Grid *gridIn, FieldProxy *fieldProxyTarget, FieldProxy *fieldProxySource)
{
    RDouble4D &fieldTarget = fieldProxyTarget->GetField_STR();
    RDouble4D &fieldSource = fieldProxySource->GetField_STR();

    fieldTarget = fieldSource;
}

FieldProxy *TransitionSolverStr::CreateFieldProxy(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);
    Range M(0, nTransitionEquation - 1);

    RDouble4D *field = new RDouble4D(I, J, K, M, fortranArray);

    FieldProxy *fieldProxy = new FieldProxy();
    fieldProxy->SetField_STR(field, true);
    return fieldProxy;
}

FieldProxy *TransitionSolverStr::CreateFlowFieldProxy(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int numberOfLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    int numberOfChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int numberOfEquation = numberOfLaminar + numberOfChemical;

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);
    Range M(0, numberOfEquation - 1);

    RDouble4D *field = new RDouble4D(I, J, K, M, fortranArray);

    FieldProxy *fieldProxy = new FieldProxy();
    fieldProxy->SetField_STR(field, true);
    return fieldProxy;
}

FieldProxy *TransitionSolverStr::GetFieldProxy(Grid *gridIn, const string &fieldName)
{
    Grid *grid = gridIn;

    RDouble4D *field = reinterpret_cast<RDouble4D *> (grid->GetDataPtr(fieldName));

    FieldProxy *fieldProxy = new FieldProxy();

    fieldProxy->SetField_STR(field);

    return fieldProxy;
}

FieldProxy *TransitionSolverStr::GetResidualProxy(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));

    FieldProxy *residualProxy = new FieldProxy();

    residualProxy->SetField_STR(&residualTransition);

    return residualProxy;
}

void TransitionSolverStr::StoreRhsByResidual(Grid *gridIn, FieldProxy *rightHandSideProxy)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));

    RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();

    rightHandSide = residualTransition;
}

void TransitionSolverStr::InitResidual(Grid *gridIn, FieldProxy *rightHandSideProxy)
{
    StructGrid *grid = StructGridCast(gridIn);

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Range M(0, nTransitionEquation - 1);
    int mst = M.first();
    int med = M.last();

    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));
    RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();

    for (int m = mst; m <= med; ++m)
    {
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    RDouble &resTransition = residualTransition(i, j, k, m);
                    resTransition = -rightHandSide(i, j, k, m);
                }
            }
        }
    }
}

void TransitionSolverStr::RecoverResidual(Grid *gridIn, FieldProxy *rightHandSideProxy)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));
    RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();

    if (grid->GetLevel() != 0)
    {
        residualTransition = rightHandSide;
    }
}

LIB_EXPORT void TransitionSolverStr::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_TransitionSolver();
    controlParameters->Init();
}

LIB_EXPORT Param_TransitionSolver *TransitionSolverStr::GetControlParameters()
{
    return static_cast <Param_TransitionSolver *> (controlParameters);
}

void TransitionSolverStr::LoadQ(Grid *gridIn, FieldProxy *qTrubProxy)
{
    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("q_transition"));
    RDouble4D &qTransitionInProxy = qTrubProxy->GetField_STR();

    qTransitionInProxy = qTransition;
}

void TransitionSolverStr::RungeKuttaResidual(Grid *gridIn, FieldProxy *dqTransitionProxy, RDouble RKCoeff)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &dqTransition = dqTransitionProxy->GetField_STR();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble3D &dt = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("dt"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    Range M(0, nTransitionEquation - 1);
    int mst = M.first();
    int med = M.last();

    for (int m = mst; m <= med; ++m)
    {
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    RDouble &dq = dqTransition(i, j, k, m);
                    dq *= dt(i, j, k) * RKCoeff;
                }
            }
        }
    }
}

//! Load rhs to residual res stored in grid.
void TransitionSolverStr::LoadResiduals(Grid *gridIn, FieldProxy *rightHandSideProxy)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();

    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Range M(0, nTransitionEquation - 1);
    int mst = M.first();
    int med = M.last();

    for (int m = mst; m <= med; ++m)
    {
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    RDouble &resTransition = residualTransition(i, j, k, m);
                    resTransition = -rightHandSide(i, j, k, m);
                }
            }
        }
    }
}

RDouble TransitionSolverStr::UnsteadyConvergence(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        return zero;
    }

    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble4D &qLaminar = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    RDouble4D &qTransitionUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n1"));
    RDouble4D &qTransitionUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n2"));
    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));

    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble *qTransition0 = new RDouble[nTransitionEquation];
    RDouble *qTransition1 = new RDouble[nTransitionEquation];
    RDouble *qTransition2 = new RDouble[nTransitionEquation];

    RDouble *qTransitionConserve0 = new RDouble[nTransitionEquation];
    RDouble *qTransitionConserve1 = new RDouble[nTransitionEquation];
    RDouble *qTransitionConserve2 = new RDouble[nTransitionEquation];

    using namespace IDX;

    RDouble summedPrimitive = zero;
    RDouble summedConserve = zero;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                for (int m = 0; m < nTransitionEquation; ++m)
                {
                    qTransition0[m] = qTransition(i, j, k, m);
                    qTransition1[m] = qTransitionUnsteadyn1(i, j, k, m);
                    qTransition2[m] = qTransitionUnsteadyn2(i, j, k, m);
                }

                RDouble rho = 1.0;

                if (nTransitionEquation >= 2)
                {
                    rho = qLaminar(i, j, k, IR);
                }

                for (int m = 0; m < nTransitionEquation; ++m)
                {
                    qTransitionConserve0[m] = rho * qTransition0[m];
                    qTransitionConserve1[m] = rho * qTransition1[m];
                    qTransitionConserve2[m] = rho * qTransition2[m];
                }

                for (int m = 0; m < nTransitionEquation; ++m)
                {
                    RDouble dqTransitionPseudo = residualTransition(i, j, k, m);           // qn+1, p+1 - qn+1, p
                    RDouble dqTransitionPhysical = qTransitionConserve0[m] - qTransitionConserve1[m];  // qn+1, p+1 - qn
                    summedPrimitive += dqTransitionPseudo * dqTransitionPseudo;
                    summedConserve += dqTransitionPhysical * dqTransitionPhysical;
                }
            }
        }
    }

    delete [] qTransition0;    qTransition0 = NULL;
    delete [] qTransition1;    qTransition1 = NULL;
    delete [] qTransition2;    qTransition2 = NULL;
    delete [] qTransitionConserve0;    qTransitionConserve0 = NULL;
    delete [] qTransitionConserve1;    qTransitionConserve1 = NULL;
    delete [] qTransitionConserve2;    qTransitionConserve2 = NULL;

    RDouble unsteadyConvergence = sqrt(ABS(summedPrimitive / (summedConserve + SMALL)));
    return unsteadyConvergence;
}

void TransitionSolverStr::UpdateUnsteadyFlow(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        return;
    }

    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    //! The volume at current timestep. If it is dualtime step method, it represents the volume at timestep of n+1.
    RDouble3D &volume = grid->GetVolume(0);
    //! The volume at timestep of n  .
    RDouble3D &volumeUnsteadyn1 = grid->GetVolume(1);
    //! The volume at timestep of n-1.
    RDouble3D &volumeUnsteadyn2 = grid->GetVolume(2);

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    RDouble4D &qTransitionUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n1"));
    RDouble4D &qTransitionUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n2"));

    RDouble4D &residualTransitionUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n1"));
    RDouble4D &residualTransitionUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n2"));
    RDouble4D &residualTransitionUnsteadyTemporary = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_tmp"));

    int nTransitionEquation = parameters->GetNTransitionEquation();

    qTransitionUnsteadyn2 = qTransitionUnsteadyn1;
    qTransitionUnsteadyn1 = qTransition;

    volumeUnsteadyn2 = volumeUnsteadyn1;
    volumeUnsteadyn1 = volume;

    //! Here the current outnstep is over, the value of the stored resTmp should be assigned to resn1 for the next outnstep.
    //! It should be noticed that residualUnsteadyTemporary only contain the inviscid and viscous flux.
    residualTransitionUnsteadyn2 = residualTransitionUnsteadyn1;
    residualTransitionUnsteadyn1 = residualTransitionUnsteadyTemporary;

    //! Statistical variables for DES simulation.
    int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
    if (ifStaticsFlowField > 0)
    {
        int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");

        if (outIterStep >= startStatisticStep)
        {
            int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
            RDouble statisticalTimePeriod = GlobalDataBase::GetDoubleParaFromDB("statisticalTimePeriod");
            RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

            RDouble coeff1 = 1.0 / (nStatisticalStep * 1.0);
            if (statisticalTimePeriod > 0.0)
            {
                coeff1 = MAX(coeff1, physicalTimeStep / statisticalTimePeriod);
            }
            RDouble coeff2 = 1.0 - coeff1;

            RDouble4D &qAverageTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("qAverage_transition"));

            for (int k = kst; k <= ked; ++k)
            {
                for (int j = jst; j <= jed; ++j)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        for (int m = 0; m < nTransitionEquation; ++m)
                        {
                            qAverageTransition(i, j, k, m) = coeff2 * qAverageTransition(i, j, k, m) + coeff1 * qTransition(i, j, k, m);
                        }
                    }
                }
            }
        }
    }
}

void TransitionSolverStr::DualTimeSource(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        return;
    }

    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble4D &qLaminar = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    RDouble4D &qTransitionUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n1"));
    RDouble4D &qTransitionUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition_unsteady_n2"));

    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));
    RDouble4D &residualTransitionUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n1"));
    RDouble4D &residualTransitionUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_n2"));
    RDouble4D &residualTransitionUnsteadyTemporary = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition_unsteady_tmp"));

    //! The volume at current timestep. If it is dualtime step method, it represents the volume at timestep of n+1.
    RDouble3D &volume = grid->GetVolume(0);
    //! The volume at timestep of n  .
    RDouble3D &volumeUnsteadyn1 = grid->GetVolume(1);
    //! the volume at timestep of n-1.
    RDouble3D &volumeUnsteadyn2 = grid->GetVolume(2);

    int nTransitionEquation = parameters->GetNTransitionEquation();

    //! The primitive variables.
    RDouble *qTransitionPrimtive0 = new RDouble[nTransitionEquation];
    RDouble *qTransitionPrimtive1 = new RDouble[nTransitionEquation];
    RDouble *qTransitionPrimtive2 = new RDouble[nTransitionEquation];

    //! The conservative variables.
    RDouble *qTransitionConserve0 = new RDouble[nTransitionEquation];
    RDouble *qTransitionConserve1 = new RDouble[nTransitionEquation];
    RDouble *qTransitionConserve2 = new RDouble[nTransitionEquation];

    //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux,
    //! used in UpdateUnsteadyFlow for the Rn in the dualtime method.
    residualTransitionUnsteadyTemporary = residualTransition;

    //! Computation of dualtime coefficients, including three coefficients for Residual of R(p), R(n) and R(n-1);
    //! and three coefficients for conservative variables of qcsv(p), qcsv(n), and qcsv(n-1).
    RDouble dualTimeCoefficient[7];
    const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
    ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);

    RDouble dualTimeResC1 = dualTimeCoefficient[0];
    RDouble dualTimeResC2 = dualTimeCoefficient[1];
    RDouble dualTimeResC3 = dualTimeCoefficient[2];
    RDouble dualTimeQC1 = dualTimeCoefficient[3];
    RDouble dualTimeQC2 = dualTimeCoefficient[4];
    RDouble dualTimeQC3 = dualTimeCoefficient[5];

    using namespace IDX;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                for (int m = 0; m < nTransitionEquation; ++m)
                {
                    qTransitionPrimtive0[m] = qTransition(i, j, k, m);
                    qTransitionPrimtive1[m] = qTransitionUnsteadyn1(i, j, k, m);
                    qTransitionPrimtive2[m] = qTransitionUnsteadyn2(i, j, k, m);
                }

                RDouble roe = 1.0;

                if (nTransitionEquation >= 2)
                {
                    roe = qLaminar(i, j, k, IR);
                }

                for (int m = 0; m < nTransitionEquation; ++m)
                {
                    qTransitionConserve0[m] = roe * qTransitionPrimtive0[m];
                    qTransitionConserve1[m] = roe * qTransitionPrimtive1[m];
                    qTransitionConserve2[m] = roe * qTransitionPrimtive2[m];
                }

                for (int m = 0; m < nTransitionEquation; ++m)
                {
                    RDouble dualSrcRes = dualTimeResC1 * residualTransition(i, j, k, m) +
                        dualTimeResC2 * residualTransitionUnsteadyn1(i, j, k, m) +
                        dualTimeResC3 * residualTransitionUnsteadyn2(i, j, k, m);

                    RDouble dualSrcQ = dualTimeQC1 * qTransitionConserve0[m] * volume(i, j, k) +
                        dualTimeQC2 * qTransitionConserve1[m] * volumeUnsteadyn1(i, j, k) +
                        dualTimeQC3 * qTransitionConserve2[m] * volumeUnsteadyn2(i, j, k);

                    residualTransition(i, j, k, m) = dualSrcRes + dualSrcQ;
                }
            }
        }
    }

    delete [] qTransitionPrimtive0;
    delete [] qTransitionPrimtive1;
    delete [] qTransitionPrimtive2;
    delete [] qTransitionConserve0;
    delete [] qTransitionConserve1;
    delete [] qTransitionConserve2;
}


void TransitionSolverStr::SourceFluxTwoEquation(Grid *gridIn)
{
    //! This routine calculates the source terms relative to k and omega equations for k-omega
    //! two-equation models.
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &qPrimitive = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &viscousLaminar = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble4D &qTurbulence = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    RDouble4D &spectrumTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_transition"));
    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));

    RDouble3D &walldist = *grid->GetWallDist();

    RDouble3D &vol = *(grid->GetCellVolume());

    Param_TransitionSolver *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble refMaNumber = parameters->GetRefMachNumber();
    string viscousName = parameters->GetViscousName();

    //! To determine if compressible correction used or not.
    int isCC = 0;
    if (GlobalDataBase::IsExist("compressibleCorrection", PHINT, 1))
    {
        isCC = GlobalDataBase::GetIntParaFromDB("compressibleCorrection");
    }

    int freeturbIntensitySRModify = 0;
    if (GlobalDataBase::IsExist("freeturbIntensitySRModify", PHINT, 1))
    {
        freeturbIntensitySRModify = GlobalDataBase::GetIntParaFromDB("freeturbIntensitySRModify");
    }
    double freeDecayXLocation = 0.0;
    if (GlobalDataBase::IsExist("freeDecayXLocation", PHDOUBLE, 1))
    {
        freeDecayXLocation = GlobalDataBase::GetDoubleParaFromDB("freeDecayXLocation");
    }

    RDouble compressibleFix = 1.0 + 0.20 * pow(refMaNumber, 2.6);
    int transitionMaFix = 1;
    if (GlobalDataBase::IsExist("transitionMaFix", PHINT, 1))
    {
        transitionMaFix = GlobalDataBase::GetIntParaFromDB("transitionMaFix");
    }
    if (transitionMaFix == 0)
    {
        compressibleFix = 1.0;
    }

    RDouble pressure = 0.0, density = 0.0, ke = 0.0, kw = 0.0;
    RDouble sij2 = 0.0, divv = 0.0, mut = 0.0;
    RDouble volume = 0.0;
    RDouble dudx0 = 0.0, dudy0 = 0.0, dudz0 = 0.0, dvdx0 = 0.0, dvdy0 = 0.0, dvdz0 = 0.0, dwdx0 = 0.0, dwdy0 = 0.0, dwdz0 = 0.0;
    RDouble s11 = 0.0, s22 = 0.0, s33 = 0.0, s12 = 0.0, s13 = 0.0, s23 = 0.0, w12 = 0.0, w13 = 0.0, w23 = 0.0;

    Range II, JJ, KK;
    GetRange(ni, nj, nk, -2, 1, II, JJ, KK);

    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    using namespace IDX;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble ca1 = parameters->Getca1();
    RDouble ca2 = parameters->Getca2();
    RDouble ce1 = parameters->Getce1();
    RDouble ce2 = parameters->Getce2();
    RDouble cct = parameters->Getcct();
    RDouble  s1 = parameters->Gets1();
    RDouble ccf = parameters->Getccf();

    RDouble3D &gamaeff = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("gamaeff"));
    RDouble3D &rescf = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("rescf"));

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                density = qPrimitive(i, j, k, IR);
                pressure = qPrimitive(i, j, k, IP);
                ke = qTurbulence(i, j, k, IKE);
                kw = qTurbulence(i, j, k, IKW);
                mut = viscousTurbulence(i, j, k);
                volume = vol(i, j, k);

                dudx0 = gradUVWTCellCenterX(i, j, k, 0);
                dudy0 = gradUVWTCellCenterY(i, j, k, 0);
                dudz0 = gradUVWTCellCenterZ(i, j, k, 0);
                dvdx0 = gradUVWTCellCenterX(i, j, k, 1);
                dvdy0 = gradUVWTCellCenterY(i, j, k, 1);
                dvdz0 = gradUVWTCellCenterZ(i, j, k, 1);
                dwdx0 = gradUVWTCellCenterX(i, j, k, 2);
                dwdy0 = gradUVWTCellCenterY(i, j, k, 2);
                dwdz0 = gradUVWTCellCenterZ(i, j, k, 2);

                s11 = dudx0;
                s22 = dvdy0;
                s33 = dwdz0;
                s12 = half * (dudy0 + dvdx0);
                s13 = half * (dudz0 + dwdx0);
                s23 = half * (dvdz0 + dwdy0);
                w12 = half * (dudy0 - dvdx0);
                w13 = half * (dudz0 - dwdx0);
                w23 = half * (dvdz0 - dwdy0);

                sij2 = two * (s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23));    //! Modulus of S.
                divv = s11 + s22 + s33;    //! Divergence of velocity.

                RDouble mul = viscousLaminar(i, j, k);
                RDouble um = qPrimitive(i, j, k, IU);
                RDouble vm = qPrimitive(i, j, k, IV);
                RDouble wm = qPrimitive(i, j, k, IW);
                RDouble intermittency = qTransition(i, j, k, IGAMA);
                RDouble Rectabar = qTransition(i, j, k, IRECT);

                RDouble wallDistance = walldist(i, j, k);

                RDouble vorx = dvdz0 - dwdy0;
                RDouble vory = dudz0 - dwdx0;
                RDouble vorz = dudy0 - dvdx0;

                RDouble vorticity = DISTANCE(vorx, vory, vorz);
                RDouble strainRate = sqrt(sij2);

                RDouble absU = MAX(DISTANCE(um, vm, wm), SMALL);

                RDouble RT = ViscosityRatio(density, mul, ke, kw, refReNumber);
                RDouble Rev = ReynoldsNumberBasedOnStrainRate(density, wallDistance, mul, strainRate, refReNumber);
                RDouble rectac = TransitionOnsetMomentumThicknessReynoldsCaliBrated(Rectabar);

                RDouble Rw = ReynoldsNumberBasedOnDissipation(density, wallDistance, mul, kw, refReNumber);
                RDouble flength = HighReynoldsCorrectionOfFlength(Rw, FlengthCaliBrated(Rectabar));
                RDouble Fonset = TransitionOnsetFunction(Rev, rectac, RT);
                RDouble Fturb = ControlFunctionFturb(RT);

                RDouble production1OfGama = ca1 * flength * density * strainRate * sqrt(intermittency * Fonset);
                RDouble production2OfGama = production1OfGama * (-ce1 * intermittency);

                RDouble productionOfGama = production1OfGama + production2OfGama;

                RDouble specGamap = MIN(0.0, 1.5 * production2OfGama + 0.5 * production1OfGama) / (density * intermittency);

                RDouble destruction2OfGama = ca2 * Fturb * density * vorticity * intermittency;
                RDouble destruction1OfGama = ce2 * intermittency * destruction2OfGama;

                RDouble destructionOfGama = destruction1OfGama + destruction2OfGama;
                RDouble specGamad = MAX(0.0, 2.0 * destruction1OfGama + destruction2OfGama) / (density * intermittency);

                residualTransition(i, j, k, IGAMA) += (productionOfGama - destructionOfGama) * volume;
                spectrumTransition(i, j, k, IGAMA) += (specGamad - specGamap) * volume;

                RDouble Fctat = 1.0, Fctatcf = 1.0;
                BlendingFunctionOfFctat(intermittency, Rw, Rectabar, vorticity, mul, density, absU, wallDistance, refReNumber, ce2, Fctat, Fctatcf);
                RDouble tscl = TimeScaleInSourceTerm(density, absU, mul, mut, refReNumber);
                RDouble TU = TurbulenceIntensity(absU, ke);
                gamaeff(i, j, k) = SeparationCorrectionOfIntermittency(intermittency, Rev, rectac, RT, Fctat, s1, refMaNumber);

                RDouble dUds = AccelerationAlongStreamline(um, vm, wm, dudx0, dudy0, dudz0, dvdx0, dvdy0, dvdz0, dwdx0, dwdy0, dwdz0);

                RDouble momentumThickness = 0.0;
                RDouble Rectat;
                const int subIter = 10;
                RDouble lamdacta0, lamdacta = 1.0;
                for (int iter = 0; iter < subIter; ++iter)
                {
                    lamdacta0 = lamdacta;
                    lamdacta = PressureGradientFunction(density, momentumThickness, mul, dUds, refReNumber);
                    lamdacta = lamdacta * compressibleFix;
                    lamdacta = MAX(-0.1, MIN(0.1, lamdacta));
                    RDouble Flamdacta = EmpiricalCorrelationOfFlamdacta(TU, lamdacta);
                    Rectat = EmpiricalCorrelationOfRectat(TU, Flamdacta);
                    momentumThickness = MomentumThickness(density, absU, mul, Rectat, refReNumber);
                    if (abs(lamdacta0 - lamdacta) <= SMALL)
                    {
                        break;
                    }
                }

                RDouble omegaStreamWise = abs((dwdy0 - dvdz0) * um + (dudz0 - dwdx0) * vm + (dvdx0 - dudy0) * wm) / absU;
                RDouble heightCrossFlow = wallDistance * omegaStreamWise / absU;
                rescf(i, j, k) = ReynoldsNumberBasedOnScf(heightCrossFlow, absU, density, mul, mut, refReNumber);

                RDouble coecommon = cct * (1.0 - Fctat) / MAX(tscl, SMALL);

                RDouble prodRe = coecommon * density * Rectat;
                RDouble destRe = coecommon * density * Rectabar;

                RDouble prodRecf = cct * density * ccf * MIN(rescf(i, j, k) - Rectabar, 0.0) * Fctatcf / tscl;
                RDouble specRe = destRe / (density * Rectabar);

                residualTransition(i, j, k, IRECT) += (prodRe + prodRecf - destRe) * volume;
                spectrumTransition(i, j, k, IRECT) += specRe * volume;
            }
        }
    }
}

void TransitionSolverStr::Diagonal(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble4D &faceVectorX = *(grid->GetFaceVectorX());
    RDouble4D &faceVectorY = *(grid->GetFaceVectorY());
    RDouble4D &faceVectorZ = *(grid->GetFaceVectorZ());
    RDouble4D &area = *(grid->GetFaceArea());
    RDouble3D &volume = *(grid->GetCellVolume());
    RDouble4D &faceNormalVelocity = *(grid->GetFaceNormalVelocity());

    RDouble4D &spectrumTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_transition"));
    RDouble4D &qLaminar = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &viscousLaminar = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTransition = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &dt = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("dt"));

    RDouble5D &matrixTransitionLeft = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_transitionl"));
    RDouble5D &matrixTransitionRight = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_transitionr"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble transitionCFLScale = parameters->GetTransitionCFLScale();
    string viscousName = parameters->GetViscousName();

    Range M(0, nTransitionEquation - 1);
    int mst = M.first();
    int med = M.last();

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble csrv = GlobalDataBase::GetDoubleParaFromDB("csrv");
    int ifLowSpeedPrecon = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");

    RDouble3D *preconCoefficient = NULL;

    if (ifLowSpeedPrecon != 0)
    {
        preconCoefficient = reinterpret_cast< RDouble3D *> (grid->GetDataPtr("preconCoefficient"));
    }

    matrixTransitionLeft = 0.0;
    matrixTransitionRight = 0.0;

    using namespace IDX;

    RDouble dualTimeSpectrumC1 = zero;
    RDouble dualTimeSpectrumC2 = one;

    //! If flow is unsteady, it needs to count the contribution of unsteady.
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble dualTimeCoefficient[7];
        const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
        dualTimeSpectrumC1 = -dualTimeCoefficient[3];
        dualTimeSpectrumC2 = dualTimeCoefficient[6];
    }

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                RDouble volumeCell = volume(i, j, k);
                RDouble dtCell = dt(i, j, k);
                RDouble specTransitionDelta = dualTimeSpectrumC1 * volumeCell + dualTimeSpectrumC2 / (transitionCFLScale * dtCell + SMALL);
                for (int m = mst; m <= med; ++m)
                {
                    RDouble &specTrans = spectrumTransition(i, j, k, m);
                    specTrans += specTransitionDelta;
                }
            }
        }
    }

    int nDim = GetDim();

    vector<RDouble> work(nTransitionEquation);

    //int il, jl, kl;
    for (int iSurface = 1; iSurface <= nDim; ++iSurface)
    {
        int il1, jl1, kl1;
        grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    int il, jl, kl;
                    grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl);

                    RDouble nx0 = faceVectorX(i, j, k, iSurface);
                    RDouble ny0 = faceVectorY(i, j, k, iSurface);
                    RDouble nz0 = faceVectorZ(i, j, k, iSurface);
                    RDouble vgn0 = faceNormalVelocity(i, j, k, iSurface) * area(i, j, k, iSurface);

                    RDouble nx1 = faceVectorX(il, jl, kl, iSurface);
                    RDouble ny1 = faceVectorY(il, jl, kl, iSurface);
                    RDouble nz1 = faceVectorZ(il, jl, kl, iSurface);
                    RDouble vgn1 = faceNormalVelocity(il, jl, kl, iSurface) * area(il, jl, kl, iSurface);

                    RDouble nx = half * (nx0 + nx1);
                    RDouble ny = half * (ny0 + ny1);
                    RDouble nz = half * (nz0 + nz1);
                    RDouble vgn = half * (vgn0 + vgn1);

                    RDouble ul = qLaminar(i, j, k, IU);
                    RDouble vl = qLaminar(i, j, k, IV);
                    RDouble wl = qLaminar(i, j, k, IW);
                    RDouble Vn = ul * nx + vl * ny + wl * nz - vgn;
                    if (ifLowSpeedPrecon != 0)
                    {
                        RDouble preconCoeff = (*preconCoefficient)(i, j, k);
                        Vn = half * Vn * (1.0 + preconCoeff);
                    }
                    RDouble absVn = ABS(Vn);

                    for (int m = M.first(); m <= M.last(); ++m)
                    {
                        RDouble &specTransition = spectrumTransition(i, j, k, m);
                        specTransition += absVn;

                        RDouble &matTransitionLeft = matrixTransitionLeft(il, jl, kl, m, iSurface);
                        RDouble &matTransitionRight = matrixTransitionRight(i, j, k, m, iSurface);
                        matTransitionLeft += half * (-Vn - absVn);
                        matTransitionRight += half * (Vn - absVn);
                    }
                }
            }
        }
    }

    RDouble dct = parameters->Getdct();
    RDouble df = parameters->Getdf();

    //int il, jl, kl;
    for (int iSurface = 1; iSurface <= nDim; ++iSurface)
    {
        int il1, jl1, kl1;
        grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

        Range M1(0, nTransitionEquation - 1);
        int mst1 = M1.first();
        int med1 = M1.last();

        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    int il, jl, kl;
                    grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl);

                    RDouble rho = qLaminar(i, j, k, IR);
                    RDouble muLaminar = viscousLaminar(i, j, k);
                    RDouble muTransition = viscousTransition(i, j, k);

                    RDouble cellVolume = volume(i, j, k) + SMALL;
                    RDouble oVolume = one / cellVolume;

                    RDouble nx0 = faceVectorX(i, j, k, iSurface);
                    RDouble ny0 = faceVectorY(i, j, k, iSurface);
                    RDouble nz0 = faceVectorZ(i, j, k, iSurface);

                    RDouble nx1 = faceVectorX(il, jl, kl, iSurface);
                    RDouble ny1 = faceVectorY(il, jl, kl, iSurface);
                    RDouble nz1 = faceVectorZ(il, jl, kl, iSurface);

                    RDouble nx = half * (nx0 + nx1);
                    RDouble ny = half * (ny0 + ny1);
                    RDouble nz = half * (nz0 + nz1);

                    RDouble ns = nx * nx + ny * ny + nz * nz;
                    RDouble ns2 = 2.0 * csrv * ns;

                    work[IGAMA] = (muLaminar + muTransition / df) / rho * ns2 * oVolume * oRefReNumber;
                    work[IRECT] = (muLaminar + muTransition) * dct / rho * ns2 * oVolume * oRefReNumber;

                    for (int m = mst1; m <= med1; ++m)
                    {
                        RDouble &specTransition = spectrumTransition(i, j, k, m);
                        specTransition += work[m];

                        RDouble &matTransitionLeft = matrixTransitionLeft(il, jl, kl, m, iSurface);
                        RDouble &matTransitionRight = matrixTransitionRight(i, j, k, m, iSurface);
                        matTransitionLeft += - half * work[m];
                        matTransitionRight += - half * work[m];
                    }
                }
            }
        }
    }
}

void TransitionSolverStr::ChangeTransitionQ(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    RDouble4D &qLaminar = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    Range M(0, nTransitionEquation - 1);
    int mst = M.first();
    int med = M.last();

    using namespace IDX;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                RDouble density = qLaminar(i, j, k, IR);
                for (int m = mst; m <= med; ++m)
                {
                    qTransition(i, j, k, m) *= density;
                }
            }
        }
    }
}

void TransitionSolverStr::RecoverTransitionQ(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    RDouble4D &qLaminar = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Range M(0, nTransitionEquation - 1);

    int mst = M.first();
    int med = M.last();

    using namespace IDX;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                RDouble rho = qLaminar(i, j, k, IR) + SMALL;
                for (int m = mst; m <= med; ++m)
                {
                    qTransition(i, j, k, m) /= rho;
                }
            }
        }
    }
}

void TransitionSolverStr::InviscidFlux(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    FieldProxy *qTransitionLeftPproxy = CreateFieldProxy(grid);
    FieldProxy *qTransitionRightPproxy = CreateFieldProxy(grid);
    FieldProxy *fluxProxy = CreateFieldProxy(grid);

    FieldProxy *qLaminarLeftProxy;
    FieldProxy *qLaminarRightProxy;

    int turbOrderStruct = 2;
    if (GlobalDataBase::IsExist("turbOrderStruct", PHINT, 1))
    {
        turbOrderStruct = GlobalDataBase::GetIntParaFromDB("turbOrderStruct");
    }
    if (turbOrderStruct == 2)
    {
        qLaminarLeftProxy = new FieldProxy();
        qLaminarRightProxy = new FieldProxy();
    }
    else
    {
        qLaminarLeftProxy = CreateFlowFieldProxy(grid);
        qLaminarRightProxy = CreateFlowFieldProxy(grid);
    }

    ChangeTransitionQ(grid);

    int nDim = GetDim();
    for (int iSurface = 1; iSurface <= nDim; ++iSurface)
    {
        if (turbOrderStruct == 2)
        {
            if (iSurface == 1)
            {
                qLaminarLeftProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1")));
                qLaminarRightProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1")));
            }
            else if (iSurface == 2)
            {
                qLaminarLeftProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2")));
                qLaminarRightProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2")));
            }
            else
            {
                qLaminarLeftProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3")));
                qLaminarRightProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3")));
            }
        }
        else if (turbOrderStruct == 1)
        {
            //! Interpolate u ,v, w, for transition equation.
            GetInvFaceQValueforTransition(grid, qLaminarLeftProxy, qLaminarRightProxy, iSurface);
        }

        //! Interpolate den*k den*w.
        GetInviscidFaceValue(grid, qTransitionLeftPproxy, qTransitionRightPproxy, iSurface);

        //! Correct the variables at the face where the BC condition is set to be wall.
        CorrectFaceVar(grid, qTransitionLeftPproxy, qTransitionRightPproxy, qLaminarLeftProxy, qLaminarRightProxy, iSurface);

        ComputeInviscidFlux(grid, qTransitionLeftPproxy, qTransitionRightPproxy, qLaminarLeftProxy, qLaminarRightProxy, fluxProxy, iSurface);

        LoadFlux(grid, fluxProxy, iSurface);
    }

    RecoverTransitionQ(grid);

    delete qTransitionLeftPproxy;
    delete qTransitionRightPproxy;

    delete qLaminarLeftProxy;
    delete qLaminarRightProxy;

    delete fluxProxy;
}

void TransitionSolverStr::GetInviscidFaceValue(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    Range M(0, nTransitionEquation - 1);

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));

    RDouble4D &qLeft = qLeftProxy->GetField_STR();
    RDouble4D &qRight = qRightProxy->GetField_STR();

    RDouble MUSCLCoefXk = GlobalDataBase::GetDoubleParaFromDB("MUSCLCoefXk");

    int il1, jl1, kl1;
    grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

    I = Range(1 - il1, ni - 1 + il1);
    J = Range(1 - jl1, nj - 1 + jl1);
    K = Range(1 - kl1, nk - 1 + kl1);
    if (nk == 1)
    {
        K.setRange(1, 1);
    }
    int turbOrderStruct = 2;
    if (GlobalDataBase::IsExist("turbOrderStruct", PHINT, 1))
    {
        turbOrderStruct = GlobalDataBase::GetIntParaFromDB("turbOrderStruct");
    }
    if (turbOrderStruct == 2)
    {
        RDouble4D *dql = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *dqr = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *dqlr = new RDouble4D(I, J, K, M, fortranArray);

        (*dql)(I, J, K, M) = qTransition(I, J, K, M) - qTransition(I - il1, J - jl1, K - kl1, M);
        (*dqr)(I, J, K, M) = qTransition(I + il1, J + jl1, K + kl1, M) - qTransition(I, J, K, M);
        (*dqlr)(I, J, K, M) = Smoothturb((*dql)(I, J, K, M), (*dqr)(I, J, K, M));

        qLeft(I, J, K, M) = qTransition(I, J, K, M) + 0.25 * (*dqlr)(I, J, K, M) * ((1.0 - MUSCLCoefXk * (*dqlr)(I, J, K, M)) * (*dql)(I, J, K, M) + (1.0 + MUSCLCoefXk * (*dqlr)(I, J, K, M)) * (*dqr)(I, J, K, M));
        qRight(I, J, K, M) = qTransition(I, J, K, M) - 0.25 * (*dqlr)(I, J, K, M) * ((1.0 - MUSCLCoefXk * (*dqlr)(I, J, K, M)) * (*dqr)(I, J, K, M) + (1.0 + MUSCLCoefXk * (*dqlr)(I, J, K, M)) * (*dql)(I, J, K, M));

        delete dql; dql = NULL;
        delete dqr; dqr = NULL;
        delete dqlr; dqlr = NULL;

    }
    else if (turbOrderStruct == 1)
    {
        qLeft(I, J, K, M) = qTransition(I, J, K, M);
        qRight(I, J, K, M) = qTransition(I, J, K, M);
    }
}

void TransitionSolverStr::GetInvFaceQValueforTransition(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range M(1, 3);

    FieldProxy *qProxy = GetFieldProxy(grid, "q");

    RDouble4D &q = qProxy->GetField_STR();
    RDouble4D &qLeft = qLeftProxy->GetField_STR();
    RDouble4D &qRight = qRightProxy->GetField_STR();

    int il1, jl1, kl1;
    grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

    I = Range(1 - il1, ni - 1 + il1);
    J = Range(1 - jl1, nj - 1 + jl1);
    K = Range(1 - kl1, nk - 1 + kl1);
    if (nk == 1)
    {
        K.setRange(1, 1);
    }

    qLeft(I, J, K, M) = q(I, J, K, M);
    qRight(I, J, K, M) = q(I, J, K, M);

    delete qProxy;
}

void TransitionSolverStr::CorrectFaceVar(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, FieldProxy *qLaminarLeftProxy, FieldProxy *qLaminarRightProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    Param_TransitionSolver *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    using namespace IDX;
    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();

    int nTransitionEquation = parameters->GetNTransitionEquation();
    Range M(0, nTransitionEquation - 1);

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++iBCRegion)
    {
        StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
        int BCType = structBC->GetBCType();

        if (!(IsWall(BCType) && viscousType > INVISCID))
        {
            continue;
        }

        int *faceDirectionIndex = structBC->GetFaceDirectionIndex();
        int iSurfOfBC = structBC->GetFaceDirection() + 1;

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        if (iSurfOfBC != iSurface)
        {
            continue;
        }

        int il, jl, kl;
        int id, jd, kd;

        GetBCFaceIDX(faceDirectionIndex, id, jd, kd);

        kst = kst + kd;
        ked = ked + kd;
        jst = jst + jd;
        jed = jed + jd;
        ist = ist + id;
        ied = ied + id;

        int di, dj, dk;

        di = ABS(faceDirectionIndex[0]);
        dj = ABS(faceDirectionIndex[1]);
        dk = ABS(faceDirectionIndex[2]);

        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    il = i - di;
                    jl = j - dj;
                    kl = k - dk;
                }
            }
        }
    }
}

void TransitionSolverStr::ComputeInviscidFlux(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, FieldProxy *qLaminarLeftProxy, FieldProxy *qLaminarRightProxy, FieldProxy *fluxProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    NNDFlux(grid, qLeftProxy, qRightProxy, qLaminarLeftProxy, qLaminarRightProxy, fluxProxy, iSurface);
}

void TransitionSolverStr::NNDFlux(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, FieldProxy *qLaminarLeftProxy, FieldProxy *qLaminarRightProxy, FieldProxy *fluxProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &qLeft = qLeftProxy->GetField_STR();
    RDouble4D &qRight = qRightProxy->GetField_STR();

    RDouble4D &qLaminarLeft = qLaminarLeftProxy->GetField_STR();
    RDouble4D &qLaminarRight = qLaminarRightProxy->GetField_STR();

    RDouble4D &flux = fluxProxy->GetField_STR();

    RDouble4D &faceNormalX = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalY = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalZ = *(grid->GetFaceNormalZ());
    RDouble4D &sarea = *(grid->GetFaceArea());
    RDouble4D &svgn = *(grid->GetFaceNormalVelocity());

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    int il1, jl1, kl1;
    grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

    int ist, ied, jst, jed, kst, ked;
    grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

    Range M(0, nTransitionEquation - 1);
    int mst = M.first();
    int med = M.last();

    using namespace IDX;

    RDouble *priml = new RDouble[nTransitionEquation];
    RDouble *primr = new RDouble[nTransitionEquation];

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                int il, jl, kl;
                grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl);

                RDouble xfn = faceNormalX(il, jl, kl, iSurface);
                RDouble yfn = faceNormalY(il, jl, kl, iSurface);
                RDouble zfn = faceNormalZ(il, jl, kl, iSurface);
                RDouble area = sarea(il, jl, kl, iSurface);
                RDouble vgn = svgn(il, jl, kl, iSurface);

                for (int m = mst; m <= med; ++m)
                {
                    priml[m] = qLeft(i, j, k, m);
                    primr[m] = qRight(il, jl, kl, m);
                }

                RDouble ul = qLaminarLeft(i, j, k, IU);
                RDouble vl = qLaminarLeft(i, j, k, IV);
                RDouble wl = qLaminarLeft(i, j, k, IW);

                RDouble ur = qLaminarRight(il, jl, kl, IU);
                RDouble vr = qLaminarRight(il, jl, kl, IV);
                RDouble wr = qLaminarRight(il, jl, kl, IW);

                RDouble vnl = xfn * ul + yfn * vl + zfn * wl - vgn;
                RDouble vnr = xfn * ur + yfn * vr + zfn * wr - vgn;

                RDouble vnla = half * (vnl + ABS(vnl));
                RDouble vnra = half * (vnr - ABS(vnr));

                for (int m = mst; m <= med; ++m)
                {
                    flux(il, jl, kl, m) = (vnla * priml[m] + vnra * primr[m]) * area;
                }
            }
        }
    }
    delete [] priml;
    delete [] primr;
}

void TransitionSolverStr::LoadFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &flux = fluxProxy->GetField_STR();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));

    int il1, jl1, kl1;
    grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    Range M(0, nTransitionEquation - 1);
    int mst = M.first();
    int med = M.last();

    for (int m = mst; m <= med; ++m)
    {
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    int il, jl, kl;
                    grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl);
                    residualTransition(i, j, k, m) -= (flux(il, jl, kl, m) - flux(i, j, k, m));
                }
            }
        }
    }
}

void TransitionSolverStr::ComputeGradientCellCenter(Grid *gridIn)
{
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &faceVectorX = *(strgrid->GetFaceVectorX());
    RDouble4D &faceVectorY = *(strgrid->GetFaceVectorY());
    RDouble4D &faceVectorZ = *(strgrid->GetFaceVectorZ());
    RDouble3D &volume = *(strgrid->GetCellVolume());

    int ndim = GetDim();
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("q_transition"));
    RDouble4D &gradTransitionCellCenterX = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradTransitionCellCenterX"));
    RDouble4D &gradTransitionCellCenterY = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradTransitionCellCenterY"));
    RDouble4D &gradTransitionCellCenterZ = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradTransitionCellCenterZ"));

    gradTransitionCellCenterX = 0.0;
    gradTransitionCellCenterY = 0.0;
    gradTransitionCellCenterZ = 0.0;

    for (int iSurface = 1; iSurface <= ndim; ++iSurface)
    {
        int il1, jl1, kl1;
        strgrid->GetNsurfIndex(il1, jl1, kl1, iSurface);

        int ist = 1;
        int ied = ni - 1 + il1;
        int jst = 1;
        int jed = nj - 1 + jl1;
        int kst = 1;
        int ked = nk - 1 + kl1;

        if (ndim == TWO_D)
        {
            ked = 1;
        }

        //! The gradient of u,v,w.
        for (int m = 0; m < nTransitionEquation; ++m)
        {
            for (int k = kst; k <= ked; ++k)
            {
                for (int j = jst; j <= jed; ++j)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        int il, jl, kl;
                        strgrid->GetLeftCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl);

                        RDouble phis = qTransition(i, j, k, m) + qTransition(il, jl, kl, m);

                        RDouble ddx = phis * faceVectorX(i, j, k, iSurface);
                        RDouble ddy = phis * faceVectorY(i, j, k, iSurface);
                        RDouble ddz = phis * faceVectorZ(i, j, k, iSurface);

                        RDouble &gradTransitionccXLeft = gradTransitionCellCenterX(i, j, k, m);
                        RDouble &gradTransitionccYLeft = gradTransitionCellCenterY(i, j, k, m);
                        RDouble &gradTransitionccZLeft = gradTransitionCellCenterZ(i, j, k, m);
                        gradTransitionccXLeft -= ddx;
                        gradTransitionccYLeft -= ddy;
                        gradTransitionccZLeft -= ddz;

                        RDouble &gradTransitionccXRight = gradTransitionCellCenterX(il, jl, kl, m);
                        RDouble &gradTransitionccYRight = gradTransitionCellCenterY(il, jl, kl, m);
                        RDouble &gradTransitionccZRight = gradTransitionCellCenterZ(il, jl, kl, m);
                        gradTransitionccXRight += ddx;
                        gradTransitionccYRight += ddy;
                        gradTransitionccZRight += ddz;
                    }
                }
            }
        }
    }

    int ist = 1;
    int ied = ni - 1;
    int jst = 1;
    int jed = nj - 1;
    int kst = 1;
    int ked = nk - 1;

    if (ndim == TWO_D) ked = 1;

    for (int m = 0; m < nTransitionEquation; ++m)
    {
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    RDouble oov = half / volume(i, j, k);
                    RDouble &gradTransitionccX = gradTransitionCellCenterX(i, j, k, m);
                    RDouble &gradTransitionccY = gradTransitionCellCenterY(i, j, k, m);
                    RDouble &gradTransitionccZ = gradTransitionCellCenterZ(i, j, k, m);
                    gradTransitionccX *= oov;
                    gradTransitionccY *= oov;
                    gradTransitionccZ *= oov;
                }
            }
        }
    }
}

void TransitionSolverStr::RotateVectorFromInterface(Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    StructGrid *grid = StructGridCast(gridIn);
    int level = grid->GetLevel();
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();
    StructBCSet *structBCSet = finestGrid->GetStructBCSet();
    RDouble4D &rotTransitiongradValueX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTransitiongradValueX"));
    RDouble4D &rotTransitiongradValueY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTransitiongradValueY"));
    RDouble4D &rotTransitiongradValueZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTransitiongradValueZ"));

    rotTransitiongradValueX = 0.0;
    rotTransitiongradValueY = 0.0;
    rotTransitiongradValueZ = 0.0;

    int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    RDouble rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle * PI / 180.0;
    if (nEquation > 0)
    {
        RDouble4D &fieldRecvTransitionY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTransitionCellCenterY"));
        RDouble4D &fieldRecvTransitionZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTransitionCellCenterZ"));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int it, jt, kt;
                finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
                finestGrid->RemapMultigridIJK(level, it, jt, kt);
                int *ibcregions = structBCSet->GetIFaceInfo();
                int ibcregion = ibcregions[iFace];
                StructBC *bcregion = structBCSet->GetBCRegion(ibcregion);
                string bcName = bcregion->GetBCName();
                if (bcName == "Periodic_up")
                {
                    for (int m = 0; m < nEquation; ++m)
                    {
                        rotTransitiongradValueY(it, jt, kt, m) = fieldRecvTransitionY(it, jt, kt, m) * cos(2 * PI - rotationAngle) - fieldRecvTransitionZ(it, jt, kt, m) * sin(2 * PI - rotationAngle);
                        rotTransitiongradValueZ(it, jt, kt, m) = fieldRecvTransitionY(it, jt, kt, m) * sin(2 * PI - rotationAngle) + fieldRecvTransitionZ(it, jt, kt, m) * cos(2 * PI - rotationAngle);
                    }
                }
                else if (bcName == "Periodic_down")
                {
                    for (int m = 0; m < nEquation; ++m)
                    {
                        rotTransitiongradValueY(it, jt, kt, m) = fieldRecvTransitionY(it, jt, kt, m) * cos(rotationAngle) - fieldRecvTransitionZ(it, jt, kt, m) * sin(rotationAngle);
                        rotTransitiongradValueZ(it, jt, kt, m) = fieldRecvTransitionY(it, jt, kt, m) * sin(rotationAngle) + fieldRecvTransitionZ(it, jt, kt, m) * cos(rotationAngle);
                    }
                }

                for (int m = 0; m < nEquation; ++m)
                {
                    fieldRecvTransitionY(it, jt, kt, m) = rotTransitiongradValueY(it, jt, kt, m);
                    fieldRecvTransitionZ(it, jt, kt, m) = rotTransitiongradValueZ(it, jt, kt, m);
                }
            }
        }
    }
}

void TransitionSolverStr::ReconGradVolWeightedwithCorrection(Grid *gridIn, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int nEquation = GetNumberOfEquations();

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));

    RDouble4D &gradTransitionX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionX"));
    RDouble4D &gradTransitionY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionY"));
    RDouble4D &gradTransitionZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionZ"));

    gradTransitionX = 0.0;
    gradTransitionY = 0.0;
    gradTransitionZ = 0.0;

    RDouble4D &gradTransitionCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionCellCenterX"));
    RDouble4D &gradTransitionCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionCellCenterY"));
    RDouble4D &gradTransitionCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionCellCenterZ"));

    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    int il1, jl1, kl1;
    grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

    Range I(1, ni - 1 + il1);
    Range J(1, nj - 1 + jl1);
    Range K(1, nk - 1 + kl1);
    if (nk == 1)
    {
        K.setRange(1, 1);
    }

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                RDouble dx = xcc(i, j, k) - xcc(i - il1, j - jl1, k - kl1);
                RDouble dy = ycc(i, j, k) - ycc(i - il1, j - jl1, k - kl1);
                RDouble dz = zcc(i, j, k) - zcc(i - il1, j - jl1, k - kl1);
                RDouble dr = sqrt(dx * dx + dy * dy + dz * dz);

                dx = dx / dr;
                dy = dy / dr;
                dz = dz / dr;

                for (int m = 0; m < nEquation; ++m)
                {
                    RDouble a = half * (gradTransitionCellCenterX(i, j, k, m) + gradTransitionCellCenterX(i - il1, j - jl1, k - kl1, m));
                    RDouble b = half * (gradTransitionCellCenterY(i, j, k, m) + gradTransitionCellCenterY(i - il1, j - jl1, k - kl1, m));
                    RDouble c = half * (gradTransitionCellCenterZ(i, j, k, m) + gradTransitionCellCenterZ(i - il1, j - jl1, k - kl1, m));
                    RDouble dqdr = (qTransition(i, j, k, m) - qTransition(i - il1, j - jl1, k - kl1, m)) / dr;
                    RDouble d = a * dx + b * dy + c * dz - dqdr;

                    RDouble &gradTranX = gradTransitionX(i, j, k, m);
                    RDouble &gradTranY = gradTransitionY(i, j, k, m);
                    RDouble &gradTranZ = gradTransitionZ(i, j, k, m);
                    gradTranX = a - d * dx;
                    gradTranY = b - d * dy;
                    gradTranZ = c - d * dz;
                }
            }
        }
    }

    GetGradientAtFace_CorrectionAtPhysicalBoundary(grid, iSurface);
}

void TransitionSolverStr::GetGradientAtFace_CorrectionAtPhysicalBoundary(Grid *gridIn, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &xfn = *(grid->GetFaceNormalX());
    RDouble4D &yfn = *(grid->GetFaceNormalY());
    RDouble4D &zfn = *(grid->GetFaceNormalZ());

    int nEquation = GetNumberOfEquations();

    RDouble4D &gradTransitionFaceX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionX"));
    RDouble4D &gradTransitionFaceY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionY"));
    RDouble4D &gradTransitionFaceZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionZ"));

    RDouble4D &gradTransitionCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionCellCenterX"));
    RDouble4D &gradTransitionCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionCellCenterY"));
    RDouble4D &gradTransitionCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionCellCenterZ"));

    StructBCSet *compositeBCRegion = grid->GetStructBCSet();
    int nBCRegion = compositeBCRegion->GetnBCRegion();

    using namespace PHENGLEI;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++iBCRegion)
    {
        StructBC *bcregion = compositeBCRegion->GetBCRegion(iBCRegion);
        int BCType = bcregion->GetBCType();
        int nsurf_bc = bcregion->GetFaceDirection() + 1;

        int ist, ied, jst, jed, kst, ked;
        bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        if (IsInterface(BCType) || nsurf_bc != iSurface) // inner P TO P BC in the current direction
        {
            continue;
        }
        else if (BCType == SYMMETRY) // symmetry BC
        {
            RDouble xfnSign;
            RDouble yfnSign;
            RDouble zfnSign;

            xfnSign = ABS(xfn(ist, jst, kst, iSurface)) - 0.1;
            yfnSign = ABS(yfn(ist, jst, kst, iSurface)) - 0.1;
            zfnSign = ABS(zfn(ist, jst, kst, iSurface)) - 0.1;

            if (xfnSign > 0.0) // for Y-Z symmetry plane
            {
                for (int k = kst; k <= ked; ++k)
                {
                    for (int j = jst; j <= jed; ++j)
                    {
                        for (int i = ist; i <= ied; ++i)
                        {
                            int ibc1, jbc1, kbc1;
                            bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                            for (int m = 0; m < nEquation; ++m)
                            {
                                gradTransitionFaceX(ibc1, jbc1, kbc1, m) = 0.0;
                                gradTransitionFaceY(ibc1, jbc1, kbc1, m) = gradTransitionCellCenterY(i, j, k, m);
                                gradTransitionFaceZ(ibc1, jbc1, kbc1, m) = gradTransitionCellCenterZ(i, j, k, m);
                            }
                        }
                    }
                }
            }

            if (yfnSign > 0.0) // for X-Z symmetry plane
            {
                for (int k = kst; k <= ked; ++k)
                {
                    for (int j = jst; j <= jed; ++j)
                    {
                        for (int i = ist; i <= ied; ++i)
                        {
                            int ibc1, jbc1, kbc1;
                            bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                            for (int m = 0; m < nEquation; ++m)
                            {
                                gradTransitionFaceX(ibc1, jbc1, kbc1, m) = gradTransitionCellCenterX(i, j, k, m);
                                gradTransitionFaceY(ibc1, jbc1, kbc1, m) = 0.0;
                                gradTransitionFaceZ(ibc1, jbc1, kbc1, m) = gradTransitionCellCenterZ(i, j, k, m);
                            }
                        }
                    }
                }
            }

            if (zfnSign > 0.0) // for X-Y symmetry plane
            {
                for (int k = kst; k <= ked; ++k)
                {
                    for (int j = jst; j <= jed; ++j)
                    {
                        for (int i = ist; i <= ied; ++i)
                        {
                            int ibc1, jbc1, kbc1;
                            bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                            for (int m = 0; m < nEquation; ++m)
                            {
                                gradTransitionFaceX(ibc1, jbc1, kbc1, m) = gradTransitionCellCenterX(i, j, k, m);
                                gradTransitionFaceY(ibc1, jbc1, kbc1, m) = gradTransitionCellCenterY(i, j, k, m);
                                gradTransitionFaceZ(ibc1, jbc1, kbc1, m) = 0.0;
                            }
                        }
                    }
                }
            }
        }
        else
        {
            for (int k = kst; k <= ked; ++k)
            {
                for (int j = jst; j <= jed; ++j)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        int ibc1, jbc1, kbc1;
                        bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                        for (int m = 0; m < nEquation; ++m)
                        {
                            gradTransitionFaceX(ibc1, jbc1, kbc1, m) = gradTransitionCellCenterX(i, j, k, m);
                            gradTransitionFaceY(ibc1, jbc1, kbc1, m) = gradTransitionCellCenterY(i, j, k, m);
                            gradTransitionFaceZ(ibc1, jbc1, kbc1, m) = gradTransitionCellCenterZ(i, j, k, m);
                        }
                    }
                }
            }
        }
    }
}

void TransitionSolverStr::ViscousFlux(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    FieldProxy *fluxProxy = CreateFieldProxy(grid);

    int nDim = GetDim();
    for (int iSurface = 1; iSurface <= nDim; ++iSurface)
    {
        ReconGradVolWeightedwithCorrection(gridIn, iSurface);

        CompVisFlux(grid, fluxProxy, iSurface);

        LoadFlux(grid, fluxProxy, iSurface);
    }

    delete fluxProxy;
}

void TransitionSolverStr::CompVisFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &gradientTransitionX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionX"));
    RDouble4D &gradientTransitionY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionY"));
    RDouble4D &gradientTransitionZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTransitionZ"));

    RDouble4D &flux = fluxProxy->GetField_STR();

    RDouble4D &faceNormalX = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalY = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalZ = *(grid->GetFaceNormalZ());
    RDouble4D &area = *(grid->GetFaceArea());

    RDouble3D &viscousLaminar = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTransition = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    Param_TransitionSolver *parameters = GetControlParameters();
    string viscousName = parameters->GetViscousName();

    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble *fvis = new RDouble[nTransitionEquation];

    int il1, jl1, kl1;
    grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

    int ist, ied, jst, jed, kst, ked;
    grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

    using namespace IDX;

    RDouble dct = parameters->Getdct();
    RDouble df = parameters->Getdf();

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                int il, jl, kl;
                grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl);

                RDouble nxs = faceNormalX(il, jl, kl, iSurface);
                RDouble nys = faceNormalY(il, jl, kl, iSurface);
                RDouble nzs = faceNormalZ(il, jl, kl, iSurface);
                RDouble nss = area(il, jl, kl, iSurface);


                RDouble muLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(il, jl, kl));
                RDouble muTransition = half * (viscousTransition(i, j, k) + viscousTransition(il, jl, kl));

                RDouble dkgamadx = gradientTransitionX(il, jl, kl, IGAMA);
                RDouble dkgamady = gradientTransitionY(il, jl, kl, IGAMA);
                RDouble dkgamadz = gradientTransitionZ(il, jl, kl, IGAMA);

                RDouble dkrectdx = gradientTransitionX(il, jl, kl, IRECT);
                RDouble dkrectdy = gradientTransitionY(il, jl, kl, IRECT);
                RDouble dkrectdz = gradientTransitionZ(il, jl, kl, IRECT);

                RDouble mlt3 = muLaminar + muTransition / df;
                RDouble mlt4 = (muLaminar + muTransition) * dct;

                RDouble dkgamadn = nxs * dkgamadx + nys * dkgamady + nzs * dkgamadz;
                RDouble dkrectdn = nxs * dkrectdx + nys * dkrectdy + nzs * dkrectdz;

                fvis[IGAMA] = -oRefReNumber * mlt3 * dkgamadn * nss;
                fvis[IRECT] = -oRefReNumber * mlt4 * dkrectdn * nss;

                for (int m = 0; m < nTransitionEquation; ++m)
                {
                    RDouble &fluxface = flux(il, jl, kl, m);
                    fluxface = fvis[m];
                }
            }
        }
    }

    delete [] fvis;    fvis = NULL;
}

void TransitionSolverStr::ZeroResidualOfSpecialCells(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int isOverset = parameters->GetIsOverLapping();
    if (!isOverset)
    {
        return;
    }

    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("res_transition"));
    StructGrid *grid = PHSPACE::StructGridCast(gridIn);
    int numberOfEquations = GetNumberOfEquations();

    int *iBlank = grid->GetCellTypeContainer();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int cellLabel = 0;
    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                if (iBlank[cellLabel] != GENERIC_COLOR)
                {
                    for (int m = 0; m < numberOfEquations; ++m)
                    {
                        residualTransition(i, j, k, m) = 0.0;
                    }
                }
                cellLabel += 1;
            }
        }
    }
}

void TransitionSolverStr::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &dqTransition = dqProxy->GetField_STR();
    RDouble4D &deltaFlux = LUplusDQ->GetField_STR();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble *dqNeighbor = new RDouble[nTransitionEquation + 1];
    RDouble *prim = new RDouble[nTransitionEquation + 1];
    RDouble *rhs0 = new RDouble[nTransitionEquation + 1];
    RDouble *dflux = new RDouble[nTransitionEquation + 1];
    RDouble *mat_abc = new RDouble[nTransitionEquation + 1];
    RDouble *dqOld = new RDouble[nTransitionEquation + 1];
    RDouble *Ux = new RDouble[nTransitionEquation + 1];

    RDouble4D &spectrumTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_transition"));
    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));
    RDouble5D &matrixTransitionLeft = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_transitionl"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int nDim = GetDim();

    Range M(0, nTransitionEquation - 1);
    int mst = M.first();
    int med = M.last();

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                for (int m = mst; m <= med; ++m)
                {
                    rhs0[m] = 0.0;
                }

                for (int m = mst; m <= med; ++m)
                {
                    //! Back up the old dq, to compute the convergence.
                    //! it is initialized by zero or res in the first sweep,
                    //! then it is updated by backward sweep, that is dq*.
                    dqOld[m] = dqTransition(i, j, k, m);

                    //! Here, the deltaFlux is computed in Upper backward sweep!
                    //! res: b      deltaFlux: Ux, which is computed in backward sweep.
                    //! the dq is not the real dq, now.
                    //! Although the 'dq' changed can not the right one, it dosen't matter, since 
                    //! the following only using the lower neighbor cells, whose 'dq' has been updated.
                    //! It is convenient for precondition transform.
                    RDouble &dq = dqTransition(i, j, k, m);
                    dq = residualTransition(i, j, k, m);

                    Ux[m] = deltaFlux(i, j, k, m);

                    //! Then reset it to zero to store the Lower forward sweep!
                    RDouble &df = deltaFlux(i, j, k, m);
                    df = 0.0;
                }

                for (int iSurface = 1; iSurface <= nDim; ++iSurface)
                {
                    int il, jl, kl, il1, jl1, kl1;
                    grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

                    grid->GetLeftCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl);

                    for (int m = mst; m <= med; ++m)
                    {
                        dqNeighbor[m] = dqTransition(il, jl, kl, m);
                    }

                    for (int m = mst; m <= med; ++m)
                    {
                        mat_abc[m] = matrixTransitionLeft(il, jl, kl, m, iSurface);
                    }

                    Transition_MxDQ(mat_abc, dqNeighbor, nTransitionEquation, dflux);

                    for (int m = 0; m < nTransitionEquation; ++m)
                    {
                        rhs0[m] += dflux[m];
                    }
                }

                for (int m = mst; m <= med; ++m)
                {
                    //! dq = { (b - Ux) - Lx } / D.
                    //! Note: the 'dq' has actually initialized by the residual at the beginning of this function.
                    //! the 'rhs0' is the total Delta Flux of the neighbors.            
                    //! rhs0: Lx    diagonal: D.
                    dqTransition(i, j, k, m) = (dqTransition(i, j, k, m) - Ux[m] - rhs0[m]) / spectrumTransition(i, j, k, m);

                    //! Store the lower forward sweep delta-flux, which will be used in the backward sweep.
                    deltaFlux(i, j, k, m) += rhs0[m];

                    sweepNormal += SQR(dqTransition(i, j, k, m) - dqOld[m]);
                }
            }
        }
    }

    delete [] dqNeighbor;
    delete [] prim;
    delete [] rhs0;
    delete [] dflux;
    delete [] mat_abc;
    delete [] dqOld;
    delete [] Ux;
}

void TransitionSolverStr::SolveLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &dqTransition = dqProxy->GetField_STR();
    RDouble4D &deltaFlux = LUplusDQ->GetField_STR();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble *dqNeighbor = new RDouble[nTransitionEquation + 1];
    RDouble *dqOld = new RDouble[nTransitionEquation + 1];
    RDouble *prim = new RDouble[nTransitionEquation + 1];
    RDouble *rhs0 = new RDouble[nTransitionEquation + 1];
    RDouble *dflux = new RDouble[nTransitionEquation + 1];
    RDouble *mat_abc = new RDouble[nTransitionEquation + 1];
    RDouble *Lx = new RDouble[nTransitionEquation + 1];

    RDouble4D &spectrumTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_transition"));
    RDouble4D &residualTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_transition"));
    RDouble5D &matrixTransitionRight = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_transitionr"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    Range M(0, nTransitionEquation - 1);

    int nDim = GetDim();

    int mst = M.first();
    int med = M.last();

    for (int k = ked; k >= kst; --k)
    {
        for (int j = jed; j >= jst; --j)
        {
            for (int i = ied; i >= ist; --i)
            {
                for (int m = mst; m <= med; ++m)
                {
                    rhs0[m] = 0.0;
                }

                for (int m = mst; m <= med; ++m)
                {
                    //! Back up the old dq, to compute the convergence.
                    //! The 'dq' is dq*, which has been updated in forward.
                    dqOld[m] = dqTransition(i, j, k, m);

                    //! The dq is not the real dq, now.
                    //! it is convenient for precondition transform.
                    dqTransition(i, j, k, m) = residualTransition(i, j, k, m);

                    Lx[m] = deltaFlux(i, j, k, m);

                    deltaFlux(i, j, k, m) = 0.0;
                }

                for (int iSurface = 1; iSurface <= nDim; ++iSurface)
                {
                    int il, jl, kl, il1, jl1, kl1;
                    grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

                    grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl);

                    for (int m = mst; m <= med; ++m)
                    {
                        dqNeighbor[m] = dqTransition(il, jl, kl, m);
                    }

                    for (int m = mst; m <= med; ++m)
                    {
                        mat_abc[m] = matrixTransitionRight(il, jl, kl, m, iSurface);
                    }

                    Transition_MxDQ(mat_abc, dqNeighbor, nTransitionEquation, dflux);

                    for (int m = mst; m <= med; ++m)
                    {
                        rhs0[m] += dflux[m];
                    }
                }

                for (int m = mst; m <= med; ++m)
                {
                    //! Note: the 'dq' has been updated by the forward sweep.
                    //! The 'rhs0' is the total Delta Flux of the neighbors.
                    //! x = {(b - LX) - Ux} / D.
                    //! rhs0: Ux.    diagonal: D.
                    dqTransition(i, j, k, m) = (dqTransition(i, j, k, m) - Lx[m] - rhs0[m]) / spectrumTransition(i, j, k, m);

                    //! Store the upper backward sweep delta-flux, which will be used in the forward sweep.
                    deltaFlux(i, j, k, m) += rhs0[m];

                    sweepNormal += SQR(dqTransition(i, j, k, m) - dqOld[m]);
                }
            }
        }
    }

    delete [] dqNeighbor;
    delete [] dqOld;
    delete [] prim;
    delete [] rhs0;
    delete [] dflux;
    delete [] mat_abc;
    delete [] Lx;
}

void TransitionSolverStr::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &qTransition = qProxy->GetField_STR();
    RDouble4D &dqTransition = dqProxy->GetField_STR();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble4D &qpmv = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    Range M(0, nTransitionEquation - 1);

    using namespace IDX;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                RDouble rho = qpmv(i, j, k, IR);
                RDouble &kgama = qTransition(i, j, k, IGAMA);
                RDouble &krect = qTransition(i, j, k, IRECT);

                RDouble dkgama = dqTransition(i, j, k, IGAMA);
                RDouble dkrect = dqTransition(i, j, k, IRECT);

                kgama += dkgama / rho;
                krect += dkrect / rho;

                kgama = MIN(kgama, 1.00);
                krect = MAX(krect, 20.0);
            }
        }
    }
}


void TransitionSolverStr::Boundary(Grid *gridIn)
{
    using namespace PHENGLEI;
    StructGrid *grid = StructGridCast(gridIn);

    StructBCSet *StructBCSet = grid->GetStructBCSet();
    int nBCRegion = StructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++iBCRegion)
    {
        StructBC *structBC = StructBCSet->GetBCRegion(iBCRegion);
        int BCType = structBC->GetBCType();

        if (IsInterface(BCType))
        {
            continue;
        }
        else if (BCType == EXTRAPOLATION)
        {
            OutFlowBC(grid, structBC);
        }
        else if (IsWall(BCType))
        {
            VisWall(grid, structBC);
        }
        else if (BCType == SYMMETRY)
        {
            SymmetryBC(grid, structBC);
        }
        else if (BCType == FARFIELD)
        {
            FarFieldBC(grid, structBC);
        }
        else if (BCType == INFLOW)
        {
            InFlowBC(grid, structBC);
        }
        else if (BCType == PRESSURE_INLET)
        {
            InFlowBC(grid, structBC);
        }
        else if (BCType == MASS_FLOW_INLET)
        {
            InFlowBC(grid, structBC);
        }
        else if (BCType == OUTFLOW)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == PRESSURE_OUTLET)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == MASS_FLOW_OUTLET)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == POLE || BCType / 10 == POLE)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == OVERSET)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == EXTERNAL_BC)
        {
            OutFlowBC(grid, structBC);
        }
        else
        {
            ostringstream oss;
            oss << "Error : Illegal BCtype ID " << BCType << endl;
            TK_Exit::ExceptionExit(oss);
        }
    }
    CornerPoint(grid);
}

void TransitionSolverStr::CornerPoint(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    FillCornerPoint3D(qTransition, ni, nj, nk, nTransitionEquation);
}

void TransitionSolverStr::InFlowBC(Grid *grid, StructBC *structBC)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    RDouble *freeStreamTransitionVar = parameters->GetFreeStreamTransitionVar();

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int it, jt, kt;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                for (int layer = 1; layer <= 2; ++layer)
                {
                    structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                    for (int m = 0; m < nTransitionEquation; ++m)
                    {
                        qTransition(it, jt, kt, m) = freeStreamTransitionVar[m];
                    }
                }
            }
        }
    }
}

void TransitionSolverStr::OutFlowBC(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int is1, js1, ks1, is2, js2, ks2;
    int it1, jt1, kt1, it2, jt2, kt2;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                structBC->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                RestrictIndex(is2, js2, ks2, ni - 1, nj - 1, nk - 1);

                structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                for (int m = 0; m < nTransitionEquation; ++m)
                {
                    qTransition(it1, jt1, kt1, m) = two * qTransition(is1, js1, ks1, m) - qTransition(is2, js2, ks2, m);
                    qTransition(it2, jt2, kt2, m) = two * qTransition(it1, jt1, kt1, m) - qTransition(is1, js1, ks1, m);
                }
            }
        }
    }
}

void TransitionSolverStr::SymmetryBC(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int is, js, ks, it, jt, kt;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                for (int layer = 1; layer <= 2; ++layer)
                {
                    structBC->GetInsideCellIndex(i, j, k, is, js, ks, layer);
                    RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);
                    structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                    for (int m = 0; m < nTransitionEquation; ++m)
                    {
                        qTransition(it, jt, kt, m) = qTransition(is, js, ks, m);
                    }
                }
            }
        }
    }
}

void TransitionSolverStr::VisWall(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    using namespace IDX;

    int is1, js1, ks1, is2, js2, ks2;
    int it1, jt1, kt1, it2, jt2, kt2;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                structBC->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                RestrictIndex(is2, js2, ks2, ni - 1, nj - 1, nk - 1);

                structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                qTransition(it1, jt1, kt1, IDX::IGAMA) = qTransition(is1, js1, ks1, IDX::IGAMA);
                qTransition(it2, jt2, kt2, IDX::IGAMA) = qTransition(it1, jt1, kt1, IDX::IGAMA);

                qTransition(it1, jt1, kt1, IDX::IRECT) = qTransition(is1, js1, ks1, IDX::IRECT);
                qTransition(it2, jt2, kt2, IDX::IRECT) = qTransition(it1, jt1, kt1, IDX::IRECT);
            }
        }
    }
}

void TransitionSolverStr::FarFieldBC(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &faceNormalX = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalY = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalZ = *(grid->GetFaceNormalZ());

    Param_TransitionSolver *parameters = GetControlParameters();
    RDouble *freeStreamTransitionVar = parameters->GetFreeStreamTransitionVar();

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    RDouble4D &qLaminar = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &gamma = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("gama"));

    int nm = GlobalDataBase::GetIntParaFromDB("nm");
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");

    using namespace IDX;

    RDouble *primitiveVariablesInfinity = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble rFarfield = primitiveVariablesInfinity[IR];
    RDouble uFarfield = primitiveVariablesInfinity[IU];
    RDouble vFarfield = primitiveVariablesInfinity[IV];
    RDouble wFarfield = primitiveVariablesInfinity[IW];
    RDouble pFarfield = primitiveVariablesInfinity[IP];

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int iSurface = structBC->GetFaceDirection() + 1;
    int leftOrRightIndex = structBC->GetFaceLeftOrRightIndex();

    RDouble *prims1 = new RDouble[nm];

    int is1, js1, ks1, is2, js2, ks2;
    int it1, jt1, kt1, it2, jt2, kt2;

    using namespace IDX;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                int in, jn, kn;
                structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);
                RDouble nxs = leftOrRightIndex * faceNormalX(in, jn, kn, iSurface);
                RDouble nys = leftOrRightIndex * faceNormalY(in, jn, kn, iSurface);
                RDouble nzs = leftOrRightIndex * faceNormalZ(in, jn, kn, iSurface);

                structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                structBC->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                RestrictIndex(is2, js2, ks2, ni - 1, nj - 1, nk - 1);

                structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                for (int m = 0; m < nm; ++m)
                {
                    prims1[m] = qLaminar(is1, js1, ks1, m);
                }

                RDouble rin = prims1[IR];
                RDouble uin = prims1[IU];
                RDouble vin = prims1[IV];
                RDouble win = prims1[IW];
                RDouble pin = prims1[IP];

                RDouble vno = nxs * uFarfield + nys * vFarfield + nzs * wFarfield;
                RDouble vni = nxs * uin + nys * vin + nzs * win;
                RDouble vei = sqrt(uin * uin + vin * vin + win * win);

                RDouble gama = gamma(is1, js1, ks1);

                RDouble coo = sqrt(ABS(refGama * pFarfield / rFarfield));
                RDouble cIN = sqrt(ABS(gama * pin / rin));

                if (vei > cIN)
                {
                    if (vni >= 0.0)
                    {
                        for (int m = 0; m < nTransitionEquation; ++m)
                        {
                            qTransition(it1, jt1, kt1, m) = qTransition(is1, js1, ks1, m);
                            qTransition(it2, jt2, kt2, m) = qTransition(is2, js2, ks2, m);
                        }
                    }
                    else
                    {
                        for (int m = 0; m < nTransitionEquation; ++m)
                        {
                            qTransition(it1, jt1, kt1, m) = freeStreamTransitionVar[m];
                            qTransition(it2, jt2, kt2, m) = freeStreamTransitionVar[m];
                        }
                    }
                    continue;
                }

                RDouble riemp = vni + 2.0 * cIN / (gama - 1.0);
                RDouble riemm = vno - 2.0 * coo / (refGama - 1.0);
                RDouble vnb = half * (riemp + riemm);

                if (vnb >= 0.0)
                {
                    //! Exit.
                    for (int m = 0; m < nTransitionEquation; ++m)
                    {
                        qTransition(it1, jt1, kt1, m) = qTransition(is1, js1, ks1, m);
                        qTransition(it2, jt2, kt2, m) = qTransition(is2, js2, ks2, m);
                    }
                }
                else
                {
                    //! Inlet.
                    for (int m = 0; m < nTransitionEquation; ++m)
                    {
                        qTransition(it1, jt1, kt1, m) = freeStreamTransitionVar[m];
                        qTransition(it2, jt2, kt2, m) = freeStreamTransitionVar[m];
                    }
                }
            }
        }
    }
    delete [] prims1;
}

void TransitionSolverStr::UploadInterfaceData(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid *grid = reinterpret_cast<StructGrid *> (GetGrid(level));
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo)
    {
        return;
    }

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble4D *qTransition = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    PHSPACE::UploadInterfaceValue(grid, qTransition, "transition::q", nTransitionEquation);
}

void TransitionSolverStr::DownloadInterfaceData(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid *grid = reinterpret_cast<StructGrid *> (GetGrid(level));
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo)
    {
        return;
    }

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble4D *qTransition = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    PHSPACE::DownloadInterfaceValue(grid, qTransition, "transition::q", nTransitionEquation);

}

void TransitionSolverStr::UploadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid *grid = StructGridCast(GetGrid(level));
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble4D *qTransition = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    PHSPACE::UploadOversetValue(grid, qTransition, "transition::q", nTransitionEquation);
}

void TransitionSolverStr::DownloadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid *grid = StructGridCast(GetGrid(level));
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble4D *qTransition = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_transition"));
    PHSPACE::DownloadOversetValue(grid, qTransition, "transition::q", nTransitionEquation);
}

void TransitionSolverStr::CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
{
    StructGrid *grid = StructGridCast(gridIn);
    int level = grid->GetLevel();
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble4D &fieldSend = fieldProxy->GetField_STR();

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
    {
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
        {
            int iFace = interfaceIndexContainerForSend[iLocalFace];
            int is, js, ks;
            finestGrid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
            finestGrid->RemapMultigridIJK(level, is, js, ks);
            for (int m = 0; m < nEquation; ++m)
            {
                PHWrite(dataContainer, fieldSend(is, js, ks, m));
            }
        }
    }
}

void TransitionSolverStr::DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
{
    StructGrid *grid = StructGridCast(gridIn);
    int level = grid->GetLevel();
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble4D &fieldRecv = fieldProxy->GetField_STR();

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
    {
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
        {
            int iFace = interfaceIndexContainerForReceive[iLocalFace];
            int it, jt, kt;
            finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
            finestGrid->RemapMultigridIJK(level, it, jt, kt);
            for (int m = 0; m < nEquation; ++m)
            {
                PHRead(dataContainer, fieldRecv(it, jt, kt, m));
            }
        }
    }
}
}