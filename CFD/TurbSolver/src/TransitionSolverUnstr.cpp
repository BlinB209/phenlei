#include "TransitionSolverUnstr.h"
#include "Geo_UnstructBC.h"
#include "IO_FileName.h"
#include "MultiGridOperation.h"
#include "PHMatrix.h"
#include "TK_Exit.h"
#include "Param_TransitionSolver.h"
#include "Residual.h"
#include "FieldProxy.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "TK_Log.h"
#include "IO_FileName.h"
#include "Transition.h"
using namespace std;

namespace PHSPACE
{

TransitionSolverUnstr::TransitionSolverUnstr()
{
    gradientTransitionField = 0;
    gradientVelocity = 0;
}

TransitionSolverUnstr::~TransitionSolverUnstr()
{
    DeAllocateGlobalVariables();

    FreeControlParameters();
}

void TransitionSolverUnstr::AllocateGlobalVar(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    if (nTransitionEquation == 0) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfTotalFace = grid->GetNTotalFace();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");
    int nTurboZone = 0;
    if (referenceFrame)
    {
        nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    }

    RDouble **qTransition = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);
    RDouble **residualTransition = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);
    RDouble **spectrumTransition = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);

    grid->UpdateDataPtr("q_transition", qTransition);    //! Unknown variable of transition model.
    grid->UpdateDataPtr("res_transition", residualTransition);    //! Residual or right-hand side.
    grid->UpdateDataPtr("spec_transition", spectrumTransition);    //! Spectral radius of transition model.

    RDouble **gradientTransitionX = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);
    RDouble **gradientTransitionY = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);
    RDouble **gradientTransitionZ = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);

    grid->UpdateDataPtr("gradTransitionX", gradientTransitionX);
    grid->UpdateDataPtr("gradTransitionY", gradientTransitionY);
    grid->UpdateDataPtr("gradTransitionZ", gradientTransitionZ);

    //! Temporary variable for periodic boundary condition.
    RDouble **rotTransitiongradValueX = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);
    RDouble **rotTransitiongradValueY = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);
    RDouble **rotTransitiongradValueZ = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);

    grid->UpdateDataPtr("rotTraTransitionitiongradValueX", rotTransitiongradValueX);
    grid->UpdateDataPtr("rotTraTransitionitiongradValueY", rotTransitiongradValueY);
    grid->UpdateDataPtr("rotTraTransitionitiongradValueZ", rotTransitiongradValueZ);

    //! allocate mixing plane variables.
    //! for multi-row turbomachinery.
    if (nTurboZone > 0)
    {
        int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");
        int nType = 2;

        //! zone with MixingIn and MixingOut.
        int nFlag = 10 * nTurboZone + 2;

        RDouble ***SpanTransFlux      = NewPointer3<RDouble>(nFlag, nTransitionEquation, nSpanSection);
        RDouble ***qTrans_Span        = NewPointer3<RDouble>(nFlag, nTransitionEquation, nSpanSection);
        RDouble ***qTransSpanNeighbor = NewPointer3<RDouble>(nFlag, nTransitionEquation, nSpanSection);
        RDouble ***dqTransSpanIn      = NewPointer3<RDouble>(nFlag, nTransitionEquation, nSpanSection);
        RDouble ***dqTransSpanEx      = NewPointer3<RDouble>(nFlag, nTransitionEquation, nSpanSection);
        RDouble ***dcTransSpan        = NewPointer3<RDouble>(nFlag, nTransitionEquation, nSpanSection);
        InitMixingPlane(SpanTransFlux, nFlag, nTransitionEquation, nSpanSection, zero);
        InitMixingPlane(qTrans_Span, nFlag, nTransitionEquation, nSpanSection, zero);
        InitMixingPlane(qTransSpanNeighbor, nFlag, nTransitionEquation, nSpanSection, zero);
        InitMixingPlane(dqTransSpanIn, nFlag, nTransitionEquation, nSpanSection, zero);
        InitMixingPlane(dqTransSpanEx, nFlag, nTransitionEquation, nSpanSection, zero);
        InitMixingPlane(dcTransSpan, nFlag, nTransitionEquation, nSpanSection, zero);
        grid->UpdateDataPtr("SpanTransFlux", SpanTransFlux);
        grid->UpdateDataPtr("qTrans_Span", qTrans_Span);
        grid->UpdateDataPtr("qTransSpanNeighbor", qTransSpanNeighbor);
        grid->UpdateDataPtr("dqTransSpanIn", dqTransSpanIn);
        grid->UpdateDataPtr("dqTransSpanEx", dqTransSpanEx);
        grid->UpdateDataPtr("dcTransSpan", dcTransSpan);
    }


    int numberOfTotalNode = grid->GetNTotalNode();
    RDouble **qTransitionNode = NewPointer2<RDouble>(nTransitionEquation, numberOfTotalNode);
    grid->UpdateDataPtr("qTransitionNode", qTransitionNode);    //! Unknown variable of transition model at node.

    int *nodeValueSliceTransition = new int[numberOfTotalNode];
    grid->UpdateDataPtr("nodeValueSliceTransition", nodeValueSliceTransition);

    int *nodeValueSliceTransitionTrade = new int[numberOfTotalNode];
    PHSPACE::SetField(nodeValueSliceTransitionTrade, 0, numberOfTotalNode);
    grid->UpdateDataPtr("nodeValueSliceTransitionTrade", nodeValueSliceTransitionTrade);

    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation)
    {
        int numberOfInterPoint = interpointInformation->GetNumberOfInterpoints();
        RDouble **qTransitionInterPoint = NewPointer2<RDouble>(nTransitionEquation, numberOfInterPoint);
        grid->UpdateDataPtr("transition::qInterpoint", qTransitionInterPoint);    //! Unknown variable of transition model at interPoint.
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **qTransitionUnsteadyN1 = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);
        RDouble **qTransitionUnsteadyN2 = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);
        RDouble **residualTransitionUnsteadyN1 = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);
        RDouble **residualTransitionUnsteadyN2 = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble **residualTransitionUnsteadyTemporary = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);

        grid->UpdateDataPtr("q_transition_unsteady_n1", qTransitionUnsteadyN1);
        grid->UpdateDataPtr("q_transition_unsteady_n2", qTransitionUnsteadyN2);
        grid->UpdateDataPtr("res_transition_unsteady_n1", residualTransitionUnsteadyN1);
        grid->UpdateDataPtr("res_transition_unsteady_n2", residualTransitionUnsteadyN2);
        grid->UpdateDataPtr("res_transition_unsteady_tmp", residualTransitionUnsteadyTemporary);

        //! Statistical variables for DES simulation.
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        if (ifStaticsFlowField > 0)
        {
            RDouble **qAverageTransition = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);
            grid->UpdateDataPtr("qAverageTransition", qAverageTransition);
        }
    }

    RDouble **matrixTransitionLeft = NewPointer2<RDouble>(nTransitionEquation, numberOfTotalFace);
    RDouble **matrixTransitionRight = NewPointer2<RDouble>(nTransitionEquation, numberOfTotalFace);

    grid->UpdateDataPtr("mat_transitionl", matrixTransitionLeft);
    grid->UpdateDataPtr("mat_transitionr", matrixTransitionRight);
}

void TransitionSolverUnstr::DeAllocateGlobalVar(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    if (nTransitionEquation == 0) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));
    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));
    RDouble **spectrumTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_transition"));

    DelPointer2(qTransition);
    DelPointer2(residualTransition);
    DelPointer2(spectrumTransition);

    RDouble **gradientTransitionX = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTransitionX"));
    RDouble **gradientTransitionY = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTransitionY"));
    RDouble **gradientTransitionZ = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTransitionZ"));

    DelPointer2(gradientTransitionX);
    DelPointer2(gradientTransitionY);
    DelPointer2(gradientTransitionZ);

    RDouble **rotTransitiongradValueX = reinterpret_cast<RDouble **> (grid->GetDataPtr("rotTransitiongradValueX"));
    RDouble **rotTransitiongradValueY = reinterpret_cast<RDouble **> (grid->GetDataPtr("rotTransitiongradValueY"));
    RDouble **rotTransitiongradValueZ = reinterpret_cast<RDouble **> (grid->GetDataPtr("rotTransitiongradValueZ"));

    DelPointer2(rotTransitiongradValueX);
    DelPointer2(rotTransitiongradValueY);
    DelPointer2(rotTransitiongradValueZ);

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");
    int nTurboZone = 0;
    if (referenceFrame)
    {
        nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    }

    //! deallocate mixing plane array.
    if (nTurboZone > 0)
    {
        int nType = 2;
        int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");

        RDouble ***SpanTransFlux      = reinterpret_cast <RDouble ***> (grid->GetDataPtr("SpanTransFlux"));
        RDouble ***qTrans_Span        = reinterpret_cast <RDouble ***> (grid->GetDataPtr("qTrans_Span"));
        RDouble ***qTransSpanNeighbor = reinterpret_cast <RDouble ***> (grid->GetDataPtr("qTransSpanNeighbor"));
        RDouble ***dqTransSpanIn      = reinterpret_cast <RDouble ***> (grid->GetDataPtr("dqTransSpanIn"));
        RDouble ***dqTransSpanEx      = reinterpret_cast <RDouble ***> (grid->GetDataPtr("dqTransSpanEx"));
        RDouble ***dcTransSpan        = reinterpret_cast <RDouble ***> (grid->GetDataPtr("dcTransSpan"));

        DelPointer3(SpanTransFlux);
        DelPointer3(qTrans_Span);
        DelPointer3(qTransSpanNeighbor);
        DelPointer3(dqTransSpanIn);
        DelPointer3(dqTransSpanEx);
        DelPointer3(dcTransSpan);
    }

    RDouble **qTransitionNode = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTransitionNode"));
    DelPointer2(qTransitionNode);

    int *nodeValueSliceTransition = reinterpret_cast<int *> (grid->GetDataPtr("nodeValueSliceTransition"));
    int *nodeValueSliceTransitionTrade = reinterpret_cast<int *> (grid->GetDataPtr("nodeValueSliceTransitionTrade"));

    delete [] nodeValueSliceTransition;
    delete [] nodeValueSliceTransitionTrade;

    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation)
    {
        RDouble **qTransitionInterPoint = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTransitionInterPoint"));
        DelPointer2(qTransitionInterPoint);
    }

    int isUnsteady = parameters->GetIsUnsteady();

    if (isUnsteady)
    {
        RDouble **qTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n1"));
        RDouble **qTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n2"));
        RDouble **residualTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n1"));
        RDouble **residualTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n2"));
        RDouble **residualTransitionUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_tmp"));

        DelPointer2(qTransitionUnsteadyN1);
        DelPointer2(qTransitionUnsteadyN2);
        DelPointer2(residualTransitionUnsteadyN1);
        DelPointer2(residualTransitionUnsteadyN2);
        DelPointer2(residualTransitionUnsteadyTemporary);

        //! Statistical variables for DES simulation.
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        if (ifStaticsFlowField > 0)
        {
            RDouble **qAverageTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTransition"));
            DelPointer2(qAverageTransition);
        }
    }

    RDouble **matrixTransitionLeft = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_transitionl"));
    RDouble **matrixTransitionRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_transitionr"));

    DelPointer2(matrixTransitionLeft);
    DelPointer2(matrixTransitionRight);
}

bool TransitionSolverUnstr::JudgeIfRestart()
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

bool TransitionSolverUnstr::JudgeIfReadAverage()
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

void TransitionSolverUnstr::DumpRestartH5(ActionKey *actkey)
{
    using namespace PHMPI;
    int currentProcessorID = GetCurrentProcessorID();

    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    int version = 1;

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    int outerIterationStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (currentProcessorID == GetServerProcessorID())
    {
        WriteData(actkey->filepos, &version, "Version");
        WriteData(actkey->filepos, &outerIterationStep, "outnstep");

        if (IsNeedStatistics())
        {
            int nStatisticalStep = 0;
            GlobalDataBase::GetData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            WriteData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
        }
    }

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));
    double **qTemporary = NewPointer2<double>(nTransitionEquation, numberOfTotal);

    for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(qTransition[iEquation][iCell]);
        }
    }

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    WriteData(grploc, qTemporary[0], "q_transition");
    WriteData(grploc, &numberOfTotalCell, "nTotalCell");

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        H5Gclose(grploc);
        DelPointer2(qTemporary);

        return;
    }

    RDouble **qTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n1"));
    for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(qTransitionUnsteadyN1[iEquation][iCell]);
        }
    }

    WriteData(grploc, qTemporary[0], "q_transition_unsteady_n1");

    RDouble **qTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n2"));
    for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(qTransitionUnsteadyN2[iEquation][iCell]);
        }
    }

    WriteData(grploc, qTemporary[0], "q_transition_unsteady_n2");

    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));
    for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(residualTransition[iEquation][iCell]);
        }
    }

    WriteData(grploc, qTemporary[0], "res_transition");

    RDouble **residualTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n1"));
    for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(residualTransitionUnsteadyN1[iEquation][iCell]);
        }
    }

    WriteData(grploc, qTemporary[0], "res_transition_unsteady_n1");

    RDouble **residualTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n2"));
    for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(residualTransitionUnsteadyN2[iEquation][iCell]);
        }
    }

    WriteData(grploc, qTemporary[0], "res_transition_unsteady_n2");

    if (IsNeedStatistics())
    {
        RDouble **qAverageTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTransition"));
        for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
            {
                qTemporary[iEquation][iCell] = static_cast<double>(qAverageTransition[iEquation][iCell]);
            }
        }

        WriteData(grploc, qTemporary[0], "qAverageTransition");
    }

    H5Gclose(grploc);
    DelPointer2(qTemporary);
}

void TransitionSolverUnstr::ReadRestartH5(ActionKey *actkey)
{
    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));
    double **qTemporary = NewPointer2<double>(nTransitionEquation, numberOfTotal);

    int outerIterationStepOfTransition = 0;
    ReadData(actkey->filepos, &outerIterationStepOfTransition, "outnstep");

    int outerIterationStepOfNS = GlobalDataBase::GetIntParaFromDB("outnstep");

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    int nTotalCellRestart = 0;
    ReadData(grploc, &nTotalCellRestart, "nTotalCell");
    if (nTotalCellRestart != numberOfTotalCell)
    {
        ostringstream erroeInfo;
        erroeInfo << " Error: the cell number in transition.dat is not equal to the cell number in grid file !" << endl;
        TK_Exit::ExceptionExit(erroeInfo.str());
    }

    ReadData(grploc, qTemporary[0], "q_transition");
    for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTransition[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **qTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n1"));
        RDouble **qTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n2"));

        RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));
        RDouble **residualTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n1"));
        RDouble **residualTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n2"));
        RDouble **qAverageTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTransition"));

        //! Start from steady flow field and reset the outer step when start from steady flow.
        int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();
        if (ifStartFromSteadyResults)
        {
            PrintToWindow("Restart from steady Transition flow field, reset outer step to be zero!\n");

            outerIterationStepOfTransition = 0;
            GlobalDataBase::UpdateData("outnstep", &outerIterationStepOfTransition, PHINT, 1);

            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
                {
                    qTransitionUnsteadyN1[m][iCell] = qTransition[m][iCell];
                    qTransitionUnsteadyN2[m][iCell] = qTransition[m][iCell];

                    residualTransition[m][iCell] = 0.0;
                    residualTransitionUnsteadyN1[m][iCell] = 0.0;
                    residualTransitionUnsteadyN2[m][iCell] = 0.0;
                }
            }
        }
        else
        {
            ReadData(grploc, qTemporary[0], "q_transition_unsteady_n1");
            for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    qTransitionUnsteadyN1[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                }
            }

            ReadData(grploc, qTemporary[0], "q_transition_unsteady_n2");
            for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    qTransitionUnsteadyN2[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                }
            }

            ReadData(grploc, qTemporary[0], "res_transition");
            for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    residualTransition[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                }
            }

            ReadData(grploc, qTemporary[0], "res_transition_unsteady_n1");
            for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    residualTransitionUnsteadyN1[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                }
            }

            ReadData(grploc, qTemporary[0], "res_transition_unsteady_n2");
            for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    residualTransitionUnsteadyN2[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                }
            }
        }

        int nStatisticalStep = 0;
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
        bool isReadAverageFlow = ifStaticsFlowField && (outerIterationStepOfTransition >= startStatisticStep);
        if (isReadAverageFlow)
        {
            if (ifStartFromSteadyResults)
            {
                for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
                {
                    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                    {
                        qAverageTransition[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                    }
                }

                GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            }
            else
            {
                ReadData(grploc, qTemporary[0], "qAverageTransition");
                for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
                {
                    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                    {
                        qAverageTransition[iEquation][iCell] = qTemporary[iEquation][iCell];
                    }
                }

                ReadData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
                GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            }
        }
    }

    H5Gclose(grploc);
    DelPointer2(qTemporary);

    CompareOutStepOfFlowfieldFile(outerIterationStepOfNS, outerIterationStepOfTransition);
}

void TransitionSolverUnstr::InitFlowAsRestart()
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    if (nTransitionEquation == 0) return;

    UnstructGrid *grid = UnstructGridCast(GetGrid(0));

    //int outerIterationStep = 0;
    //GlobalDataBase::UpdateData("outnstep", &outerIterationStep, PHINT, 1);

    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    RDouble *freeStreamTransitionVar = parameters->GetFreeStreamTransitionVar();
    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTransition[m][iCell] = freeStreamTransitionVar[m];
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    RDouble **qTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n1"));
    RDouble **qTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n2"));

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTransitionUnsteadyN1[m][iCell] = qTransition[m][iCell];
            qTransitionUnsteadyN2[m][iCell] = qTransition[m][iCell];
        }
    }

    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));
    RDouble **residualTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n1"));
    RDouble **residualTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n2"));
    RDouble **residualTransitionUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_tmp"));

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTransition[m][iCell] = 0.0;
            residualTransitionUnsteadyN1[m][iCell] = 0.0;
            residualTransitionUnsteadyN2[m][iCell] = 0.0;
            residualTransitionUnsteadyTemporary[m][iCell] = 0.0;
        }
    }
}

void TransitionSolverUnstr::ZeroResiduals(Grid *gridIn)
{
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(gridIn);

    int numberOfGeneralCells = grid->GetNTotalCell() + grid->GetNBoundFace();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **residualTransition = reinterpret_cast<RDouble **>(grid->GetDataPtr("res_transition"));

    for (int iEquation = 0; iEquation < nTransitionEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfGeneralCells; ++ iCell)
        {
            residualTransition[iEquation][iCell] = 0.0;
        }
    }
    return;
}

void TransitionSolverUnstr::InitSpectrum(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    if (nTransitionEquation == 0) return;

    RDouble transitionCFLScale = parameters->GetTransitionCFLScale();

    RDouble *vol = grid->GetCellVolume();
    RDouble *dt = reinterpret_cast<RDouble *> (grid->GetDataPtr("dt"));
    RDouble **spectrumTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_transition"));

    RDouble dualTimeSpectrumC1 = zero;
    RDouble dualTimeSpectrumC2 = one;

    //! If flow is unsteady, need to count the contribution of unsteady.
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble dualTimeCoefficient[7];
        const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
        dualTimeSpectrumC1 = - dualTimeCoefficient[3];
        dualTimeSpectrumC2 = dualTimeCoefficient[6];
    }

    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
    {
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            spectrumTransition[m][iCell] = 0.0;
        }
    }

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            spectrumTransition[m][iCell] = dualTimeSpectrumC1 * vol[iCell] + dualTimeSpectrumC2 / (transitionCFLScale * dt[iCell] + SMALL);
        }
    }
}

void TransitionSolverUnstr::CorrectFineGrid(Grid *fineGrid, Grid *coarseGrid)
{
    int numberOfTotalCellOnFineGrid = UnstructGridCast(fineGrid)->GetNTotalCell();
    int *cell2CoarseGridCell = UnstructGridCast(fineGrid)->GetCell2CoarseGridCell();

    RDouble **qTransitionFineGrid = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("q_transition"));
    RDouble **qTransitionCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("q_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotalCellOnFineGrid; ++ iCell)
        {
            qTransitionFineGrid[m][iCell] += qTransitionCoarseGrid[m][cell2CoarseGridCell[iCell]];
        }
    }
}

void TransitionSolverUnstr::InterpolatFineGrid(Grid *fineGridIn, Grid *coarseGridIn)
{
    UnstructGrid *fineGrid = UnstructGridCast(fineGridIn);
    UnstructGrid *coarseGrid = UnstructGridCast(coarseGridIn);

    RDouble **qTransitionFineGrid = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("q_transition"));
    RDouble **qTransitionCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("q_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        InterpolatQ(fineGrid, qTransitionFineGrid[m], coarseGrid, qTransitionCoarseGrid[m]);
    }
}

void TransitionSolverUnstr::PutCorrectionBack(Grid *gridIn, FieldProxy *qProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble **q = qProxy->GetField_UNS();
    RDouble **qOld = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            qOld[m][iCell] += q[m][iCell];
        }
    }
}

void TransitionSolverUnstr::RestrictDefect(Grid *fineGridIn, Grid *coarseGridIn)
{
    UnstructGrid *fineGrid = UnstructGridCast(fineGridIn);
    UnstructGrid *coarseGrid = UnstructGridCast(coarseGridIn);

    int numberOfTotalCellOnFineGrid = fineGrid->GetNTotalCell();
    int *cell2CoarseGridCell = fineGrid->GetCell2CoarseGridCell();

    ZeroResiduals(coarseGrid);

    RDouble **residualTransitionCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("res_transition"));
    RDouble **residualTransitionFineGrid = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("res_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotalCellOnFineGrid; ++ iCell)
        {
            residualTransitionCoarseGrid[m][cell2CoarseGridCell[iCell]] -= residualTransitionFineGrid[m][iCell];
        }
    }
}

void TransitionSolverUnstr::LoadQ(Grid *gridIn, FieldProxy *qProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble **qqTransition = qProxy->GetField_UNS();

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            qqTransition[m][iCell] = qTransition[m][iCell];
        }
    }
}

void TransitionSolverUnstr::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble **qPrimitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTransition = qProxy->GetField_UNS();
    RDouble **dqTransition = dqProxy->GetField_UNS();

    using namespace IDX;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        RDouble rho = qPrimitiveVariable[IR][iCell];
        RDouble &kgama = qTransition[IGAMA][iCell];
        RDouble &krect = qTransition[IRECT][iCell];

        RDouble dkgama = dqTransition[IGAMA][iCell];
        RDouble dkrect = dqTransition[IRECT][iCell];

        kgama += dkgama / rho;
        krect += dkrect / rho;

        kgama = MIN(kgama, 1.00);
        krect = MAX(krect, 20.0);
    }
}

void TransitionSolverUnstr::SMoothTransitionPoint(Grid *gridIn, RDouble *primitiveVariable, int i)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    int le, re, ie, face;
    RDouble volumeSum = zero;
    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        primitiveVariable[m] = zero;
    }

    for (int j = 0; j < faceNumberOfEachCell[i]; ++ j)
    {
        face = cell2Face[i][j];
        le = leftCellOfFace[face];
        re = rightCellOfFace[face];

        if (le != i)
        {
            ie = le;
        }
        else
        {
            ie = re;
        }

        if (ie >= numberOfTotalCell) continue;

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            primitiveVariable[m] += qTransition[m][ie];
        }
        volumeSum += 1.0;
    }

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        primitiveVariable[m] /= volumeSum;
    }
}

void TransitionSolverUnstr::SMoothTransition(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *volume = grid->GetCellVolume();

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **tmporary = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
    {
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            qTransition[m][iCell] *= volume[iCell];
            tmporary[m][iCell] = qTransition[m][iCell];
        }
    }

    RDouble volumeSum;

    int le, re, ie, face;
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        volumeSum = volume[iCell];
        for (int j = 0; j < faceNumberOfEachCell[iCell]; ++ j)
        {
            face = cell2Face[iCell][j];
            le = leftCellOfFace[face];
            re = rightCellOfFace[face];

            if (le != iCell)
            {
                ie = le;
            }
            else
            {
                ie = re;
            }

            if (ie >= numberOfTotalCell) continue;

            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                qTransition[m][iCell] += tmporary[m][ie];
            }
            volumeSum += volume[ie];
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            qTransition[m][iCell] /= volumeSum;
        }
    }

    DelPointer2(tmporary);
}

void TransitionSolverUnstr::RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coefficient)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble **dqTransition = dqProxy->GetField_UNS();

    RDouble *dt = reinterpret_cast<RDouble *> (grid->GetDataPtr("dt"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            dqTransition[m][iCell] = dt[iCell] * coefficient * dqTransition[m][iCell];
        }
    }
}

void TransitionSolverUnstr::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **dq = dqProxy->GetField_UNS();
    int *iBlank = grid->GetBlankIndex();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble *df = new RDouble[nTransitionEquation];
    RDouble *rhs0 = new RDouble[nTransitionEquation];
    RDouble *dqOld = new RDouble[nTransitionEquation];

    RDouble **spectrumTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_transition"));
    RDouble **matrixTransitionLeft = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_transitionl"));
    RDouble **matrixTransitionRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_transitionr"));
    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();
    RDouble *matrixUx = new RDouble[nTransitionEquation];

    //! Now the Forward Sweep
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        if (iBlank[iCell] <= 0)
        {
            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                dq[m][iCell] = 0.0;
            }
            continue;
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            rhs0[m] = 0.0;
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            //! Back up the old dq, to compute the convergency.
            dqOld[m] = dq[m][iCell];

            //! Here, the deltaFlux is computed in Upper backward sweep!
            //! res: b      deltaFlux: Ux, which is computed in backward sweep.
            //! the dq is not the real dq, now.
            //! Although the 'dq' changed can not the right one, it dosen't matter, since 
            //! the following only using the lower neighbor cells, whose 'dq' has been updated.
            //! It is convenient for precondition transform.
            dq[m][iCell] = residualTransition[m][iCell];
            matrixUx[m] = deltaFlux[m][iCell];

            //! Then reset it to zero to store the Lower forward sweep!
            deltaFlux[m][iCell] = 0.0;
        }

        for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++ jFace)
        {
            int face = cell2Face[iCell][jFace];
            int le = leftCellOfFace[face];
            int re = rightCellOfFace[face];

            //! One of le and re must be cell itself.
            if (le > iCell || re > iCell) continue;

            if (re == iCell)
            {
                //! Now its neighboring cell belongs to lower triangular
                SWAP(le, re);

                for (int m = 0; m < nTransitionEquation; ++ m)
                {
                    df[m] = matrixTransitionLeft[m][face] * dq[m][re];
                }
            }
            else
            {
                for (int m = 0; m < nTransitionEquation; ++ m)
                {
                    df[m] = matrixTransitionRight[m][face] * dq[m][re];
                }
            }

            //! Add Flux together
            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            //! dq = { (b - Ux) - Lx } / D.
            //! Note: the 'dq' has actually initialized by the residual at the beginning of this function.
            //! the 'rhs0' is the total Delta Flux of the neighbors.            
            //! rhs0: Lx    diagonal: D.
            dq[m][iCell] = (dq[m][iCell] - matrixUx[m] - rhs0[m]) / spectrumTransition[m][iCell];

            //! Store the lower forward sweep delta-flux, which will be used in the backward sweep.
            deltaFlux[m][iCell] += rhs0[m];

            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    }

    delete [] rhs0;
    delete [] df;
    delete [] dqOld;
    delete [] matrixUx;
}

void TransitionSolverUnstr::SolveLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **dq = dqProxy->GetField_UNS();
    int *iBlank = grid->GetBlankIndex();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble *df = new RDouble[nTransitionEquation];
    RDouble *rhs0 = new RDouble[nTransitionEquation];
    RDouble *dqOld = new RDouble[nTransitionEquation];

    RDouble **spectrumTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_transition"));
    RDouble **matrixTransitionLeft = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_transitionl"));
    RDouble **matrixTransitionRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_transitionr"));
    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();
    RDouble *matrixLx = new RDouble[nTransitionEquation];

    // Backward Sweep
    for (int iCell = numberOfTotalCell - 1; iCell >= 0; -- iCell)
    {
        if (iBlank[iCell] <= 0)
        {
            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                dq[m][iCell] = 0.0;
            }
            continue;
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            rhs0[m] = 0.0;
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            //! Back up the old dq, to compute the convergence.
            //! the 'dq' is dq*, which has been updated in forward.
            dqOld[m] = dq[m][iCell];

            //! the dq is not the real dq, now.
            //! it is convenient for precondition transform.
            dq[m][iCell] = residualTransition[m][iCell];

            matrixLx[m] = deltaFlux[m][iCell];

            deltaFlux[m][iCell] = 0.0;
        }

        for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++ jFace)
        {
            int face = cell2Face[iCell][jFace];
            int le = leftCellOfFace[face];
            int re = rightCellOfFace[face];

            //! One of le and re must be cell itself.
            if (le < iCell || re < iCell) continue;

            if (re == iCell)
            {
                //! Now its neighboring cell belongs to upper triangular
                SWAP(le, re);

                for (int m = 0; m < nTransitionEquation; ++ m)
                {
                    df[m] = matrixTransitionLeft[m][face] * dq[m][re];
                }
            }
            else
            {
                for (int m = 0; m < nTransitionEquation; ++ m)
                {
                    df[m] = matrixTransitionRight[m][face] * dq[m][re];
                }
            }

            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            //! Note: the 'dq' has been updated by the forward sweep.
            //! the 'rhs0' is the total Delta Flux of the neighbors.
            //! x = {(b - LX) - Ux} / D.
            //! rhs0: Ux.    diagonal: D.
            dq[m][iCell] = (dq[m][iCell] - matrixLx[m] - rhs0[m]) / spectrumTransition[m][iCell];

            //! Store the upper backward sweep delta-flux, which will be used in the forward sweep.
            deltaFlux[m][iCell] += rhs0[m];

            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    }

    delete [] rhs0;
    delete [] df;
    delete [] dqOld;
    delete [] matrixLx;
}

// Bell 20130401 add
void TransitionSolverUnstr::SetGhostDQLUSGS(Grid *gridIn, RDouble **dq)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    if (nTransitionEquation != 1)
    {
        cout << "SetGhostDQLUSGS function is turned off !\n";
        return;
    }

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *rightCellOfFace = grid->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (IsInterface(bcType))
        {
            continue;
        }

        if (IsWall(bcType))
        {
            SetWallBCGhostDQLUSGS(grid, bcRegion, dq);
        }
        //FARFIELD, SYMMETRY, INFLOW and OUTFLOW can use the same module.
        else if (bcType == PHENGLEI::FARFIELD || bcType == PHENGLEI::SYMMETRY
            || bcType == PHENGLEI::INFLOW || bcType == PHENGLEI::OUTFLOW)
        {
            SetFarfieldBCGhostDQLUSGS(grid, bcRegion, dq);
        }
        else
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int re = rightCellOfFace[iFace];
                for (int m = 0; m < nTransitionEquation; ++ m)
                {
                    dq[m][re] = 0.0;
                }
            }
        }
    }
}

void TransitionSolverUnstr::SetWallBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        int iFace = *iter;
        int re = rightCellOfFace[iFace];
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            dq[m][re] = 0.0;
        }
    }
}

void TransitionSolverUnstr::SetFarfieldBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            dq[m][re] = dq[m][le];
        }
    }
}

void TransitionSolverUnstr::LoadResiduals(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **rhs = rhsProxy->GetField_UNS();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTransition[m][iCell] = - rhs[m][iCell];
        }
    }
}

RDouble TransitionSolverUnstr::UnsteadyConvergence(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return zero;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble **qPrimitive = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));
    RDouble **qTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n1"));
    RDouble **qTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n2"));
    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));

    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble *primitive0 = new RDouble[nTransitionEquation];
    RDouble *primitive1 = new RDouble[nTransitionEquation];
    RDouble *primitive2 = new RDouble[nTransitionEquation];

    RDouble *qConservative0 = new RDouble[nTransitionEquation];
    RDouble *qConservative1 = new RDouble[nTransitionEquation];
    RDouble *qConservative2 = new RDouble[nTransitionEquation];

    RDouble sum1 = zero;
    RDouble sum2 = zero;

    using namespace IDX;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            primitive0[m] = qTransition[m][iCell];
            primitive1[m] = qTransitionUnsteadyN1[m][iCell];
            primitive2[m] = qTransitionUnsteadyN2[m][iCell];
        }

        RDouble coefficient = 1.0;

        if (nTransitionEquation >= 2)
        {
            coefficient = qPrimitive[IR][iCell];
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            qConservative0[m] = coefficient * primitive0[m];
            qConservative1[m] = coefficient * primitive1[m];
            qConservative2[m] = coefficient * primitive2[m];
        }

        //! nl is needed.
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            //! Here, res must be transformed into dqthat is, res is not equal to rhs but dq.
            RDouble dqp = residualTransition[m][iCell];        // qn+1,p+1 - qn+1,p
            RDouble dqn = qConservative0[m] - qConservative1[m];  // qn+1,p+1 - qn
            sum1 += dqp * dqp;
            sum2 += dqn * dqn;
        }
    }

    delete [] primitive0;
    delete [] primitive1;
    delete [] primitive2;
    delete [] qConservative0;
    delete [] qConservative1;
    delete [] qConservative2;

    RDouble cvg = sqrt(ABS(sum1 / (sum2 + SMALL)));
    return cvg;
}

void TransitionSolverUnstr::UpdateUnsteadyFlow(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));
    RDouble **qTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n1"));
    RDouble **qTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n2"));

    RDouble **residualTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n1"));
    RDouble **residualTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n2"));
    RDouble **residualTransitionUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_tmp"));

    int nTransitionEquation = parameters->GetNTransitionEquation();

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTransitionUnsteadyN2[m][iCell] = qTransitionUnsteadyN1[m][iCell];
        }
    }

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTransitionUnsteadyN1[m][iCell] = qTransition[m][iCell];
        }
    }

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTransitionUnsteadyN2[m][iCell] = residualTransitionUnsteadyN1[m][iCell];
        }
    }

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTransitionUnsteadyN1[m][iCell] = residualTransitionUnsteadyTemporary[m][iCell]; //! Here the current outerIterationStep is over, the value of the stored resTmp should be assigned to resn1 for the next outerIterationStep. It should be noticed that resTmp only contain the inviscid and viscous flux.
        }
    }

    //! Statistical variables for DES simulation.
    int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
    if (ifStaticsFlowField > 0)
    {
        int outerIterationStep = GlobalDataBase::GetIntParaFromDB("outnstep");
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");

        if (outerIterationStep >= startStatisticStep)
        {
            int numberOfStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
            RDouble statisticalTimePeriod = GlobalDataBase::GetDoubleParaFromDB("statisticalTimePeriod");
            RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

            RDouble c1 = 1.0 / (numberOfStatisticalStep * 1.0);
            if (statisticalTimePeriod > 0.0)
            {
                c1 = MAX(c1, physicalTimeStep / statisticalTimePeriod);
            }
            RDouble c2 = 1.0 - c1;

            RDouble **qAverageTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTransition"));
            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    qAverageTransition[m][iCell] = c2 * qAverageTransition[m][iCell] + c1 * qTransition[m][iCell];
                }
            }
        }
    }
}

void TransitionSolverUnstr::DualTimeSource(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble **qPrimitive = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));
    RDouble **qTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n1"));
    RDouble **qTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition_unsteady_n2"));

    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));
    RDouble **residualTransitionUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n1"));
    RDouble **residualTransitionUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_n2"));
    RDouble **residualTransitionUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition_unsteady_tmp"));

    RDouble *volume = grid->GetCellVolume();
    RDouble *volumeN1 = grid->GetCellVolume();
    RDouble *volumeN2 = grid->GetCellVolume();

    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble *primitive0 = new RDouble[nTransitionEquation];
    RDouble *primitive1 = new RDouble[nTransitionEquation];
    RDouble *primitive2 = new RDouble[nTransitionEquation];

    RDouble *qConservative0 = new RDouble[nTransitionEquation];
    RDouble *qConservative1 = new RDouble[nTransitionEquation];
    RDouble *qConservative2 = new RDouble[nTransitionEquation];

    //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            residualTransitionUnsteadyTemporary[m][iCell] = residualTransition[m][iCell];
        }
    }

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

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            primitive0[m] = qTransition[m][iCell];
            primitive1[m] = qTransitionUnsteadyN1[m][iCell];
            primitive2[m] = qTransitionUnsteadyN2[m][iCell];
        }

        RDouble coefficient = 1.0;

        if (nTransitionEquation >= 2)
        {
            coefficient = qPrimitive[IR][iCell];
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            qConservative0[m] = coefficient * primitive0[m];
            qConservative1[m] = coefficient * primitive1[m];
            qConservative2[m] = coefficient * primitive2[m];
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            RDouble dualSrcRes = dualTimeResC1 * residualTransition[m][iCell] +
                dualTimeResC2 * residualTransitionUnsteadyN1[m][iCell] +
                dualTimeResC3 * residualTransitionUnsteadyN2[m][iCell];

            RDouble dualSrcQ = dualTimeQC1 * qConservative0[m] * volume[iCell] +
                dualTimeQC2 * qConservative1[m] * volumeN1[iCell] +
                dualTimeQC3 * qConservative2[m] * volumeN2[iCell];

            residualTransition[m][iCell] = dualSrcRes + dualSrcQ;
        }
    }

    delete [] primitive0;
    delete [] primitive1;
    delete [] primitive2;
    delete [] qConservative0;
    delete [] qConservative1;
    delete [] qConservative2;
}

void TransitionSolverUnstr::FillField(Grid *gridIn, FieldProxy *fieldProxy, RDouble value)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **field = fieldProxy->GetField_UNS();

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            field[m][iCell] = value;
        }
    }
}

void TransitionSolverUnstr::FillField(Grid *gridIn, FieldProxy *field1Proxy, FieldProxy *field2Proxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **field1 = field1Proxy->GetField_UNS();
    RDouble **field2 = field2Proxy->GetField_UNS();

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            field1[m][iCell] = field2[m][iCell];
        }
    }
}

FieldProxy *TransitionSolverUnstr::CreateFieldProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **field = NewPointer2<RDouble>(nTransitionEquation, numberOfTotal);

    FieldProxy *fieldProxy = new FieldProxy();

    fieldProxy->SetField_UNS(field, true);

    return fieldProxy;
}

FieldProxy *TransitionSolverUnstr::GetFieldProxy(Grid *gridIn, const string &fieldName)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **field = reinterpret_cast<RDouble **> (grid->GetDataPtr(fieldName));

    FieldProxy *fieldProxy = new FieldProxy();

    fieldProxy->SetField_UNS(field);

    return fieldProxy;
}

FieldProxy *TransitionSolverUnstr::GetResidualProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));

    FieldProxy *residualTransitionProxy = new FieldProxy();

    residualTransitionProxy->SetField_UNS(residualTransition);

    return residualTransitionProxy;
}

void TransitionSolverUnstr::RecoverResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));

    if (grid->GetLevel() != 0)
    {
        RDouble **rhs = rhsProxy->GetField_UNS();
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
            {
                residualTransition[m][iCell] = rhs[m][iCell];
            }
        }
    }
}

void TransitionSolverUnstr::StoreRhsByResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));
    RDouble **rhs = rhsProxy->GetField_UNS();

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            rhs[m][iCell] = residualTransition[m][iCell];
        }
    }
}

void TransitionSolverUnstr::InitResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));
    RDouble **rhs = rhsProxy->GetField_UNS();

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTransition[m][iCell] = - rhs[m][iCell];
        }
    }
}

void TransitionSolverUnstr::RestrictAllQ(Grid *fineGridIn, Grid *coarseGridIn)
{
    UnstructGrid *fineGrid = UnstructGridCast(fineGridIn);
    UnstructGrid *coarseGrid = UnstructGridCast(coarseGridIn);

    RDouble **qTransitionFineGrid = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("q_transition"));
    RDouble **qTransitionCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("q_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        RestrictQ(fineGrid, qTransitionFineGrid[m], coarseGrid, qTransitionCoarseGrid[m]);
    }
}

void TransitionSolverUnstr::LoadFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **flux = faceProxy->GetFlux();

    //! Determine if there are boundary faces.
    int nMid = localStart;
    if (localEnd <= numberOfBoundaryFace)
    {
        //! If all boundary faces
        nMid = localEnd;
    }
    else if (localStart < numberOfBoundaryFace)
    {
        //! Part of them are boundary faces
        nMid = numberOfBoundaryFace;
    }

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));
    using namespace IDX;

    //! For boundary faces, remember re is ghost cell
    for (int i = localStart; i < nMid; ++ i)
    {
        int le, j;
        le = leftCellOfFace[i];
        j = i - localStart;
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            residualTransition[m][le] -= flux[m][j];
        }
    }

    //! Interior faces
    for (int i = nMid; i < localEnd; ++ i)
    {
        int le, re, j;
        le = leftCellOfFace[i];
        re = rightCellOfFace[i];
        j = i - localStart;

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            residualTransition[m][le] -= flux[m][j];
            residualTransition[m][re] += flux[m][j];
        }
    }
}

void TransitionSolverUnstr::SourceFluxTwoEquation(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **primitiveVariables = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));
    RDouble **spectrumTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_transition"));
    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));

    RDouble *wallDistance = grid->GetWallDist();
    RDouble *vol = grid->GetCellVolume();

    Param_TransitionSolver *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();

    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    RDouble density, ke, kw;
    RDouble sij2, mut;
    RDouble volume;

    RDouble s11, s22, s33, s12, s13, s23, w12, w13, w23;

    RDouble *dudx = new RDouble[numberOfTotal]();
    RDouble *dudy = new RDouble[numberOfTotal]();
    RDouble *dudz = new RDouble[numberOfTotal]();

    RDouble *dvdx = new RDouble[numberOfTotal]();
    RDouble *dvdy = new RDouble[numberOfTotal]();
    RDouble *dvdz = new RDouble[numberOfTotal]();

    RDouble *dwdx = new RDouble[numberOfTotal]();
    RDouble *dwdy = new RDouble[numberOfTotal]();
    RDouble *dwdz = new RDouble[numberOfTotal]();

    using namespace IDX;
    grid->CompGradient(primitiveVariables[IU], dudx, dudy, dudz);
    grid->CompGradient(primitiveVariables[IV], dvdx, dvdy, dvdz);
    grid->CompGradient(primitiveVariables[IW], dwdx, dwdy, dwdz);

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

    using namespace IDX;
    RDouble refMaNumber = parameters->GetRefMachNumber();

    RDouble compressibleFix = 1.0 + 0.40 * pow(refMaNumber, 1.75);

    RDouble ca1 = parameters->Getca1();
    RDouble ca2 = parameters->Getca2();
    RDouble ce1 = parameters->Getce1();
    RDouble ce2 = parameters->Getce2();
    RDouble cct = parameters->Getcct();
    RDouble s1 = parameters->Gets1();
    RDouble ccf = parameters->Getccf();

    RDouble *gamaeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("gamaeff"));

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        density = primitiveVariables[IR][iCell];
        ke = qTurbulence[IKE][iCell];
        kw = qTurbulence[IKW][iCell];
        mut = viscousTurbulence[iCell];
        volume = vol[iCell];

        s11 = dudx[iCell];
        s22 = dvdy[iCell];
        s33 = dwdz[iCell];
        s12 = half * (dudy[iCell] + dvdx[iCell]);
        s13 = half * (dudz[iCell] + dwdx[iCell]);
        s23 = half * (dvdz[iCell] + dwdy[iCell]);
        w12 = half * (dudy[iCell] - dvdx[iCell]);
        w13 = half * (dudz[iCell] - dwdx[iCell]);
        w23 = half * (dvdz[iCell] - dwdy[iCell]);

        sij2 = two * (s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23));

        RDouble mul = viscousLaminar[iCell];
        RDouble um = primitiveVariables[IU][iCell];
        RDouble vm = primitiveVariables[IV][iCell];
        RDouble wm = primitiveVariables[IW][iCell];
        RDouble intermittency = qTransition[IGAMA][iCell];
        RDouble Rectabar = qTransition[IRECT][iCell];

        RDouble wallDist = wallDistance[iCell];

        RDouble vorx = dvdz[iCell] - dwdy[iCell];
        RDouble vory = dudz[iCell] - dwdx[iCell];
        RDouble vorz = dudy[iCell] - dvdx[iCell];

        RDouble vorticity = DISTANCE(vorx, vory, vorz);
        RDouble strainRate = sqrt(sij2);

        RDouble absU = MAX(DISTANCE(um, vm, wm), SMALL);

        RDouble RT = ViscosityRatio(density, mul, ke, kw, refReNumber);
        RDouble Rev = ReynoldsNumberBasedOnStrainRate(density, wallDist, mul, strainRate, refReNumber);
        RDouble Rectac = TransitionOnsetMomentumThicknessReynolds(Rectabar);

        RDouble Rw = ReynoldsNumberBasedOnDissipation(density, wallDist, mul, kw, refReNumber);
        RDouble Flength = HighReynoldsCorrectionOfFlength(Rw, FlengthGivenByLangtry(Rectabar));

        RDouble Fonset = TransitionOnsetFunction(Rev, Rectac, RT);
        RDouble Fturb = ControlFunctionFturb(RT);

        RDouble production1OfGama = ca1 * Flength * density * strainRate * sqrt(intermittency * Fonset);
        RDouble production2OfGama = production1OfGama * (-ce1 * intermittency);

        RDouble productionOfGama = production1OfGama + production2OfGama;

        RDouble specGamap = MIN(0.0, 1.5 * production2OfGama + 0.5 * production1OfGama) / (density * intermittency);

        RDouble destruction2OfGama = -ca2 * Fturb * density * vorticity * intermittency;
        RDouble destruction1OfGama = -ce2 * intermittency * destruction2OfGama;

        RDouble destructionOfGama = destruction1OfGama + destruction2OfGama;
        RDouble specGamad = MAX(0.0, 2.0 * destruction1OfGama + destruction2OfGama) / (density * intermittency);

        residualTransition[IGAMA][iCell] += (productionOfGama - destructionOfGama) * volume;
        spectrumTransition[IGAMA][iCell] += (specGamad - specGamap) * volume;

        RDouble Fctat = 1.0, Fctatcf = 1.0;
        BlendingFunctionOfFctat(intermittency, Rw, Rectabar, vorticity, mul, density, absU, wallDist, refReNumber, ce2, Fctat, Fctatcf);
        RDouble tscl = TimeScaleInSourceTerm(density, absU, mul, mut, refReNumber);
        RDouble TU = TurbulenceIntensity(absU, ke);
        gamaeff[iCell] = SeparationCorrectionOfIntermittency(intermittency, Rev, Rectac, RT, Fctat, s1, refMaNumber);

        RDouble dUds = AccelerationAlongStreamline(um, vm, wm, dudx[iCell], dudy[iCell], dudz[iCell], dvdx[iCell], dvdy[iCell], dvdz[iCell], dwdx[iCell], dwdy[iCell], dwdz[iCell]);

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
        RDouble omegaStreamWise = abs((dwdy[iCell] - dvdz[iCell]) * um + (dudz[iCell] - dwdx[iCell]) * vm + (dvdx[iCell] - dudy[iCell]) * wm) / absU;
        RDouble heightCrossFlow = wallDist * omegaStreamWise / absU;
        RDouble Rescf = ReynoldsNumberBasedOnScf(heightCrossFlow, absU, density, mul, mut, refReNumber);

        RDouble coecommon = cct * (1.0 - Fctat) / MAX(tscl, SMALL);

        RDouble prodRe = coecommon * density * Rectat;
        RDouble destRe = coecommon * density * Rectabar;

        RDouble prodRecf = cct * density * ccf * MIN(Rescf - Rectabar, 0.0) * Fctatcf / tscl;
        RDouble specRe = destRe / (density * Rectabar);

        residualTransition[IRECT][iCell] += (prodRe + prodRecf - destRe) * volume;
        spectrumTransition[IRECT][iCell] += specRe * volume;
    }

    delete [] dudx;
    delete [] dudy;
    delete [] dudz;

    delete [] dvdx;
    delete [] dvdy;
    delete [] dvdz;

    delete [] dwdx;
    delete [] dwdy;
    delete [] dwdz;
}

FaceProxy *TransitionSolverUnstr::CreateFaceProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    FaceProxy *faceProxy = CreateFaceProxyTransition(grid);
    faceProxy->SetNext(CreateFaceProxyNS(grid));

    return faceProxy;
}

FaceProxy *TransitionSolverUnstr::CreateFaceProxyNS(Grid *gridIn)
{
    int nl;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);

    int nchem;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);

    int numberOfEquation = nl + nchem;

    FaceProxy *faceProxy = new FaceProxy();
    faceProxy->Create(SEGCTION_LENGTH, numberOfEquation);

    return faceProxy;
}

FaceProxy *TransitionSolverUnstr::CreateFaceProxyTransition(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    FaceProxy *faceProxy = new FaceProxy();
    faceProxy->Create(SEGCTION_LENGTH, nTransitionEquation);

    return faceProxy;
}

TransitionFaceValue *TransitionSolverUnstr::CreateTransitionFaceValue(Grid *gridIn)
{
    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    TransitionFaceValue *facevar = new TransitionFaceValue(nTransitionEquation, SEGCTION_LENGTH);
    return facevar;
}

void TransitionSolverUnstr::ComputeQTransitionNodeValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalNode = grid->GetNTotalNode();
    int numberOfTotalFace = grid->GetNTotalFace();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **qTransitionNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qTransitionNode"));
    RDouble **qTransition = reinterpret_cast<RDouble **>(grid->GetDataPtr("q_transition"));
    int *nodeValueSliceTransition = reinterpret_cast<int *>(grid->GetDataPtr("nodeValueSliceTransition"));

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *face2node = grid->GetFace2Node();
    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iNode = 0; iNode < numberOfTotalNode; ++ iNode)
    {
        nodeValueSliceTransition[iNode] = 0;

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            qTransitionNode[m][iNode] = 0;
        }
    }

    //! Boundary faces
    int nodePosition = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];

            if (bcType != PHENGLEI::INTERFACE)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    //! From left
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransitionNode[m][point] += qTransition[m][le];
                    }
                    nodeValueSliceTransition[point] += 1;

                    //! From right
                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransitionNode[m][point] += qTransition[m][re];
                    }
                    nodeValueSliceTransition[point] += 1;
                }
            }
            else
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];
                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransition[m][point] += qTransition[m][le];
                    }
                    nodeValueSliceTransition[point] += 1;
                }
            }
            nodePosition += nodeNumberOfEachFace[iFace];
        }
    }

    //! Interior faces
    for (int iFace = numberOfBoundaryFace; iFace < numberOfTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
        {
            //! From left
            int point = face2node[nodePosition + jNode];
            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                qTransitionNode[m][point] += qTransition[m][le];
            }
            ++ nodeValueSliceTransition[point];

            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                qTransitionNode[m][point] += qTransition[m][re];
            }
            ++ nodeValueSliceTransition[point];
        }
        nodePosition += nodeNumberOfEachFace[iFace];
    }

    nodePosition = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType == PHENGLEI::SYMMETRY)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransitionNode[m][point] = 0;
                    }
                    nodeValueSliceTransition[point] = 0;
                }
            }
            nodePosition += nodeNumberOfEachFace[iFace];
        }
    }

    nodePosition = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            int le = leftCellOfFace[iFace];
            int re = iFace + numberOfTotalCell;

            if (bcType == PHENGLEI::SYMMETRY)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransitionNode[m][point] += qTransition[m][le];
                    }
                    ++ nodeValueSliceTransition[point];

                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransitionNode[m][point] += qTransition[m][re];
                    }
                    ++ nodeValueSliceTransition[point];
                }
            }
            nodePosition += nodeNumberOfEachFace[iFace];
        }
    }

    nodePosition = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType == PHENGLEI::FARFIELD)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransitionNode[m][point] = 0;
                    }
                    nodeValueSliceTransition[point] = 0;
                }
            }
            nodePosition += nodeNumberOfEachFace[iFace];
        }
    }

    nodePosition = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            int le = leftCellOfFace[iFace];
            int re = iFace + numberOfTotalCell;

            if (bcType == PHENGLEI::FARFIELD)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransitionNode[m][point] += qTransition[m][le];
                    }
                    ++ nodeValueSliceTransition[point];

                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransitionNode[m][point] += qTransition[m][re];
                    }
                    ++ nodeValueSliceTransition[point];
                }
            }
            nodePosition += nodeNumberOfEachFace[iFace];
        }
    }

    nodePosition = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransitionNode[m][point] = 0;
                    }
                    nodeValueSliceTransition[point] = 0;
                }
            }
            nodePosition += nodeNumberOfEachFace[iFace];
        }
    }

    nodePosition = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            int le = leftCellOfFace[iFace];
            int re = iFace + numberOfTotalCell;
            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransitionNode[m][point] += qTransition[m][le];
                    }
                    ++ nodeValueSliceTransition[point];

                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        qTransitionNode[m][point] += qTransition[m][re];
                    }
                    ++ nodeValueSliceTransition[point];
                }
            }
            nodePosition += nodeNumberOfEachFace[iFace];
        }
    }

    //! Change the number of points divided by the interpoint to the number of points before partition,because it is a face loop,each point is calculated twice,so the number of points needs to be multiplied by 2.
    InterpointInformation *interPointInfor = grid->GetInterpointInfo();
    if (interPointInfor)
    {
        int numberOfInterPoint = interPointInfor->GetNumberOfInterpoints();
        int *labelOfInterPoint = interPointInfor->GetLabelOfInterPoint();
        int *interPoint2GlobalPoint = interPointInfor->GetInterPoint2GlobalPoint();

        int globalPoint = 0;
        for (int iPoint = 0; iPoint < numberOfInterPoint; ++ iPoint)
        {
            globalPoint = interPoint2GlobalPoint[iPoint];
            if (labelOfInterPoint[iPoint] != 1)
            {
                for (int m = 0; m < nTransitionEquation; ++ m)
                {
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    qTransitionNode[m][globalPoint] = 0;
                }
            }
        }
        RDouble **qTransitionInterPoint = reinterpret_cast<RDouble **> (grid->GetDataPtr("transition::qInterpoint"));
        for (int iPoint = 0; iPoint < numberOfInterPoint; ++ iPoint)
        {
            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                qTransitionInterPoint[m][iPoint] = 0;
            }
        }
    }
}

void TransitionSolverUnstr::ModifyQTransitionNodeValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalNode = grid->GetNTotalNode();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **qTransitionNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qTransitionNode"));
    int *nodeValueSliceTransition = reinterpret_cast<int *>(grid->GetDataPtr("nodeValueSliceTransition"));

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nodeValueSliceTransition[iNode] != 0)
        {
            for (int m = 0; m < nTransitionEquation; m++)
            {
                qTransitionNode[m][iNode] /= nodeValueSliceTransition[iNode];
            }
        }
        else
        {
            for (int m = 0; m < nTransitionEquation; m++)
            {
                qTransitionNode[m][iNode] = 0.0;
            }
        }
    }
}

void TransitionSolverUnstr::ComputeFaceWeight(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *deltL = faceProxy->GetWeightL();
    RDouble *deltR = faceProxy->GetWeightR();

    grid->FaceWeight(deltL, deltR, localStart, localEnd);
}

void TransitionSolverUnstr::ViscousFlux(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalFace = grid->GetNTotalFace();

    FaceProxy *faceProxy = CreateFaceProxy(grid);
    faceProxy->SetTransitionFaceValue(CreateTransitionFaceValue(grid));

    int localStart, localEnd;
    localStart = 0;
    do
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > numberOfTotalFace) localEnd = numberOfTotalFace;

        ComputeFaceWeight(grid, faceProxy, localStart, localEnd);

        GetVisFaceValue(grid, faceProxy, localStart, localEnd);

        ComputeVisflux(grid, faceProxy, localStart, localEnd);

        //LoadFlux(grid, faceProxy, localStart, localEnd);

        localStart = localEnd;
    } while (localStart < numberOfTotalFace);

    delete faceProxy;
}

void TransitionSolverUnstr::GetVisFaceValue(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    TransitionFaceValue *faceVariable = faceProxy->GetTransitionFaceValue();

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));

    Param_TransitionSolver *parameters = GetControlParameters();

    RDouble *mul = faceVariable->mul;
    RDouble *mut = faceVariable->mut;
    RDouble **mlt = faceVariable->mlt;

    int le, re, jFace;
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        le = leftCellOfFace[iFace];
        re = rightCellOfFace[iFace];
        if (iFace < numberOfBoundaryFace) re = iFace + numberOfTotalCell;
        jFace = iFace - localStart;

        mul[jFace] = half * (viscousLaminar[le] + viscousLaminar[re]);
        mut[jFace] = half * (viscousTurbulence[le] + viscousTurbulence[re]);
    }

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int nMid = numberOfBoundaryFace;
    if (localEnd < numberOfBoundaryFace) nMid = localEnd;

    if (localStart < numberOfBoundaryFace)
    {
        for (int iFace = localStart; iFace < nMid; ++ iFace)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
            int bcType = bcRegion->GetBCType();

            jFace = iFace - localStart;

            le = leftCellOfFace[iFace];

            if (IsInterface(bcType))
            {
                continue;
            }

            if (IsWall(bcType))
            {
                mul[jFace] = viscousLaminar[le];
                mut[jFace] = 0.0;
            }
            else
            {
                mul[jFace] = viscousLaminar[le];
                mut[jFace] = viscousTurbulence[le];
            }
        }
    }

    using namespace IDX;

    RDouble dct = parameters->Getdct();
    RDouble df = parameters->Getdf();

    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        le = leftCellOfFace[iFace];
        re = rightCellOfFace[iFace];
        if (iFace < numberOfBoundaryFace) re = iFace + numberOfTotalCell;
        jFace = iFace - localStart;

        mlt[IGAMA][jFace] = mul[jFace] + mut[jFace] / df;
        mlt[IRECT][jFace] = (mul[jFace] + mut[jFace]) * dct;
    }
}

//! Bell 20130513 mod.
void TransitionSolverUnstr::ComputeVisflux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int numberOfBoundaryFace = grid->GetNBoundFace();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    RDouble **flux = faceProxy->GetFlux();
    TransitionFaceValue *faceVariable = faceProxy->GetTransitionFaceValue();

    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();

    RDouble *cellCenterX = grid->GetCellCenterX();
    RDouble *cellCenterY = grid->GetCellCenterY();
    RDouble *cellCenterZ = grid->GetCellCenterZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble **mlt = faceVariable->mlt;

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));

    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();

    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble *dfd1 = new RDouble[nTransitionEquation];
    RDouble *dfd2 = new RDouble[nTransitionEquation];
    RDouble *dfdn = new RDouble[nTransitionEquation];

    RDouble *f1 = new RDouble[nTransitionEquation];
    RDouble *f2 = new RDouble[nTransitionEquation];
    RDouble *fMid = new RDouble[nTransitionEquation];

    int le, re, jFace;
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        le = leftCellOfFace[iFace];
        re = rightCellOfFace[iFace];
        jFace = iFace - localStart;

        RDouble nxs = xfn[iFace];
        RDouble nys = yfn[iFace];
        RDouble nzs = zfn[iFace];

        RDouble dxL = cellCenterX[le] - xfc[iFace];
        RDouble dyL = cellCenterY[le] - yfc[iFace];
        RDouble dzL = cellCenterZ[le] - zfc[iFace];

        RDouble dxR = cellCenterX[re] - xfc[iFace];
        RDouble dyR = cellCenterY[re] - yfc[iFace];
        RDouble dzR = cellCenterZ[re] - zfc[iFace];

        RDouble dL = nxs * dxL + nys * dyL + nzs * dzL;
        RDouble dR = nxs * dxR + nys * dyR + nzs * dzR;

        RDouble distTemp = - dL / (DISTANCE(dxL, dyL, dzL) + SMALL);
        distTemp = MIN(distTemp, 1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle1 = asin(distTemp) * 180.0 / PI;

        distTemp = dR / (DISTANCE(dxR, dyR, dzR) + SMALL);
        distTemp = MIN(distTemp, 1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle2 = asin(distTemp) * 180.0 / PI;

        //! quantities at points 1 and 2
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            f1[m] = qTransition[m][le];
            f2[m] = qTransition[m][re];
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            fMid[m] = half * (f1[m] + f2[m]);
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            dfdn[m] = 0.0;
        }

        if (angle1 > 0.0 && angle2 > 0.0 && ABS(dL) > TINY &&ABS(dR) > TINY)
        {
            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                dfd1[m] = (f1[m] - fMid[m]) / dL;
            }

            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                dfd2[m] = (f2[m] - fMid[m]) / dR;
            }

            RDouble dtmp = dL * dL + dR * dR;
            RDouble weightL = dL * dL / dtmp;
            RDouble weightR = 1.0 - weightL;

            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                dfdn[m] = weightL * dfd1[m] + weightR * dfd2[m];
            }

            if (iFace < numberOfBoundaryFace)
            {
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
                int bcType = bcRegion->GetBCType();

                if (bcType != PHENGLEI::INTERFACE && bcType != PHENGLEI::SYMMETRY)
                {
                    RDouble normalDist = (cellCenterX[re] - cellCenterX[le]) * nxs +
                        (cellCenterY[re] - cellCenterY[le]) * nys +
                        (cellCenterZ[re] - cellCenterZ[le]) * nzs;

                    for (int m = 0; m < nTransitionEquation; ++ m)
                    {
                        dfdn[m] = (f2[m] - f1[m]) / normalDist;
                    }
                }
            }
        }


        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            flux[m][jFace] = - oRefReNumber * mlt[m][jFace] * dfdn[m] * area[iFace];
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            residualTransition[m][le] -= flux[m][jFace];

            if (iFace >= numberOfBoundaryFace)
            {
                residualTransition[m][re] += flux[m][jFace];
            }
        }

    }

    delete [] dfd1;
    delete [] dfd2;
    delete [] dfdn;

    delete [] f1;
    delete [] f2;
    delete [] fMid;
}

void TransitionSolverUnstr::Diagonal(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalFace = grid->GetNTotalFace();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));

    RDouble **spectrumTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_transition"));
    RDouble **matrixTransitionLeft = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_transitionl"));
    RDouble **matrixTransitionRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_transitionr"));

    RDouble *vol = grid->GetCellVolume();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vgn = grid->GetFaceNormalVelocity();

    RDouble uL, vL, wL, vnL, absVnL;
    RDouble uR, vR, wR, vnR, absVnR;

    RDouble vn, absVn, ns2, rho, mul, mut, volume, oVolume;

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    vector<RDouble> work(nTransitionEquation);

    using namespace IDX;

    for (int iFace = 0; iFace < numberOfTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        //! Calculate velocity at the cell interface.
        uL = q[IU][le];
        vL = q[IV][le];
        wL = q[IW][le];
        vnL = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];

        absVnL = ABS(vnL);

        uR = q[IU][re];
        vR = q[IV][re];
        wR = q[IW][re];
        vnR = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR - vgn[iFace];

        absVnR = ABS(vnR);

        vn = half * (vnL + vnR);
        absVn = ABS(vn);

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            spectrumTransition[m][le] += area[iFace] * half * absVn;
            spectrumTransition[m][re] += area[iFace] * half * absVn;

            matrixTransitionRight[m][iFace] = area[iFace] * half * (vnR - absVn);
            matrixTransitionLeft[m][iFace] = area[iFace] * half * (- vnL - absVn);
        }
    }

    RDouble dct = parameters->Getdct();
    RDouble df = parameters->Getdf();

    for (int iFace = 0; iFace < numberOfTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        rho = half * (q[IR][le] + q[IR][re]);
        mul = half * (viscousLaminar[le] + viscousLaminar[re]);
        mut = half * (viscousTurbulence[le] + viscousTurbulence[re]);


        volume = half * (vol[le] + vol[re]);
        oVolume = one / volume;

        ns2 = area[iFace] * area[iFace];

        work[IGAMA] = (mul + mut / dct) / rho * ns2 * oVolume;
        work[IRECT] = (mul + mut) * df / rho * ns2 * oVolume;

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            spectrumTransition[m][le] += oRefReNumber * work[m];
            spectrumTransition[m][re] += oRefReNumber * work[m];

            matrixTransitionLeft[m][iFace] += - oRefReNumber * work[m];
            matrixTransitionRight[m][iFace] += - oRefReNumber * work[m];
        }
    }
}

void TransitionSolverUnstr::GetQlQr(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    GetQlQrTransition(grid, faceProxy, localStart, localEnd);

    FaceProxy *faceProxyNS = faceProxy->GetNext();
    GetQlQrNS(grid, faceProxyNS, localStart, localEnd);

}

void TransitionSolverUnstr::GetQlQrTransition(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **qTransitionVarables = reinterpret_cast <RDouble **> (grid->GetDataPtr("q_transition"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le, re, jFace;
        le = leftCellOfFace[iFace];
        re = rightCellOfFace[iFace];
        jFace = iFace - localStart;

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            qL[m][jFace] = qTransitionVarables[m][le];
            qR[m][jFace] = qTransitionVarables[m][re];
        }
    }

    //! Boundary fix.
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int nBoundFace = grid->GetNBoundFace();
    int nMid = 0;

    //! Check if there are boundary faces. If no, return.
    if (localStart >= nBoundFace)
    {
        return;
    }

    if (localEnd <= nBoundFace)
    {
        //! If they are all boundary faces.
        nMid = localEnd;
    }
    else
    {
        //! Part of them are boundary faces.
        nMid = nBoundFace;
    }

    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        int jFace = iFace - localStart;

        if (bcType == PHENGLEI::INTERFACE)
        {
            continue;
        }

        if (bcType == PHENGLEI::SYMMETRY)
        {
            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                qR[m][jFace] = qL[m][jFace];
            }
        }
        else
        {
            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                RDouble temp = half * (qR[m][jFace] + qL[m][jFace]);
                qR[m][jFace] = temp;
                qL[m][jFace] = temp;
            }
        }
    }
}

void TransitionSolverUnstr::GetQlQrNS(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    int nl = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation = nl + nchem;

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le, re, jFace;
        le = leftCellOfFace[iFace];
        re = rightCellOfFace[iFace];
        jFace = iFace - localStart;

        for (int m = 0; m < nEquation; ++ m)
        {
            qL[m][jFace] = q[m][le];
            qR[m][jFace] = q[m][re];
        }
    }

    //! Boundary fix.
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int nBoundFace = grid->GetNBoundFace();
    int nMid = 0;

    //! Check if there are boundary faces. If no, return.
    if (localStart >= nBoundFace)
    {
        return;
    }

    if (localEnd <= nBoundFace)
    {
        //! If they are all boundary faces.
        nMid = localEnd;
    }
    else
    {
        //! Part of them are boundary faces.
        nMid = nBoundFace;
    }

    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        int jFace = iFace - localStart;

        if (bcType == PHENGLEI::INTERFACE)
        {
            continue;
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            RDouble temp = half * (qR[m][jFace] + qL[m][jFace]);
            qR[m][jFace] = temp;
            qL[m][jFace] = temp;
        }
    }
}

void TransitionSolverUnstr::InviscidFlux(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalFace = grid->GetNTotalFace();

    FaceProxy *faceProxy = CreateFaceProxy(grid);

    int localStart, localEnd;

    localStart = 0;
    do
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > numberOfTotalFace) localEnd = numberOfTotalFace;

        GetQlQr(grid, faceProxy, localStart, localEnd);

        ComputeInviscidFlux(grid, faceProxy, localStart, localEnd);

        LoadFlux(grid, faceProxy, localStart, localEnd);

        localStart = localEnd;

    } while (localStart < numberOfTotalFace);

    delete faceProxy;
}

void TransitionSolverUnstr::ComputeInviscidFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    NNDFlux(grid, faceProxy, localStart, localEnd);
}

void TransitionSolverUnstr::NNDFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *vgn = grid->GetFaceNormalVelocity();
    RDouble *area = grid->GetFaceArea();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();
    RDouble **flux = faceProxy->GetFlux();

    FaceProxy *faceProxyNS = faceProxy->GetNext();

    RDouble **qPrimitiveVariableL = faceProxyNS->GetQL();
    RDouble **qPrimitiveVariableR = faceProxyNS->GetQR();

    RDouble *primitiveVariableL = new RDouble[nTransitionEquation];
    RDouble *primitiveVariableR = new RDouble[nTransitionEquation];
    RDouble *f = new RDouble[nTransitionEquation];

    RDouble rL, uL, vL, wL, vnL, vnLL, rR, uR, vR, wR, vnR, vnRR;
    RDouble rhoFlag = 0, rhoCoefficientL, rhoCoefficientR;

    using namespace IDX;

    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le, re, jFace;
        jFace = iFace - localStart;
        le = leftCellOfFace[iFace];
        re = rightCellOfFace[iFace];

        rL = qPrimitiveVariableL[IR][jFace];
        uL = qPrimitiveVariableL[IU][jFace];
        vL = qPrimitiveVariableL[IV][jFace];
        wL = qPrimitiveVariableL[IW][jFace];

        rR = qPrimitiveVariableR[IR][jFace];
        uR = qPrimitiveVariableR[IU][jFace];
        vR = qPrimitiveVariableR[IV][jFace];
        wR = qPrimitiveVariableR[IW][jFace];

        vnL = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];
        vnR = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR - vgn[iFace];

        vnLL = half * (vnL + ABS(vnL));
        vnRR = half * (vnR - ABS(vnR));

        rhoCoefficientL = rhoFlag * one + (one - rhoFlag) * rL;
        rhoCoefficientR = rhoFlag * one + (one - rhoFlag) * rR;

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            primitiveVariableL[m] = rhoCoefficientL * qL[m][jFace];
            primitiveVariableR[m] = rhoCoefficientR * qR[m][jFace];
        }

        //! Consider vgn here!
        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            f[m] = (vnLL * primitiveVariableL[m] + vnRR * primitiveVariableR[m]);
        }

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            flux[m][jFace] = f[m] * area[iFace];
        }
    }

    delete [] primitiveVariableL;
    delete [] primitiveVariableR;
    delete [] f;
}

void TransitionSolverUnstr::GetResidual(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    RDouble **residualTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_transition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    PHSPACE::ComputeResidualonGrid(grid, actkey, residualTransition, nTransitionEquation);
}

void TransitionSolverUnstr::CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation) return;

    int iNeighborZone = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->MoveToBegin();

    RDouble **dq = fieldProxy->GetField_UNS();

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
    {
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int sourceCell;
            int iFace = interfaceIndexContainerForSend[iLocalFace];
            grid->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
            for (int m = 0; m < nEquation; ++ m)
            {
                PHWrite(dataContainer, dq[m][sourceCell]);
            }
        }
    }
}

void TransitionSolverUnstr::DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation) return;

    int iNeighborZone = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();
    RDouble **dq = fieldProxy->GetField_UNS();

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
    {
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForReceive[iLocalFace];
            int targetCell;
            grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);

            for (int m = 0; m < nEquation; ++ m)
            {
                PHRead(dataContainer, dq[m][targetCell]);
            }
        }
    }
}

void TransitionSolverUnstr::Boundary(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    using namespace PHENGLEI;
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();
   
    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");

    //! Parallel
    int currentZone = grid->GetOrdinaryGridIndex();
    //! Serial
    if (currentZone == -1)
    {
        currentZone = grid->GetZoneID();
    }
    int MixingPlaneNumber = 0;

    string MixingIn, MixingOut;
    if (referenceFrame == ROTATIONAL_FRAME)
    {
        int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
        string MixingPlane[100];
        GlobalDataBase::GetData("MixingPlane", &MixingPlane, PHSTRING, 2 * nTurboZone - 2);

        if (currentZone == 0)
        {
            MixingOut = MixingPlane[2 * currentZone];
        }
        else if (currentZone == nTurboZone)
        {
            MixingIn = MixingPlane[2 * currentZone - 1];
        }
        else
        {
            MixingIn = MixingPlane[2 * currentZone - 1];
            MixingOut = MixingPlane[2 * currentZone];
        }
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        if (bcName == MixingIn || bcName == MixingOut)
        {
            continue;
        }
        if (IsInterface(bcType))
        {
            continue;
        }
        if (bcType == EXTRAPOLATION || bcType == SYMMETRY || bcType == OUTFLOW || bcType == PRESSURE_OUTLET || bcType == MASS_FLOW_OUTLET)
        {
            OutflowBCRegion(grid, bcRegion);
        }
        else if (bcType == INFLOW || bcType == PRESSURE_INLET || bcType == MASS_FLOW_INLET)
        {
            InflowBCRegion(grid, bcRegion);
        }
        else if (IsWall(bcType))
        {
            WallBCRegion(grid, bcRegion);
        }
        else if (bcType == FARFIELD)
        {
            FarfieldBCRegion(grid, bcRegion);
        }
        else
        {
            OutflowBCRegion(grid, bcRegion);
        }
    }
}

void TransitionSolverUnstr::OutflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        // iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            qTransition[m][re] = qTransition[m][le];
        }
    }
}

void TransitionSolverUnstr::InflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble *freeStreamTransitionVar = parameters->GetFreeStreamTransitionVar();
    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        // iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int re = rightCellOfFace[iFace];

        for (int m = 0; m < nTransitionEquation; ++ m)
        {
            qTransition[m][re] = freeStreamTransitionVar[m];
        }
    }
}

void TransitionSolverUnstr::WallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));

    using namespace IDX;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        // iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        qTransition[IGAMA][re] = qTransition[IGAMA][le];
        qTransition[IRECT][re] = qTransition[IRECT][le];
    }
}

void TransitionSolverUnstr::FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *nxs = grid->GetFaceNormalX();
    RDouble *nys = grid->GetFaceNormalY();
    RDouble *nzs = grid->GetFaceNormalZ();

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble refGama = parameters->GetRefGama();
    RDouble *freeStreamTransitionVar = parameters->GetFreeStreamTransitionVar();

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));
    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *gamma = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));

    using namespace IDX;
    RDouble *primitiveVariableInflow = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble roo = primitiveVariableInflow[IR];
    RDouble uoo = primitiveVariableInflow[IU];
    RDouble voo = primitiveVariableInflow[IV];
    RDouble woo = primitiveVariableInflow[IW];
    RDouble poo = primitiveVariableInflow[IP];

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        // iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        RDouble rin, uin, vin, win, pin;
        RDouble vno, vni, vei, coo, cIN;
        RDouble riemp, riemm, cb, vnb;

        rin = q[IR][le];
        uin = q[IU][le];
        vin = q[IV][le];
        win = q[IW][le];
        pin = q[IP][le];

        vno = nxs[iFace] * uoo + nys[iFace] * voo + nzs[iFace] * woo;
        vni = nxs[iFace] * uin + nys[iFace] * vin + nzs[iFace] * win;
        vei = sqrt(uin * uin + vin * vin + win * win);

        RDouble gama = gamma[le];

        coo = sqrt(ABS(refGama * poo / roo));
        cIN = sqrt(ABS(gama * pin / rin));

        //! supersonic
        if (vei > cIN)
        {
            if (vni >= 0.0)
            {
                //! exit
                for (int m = 0; m < nTransitionEquation; ++ m)
                {
                    qTransition[m][re] = qTransition[m][le];
                }
            }
            else
            {
                //! inlet
                for (int m = 0; m < nTransitionEquation; ++ m)
                {
                    qTransition[m][re] = freeStreamTransitionVar[m];
                }
            }
            continue;
        }

        //! subsonic
        riemp = vni + 2.0 * cIN / (gama - 1.0);
        riemm = vno - 2.0 * coo / (refGama - 1.0);
        vnb = half * (riemp + riemm);
        cb = fourth * (riemp - riemm) * gama;

        if (vnb >= 0.0)
        {
            //! exit
            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                qTransition[m][re] = qTransition[m][le];
            }
        }
        else
        {
            //! inlet
            for (int m = 0; m < nTransitionEquation; ++ m)
            {
                qTransition[m][re] = freeStreamTransitionVar[m];
            }
        }

    }
}

void TransitionSolverUnstr::InitCGrid(Grid *fineGridIn, Grid *coarseGridIn)
{
    UnstructGrid *fineGrid = UnstructGridCast(fineGridIn);
    UnstructGrid *coarseGrid = UnstructGridCast(coarseGridIn);

    int *cell2CoarseGridCell = fineGrid->GetCell2CoarseGridCell();
    int numberOfTotalCellOnFineGrid = fineGrid->GetNTotalCell();
    int numberOfTotalCellOnCoarseGrid = coarseGrid->GetNTotalCell();

    RDouble *volumeOfCoarseGrid = coarseGrid->GetCellVolume();
    RDouble *volumeOfFineGrid = fineGrid->GetCellVolume();

    RDouble *walldistanceOfFineGrid = fineGrid->GetWallDist();

    if (!coarseGrid->GetWallDist())
    {
        coarseGrid->AllocateWalldist();
    }
    RDouble *walldistanceOfCoarseGrid = coarseGrid->GetWallDist();

    for (int iCell = 0; iCell < numberOfTotalCellOnCoarseGrid; ++ iCell)
    {
        walldistanceOfCoarseGrid[iCell] = 0.0;
    }

    for (int iCell = 0; iCell < numberOfTotalCellOnFineGrid; ++ iCell)
    {
        walldistanceOfCoarseGrid[cell2CoarseGridCell[iCell]] += volumeOfFineGrid[iCell] * walldistanceOfFineGrid[iCell];
    }

    for (int iCell = 0; iCell < numberOfTotalCellOnCoarseGrid; ++ iCell)
    {
        walldistanceOfCoarseGrid[iCell] /= volumeOfCoarseGrid[iCell];
    }
}

void TransitionSolverUnstr::UploadInterfaceData(ActionKey *actkey)
{
    UploadInterfaceValue(actkey);
}

void TransitionSolverUnstr::UploadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (interfaceInformation == 0) return;

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));
    PHSPACE::UploadInterfaceValue(grid, qTransition, "transition::q", nTransitionEquation);
}

void TransitionSolverUnstr::UploadInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterpointInformation *interPointInformation = grid->GetInterpointInfo();
    if (interPointInformation == 0)
    {
        return;
    }

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble **qTransitionNode = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTransitionNode"));
    PHSPACE::UploadInterpointValue(grid, qTransitionNode, "qTransitionNode", nTransitionEquation);
}

void TransitionSolverUnstr::DownloadInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterpointInformation *interPointInformation = grid->GetInterpointInfo();
    if (interPointInformation == 0)
    {
        return;
    }

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble **qTransitionNode = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTransitionNode"));
    PHSPACE::DownloadInterpointValue(grid, qTransitionNode, "qTransitionNode", nTransitionEquation);
}

void TransitionSolverUnstr::CommunicationInterpointWeight(ActionKey *actkey)
{
    int level = actkey->level;

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int myid = GetCurrentProcessorID();
        int tag = iZone + level * nZones;
        int send_proc = GetZoneProcessorID(iZone);

        ZoneNeighbor *zoneNeighborForPoint = zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);
        std::size_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();

        for (std::size_t iNeighborZone = 0; iNeighborZone < numberOfNeighbor; ++ iNeighborZone)
        {
            int neighborZoneIndex = zoneNeighborForPoint->GetZoneIndexOfNeighbor(iNeighborZone);
            int recv_proc = GetZoneProcessorID(neighborZoneIndex);

            DataContainer *cdata = new DataContainer();

            if (myid == send_proc)
            {
                UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(iZone, level));
                InterpointInformation *interpointInformation = grid->GetInterpointInfo();
                int *nodeValueSliceTransition = reinterpret_cast<int *>(grid->GetDataPtr("nodeValueSliceTransition"));

                cdata->MoveToBegin();

                int iNeighbor = interpointInformation->FindIthNeighbor(neighborZoneIndex);
                int numberOfInterpointsForNeighbor = interpointInformation->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForSend = interpointInformation->GetPointIndexForSend(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                int valueWeight;
                int *qtmp = new int[numberOfInterpointsForNeighbor];
                for (int ipointLocal = 0; ipointLocal < numberOfInterpointsForNeighbor; ++ ipointLocal)
                {
                    iPoint = pointIndexForSend[ipointLocal];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    valueWeight = nodeValueSliceTransition[globalPoint];

                    qtmp[ipointLocal] = valueWeight;
                }
                cdata->Write(qtmp, numberOfInterpointsForNeighbor * sizeof(int));
                delete [] qtmp;     qtmp = nullptr;
            }

            PH_Trade(cdata, send_proc, recv_proc, tag);

            if (myid == recv_proc)
            {
                UnstructGrid *gridNeighbor = UnstructGridCast(PHSPACE::GetGrid(neighborZoneIndex, level));
                InterpointInformation *interpointInformationNeighbor = gridNeighbor->GetInterpointInfo();
                int *nodeValueSliceTransitionTrade = reinterpret_cast<int *>(gridNeighbor->GetDataPtr("nodeValueSliceTransitionTrade"));

                cdata->MoveToBegin();

                int iNeighbor = interpointInformationNeighbor->FindIthNeighbor(iZone);
                int numberOfInterpointsForNeighbor = interpointInformationNeighbor->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForReceive = interpointInformationNeighbor->GetPointIndexForReceive(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformationNeighbor->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                int valueWeight;
                for (int j = 0; j < numberOfInterpointsForNeighbor; ++ j)
                {
                    iPoint = pointIndexForReceive[j];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    cdata->Read(&valueWeight, sizeof(int));

                    nodeValueSliceTransitionTrade[globalPoint] += valueWeight;
                }
            }

            delete cdata;    cdata = nullptr;
        }
    }
}

void TransitionSolverUnstr::DownloadInterpointWeight(ActionKey *actkey)
{
    int level = actkey->level;

    using namespace PHMPI;

    int myid = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int currentZoneProc = GetZoneProcessorID(iZone);

        if (myid != currentZoneProc)
        {
            continue;
        }

        UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(iZone, level));
        int nTotalNode = grid->GetNTotalNode();
        int *nodeValueSliceTransition = reinterpret_cast<int *>(grid->GetDataPtr("nodeValueSliceTransition"));
        int *nodeValueSliceTransitionTrade = reinterpret_cast<int *>(grid->GetDataPtr("nodeValueSliceTransitionTrade"));

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            nodeValueSliceTransition[iNode] += nodeValueSliceTransitionTrade[iNode];
        }

        PHSPACE::SetField(nodeValueSliceTransitionTrade, 0, nTotalNode);
    }
}

// Bell 20120910 add
void TransitionSolverUnstr::DownloadInterfaceData(ActionKey *actkey)
{
    DownloadInterfaceValue(actkey);
}
//

void TransitionSolverUnstr::DownloadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (interfaceInformation == 0) return;

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();
    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));
    PHSPACE::DownloadInterfaceValue(grid, qTransition, "transition::q", nTransitionEquation);
}

//! need further compllished!!!
void TransitionSolverUnstr::UploadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));
    PHSPACE::UploadOversetValue(grid, qTransition, "transition::q", nTransitionEquation);
}

//! need further compllished!!!
void TransitionSolverUnstr::DownloadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    RDouble **qTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_transition"));

    PHSPACE::DownloadOversetValue(grid, qTransition, "transition::q", nTransitionEquation);
}

void TransitionSolverUnstr::ReadStatisticalFlow(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **qAverageTransition = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTransition"));

    Param_TransitionSolver *parameters = GetControlParameters();
    int nTransitionEquation = parameters->GetNTransitionEquation();

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    int numberOfStatisticalStep = 0;
    cdata->Read(&numberOfStatisticalStep, sizeof(int));

    GlobalDataBase::UpdateData("nStatisticalStep", &numberOfStatisticalStep, PHINT, 1);

    for (int m = 0; m < nTransitionEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            cdata->Read(&qAverageTransition[m][iCell], sizeof(RDouble));
        }
    }
}

void TransitionSolverUnstr::ComputeGradient(Grid *gridIn)
{
    GetGradientField(gridIn);
}

void TransitionSolverUnstr::InitMixingPlane(RDouble ***MixingPlaneVar, int Dim1, int Dim2, int Dim3, RDouble Value)
{
    for (int iDim1 = 0; iDim1 < Dim1; iDim1++)
    {
        for (int iDim2 = 0; iDim2 < Dim2; iDim2++)
        {
            for (int iDim3 = 0; iDim3 < Dim3; iDim3++)
            {
                MixingPlaneVar[iDim1][iDim2][iDim3] = Value;
            }
        }
    }
}

//! need to input q_average, which means flux averaged primitive variable of mixing plane.
//! this function should be called in boundary condition part.
void TransitionSolverUnstr::AverageMixingPlane(Grid *grid)
{
    UnstructGrid *gridUnstruct   = UnstructGridCast(grid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_TransitionSolver *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1   = refGama - 1;

    RDouble refDensity  = parameters->GetRefDimensionalDensity();
    RDouble refPressure = parameters->GetRefDimensionalPressure();

    RDouble **q           = reinterpret_cast<RDouble **> (gridUnstruct->GetDataPtr("q"));
    RDouble **qTransition = reinterpret_cast<RDouble **> (gridUnstruct->GetDataPtr("q_transition"));

    RDouble **r_ori    = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Face"));
    RDouble **r_target = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Span"));

    RDouble **SpanArea = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("SpanArea"));

    //! nSpan is the number of spanwise point to do data average
    //! also interpolate data on these points.
    //! nSpanSection
    int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");
    RDouble nSpan = (RDouble)(nSpanSection);

    RDouble Area;

    RDouble *nxs, *nys, *nzs, *ns;

    nxs = gridUnstruct->GetFaceNormalX();
    nys = gridUnstruct->GetFaceNormalY();
    nzs = gridUnstruct->GetFaceNormalZ();
    ns  = gridUnstruct->GetFaceArea();

    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    string MixingPlane[100];
    GlobalDataBase::GetData("MixingPlane", &MixingPlane, PHSTRING, 2 * nTurboZone - 2);

    //r_ori should be an array that contains coordinate of each node.
    //the range of this array is equal to the face number of this boundary.
    RDouble r_min, r_max;

    r_min = 1e6;
    r_max = 0;

    RDouble *xfc = gridUnstruct->GetFaceCenterX();
    RDouble *yfc = gridUnstruct->GetFaceCenterY();
    RDouble *zfc = gridUnstruct->GetFaceCenterZ();

    RDouble *gama = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama"));

    int **face2nodeArray = gridUnstruct->GetFace2NodeArray();
    int * node_number_of_each_face = gridUnstruct->GetNodeNumberOfEachFace();

    int *leftCellofFace  = gridUnstruct->GetLeftCellOfFace();
    int *rightCellofFace = gridUnstruct->GetRightCellOfFace();
    int nTotalCell       = gridUnstruct->GetNTotalCell();

    int nBCRegion = unstructBCSet->GetnBCRegion();

    RDouble PeriodicRotationAngle[100];
    GlobalDataBase::GetData("PeriodicRotationAngle", &PeriodicRotationAngle, PHDOUBLE, nTurboZone);

    using namespace IDX;

    RDouble ***SpanFlux      = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanFlux"));
    RDouble ***SpanTransFlux = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanTransFlux"));
    int FlowDirection = -1;

    int nEquation = GetNumberOfEquations();

    RDouble *vgn = gridUnstruct->GetFaceNormalVelocity();

    RDouble *primi = new RDouble[nEquation]();
    RDouble *primo = new RDouble[nEquation]();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        //! find mixing plane and average data.
        //! Parallel
        int iTurboZone = gridUnstruct->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = gridUnstruct->GetZoneID();
        }

        RDouble rotationAngle = PeriodicRotationAngle[iTurboZone] / 180 * PI;

        string MixingPlaneIn, MixingPlaneOut;
        int SourceFlag;

        if (iTurboZone == 0)
        {
            MixingPlaneOut = MixingPlane[0];
        }
        else if (iTurboZone == nTurboZone)
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
        }
        else
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
            MixingPlaneOut = MixingPlane[2 * iTurboZone];
        }

        if (bcName == MixingPlaneIn || bcName == MixingPlaneOut)
        {
            if (bcName == MixingPlaneIn)
            {
                FlowDirection = 0;
                SourceFlag = 1 + (iTurboZone - 1) * 10;
            }
            if (bcName == MixingPlaneOut)
            {
                FlowDirection = 1;
                SourceFlag = 0 + (iTurboZone + 1) * 10;
            }
            int MixingFlag = FlowDirection + iTurboZone * 10;

            //! preprocess of mixing plane.
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;

                int le = leftCellofFace[iFace];
                int re = rightCellofFace[iFace];

                //! need to compute total area, total flux of this bc region.
                //also need to compute flux of each spanwise.
                //spanwise can be defined based on mesh or manually.
                r_ori[MixingFlag][iFace] = sqrt(yfc[iFace] * yfc[iFace] + zfc[iFace] * zfc[iFace]);

                r_min = MIN(r_ori[MixingFlag][iFace], r_min);
                r_max = MAX(r_ori[MixingFlag][iFace], r_max);
            }

            //! get target radius.
            for (int iSpan = 0; iSpan <= nSpanSection; iSpan++)
            {
                r_target[MixingFlag][iSpan] = r_min + iSpan * (r_max - r_min) / nSpan;
            }

            //! compute average flux.
            for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
            {
                RDouble Flux[5]      = { 0, 0, 0, 0, 0 };
                RDouble TransFlux[5] = { 0, 0, 0, 0, 0 };

                RDouble TotalArea = 0.0;

                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    //! iFace is the face number in the set of faceIndex.
                    int iFace = *iter;

                    int le = leftCellofFace[iFace];
                    int re = rightCellofFace[iFace];

                    Area = ns[iFace];

                    RDouble r_Face = sqrt(yfc[iFace] * yfc[iFace] + zfc[iFace] * zfc[iFace]);

                    RDouble Velocity_ti = (q[IV][le] * zfc[iFace] - q[IW][le] * yfc[iFace]) / r_Face;
                    RDouble Velocity_ri = (q[IV][le] * yfc[iFace] + q[IW][le] * zfc[iFace]) / r_Face;

                    RDouble Velocity_to = (q[IV][re] * zfc[iFace] - q[IW][re] * yfc[iFace]) / r_Face;
                    RDouble Velocity_ro = (q[IV][re] * yfc[iFace] + q[IW][re] * zfc[iFace]) / r_Face;

                    RDouble nrs = (nys[iFace] * yfc[iFace] + nzs[iFace] * zfc[iFace]) / r_Face;

                    RDouble vni = q[IU][le] * nxs[iFace] + Velocity_ri * nrs - vgn[iFace];
                    RDouble vno = q[IU][re] * nxs[iFace] + Velocity_ro * nrs - vgn[iFace];
                    //! find face center that locate in target span.
                    //! compute flux at each span.
                    if (r_Face > r_target[MixingFlag][iSpan] && r_Face <= r_target[MixingFlag][iSpan + 1])
                    {
                        TotalArea += Area;
                        for (int m = 0; m < nEquation; m++)
                        {
                            TransFlux[m] += q[IR][le] * vni * qTransition[m][le] * ns[iFace];
                        }
                    }
                }

                SpanArea[MixingFlag][iSpan] = TotalArea;

                for (int m = 0; m < nEquation; m++)
                {
                    SpanTransFlux[MixingFlag][m][iSpan] = TransFlux[m] / SpanArea[MixingFlag][iSpan];
                }
            }
        }
    }
}

//! SpanFlux from AverageMixingPlane should be input data.
//! Here we will obtain averaged primitive variable of each span from SpanFlux.
//! The averaged primitive variable will be treated as outlet data of upstream and inlet data of downsteam.
void TransitionSolverUnstr::MixingPlaneDataTransfer(Grid *grid, Grid *NeighborGrid)
{
    UnstructGrid *gridUnstruct   = UnstructGridCast(grid);
    UnstructGrid *gridSource     = UnstructGridCast(NeighborGrid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_TransitionSolver *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1 = refGama - 1;

    RDouble refDensity  = parameters->GetRefDimensionalDensity();
    RDouble refPressure = parameters->GetRefDimensionalPressure();

    RDouble **q           = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q"));
    RDouble **qTransition = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q_transition"));

    RDouble **r_target = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Span"));

    RDouble **SpanArea = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("SpanArea"));
    //! nSpan is the number of spanwise point to do data average
    //! also interpolate data on these points.
    int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");

    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    string MixingPlane[100];
    GlobalDataBase::GetData("MixingPlane", &MixingPlane, PHSTRING, 2 * nTurboZone - 2);

    RDouble *xfc = gridUnstruct->GetFaceCenterX();
    RDouble *yfc = gridUnstruct->GetFaceCenterY();
    RDouble *zfc = gridUnstruct->GetFaceCenterZ();

    int nBCRegion = unstructBCSet->GetnBCRegion();

    //! qSpanNeighbor: record flow variable in neighbor zone at each span.
    RDouble ***qTrans_Span   = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTrans_Span"));
    RDouble ***SourceFlux    = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanFlux"));
    RDouble ***SpanTransFlux = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanTransFlux"));
    RDouble ***qSpanNeighbor = reinterpret_cast< RDouble *** > (gridSource->GetDataPtr("qTransSpanNeighbor"));

    using namespace IDX;

    int FlowDirection = -1;

    int nEquation = GetNumberOfEquations();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        //! Parallel
        int iTurboZone = gridUnstruct->GetOrdinaryGridIndex();
        int SourceZone = gridSource->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = gridUnstruct->GetZoneID();
            SourceZone = gridSource->GetZoneID();
        }

        string MixingPlaneIn, MixingPlaneOut;

        if (iTurboZone == 0)
        {
            MixingPlaneOut = MixingPlane[0];
        }
        else if (iTurboZone == nTurboZone)
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
        }
        else
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
            MixingPlaneOut = MixingPlane[2 * iTurboZone];
        }

        //! TargetMixingFlag is Mixingplane in current zone, which need neighbor zone data.
        //! SourceMixingFlag is Mixingplane in neighbor zone, which provide the data. 
        int MixingFlag, TargetFlag;
        if (bcName == MixingPlaneIn || bcName == MixingPlaneOut)
        {
            //! downstream to upstream: transfer data from MixingIn to MixingOut.
            if (SourceZone == iTurboZone + 1)
            {
                if (bcName == MixingPlaneIn) continue;
                //! upstream MixingIn.
                TargetFlag = 0 + SourceZone * 10;
                //! current MixingOut.
                MixingFlag = 1 + iTurboZone * 10;
            }
            //! upstream to downstream: transfer data from MixingOut to MixingIn.
            else if (SourceZone == iTurboZone - 1)
            {
                if (bcName == MixingPlaneOut) continue;
                TargetFlag = 1 + SourceZone * 10;
                MixingFlag = 0 + iTurboZone * 10;
            }

            //! compute flow variables.
            for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
            {
                for (int m = 0; m < nEquation; m++)
                {
                    qTrans_Span[MixingFlag][m][iSpan] = SpanTransFlux[MixingFlag][m][iSpan] / SourceFlux[MixingFlag][IR][iSpan];
                }

                for (int m = 0; m < nEquation; m++)
                {
                    qSpanNeighbor[TargetFlag][m][iSpan] = qTrans_Span[MixingFlag][m][iSpan];
                }
            }
        }
    }
}

void TransitionSolverUnstr::NonReflective(Grid *grid, Grid *NeighborGrid)
{
    UnstructGrid *gridUnstruct = UnstructGridCast(grid);
    UnstructGrid *gridTarget = UnstructGridCast(NeighborGrid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_TransitionSolver *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1 = refGama - 1;

    RDouble refDensity = parameters->GetRefDimensionalDensity();
    RDouble refPressure = parameters->GetRefDimensionalPressure();

    RDouble **q = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q"));
    RDouble **qTransition = reinterpret_cast< RDouble **> (gridUnstruct->GetDataPtr("q_transition"));

    RDouble **r_ori = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Face"));
    RDouble **r_target = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Span"));

    RDouble ***qTrans_Span = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTrans_Span"));

    RDouble **SpanArea = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("SpanArea"));
    //! nSpan is the number of spanwise point to do data average
    //! also interpolate data on these points.
    int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");

    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    string MixingPlane[100];
    GlobalDataBase::GetData("MixingPlane", &MixingPlane, PHSTRING, 2 * nTurboZone - 2);

    RDouble *xcc = gridUnstruct->GetCellCenterX();
    RDouble *ycc = gridUnstruct->GetCellCenterY();
    RDouble *zcc = gridUnstruct->GetCellCenterZ();

    RDouble *xfc = gridUnstruct->GetFaceCenterX();
    RDouble *yfc = gridUnstruct->GetFaceCenterY();
    RDouble *zfc = gridUnstruct->GetFaceCenterZ();

    RDouble *ns = gridUnstruct->GetFaceArea();

    int *leftCellofFace  = gridUnstruct->GetLeftCellOfFace();
    int *rightCellofFace = gridUnstruct->GetRightCellOfFace();

    int nBCRegion = unstructBCSet->GetnBCRegion();

    RDouble **nxsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("nxsSpan"));
    RDouble **ntsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("ntsSpan"));
    RDouble **nrsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("nrsSpan"));

    //! qSpanNeighbor: record flow variable in neighbor zone at each span.
    RDouble ***q_SpanCurrent = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("q_Span"));
    RDouble ***q_SpanTarget  = reinterpret_cast< RDouble *** > (gridTarget->GetDataPtr("q_Span"));

    RDouble ***qTrans_SpanCurrent = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTrans_Span"));
    RDouble ***qTransSpanNeighbor = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTransSpanNeighbor"));
    RDouble ***SpanTransFlux      = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanTransFlux"));

    RDouble ***qTrans_SpanTarget = reinterpret_cast< RDouble *** > (gridTarget->GetDataPtr("qTrans_Span"));

    RDouble ***dqTransSpanIn = reinterpret_cast <RDouble ***> (gridUnstruct->GetDataPtr("dqTransSpanIn"));
    RDouble ***dqTransSpanEx = reinterpret_cast <RDouble ***> (gridUnstruct->GetDataPtr("dqTransSpanEx"));
    RDouble ***dcTransSpan   = reinterpret_cast <RDouble ***> (gridUnstruct->GetDataPtr("dcTransSpan"));

    using namespace IDX;

    int FlowDirection = -1;

    int nEquation = GetNumberOfEquations();

    RDouble *ExtraDqTrans = new RDouble[nEquation];
    RDouble *deltaQTrans  = new RDouble[nEquation];
    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion++)
    {
        UnstructBC *bcRegion   = unstructBCSet->GetBCRegion(iBCRegion);
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        int bcType    = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        //! Parallel
        int iTurboZone = gridUnstruct->GetOrdinaryGridIndex();
        int TargetZone = gridTarget->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = gridUnstruct->GetZoneID();
            TargetZone = gridTarget->GetZoneID();
        }

        string MixingPlaneIn, MixingPlaneOut;

        if (iTurboZone == 0)
        {
            MixingPlaneOut = MixingPlane[0];
        }
        else if (iTurboZone == nTurboZone)
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
        }
        else
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
            MixingPlaneOut = MixingPlane[2 * iTurboZone];
        }

        //! TargetMixingFlag is Mixingplane in current zone, which need neighbor zone data.
        //! SourceMixingFlag is Mixingplane in neighbor zone, which provide the data. 
        int MixingFlag, TargetFlag;
        if (bcName == MixingPlaneIn || bcName == MixingPlaneOut)
        {
            //! downstream to upstream: transfer data from MixingIn to MixingOut.
            if (TargetZone == iTurboZone + 1)
            {
                if (bcName == MixingPlaneIn) continue;
                //! upstream MixingIn.
                TargetFlag = 0 + TargetZone * 10;
                //! current MixingOut.
                MixingFlag = 1 + iTurboZone * 10;
            }
            //! upstream to downstream: transfer data from MixingOut to MixingIn.
            else if (TargetZone == iTurboZone - 1)
            {
                if (bcName == MixingPlaneOut) continue;
                TargetFlag = 1 + TargetZone * 10;
                MixingFlag = 0 + iTurboZone * 10;
            }
            if (bcName == MixingPlaneIn || bcName == MixingPlaneOut)
            {
                //! compute variable perturbation
                for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
                {
                    for (int m = 0; m < nEquation; m++)
                    {
                        ExtraDqTrans[m] = 0.0;
                    }

                    //! interpolation
                    for (int m = 0; m < nEquation; m++)
                    {
                        dqTransSpanIn[MixingFlag][m][iSpan] = qTrans_SpanTarget[TargetFlag][m][iSpan] - qTrans_SpanCurrent[MixingFlag][m][iSpan];
                    }
                    //! extrapolation
                    vector<int> *faceIndex = bcRegion->GetFaceIndex();
                    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                    {
                        int iFace = *iter;

                        int re = rightCellofFace[iFace];
                        int le = leftCellofFace[iFace];

                        RDouble r_Face = sqrt(yfc[iFace] * yfc[iFace] + zfc[iFace] * zfc[iFace]);

                        //! extrapolating q[m][le] and q on the interface to obtain ExtraDq.
                        if (r_Face > r_target[MixingFlag][iSpan] && r_Face <= r_target[MixingFlag][iSpan + 1])
                        {
                            for (int m = 0; m < nEquation; m++)
                            {
                                ExtraDqTrans[m] = half * (qTrans_Span[MixingFlag][m][iSpan] - qTransition[m][le]) * ns[iFace];
                            }
                        }
                    }

                    for (int m = 0; m < nEquation; m++)
                    {
                        dqTransSpanEx[MixingFlag][m][iSpan] = ExtraDqTrans[m] / SpanArea[MixingFlag][iSpan];
                    }
                }
                for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
                {
                    //! after data transfer in NSSlover, IU, IV, IW refers to vx, vtheta, vr.
                    RDouble vn = nxsSpan[MixingFlag][iSpan] * q_SpanCurrent[MixingFlag][IU][iSpan]
                        + nrsSpan[MixingFlag][iSpan] * q_SpanCurrent[MixingFlag][IW][iSpan];

                    RDouble avgRho = half * (q_SpanCurrent[MixingFlag][IR][iSpan] + q_SpanTarget[TargetFlag][IR][iSpan]);
                    RDouble avgP = half * (q_SpanCurrent[MixingFlag][IP][iSpan] + q_SpanTarget[TargetFlag][IP][iSpan]);

                    RDouble c = sqrt(refGama * (avgP / avgRho));

                    RDouble rhoc = half * (q_SpanCurrent[MixingFlag][IR][iSpan] + q_SpanTarget[TargetFlag][IR][iSpan]) * c;

                    if (vn > 0)
                    {
                        for (int m = 0; m < nEquation; m++)
                        {
                            dcTransSpan[MixingFlag][m][iSpan] = rhoc * dqTransSpanEx[MixingFlag][m][iSpan];
                        }
                    }
                    else
                    {
                        for (int m = 0; m < nEquation; m++)
                        {
                            dcTransSpan[MixingFlag][m][iSpan] = rhoc * dqTransSpanIn[MixingFlag][m][iSpan];
                        }
                    }

                    for (int m = 0; m < nEquation; m++)
                    {
                        deltaQTrans[m] = dcTransSpan[MixingFlag][m][iSpan] / rhoc;
                    }

                    for (int m = 0; m < nEquation; m++)
                    {
                        qTransSpanNeighbor[MixingFlag][m][iSpan] = qTrans_SpanCurrent[MixingFlag][m][iSpan] + deltaQTrans[m];
                    }
                }
            }
        }
    }
    delete [] ExtraDqTrans;
    delete [] deltaQTrans;
}

//! function to set radial profile of mixingin and mixingout.
//! treat mixingin as a kind of inlet and treat mixing as a kind of outlet.
void TransitionSolverUnstr::SetMixingPlaneData(Grid *grid)
{
    UnstructGrid *gridUnstruct = UnstructGridCast(grid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_TransitionSolver *parameters = GetControlParameters();
    RDouble refPressure = parameters->GetRefDimensionalPressure();
    RDouble refGama     = parameters->GetRefGama();
    RDouble gama1       = refGama - 1;

    RDouble **q           = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q"));
    RDouble **qTransition = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q_transition"));

    RDouble **r_ori    = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Face"));
    RDouble **r_target = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Span"));

    //! nSpan is the number of spanwise point to do data average
    //! also interpolate data on these points.
    //! nSpan is set to 50 temporary.
    int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");
    RDouble nSpan = (RDouble)(nSpanSection);

    RDouble *nxs, *nys, *nzs, *ns;

    nxs = gridUnstruct->GetFaceNormalX();
    nys = gridUnstruct->GetFaceNormalY();
    nzs = gridUnstruct->GetFaceNormalZ();
    ns = gridUnstruct->GetFaceArea();

    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    string MixingPlane[100];
    GlobalDataBase::GetData("MixingPlane", &MixingPlane, PHSTRING, 2 * nTurboZone - 2);

    RDouble *xfc = gridUnstruct->GetFaceCenterX();
    RDouble *yfc = gridUnstruct->GetFaceCenterY();
    RDouble *zfc = gridUnstruct->GetFaceCenterZ();

    RDouble *xcc = gridUnstruct->GetCellCenterX();
    RDouble *ycc = gridUnstruct->GetCellCenterY();
    RDouble *zcc = gridUnstruct->GetCellCenterZ();

    int **face2nodeArray = gridUnstruct->GetFace2NodeArray();
    int * node_number_of_each_face = gridUnstruct->GetNodeNumberOfEachFace();

    int *leftCellofFace = gridUnstruct->GetLeftCellOfFace();
    int *rightCellofFace = gridUnstruct->GetRightCellOfFace();
    int nTotalCell = gridUnstruct->GetNTotalCell();

    int nBCRegion = unstructBCSet->GetnBCRegion();

    int nEquation = GetNumberOfEquations();

    using namespace IDX;

    RDouble ***qTrans_Span   = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTrans_Span"));
    RDouble ***qSpanNeighbor = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTransSpanNeighbor"));

    int FlowDirection = -1;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        //! Parallel
        int iTurboZone = gridUnstruct->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = gridUnstruct->GetZoneID();
        }

        //! use new bcType to define mixing plane.
        //! if (bcType == -1)
        string MixingPlaneIn = "";
        string MixingPlaneOut = "";

        if (iTurboZone == 0)
        {
            MixingPlaneOut = MixingPlane[0];
        }
        else if (iTurboZone == nTurboZone)
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
        }
        else
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
            MixingPlaneOut = MixingPlane[2 * iTurboZone];
        }

        if (bcName == MixingPlaneIn || bcName == MixingPlaneOut)
        {
            if (bcName == MixingPlaneIn)
            {
                FlowDirection = 0;
            }
            if (bcName == MixingPlaneOut)
            {
                FlowDirection = 1;
            }
            int MixingFlag = FlowDirection + iTurboZone * 10;

            for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
            {
                for (int m = 0; m < nEquation; m++)
                {
                    qTrans_Span[MixingFlag][m][iSpan] = qSpanNeighbor[MixingFlag][m][iSpan];
                }
            }

            for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
            {
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    //! iFace is the face number in the set of faceIndex.
                    int iFace = *iter;

                    int le = leftCellofFace[iFace];
                    int re = rightCellofFace[iFace];
                    //! find face center that locate in target span.
                    //! compute variable at each span.
                    if (r_ori[MixingFlag][iFace] > r_target[MixingFlag][iSpan] && r_ori[MixingFlag][iFace] <= r_target[MixingFlag][iSpan + 1])
                    {
                        for (int m = 0; m < nEquation; m++)
                        {
                            qTransition[m][re] = qTrans_Span[MixingFlag][m][iSpan];
                        }
                    }
                }
            }
        }
    }
}

void TransitionSolverUnstr::RotateVectorFromInterface(Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    RDouble **rotTransitiongradValueX = reinterpret_cast <RDouble **> (grid->GetDataPtr("rotTransitiongradValueX"));
    RDouble **rotTransitiongradValueY = reinterpret_cast <RDouble **> (grid->GetDataPtr("rotTransitiongradValueY"));
    RDouble **rotTransitiongradValueZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("rotTransitiongradValueZ"));

    int iNeighborZone = interfaceInformation->FindIthNeighbor(neighborZoneIndex);

    if (iNeighborZone = -1)
    {
        return;
    }

    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    int* interFace2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();

    RDouble rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle * PI / 180.0;

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");

    if (referenceFrame == ROTATIONAL_FRAME)
    {
        int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
        RDouble PeriodicRotationAngle[100];
        GlobalDataBase::GetData("PeriodicRotationAngle", &PeriodicRotationAngle, PHDOUBLE, nTurboZone);

        //! Parallel
        int iTurboZone = grid->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = grid->GetZoneID();
        }
        PeriodicRotationAngle[iTurboZone] = PeriodicRotationAngle[iTurboZone] * PI / 180.0;
        rotationAngle = PeriodicRotationAngle[iTurboZone];
    }

    if (nEquation > 0)
    {
        RDouble **fieldRecvTransitionY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTransitionY"));
        RDouble **fieldRecvTransitionZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTransitionZ"));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
            {
                int t1;
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                grid->GetTargetIndex(iFace, 1, t1);

                int iBFace = interFace2BoundaryFace[iFace];
                UnstructBC *bcregion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iBFace]);
                string bcName = bcregion->GetBCName();

                if (referenceFrame == ROTATIONAL_FRAME)
                {
                    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
                    string Periodic_Name[100];
                    GlobalDataBase::GetData("Periodic_Name", &Periodic_Name, PHSTRING, 2 * nTurboZone);

                    //! Parallel
                    int iTurboZone = grid->GetOrdinaryGridIndex();
                    //! Serial
                    if (iTurboZone == -1)
                    {
                        iTurboZone = grid->GetZoneID();
                    }

                    //! need loop over bcregions
                    if (bcName == Periodic_Name[2 * iTurboZone])
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {
                            rotTransitiongradValueY[m][t1] = fieldRecvTransitionY[m][t1] * cos(2 * PI - rotationAngle) - fieldRecvTransitionZ[m][t1] * sin(2 * PI - rotationAngle);
                            rotTransitiongradValueZ[m][t1] = fieldRecvTransitionY[m][t1] * sin(2 * PI - rotationAngle) + fieldRecvTransitionZ[m][t1] * cos(2 * PI - rotationAngle);
                        }

                        for (int m = 0; m < nEquation; ++m)
                        {
                            fieldRecvTransitionY[m][t1] = rotTransitiongradValueY[m][t1];
                            fieldRecvTransitionZ[m][t1] = rotTransitiongradValueZ[m][t1];
                        }
                    }
                    else if (bcName == Periodic_Name[2 * iTurboZone + 1])
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {
                            rotTransitiongradValueY[m][t1] = fieldRecvTransitionY[m][t1] * cos(rotationAngle) - fieldRecvTransitionZ[m][t1] * sin(rotationAngle);
                            rotTransitiongradValueZ[m][t1] = fieldRecvTransitionY[m][t1] * sin(rotationAngle) + fieldRecvTransitionZ[m][t1] * cos(rotationAngle);
                        }

                        for (int m = 0; m < nEquation; ++m)
                        {
                            fieldRecvTransitionY[m][t1] = rotTransitiongradValueY[m][t1];
                            fieldRecvTransitionZ[m][t1] = rotTransitiongradValueZ[m][t1];
                        }
                    }
                }

                else
                {
                    //! need loop over bcregions
                    if (bcName == "Periodic_up")
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {
                            rotTransitiongradValueY[m][t1] = fieldRecvTransitionY[m][t1] * cos(2 * PI - rotationAngle) - fieldRecvTransitionZ[m][t1] * sin(2 * PI - rotationAngle);
                            rotTransitiongradValueZ[m][t1] = fieldRecvTransitionY[m][t1] * sin(2 * PI - rotationAngle) + fieldRecvTransitionZ[m][t1] * cos(2 * PI - rotationAngle);
                        }

                        for (int m = 0; m < nEquation; ++m)
                        {
                            fieldRecvTransitionY[m][t1] = rotTransitiongradValueY[m][t1];
                            fieldRecvTransitionZ[m][t1] = rotTransitiongradValueZ[m][t1];
                        }
                    }
                    else if (bcName == "Periodic_down")
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {
                            rotTransitiongradValueY[m][t1] = fieldRecvTransitionY[m][t1] * cos(rotationAngle) - fieldRecvTransitionZ[m][t1] * sin(rotationAngle);
                            rotTransitiongradValueZ[m][t1] = fieldRecvTransitionY[m][t1] * sin(rotationAngle) + fieldRecvTransitionZ[m][t1] * cos(rotationAngle);
                        }

                        for (int m = 0; m < nEquation; ++m)
                        {
                            fieldRecvTransitionY[m][t1] = rotTransitiongradValueY[m][t1];
                            fieldRecvTransitionZ[m][t1] = rotTransitiongradValueZ[m][t1];
                        }
                    }
                }
            }
        }
    }
}

void TransitionSolverUnstr::GetGradientField(Grid *gridIn)
{
    UnstructGrid *gridUnstruct = UnstructGridCast(gridIn);
    int numberOfEquations = this->GetNumberOfEquations();

    RDouble **gradientTransitionX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTransitionX"));
    RDouble **gradientTransitionY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTransitionY"));
    RDouble **gradientTransitionZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTransitionZ"));
    RDouble **qTransitionNode = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("qTransitionNode"));
    RDouble **qTransition = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("q_transition"));

    for (int m = 0; m < numberOfEquations; ++ m)
    {
        string gradientName = GlobalDataBase::GetStrParaFromDB("gradientName");
        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGNode_NEW(qTransition[m], qTransitionNode[m], gradientTransitionX[m], gradientTransitionY[m], gradientTransitionZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(qTransition[m], gradientTransitionX[m], gradientTransitionY[m], gradientTransitionZ[m]);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(qTransition[m], gradientTransitionX[m], gradientTransitionY[m], gradientTransitionZ[m]);
        }
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(qTransition[m], gradientTransitionX[m], gradientTransitionY[m], gradientTransitionZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(qTransition[m], gradientTransitionX[m], gradientTransitionY[m], gradientTransitionZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(qTransition[m], gradientTransitionX[m], gradientTransitionY[m], gradientTransitionZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(qTransition[m], gradientTransitionX[m], gradientTransitionY[m], gradientTransitionZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }
}

LIB_EXPORT void TransitionSolverUnstr::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_TransitionSolver();
    controlParameters->Init();

    //! Compute machin zero. reference from cfl3d.
    machZeroExp = 1;
    RDouble add, x11;
    RDouble compare = 1.0;
    for (int iCount = 0; iCount < 20; ++ iCount)
    {
        add = 1.0;
        for (int n = 0; n < iCount + 1; ++ n)
        {
            add = add * 0.1;
        }

        x11 = compare + add;
        if (x11 == compare)
        {
            machZeroExp = iCount + 1;
            break;
        }
    }
}

LIB_EXPORT Param_TransitionSolver *TransitionSolverUnstr::GetControlParameters()
{
    return static_cast<Param_TransitionSolver *>(controlParameters);
}

}