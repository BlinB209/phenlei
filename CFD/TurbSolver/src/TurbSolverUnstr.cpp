#include "TurbSolverUnstr.h"
#include "Geo_UnstructBC.h"
#include "IO_FileName.h"
#include "MultiGridOperation.h"
#include "PHMatrix.h"
#include "TK_Exit.h"
#include "Param_TurbSolverUnstruct.h"
#include "Residual.h"
#include "FieldProxy.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "TK_Log.h"
#include "IO_FileName.h"
#include "Transition.h"
#include "Flux_Inviscid.h"
#include "Gas.h"
#ifdef USE_GMRESSOLVER
#include "Sacado.hpp"
#endif
using namespace std;

namespace PHSPACE
{
#ifdef USE_GMRESSOLVER
typedef Sacado::Fad::DFad<RDouble> ADReal;
#endif
TurbSolverUnstr::TurbSolverUnstr()
{
    DESLength = 0;
    gradientTurbField = 0;
    gradientVelocity = 0;
}

TurbSolverUnstr::~TurbSolverUnstr()
{
    DeAllocateGlobalVariables();

    delete []DESLength; 
    DESLength = 0;

    FreeControlParameters();
}

void TurbSolverUnstr::AllocateGlobalVar(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int viscousType = parameters->GetViscousType();
    if (nTurbulenceEquation == 0) return;

    int transitionType = parameters->GetTransitionType();

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotalFace    = grid->GetNTotalFace();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **qTurbulence        = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
    RDouble **residualTurbulence = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
    RDouble **spectrumTurbulence = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

    grid->UpdateDataPtr("q_turb"   , qTurbulence      );    //! Unknown variable of turbulent model.
    grid->UpdateDataPtr("res_turb" , residualTurbulence);    //! Residual or right-hand side.
    grid->UpdateDataPtr("spec_turb", spectrumTurbulence);    //! Spectral radius of turbulent model.

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");
    int nTurboZone = 0;
    if(referenceFrame)
    {
        nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    }

#ifdef USE_GMRESSOLVER
    //! GMRESturb
    int originaltscheme = GlobalDataBase::GetIntParaFromDB("OriginalTscheme");

    if( originaltscheme == GMRES )
    {

        vector<int> AI = grid->GetJacobianAI4GMRES();
        int TotalSize = AI[numberOfTotal];
        RDouble** dRdq_turb = NewPointer2<RDouble>(nTurbulenceEquation, nTurbulenceEquation * TotalSize);
        grid->UpdateDataPtr("dRdq_turb", dRdq_turb);
        PHSPACE::SetField(dRdq_turb, nTurbulenceEquation, nTurbulenceEquation * TotalSize, 0.0);

        // GMRESCoupled1st
        vector<int> AI1st   = grid->GetJacobianAI1st4GMRES(); 
        int JacOrder        = grid->GetJacobianOrder();
        if( JacOrder == 2 )
        {
            int TotalSize1st            = AI1st[numberOfTotal];
            RDouble** dRdq_turb1st      = NewPointer2<RDouble>(nTurbulenceEquation, nTurbulenceEquation * TotalSize1st);
            grid->UpdateDataPtr("dRdq_turb1st", dRdq_turb1st);
            PHSPACE::SetField(dRdq_turb1st, nTurbulenceEquation, nTurbulenceEquation * TotalSize1st, 0.0);
        }
        else
        {
            RDouble** dRdq_turb1st = NewPointer2<RDouble>(1,1);
            grid->UpdateDataPtr("dRdq_turb1st",dRdq_turb1st);
        }

    }


    RDouble **dDdP_turb = NewPointer2<RDouble>(nTurbulenceEquation, nTurbulenceEquation*numberOfBoundaryFace);
    grid->UpdateDataPtr("dDdP_turb", dDdP_turb);
    PHSPACE::SetField(dDdP_turb, nTurbulenceEquation, nTurbulenceEquation*numberOfBoundaryFace, 0.0);
#endif

    RDouble **gradientTurbulenceX = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
    RDouble **gradientTurbulenceY = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
    RDouble **gradientTurbulenceZ = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

    grid->UpdateDataPtr("gradTurbulenceX", gradientTurbulenceX);
    grid->UpdateDataPtr("gradTurbulenceY", gradientTurbulenceY);
    grid->UpdateDataPtr("gradTurbulenceZ", gradientTurbulenceZ);

    //! Temporary variable for periodic boundary condition.
    RDouble **rotNSgradValueX = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
    RDouble **rotNSgradValueY = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
    RDouble **rotNSgradValueZ = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
    grid->UpdateDataPtr("rotTurbgradValueX", rotNSgradValueX);
    grid->UpdateDataPtr("rotTurbgradValueY", rotNSgradValueY);
    grid->UpdateDataPtr("rotTurbgradValueZ", rotNSgradValueZ);

    //! allocate mixing plane variables.
    //! for multi-row turbomachinery.
    if (nTurboZone > 0)
    {
        int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");
        int nType = 2;

        //! zone with MixingIn and MixingOut.
        int nFlag = 10 * nTurboZone + 2;

        RDouble ***SpanTurbFlux      = NewPointer3<RDouble>(nFlag, nTurbulenceEquation, nSpanSection);
        RDouble ***qTurb_Span        = NewPointer3<RDouble>(nFlag, nTurbulenceEquation, nSpanSection);
        RDouble ***qTurbSpanNeighbor = NewPointer3<RDouble>(nFlag, nTurbulenceEquation, nSpanSection);
        RDouble ***dqTurbSpanIn      = NewPointer3<RDouble>(nFlag, nTurbulenceEquation, nSpanSection);
        RDouble ***dqTurbSpanEx      = NewPointer3<RDouble>(nFlag, nTurbulenceEquation, nSpanSection);
        RDouble ***dcTurbSpan        = NewPointer3<RDouble>(nFlag, nTurbulenceEquation, nSpanSection);
        InitMixingPlane(SpanTurbFlux, nFlag, nTurbulenceEquation, nSpanSection, zero);
        InitMixingPlane(qTurb_Span, nFlag, nTurbulenceEquation, nSpanSection, zero);
        InitMixingPlane(qTurbSpanNeighbor, nFlag, nTurbulenceEquation, nSpanSection, zero);
        InitMixingPlane(dqTurbSpanIn, nFlag, nTurbulenceEquation, nSpanSection, zero);
        InitMixingPlane(dqTurbSpanEx, nFlag, nTurbulenceEquation, nSpanSection, zero);
        InitMixingPlane(dcTurbSpan, nFlag, nTurbulenceEquation, nSpanSection, zero);
        grid->UpdateDataPtr("SpanTurbFlux", SpanTurbFlux);
        grid->UpdateDataPtr("qTurb_Span", qTurb_Span);
        grid->UpdateDataPtr("qTurbSpanNeighbor", qTurbSpanNeighbor);
        grid->UpdateDataPtr("dqTurbSpanIn", dqTurbSpanIn);
        grid->UpdateDataPtr("dqTurbSpanEx", dqTurbSpanEx);
        grid->UpdateDataPtr("dcTurbSpan", dcTurbSpan);
    }

    int numberOfTotalNode = grid->GetNTotalNode();
    RDouble **qTurbulenceNode = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotalNode);
    grid->UpdateDataPtr("qTurbNode", qTurbulenceNode);    //! Unknown variable of turbulent model at node.

    int *nodeValueSliceTurb = new int[numberOfTotalNode];
    grid->UpdateDataPtr("nodeValueSliceTurb", nodeValueSliceTurb);

    int *nodeValueSliceTurbTrade = new int[numberOfTotalNode];
    PHSPACE::SetField(nodeValueSliceTurbTrade, 0, numberOfTotalNode);
    grid->UpdateDataPtr("nodeValueSliceTurbTrade", nodeValueSliceTurbTrade);

    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation)
    {
        int numberOfInterPoint = interpointInformation->GetNumberOfInterpoints();
        RDouble **qTurbulenceInterPoint = NewPointer2<RDouble>(nTurbulenceEquation, numberOfInterPoint);
        grid->UpdateDataPtr("turb::qInterpoint", qTurbulenceInterPoint);    //! Unknown variable of turbulent model at interPoint.
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **qTurbulenceUnsteadyN1 = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
        RDouble **qTurbulenceUnsteadyN2 = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
        RDouble **residualTurbulenceUnsteadyN1 = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
        RDouble **residualTurbulenceUnsteadyN2 = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble **residualTurbulenceUnsteadyTemporary = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

        grid->UpdateDataPtr("q_turb_unsteady_n1", qTurbulenceUnsteadyN1);
        grid->UpdateDataPtr("q_turb_unsteady_n2", qTurbulenceUnsteadyN2);
        grid->UpdateDataPtr("res_turb_unsteady_n1", residualTurbulenceUnsteadyN1);
        grid->UpdateDataPtr("res_turb_unsteady_n2", residualTurbulenceUnsteadyN2);
        grid->UpdateDataPtr("res_turb_unsteady_tmp", residualTurbulenceUnsteadyTemporary);

        //! Statistical variables for DES simulation.
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        if (ifStaticsFlowField > 0)
        {
            RDouble **qAverageTurbulence = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
            grid->UpdateDataPtr("qAverageTurb", qAverageTurbulence);
        }
    }

    if (viscousType == TWO_EQU)
    {
        RDouble *cross = new RDouble[numberOfTotal];
        grid->UpdateDataPtr("cross", cross);    //! For the crossing term in two equation turbulence model.
        RDouble *blend = new RDouble[numberOfTotal];
        grid->UpdateDataPtr("blend", blend);    //! For blending function F1.
        RDouble *SST_F2 = new RDouble[numberOfTotal];
        grid->UpdateDataPtr("SST_F2", SST_F2);    //! For blending function F2.
        if (transitionType == IREGAMA)
        {
            RDouble *SpSdRatio = new RDouble[numberOfTotal];
            grid->UpdateDataPtr("SpSdRatio", SpSdRatio);

            RDouble *gamaeff = new RDouble[numberOfTotal];
            grid->UpdateDataPtr("gamaeff", gamaeff);
        }
    }

    int neasm    = GlobalDataBase::GetIntParaFromDB("neasm");

    if (viscousType >= TWO_EQU && neasm > 0)
    {
        int naniss = 5;
        RDouble **aniss = NewPointer2<RDouble>(naniss, numberOfTotal);

        grid->UpdateDataPtr("aniss", aniss);

        RDouble *sengy = new RDouble[numberOfTotal];

        grid->UpdateDataPtr("sengy", sengy);
    }

    RDouble **matrixTurbulenceLeft  = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotalFace);
    RDouble **matrixTurbulenceRight = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotalFace);

    grid->UpdateDataPtr("mat_turbl", matrixTurbulenceLeft);
    grid->UpdateDataPtr("mat_turbr", matrixTurbulenceRight);
}

void TurbSolverUnstr::DeAllocateGlobalVar(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int viscousType = parameters->GetViscousType();
    if (nTurbulenceEquation == 0) return;

    int transitionType = parameters->GetTransitionType();

    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **qTurbulence        = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"  ));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));

    DelPointer2(qTurbulence      );
    DelPointer2(residualTurbulence);
    DelPointer2(spectrumTurbulence);

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");
    int nTurboZone = 0;
    if (referenceFrame)
    {
        nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    }

#ifdef USE_GMRESSOLVER
    // GMRESCoupled1st
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");

    if( tscheme == GMRES )
    {
        RDouble** dRdq_turb     = reinterpret_cast<RDouble**> (grid->GetDataPtr("dRdq_turb"));
        RDouble** dRdq_turb1st  = reinterpret_cast<RDouble**> (grid->GetDataPtr("dRdq_turb1st"));
        DelPointer2(dRdq_turb      );
        DelPointer2(dRdq_turb1st);
    }
    
    RDouble** dDdP_turb     = reinterpret_cast<RDouble**> (grid->GetDataPtr("dDdP_turb"));
    DelPointer2(dDdP_turb);
#endif

    RDouble **gradientTurbulenceX = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTurbulenceX"));
    RDouble **gradientTurbulenceY = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTurbulenceY"));
    RDouble **gradientTurbulenceZ = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTurbulenceZ"));

    DelPointer2(gradientTurbulenceX);
    DelPointer2(gradientTurbulenceY);
    DelPointer2(gradientTurbulenceZ);

    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTurbNode"));
    DelPointer2(qTurbulenceNode);

    int *nodeValueSliceTurb = reinterpret_cast<int *> (grid->GetDataPtr("nodeValueSliceTurb"));
    int *nodeValueSliceTurbTrade = reinterpret_cast<int *> (grid->GetDataPtr("nodeValueSliceTurbTrade"));

    delete [] nodeValueSliceTurb;
    delete [] nodeValueSliceTurbTrade;

    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation)
    {
        RDouble **qTurbulenceInterPoint = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTurbulenceInterPoint"));
        DelPointer2(qTurbulenceInterPoint);
    }

    int isUnsteady = parameters->GetIsUnsteady();

    if (isUnsteady)
    {
        RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
        RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));
        RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
        RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
        RDouble **residualTurbulenceUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_tmp"));

        DelPointer2(qTurbulenceUnsteadyN1);
        DelPointer2(qTurbulenceUnsteadyN2);
        DelPointer2(residualTurbulenceUnsteadyN1);
        DelPointer2(residualTurbulenceUnsteadyN2);
        DelPointer2(residualTurbulenceUnsteadyTemporary);

        //! Statistical variables for DES simulation.
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        if (ifStaticsFlowField > 0)
        {
            RDouble **qAverageTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));
            DelPointer2(qAverageTurbulence);
        }
    }

    //! deallocate mixing plane array.
    if (nTurboZone > 0)
    {
        int nType = 2;
        int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");

        RDouble ***SpanTurbFlux      = reinterpret_cast <RDouble ***> (grid->GetDataPtr("SpanTurbFlux"));
        RDouble ***qTurb_Span        = reinterpret_cast <RDouble ***> (grid->GetDataPtr("qTurb_Span"));
        RDouble ***qTurbSpanNeighbor = reinterpret_cast <RDouble ***> (grid->GetDataPtr("qTurbSpanNeighbor"));
        RDouble ***dqTurbSpanIn      = reinterpret_cast <RDouble ***> (grid->GetDataPtr("dqTurbSpanIn"));
        RDouble ***dqTurbSpanEx      = reinterpret_cast <RDouble ***> (grid->GetDataPtr("dqTurbSpanEx"));
        RDouble ***dcTurbSpan        = reinterpret_cast <RDouble ***> (grid->GetDataPtr("dcTurbSpan"));

        DelPointer3(SpanTurbFlux);
        DelPointer3(qTurb_Span);
        DelPointer3(qTurbSpanNeighbor);
        DelPointer3(dqTurbSpanIn);
        DelPointer3(dqTurbSpanEx);
        DelPointer3(dcTurbSpan);
    }

    if (viscousType == TWO_EQU)
    {
        RDouble *cross  = reinterpret_cast<RDouble *> (grid->GetDataPtr("cross"));
        delete []cross;
        RDouble *blend  = reinterpret_cast<RDouble *> (grid->GetDataPtr("blend"));
        delete []blend;
        RDouble *SST_F2 = reinterpret_cast<RDouble *> (grid->GetDataPtr("SST_F2"));
        delete []SST_F2;
        if (transitionType == IREGAMA)
        {
            RDouble *SpSdRatio = reinterpret_cast<RDouble *> (grid->GetDataPtr("SpSdRatio"));
            delete [] SpSdRatio;
            RDouble *gamaeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("gamaeff"));
            delete [] gamaeff;
        }
    }

    int neasm    = GlobalDataBase::GetIntParaFromDB("neasm");

    if (viscousType >= TWO_EQU && neasm > 0)
    {
        RDouble **aniss = reinterpret_cast<RDouble **> (grid->GetDataPtr("aniss"));
        DelPointer2(aniss);

        RDouble *sengy = reinterpret_cast<RDouble * > (grid->GetDataPtr("sengy"));
        delete []sengy;
    }

    RDouble **matrixTurbulenceLeft  = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));

    DelPointer2(matrixTurbulenceLeft);
    DelPointer2(matrixTurbulenceRight);
}

bool TurbSolverUnstr::JudgeIfRestart()
{
    string turbFile = ".\results\turb.dat";
    GlobalDataBase::GetData("turbfile", &turbFile, PHSTRING, 1);

    if (PHMPI::IsParallelRun())
    {
        turbFile = PHSPACE::AddSymbolToFileName(turbFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(turbFile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

bool TurbSolverUnstr::JudgeIfReadAverage()
{
    string turbFile = ".\results\turb.dat";
    GlobalDataBase::GetData("turbfile", &turbFile, PHSTRING, 1);
    string averageTurbFile = AddSymbolToFileName(turbFile, "_", "Average");

    if (PHMPI::IsParallelRun())
    {
        averageTurbFile = PHSPACE::AddSymbolToFileName(averageTurbFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(averageTurbFile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

void TurbSolverUnstr::DumpRestart(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int outerIterationStep = GlobalDataBase::GetIntParaFromDB("outnstep");

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();
    cdata->Write(&outerIterationStep, sizeof(int));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->Write(qTurbulence[m], sizeof(RDouble) * numberOfTotal);
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->Write(qTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->Write(qTurbulenceUnsteadyN2[m], sizeof(RDouble) * numberOfTotal);
    }

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
    RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->Write(residualTurbulence[m], sizeof(RDouble) * numberOfTotal);
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->Write(residualTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->Write(residualTurbulenceUnsteadyN2[m], sizeof(RDouble) * numberOfTotal);
    }
}

void TurbSolverUnstr::DumpRestartH5(ActionKey *actkey)
{
    using namespace PHMPI;
    int currentProcessorID = GetCurrentProcessorID();

    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    int version = 1;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
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

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    double **qTemporary  = NewPointer2<double>(nTurbulenceEquation, numberOfTotal);

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(qTurbulence[iEquation][iCell]);
        }
    }
    WriteData(grploc, qTemporary[0], "q_turb");
    WriteData(grploc, &numberOfTotalCell, "nTotalCell");

    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    double *qTemporary2  = new double[numberOfTotal];
    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
    {
        qTemporary2[iCell] = static_cast<double>(viscousTurbulence[iCell]);
    }

    WriteData(grploc, qTemporary2, "visturb");

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        H5Gclose(grploc);
        DelPointer2(qTemporary);

        return;
    }

    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(qTurbulenceUnsteadyN1[iEquation][iCell]);
        }
    }

    WriteData(grploc, qTemporary[0], "q_turb_unsteady_n1");

    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(qTurbulenceUnsteadyN2[iEquation][iCell]);
        }
    }

    WriteData(grploc, qTemporary[0], "q_turb_unsteady_n2");

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(residualTurbulence[iEquation][iCell]);
        }
    }

    WriteData(grploc, qTemporary[0], "res_turb");

    RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(residualTurbulenceUnsteadyN1[iEquation][iCell]);
        }
    }

    WriteData(grploc, qTemporary[0], "res_turb_unsteady_n1");

    RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = static_cast<double>(residualTurbulenceUnsteadyN2[iEquation][iCell]);
        }
    }

    WriteData(grploc, qTemporary[0], "res_turb_unsteady_n2");

    if (IsNeedStatistics())
    {
        RDouble **qAverageTurb = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));
        for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
            {
                qTemporary[iEquation][iCell] = static_cast<double>(qAverageTurb[iEquation][iCell]);
            }
        }

        WriteData(grploc, qTemporary[0], "qAverageTurb");
    }

    H5Gclose(grploc);
    DelPointer2(qTemporary);
    delete qTemporary2;    qTemporary2 = nullptr;
}

void TurbSolverUnstr::ReadRestartH5(ActionKey *actkey)
{
    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulence = reinterpret_cast<RDouble **>(grid->GetDataPtr("q_turb"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    double **qTemporary = NewPointer2<double>(nTurbulenceEquation, numberOfTotal);

    int outerIterationStepOfTurbulence = 0;
    ReadData(actkey->filepos, &outerIterationStepOfTurbulence, "outnstep");

    int outerIterationStepOfNS = GlobalDataBase::GetIntParaFromDB("outnstep");

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    int nTotalCellRestart = 0;
    ReadData(grploc, &nTotalCellRestart, "nTotalCell");
    if(nTotalCellRestart != numberOfTotalCell)
    {
        ostringstream erroeInfo;
        erroeInfo << " Error: the cell number in turb.dat is not equal to the cell number in grid file !" << endl;
        TK_Exit::ExceptionExit(erroeInfo.str());
    }

    ReadData(grploc, qTemporary[0], "q_turb");
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulence[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
        }
    }
    if (CheckDataExist(grploc, "visturb"))
    {
        double *qTemporary1d  = new double[numberOfTotal];
        ReadData(grploc, qTemporary1d, "visturb");
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            viscousTurbulence[iCell] = static_cast<RDouble>(qTemporary1d[iCell]);
        }
        delete qTemporary1d;
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
        RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

        RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
        RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
        RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
        RDouble **qAverageTurb = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));

        //! Start from steady flow field and reset the outer step when start from steady flow.
        int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();
        if (ifStartFromSteadyResults)
        {
            PrintToWindow("Restart from steady Turbulent flow field, reset outer step to be zero!\n");

            outerIterationStepOfTurbulence = 0;
            GlobalDataBase::UpdateData("outnstep", &outerIterationStepOfTurbulence, PHINT, 1);

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
                {
                    qTurbulenceUnsteadyN1[m][iCell] = qTurbulence[m][iCell];
                    qTurbulenceUnsteadyN2[m][iCell] = qTurbulence[m][iCell];

                    residualTurbulence[m][iCell] = 0.0;
                    residualTurbulenceUnsteadyN1[m][iCell] = 0.0;
                    residualTurbulenceUnsteadyN2[m][iCell] = 0.0;
                }
            }
        }
        else
        {
            ReadData(grploc, qTemporary[0], "q_turb_unsteady_n1");
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    qTurbulenceUnsteadyN1[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                }
            }

            ReadData(grploc, qTemporary[0], "q_turb_unsteady_n2");
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    qTurbulenceUnsteadyN2[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                }
            }

            ReadData(grploc, qTemporary[0], "res_turb");
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    residualTurbulence[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                }
            }

            ReadData(grploc, qTemporary[0], "res_turb_unsteady_n1");
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    residualTurbulenceUnsteadyN1[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                }
            }

            ReadData(grploc, qTemporary[0], "res_turb_unsteady_n2");
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    residualTurbulenceUnsteadyN2[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                }
            }
        }

        int nStatisticalStep = 0;
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
        bool isReadAverageFlow = ifStaticsFlowField && (outerIterationStepOfTurbulence >= startStatisticStep);
        if (isReadAverageFlow)
        {
            if (ifStartFromSteadyResults)
            {
                for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
                {
                    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                    {
                        qAverageTurb[iEquation][iCell] = static_cast<RDouble>(qTemporary[iEquation][iCell]);
                    }
                }

                GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            }
            else
            {
                ReadData(grploc, qTemporary[0], "qAverageTurb");
                for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
                {
                    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                    {
                        qAverageTurb[iEquation][iCell] = qTemporary[iEquation][iCell];
                    }
                }

                ReadData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
                GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            }
        }
    }

    H5Gclose(grploc);
    DelPointer2(qTemporary);

    CompareOutStepOfFlowfieldFile(outerIterationStepOfNS, outerIterationStepOfTurbulence);
}

void TurbSolverUnstr::ReadRestart(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();
    int outerIterationStepOfTurbulence = 0;
    cdata->Read(&outerIterationStepOfTurbulence, sizeof(int));

    int outerIterationStepOfNS = GlobalDataBase::GetIntParaFromDB("outnstep");

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            cdata->Read(&qTurbulence[m][iCell], sizeof(RDouble));
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
        RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

        RDouble **residualTurbulence           = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
        RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
        RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));

        int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();

        if (ifStartFromSteadyResults)
        {
            PrintToWindow("Restart from steady Turbulent flow field, reset outer step to be zero!\n");

            //! Start from steady flow field and reset the outer step when start from steady flow.
            outerIterationStepOfTurbulence = 0;
            GlobalDataBase::UpdateData("outnstep", &outerIterationStepOfTurbulence, PHINT, 1);

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
                {
                    qTurbulenceUnsteadyN1[m][iCell] = qTurbulence[m][iCell];
                    qTurbulenceUnsteadyN2[m][iCell] = qTurbulence[m][iCell];

                    residualTurbulence          [m][iCell] = 0.0;
                    residualTurbulenceUnsteadyN1[m][iCell] = 0.0;
                    residualTurbulenceUnsteadyN2[m][iCell] = 0.0;
                }
            }
        }
        else
        {
            //! Start from unsteady flow field.
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                cdata->Read(qTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
            }
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                cdata->Read(qTurbulenceUnsteadyN2[m], sizeof(RDouble) * numberOfTotal);
            }

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                cdata->Read(residualTurbulence[m], sizeof(RDouble) * numberOfTotal);
            }
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                cdata->Read(residualTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
            }
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                cdata->Read(residualTurbulenceUnsteadyN2[m], sizeof(RDouble) * numberOfTotal);
            }
        }
    }

    CompareOutStepOfFlowfieldFile(outerIterationStepOfNS, outerIterationStepOfTurbulence);
}

void TurbSolverUnstr::InitFlowAsRestart()
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    if (nTurbulenceEquation == 0) return;

    UnstructGrid *grid = UnstructGridCast(GetGrid(0));

    //int outerIterationStep = 0;
    //GlobalDataBase::UpdateData("outnstep", &outerIterationStep, PHINT, 1);

    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulence[m][iCell] = freeStreamTurbVar[m];
        }
    }

    RDouble *SpSdRatio = reinterpret_cast<RDouble *> (grid->GetDataPtr("SpSdRatio"));
    RDouble *gamaeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("gamaeff"));
    int transitionType = parameters->GetTransitionType();
    if (transitionType == IREGAMA)
    {
        SetField(SpSdRatio, 1.0, numberOfTotal);
        SetField(gamaeff, 0.0, numberOfTotal);
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulenceUnsteadyN1[m][iCell] = qTurbulence[m][iCell];
        }
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulenceUnsteadyN2[m][iCell] = qTurbulence[m][iCell];
        }
    }

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
    RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
    RDouble **residualTurbulenceUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_tmp"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTurbulence[m][iCell] = 0.0;
            residualTurbulenceUnsteadyN1[m][iCell] = 0.0;
            residualTurbulenceUnsteadyN2[m][iCell] = 0.0;
            residualTurbulenceUnsteadyTemporary[m][iCell] = 0.0;
        }
    }
}

void TurbSolverUnstr::ZeroResiduals(Grid *gridIn)
{
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(gridIn);

    int numberOfGeneralCells = grid->GetNTotalCell() + grid->GetNBoundFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **>(grid->GetDataPtr("res_turb"));

    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfGeneralCells; ++ iCell)
        {
            residualTurbulence[iEquation][iCell] = 0.0;
        }
    }

#ifdef USE_GMRESSOLVER
    //! GMRESturb zero residual

    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");

    if( tscheme == GMRES )
    {
        RDouble** dRdq_turb = reinterpret_cast<RDouble**> (grid->GetDataPtr("dRdq_turb"));
        //! GMRESturb CSR
        vector<int> AI = grid->GetJacobianAI4GMRES();
        PHSPACE::SetField(dRdq_turb, nTurbulenceEquation, nTurbulenceEquation * AI[numberOfGeneralCells], 0.0);
        //! GMRESCoupled
        int nLaminar = grid->GetBlockSizeOfJacobianMatrix();
        nLaminar -= nTurbulenceEquation;
        RDouble **dRdqCoupledTerm_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dRdqCoupledTerm_turb"));
        PHSPACE::SetField(dRdqCoupledTerm_turb, nTurbulenceEquation, nLaminar * AI[numberOfGeneralCells], 0.0);
    }
#endif

    return;
}

void TurbSolverUnstr::InitSpectrum(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    if (nTurbulenceEquation == 0) return;

    RDouble turbCFLScale = parameters->GetTurbCFLScale();

    RDouble *vol  = grid->GetCellVolume();
    RDouble *dt = reinterpret_cast<RDouble *> (grid->GetDataPtr("dt"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));

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
        dualTimeSpectrumC2 =   dualTimeCoefficient[6];
    }

    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            spectrumTurbulence[m][iCell] = 0.0;
        }
    }

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            spectrumTurbulence[m][iCell] = dualTimeSpectrumC1 * vol[iCell] + dualTimeSpectrumC2 / (turbCFLScale * dt[iCell] + SMALL);
        }
    }
}

void TurbSolverUnstr::CorrectFineGrid(Grid *fineGrid, Grid *coarseGrid)
{
    int numberOfTotalCellOnFineGrid = UnstructGridCast(fineGrid)->GetNTotalCell();
    int *cell2CoarseGridCell = UnstructGridCast(fineGrid)->GetCell2CoarseGridCell();

    RDouble **qTurbulenceFineGrid = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell <numberOfTotalCellOnFineGrid; ++ iCell)
        {
            qTurbulenceFineGrid[m][iCell] += qTurbulenceCoarseGrid[m][cell2CoarseGridCell[iCell]];
        }
    }
}

void TurbSolverUnstr::InterpolatFineGrid(Grid *fineGridIn, Grid *coarseGridIn)
{
    UnstructGrid *fineGrid = UnstructGridCast(fineGridIn);
    UnstructGrid *coarseGrid = UnstructGridCast(coarseGridIn);

    RDouble **qTurbulenceFineGrid = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        InterpolatQ(fineGrid, qTurbulenceFineGrid[m], coarseGrid, qTurbulenceCoarseGrid[m]);
    }

    RDouble *viscousOfFineGrid = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("vist"));
    RDouble *viscousOfCoarseGrid = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("vist"));
    InterpolatQ(fineGrid, viscousOfFineGrid, coarseGrid, viscousOfCoarseGrid);
}

void TurbSolverUnstr::PutCorrectionBack(Grid *gridIn, FieldProxy *qProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **q = qProxy->GetField_UNS();
    RDouble **qOld = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qOld[m][iCell] += q[m][iCell];
        }
    }
}

void TurbSolverUnstr::RestrictDefect(Grid *fineGridIn, Grid *coarseGridIn)
{
    UnstructGrid *fineGrid = UnstructGridCast(fineGridIn);
    UnstructGrid *coarseGrid = UnstructGridCast(coarseGridIn);

    int numberOfTotalCellOnFineGrid = fineGrid->GetNTotalCell();
    int *cell2CoarseGridCell = fineGrid->GetCell2CoarseGridCell();

    ZeroResiduals(coarseGrid);

    RDouble **residualTuebulenceCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("res_turb"));
    RDouble **residualTuebulenceFineGrid = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("res_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell <numberOfTotalCellOnFineGrid; ++ iCell)
        {
            residualTuebulenceCoarseGrid[m][cell2CoarseGridCell[iCell]] -= residualTuebulenceFineGrid[m][iCell];
        }
    }
}

void TurbSolverUnstr::LoadQ(Grid *gridIn, FieldProxy *qProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **qqTurbulence = qProxy->GetField_UNS();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qqTurbulence[m][iCell] = qTurbulence[m][iCell];
        }
    }
}

void TurbSolverUnstr::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{
    UnstructGrid *grid  = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble **qPrimitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence  = qProxy->GetField_UNS();
    RDouble **dqTurbulence = dqProxy->GetField_UNS();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int viscousType = parameters->GetViscousType();
    int DESType = parameters->GetDESType();

    using namespace IDX;

    RDouble *primitiveVariable = new RDouble [nTurbulenceEquation];

    if (viscousType == ONE_EQU)
    {
        //! Spalart-allmaras model.

        RDouble nuoo  = 1.0;    //!Nondimensional inflow laminar viscosity / rho.
        RDouble nuMin = 1.0e-5 * nuoo;
        RDouble nuMax = 1.0e10 * nuoo;
        int newNstep = GlobalDataBase::GetIntParaFromDB("newnstep");
        int qTurbLimit = GlobalDataBase::GetIntParaFromDB("qTurbLimit");
        if(newNstep < 1000)
        {
            nuMax = 1.0e5 * nuoo;
        }

        if(DESType == RANS && qTurbLimit == 0)
        {
            for (int iCell = 0; iCell < numberOfTotalCell; ++iCell)
            {
                qTurbulence[ISA][iCell] += dqTurbulence[ISA][iCell];
            }
        }
        else
        {
            for (int iCell = 0; iCell < numberOfTotalCell; ++iCell)
            {
                qTurbulence[ISA][iCell] += dqTurbulence[ISA][iCell];
                qTurbulence[ISA][iCell] = MAX(qTurbulence[ISA][iCell], nuMin);
                qTurbulence[ISA][iCell] = MIN(qTurbulence[ISA][iCell], nuMax);
            }
        }
    }
    else if (viscousType == TWO_EQU)
    {
        RDouble *freeStreamTurbVar = parameters->GetFreeStreamTurbVar();

        RDouble kelim = 1.0e-5 * freeStreamTurbVar[IKE];
        RDouble kwlim = 1.0e-5 * freeStreamTurbVar[IKW];

        for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
        {
            RDouble ke,kw,dke,dkw,rho;
            rho = qPrimitiveVariable[IR][iCell];
            ke  = qTurbulence[IKE][iCell];
            kw  = qTurbulence[IKW][iCell]; 
            dke = dqTurbulence[IKE][iCell];
            dkw = dqTurbulence[IKW][iCell];
            ke += dke/rho;
            kw += dkw/rho;
            ke = MAX(ke,kelim);
            kw = MAX(kw,kwlim);
            qTurbulence[IKE][iCell] = ke;
            qTurbulence[IKW][iCell] = kw;
        }
    }

    delete [] primitiveVariable;    primitiveVariable = nullptr;
}

void TurbSolverUnstr::SMoothTurbulencePoint(Grid *gridIn, RDouble *primitiveVariable, int i)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    int le,re,ie,face;
    RDouble volumeSum = zero;
    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        primitiveVariable[m] = zero;
    }

    for (int j = 0; j < faceNumberOfEachCell[i]; ++ j)
    {
        face = cell2Face[i][j];
        le   = leftCellOfFace[face];
        re   = rightCellOfFace[face];

        if (le != i)
        {
            ie = le;
        }
        else
        {
            ie = re;
        }

        if (ie >= numberOfTotalCell) continue;

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            primitiveVariable[m] += qTurbulence[m][ie];
        }
        volumeSum += 1.0;
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        primitiveVariable[m] /= volumeSum;
    }
}

void TurbSolverUnstr::SMoothTurbulence(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace  = grid->GetRightCellOfFace();

    RDouble *volume = grid->GetCellVolume() ;

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **tmporary = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qTurbulence[m][iCell] *= volume[iCell];
            tmporary[m][iCell]     = qTurbulence[m][iCell];
        }
    }

    RDouble volumeSum;

    int le,re,ie,face;
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        volumeSum = volume[iCell];
        for (int j = 0; j < faceNumberOfEachCell[iCell]; ++ j)
        {
            face = cell2Face[iCell][j];
            le    = leftCellOfFace [face];
            re    = rightCellOfFace[face];
            
            if (le != iCell)
            {
                ie = le;
            }
            else
            {
                ie = re;
            }

            if (ie >= numberOfTotalCell) continue;

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qTurbulence[m][iCell] += tmporary[m][ie];
            }
            volumeSum += volume[ie];
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qTurbulence[m][iCell] /= volumeSum;
        }
    }

    DelPointer2(tmporary);
}

void TurbSolverUnstr::RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coefficient)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble **dqTurbulence = dqProxy->GetField_UNS();

    RDouble *dt = reinterpret_cast<RDouble *> (grid->GetDataPtr("dt"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            dqTurbulence[m][iCell] = dt[iCell] * coefficient * dqTurbulence[m][iCell];
        }
    }
}

void TurbSolverUnstr::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **dq = dqProxy->GetField_UNS();
    int *iBlank = grid->GetBlankIndex();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble *df    = new RDouble [nTurbulenceEquation];
    RDouble *rhs0  = new RDouble [nTurbulenceEquation];
    RDouble *dqOld = new RDouble [nTurbulenceEquation];

    RDouble **spectrumTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft  = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));
    RDouble **residualTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();
    RDouble *matrixUx = new RDouble [nTurbulenceEquation];

    //! Now the Forward Sweep
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        if (iBlank[iCell] != ACTIVE)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dq[m][iCell] = 0.0;
            }
            continue;
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            rhs0[m] = 0.0;
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            //! Back up the old dq, to compute the convergency.
            dqOld[m] = dq[m][iCell];

            //! Here, the deltaFlux is computed in Upper backward sweep!
            //! res: b      deltaFlux: Ux, which is computed in backward sweep.
            //! the dq is not the real dq, now.
            //! Although the 'dq' changed can not the right one, it dosen't matter, since 
            //! the following only using the lower neighbor cells, whose 'dq' has been updated.
            //! It is convenient for precondition transform.
            dq[m][iCell] = residualTurbulence[m][iCell];
            matrixUx[m]  = deltaFlux[m][iCell];

            //! Then reset it to zero to store the Lower forward sweep!
            deltaFlux[m][iCell] = 0.0;
        }

        for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++ jFace)
        {
            int face = cell2Face[iCell][jFace];
            int le   = leftCellOfFace [face];
            int re   = rightCellOfFace [face];
            
            if (iBlank[le] == INTERPOLATION || iBlank[re] == INTERPOLATION) continue; 

            //! One of le and re must be cell itself.
            if (le > iCell || re > iCell) continue;

            if (re == iCell)
            {
                //! Now its neighboring cell belongs to lower triangular
                SWAP(le, re);

                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    df[m] = matrixTurbulenceLeft[m][face] * dq[m][re];
                }
            }
            else
            {
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    df[m] = matrixTurbulenceRight[m][face] * dq[m][re];
                }
            }

            //! Add Flux together
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            //! dq = { (b - Ux) - Lx } / D.
            //! Note: the 'dq' has actually initialized by the residual at the beginning of this function.
            //! the 'rhs0' is the total Delta Flux of the neighbors.
            //! rhs0: Lx    diagonal: D.
            dq[m][iCell] = (dq[m][iCell] - matrixUx[m] - rhs0[m]) / spectrumTurbulence[m][iCell];

            //! Store the lower forward sweep delta-flux, which will be used in the backward sweep.
            deltaFlux[m][iCell] += rhs0[m];

            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    }

    delete [] rhs0;    rhs0 = nullptr;
    delete [] df;    df = nullptr;
    delete [] dqOld;    dqOld = nullptr;
    delete [] matrixUx;    matrixUx = nullptr;
}

void TurbSolverUnstr::SolveLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **dq = dqProxy->GetField_UNS();
    int *iBlank = grid->GetBlankIndex();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble *df    = new RDouble [nTurbulenceEquation];
    RDouble *rhs0  = new RDouble [nTurbulenceEquation];
    RDouble *dqOld = new RDouble [nTurbulenceEquation];

    RDouble **spectrumTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft  = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));
    RDouble **residualTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();
    RDouble *matrixLx = new RDouble [nTurbulenceEquation];

    //! Backward Sweep
    for (int iCell = numberOfTotalCell - 1; iCell >= 0; -- iCell)
    {
        if (iBlank[iCell] != ACTIVE)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dq[m][iCell] = 0.0;
            }
            continue;
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            rhs0[m] = 0.0;
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            //! Back up the old dq, to compute the convergence.
            //! the 'dq' is dq*, which has been updated in forward.
            dqOld[m] = dq[m][iCell];

            //! the dq is not the real dq, now.
            //! it is convenient for precondition transform.
            dq[m][iCell] = residualTurbulence[m][iCell];

            matrixLx[m] = deltaFlux[m][iCell];

            deltaFlux[m][iCell] = 0.0;
        }

        for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++ jFace)
        {
            int face  = cell2Face[iCell][jFace];
            int le    = leftCellOfFace  [face];
            int re    = rightCellOfFace [face];

            if (iBlank[le] == INTERPOLATION || iBlank[re] == INTERPOLATION) continue; 

            //! One of le and re must be cell itself.
            if (le < iCell || re < iCell) continue;

            if (re == iCell)
            {
                //! Now its neighboring cell belongs to upper triangular
                SWAP(le, re);

                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    df[m] = matrixTurbulenceLeft[m][face] * dq[m][re];
                }
            }
            else
            {
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    df[m] = matrixTurbulenceRight[m][face] * dq[m][re];
                }
            }

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            //! Note: the 'dq' has been updated by the forward sweep.
            //! the 'rhs0' is the total Delta Flux of the neighbors.
            //! x = {(b - LX) - Ux} / D.
            //! rhs0: Ux.    diagonal: D.
            dq[m][iCell] = (dq[m][iCell] - matrixLx[m] - rhs0[m]) / spectrumTurbulence[m][iCell];

            //! Store the upper backward sweep delta-flux, which will be used in the forward sweep.
            deltaFlux[m][iCell] += rhs0[m];

            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    }

    delete [] rhs0;    rhs0 = nullptr;
    delete [] df;    df = nullptr;
    delete [] dqOld;    dqOld = nullptr;
    delete [] matrixLx;    matrixLx = nullptr;
}

//! Bell 20130401 add
void TurbSolverUnstr::SetGhostDQLUSGS(Grid *gridIn, RDouble **dq)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    if (nTurbulenceEquation != 1)
    {
        cout << "SetGhostDQLUSGS function is turned off !\n";
        return;
    }

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *rightCellOfFace  = grid->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (IsInterface(bcType) || bcType == PHENGLEI::OVERSET)
        {
            continue;
        }

        if (IsWall(bcType))
        {
            SetWallBCGhostDQLUSGS(grid, bcRegion, dq);
        }
        //! FARFIELD, SYMMETRY, INFLOW and OUTFLOW can use the same module.
        else if (bcType == PHENGLEI::FARFIELD || bcType == PHENGLEI::SYMMETRY 
            || bcType == PHENGLEI::INFLOW   || bcType == PHENGLEI::OUTFLOW)
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
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    dq[m][re] = 0.0; 
                }
            }
        }
    }
}

void TurbSolverUnstr::SetWallBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *rightCellOfFace  = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        int iFace = *iter;
        int re = rightCellOfFace[iFace];
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            dq[m][re] = 0.0;
        }
    }
}

void TurbSolverUnstr::SetFarfieldBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace  = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            dq[m][re] = dq[m][le]; 
        }
    }
}


void TurbSolverUnstr::ModifyResiduals(Grid *grid)
{
    ModifyResidualOnWallBoundary(grid);
}

void TurbSolverUnstr::LoadResiduals(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **rhs = rhsProxy->GetField_UNS();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTurbulence[m][iCell] = - rhs[m][iCell];
        }
    }
}

RDouble TurbSolverUnstr::UnsteadyConvergence(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return zero;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **qPrimitive  = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble *primitive0 = new RDouble[nTurbulenceEquation];
    RDouble *primitive1 = new RDouble[nTurbulenceEquation];
    RDouble *primitive2 = new RDouble[nTurbulenceEquation];

    RDouble *qConservative0 = new RDouble[nTurbulenceEquation];
    RDouble *qConservative1 = new RDouble[nTurbulenceEquation];
    RDouble *qConservative2 = new RDouble[nTurbulenceEquation];

    RDouble sum1 = zero;
    RDouble sum2 = zero;

    using namespace IDX;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            primitive0[m] = qTurbulence[m][iCell];
            primitive1[m] = qTurbulenceUnsteadyN1[m][iCell];
            primitive2[m] = qTurbulenceUnsteadyN2[m][iCell];
        }

        RDouble coefficient = 1.0;

        if (nTurbulenceEquation >= 2)
        {
            coefficient = qPrimitive[IR][iCell];
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qConservative0[m] = coefficient * primitive0[m];
            qConservative1[m] = coefficient * primitive1[m];
            qConservative2[m] = coefficient * primitive2[m];
        }

        //! nl is needed.
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            //! Here, res must be transformed into dq, res is not equal to rhs but dq.
            RDouble dqp = residualTurbulence[m][iCell];        // qn+1,p+1 - qn+1,p
            RDouble dqn = qConservative0[m] - qConservative1[m];  // qn+1,p+1 - qn
            sum1 += dqp * dqp;
            sum2 += dqn * dqn;
        }
    }

    delete []primitive0;
    delete []primitive1;
    delete []primitive2;
    delete []qConservative0;
    delete []qConservative1;
    delete []qConservative2;

    RDouble cvg = sqrt(ABS(sum1 / (sum2 + SMALL)));
    return cvg;
}

void TurbSolverUnstr::UpdateUnsteadyFlow(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

    RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
    RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
    RDouble **residualTurbulenceUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_tmp"));

    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulenceUnsteadyN2[m][iCell] = qTurbulenceUnsteadyN1[m][iCell];
        }
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulenceUnsteadyN1[m][iCell] = qTurbulence[m][iCell];
        }
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTurbulenceUnsteadyN2[m][iCell] = residualTurbulenceUnsteadyN1[m][iCell];
        }
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTurbulenceUnsteadyN1[m][iCell] = residualTurbulenceUnsteadyTemporary[m][iCell]; //! Here the current outerIterationStep is over, the value of the stored resTmp should be assigned to resn1 for the next outerIterationStep. It should be noticed that resTmp only contain the inviscid and viscous flux.
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
            if(statisticalTimePeriod > 0.0)
            {
                c1 = MAX(c1, physicalTimeStep / statisticalTimePeriod);
            }
            RDouble c2 = 1.0 - c1;

            RDouble **qAverageTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    //qAverageTurbulence[m][iCell] = c1 * qAverageTurbulence[m][iCell] + c2 * qTurbulence[m][iCell];
                    qAverageTurbulence[m][iCell] = c2 * qAverageTurbulence[m][iCell] + c1 * qTurbulence[m][iCell];
                }
            }
        }
    }
}

void TurbSolverUnstr::DualTimeSource(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **qPrimitive  = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
    RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
    RDouble **residualTurbulenceUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_tmp"));

    RDouble *volume   = grid->GetCellVolume();
    RDouble *volumeN1 = grid->GetCellVolume();
    RDouble *volumeN2 = grid->GetCellVolume();

    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble *primitive0 = new RDouble[nTurbulenceEquation];
    RDouble *primitive1 = new RDouble[nTurbulenceEquation];
    RDouble *primitive2 = new RDouble[nTurbulenceEquation];

    RDouble *qConservative0 = new RDouble[nTurbulenceEquation];
    RDouble *qConservative1 = new RDouble[nTurbulenceEquation];
    RDouble *qConservative2 = new RDouble[nTurbulenceEquation];

    //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            residualTurbulenceUnsteadyTemporary[m][iCell] = residualTurbulence[m][iCell];
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
    RDouble dualTimeQC1   = dualTimeCoefficient[3];
    RDouble dualTimeQC2   = dualTimeCoefficient[4];
    RDouble dualTimeQC3   = dualTimeCoefficient[5];

    using namespace IDX;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            primitive0[m] = qTurbulence[m][iCell];
            primitive1[m] = qTurbulenceUnsteadyN1[m][iCell];
            primitive2[m] = qTurbulenceUnsteadyN2[m][iCell];
        }

        RDouble coefficient = 1.0;

        if (nTurbulenceEquation >= 2)
        {
            coefficient  = qPrimitive[IR][iCell];
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qConservative0[m] = coefficient * primitive0[m];
            qConservative1[m] = coefficient * primitive1[m];
            qConservative2[m] = coefficient * primitive2[m];
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            RDouble dualSrcRes = dualTimeResC1 * residualTurbulence[m][iCell] +
                                 dualTimeResC2 * residualTurbulenceUnsteadyN1[m][iCell] +
                                 dualTimeResC3 * residualTurbulenceUnsteadyN2[m][iCell];

            RDouble dualSrcQ   = dualTimeQC1 * qConservative0[m]* volume  [iCell] + 
                                 dualTimeQC2 * qConservative1[m]* volumeN1[iCell] + 
                                 dualTimeQC3 * qConservative2[m]* volumeN2[iCell];

            residualTurbulence[m][iCell] = dualSrcRes + dualSrcQ;
        }
    }

    delete [] primitive0;    primitive0 = nullptr;
    delete [] primitive1;    primitive1 = nullptr;
    delete [] primitive2;    primitive2 = nullptr;
    delete [] qConservative0;    qConservative0 = nullptr;
    delete [] qConservative1;    qConservative1 = nullptr;
    delete [] qConservative2;    qConservative2 = nullptr;
}

void TurbSolverUnstr::FillField(Grid *gridIn, FieldProxy *fieldProxy, RDouble value)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **field = fieldProxy->GetField_UNS();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            field[m][iCell] = value;
        }
    }
}

void TurbSolverUnstr::FillField(Grid *gridIn, FieldProxy *field1Proxy, FieldProxy *field2Proxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **field1 = field1Proxy->GetField_UNS();
    RDouble **field2 = field2Proxy->GetField_UNS();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            field1[m][iCell] = field2[m][iCell];
        }
    }
}

FieldProxy * TurbSolverUnstr::CreateFieldProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **field = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

    FieldProxy *fieldProxy = new FieldProxy();

    fieldProxy->SetField_UNS(field, true);

    return fieldProxy;
}

FieldProxy * TurbSolverUnstr::GetFieldProxy(Grid *gridIn, const string &fieldName)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **field = reinterpret_cast<RDouble **> (grid->GetDataPtr(fieldName));

    FieldProxy *fieldProxy = new FieldProxy();

    fieldProxy->SetField_UNS(field);

    return fieldProxy;
}

FieldProxy * TurbSolverUnstr::GetResidualProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    FieldProxy *residualTurbulenceProxy = new FieldProxy();

    residualTurbulenceProxy->SetField_UNS(residualTurbulence);

    return residualTurbulenceProxy;
}

void TurbSolverUnstr::RecoverResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    if (grid->GetLevel() != 0) 
    {
        RDouble **rhs = rhsProxy->GetField_UNS();
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
            {
                residualTurbulence[m][iCell] = rhs[m][iCell];
            }
        }
    }
}

void TurbSolverUnstr::StoreRhsByResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **rhs = rhsProxy->GetField_UNS();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            rhs[m][iCell] = residualTurbulence[m][iCell];
        }
    }
}

void TurbSolverUnstr::InitResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **rhs = rhsProxy->GetField_UNS();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTurbulence[m][iCell] = - rhs[m][iCell];
        }
    }
}

void TurbSolverUnstr::RestrictAllQ(Grid *fineGridIn, Grid *coarseGridIn)
{
    UnstructGrid *fineGrid   = UnstructGridCast(fineGridIn);
    UnstructGrid *coarseGrid = UnstructGridCast(coarseGridIn);

    RDouble **qTurbulenceFineGrid   = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        RestrictQ(fineGrid, qTurbulenceFineGrid[m], coarseGrid, qTurbulenceCoarseGrid[m]);
    }

    RDouble *viscousOfFineGrid   = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("vist"));
    RDouble *viscousOfCoarseGrid = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("vist"));
    RestrictQ(fineGrid, viscousOfFineGrid, coarseGrid, viscousOfCoarseGrid);
    //coarseGrid->SetGhostCellExceptInterface(cvist);
    coarseGrid->SetGhostCell(viscousOfCoarseGrid);
}

// Bell 20130513 mod
void TurbSolverUnstr::LoadFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **flux = faceProxy->GetFlux();

    //! Determine if there are boundary faces.
    int nMid  = localStart;
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

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    using namespace IDX;

    //! For boundary faces, remember re is ghost cell
    for (int i = localStart; i < nMid; ++ i)
    {
        int le,j;
        le = leftCellOfFace[i];
        j = i - localStart;
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            residualTurbulence[m][le] -= flux[m][j];
        }
    }

    //! Interior faces
    for (int i = nMid; i <localEnd; ++ i)
    {
        int le,re,j;
        le = leftCellOfFace[i];
        re = rightCellOfFace[i];
        j  = i - localStart;

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            residualTurbulence[m][le] -= flux[m][j];
            residualTurbulence[m][re] += flux[m][j];
        }
    }
}

void TurbSolverUnstr::ModifyResidualOnWallBoundary(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int modifyTurbulenceResidual = GlobalDataBase::GetIntParaFromDB("mod_turb_res");

    //! modify the residuals for the cells next to the wall
    if (modifyTurbulenceResidual != 1) return;

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    using namespace PHENGLEI;
    if (viscousType == ONE_EQU)
    {
        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegion = unstructBCSet->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            int bcType = bcRegion->GetBCType();

            if (IsWall(bcType))
            {
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter ++)
                {
                    //! iFace is the face number in the set of faceIndex.
                    int iFace = *iter;
                    int le = leftCellOfFace[iFace];
                    int re = rightCellOfFace[iFace];
                    residualTurbulence[0][le] = 0.0;
                    residualTurbulence[0][re] = 0.0;    //! Bell 20130508 add
                }
            }
        }
    }
}

void TurbSolverUnstr::ResetWallScalar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    using namespace PHENGLEI;
    if (viscousType == ONE_EQU)
    {
        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegion = unstructBCSet->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            int bcType = bcRegion->GetBCType();

            if (IsWall(bcType))
            {
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    //! iFace is the face number in the set of faceIndex.
                    int iFace = *iter;
                    int le = leftCellOfFace[iFace];
                    qTurbulence[0][le] = 0.0;
                }
            }
        }
    }
}

void TurbSolverUnstr::SourceFluxOneEquation(Grid *grid)
{
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");

    ComputeLengthScaleofOneEquationModel(grid);

    if (TurbEquationFormation() == FORMATION_SA_CFL3D)
    {
        SourceFluxOneEquationOriginalCFL3D(grid);
    }
    else
    {
        if (tscheme != GMRES)
        {
            SourceFluxOneEquationOriginal(grid);
        }
#ifdef USE_GMRESSOLVER
        else
        {
            GMRES_SourceFluxOneEquationOriginal(grid);
        }
#endif
    }
}

#ifdef USE_GMRESSOLVER
//! GMRESturb
void TurbSolverUnstr::GMRES_SourceFluxOneEquationOriginal(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **gradientPrimtiveVariableX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradientPrimtiveVariableY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradientPrimtiveVariableZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble *dudx = gradientPrimtiveVariableX[1];
    RDouble *dudy = gradientPrimtiveVariableY[1];
    RDouble *dudz = gradientPrimtiveVariableZ[1];

    RDouble *dvdx = gradientPrimtiveVariableX[2];
    RDouble *dvdy = gradientPrimtiveVariableY[2];
    RDouble *dvdz = gradientPrimtiveVariableZ[2];

    RDouble *dwdx = gradientPrimtiveVariableX[3];
    RDouble *dwdy = gradientPrimtiveVariableY[3];
    RDouble *dwdz = gradientPrimtiveVariableZ[3];

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceZ"));

    RDouble **q         = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));

    RDouble **dRdq_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dRdq_turb"));
    RDouble **dDdP_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP_turb"));
    RDouble **dRdqCoupledTerm_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dRdqCoupledTerm_turb"));
    RDouble **dDdP = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP"));
    vector<int> AI = grid->GetJacobianAI4GMRES();
    vector<int> AJ = grid->GetJacobianAJ4GMRES();

    RDouble *lengthScale = this->GetLengthScale(gridIn);
    RDouble *nxs,*nys,*nzs,*ns,*vol;

    nxs = grid->GetFaceNormalX();
    nys = grid->GetFaceNormalY();
    nzs = grid->GetFaceNormalZ();
    ns  = grid->GetFaceArea();
    vol = grid->GetCellVolume();

    int **cell2face = grid->GetCell2Face();
    vector<int> * c2c = grid->GetCell2Cell();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    using namespace IDX;

    //! Work    = temporary storage array, contains on input the
    //!          velocity gradients:
    //!          work(n,1) = ui,i
    //!          work(n,2) = u1,1**2 + u2,2**2 + u3,3**2
    //!          work(n,3) = (u1,2+u2,1)**2 + (u1,3+u3,1)**2 + (u2,3+u3,2)**2
    //!          work(n,4) = sqrt((u3,2-u2,3)**2 + (u1,3-u3,1)**2
    //!                                           + (u2,1-u1,2)**2)
    //!  calculation of dnu/dxi
    //!  dnu/dx in work(.,5), dnu/dy in work(.,6), dnu/dz in work(.,7)

    //! Von karman constant
    const RDouble karm   = 0.41;
    const RDouble cb1    = 0.1355;
    const RDouble cb2    = 0.622;
    const RDouble cv1    = 7.1;
    const RDouble cw2    = 0.3;
    const RDouble cw3    = 2.0;
    const RDouble sigma  = 2.0/3.0;
    const RDouble sdilim = 1.0e9;

    RDouble kap2,oKap2;
    kap2  = karm * karm;
    oKap2 = one / kap2;

    RDouble osigma = 1.0 / sigma;
    RDouble cw1    = cb1 * oKap2 + (one + cb2) / sigma;
    RDouble cv13   = cv1 * cv1 * cv1;
    RDouble cw36   = pow(cw3,6);
    RDouble or6    = one / six;

    RDouble cb2s, cw1Kap2;
    RDouble distance2, oDistance2;
    RDouble volume, ftrans;
    RDouble saC2, saC3;

    ADReal nueSA, gradx, grady, gradz, cmachfix, xsi, xsi2, xsi3, fv1, fv2, ft2, nuetRs;
    ADReal r, r5, g, g6, fw, grd2, rs, std, ostd, sBar, omega, olam;
    ADReal production, diffsource, wallDest, source;
    ADReal prod, dest, prodp, destp, prde, prdep, part1, part2;
    ADReal dsdnu, dfv2dk, dfwdg, dgdr, drdnu, dfwdnu, temp;

    ADReal dudy_ad, dudz_ad;
    ADReal dvdx_ad, dvdz_ad;
    ADReal dwdx_ad, dwdy_ad;
    ADReal rho, pressure, mul;

    RDouble coefficientofStateEquation = GAS_SPACE::gas->GetCoefficientOfStateEquation();
    RDouble nonDimensionalSutherlandTemperature = GlobalDataBase::GetDoubleParaFromDB("tsuth");
    RDouble viscousLaminarMin;
    GlobalDataBase::GetData("visl_min", &viscousLaminarMin, PHDOUBLE, 1);

    int nLaminar = grid->GetBlockSizeOfJacobianMatrix();
    nLaminar--;

    //! GMRESCoupled
    RDouble *gradperturb = new RDouble[nLaminar];
    for (int m = 0; m < nLaminar; m++)
    {
        gradperturb[m] = 0.0;
    }

    int colidx, colidx_bc;

    saC2 = 0.7;
    saC3 = 0.9;

    cb2s = osigma * cb2;
    cw1Kap2 = cw1 * kap2;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        //! AD info
        nueSA    = qTurbulence[ISA][iCell];
        dudy_ad  = dudy[iCell];
        dudz_ad  = dudz[iCell];
        dvdx_ad  = dvdx[iCell];
        dvdz_ad  = dvdz[iCell];
        dwdx_ad  = dwdx[iCell];
        dwdy_ad  = dwdy[iCell];
        rho      = q[IR][iCell];
        pressure = q[IP][iCell];

        nueSA.diff(0, 9);
        dudy_ad.diff(1, 9);
        dudz_ad.diff(2, 9);
        dvdx_ad.diff(3, 9);
        dvdz_ad.diff(4, 9);
        dwdx_ad.diff(5, 9);
        dwdy_ad.diff(6, 9);
        rho.diff(7, 9);
        pressure.diff(8, 9);

        //! compute viscousLaminar (mul)
        ADReal t = pressure / (coefficientofStateEquation * rho);
        mul = t * sqrt(t) * (1.0 + nonDimensionalSutherlandTemperature) / (t + nonDimensionalSutherlandTemperature);
        mul = max(viscousLaminarMin, mul);

        //! Omega: magnitude of the vorticity.
        // omega = sqrt(SQR(dwdy[iCell] - dvdz[iCell]) +
        //             SQR(dudz[iCell] - dwdx[iCell]) +
        //             SQR(dvdx[iCell] - dudy[iCell]));
        omega = pow(dwdy_ad - dvdz_ad, 2) + pow(dudz_ad - dwdx_ad, 2) + pow(dvdx_ad - dudy_ad,2);
        omega = pow(omega,0.5);
        

        volume = vol[iCell];
        cmachfix = nueSA * (dudx[iCell] + dvdy[iCell] + dwdz[iCell]);
        cmachfix = 0.0;
        olam     = rho / (mul + SMALL);
        distance2 = lengthScale[iCell] * lengthScale[iCell];
        oDistance2 = one / distance2;

        xsi      = nueSA * olam + SMALL;
        xsi3     = xsi * xsi * xsi;
        fv1      = xsi3 / (xsi3 + cv13);
        fv2      = one - xsi / (one + xsi * fv1);

        //! Here, absorb oRefReNumber into rs.
        rs     = nueSA * oKap2 * oDistance2 * oRefReNumber;
        nuetRs = nueSA * rs;

        //! sBar = fv2 * nueSA / (k2 * distance2).
        sBar    = rs * fv2;

        if (sBar >= - saC2 * omega)
        {
            std = omega + sBar;
        }
        else
        {
            std = omega + omega * (saC2 * saC2 * omega + saC3 * sBar) / ((saC3 - 2 *saC2) * omega - sBar);
        }

        ostd    = one / (std + SMALL);
        r       = rs * ostd;

        r    = min(r, static_cast<RDouble>(ten));
        r5   = pow(r, 5);
        g    = r + cw2 * (r5 * r - r);
        g6   = pow(g, 6);
        fw   = g * pow((static_cast<RDouble>(one) + cw36) / (g6 + cw36), or6);

        //! Calculate "negative" production term in order to assure a zero.
        //! Solution in the laminar region.

        ftrans = zero;
        xsi2   = xsi * xsi;
        ft2    = 0.0;       //! Do Not consider transition here!
        grd2 = SQR(dqdx[ISA][iCell]) + SQR(dqdy[ISA][iCell]) + SQR(dqdz[ISA][iCell]);

        //! cb1  * omega * nueSA + cb1  * nueSA_rs * fv2.
        production  = cb1 * (one - ft2) * std * nueSA;       

        //! cb2s * rho * gradnue2.
        diffsource  = cb2s * grd2 * oRefReNumber;           

        //! ft2 = 0.0, cw1Kap2 = Cw1 * Kap2, so,
        //! cw1Kap2 * fw * nuetRs = (Cw1 * Kap2 * fw * nueSA) * (nueSA2 * oKap2 * oDistance2)
        //!                       = Cw1 * fw * (nueSA/d)^2
        wallDest   = (cw1Kap2 * fw - cb1 * ft2) * nuetRs; //! The Re has been included in nuetRs.

        diffsource = min(sdilim * production, diffsource);
        source     = production - wallDest;

        if(TurbEquationFormation() == FORMATION_SA_NSMB)
        {
            //source  += diffsource;                        //! Bell 20130513 delete for modsa.
        }
        else if(TurbEquationFormation() == FORMATION_ORGSA_SST)
        {
            source  += diffsource;
        }

        source *= volume;

        residualTurbulence[ISA][iCell] += source.val();

        // assemable Jacobian matrix, this is 1eq-sa, so only scalar production
        colidx = GetIndexOfBlockJacobianMatrix(AI,AJ,1,iCell,iCell);
        dRdq_turb[ISA][colidx] -= source.dx(0);

        //! GMRESCoupled
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, iCell, iCell);
        dRdqCoupledTerm_turb[ISA][colidx + IR] -= source.dx(7);
        dRdqCoupledTerm_turb[ISA][colidx + IP] -= source.dx(8);

        // find neighbors
        for (int j = 0; j < faceNumberOfEachCell[iCell]; j++)
        {
            int face = cell2face[iCell][j];
            int le = leftCellofFace[face];
            int re = rightCellofFace[face];
            gradperturb[IU] = (nys[face] * source.dx(1) + nzs[face] * source.dx(2)) / volume * 0.5 * ns[face];
            gradperturb[IV] = (nxs[face] * source.dx(3) + nzs[face] * source.dx(4)) / volume * 0.5 * ns[face];
            gradperturb[IW] = (nxs[face] * source.dx(5) + nys[face] * source.dx(6)) / volume * 0.5 * ns[face];
            //! if neighbor cell sits at the left cell of face
            if(iCell == re)
            {
                // perturbation of the left cell, influences the right cell
                // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, 1, iCell, le);
                // dRdq_turb[ISA][colidx] -= (-nxs[face] * source.dx(1) - nys[face] * source.dx(2) - nzs[face] * source.dx(3)) / volume * 0.5 * ns[face];
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, iCell, le);
                dRdqCoupledTerm_turb[ISA][colidx + IU] += gradperturb[IU];
                dRdqCoupledTerm_turb[ISA][colidx + IV] += gradperturb[IV];
                dRdqCoupledTerm_turb[ISA][colidx + IW] += gradperturb[IW];
            }
            //! if neighbor cell sits at the right cell of face
            else
            {
                //! not boundary face
                if(re < numberOfTotalCell)
                {
                    // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, 1, iCell, re);
                    // dRdq_turb[ISA][colidx] -= (nxs[face] * source.dx(1) + nys[face] * source.dx(2) + nzs[face] * source.dx(3)) / volume * 0.5 *  ns[face];
                    colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, iCell, re);
                    dRdqCoupledTerm_turb[ISA][colidx + IU] -= gradperturb[IU];
                    dRdqCoupledTerm_turb[ISA][colidx + IV] -= gradperturb[IV];
                    dRdqCoupledTerm_turb[ISA][colidx + IW] -= gradperturb[IW];
                }
                else
                {
                    if (interfaceInfo != nullptr && interfaceInfo->MatchedLocalPhysicalCellIndex(iCell) >= -1)
                    {
                        // it is for the interface boundary
                        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, re, iCell);
                        dRdqCoupledTerm_turb[ISA][colidx + IU] += gradperturb[IU];
                        dRdqCoupledTerm_turb[ISA][colidx + IV] += gradperturb[IV];
                        dRdqCoupledTerm_turb[ISA][colidx + IW] += gradperturb[IW];
                    }
                    else
                    {
                        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, iCell, iCell);
                        colidx_bc = (re - numberOfTotalCell) * nLaminar;
                        for (int indexI = 0; indexI < nLaminar; indexI++)
                        {
                            for (int indexJ = 0; indexJ < nLaminar; indexJ++)
                            {
                                dRdqCoupledTerm_turb[ISA][colidx + indexI] -= gradperturb[indexJ] * dDdP[indexJ][colidx_bc + indexI];
                            }
                        }
                    }
                }
            }
        }

        //! Spectrum.
        dfv2dk = (three * cv13 * pow(fv1/(xsi+SMALL), 2) - one) / pow(one + xsi * fv1, 2);
        dsdnu  = oRefReNumber * oKap2 * oDistance2 * (fv2 + xsi * dfv2dk);
        dfwdg   = fw / (g + SMALL) * (one - g6 / (g6 + cw36));
        dgdr    = one + cw2 * (six * r5 - one);
        drdnu   = oRefReNumber * oKap2 * oDistance2 * ostd * (one - nueSA * ostd * dsdnu);
        dfwdnu  = dfwdg * dgdr * drdnu;

        prod = cb1 * (one - ft2) * std;
        dest = (cw1*fw - cb1*ft2*oKap2) * nueSA * oDistance2 * oRefReNumber;

        prodp = cb1 * ((one - ft2) * dsdnu);
        destp = (cw1 * (fw + dfwdnu * nueSA) - cb1 * ft2 * oKap2) * oDistance2 * oRefReNumber;

        prde  = prod  - dest;
        prdep = prodp - destp;

        part1 = - half * (prde  - ABS(prde));
        part2 = - half * (prdep - ABS(prdep)) * nueSA;

        temp = (part1 + part2) * volume; 
        spectrumTurbulence[ISA][iCell] += temp.val(); 
    }

    delete []gradperturb;
}
#endif

//! Bell 20130513 mod
void TurbSolverUnstr::SourceFluxOneEquationOriginal(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble **gradientPrimtiveVariableX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradientPrimtiveVariableY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradientPrimtiveVariableZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble *dudx = gradientPrimtiveVariableX[IDX::IU];
    RDouble *dudy = gradientPrimtiveVariableY[IDX::IU];
    RDouble *dudz = gradientPrimtiveVariableZ[IDX::IU];

    RDouble *dvdx = gradientPrimtiveVariableX[IDX::IV];
    RDouble *dvdy = gradientPrimtiveVariableY[IDX::IV];
    RDouble *dvdz = gradientPrimtiveVariableZ[IDX::IV];

    RDouble *dwdx = gradientPrimtiveVariableX[IDX::IW];
    RDouble *dwdy = gradientPrimtiveVariableY[IDX::IW];
    RDouble *dwdz = gradientPrimtiveVariableZ[IDX::IW];

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceZ"));

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));

    RDouble *lengthScale = this->GetLengthScale(gridIn);
    RDouble *nxs, *nys, *nzs, *ns, *vol;

    nxs = grid->GetFaceNormalX();
    nys = grid->GetFaceNormalY();
    nzs = grid->GetFaceNormalZ();
    ns  = grid->GetFaceArea();
    vol = grid->GetCellVolume();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    using namespace IDX;

    //! Work   = temporary storage array, contains on input the
    //!          velocity gradients:
    //!          work(n,1) = ui,i
    //!          work(n,2) = u1,1**2 + u2,2**2 + u3,3**2
    //!          work(n,3) = (u1,2+u2,1)**2 + (u1,3+u3,1)**2 + (u2,3+u3,2)**2
    //!          work(n,4) = sqrt((u3,2-u2,3)**2 + (u1,3-u3,1)**2
    //!                                           + (u2,1-u1,2)**2)
    //!  calculation of dnu/dxi
    //!  dnu/dx in work(.,5), dnu/dy in work(.,6), dnu/dz in work(.,7)

    //! Von karman constant
    const RDouble karm   = 0.41;
    const RDouble cb1    = 0.1355;
    const RDouble cb2    = 0.622;
    const RDouble cv1    = 7.1;
    const RDouble cw2    = 0.3;
    const RDouble cw3    = 2.0;
    const RDouble sigma  = 2.0/3.0;
    const RDouble sdilim = 1.0e9;

    RDouble kap2,oKap2;
    kap2  = karm * karm;
    oKap2 = one / kap2;

    RDouble osigma = 1.0 / sigma;
    RDouble cw1    = cb1 * oKap2 + (one + cb2) / sigma;
    RDouble cv13   = cv1 * cv1 * cv1;
    RDouble cw36   = pow(cw3,6);
    RDouble or6    = one / six;
    RDouble ct3    = 1.2;

    RDouble production, diffsource, wallDest, source;
    RDouble cb2s, cw1Kap2, rs, std, ostd, nuetRs;
    RDouble prod, dest, prodp, destp, prde, prdep, part1, part2;
    RDouble distance2, oDistance2, xsi, xsi2, xsi3, fv1, fv2, ft2;
    RDouble r, r5, g, g6, fw, grd2;
    RDouble dsdnu, dfv2dk, dfwdg, dgdr, drdnu, dfwdnu;
    RDouble nueSA, volume, ftrans, olam, cmachfix;
    RDouble sBar, omega;
    RDouble saC2, saC3;

    saC2 = 0.7;
    saC3 = 0.9;

    cb2s = osigma * cb2;
    cw1Kap2 = cw1 * kap2;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        nueSA = qTurbulence[ISA][iCell];

        //! Omega: magnitude of the vorticity.
        omega = sqrt(SQR(dwdy[iCell] - dvdz[iCell]) +
                    SQR(dudz[iCell] - dwdx[iCell]) +
                    SQR(dvdx[iCell] - dudy[iCell]));

        volume   = vol[iCell];
        cmachfix = nueSA * (dudx[iCell] + dvdy[iCell] + dwdz[iCell]);
        cmachfix = 0.0;
        olam     = q[IR][iCell] / (viscousLaminar[iCell] + SMALL);
        distance2= lengthScale[iCell] * lengthScale[iCell];
        oDistance2 = one / distance2;

        xsi  = nueSA * olam + SMALL;
        xsi3 = xsi * xsi * xsi;
        fv1  = xsi3 / (xsi3 + cv13);
        fv2  = one - xsi / (one + xsi * fv1);

        //! Here, absorb oRefReNumber into rs.
        rs     = nueSA * oKap2 * oDistance2 * oRefReNumber;
        nuetRs = nueSA * rs;

        //! sBar = fv2 * nueSA / (k2 * distance2).
        sBar    = rs * fv2;

        if (sBar >= - saC2 * omega)
        {
            std = omega + sBar;
        }
        else
        {
            std = omega + omega * (saC2 * saC2 * omega + saC3 * sBar) / ((saC3 - 2 *saC2) * omega - sBar);
        }

        //std     = MAX(str + rs * fv2, 0.3*str);
        ostd = one / (std + SMALL);
        r    = rs * ostd;

        r  = MIN(r, static_cast<RDouble>(ten));
        r5 = pow(r, 5);
        g  = r + cw2 * (r5 * r - r);
        g6 = pow(g, 6);
        fw = g * pow((static_cast<RDouble>(one) + cw36) / (g6 + cw36), or6);

        //! Calculate "negative" production term in order to assure a zero.
        //! Solution in the laminar region.

        ftrans = zero;
        xsi2   = xsi * xsi;
        ft2    = 0.0;       //! Do Not consider transition here!
        //ft2    = ct3 * exp(- ct4 * xsi2) * ftrans;
        grd2 = SQR(dqdx[ISA][iCell]) + SQR(dqdy[ISA][iCell]) + SQR(dqdz[ISA][iCell]);

        //! cb1  * omega * nueSA + cb1  * nueSA_rs * fv2.
        production = cb1 * (one - ft2) * std * nueSA;

        //! cb2s * rho * gradnue2.
        diffsource = cb2s * grd2 * oRefReNumber;

        //! ft2 = 0.0, cw1Kap2 = Cw1 * Kap2, so,
        //! cw1Kap2 * fw * nuetRs = (Cw1 * Kap2 * fw * nueSA) * (nueSA2 * oKap2 * oDistance2)
        //!                       = Cw1 * fw * (nueSA/d)^2
        wallDest = (cw1Kap2 * fw - cb1 * ft2) * nuetRs;    //! The Re has been included in nuetRs.

        RDouble destrAlpha = 50.0;
        if(nueSA < 0.0)
        {
            production = cb1 * (one - ct3) * omega * nueSA;
            wallDest   = - destrAlpha * cw1Kap2 * nuetRs;
        }

        diffsource = MIN(sdilim * production, diffsource);
        source     = production - wallDest;

        if(TurbEquationFormation() == FORMATION_SA_NSMB)
        {
            //source  += diffsource;    //! Bell 20130513 delete for modsa.
        }
        else if(TurbEquationFormation() == FORMATION_ORGSA_SST)
        {
            source += diffsource;
        }

        residualTurbulence[ISA][iCell] += (source) * volume;

        //! Spectrum.
        dfv2dk = (three * cv13 * pow(fv1/(xsi+SMALL), 2) - one) / pow(one + xsi * fv1, 2);
        dsdnu  = oRefReNumber * oKap2 * oDistance2 * (fv2 + xsi * dfv2dk);
        dfwdg  = fw / (g + SMALL) * (one - g6 / (g6 + cw36));
        dgdr   = one + cw2 * (six * r5 - one);
        drdnu  = oRefReNumber * oKap2 * oDistance2 * ostd * (one - nueSA * ostd * dsdnu);
        dfwdnu = dfwdg * dgdr * drdnu;

        //dft2dnu = - two * ct3 * ct4 * xsi * exp (- ct4 * xsi2) * olam * ftrans;

        prod = cb1 * (one - ft2) * std;
        dest = (cw1*fw - cb1*ft2*oKap2) * nueSA * oDistance2 * oRefReNumber;

        prodp = cb1 * ((one - ft2) * dsdnu);
        destp = (cw1 * (fw + dfwdnu * nueSA) - cb1 * ft2 * oKap2) * oDistance2 * oRefReNumber;

        if(nueSA < 0.0)
        {
            prod = cb1 * (one - ct3) * omega;
            dest = - destrAlpha * cw1 * nueSA * oDistance2 * oRefReNumber;

            prodp = cb1 * ((one - ct3) * dsdnu);
            destp = -destrAlpha * cw1 * oDistance2 * oRefReNumber;
        }

        prde  = prod  - dest;
        prdep = prodp - destp;

        part1 = - half * (prde  - ABS(prde));
        part2 = - half * (prdep - ABS(prdep)) * nueSA;

        spectrumTurbulence[ISA][iCell] += (part1 + part2) * volume;
    }
}

void TurbSolverUnstr::SourceFluxOneEquationOriginalCFL3D(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **gradientPrimtiveVariableX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradientPrimtiveVariableY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradientPrimtiveVariableZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble *drhodx = gradientPrimtiveVariableX[IDX::IR];
    RDouble *drhody = gradientPrimtiveVariableY[IDX::IR];
    RDouble *drhodz = gradientPrimtiveVariableZ[IDX::IR];

    RDouble *dudy = gradientPrimtiveVariableY[IDX::IU];
    RDouble *dudz = gradientPrimtiveVariableZ[IDX::IU];

    RDouble *dvdx = gradientPrimtiveVariableX[IDX::IV];
    RDouble *dvdz = gradientPrimtiveVariableZ[IDX::IV];

    RDouble *dwdx = gradientPrimtiveVariableX[IDX::IW];
    RDouble *dwdy = gradientPrimtiveVariableY[IDX::IW];

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceZ"));

    RDouble **q         = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));

    RDouble *lengthScale = this->GetLengthScale(gridIn);
    RDouble *nxs,*nys,*nzs,*ns,*vol;

    nxs = grid->GetFaceNormalX();
    nys = grid->GetFaceNormalY();
    nzs = grid->GetFaceNormalZ();
    ns  = grid->GetFaceArea();
    vol = grid->GetCellVolume();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    using namespace IDX;

    //! Von karman constant
    const RDouble karm   = 0.41;
    const RDouble cb1    = 0.1355;
    const RDouble cb2    = 0.622;
    const RDouble cv1    = 7.1;
    const RDouble cw2    = 0.3;
    const RDouble cw3    = 2.0;
    const RDouble sigma  = 2.0/3.0;

    RDouble kap2,oKap2;
    kap2  = karm * karm;
    oKap2 = one / kap2;

    RDouble osigma = 1.0 / sigma;
    RDouble cw1    = cb1 * oKap2 + (one + cb2) / sigma;
    RDouble cv13   = cv1 * cv1 * cv1;
    RDouble cw36   = pow(cw3, 6.0);
    RDouble or6    = one / six;

    RDouble ct3 = 1.2;
    RDouble ct4 = 0.5; 

    RDouble production, wallDest, source;
    RDouble cb2s, cw1Kap2, S, oS;
    RDouble distance2, oDistance2, xsi, xsi3, fv1, fv2, ft2;
    RDouble r, r5, g, g6, fw;
    RDouble nue, nueSA, volume, ftrans, rho, vis;
    RDouble sBar, omega;
    RDouble saC2, saC3;

    saC2 = 0.7;
    saC3 = 0.9;

    cb2s = osigma * cb2;
    cw1Kap2 = cw1 * kap2;

    RDouble machinZero = pow(10.0, -machZeroExp + 1);

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        nueSA = qTurbulence[ISA][iCell];
        rho   = q[IR][iCell];
        vis   = viscousLaminar[iCell];
        nue   = vis / (rho + SMALL);        

        //! Omega: magnitude of the vorticity.
        omega = sqrt(SQR(dwdy[iCell] - dvdz[iCell]) +
                     SQR(dudz[iCell] - dwdx[iCell]) +
                     SQR(dvdx[iCell] - dudy[iCell]));   

        volume   = vol[iCell];
  
        distance2= lengthScale[iCell] * lengthScale[iCell];
        oDistance2 = one / distance2;

        xsi      = nueSA / nue;
        xsi3     = xsi * xsi * xsi;
        fv1      = xsi3 / (xsi3 + cv13);
        fv2      = one - xsi / (one + xsi * fv1);

        //! Important, Reynolds will exist when nueSA / ( k2 * distance2) appear,
        //! so include the oRefReNumber for coding convenient.
        //! Here, absorb oRefReNumber into rs(term0): nueSA / (k2 * distance2).
        //! oRefReNumber is used both in S and r, so it is included into term0.
        RDouble term0 = nueSA * oKap2 * oDistance2 * oRefReNumber;        

        //!sBar = fv2 * nueSA / (k2 * distance2).
        sBar    = term0 * fv2;

        if (sBar >= - saC2 * omega)
        {
            S = omega + sBar;
        }
        else
        {
            S = omega + omega * (saC2 * saC2 * omega + saC3 * sBar) / ((saC3 - 2 *saC2) * omega - sBar);
        }
        oS    = one / (S + SMALL);

        r     = term0 * oS;
        r    = MIN(r, static_cast<RDouble>(ten));

        r5   = pow(r, 5);
        g    = r + cw2 * (r5 * r - r);
        g    = MAX(g, machinZero);

        g6   = pow(g, 6);
        fw   = g * pow((static_cast<RDouble>(one) + cw36) / (g6 + cw36), or6);

        //! ftrans: transition or not -- 1/0.
        //! More work should be done to compare the influence of ft2.
        ftrans = zero;  // one.
        ft2    = ct3 * exp(- ct4 * xsi * xsi);
        ft2   *= ftrans;

        //! Using CFL3D terms.
        //! CFL3D user's manual v5.0, and, open source code.
        //! Here, consider omega only, since the sBar has been moved into 
        //! the wallDest term because they all have the Reynolds.
        RDouble term1_CFL3D = cb1 * (one - ft2) * omega;
        RDouble term2_CFL3D = cb1 * ((one - ft2) * fv2 + ft2) * oKap2 - cw1 * fw;

        RDouble destrAlpha = 50.0;
        if(nueSA < 0.0)
        {
            term1_CFL3D = cb1 * ( one - ct3) * omega;
            term2_CFL3D = destrAlpha * cw1;
        }

        production  = rho * term1_CFL3D * nueSA;     
        wallDest    = rho * term2_CFL3D * nueSA * nueSA * oDistance2 * oRefReNumber;
        source      = production + wallDest;

        //! Conservative term: (1+cb2) * Grad(Rho) * Grad(nue) / sigma.
        RDouble gradRhoNue = drhodx[iCell] * dqdx[ISA][iCell] + 
                             drhody[iCell] * dqdy[ISA][iCell] +
                             drhodz[iCell] * dqdz[ISA][iCell];
        RDouble nuePositive  = PositiveSA(nue, vis, rho);
        RDouble conserveTerm = (1 + cb2) * osigma * (vis / rho + nuePositive) * gradRhoNue * oRefReNumber;
        source -= conserveTerm;

        residualTurbulence[ISA][iCell] += source * volume;
        
        //! Spectrum, reference from CFL3D codes.
        RDouble derivativePart1 = 2.0 * term2_CFL3D * nueSA * oDistance2;
        RDouble dfv1  = (fv1 - fv1 * fv1) * 3.0 / nueSA;
        RDouble dfv2  = (fv2 - 1.0) / nueSA + (1.0 - fv2) * (1.0 - fv2) * (fv1 / nueSA + dfv1);
        RDouble dft2  = -(2.0 * ct4 * nueSA / (nue * nue)) * ft2;
        RDouble drr   = r / nueSA - r * r * (fv2 / nueSA + dfv2); 
        RDouble dgg   = (1.0 - cw2 + 6.0 * cw2 * pow(r, 5.0)) * drr;

        g  = MAX(g, 10.0 * machinZero); //! Fix for single prec. via CFL3D via TLNS3D via FUN3D vi Tim Barth
        g6 = pow(g, 6);

        RDouble dfw   = pow((1.0 + cw36) / (g6 + cw36), or6) - 
                       (pow(1.0 + cw36, or6) / (pow(g6 + cw36,(7.0 / 6.0)))) * g6;      
        dfw *= dgg;  
        RDouble derivativePart2  = oDistance2 * nueSA * nueSA * 
            (cb1 * oKap2 * (dfv2 - ft2 * dfv2 - fv2 * dft2 + dft2)-cw1 * dfw);

        if(nueSA < 0.0)
        {
            derivativePart2 = 0.0;
        }

        //! The source spectrum terms lie on left side of equation, so minus it.
        spectrumTurbulence[ISA][iCell] -= (derivativePart1 + derivativePart2) * oRefReNumber * volume;

        //! The following term is negatived. it may delay the converge, although may improve robust.
        spectrumTurbulence[ISA][iCell] -= term1_CFL3D * volume;
    }
}

void TurbSolverUnstr::SourceFluxTwoEquation(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **primitiveVariables = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    RDouble *xc   = grid->GetCellCenterX();
    RDouble *vol  = grid->GetCellVolume();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble refReNumber  = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    string viscousName = parameters->GetViscousName();

    int neasm = -1;
    GlobalDataBase::GetData("neasm",&neasm, PHINT, 1);

    int SSTProductType = parameters->GetSSTProductType();

    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble cblend, density, ke, kw;
    RDouble sij2, divv, prodk, dissk, prodw, dissw, cdkww, crss, mut;
    RDouble srck, srcw, fskn, fswn, diak, diaw, oork, oorw, volume;
    /*RDouble turbXk, turbXw, turbXk2, turbFk2, turbFw;*/
    RDouble s11, s22, s33, s12, s13, s23, w12, w13, w23;
    RDouble vort2;

    RDouble turbFbetas, turbFbeta, turbBeta, turbAlphaw/*, turbSigd*/;

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

    RDouble **aniss = 0;
    RDouble *sengy  = 0;

    if (neasm > 0)
    {
        aniss = reinterpret_cast<RDouble **> (grid->GetDataPtr("aniss"));
        sengy = reinterpret_cast<RDouble * > (grid->GetDataPtr("sengy"));
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

    using namespace IDX;
    RDouble refMaNumber = parameters->GetRefMachNumber();
    RDouble SpSdlimit = MIN(1.0 + 0.00075 * refMaNumber * refMaNumber * refMaNumber * refMaNumber, 1.2);

    int transitionType = parameters->GetTransitionType();

    RDouble *freeStreamTurbVar = parameters->GetFreeStreamTurbVar();

    //! K-OMEGA MODELS
    if (viscousName.substr(0, 6) == "2eq-kw")
    {
        //! MENTER'S K-OMEGA SST MODEL
        if (viscousName.substr(0, 17) == "2eq-kw-menter-sst")
        {
            //turbulence model constants
            RDouble SSTProductLimit = parameters->GetSSTProductLimit();
            GlobalDataBase::GetData("turb_fbetas", &turbFbetas, PHDOUBLE, 1);
            GlobalDataBase::GetData("turb_fbeta" , &turbFbeta , PHDOUBLE, 1);
            RDouble SST_beta1 = parameters->GetSST_beta1();
            RDouble SST_beta2 = parameters->GetSST_beta2();
            RDouble SST_betaStar = parameters->GetSST_betaStar();
            RDouble KW_sigmaW2 = parameters->GetKW_sigmaW2();
            RDouble SST_alphaw1 = parameters->GetSST_alphaw1();
            RDouble SST_alphaw2 = parameters->GetSST_alphaw2();

            RDouble *cross = reinterpret_cast<RDouble  *> (grid->GetDataPtr("cross"));
            RDouble *blend = reinterpret_cast<RDouble  *> (grid->GetDataPtr("blend"));
            RDouble *SpSdRatio = NULL, *gamaeff = NULL;
            if (transitionType == IREGAMA)
            {
                SpSdRatio = reinterpret_cast<RDouble *> (grid->GetDataPtr("SpSdRatio"));
                gamaeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("gamaeff"));
            }

            for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
            {
                cblend  = blend[iCell];
                crss    = cross[iCell];
                density = primitiveVariables[IR][iCell];
                ke      = qTurbulence[IKE][iCell];
                kw      = qTurbulence[IKW][iCell];
                mut     = viscousTurbulence[iCell];
                volume  = vol[iCell];

                turbBeta   = cblend * SST_beta1   + (1.0 - cblend) * SST_beta2;
                turbAlphaw = cblend * SST_alphaw1 + (1.0 - cblend) * SST_alphaw2;

                s11 = dudx[iCell];
                s22 = dvdy[iCell];
                s33 = dwdz[iCell];
                s12 = half * (dudy[iCell] + dvdx[iCell]);
                s13 = half * (dudz[iCell] + dwdx[iCell]);
                s23 = half * (dvdz[iCell] + dwdy[iCell]);
                w12 = half * (dudy[iCell] - dvdx[iCell]);
                w13 = half * (dudz[iCell] - dwdx[iCell]);
                w23 = half * (dvdz[iCell] - dwdy[iCell]);
                
                sij2 = two * (s11*s11 + s22*s22 + s33*s33 + two*(s12*s12 + s13*s13 + s23*s23));
                divv = s11 + s22 + s33;

                if (neasm < 0)
                {
                    if (SSTProductType == 0)
                    {
                        //! Boussinesq approximation,full production.
                        prodk = mut * (sij2 - two3rd * divv * divv) * oRefReNumber - two3rd * density * ke * divv;
                    }
                    else
                    {
                        vort2 = four * (w12 * w12 + w13 * w13 + w23 * w23);
                        prodk = mut * (vort2) * oRefReNumber - two3rd * density * ke * divv;
                    }
                }
                dissk = turbFbetas * SST_betaStar * density * ke * kw * refReNumber;

                //! Limit production
                prodk = MIN(prodk, SSTProductLimit * dissk);
                prodw = turbAlphaw * density / (mut + SMALL) * prodk;
                dissw = turbFbeta * turbBeta * density * kw * kw * refReNumber;
                cdkww = two * (one - cblend) * density * KW_sigmaW2 * crss / (kw + SMALL) * oRefReNumber;
                //! Re-gama transition model
                if (transitionType == IREGAMA && freeturbIntensitySRModify == 1 && xc[iCell] < freeDecayXLocation)
                {
                    //! S-R Modify
                    RDouble keamb = freeStreamTurbVar[IKE];
                    RDouble kwamb = freeStreamTurbVar[IKW];
                    prodk += turbFbetas * SST_betaStar * density * keamb * kwamb * refReNumber;
                    prodw += turbFbeta * turbBeta * density * kwamb * kwamb * refReNumber;
                }

                //! Re-gama transition model
                if (transitionType == IREGAMA)
                {
                    RDouble gmeff = gamaeff[iCell];
                    prodk = CorrectionOfProductionInKEquation(gmeff, prodk);
                    dissk = CorrectionOfDestructionInKEquation(gmeff, dissk);
                    SpSdRatio[iCell] = MIN(MAX(prodk / (dissk + SMALL), 1.0), SpSdlimit);
                }

                srck  = prodk - dissk;
                srcw  = prodw - dissw + cdkww;

                residualTurbulence[IKE][iCell] += (srck) * volume;
                residualTurbulence[IKW][iCell] += (srcw) * volume;

                if (neasm > 0)
                {
                    sengy[iCell] = - (srck) * volume;
                }

                diak = - two * turbFbetas * SST_betaStar * kw;
                diaw = - two * turbFbeta  * turbBeta  * kw - ABS(cdkww) / (density * kw + SMALL);

                fskn = half * (srck - ABS(srck));
                fswn = half * (srcw - ABS(srcw));

                oork = 1.0 / (density * ke + SMALL);
                oorw = 1.0 / (density * kw + SMALL);

                spectrumTurbulence[IKE][iCell] += - (diak + fskn * oork) * volume;
                spectrumTurbulence[IKW][iCell] += - (diaw + fswn * oorw) * volume;
            }
        }
    }

    delete []dudx;
    delete []dudy;
    delete []dudz;

    delete []dvdx;
    delete []dvdy;
    delete []dvdz;

    delete []dwdx;
    delete []dwdy;
    delete []dwdz;
}

void TurbSolverUnstr::Crossing(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble *dkedx = new RDouble [numberOfTotal];
    RDouble *dkedy = new RDouble [numberOfTotal];
    RDouble *dkedz = new RDouble [numberOfTotal];
    RDouble *dkwdx = new RDouble [numberOfTotal];
    RDouble *dkwdy = new RDouble [numberOfTotal];
    RDouble *dkwdz = new RDouble [numberOfTotal];

    RDouble *cross = reinterpret_cast<RDouble  *> (grid->GetDataPtr("cross"));

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    grid->CompGradient(qTurbulence[0], dkedx, dkedy, dkedz);
    grid->CompGradient(qTurbulence[1], dkwdx, dkwdy, dkwdz);

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        cross[iCell] = dkedx[iCell] * dkwdx[iCell] + dkedy[iCell] * dkwdy[iCell] + dkedz[iCell] * dkwdz[iCell];
    }
    delete []dkedx;
    delete []dkedy;
    delete []dkedz;

    delete []dkwdx;
    delete []dkwdy;
    delete []dkwdz;
}

void TurbSolverUnstr::Blending(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar   = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *cross = reinterpret_cast<RDouble  *> (grid->GetDataPtr("cross"));
    RDouble *blend = reinterpret_cast<RDouble  *> (grid->GetDataPtr("blend"));

    RDouble *wallDistance = grid->GetWallDist();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble oReynolds = 1.0 / refReNumber;
    RDouble oReynoldsSqr = oReynolds * oReynolds;
    int transitionType = parameters->GetTransitionType();

    RDouble betas = 0.09;
    RDouble sigw2 = 0.856;

    using namespace IDX;

    //! Compute maximum cross diffusion term across flowfield
    RDouble cdkwmax, crossDiff;
    cdkwmax = 0.0;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        RDouble rho  = q[IR][iCell];
        RDouble kw   = qTurbulence[IKW][iCell];
        RDouble crss = cross[iCell];

        crossDiff = 2.0 * (rho) * crss * sigw2 / (kw + SMALL);
    }

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        RDouble rho   = q[IR][iCell];
        RDouble ke    = qTurbulence[IKE][iCell];
        RDouble kw    = qTurbulence[IKW][iCell];
        RDouble crss  = cross[iCell];
        RDouble visL  = viscousLaminar[iCell];
        RDouble dist  = wallDistance[iCell];
        RDouble dist2 = dist * dist;
        RDouble cdkw, term1, term2, term3, term4, arg1;

        //! calculate cd_kw
        crossDiff = 2.0 * (rho) * sigw2 * crss / (kw + SMALL);

        //! Original Menter CD_kw calculation
        cdkw = MAX(crossDiff, 1.0e-10);

        //! calculate arg1
        term1 = sqrt(ABS(ke)) * oReynolds / (betas * dist * kw + SMALL);
        term2 = 500.0 * visL * oReynoldsSqr / (rho * kw * dist2 + SMALL);
        term3 = MAX(term1, term2);
        term4 = 4.0 * rho * sigw2 * ke / (cdkw * dist2 + SMALL);
        arg1  = MIN(term3, term4);

        //! Compute the blending function fbsl
        blend[iCell] = tanh(arg1 * arg1 * arg1 * arg1);

        if (transitionType == IREGAMA)
        {
            CorrectionOfBlendingFunctionInSST(rho, dist, visL, ke, refReNumber, blend[iCell]);
        }
    }

    grid->SetGhostCellExceptInterface(blend);
}

FaceProxy * TurbSolverUnstr::CreateFaceProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast (gridIn);

    FaceProxy *faceProxy = CreateFaceProxyTurbulence(grid);
    faceProxy->SetNext(CreateFaceProxyNS(grid));

    return faceProxy;
}

FaceProxy * TurbSolverUnstr::CreateFaceProxyNS(Grid *gridIn)
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

FaceProxy * TurbSolverUnstr::CreateFaceProxyTurbulence(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    FaceProxy *faceProxy = new FaceProxy();
    faceProxy->Create(SEGCTION_LENGTH, nTurbulenceEquation);

    return faceProxy;
}

TurbFaceValue * TurbSolverUnstr::CreateTurbFaceValue(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    TurbFaceValue *facevar = new TurbFaceValue(nTurbulenceEquation, SEGCTION_LENGTH);
    return facevar;
}

void TurbSolverUnstr::ComputeQTurbNodeValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalNode = grid->GetNTotalNode();
    int numberOfTotalFace = grid->GetNTotalFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qTurbNode"));
    RDouble     **qTurbulence = reinterpret_cast<RDouble **>(grid->GetDataPtr("q_turb"));
    int   *nodeValueSliceTurb = reinterpret_cast<int      *>(grid->GetDataPtr("nodeValueSliceTurb"));

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *face2node = grid->GetFace2Node();
    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
 
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iNode = 0; iNode < numberOfTotalNode; ++ iNode)
    {
        nodeValueSliceTurb[iNode] = 0;

        for(int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qTurbulenceNode[m][iNode] = 0; 
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
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            int le = leftCellOfFace [iFace];
            int re = rightCellOfFace[iFace];

            if (bcType != PHENGLEI::INTERFACE && bcType != PHENGLEI::OVERSET)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    //! From left
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] += qTurbulence[m][le];
                    }            
                    nodeValueSliceTurb[point] += 1;

                    //! From right
                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] += qTurbulence[m][re];
                    }
                    nodeValueSliceTurb[point] += 1;
                }
            }
            else
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];
                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] += qTurbulence[m][le];
                    }
                    nodeValueSliceTurb[point] += 1;
                }
            }
            nodePosition += nodeNumberOfEachFace[iFace];
        }
    }

    //! Interior faces
    for (int iFace = numberOfBoundaryFace; iFace < numberOfTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];
        for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
        {
            //! From left
            int point = face2node[nodePosition + jNode];
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qTurbulenceNode[m][point] += qTurbulence[m][le];
            }
            ++ nodeValueSliceTurb[point];

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qTurbulenceNode[m][point] += qTurbulence[m][re];
            }
            ++ nodeValueSliceTurb[point];
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
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType == PHENGLEI::SYMMETRY)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] = 0;
                    }            
                    nodeValueSliceTurb[point] = 0;
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
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            int le = leftCellOfFace[iFace];
            int re = iFace + numberOfTotalCell;

            if (bcType == PHENGLEI::SYMMETRY)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] += qTurbulence[m][le];
                    }            
                    ++ nodeValueSliceTurb[point];

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] += qTurbulence[m][re];
                    }
                    ++ nodeValueSliceTurb[point];
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
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType == PHENGLEI::FARFIELD)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] = 0;
                    }            
                    nodeValueSliceTurb[point] = 0;
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
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            int le = leftCellOfFace[iFace];
            int re = iFace + numberOfTotalCell;

            if (bcType == PHENGLEI::FARFIELD)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] += qTurbulence[m][le];
                    }            
                    ++ nodeValueSliceTurb[point];

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] += qTurbulence[m][re];
                    }
                    ++ nodeValueSliceTurb[point];
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
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] = 0;
                    }            
                    nodeValueSliceTurb[point] = 0;
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
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            int le = leftCellOfFace[iFace];
            int re = iFace + numberOfTotalCell;
            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2node[nodePosition + jNode];

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] += qTurbulence[m][le];
                    }            
                    ++ nodeValueSliceTurb[point];

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulenceNode[m][point] += qTurbulence[m][re];
                    }
                    ++ nodeValueSliceTurb[point];
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
                for (int m = 0; m < nTurbulenceEquation;++ m)
                {    
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    qTurbulenceNode[m][globalPoint] = 0;
                }
            }
        }
        RDouble **qTurbulenceInterPoint = reinterpret_cast<RDouble **> (grid->GetDataPtr("turb::qInterpoint"));
        for (int iPoint = 0; iPoint < numberOfInterPoint; ++ iPoint)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qTurbulenceInterPoint[m][iPoint] = 0;
            }
        }
    }
}

void TurbSolverUnstr::ModifyQTurbNodeValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalNode = grid->GetNTotalNode();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qTurbNode"));
    int   *nodeValueSliceTurb = reinterpret_cast<int      *>(grid->GetDataPtr("nodeValueSliceTurb"));

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nodeValueSliceTurb[iNode] != 0)
        {
            for (int m = 0; m < nTurbulenceEquation; m++)
            {
                qTurbulenceNode[m][iNode] /= nodeValueSliceTurb[iNode];
            }
        }
        else
        {
            for (int m = 0; m < nTurbulenceEquation; m++)
            {
                qTurbulenceNode[m][iNode] = 0.0;
            }
        }
    }
}

void TurbSolverUnstr::ComputeFaceWeight(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *deltL = faceProxy->GetWeightL();
    RDouble *deltR = faceProxy->GetWeightR();

    grid->FaceWeight(deltL, deltR, localStart, localEnd);
}

void TurbSolverUnstr::ViscousFlux(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalFace = grid->GetNTotalFace();
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    FaceProxy *faceProxy = CreateFaceProxy(grid);
    faceProxy->SetTurbFaceValue(CreateTurbFaceValue(grid));

    int localStart,localEnd;
    localStart = 0;
    do
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > numberOfTotalFace) localEnd = numberOfTotalFace;

        ComputeFaceWeight(grid, faceProxy, localStart, localEnd);

        // GetVisFaceValue(grid, faceProxy, localStart, localEnd);

        if(tscheme != GMRES)
        {
            GetVisFaceValue(grid, faceProxy, localStart, localEnd);

            ComputeVisflux(grid, faceProxy, localStart, localEnd);
        }
#ifdef USE_GMRESSOLVER
        else
        {
            GMRES_ComputeVisflux(grid, faceProxy, localStart, localEnd);
        }
#endif
        localStart = localEnd;
    } while (localStart < numberOfTotalFace);

    delete faceProxy;
}

void TurbSolverUnstr::GetVisFaceValue(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    //! GMRESturb
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    TurbFaceValue *faceVariable = faceProxy->GetTurbFaceValue();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    string viscousName = parameters->GetViscousName();

    RDouble *mul  = faceVariable->mul;
    RDouble *mut  = faceVariable->mut;
    RDouble **mlt = faceVariable->mlt;

    int le,re,jFace;
    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        if (iFace < numberOfBoundaryFace) re = iFace + numberOfTotalCell;
        jFace  = iFace - localStart;

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

            if (IsInterface(bcType) || bcType == PHENGLEI::OVERSET)
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

    if (viscousType == ONE_EQU)
    {
        RDouble KW_sigma = parameters->GetKW_sigma();
        RDouble oSigma = 1.0 / KW_sigma;
        RDouble **qPrimitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

        for (int iFace = localStart; iFace <localEnd; ++ iFace)
        {
            le = leftCellOfFace [iFace];
            re = rightCellOfFace[iFace];
            jFace = iFace - localStart;

            RDouble nueLeft = qTurbulence[ISA][le];
            RDouble visLeft = viscousLaminar[le];
            RDouble rhoLeft = qPrimitiveVariable[IR][le]; 

            RDouble nueRight= qTurbulence[ISA][re]; 
            RDouble visRight= viscousLaminar[re];
            RDouble rhoRight= qPrimitiveVariable[IR][re];

            if(iFace < numberOfBoundaryFace)
            {
                //! move to the upper BC loop future.
                if(nueLeft < 0.0)
                {
                    nueLeft  = PositiveSA(nueLeft, visLeft, rhoLeft);

                    int bcType = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace])->GetBCType();
                    if(IsWall(bcType))
                    {
                        nueRight = -nueLeft;
                    }
                }
            }
            else
            {
                nueLeft  = PositiveSA(nueLeft, visLeft, rhoLeft);
                nueRight = PositiveSA(nueRight, visRight, rhoRight);
            }

            mlt[ISA][jFace] = oSigma * half * (visLeft / rhoLeft + visRight / rhoRight
                                             + nueLeft + nueRight );
        }    //! iFace loop
    }
    else
    {
        if (viscousName.substr(0, 13) == "2eq-kw-menter" || viscousName.substr(0, 12) == "easm-kw-2005")
        {
            RDouble *blend = reinterpret_cast<RDouble  *> (grid->GetDataPtr("blend"));

            RDouble turbSigk,turbSigw;

            RDouble KW_sigmaK1 = parameters->GetKW_sigmaK1();
            RDouble KW_sigmaW1 = parameters->GetKW_sigmaW1();
            RDouble KW_sigmaK2 = parameters->GetKW_sigmaK2();
            RDouble KW_sigmaW2 = parameters->GetKW_sigmaW2();

            for (int iFace = localStart; iFace <localEnd; ++ iFace)
            {
                le = leftCellOfFace [iFace];
                re = rightCellOfFace[iFace];
                if (iFace < numberOfBoundaryFace) re = iFace + numberOfTotalCell;
                jFace  = iFace - localStart;

                //RDouble cblend = half * (blend[le] + blend[re]);
                RDouble cblend = blend[le];
                turbSigk = cblend * KW_sigmaK1 + (1.0 - cblend) * KW_sigmaK2;
                turbSigw = cblend * KW_sigmaW1 + (1.0 - cblend) * KW_sigmaW2;

                mlt[IKE][jFace] = mul[jFace] + mut[jFace] * turbSigk;
                mlt[IKW][jFace] = mul[jFace] + mut[jFace] * turbSigw;
            }
        }
        else
        {
            RDouble turbSigk, turbSigw;

            GlobalDataBase::GetData("turbSigk", &turbSigk, PHDOUBLE, 1);
            GlobalDataBase::GetData("turbSigw", &turbSigw, PHDOUBLE, 1);

            for (int iFace = localStart; iFace <localEnd; ++ iFace)
            {
                le = leftCellOfFace [iFace];
                re = rightCellOfFace[iFace];
                if (iFace < numberOfBoundaryFace) re = iFace + numberOfTotalCell;
                jFace  = iFace - localStart;

                mlt[IKE][jFace] = mul[jFace] + mut[jFace] * turbSigk;
                mlt[IKW][jFace] = mul[jFace] + mut[jFace] * turbSigw;
            }
        }
    }
}

#ifdef USE_GMRESSOLVER
//! GMRESCoupled
template <typename T>
void TurbSolverUnstr::GMRES_GetVisFaceValue(UnstructBCSet *unstructBCSet, int iFace, int nBoundFace, int nTurbulenceEquation, RDouble KW_sigma, T rhol, T rhor, T pressurel, T pressurer, T* turbQl, T* turbQr, T& mul, T& mlt, T& viscousLaminarl, T& viscousLaminarr)
{
    RDouble coefficientOfStateEquation = GAS_SPACE::gas->GetCoefficientOfStateEquation();
    RDouble nonDimensionalSutherlandTemperature = GlobalDataBase::GetDoubleParaFromDB("tsuth");
    RDouble viscousLaminarMin;
    GlobalDataBase::GetData("visl_min", &viscousLaminarMin, PHDOUBLE, 1);

    //! compute viscousLaminar (mul)
    T tl = pressurel / (coefficientOfStateEquation * rhol);
    T tr = pressurer / (coefficientOfStateEquation * rhor);
    viscousLaminarl = tl * sqrt(tl) * (1.0 + nonDimensionalSutherlandTemperature) / (tl + nonDimensionalSutherlandTemperature);
    viscousLaminarr = tr * sqrt(tr) * (1.0 + nonDimensionalSutherlandTemperature) / (tr + nonDimensionalSutherlandTemperature);
    viscousLaminarl = max(viscousLaminarMin, viscousLaminarl);
    viscousLaminarr = max(viscousLaminarMin, viscousLaminarr);
    mul = half * (viscousLaminarl + viscousLaminarr);

    if(iFace < nBoundFace)
    {
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();
        if(!IsInterface(bcType))
        {
            mul = viscousLaminarl;
        }
    }

    RDouble oSigma = 1. / KW_sigma;
    T oRho = two / (rhol + rhor + SMALL);
    mlt = oSigma * (mul * oRho + half * (turbQl[IDX::ISA] + turbQr[IDX::ISA]));

    return;
}

//! GMRESturb
void TurbSolverUnstr::GMRES_ComputeVisflux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    RDouble **flux  = faceProxy->GetFlux();
    TurbFaceValue *faceVariable = faceProxy->GetTurbFaceValue();
    
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    
    RDouble *cellCenterX = grid->GetCellCenterX();
    RDouble *cellCenterY = grid->GetCellCenterY();
    RDouble *cellCenterZ = grid->GetCellCenterZ();

    RDouble *xfc  = grid->GetFaceCenterX();
    RDouble *yfc  = grid->GetFaceCenterY();
    RDouble *zfc  = grid->GetFaceCenterZ();

    // RDouble **mlt  = faceVariable->mlt;

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    // RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble **qPrimitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    RDouble **dDdP_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP_turb"));
    RDouble **dRdq_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dRdq_turb"));
    RDouble **dDdP = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP"));
    RDouble **dRdqCoupledTerm_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dRdqCoupledTerm_turb"));
    vector<int> AI = grid->GetJacobianAI4GMRES();
    vector<int> AJ = grid->GetJacobianAJ4GMRES();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();

    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int nLaminar = grid->GetBlockSizeOfJacobianMatrix();
    nLaminar -= nTurbulenceEquation;

    RDouble KW_sigma = parameters->GetKW_sigma();

    const RDouble cb2    = 0.622;
    const RDouble sigma  = 2.0/3.0;
    const RDouble osigma = 1.0 / sigma;

    ADReal rhol, pressurel;
    ADReal rhor, pressurer;
    ADReal viscousLaminarl, viscousLaminarr;
    ADReal mlt, mul;

    ADReal *dfd1  = new ADReal[nTurbulenceEquation];
    ADReal *dfd2  = new ADReal[nTurbulenceEquation];
    ADReal *dfdn  = new ADReal[nTurbulenceEquation];

    ADReal *f1    = new ADReal[nTurbulenceEquation];
    ADReal *f2    = new ADReal[nTurbulenceEquation];
    ADReal *fMid  = new ADReal[nTurbulenceEquation];

    // ADReal *ADmlt  = new ADReal[nTurbulenceEquation];
    ADReal *turbQl = new ADReal[nTurbulenceEquation];
    ADReal *turbQr = new ADReal[nTurbulenceEquation];
    ADReal *fluxt  = new ADReal[nTurbulenceEquation];
    ADReal *fluxtr = new ADReal[nTurbulenceEquation];

    RDouble dfluxdturbQl[nTurbulenceEquation][nTurbulenceEquation];
    RDouble dfluxdturbQr[nTurbulenceEquation][nTurbulenceEquation];
    // RDouble dfluxdmlt[nTurbulenceEquation][nTurbulenceEquation];

    RDouble dDdPlocal[nTurbulenceEquation][nTurbulenceEquation];
    RDouble dDdplocal_ns[nLaminar][nLaminar];
    // RDouble tmpmatrix[nTurbulenceEquation][nTurbulenceEquation];

    using namespace IDX;

    int le, re, jFace, colidx;
    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        jFace  = iFace - localStart;

        //! AD info
        for (int m = 0; m < nTurbulenceEquation; m++)
        {
            turbQl[m] = qTurbulence[m][le];
            turbQr[m] = qTurbulence[m][re];

            turbQl[m].diff(m, 6);
            turbQr[m].diff(m + 3, 6);
        }
        rhol      = qPrimitiveVariable[IR][le];
        rhor      = qPrimitiveVariable[IR][re];
        pressurel = qPrimitiveVariable[IP][le];
        pressurer = qPrimitiveVariable[IP][re];
        rhol.diff(1, 6);
        pressurel.diff(2, 6);
        rhor.diff(4, 6);
        pressurer.diff(5, 6);

        //! compute mul and mlt
        GMRES_GetVisFaceValue(unstructBCSet, iFace, numberOfBoundaryFace, nTurbulenceEquation, KW_sigma, rhol, rhor, pressurel, pressurer, turbQl, turbQr, mul, mlt, viscousLaminarl, viscousLaminarr);

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
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle1 = asin(distTemp) * 180.0 / PI;

        distTemp = dR / (DISTANCE(dxR, dyR, dzR) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle2 = asin(distTemp) * 180.0 / PI;

        //! Quantities at points 1 and 2
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            f1[m] = turbQl[m];
            f2[m] = turbQr[m];
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            fMid[m] = half * (f1[m] + f2[m]);
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            dfdn[m] = 0.0;
        }

        if (angle1 > 0.0 && angle2 > 0.0 && ABS(dL) > TINY && ABS(dR) > TINY)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dfd1[m] = (f1[m] - fMid[m]) / dL;
            }

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dfd2[m] = (f2[m] - fMid[m]) / dR;
            }

            RDouble dtmp = dL * dL + dR * dR;
            RDouble weightL = dL * dL / dtmp;
            RDouble weightR = 1.0 - weightL;

            for (int m = 0; m < nTurbulenceEquation; ++ m)
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

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        dfdn[m] = (f2[m] - f1[m]) / normalDist;
                    }
                }
            }
        }
        
 if(TurbEquationFormation() == FORMATION_SA_NSMB)
        {
            int m = IDX::ISA;

            ADReal flux1 = - oRefReNumber * (1.0 + cb2) * mlt * dfdn[m] * area[iFace];
            ADReal flux2 =  oRefReNumber * cb2 * osigma * dfdn[m] * area[iFace];
            fluxt[m] = flux1 + flux2 * (viscousLaminarl / rhol + turbQl[m]);

            residualTurbulence[m][le] -= fluxt[m].val();

            dfluxdturbQl[m][m] = fluxt[m].dx(0);
            dfluxdturbQr[m][m] = fluxt[m].dx(3);

            //! Jacobian components for residual of the left cell
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, le, le);
            dRdq_turb[m][colidx] += dfluxdturbQl[m][m];
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, le, le);
            dRdqCoupledTerm_turb[ISA][colidx + IR] += fluxt[m].dx(1);
            dRdqCoupledTerm_turb[ISA][colidx + IP] += fluxt[m].dx(2);


            if(iFace < numberOfBoundaryFace)
            {
                UnstructBC* bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
                int bcType = bcRegion->GetBCType();
                
                if(bcType != PHENGLEI::INTERFACE)
                {
                    int indexre = (re - nTotalCell) * nTurbulenceEquation;
                    dDdPlocal[m][m] = dDdP_turb[m][indexre + m];
                    indexre = (re - nTotalCell) * nLaminar; //! index for coupled term

                    //! firstly, Jacobian components for residual of the right cell
                    // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, le, re);
                    // dRdq_turb[m][colidx] += dfluxdturbQr[m][m];
                    // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, le, re);
                    // dRdqCoupledTerm_turb[ISA][colidx + IR] += fluxt[m].dx(4);
                    // dRdqCoupledTerm_turb[ISA][colidx + IP] += fluxt[m].dx(5);

                    // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, re, re);
                    // dRdq_turb[m][colidx] -= dfluxdturbQr[m][m];
                    // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, re, re);
                    // dRdqCoupledTerm_turb[ISA][colidx + IR] -= fluxt[m].dx(4);
                    // dRdqCoupledTerm_turb[ISA][colidx + IP] -= fluxt[m].dx(5);

                    colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, le, le);
                    dRdq_turb[m][colidx] += dfluxdturbQr[m][m] * dDdPlocal[m][m];
                    colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, le, le);
                    dRdqCoupledTerm_turb[ISA][colidx + IR] += (fluxt[ISA].dx(4) * dDdP[IR][indexre + IR] + fluxt[ISA].dx(5) * dDdP[IP][indexre + IR]);
                    dRdqCoupledTerm_turb[ISA][colidx + IU] += (fluxt[ISA].dx(4) * dDdP[IR][indexre + IU] + fluxt[ISA].dx(5) * dDdP[IP][indexre + IU]);
                    dRdqCoupledTerm_turb[ISA][colidx + IV] += (fluxt[ISA].dx(4) * dDdP[IR][indexre + IV] + fluxt[ISA].dx(5) * dDdP[IP][indexre + IV]);
                    dRdqCoupledTerm_turb[ISA][colidx + IW] += (fluxt[ISA].dx(4) * dDdP[IR][indexre + IW] + fluxt[ISA].dx(5) * dDdP[IP][indexre + IW]);
                    dRdqCoupledTerm_turb[ISA][colidx + IP] += (fluxt[ISA].dx(4) * dDdP[IR][indexre + IP] + fluxt[ISA].dx(5) * dDdP[IP][indexre + IP]);

                    // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, re, le);
                    // dRdq_turb[m][colidx] -= dfluxdturbQr[m][m] * dDdPlocal[m][m];
                    // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, re, le);
                    // dRdqCoupledTerm_turb[ISA][colidx + IR] -= (fluxt[ISA].dx(4) * dDdP[IR][indexre + IR] + fluxt[ISA].dx(5) * dDdP[IP][indexre + IR]);
                    // dRdqCoupledTerm_turb[ISA][colidx + IU] -= (fluxt[ISA].dx(4) * dDdP[IR][indexre + IU] + fluxt[ISA].dx(5) * dDdP[IP][indexre + IU]);
                    // dRdqCoupledTerm_turb[ISA][colidx + IV] -= (fluxt[ISA].dx(4) * dDdP[IR][indexre + IV] + fluxt[ISA].dx(5) * dDdP[IP][indexre + IV]);
                    // dRdqCoupledTerm_turb[ISA][colidx + IW] -= (fluxt[ISA].dx(4) * dDdP[IR][indexre + IW] + fluxt[ISA].dx(5) * dDdP[IP][indexre + IW]);
                    // dRdqCoupledTerm_turb[ISA][colidx + IP] -= (fluxt[ISA].dx(4) * dDdP[IR][indexre + IP] + fluxt[ISA].dx(5) * dDdP[IP][indexre + IP]);
                }
                else
                {
                    fluxtr[m] = flux1 + flux2 * (viscousLaminarr / rhor + turbQr[m]);

                    residualTurbulence[m][re] += fluxtr[m].val();

                    dfluxdturbQl[m][m] = fluxtr[m].dx(0);
                    dfluxdturbQr[m][m] = fluxtr[m].dx(3);

                    //! Jacobian components for residual of the right cell
                    colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, re, le);
                    dRdq_turb[m][colidx] -= dfluxdturbQl[m][m];
                    colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, re, le);
                    dRdqCoupledTerm_turb[ISA][colidx + IR] -= fluxt[m].dx(1);
                    dRdqCoupledTerm_turb[ISA][colidx + IP] -= fluxt[m].dx(2);

                    // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, re, re);
                    // dRdq_turb[m][colidx] -= dfluxdturbQr[m][m];
                    // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, re, re);
                    // dRdqCoupledTerm_turb[ISA][colidx + IR] -= fluxt[m].dx(4);
                    // dRdqCoupledTerm_turb[ISA][colidx + IP] -= fluxt[m].dx(5);  
                }
            }
            else
            {
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, le, re);
                dRdq_turb[m][colidx] += dfluxdturbQr[m][m];
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, le, re);
                dRdqCoupledTerm_turb[ISA][colidx + IR] += fluxt[m].dx(4);
                dRdqCoupledTerm_turb[ISA][colidx + IP] += fluxt[m].dx(5);    

                fluxtr[m] = flux1 + flux2 * (viscousLaminarr / rhor + turbQr[m]);

                residualTurbulence[m][re] += fluxtr[m].val();

                dfluxdturbQl[m][m] = fluxtr[m].dx(0);
                dfluxdturbQr[m][m] = fluxtr[m].dx(3);

                //! Jacobian components for residual of the right cell
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, re, le);
                dRdq_turb[m][colidx] -= dfluxdturbQl[m][m];
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, re, le);
                dRdqCoupledTerm_turb[ISA][colidx + IR] -= fluxt[m].dx(1);
                dRdqCoupledTerm_turb[ISA][colidx + IP] -= fluxt[m].dx(2);

                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, re, re);
                dRdq_turb[m][colidx] -= dfluxdturbQr[m][m];
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, re, re);
                dRdqCoupledTerm_turb[ISA][colidx + IR] -= fluxt[m].dx(4);
                dRdqCoupledTerm_turb[ISA][colidx + IP] -= fluxt[m].dx(5);  
            }
        }
        else 
        {
            printf("TurbEquationFormation(): %d", TurbEquationFormation());
            TK_Exit::ExceptionExit("GMRES RANS solver supports 1eq-sa-Original so far");
        }
    }

    delete []dfd1;
    delete []dfd2;
    delete []dfdn;

    delete []f1;
    delete []f2;
    delete []fMid;

    delete []turbQl;
    delete []turbQr;
    delete []fluxt;
    delete []fluxtr;
}

#endif
//! Bell 20130513 mod.
void TurbSolverUnstr::ComputeVisflux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int numberOfBoundaryFace = grid->GetNBoundFace();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    RDouble **flux = faceProxy->GetFlux();
    TurbFaceValue *faceVariable = faceProxy->GetTurbFaceValue();
    
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    
    RDouble *cellCenterX = grid->GetCellCenterX();
    RDouble *cellCenterY = grid->GetCellCenterY();
    RDouble *cellCenterZ = grid->GetCellCenterZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble **mlt = faceVariable->mlt;

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble **qPrimitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble skewnessAngle = parameters->GetSkewnessAngle();
    RDouble turbSkewnessAngle = parameters->GetTurbSkewnessAngle();

    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    const RDouble cb2    = 0.622;
    const RDouble sigma  = 2.0/3.0;
    const RDouble osigma = 1.0 / sigma;

    RDouble *dfd1 = new RDouble[nTurbulenceEquation];
    RDouble *dfd2 = new RDouble[nTurbulenceEquation];
    RDouble *dfdn = new RDouble[nTurbulenceEquation];

    RDouble *f1   = new RDouble[nTurbulenceEquation];
    RDouble *f2   = new RDouble[nTurbulenceEquation];
    RDouble *fMid = new RDouble[nTurbulenceEquation];

    int le, re, jFace;
    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        le = leftCellOfFace [iFace];
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
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle1 = asin(distTemp) * 180.0 / PI;

        distTemp = dR / (DISTANCE(dxR, dyR, dzR) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle2 = asin(distTemp) * 180.0 / PI;

        //! quantities at points 1 and 2
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            f1[m] = qTurbulence[m][le];
            f2[m] = qTurbulence[m][re];
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            fMid[m] = half * (f1[m] + f2[m]);
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            dfdn[m] = 0.0;
        }

        if (angle1 > turbSkewnessAngle && angle2 > turbSkewnessAngle && ABS(dL) > TINY && ABS(dR) > TINY)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dfd1[m] = (f1[m] - fMid[m]) / dL;
            }

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dfd2[m] = (f2[m] - fMid[m]) / dR;
            }

            RDouble dtmp = dL * dL + dR * dR;
            RDouble weightL = dL * dL / dtmp;
            RDouble weightR = 1.0 - weightL;

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dfdn[m] = weightL * dfd1[m] + weightR * dfd2[m];
            }

            if (iFace < numberOfBoundaryFace)
            {
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
                int bcType = bcRegion->GetBCType();

                if (bcType != PHENGLEI::INTERFACE && bcType != PHENGLEI::SYMMETRY && bcType != PHENGLEI::OVERSET)
                {
                    RDouble normalDist = (cellCenterX[re] - cellCenterX[le]) * nxs +
                                         (cellCenterY[re] - cellCenterY[le]) * nys +
                                         (cellCenterZ[re] - cellCenterZ[le]) * nzs;

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        dfdn[m] = (f2[m] - f1[m]) / normalDist;
                    }
                }
            }
        }

        if(TurbEquationFormation() == FORMATION_SA_CFL3D)
        {
            int m = IDX::ISA;

            RDouble nueLeft = qTurbulence[IDX::ISA][le];
            RDouble visLeft = viscousLaminar[le];
            RDouble rhoLeft = qPrimitiveVariable[IDX::IR][le]; 
            RDouble nueRight= qTurbulence[IDX::ISA][re]; 
            RDouble visRight= viscousLaminar[re];
            RDouble rhoRight= qPrimitiveVariable[IDX::IR][re];
            nueLeft  = PositiveSA(nueLeft, visLeft, rhoLeft);
            nueRight = PositiveSA(nueRight, visRight, rhoRight);
            if(iFace < numberOfBoundaryFace)
            {
                int bcType = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace])->GetBCType();
                if(IsWall(bcType))
                {
                    nueRight = -nueLeft;
                }
            }

            RDouble viscousFace = half * (viscousLaminar[le] + viscousLaminar[re]);
            RDouble rhoFace     = half * (qPrimitiveVariable[IDX::IR][le] + qPrimitiveVariable[IDX::IR][re]);
            RDouble nueFace     = viscousFace / rhoFace;
            RDouble nueSAFace   = half * (nueLeft + nueRight);

            //! Term1 && term2.
            RDouble flux1 = - oRefReNumber * osigma * (nueFace + (1.0 + cb2) * nueSAFace) * dfdn[m] * area[iFace];
            RDouble flux2 =   oRefReNumber * osigma * cb2 * dfdn[m] * area[iFace];

            residualTurbulence[m][le] -= flux1;
            residualTurbulence[m][le] -= flux2 * qTurbulence[IDX::ISA][le];

            if(iFace >= numberOfBoundaryFace)
            {
                residualTurbulence[m][re] += flux1;
                residualTurbulence[m][re] += flux2 * qTurbulence[IDX::ISA][re];
            }
        }
        else if(TurbEquationFormation() == FORMATION_SA_NSMB)
        {
            int m = IDX::ISA;

            RDouble flux1 = - oRefReNumber * (1.0 + cb2) * mlt[m][jFace] * dfdn[m] * area[iFace];
            RDouble flux2 =  oRefReNumber * cb2 * osigma * dfdn[m] * area[iFace];

            RDouble nueLeft = qTurbulence[IDX::ISA][le];
            RDouble visLeft = viscousLaminar[le];
            RDouble rhoLeft = qPrimitiveVariable[IDX::IR][le];
            RDouble nueRight= qTurbulence[IDX::ISA][re];
            RDouble visRight= viscousLaminar[re];
            RDouble rhoRight= qPrimitiveVariable[IDX::IR][re];
            nueLeft  = PositiveSA(nueLeft, visLeft, rhoLeft);
            nueRight = PositiveSA(nueRight, visRight, rhoRight);

            residualTurbulence[m][le] -= flux1;
            residualTurbulence[m][le] -= flux2 * (viscousLaminar[le] / qPrimitiveVariable[IDX::IR][le] + nueLeft);

            if(iFace >= numberOfBoundaryFace)
            {
                residualTurbulence[m][re] += flux1;
                residualTurbulence[m][re] += flux2 * (viscousLaminar[re] / qPrimitiveVariable[IDX::IR][re] + nueRight);
            }
        }
        else if(TurbEquationFormation() == FORMATION_ORGSA_SST)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                flux[m][jFace] = - oRefReNumber * mlt[m][jFace] * dfdn[m] * area[iFace];
            }

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                residualTurbulence[m][le] -= flux[m][jFace];

                if(iFace >= numberOfBoundaryFace)
                {
                    residualTurbulence[m][re] += flux[m][jFace];
                }
            }
        }
    }

    delete [] dfd1;    dfd1 = nullptr;
    delete [] dfd2;    dfd2 = nullptr;
    delete [] dfdn;    dfdn = nullptr;

    delete [] f1;    f1 = nullptr;
    delete [] f2;    f2 = nullptr;
    delete [] fMid;    fMid = nullptr;
}

void TurbSolverUnstr::Diagonal(Grid *grid)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    if (viscousType == ONE_EQU)
    {
        SpectrumRadiusOfOneEquation(grid);
    }
    else if (viscousType == TWO_EQU)
    {
        SpectrumRadiusOfTwoEquation(grid);
    }
}

void TurbSolverUnstr::SpectrumRadiusOfOneEquation(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalFace    = grid->GetNTotalFace();
    int numberOfBoundaryFace = grid->GetNBoundFace();

    RDouble **q         = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();

    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    //! Initialize some constants
    RDouble sigma  = 2.0 / 3.0;
    RDouble osigma = one / sigma;
    const RDouble cb2 = 0.622;

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *vol  = grid->GetCellVolume();
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vgn  = grid->GetFaceNormalVelocity();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    using namespace IDX;
    //! rhoFlag, is used to compatible both original (0) and conservative (1) variable form.
    int rhoFlag = parameters->UsingConservativeForm();

    RDouble invJacobianL, invJacobianR;
    for (int iFace = 0; iFace < numberOfTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];

        //! Consider vgn here!

        //! Inviscid Jacobian.
        RDouble rhoL = q[IR][le];
        RDouble uL = q[IU][le];
        RDouble vL = q[IV][le];
        RDouble wL = q[IW][le];
        //RDouble vnL = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL;
        RDouble vnL = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];
        RDouble rhoCoefficientL = rhoFlag * one +  (1 - rhoFlag) * rhoL;
        if (vnL > 0.0)
        {
            invJacobianL = rhoCoefficientL * vnL * area[iFace];
        }
        else
        {
            invJacobianL = 0.0;
        }

        RDouble rhoR = q[IR][re];
        RDouble uR = q[IU][re];
        RDouble vR = q[IV][re];
        RDouble wR = q[IW][re];
        //RDouble vnR = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR;
        RDouble vnR = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR - vgn[iFace];

        RDouble rhoCoefficientR = rhoFlag * one +  (1 - rhoFlag) * rhoR;
        if (vnR < 0.0)
        {
            invJacobianR = rhoCoefficientR * vnR * area[iFace];
        }
        else
        {
            invJacobianR = 0.0;
        }

        spectrumTurbulence[ISA][le] += invJacobianL;
        spectrumTurbulence[ISA][re] -= invJacobianR;

        matrixTurbulenceRight[ISA][iFace] = invJacobianR;
        matrixTurbulenceLeft[ISA][iFace]  = -invJacobianL;

        //! Viscous Jacobian.
        RDouble orL = 1.0 / (rhoL + SMALL);
        RDouble orR = 1.0 / (rhoR + SMALL);

        RDouble nueL = qTurbulence[ISA][le];
        RDouble nueR = qTurbulence[ISA][re];

        nueL = PositiveSA(nueL, viscousLaminar[le], rhoL);
        nueR = PositiveSA(nueR, viscousLaminar[re], rhoR);
        if(iFace < numberOfBoundaryFace)
        {
            int bcType = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace])->GetBCType();
            if(IsWall(bcType))
            {
                nueR = -nueL;
            }
        }

        RDouble nueffL = viscousLaminar[le] * orL + ABS(nueL);
        RDouble nueffR = viscousLaminar[re] * orR + ABS(nueR);

        RDouble nueff   = half * (nueffL + nueffR);

        RDouble volume  = half * (vol[le] + vol[re]);
        RDouble oVolume = one / volume;

        RDouble ns2     = area[iFace] * area[iFace];
        RDouble lengthScale = ns2 * oVolume;

        //! Since the accuracy of the viscous spectrum radius has little effect on time advancing,
        //! the approximate method is used here.
        if (TurbEquationFormation() == FORMATION_SA_CFL3D)
        {
            RDouble viscousFace = half * (viscousLaminar[le] + viscousLaminar[re]);
            RDouble rhoFace     = half * (q[IDX::IR][le] + q[IDX::IR][re]);
            RDouble nueFace     = viscousFace / rhoFace;
            RDouble nueSAFace   = half * (nueL + nueR);

            //! Term1 && term2.
            RDouble term1 =   oRefReNumber * osigma * (1.0 + cb2) * (nueFace + nueSAFace) * lengthScale;
            RDouble term2 = - oRefReNumber * osigma * cb2 * lengthScale;

            RDouble jocabianL = (term1 + term2 * nueL);
            RDouble jocabianR = (term1 + term2 * nueR);

            //! For iCell itself, Jacobian is negative, move it to left side, neg-neg is positive now.
            spectrumTurbulence[ISA][le] += jocabianL;
            spectrumTurbulence[ISA][re] += jocabianR;

            matrixTurbulenceRight[ISA][iFace] -= jocabianR;    //! Right cell of iFace.
            matrixTurbulenceLeft[ISA][iFace]  -= jocabianL;    //! Left cell of iFace.
        }
        else if (TurbEquationFormation() == FORMATION_SA_NSMB)
        {
            spectrumTurbulence[ISA][le] += oRefReNumber * oVolume * (1.0 + cb2) * (osigma * nueff * ns2);
            spectrumTurbulence[ISA][re] += oRefReNumber * oVolume * (1.0 + cb2) * (osigma * nueff * ns2);

            matrixTurbulenceLeft[ISA][iFace]  += - oRefReNumber * oVolume * (1.0 + cb2) * (osigma * nueff * ns2);
            matrixTurbulenceRight[ISA][iFace]  += - oRefReNumber * oVolume * (1.0 + cb2) * (osigma * nueff * ns2);

            //! �����nsmb��2��ò��ҲӦ�üӵ�spectrumTurbulence��
            RDouble coeffLe, coeffRe;
            coeffLe = cb2 * osigma * nueffL;
            coeffRe = cb2 * osigma * nueffR;
            matrixTurbulenceLeft[ISA][iFace] += oRefReNumber * oVolume * ns2 * coeffLe;
            matrixTurbulenceRight[ISA][iFace] += oRefReNumber * oVolume * ns2 * coeffRe;
        }
        else if (TurbEquationFormation() == FORMATION_ORGSA_SST)
        {
            spectrumTurbulence[ISA][le] += oRefReNumber * oVolume * (osigma * nueff * ns2);
            spectrumTurbulence[ISA][re] += oRefReNumber * oVolume * (osigma * nueff * ns2);

            matrixTurbulenceLeft[ISA][iFace] += - oRefReNumber * oVolume * (osigma * nueff * ns2);
            matrixTurbulenceRight[ISA][iFace] += - oRefReNumber * oVolume * (osigma * nueff * ns2);
        }
    }
}

void TurbSolverUnstr::SpectrumRadiusOfTwoEquation(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalFace = grid->GetNTotalFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));
    
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));

    RDouble *vol  = grid->GetCellVolume();
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vgn  = grid->GetFaceNormalVelocity();

    RDouble uL, vL, wL, vnL, absVnL;
    RDouble uR, vR, wR, vnR, absVnR;

    RDouble vn, absVn, ns2, rho, mul, mut, cblend, volume, oVolume;

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    vector<RDouble> work(nTurbulenceEquation);

    using namespace IDX;

    for (int iFace = 0; iFace < numberOfTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];

        //! Calculate velocity at the cell interface.
        uL = q[IU][le];
        vL = q[IV][le];
        wL = q[IW][le];
        //vnL = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL;
        vnL = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];

        absVnL = ABS(vnL);

        uR = q[IU][re];
        vR = q[IV][re];
        wR = q[IW][re];
        //vnR = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR;
        vnR = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR - vgn[iFace];

        absVnR = ABS(vnR);

        vn = half * (vnL + vnR);
        absVn = ABS(vn);

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            spectrumTurbulence[m][le] += area[iFace] * half * absVn;
            spectrumTurbulence[m][re] += area[iFace] * half * absVn;

            matrixTurbulenceRight[m][iFace]  = area[iFace] * half * ( vnR - absVn);
            matrixTurbulenceLeft[m][iFace]   = area[iFace] * half * (- vnL - absVn);
        }
    }

    RDouble *blend = reinterpret_cast<RDouble  *> (grid->GetDataPtr("blend"));

    RDouble turbSigk, turbSigw;

    RDouble KW_sigmaK1 = parameters->GetKW_sigmaK1();
    RDouble KW_sigmaW1 = parameters->GetKW_sigmaW1();
    RDouble KW_sigmaK2 = parameters->GetKW_sigmaK2();
    RDouble KW_sigmaW2 = parameters->GetKW_sigmaW2();

    for (int iFace = 0; iFace < numberOfTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];
        rho = half * (q[IR][le] + q[IR][re]);
        mul = half * (viscousLaminar[le] + viscousLaminar[re]);
        mut = half * (viscousTurbulence[le] + viscousTurbulence[re]);

        cblend    = half * (blend[le] + blend[re]);
        turbSigk = cblend * KW_sigmaK1 + (1.0 - cblend) * KW_sigmaK2;
        turbSigw = cblend * KW_sigmaW1 + (1.0 - cblend) * KW_sigmaW2;

        volume  = half * (vol[le] + vol[re]);
        oVolume = one / volume;

        ns2     = area[iFace] * area[iFace];
        work[IKE] = (mul + mut * turbSigk) / rho * ns2 * oVolume;
        work[IKW] = (mul + mut * turbSigw) / rho * ns2 * oVolume;

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            spectrumTurbulence[m][le] += oRefReNumber * work[m];
            spectrumTurbulence[m][re] += oRefReNumber * work[m];

            matrixTurbulenceLeft[m][iFace]  += - oRefReNumber * work[m];
            matrixTurbulenceRight[m][iFace] += - oRefReNumber * work[m];
        }
    }
}

void TurbSolverUnstr::GetQlQr(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    GetQlQrTurbulence(grid, faceProxy, localStart, localEnd);

    FaceProxy *faceProxyNS = faceProxy->GetNext();
    GetQlQrNS(grid, faceProxyNS, localStart, localEnd);

}

void TurbSolverUnstr::GetQlQrTurbulence(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulenceVarables = reinterpret_cast <RDouble **> (grid->GetDataPtr("q_turb"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    //! GMRESturb
    faceProxy->SetlocalStart(localStart);
    faceProxy->SetlocalEnd(localEnd);

    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        int le, re, jFace;
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        jFace = iFace - localStart;

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qL[m][jFace] = qTurbulenceVarables[m][le];
            qR[m][jFace] = qTurbulenceVarables[m][re];
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

        if (bcType == PHENGLEI::INTERFACE || bcType == PHENGLEI::OVERSET)
        {
            continue;
        }

        if (bcType == PHENGLEI::SYMMETRY)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qR[m][jFace] = qL[m][jFace];
            }
        }
        else
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                RDouble temp = half * (qR[m][jFace] + qL[m][jFace]);
                qR[m][jFace] = temp;
                qL[m][jFace] = temp;
            }
        }
    }
}

void TurbSolverUnstr::GetQlQrNS(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation = nl + nchem;
   
    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        int le, re, jFace;
        le = leftCellOfFace [iFace];
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

        if (bcType == PHENGLEI::INTERFACE || bcType == PHENGLEI::OVERSET)
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
#ifdef USE_GMRESSOLVER
//! GMRESCoupled
template <typename T>
void TurbSolverUnstr::FixBoundaryQlQr(UnstructBCSet *unstructBCSet, int iFace, int nTurbulenceEquation, int nEquation, T *primQl, T *primQr, T *turbQl, T *turbQr)
{
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
    int bcType = bcRegion->GetBCType();

    if(bcType == PHENGLEI::INTERFACE)
    {
        return;
    }
    
    if(bcType == PHENGLEI::SYMMETRY)
    {
        //! Turbulence
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            turbQr[m] = turbQl[m];
        }
        //! NS
        for (int m = 0; m < nEquation; ++ m)
        {
            T temp = (primQl[m] + primQr[m]) * half;
            primQl[m] = temp;
            primQr[m] = temp;
        }
    }
    else
    {
        //! Turbulence
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            T temp = (turbQl[m] + turbQr[m]) * half;
            turbQl[m] = temp;
            turbQr[m] = temp;
        }
        //! NS
        for (int m = 0; m < nEquation; ++ m)
        {
            T temp = (primQl[m] + primQr[m]) * half;
            primQl[m] = temp;
            primQr[m] = temp;
        }
    }
} 
#endif

void TurbSolverUnstr::InviscidFlux(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalFace = grid->GetNTotalFace();

    FaceProxy *faceProxy = CreateFaceProxy(grid);

    //! GMRESCoupled
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");

    int localStart, localEnd;

    localStart = 0;
    do
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > numberOfTotalFace) localEnd = numberOfTotalFace;

        if(tscheme != GMRES)
        {
            GetQlQr(grid, faceProxy, localStart, localEnd);
        }

        ComputeInviscidFlux(grid, faceProxy, localStart, localEnd);

        LoadFlux(grid, faceProxy, localStart, localEnd);

        localStart = localEnd;

    } while (localStart < numberOfTotalFace);

    delete faceProxy;
}

void TurbSolverUnstr::ComputeInviscidFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");

    if(tscheme != GMRES)
    {
        NNDFlux(grid, faceProxy, localStart, localEnd);
    }
#ifdef USE_GMRESSOLVER  
    else
    {
        GMRES_NNDFlux(grid, faceProxy, localStart, localEnd);
    }
 #endif
}
#ifdef USE_GMRESSOLVER  
void TurbSolverUnstr::GMRES_NNDFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    int nBoundFace = grid->GetNBoundFace();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *vgn  = grid->GetFaceNormalVelocity();
    RDouble *area = grid->GetFaceArea();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int viscousType = parameters->GetViscousType();

    int subJacMatSize = grid->GetBlockSizeOfJacobianMatrix();
    int nLaminar = subJacMatSize - nTurbulenceEquation;
    
    // RDouble **qL = faceProxy->GetQL();
    // RDouble **qR   = faceProxy->GetQR();
    RDouble **flux = faceProxy->GetFlux();

    int nMid = 0;
    int colidx;

    //! Check if there are boundary faces, if no 
    if (localStart >= nBoundFace)
    {
        nMid = localStart;
    }
    else if (localEnd <= nBoundFace)
    {
        //! If they are all boundary faces.
        nMid = localEnd;
    }
    else
    {
        //! Part of them are boundary faces.
        nMid = nBoundFace;
    }

    FaceProxy *faceProxyNS = faceProxy->GetNext();

    // RDouble **qPrimitiveVariableL = faceProxyNS->GetQL();
    // RDouble **qPrimitiveVariableR = faceProxyNS->GetQR();
    RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));
    RDouble **q_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("q_turb"));


    ADReal *primitiveVariableL = new ADReal[nTurbulenceEquation];
    ADReal *primitiveVariableR = new ADReal[nTurbulenceEquation];
    ADReal *turbQl             = new ADReal[nTurbulenceEquation];
    ADReal *turbQr             = new ADReal[nTurbulenceEquation];
    ADReal *primQl             = new ADReal[nLaminar];
    ADReal *primQr             = new ADReal[nLaminar];
    ADReal *turbQle            = new ADReal[nTurbulenceEquation];
    ADReal *turbQre            = new ADReal[nTurbulenceEquation];
    ADReal *primQle            = new ADReal[nLaminar];
    ADReal *primQre            = new ADReal[nLaminar];
    ADReal *f = new ADReal[nTurbulenceEquation];

    RDouble **dRdq_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dRdq_turb"));
    RDouble **dDdP_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP_turb"));
    RDouble **dDdP = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP"));
    RDouble **dRdqCoupledTerm_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dRdqCoupledTerm_turb"));
    RDouble dFluxdturbQl[nTurbulenceEquation][nTurbulenceEquation];
    RDouble dFluxdturbQr[nTurbulenceEquation][nTurbulenceEquation];
    RDouble dDdPlocal[nTurbulenceEquation][nTurbulenceEquation];
    RDouble dFluxdpBC[nTurbulenceEquation][nTurbulenceEquation];
    RDouble dDdPlocal_ns[nLaminar][nLaminar];

    vector<int> AI = grid->GetJacobianAI4GMRES();
    vector<int> AJ = grid->GetJacobianAJ4GMRES();
    int nTotalCell = grid->GetNTotalCell();

    ADReal rL,uL,vL,wL,vnL,vnLL,rR,uR,vR,wR,vnR,vnRR;
    ADReal rhoFlag = 0, rhoCoefficientL, rhoCoefficientR;

    using namespace IDX;

    if (viscousType == ONE_EQU)
    {
        //! rhoFlag, is used to compatible both original and conversative variables.
        //! For SA equation, using original variable; else, using conversative variables.
        rhoFlag = 1;
    }

    // boundary face
    for (int iFace = localStart; iFace <nMid; ++ iFace)
    {
        int le, re, jFace;
        jFace = iFace - localStart;
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        // AD init
        for (int m = 0; m < nTurbulenceEquation; m++)
        {
            turbQle[m] = q_turb[m][le];
            turbQre[m] = q_turb[m][re];
            turbQle[m].diff(m, 2 * subJacMatSize);
            turbQre[m].diff(m + subJacMatSize, 2 * subJacMatSize);
        }
        for (int m = 0; m < nLaminar; m++)
        {
            primQle[m] = q[m][le];
            primQre[m] = q[m][re];
            primQle[m].diff(nTurbulenceEquation + m, 2 * subJacMatSize);
            primQre[m].diff(nTurbulenceEquation + m + subJacMatSize, 2 * subJacMatSize);
        }

        //! GMRESCoupled get QlQr
        for (int m = 0; m < nTurbulenceEquation; m++)
        {
            turbQl[m] = turbQle[m];
            turbQr[m] = turbQre[m];
        }
        for (int m = 0; m < nLaminar; m++)
        {
            primQl[m] = primQle[m];
            primQr[m] = primQre[m];
        }
        FixBoundaryQlQr(unstructBCSet, iFace, nTurbulenceEquation, nLaminar, primQl, primQr, turbQl, turbQr);

        rL = primQl[IR];
        uL = primQl[IU];
        vL = primQl[IV];
        wL = primQl[IW];

        rR = primQr[IR];
        uR = primQr[IU];
        vR = primQr[IV];
        wR = primQr[IW];

        vnL  = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];
        vnR  = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR - vgn[iFace];

        vnLL = half * (vnL + ABS(vnL));
        vnRR = half * (vnR - ABS(vnR));

        rhoCoefficientL = rhoFlag * one + (one - rhoFlag) * rL;
        rhoCoefficientR = rhoFlag * one + (one - rhoFlag) * rR;

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            primitiveVariableL[m] = rhoCoefficientL * turbQl[m];
            primitiveVariableR[m] = rhoCoefficientR * turbQr[m];
        }

        //! Consider vgn here!
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            f[m] = (vnLL * primitiveVariableL[m] + vnRR * primitiveVariableR[m]);
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            f[m] *= area[iFace];
            flux[m][jFace] = f[m].val();            
            for (int n = 0; n < nTurbulenceEquation; n++)
            {
                dFluxdturbQl[m][n] = f[m].dx(n);
                dFluxdturbQr[m][n] = f[m].dx(subJacMatSize + n);
            }
        }

        // //! consider boundary
        // int indexre = (re - nTotalCell) * nTurbulenceEquation;
        // for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        // {
        //     for (int indexJ = 0; indexJ < nTurbulenceEquation; indexJ++)
        //     {
        //         dDdPlocal[indexI][indexJ] = dDdP_turb[indexI][indexre + indexJ];
        //     }
        // }
        // indexre = (re - nTotalCell) * nLaminar;
        // for (int indexI = 0; indexI < nLaminar; indexI++)
        // {
        //     for (int indexJ = 0; indexJ < nLaminar; indexJ++)
        //     {
        //         dDdPlocal_ns[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
        //     }
        // }
        

        //! for boundary condition except symmetry and interface
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        // //! turbulence
        // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, re, le);
        // for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        // {
        //     for (int indexJ = 0; indexJ < nTurbulenceEquation; indexJ++)
        //     {
        //         dRdq_turb[indexI][colidx + indexJ] -= dFluxdturbQl[indexI][indexJ];
        //         //! consider boundary
        //         for (int indexK = 0; indexK < nTurbulenceEquation; indexK++)
        //         {
        //             dRdq_turb[indexI][colidx + indexJ] -= dFluxdturbQr[indexI][indexK] * dDdPlocal[indexK][indexJ];                           
        //         }
        //     }
        // }
        // //! NS
        // colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, re, le);
        // for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        // {
        //     for (int indexJ = 0; indexJ < nLaminar; indexJ++)
        //     {
        //         dRdqCoupledTerm_turb[indexI][colidx + indexJ] -= f[indexI].dx(nTurbulenceEquation + indexJ);
        //         //! consider boundary
        //         for (int indexK = 0; indexK < nLaminar; indexK++)
        //         {
        //             dRdqCoupledTerm_turb[indexI][colidx + indexJ] -= f[indexI].dx(nTurbulenceEquation + subJacMatSize + indexK) * dDdPlocal_ns[indexK][indexJ];
        //         }
        //     }
        // }

        if (bcType != PHENGLEI::INTERFACE)
        {
            //! consider boundary
            int indexre = (re - nTotalCell) * nTurbulenceEquation;
            for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nTurbulenceEquation; indexJ++)
                {
                    dDdPlocal[indexI][indexJ] = dDdP_turb[indexI][indexre + indexJ];
                }
            }
            indexre = (re - nTotalCell) * nLaminar;
            for (int indexI = 0; indexI < nLaminar; indexI++)
            {
                for (int indexJ = 0; indexJ < nLaminar; indexJ++)
                {
                    dDdPlocal_ns[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
                }
            }
        
        //! turbulence
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, le, le);
        for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        {
            for (int indexJ = 0; indexJ < nTurbulenceEquation; indexJ++)
            {
                dRdq_turb[indexI][colidx + indexJ] += dFluxdturbQl[indexI][indexJ];
                //! consider boundary
                for (int indexK = 0; indexK < nTurbulenceEquation; indexK++)
                {
                    dRdq_turb[indexI][colidx + indexJ] += dFluxdturbQr[indexI][indexK] * dDdPlocal[indexK][indexJ];
                }
            }
        }
        //! NS
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, le, le);
        for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        {
            for (int indexJ = 0; indexJ < nLaminar; indexJ++)
            {
                dRdqCoupledTerm_turb[indexI][colidx + indexJ] += f[indexI].dx(nTurbulenceEquation + indexJ);
                //! consider boundary
                for (int indexK = 0; indexK < nLaminar; indexK++)
                {
                    dRdqCoupledTerm_turb[indexI][colidx + indexJ] += f[indexI].dx(nTurbulenceEquation + subJacMatSize + indexK) * dDdPlocal_ns[indexK][indexJ];
                }            
            }
        }
        }
        else
        {
            //! turbulence
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, re, le);
            for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nTurbulenceEquation; indexJ++)
                {
                    dRdq_turb[indexI][colidx + indexJ] -= dFluxdturbQl[indexI][indexJ];
                }
            }
            //! NS
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, re, le);
            for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nLaminar; indexJ++)
                {
                    dRdqCoupledTerm_turb[indexI][colidx + indexJ] -= f[indexI].dx(nTurbulenceEquation + indexJ);
                }
            }
        
            //! turbulence
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, le, le);
            for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nTurbulenceEquation; indexJ++)
                {
                    dRdq_turb[indexI][colidx + indexJ] += dFluxdturbQl[indexI][indexJ];
                }
            }
            //! NS
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, le, le);
            for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nLaminar; indexJ++)
                {
                    dRdqCoupledTerm_turb[indexI][colidx + indexJ] += f[indexI].dx(nTurbulenceEquation + indexJ);
                }
            }
        }
    }
    
    //! interior face
    for (int iFace = nMid; iFace <localEnd; ++ iFace)
    {
        int le, re, jFace;
        jFace = iFace - localStart;
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        // AD init
        for (int m = 0; m < nTurbulenceEquation; m++)
        {
            turbQle[m] = q_turb[m][le];
            turbQre[m] = q_turb[m][re];
            turbQle[m].diff(m, 2 * subJacMatSize);
            turbQre[m].diff(m + subJacMatSize, 2 * subJacMatSize);
        }
        for (int m = 0; m < nLaminar; m++)
        {
            primQle[m] = q[m][le];
            primQre[m] = q[m][re];
            primQle[m].diff(nTurbulenceEquation + m, 2 * subJacMatSize);
            primQre[m].diff(nTurbulenceEquation + m + subJacMatSize, 2 * subJacMatSize);
        }

        //! GMRESCoupled get QlQr
        for (int m = 0; m < nTurbulenceEquation; m++)
        {
            turbQl[m] = turbQle[m];
            turbQr[m] = turbQre[m];
        }
        for (int m = 0; m < nLaminar; m++)
        {
            primQl[m] = primQle[m];
            primQr[m] = primQre[m];
        }

        rL = primQl[IR];
        uL = primQl[IU];
        vL = primQl[IV];
        wL = primQl[IW];

        rR = primQr[IR];
        uR = primQr[IU];
        vR = primQr[IV];
        wR = primQr[IW];

        vnL  = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];
        vnR  = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR - vgn[iFace];

        vnLL = half * (vnL + ABS(vnL));
        vnRR = half * (vnR - ABS(vnR));

        rhoCoefficientL = rhoFlag * one + (one - rhoFlag) * rL;
        rhoCoefficientR = rhoFlag * one + (one - rhoFlag) * rR;

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            primitiveVariableL[m] = rhoCoefficientL * turbQl[m];
            primitiveVariableR[m] = rhoCoefficientR * turbQr[m];
        }

        //! Consider vgn here!
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            f[m] = (vnLL * primitiveVariableL[m] + vnRR * primitiveVariableR[m]);
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            f[m] *= area[iFace];
            flux[m][jFace] = f[m].val();            
            for (int n = 0; n < nTurbulenceEquation; n++)
            {
                dFluxdturbQl[m][n] = f[m].dx(n);
                dFluxdturbQr[m][n] = f[m].dx(subJacMatSize + n);
            }
        }

        //! turbulence
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, re, le);
        for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        {
            for (int indexJ = 0; indexJ < nTurbulenceEquation; indexJ++)
            {
                dRdq_turb[indexI][colidx + indexJ] -= dFluxdturbQl[indexI][indexJ];
            }
        }
        //! NS
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, re, le);
        for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        {
            for (int indexJ = 0; indexJ < nLaminar; indexJ++)
            {
                dRdqCoupledTerm_turb[indexI][colidx + indexJ] -= f[indexI].dx(nTurbulenceEquation + indexJ);
            }
        }
        
        //! turbulence
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, le, le);
        for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        {
            for (int indexJ = 0; indexJ < nTurbulenceEquation; indexJ++)
            {
                dRdq_turb[indexI][colidx + indexJ] += dFluxdturbQl[indexI][indexJ];
            }
        }
        //! NS
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, le, le);
        for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        {
            for (int indexJ = 0; indexJ < nLaminar; indexJ++)
            {
                dRdqCoupledTerm_turb[indexI][colidx + indexJ] += f[indexI].dx(nTurbulenceEquation + indexJ);
            }
        }

        //! turbulence
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, re, re);
        for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        {
            for (int indexJ = 0; indexJ < nTurbulenceEquation; indexJ++)
            {
                dRdq_turb[indexI][colidx + indexJ] -= dFluxdturbQr[indexI][indexJ];
            }
        }
        //! NS
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, re, re);
        for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        {
            for (int indexJ = 0; indexJ < nLaminar; indexJ++)
            {
                dRdqCoupledTerm_turb[indexI][colidx + indexJ] -= f[indexI].dx(nTurbulenceEquation + subJacMatSize + indexJ);
            }
        }

        //! turbulence
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nTurbulenceEquation, le, re);
        for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        {
            for (int indexJ = 0; indexJ < nTurbulenceEquation; indexJ++)
            {
                dRdq_turb[indexI][colidx + indexJ] += dFluxdturbQr[indexI][indexJ];
            }
        }
        //! NS
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nLaminar, le, re);
        for (int indexI = 0; indexI < nTurbulenceEquation; indexI++)
        {
            for (int indexJ = 0; indexJ < nLaminar; indexJ++)
            {
                dRdqCoupledTerm_turb[indexI][colidx + indexJ] += f[indexI].dx(nTurbulenceEquation + subJacMatSize + indexJ);
            }
        }
    }
   
    delete []primitiveVariableL;
    delete []primitiveVariableR;
    delete []f;
    delete []turbQl;
    delete []turbQr;
    delete []turbQle;
    delete []turbQre;
    delete []primQl;
    delete []primQr;
    delete []primQle;
    delete []primQre;
}
 #endif
void TurbSolverUnstr::NNDFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *vgn  = grid->GetFaceNormalVelocity();
    RDouble *area = grid->GetFaceArea();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qL   = faceProxy->GetQL();
    RDouble **qR   = faceProxy->GetQR();
    RDouble **flux = faceProxy->GetFlux();

    FaceProxy *faceProxyNS = faceProxy->GetNext();

    RDouble **qPrimitiveVariableL = faceProxyNS->GetQL();
    RDouble **qPrimitiveVariableR = faceProxyNS->GetQR();

    RDouble *qTurbL = new RDouble[nTurbulenceEquation];
    RDouble *qTurbR = new RDouble[nTurbulenceEquation];
    RDouble *f = new RDouble[nTurbulenceEquation];

    RDouble rL, uL, vL, wL, vnL, vnLL, rR, uR, vR, wR, vnR, vnRR;
    RDouble rhoCoefficientL, rhoCoefficientR;

    using namespace IDX;
    //! rhoFlag, is used to compatible both original (0) and conservative (1) variable form.
    int rhoFlag = parameters->UsingConservativeForm();

    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        int le, re, jFace;
        jFace = iFace - localStart;
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        rL = qPrimitiveVariableL[IR][jFace];
        uL = qPrimitiveVariableL[IU][jFace];
        vL = qPrimitiveVariableL[IV][jFace];
        wL = qPrimitiveVariableL[IW][jFace];

        rR = qPrimitiveVariableR[IR][jFace];
        uR = qPrimitiveVariableR[IU][jFace];
        vR = qPrimitiveVariableR[IV][jFace];
        wR = qPrimitiveVariableR[IW][jFace];

        vnL  = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];
        vnR  = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR - vgn[iFace];

        vnLL = half * (vnL + ABS(vnL));
        vnRR = half * (vnR - ABS(vnR));

        rhoCoefficientL = (1 - rhoFlag) * one + rhoFlag * rL;
        rhoCoefficientR = (1 - rhoFlag) * one + rhoFlag * rR;

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qTurbL[m] = rhoCoefficientL * qL[m][jFace];
            qTurbR[m] = rhoCoefficientR * qR[m][jFace];
        }

        //! Consider vgn here!
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            f[m] = (vnLL * qTurbL[m] + vnRR * qTurbR[m]);
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            flux[m][jFace] = f[m] * area[iFace];
        }
    }

    DelPointer(qTurbL);
    DelPointer(qTurbR);
    DelPointer(f);
}

void TurbSolverUnstr::GetResidual(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    RDouble **residualTurbulence  = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    PHSPACE::ComputeResidualonGrid(grid, actkey, residualTurbulence, nTurbulenceEquation);
}

void TurbSolverUnstr::CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->MoveToBegin();
    if (0 == nEquation)
    {
        RDouble *fieldSend = reinterpret_cast <RDouble *> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int sourceCell;
                int iFace = interfaceIndexContainerForSend[iLocalFace];
                grid->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
                PHWrite(dataContainer, fieldSend[sourceCell]);
            }
        }
    }
    else
    {
        RDouble **fieldSend = reinterpret_cast <RDouble **> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int sourceCell;
                int iFace = interfaceIndexContainerForSend[iLocalFace];
                grid->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHWrite(dataContainer, fieldSend[m][sourceCell]);
                }
            }
        }
    }
}

void TurbSolverUnstr::DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    
    dataContainer->MoveToBegin();

    if (0 == nEquation)
    {
        RDouble *fieldReceive = reinterpret_cast <RDouble *> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int targetCell;
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);
                PHRead(dataContainer, fieldReceive[targetCell]);
            }
        }
    }
    else
    {
        RDouble **fieldReceive = reinterpret_cast <RDouble **> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int targetCell;
                grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHRead(dataContainer, fieldReceive[m][targetCell]);
                }
            }
        }
    }
}

void TurbSolverUnstr::Boundary(Grid *gridIn)
{
    using namespace PHENGLEI;
    UnstructGrid *grid = UnstructGridCast(gridIn);
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int isOverset = parameters->GetIsOverLapping();

    if (isOverset)
    {
        ReSetOversetBoundary(grid);
    }

    int wallFunctionType = parameters->GetWallFunctionType();
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

        if (IsInterface(bcType) || bcType == PHENGLEI::OVERSET)
        {
            continue;
        }
        if (bcType == EXTRAPOLATION || bcType == SYMMETRY || bcType == OUTFLOW || bcType == PRESSURE_OUTLET || bcType == MASS_FLOW_OUTLET)
        {
            OutflowBCRegion(grid, bcRegion);
        }
        else if(bcType == INFLOW || bcType == PRESSURE_INLET || bcType == MASS_FLOW_INLET)
        {
            InflowBCRegion(grid, bcRegion);
        }
        else if(IsWall(bcType))
        {
            if (wallFunctionType == WALLFUNCTION::NONE)
            {
               WallBCRegion(grid, bcRegion);
            }
            else if (wallFunctionType == WALLFUNCTION::STANDARD)
            {
                VisWallWithWallFunctionStandard(grid, bcRegion);
            }
            else if (wallFunctionType == WALLFUNCTION::PAB3D)
            {
                VisWallWithWallFunctionPAB3D(grid, bcRegion);
            }
            else
            {
                TK_Exit::UnexpectedVarValue("wallFunctionType", wallFunctionType);
            }
        }
        else if (bcType == FARFIELD)
        {
            FarfieldBCRegion(grid, bcRegion);
        }
        else
        {
            TK_Exit::ExceptionExit("Error: this boundary type does not exist!\n");
        }
    }
}

void TurbSolverUnstr::OutflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        viscousTurbulence[re] = viscousTurbulence[le];
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qTurbulence[m][re] = qTurbulence[m][le];
        }
    }
}

//! Overlapping boundary : The flow field variables are interpolated from the interpolating cell to the invalid cell 
void TurbSolverUnstr::ReSetOversetBoundary(Grid* gridIn)
{
    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");
    if (isOversetSlip)
    {
        return;
    }
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    vector<int> *c2c = grid->GetCell2Cell();
    int *iBlank = grid->GetBlankIndex();

    //! Get the number of equations.
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));

    for (int iCell = 0; iCell < nTotalCell; iCell++)
    {
        if (iBlank[iCell] == INTERPOLATION)
        {
            int nNeighbor = static_cast<int>(c2c[iCell].size());
            for (int j = 0; j < nNeighbor; ++j)
            {
                int neighborCell = c2c[iCell][j];
                if (neighborCell < nTotalCell && iBlank[neighborCell] == INACTIVE)
                {
                    viscousTurbulence[neighborCell] = viscousTurbulence[iCell];
                    for (int m = 0; m < nTurbulenceEquation; ++m)
                    {
                        qTurbulence[m][neighborCell] = qTurbulence[m][iCell];
                    }
                }
            }
        }
    }
}

void TurbSolverUnstr::InflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();    
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int re = rightCellOfFace[iFace];

        viscousTurbulence[re] = freeStreamViscosity;
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qTurbulence[m][re] = freeStreamTurbVar[m];
        }
    }
}

void TurbSolverUnstr::WallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *wallDistance = grid->GetWallDist();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **q           = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar    = reinterpret_cast<RDouble  *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));

    using namespace IDX;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    if(tscheme != GMRES)
    {
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];

            viscousTurbulence[re] = -viscousTurbulence[le];
            if (viscousType == ONE_EQU)
            {
                qTurbulence[ISA][re] = - qTurbulence[ISA][le];
            }
            else if (viscousType == TWO_EQU)
            {
                RDouble beta1 = 0.075;
                RDouble rho   = q[IR][le];
                RDouble dist  = wallDistance[le];
                RDouble dist2 = dist * dist;

                RDouble kwWall = 60.0 * viscousLaminar[le] / (rho * beta1 * dist2 * refReNumber * refReNumber);

                qTurbulence[IKE][re] = -qTurbulence[IKE][le];
                //qTurbulence[IKW][re] = MAX(two * kw_wall - qTurbulence[IKW][le],SMALL);
                qTurbulence[IKW][re] = two * kwWall - qTurbulence[IKW][le];
            }
        }
    }
    else
    {
#ifdef USE_GMRESSOLVER
        //! GMRESturb Allocation
        int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
        int nTotalCell = grid->GetNTotalCell();
        ADReal *turbQl = new ADReal[nTurbulenceEquation]();
        ADReal *turbQr = new ADReal[nTurbulenceEquation]();
        RDouble **dDdP_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP_turb"));

        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];

            for (int n = 0; n < nTurbulenceEquation; n++)
            {
                turbQl[n] = qTurbulence[n][le];
                turbQl[n].diff(n, nTurbulenceEquation);
            }

            viscousTurbulence[re] = -viscousTurbulence[le];
            if (viscousType == ONE_EQU)
            {
                for (int n = 0; n < nTurbulenceEquation; n++)
                {
                    turbQr[n] = -turbQl[n];
                }

                int idx = (re - nTotalCell) * nTurbulenceEquation;        
                for (int m = 0; m < nTurbulenceEquation; m++)
                {
                    qTurbulence[m][re] = turbQr[m].val();
                    for (int n = 0; n < nTurbulenceEquation; n++)
                    {
                        dDdP_turb[m][idx + n] = turbQr[m].dx(n);
                    }
                }
            }
            else
            {
                TK_Exit::ExceptionExit("GMRES for 1eq-sa only so far");
            }
        }

        delete[] turbQl; turbQl = nullptr;
        delete[] turbQr; turbQr = nullptr;
#endif
    }
}

#ifdef USE_GMRESSOLVER
//! GMRESturb
void TurbSolverUnstr::GMRES_WallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *wallDistance = grid->GetWallDist();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    RDouble refReNumber = parameters->GetRefReNumber();
    int transitionType = parameters->GetTransitionType();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble  *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));

    //! GMRESturb Allocation
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int nTotalCell = grid->GetNTotalCell();
    ADReal *turbQl = new ADReal[nTurbulenceEquation]();
    ADReal *turbQr = new ADReal[nTurbulenceEquation]();
    RDouble **dDdP_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP_turb"));

    using namespace IDX;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        for (int n = 0; n < nTurbulenceEquation; n++)
        {
            turbQl[n] = qTurbulence[n][le];
            turbQl[n].diff(n, nTurbulenceEquation);
        }

        viscousTurbulence[re] = -viscousTurbulence[le];
        if (viscousType == ONE_EQU)
        {
            for (int n = 0; n < nTurbulenceEquation; n++)
            {
                turbQr[n] = -turbQl[n];
            }

            int idx = (re - nTotalCell) * nTurbulenceEquation;        
            for (int m = 0; m < nTurbulenceEquation; m++)
            {
                qTurbulence[m][re] = turbQr[m].val();
                for (int n = 0; n < nTurbulenceEquation; n++)
                {
                    dDdP_turb[m][idx + n] = turbQr[m].dx(n);
                }
            }
            /* printf("turbulence wall bc: %d\n", le);
            for (int m = 0; m < nTurbulenceEquation; m++)
            {
            for (int n = 0; n < nTurbulenceEquation; n++)
            {
            int index;
            index = (re - nTotalCell) * nTurbulenceEquation;
            printf("%lf  ", dDdP_turb[m][index + n]);
            }
            printf("\n");
            }
            printf("-------------------------------------------------\n"); */
        }
        else
        {
            TK_Exit::ExceptionExit("GMRES for 1eq-sa only so far");
        }
    }

    delete[] turbQl; turbQl = nullptr;
    delete[] turbQr; turbQr = nullptr;
}
#endif

void TurbSolverUnstr::VisWallWithWallFunctionStandard(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *faceNormalX = grid->GetFaceNormalX();
    RDouble *faceNormalY = grid->GetFaceNormalY();
    RDouble *faceNormalZ = grid->GetFaceNormalZ();
    RDouble *wallDistance = grid->GetWallDist();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble **qTurbulence = reinterpret_cast<RDouble**> (grid->GetDataPtr("q_turb"));
    RDouble **qLaminar = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble*> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble*> (grid->GetDataPtr("vist"));

    RDouble prandtlLaminar = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("prt");

    RDouble jayatillekeP = 9.24 * (pow((prandtlLaminar / prandtlTurbulence), 0.75) - 1.0)
        * (1.0 + 0.28 * exp(-0.007 * prandtlLaminar / prandtlTurbulence));
    RDouble refGama = parameters->GetRefGama();
    RDouble refMachNumber = parameters->GetRefMachNumber();

    const RDouble E_ = 9.793;
    const RDouble beta1 = 0.075;
    const RDouble yPluslimit = 11.225;
    const RDouble kappa1 = 0.4178;
    const RDouble CMu = 0.09;

    using namespace IDX;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        RDouble densityWall = qLaminar[IR][le];
        RDouble up = qLaminar[IU][le];
        RDouble vp = qLaminar[IV][le];
        RDouble wp = qLaminar[IW][le];
        RDouble laminarViscosity = viscousLaminar[le];
        RDouble turbulentViscosity = viscousTurbulence[le];

        RDouble xfn = faceNormalX[iFace];
        RDouble yfn = faceNormalY[iFace];
        RDouble zfn = faceNormalZ[iFace];

        RDouble vpn = up * xfn + vp * yfn + wp * zfn;
        RDouble vpt = sqrt(up * up + vp * vp + wp * wp - vpn * vpn);
        RDouble small = 1.0e-15;
        RDouble uvw = max(vpt, small);

        RDouble wallDist = wallDistance[le];

        //! Obtain the reynold number.
        RDouble refReNumber = parameters->GetRefReNumber();

        RDouble taow = (laminarViscosity + turbulentViscosity) * uvw / wallDist / refReNumber;
        RDouble utao = sqrt(taow / densityWall);

        RDouble yPlus = refReNumber * wallDist * densityWall * utao / laminarViscosity;
        RDouble epslon = 1.0;
        RDouble viscousTurbulenceWall;
        RDouble tempturePlus = 1.0;
        int iterStep = 0;
        while (epslon > 0.001 && ++iterStep < 10)
        {
            RDouble utemp = utao;

            if (yPlus > yPluslimit)
            {
                utao = kappa1 * uvw / (log(E_ * yPlus));
                tempturePlus = prandtlTurbulence * (log(E_ * yPlus) / kappa1 + jayatillekeP);
            }
            else
            {
                utao = uvw / yPlus;
                tempturePlus = prandtlLaminar * yPlus;
            }
            yPlus = refReNumber * wallDist * densityWall * utao / laminarViscosity;
            taow = densityWall * utao * utao;

            epslon = abs(utao - utemp) / (abs(utemp) + 1.0e-20);
        }

        viscousTurbulenceWall = laminarViscosity * (taow / (uvw + small) / laminarViscosity * wallDist * refReNumber - 1.0);
        viscousTurbulenceWall = MAX(viscousTurbulenceWall, 0.0);
        viscousTurbulence[re] = 2.0 * viscousTurbulenceWall - viscousTurbulence[le];

        if (viscousType == ONE_EQU)
        {
            qTurbulence[ISA][re] = -qTurbulence[ISA][le];
        }
        else if (viscousType == TWO_EQU)
        {
            RDouble distanceSquare = wallDist * wallDist;

            RDouble omgi = 6.0 * laminarViscosity / (densityWall * beta1 * distanceSquare * refReNumber * refReNumber);
            RDouble omgo = utao / (0.126 * wallDist * refReNumber);
            RDouble kwWall = sqrt(omgi * omgi + omgo * omgo);

            if (yPlus > 11.225)
            {
                qTurbulence[IKW][le] = kwWall;
                qTurbulence[IKW][re] = 20.0 * kwWall - qTurbulence[IKW][le];

                qTurbulence[IKE][le] = qTurbulence[IKW][le] * viscousTurbulence[le] / densityWall;
                qTurbulence[IKE][re] = -qTurbulence[IKE][le];
            }
            else
            {
                omgi = MAX(omgi, 0.0);
                qTurbulence[IKW][re] = 20.0 * omgi - qTurbulence[IKW][le];
                qTurbulence[IKE][re] = -qTurbulence[IKE][le];
            }
        }
    }
}

void TurbSolverUnstr::VisWallWithWallFunctionPAB3D(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *faceNormalX = grid->GetFaceNormalX();
    RDouble *faceNormalY = grid->GetFaceNormalY();
    RDouble *faceNormalZ = grid->GetFaceNormalZ();
    RDouble *wallDistance = grid->GetWallDist();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble **qTurbulence = reinterpret_cast<RDouble**> (grid->GetDataPtr("q_turb"));
    RDouble **qLaminar = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble*> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble*> (grid->GetDataPtr("vist"));
    RDouble **heatTransferCoeff = reinterpret_cast<RDouble**> (grid->GetDataPtr("heatTransferCoeff"));

    RDouble prandtlLaminar = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("prt");

    RDouble paba10[7] = { 2.354039, 0.1179840, -4.2899192e-04, 2.0404148e-06,-5.1775775e-09, 6.2687308e-12, -2.916958e-15 };
    RDouble paba11[5] = { 5.777191, 6.8756983e-02, -7.1582745e-06, 1.5594904e-09, -1.4865778e-13 };
    RDouble paba12[5] = { 31.08654, 5.0429072e-02, -2.0072314e-8 };
    RDouble beta1 = 0.075;

    using namespace IDX;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        RDouble rhoInside = qLaminar[IR][le];
        RDouble rhoGhost = qLaminar[IR][re];
        RDouble densityWall = half * (rhoInside + rhoGhost);
        RDouble up = qLaminar[IU][le];
        RDouble vp = qLaminar[IV][le];
        RDouble wp = qLaminar[IW][le];
        RDouble laminarViscosity = viscousLaminar[le];


        RDouble xfn = faceNormalX[iFace];
        RDouble yfn = faceNormalY[iFace];
        RDouble zfn = faceNormalZ[iFace];

        RDouble vpn = up * xfn + vp * yfn + wp * zfn;
        RDouble vpt = sqrt(up * up + vp * vp + wp * wp - vpn * vpn);
        RDouble small = 1.0e-15;
        RDouble uvw = max(vpt, small);

        RDouble wallDist = wallDistance[le];
        RDouble refReNumber = parameters->GetRefReNumber(); //! Obtain the reynold number.

        RDouble rc = rhoInside * refReNumber * wallDist * uvw / laminarViscosity;
        RDouble xnplus = 0.0;

        //! The original one here is yPlus.
        if (rc <= 20.24)
        {
            xnplus = sqrt(rc);
        }
        else if (rc <= 435.0)
        {
            xnplus = paba10[0] + paba10[1] * rc + paba10[2] * rc * rc + paba10[3] * rc * rc * rc + paba10[4] * rc * rc * rc * rc + paba10[5] * rc * rc * rc * rc * rc + paba10[6] * rc * rc * rc * rc * rc * rc;
        }
        else if (rc <= 4000.0)
        {
            xnplus = paba11[0] + paba11[1] * rc + paba11[2] * rc * rc + paba11[3] * rc * rc * rc + paba11[4] * rc * rc * rc * rc;
        }
        else
        {
            xnplus = paba12[0] + paba12[1] * rc + paba12[2] * rc * rc;
        }
        RDouble xnplussav = xnplus;

        //! Newton iteration to solve for nplus, assuming it is in log region:
        for (int iterStep = 0; ; ++iterStep)
        {
            RDouble f = rc / xnplus - 2.44 * (log(xnplus)) - 5.2;
            RDouble dfdn = -rc / (xnplus * xnplus) - 2.44 / xnplus;
            RDouble delta = -f / dfdn;
            xnplus = fabs(xnplus + delta);
            if (iterStep > 10)
            {
                //! Revert back to approx series soln if Newton iteration fails.
                xnplus = xnplussav;
                break;
            }
            else if (abs(delta) < 1.e-3)
            {
                break;
            }
        }

        RDouble viscousTurbulenceWall = laminarViscosity * (laminarViscosity * xnplus * xnplus / (wallDist * rhoInside * uvw * refReNumber) - 1.0);
        viscousTurbulenceWall = MAX(viscousTurbulenceWall, 0.0);
        viscousTurbulence[re] = 2.0 * viscousTurbulenceWall - viscousTurbulence[le];
        if (viscousType == ONE_EQU)
        {
            qTurbulence[ISA][re] = -qTurbulence[ISA][le];
        }
        else if (viscousType == TWO_EQU)
        {
            RDouble distanceSquare = wallDist * wallDist;

            RDouble omgi = 6.0 * laminarViscosity / (rhoInside * beta1 * distanceSquare * refReNumber * refReNumber);
            RDouble taow = laminarViscosity * uvw / wallDist / refReNumber;
            RDouble utao = sqrt(taow / rhoInside);
            RDouble omgo = utao / (0.126 * wallDist * refReNumber);
            RDouble kwWall = sqrt(omgi * omgi + omgo * omgo);

            qTurbulence[IKW][re] = 20.0 * kwWall - qTurbulence[IKW][le];

            qTurbulence[IKE][re] = qTurbulence[IKW][re] * viscousTurbulence[re] / rhoGhost;
        }
    }
}

void TurbSolverUnstr::FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *nxs = grid->GetFaceNormalX();
    RDouble *nys = grid->GetFaceNormalY();
    RDouble *nzs = grid->GetFaceNormalZ();
    RDouble *vgn = grid->GetFaceNormalVelocity();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble refGama = parameters->GetRefGama();
    RDouble *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **q           = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));
    RDouble *gamma = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));

    using namespace IDX;
    RDouble *primitiveVariableInflow = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble roo = primitiveVariableInflow[IR];
    RDouble uoo = primitiveVariableInflow[IU];
    RDouble voo = primitiveVariableInflow[IV];
    RDouble woo = primitiveVariableInflow[IW];
    RDouble poo = primitiveVariableInflow[IP];

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    if (tscheme != GMRES)
    {
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            //! iFace is the face number in the set of faceIndex.
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

        vno = nxs[iFace] * uoo + nys[iFace] * voo + nzs[iFace] * woo - vgn[iFace];
        vni = nxs[iFace] * uin + nys[iFace] * vin + nzs[iFace] * win - vgn[iFace];
            vei = sqrt(uin * uin + vin * vin + win * win);

            RDouble gama = gamma[le];

            coo = sqrt(ABS(refGama * poo / roo));
            cIN = sqrt(ABS(gama  * pin / rin));

            //! supersonic
            if (vei > cIN)
            {
                if (vni >= 0.0)
                {
                    //! exit
                    viscousTurbulence[re] = viscousTurbulence[le];
                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulence[m][re] = qTurbulence[m][le];
                    }
                }
                else
                {
                    //! inlet
                    viscousTurbulence[re] = freeStreamViscosity;
                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        qTurbulence[m][re] = freeStreamTurbVar[m];
                    }
                }
                continue;
            }

            //! subsonic
            riemp = vni + 2.0 * cIN / (gama  - 1.0);
            riemm = vno - 2.0 * coo / (refGama - 1.0);
            vnb   = half   * (riemp + riemm);
            cb    = fourth * (riemp - riemm) * gama;

            if (vnb >= 0.0)
            {
                //! exit
                viscousTurbulence[re] = viscousTurbulence[le];
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulence[m][re] = qTurbulence[m][le];
                }
            }
            else
            {
                //! inlet
                viscousTurbulence[re] = freeStreamViscosity;
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulence[m][re] = freeStreamTurbVar[m];
                }
            }

        }
    }
    else
    {
#ifdef USE_GMRESSOLVER
        //! GMRESturb
        int nTotalCell = grid->GetNTotalCell();
        ADReal *turbQl = new ADReal[nTurbulenceEquation]();
        ADReal *turbQr = new ADReal[nTurbulenceEquation]();
        RDouble **dDdP_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP_turb"));

        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];

            for (int m = 0; m < nTurbulenceEquation; m++)
            {
                turbQl[m] = qTurbulence[m][le];
                turbQl[m].diff(m, nTurbulenceEquation);
            }


            RDouble rin,uin,vin,win,pin;
            RDouble vno,vni,vei,coo,cin;
            RDouble riemp,riemm,cb,vnb;

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
            cin = sqrt(ABS(gama  * pin / rin));

            //! Supersonic
            if (vei > cin)
            {
                if (vni >= 0.0)
                {
                    //! Exit
                    viscousTurbulence[re] = viscousTurbulence[le];
                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        turbQr[m] = turbQl[m];
                    }
                }
                else
                {
                    //! Inlet
                    viscousTurbulence[re] = freeStreamViscosity;
                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        turbQr[m] = freeStreamTurbVar[m];
                    }
                }
                continue;
            }

            //! Subsonic
            riemp = vni + 2.0 * cin / (gama  - 1.0);
            riemm = vno - 2.0 * coo / (refGama - 1.0);
            vnb   = half   * (riemp + riemm);
            cb    = fourth * (riemp - riemm) * gama;

            if (vnb >= 0.0)
            {
                //! Exit
                viscousTurbulence[re] = viscousTurbulence[le];
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    turbQr[m] = turbQl[m];
                }
            }
            else
            {
                //! Inlet
                viscousTurbulence[re] = freeStreamViscosity;
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    turbQr[m] = freeStreamTurbVar[m];
                }
            }

            int idx = (re - nTotalCell) * nTurbulenceEquation;        
            for (int m = 0; m < nTurbulenceEquation; m++)
            {
                qTurbulence[m][re] = turbQr[m].val();
                for (int n = 0; n < nTurbulenceEquation; n++)
                {
                    dDdP_turb[m][idx + n] = turbQr[m].dx(n);
                }
            }
        }
#endif
    }
}

#ifdef USE_GMRESSOLVER
//! GMRESturb
void TurbSolverUnstr::GMRES_FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble  *nxs  = grid->GetFaceNormalX();
    RDouble  *nys  = grid->GetFaceNormalY();
    RDouble  *nzs  = grid->GetFaceNormalZ();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble refGama = parameters->GetRefGama();
    RDouble *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));
    RDouble *gamma = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));

    //! GMRESturb
    int nTotalCell = grid->GetNTotalCell();
    ADReal *turbQl = new ADReal[nTurbulenceEquation]();
    ADReal *turbQr = new ADReal[nTurbulenceEquation]();
    RDouble **dDdP_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP_turb"));

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
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        for (int m = 0; m < nTurbulenceEquation; m++)
        {
            turbQl[m] = qTurbulence[m][le];
            turbQl[m].diff(m, nTurbulenceEquation);
        }
        

        RDouble rin,uin,vin,win,pin;
        RDouble vno,vni,vei,coo,cin;
        RDouble riemp,riemm,cb,vnb;

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
        cin = sqrt(ABS(gama  * pin / rin));

        //! Supersonic
        if (vei > cin)
        {
            if (vni >= 0.0)
            {
                //! Exit
                viscousTurbulence[re] = viscousTurbulence[le];
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    turbQr[m] = turbQl[m];
                }
            }
            else
            {
                //! Inlet
                viscousTurbulence[re] = freeStreamViscosity;
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    turbQr[m] = freeStreamTurbVar[m];
                }
            }
            continue;
        }

        //! Subsonic
        riemp = vni + 2.0 * cin / (gama  - 1.0);
        riemm = vno - 2.0 * coo / (refGama - 1.0);
        vnb   = half   * (riemp + riemm);
        cb    = fourth * (riemp - riemm) * gama;

        if (vnb >= 0.0)
        {
            //! Exit
            viscousTurbulence[re] = viscousTurbulence[le];
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                turbQr[m] = turbQl[m];
            }
        }
        else
        {
            //! Inlet
            viscousTurbulence[re] = freeStreamViscosity;
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                turbQr[m] = freeStreamTurbVar[m];
            }
        }

        int idx = (re - nTotalCell) * nTurbulenceEquation;        
        for (int m = 0; m < nTurbulenceEquation; m++)
        {
            qTurbulence[m][re] = turbQr[m].val();
            for (int n = 0; n < nTurbulenceEquation; n++)
            {
                dDdP_turb[m][idx + n] = turbQr[m].dx(n);
            }
        }
        /* printf("turbulence far bc: %d\n", le);
        for (int m = 0; m < nTurbulenceEquation; m++)
        {
            for (int n = 0; n < nTurbulenceEquation; n++)
            {
                int index;
                index = (re - nTotalCell) * nTurbulenceEquation;
                printf("%lf  ", dDdP_turb[m][index + n]);
            }
            printf("\n");
        }
        printf("-------------------------------------------------\n"); */
    }
}
#endif
void TurbSolverUnstr::ComputeViscousCoeff(Grid *grid)
{
    Viscosity(grid);
}

void TurbSolverUnstr::Viscosity(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    RDouble *viscousLaminar = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));

    RDouble *SST_F2 = reinterpret_cast<RDouble *> (grid->GetDataPtr("SST_F2"));

    RDouble *wallDistance = grid->GetWallDist();

    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfBoundaryFace + numberOfTotalCell;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();
    int monitorViscousTurbulenceMax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");
    string viscousName  = parameters->GetViscousName();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble reynoldsSquare = refReNumber * refReNumber;

    RDouble viscousTurbulenceMax = 0.0;
    RDouble viscousTurbulenceMin = 1.0e30;

    RDouble part1, part2, arg2, f2;
    RDouble ld, ld3, fv1;
    RDouble rho, ke, kw, ds, ds2;

    int iMax = 0;
    int iMin = 0;

    using namespace IDX;

    if (viscousName.substr(0, 6) == "1eq-sa")
    {
        RDouble SA_cv1_cube = parameters->GetSA_cv1_cube();

        for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
        {
            RDouble nue = qTurbulence[0][iCell];
            if(nue > 0.0)
            {
                ld  = nue * q[IR][iCell] / (viscousLaminar[iCell] + SMALL);
                ld3 = ld * ld * ld;
                fv1 = ld3 / (ld3 + SA_cv1_cube);
                viscousTurbulence[iCell] = q[IR][iCell] * fv1 * nue;
            }
            else
            {
                viscousTurbulence[iCell] = 0.0;
            }

            if (iCell < numberOfTotalCell && viscousTurbulenceMax < viscousTurbulence[iCell])    //! Bell 20130527 add iCell <numberOfTotalCell
            {
                viscousTurbulenceMax = viscousTurbulence[iCell];
                iMax = iCell;
            }
            if (iCell < numberOfTotalCell && viscousTurbulenceMin > viscousTurbulence[iCell])    //! Bell 20130527 add iCell <numberOfTotalCell
            {
                viscousTurbulenceMin = viscousTurbulence[iCell];
                iMin = iCell;
            }
        }
    }
    else if (viscousName.substr(0, 6) == "2eq-kw")
    {
        if (viscousName.substr(0, 17) == "2eq-kw-menter-sst")
        {
            //! Original Menter k-w-SST model
            //! Important: no repalce! otherwise it can not repeat.
            RDouble SST_a1 = parameters->GetSST_a1();
            RDouble SST_betaStar = parameters->GetSST_betaStar();

            RDouble *dudx = new RDouble [numberOfTotal];
            RDouble *dudy = new RDouble [numberOfTotal];
            RDouble *dudz = new RDouble [numberOfTotal];

            RDouble *dvdx = new RDouble [numberOfTotal];
            RDouble *dvdy = new RDouble [numberOfTotal];
            RDouble *dvdz = new RDouble [numberOfTotal];

            RDouble *dwdx = new RDouble [numberOfTotal];
            RDouble *dwdy = new RDouble [numberOfTotal];
            RDouble *dwdz = new RDouble [numberOfTotal];

            grid->CompGradient(q[IU], dudx, dudy, dudz);
            grid->CompGradient(q[IV], dvdx, dvdy, dvdz);
            grid->CompGradient(q[IW], dwdx, dwdy, dwdz);

            RDouble *SpSdRatio = reinterpret_cast<RDouble *> (grid->GetDataPtr("SpSdRatio"));

            for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
            {
                rho = ABS(q[IR][iCell]) + SMALL;
                ke  = qTurbulence[IKE][iCell];
                kw  = qTurbulence[IKW][iCell];
                ds  = wallDistance[iCell];
                ds2 = ds * ds;

                RDouble s11 = dudx[iCell];
                RDouble s22 = dvdy[iCell];
                RDouble s33 = dwdz[iCell];
                RDouble s12 = half * (dudy[iCell] + dvdx[iCell]);
                RDouble s13 = half * (dudz[iCell] + dwdx[iCell]);
                RDouble s23 = half * (dvdz[iCell] + dwdy[iCell]);

                RDouble sij2 = two * (PHSPACE::SQR(s11, s22, s33) + two * PHSPACE::SQR(s12, s13, s23));

                part1 = 2.0 * sqrt(ke) / (SST_betaStar * kw * ds * refReNumber);
                part2 = 500.0 * viscousLaminar[iCell] / (rho * kw * ds2 * reynoldsSquare);
                arg2  = MAX(part1, part2);
                f2    = tanh(arg2 * arg2);
                SST_F2[iCell] = f2;

                //viscousTurbulence[iCell] = rho * ke / MAX(kw, vort * f2 / turbA1) * refReNumber;
                viscousTurbulence[iCell] = rho * ke / MAX(kw, sqrt(sij2) * f2 / (SST_a1 * refReNumber));

                int transitionType = parameters->GetTransitionType();
                if (transitionType == IREGAMA)
                {
                    viscousTurbulence[iCell] = viscousTurbulence[iCell] * sqrt(SpSdRatio[iCell]);
                }

                viscousTurbulence[iCell] = MIN(eddyViscosityLimit, viscousTurbulence[iCell]); //viscousLaminar[iCell] * coef_kvist * 

                if (viscousTurbulenceMax < viscousTurbulence[iCell]) 
                {
                    viscousTurbulenceMax = viscousTurbulence[iCell];
                    iMax = iCell;
                }
                if (viscousTurbulenceMin > viscousTurbulence[iCell])
                {
                    viscousTurbulenceMin = viscousTurbulence[iCell];
                    iMin = iCell;
                }
            }

            delete [] dudx;    dudx = nullptr;
            delete [] dudy;    dudy = nullptr;
            delete [] dudz;    dudz = nullptr;

            delete [] dvdx;    dvdx = nullptr;
            delete [] dvdy;    dvdy = nullptr;
            delete [] dvdz;    dvdz = nullptr;

            delete [] dwdx;    dwdx = nullptr;
            delete [] dwdy;    dwdy = nullptr;
            delete [] dwdz;    dwdz = nullptr;
        }
    }

    RDouble vistMax = viscousTurbulenceMax;
    RDouble vistMin = viscousTurbulenceMin;

    grid->UpdateData("vist_max", &vistMax, PHDOUBLE, 1);
    grid->UpdateData("vist_min", &vistMin, PHDOUBLE, 1);

    if (monitorViscousTurbulenceMax)
    {
        int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
        int intervalStepRes = GlobalDataBase::GetIntParaFromDB("intervalStepRes");
        if (outnstep % intervalStepRes == 0)
        {
        grid->CompareMaxMinValue(viscousTurbulenceMax, 1);
        grid->CompareMaxMinValue(viscousTurbulenceMin, 2);

        ostringstream oss;
        cout << setprecision(6);
        cout << setiosflags(ios::scientific);
        oss  << "Max && min viscosity: " << viscousTurbulenceMax << " " << viscousTurbulenceMin << " " << "\n";
        PrintToWindow(oss);
        }
    }
}

void TurbSolverUnstr::InitCGrid(Grid *fineGridIn, Grid *coarseGridIn)
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

    for (int iCell = 0; iCell <numberOfTotalCellOnCoarseGrid; ++ iCell)
    {
        walldistanceOfCoarseGrid[iCell] = 0.0;
    }

    for (int iCell = 0; iCell <numberOfTotalCellOnFineGrid; ++ iCell)
    {
        walldistanceOfCoarseGrid[cell2CoarseGridCell[iCell]] += volumeOfFineGrid[iCell] * walldistanceOfFineGrid[iCell];
    }

    for (int iCell = 0; iCell <numberOfTotalCellOnCoarseGrid; ++ iCell)
    {
        walldistanceOfCoarseGrid[iCell] /= volumeOfCoarseGrid[iCell];
    }
}

void TurbSolverUnstr::UploadInterfaceData(ActionKey *actkey)
{
    UploadInterfaceValue(actkey);
}

void TurbSolverUnstr::UploadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (interfaceInformation == 0) return;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    PHSPACE::UploadInterfaceValue(grid, qTurbulence, "turb::q", nTurbulenceEquation);

    RDouble *viscousTurbulence = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    PHSPACE::UploadInterfaceValue(grid, viscousTurbulence, "vist");
}

void TurbSolverUnstr::UploadInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterpointInformation *interPointInformation = grid->GetInterpointInfo();
    if (interPointInformation == 0)
    {
        return;
    }

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTurbNode"));
    PHSPACE::UploadInterpointValue(grid, qTurbulenceNode, "qTurbNode", nTurbulenceEquation);
}

void TurbSolverUnstr::DownloadInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterpointInformation *interPointInformation = grid->GetInterpointInfo();
    if (interPointInformation == 0)
    {
        return;
    }

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTurbNode"));
    PHSPACE::DownloadInterpointValue(grid, qTurbulenceNode, "qTurbNode", nTurbulenceEquation);
}

void TurbSolverUnstr::CommunicationInterpointWeight(ActionKey *actkey)
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
                int *nodeValueSliceTurb = reinterpret_cast<int *>(grid->GetDataPtr("nodeValueSliceTurb"));

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
                    valueWeight = nodeValueSliceTurb[globalPoint];

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
                int *nodeValueSliceTurbTrade = reinterpret_cast<int *>(gridNeighbor->GetDataPtr("nodeValueSliceTurbTrade"));

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

                    nodeValueSliceTurbTrade[globalPoint] += valueWeight;
                }
            }

            delete cdata;    cdata = nullptr;
        }
    }
}

void TurbSolverUnstr::DownloadInterpointWeight(ActionKey *actkey)
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
        int *nodeValueSliceTurb = reinterpret_cast<int *>(grid->GetDataPtr("nodeValueSliceTurb"));
        int *nodeValueSliceTurbTrade = reinterpret_cast<int *>(grid->GetDataPtr("nodeValueSliceTurbTrade"));

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            nodeValueSliceTurb[iNode] += nodeValueSliceTurbTrade[iNode];
        }

        PHSPACE::SetField(nodeValueSliceTurbTrade, 0, nTotalNode);
    }
}

//! Bell 20120910 add
void TurbSolverUnstr::DownloadInterfaceData(ActionKey *actkey)
{
    DownloadInterfaceValue(actkey);
}
//

void TurbSolverUnstr::DownloadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (interfaceInformation == 0) return;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    PHSPACE::DownloadInterfaceValue(grid, qTurbulence, "turb::q", nTurbulenceEquation); 

    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    PHSPACE::DownloadInterfaceValue(grid, viscousTurbulence, "vist");
}

//! need further compllished!!!
void TurbSolverUnstr::UploadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int isSolve = GlobalDataBase::GetIntParaFromDB("isSolve");
    if (!isSolve)
    {
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    PHSPACE::UploadOversetValue(grid, qTurbulence, "turb::q", nTurbulenceEquation);

        RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
        PHSPACE::UploadOversetValue(grid, viscousTurbulence, "vist");
}
    else
    {
        RDouble **gradientTurbulenceX = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTurbulenceX"));
        RDouble **gradientTurbulenceY = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTurbulenceY"));
        RDouble **gradientTurbulenceZ = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTurbulenceZ"));

        PHSPACE::UploadOversetValue(grid, gradientTurbulenceX, "turb::dqdx", nTurbulenceEquation);
        PHSPACE::UploadOversetValue(grid, gradientTurbulenceY, "turb::dqdy", nTurbulenceEquation);
        PHSPACE::UploadOversetValue(grid, gradientTurbulenceZ, "turb::dqdz", nTurbulenceEquation);
    }
}

//! need further compllished!!!
void TurbSolverUnstr::DownloadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int isSolve = GlobalDataBase::GetIntParaFromDB("isSolve");
    if (!isSolve)
    {
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    PHSPACE::DownloadOversetValue(grid, qTurbulence, "turb::q", nTurbulenceEquation);

        RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
        PHSPACE::DownloadOversetValue(grid, viscousTurbulence, "vist");
}
    else
    {
        RDouble **gradientTurbulenceX = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTurbulenceX"));
        RDouble **gradientTurbulenceY = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTurbulenceY"));
        RDouble **gradientTurbulenceZ = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTurbulenceZ"));

        PHSPACE::DownloadOversetValue(grid, gradientTurbulenceX, "turb::dqdx", nTurbulenceEquation);
        PHSPACE::DownloadOversetValue(grid, gradientTurbulenceY, "turb::dqdy", nTurbulenceEquation);
        PHSPACE::DownloadOversetValue(grid, gradientTurbulenceZ, "turb::dqdz", nTurbulenceEquation);
    }
}

//! Bell 20120910 add
void TurbSolverUnstr::DumpInterfaceInfoBell()
{
    int level = 0;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));
    
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation) return;

    fstream file;
    std::ostringstream oss;
    oss << AbsoluteFilename("interface") << grid->GetZoneID() << ".dat";
    file.open(oss.str().c_str(), ios_base::out);
    file.setf(ios::scientific);
    file.precision(8);
    int seprate = 20;

    int numberOfInterFace = interfaceInformation->GetNIFace();
    if (numberOfInterFace == 0) return;

    int *interFace2ZoneID       = interfaceInformation->GetInterFace2ZoneID();
    int *interFace2InterFaceID  = interfaceInformation->GetInterFace2InterFaceID();
    int *interface2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();

    int numberOfBoundaryFace = grid->GetNBoundFace();

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));

    RDouble **Q1;
    RDouble **Q2;

    RDouble **visLaminar = new RDouble *[1];
    RDouble **VisTurbulence = new RDouble *[1];
    visLaminar[0] = viscousLaminar;
    VisTurbulence[0] = viscousTurbulence;
    Q1 = visLaminar;
    Q2 = VisTurbulence;

    file << " numberOfInterFace = " << numberOfInterFace << " numberOfBoundaryFace = " << numberOfBoundaryFace << "numberOfTotalCell = " << grid->GetNTotalCell() << "\n";
    file << "i" << setw(seprate) << "nb_zone" << setw(seprate) << "le" << setw(seprate) << "re" << setw(seprate) << "tle" << setw(seprate)  << "tre" << setw(seprate)
         << "qTurbulence[L]" << setw(seprate) << "qTurbulence[R]"
         << endl;
    
    for (int iFace = 0; iFace < numberOfInterFace; ++ iFace)
    {
        int targetZone = interFace2ZoneID[iFace];
        int targetInterface = interFace2InterFaceID[iFace];
        UnstructGrid *targetGrid = UnstructGridCast(PHSPACE::GetGrid(targetZone,0));

        if(!targetGrid) return;

        InterfaceInfo *targetInterfaceInformation = targetGrid->GetInterfaceInfo();
       
        int *targetInterface2BoundaryFace = targetInterfaceInformation->GetInterFace2BoundaryFace();

        int *targetLeftCellOfFace = targetGrid->GetLeftCellOfFace();
        int *targetRightCellOfFace = targetGrid->GetRightCellOfFace();

        int le = leftCellOfFace [interface2BoundaryFace[iFace]];
        int re = rightCellOfFace[interface2BoundaryFace[iFace]];

        int tle = targetLeftCellOfFace [targetInterface2BoundaryFace[targetInterface]];
        int tre = targetRightCellOfFace[targetInterface2BoundaryFace[targetInterface]];

        file << iFace << setw(seprate) << targetZone << setw(seprate) << le << setw(seprate) << re << setw(seprate) << tle << setw(seprate) << tre << setw(seprate)
             << Q1[0][le]<< setw(seprate) << Q1[0][re] << setw(seprate)
             << Q2[0][le]<< setw(seprate) << Q2[0][re]
             << endl;

    }
    delete [] visLaminar[1];    visLaminar[1] = NULL;
    delete [] VisTurbulence[1];    VisTurbulence[1] = NULL;
    file.close();
    file.clear();
}
//*/

void TurbSolverUnstr::ComputeLengthScaleofOneEquationModel(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();

    if (DESType == DES)
    {
        ComputeDESLength(gridIn);
    }
    else if (DESType == DDES)
    {
        ComputeDDESLength(gridIn);
    }
    else if (DESType == IDDES)
    {
        ComputeIDDESLength(gridIn);
    }
}

void TurbSolverUnstr::ComputeDESLength(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    RDouble *wallDistance = grid->GetWallDist();

    RDouble *DESLength = this->GetDESLength();
    if(!DESLength)
    {
        DESLength = new RDouble [numberOfTotalCell];
    }
    
    RDouble *lengthOfLES = new RDouble [numberOfTotalCell];
    RDouble *lowReynoldsNumberCorrection = new RDouble [numberOfTotalCell];

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();
    ComputeLESLengthofSA(gridIn, lengthOfLES, lowReynoldsNumberCorrection, DESType);

    RDouble lengthRANS, lengthLES;
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
    {
        lengthRANS = wallDistance[iCell];
        lengthLES  = lengthOfLES[iCell];

        DESLength[iCell] = MIN(lengthRANS, lengthLES);
    }

    this->SetDESLength(DESLength);
    
    delete []lengthOfLES;
    delete []lowReynoldsNumberCorrection;
}

void TurbSolverUnstr::ComputeDDESLength(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    RDouble *wallDist = grid->GetWallDist();

    RDouble *DDESLength = this->GetDESLength();
    if(!DDESLength)
    {
        DDESLength = new RDouble [numberOfTotalCell];
    }

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();

    RDouble *lengthOfLES = new RDouble [numberOfTotalCell];
    RDouble *lowReynoldsNumberCorrection = new RDouble [numberOfTotalCell];
    ComputeLESLengthofSA(gridIn, lengthOfLES, lowReynoldsNumberCorrection, DESType);

    RDouble **q     = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar  = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));

    RDouble **gradientPrimtiveVariableX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradientPrimtiveVariableY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradientPrimtiveVariableZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble *dudx = gradientPrimtiveVariableX[1];
    RDouble *dudy = gradientPrimtiveVariableY[1];
    RDouble *dudz = gradientPrimtiveVariableZ[1];

    RDouble *dvdx = gradientPrimtiveVariableX[2];
    RDouble *dvdy = gradientPrimtiveVariableY[2];
    RDouble *dvdz = gradientPrimtiveVariableZ[2];

    RDouble *dwdx = gradientPrimtiveVariableX[3];
    RDouble *dwdy = gradientPrimtiveVariableY[3];
    RDouble *dwdz = gradientPrimtiveVariableZ[3];

    RDouble oRefReNumber = parameters->GetoRefReNumber();
    const RDouble karm   = 0.41;
    const RDouble kap2   = karm * karm;

    using namespace IDX;
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
    {
        RDouble laminarViscosity = viscousLaminar[iCell];
        RDouble turbulentViscosity = viscousTurbulence[iCell];
        RDouble density = q[IR][iCell];

        RDouble wallDistance  = wallDist[iCell];
        RDouble wallDistance2 = SQR(wallDistance);
        RDouble lengthRans    = wallDistance;
        RDouble lengthLes     = lengthOfLES[iCell];

        RDouble sumOfGradient2 = SQR(dudx[iCell], dudy[iCell], dudz[iCell]) + 
            SQR(dvdx[iCell], dvdy[iCell], dvdz[iCell]) + 
            SQR(dwdx[iCell], dwdy[iCell], dwdz[iCell]);

        RDouble nueff = oRefReNumber * (laminarViscosity + turbulentViscosity) / density;
        RDouble rd    = nueff / (kap2 * wallDistance2 * MAX(sqrt(sumOfGradient2), 1.0e-20));
        RDouble fd    = 1.0 - tanh(POWER3(8.0 * rd));

        DDESLength[iCell] = lengthRans - fd * MAX(0.0, (lengthRans - lengthLes));
    }

    this->SetDESLength(DDESLength);

    delete [] lengthOfLES;
    delete [] lowReynoldsNumberCorrection;
}

void TurbSolverUnstr::ComputeIDDESLength(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    RDouble *wallDist = grid->GetWallDist();
    RDouble *largestLocalGridLength = grid->GetLargestLocalGridLength();

    RDouble *IDDESLength = this->GetDESLength();
    if(!IDDESLength)
    {
        IDDESLength = new RDouble [numberOfTotalCell];
    }

    RDouble *lengthOfLES = new RDouble [numberOfTotalCell];
    RDouble *lowReynoldsNumberCorrection = new RDouble [numberOfTotalCell];

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();
    ComputeLESLengthofSA(gridIn, lengthOfLES, lowReynoldsNumberCorrection, DESType);

    RDouble **q     = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar  = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));

    RDouble **gradientPrimtiveVariableX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradientPrimtiveVariableY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradientPrimtiveVariableZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble *dudx = gradientPrimtiveVariableX[1];
    RDouble *dudy = gradientPrimtiveVariableY[1];
    RDouble *dudz = gradientPrimtiveVariableZ[1];

    RDouble *dvdx = gradientPrimtiveVariableX[2];
    RDouble *dvdy = gradientPrimtiveVariableY[2];
    RDouble *dvdz = gradientPrimtiveVariableZ[2];

    RDouble *dwdx = gradientPrimtiveVariableX[3];
    RDouble *dwdy = gradientPrimtiveVariableY[3];
    RDouble *dwdz = gradientPrimtiveVariableZ[3];

    RDouble oRefReNumber = parameters->GetoRefReNumber();
    const RDouble karm = 0.41;
    const RDouble kap2 = karm * karm;

    using namespace IDX;
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
    {
        RDouble laminarViscosity        = viscousLaminar[iCell];
        RDouble turbulentViscosity      = viscousTurbulence[iCell];
        RDouble density                 = q[IR][iCell];

        RDouble wallDistance            = wallDist[iCell];
        RDouble wallDistance2           = SQR(wallDistance);
        RDouble lengthRans              = wallDistance;
        RDouble lengthLes               = lengthOfLES[iCell];
        RDouble largestLocalGridSpacing = largestLocalGridLength[iCell];

        RDouble alf  = 0.25 - wallDistance / largestLocalGridSpacing;
        RDouble fb   = MIN(2.0 * exp(- 9.0 * alf * alf),  1.0);

        RDouble sumOfGradient2 = SQR(dudx[iCell], dudy[iCell], dudz[iCell]) + 
                                SQR(dvdx[iCell], dvdy[iCell], dvdz[iCell]) + 
                                SQR(dwdx[iCell], dwdy[iCell], dwdz[iCell]);

        RDouble nuet = oRefReNumber * turbulentViscosity / density;
        RDouble rdt  = nuet / (kap2 * wallDistance2 * MAX(sqrt(sumOfGradient2), 1.0e-10));

        RDouble fdt  = 1.0 - tanh(POWER3(8.0 * rdt));
        RDouble fdb  = MAX((1.0 - fdt), fb);

        RDouble ct   = 1.63;
        RDouble cL   = 3.55;
        RDouble nuel = oRefReNumber * laminarViscosity / density;
        RDouble rdl  = nuel / (kap2 * wallDistance2 * MAX(sqrt(sumOfGradient2), 1.0e-10));

        RDouble fl   = tanh(pow(cL * cL * rdl, 10.0));
        RDouble ft   = tanh(POWER3(ct * ct * rdt));

        RDouble fe2  = 1.0 - MAX(ft, fl);

        RDouble fe1 = 0.0;
        if (alf < 0.0)
        {
            fe1 = 2.0 * exp(- 9.0 * alf * alf);
        }
        else
        {
            fe1 = 2.0 * exp(- 11.09 * alf * alf);
        }

        RDouble fe = fe2 * MAX((fe1 - 1.0), 0.0) * lowReynoldsNumberCorrection[iCell];

        IDDESLength[iCell]  =  fdb * (1.0 + fe) * lengthRans  +  (1.0 - fdb) * lengthLes;
    }

    this->SetDESLength(IDDESLength);

    delete []lengthOfLES;
    delete []lowReynoldsNumberCorrection;
}

RDouble * TurbSolverUnstr::GetLengthScale(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();

    if (DESType != 0)
    {
        return this->GetDESLength();
    }
    else
    {
        return grid->GetWallDist();
    }
}

void TurbSolverUnstr::ReadStatisticalFlow(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **qAverageTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    int numberOfStatisticalStep = 0;
    cdata->Read(&numberOfStatisticalStep, sizeof(int));

    GlobalDataBase::UpdateData("nStatisticalStep", &numberOfStatisticalStep, PHINT, 1);

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            cdata->Read(&qAverageTurbulence[m][iCell], sizeof(RDouble));
        }
    }
}

void ComputeLESLengthofSA(Grid *gridIn, RDouble *lengthOfLES, RDouble *lowReynoldsNumberCorrection, int DESType)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    //! Compute LES length.
    RDouble CoeffOfDES = 0.65;

    ComputeLowReynoldsNumberCorrection(gridIn, lowReynoldsNumberCorrection);

    if (DESType == IDDES)
    {
        RDouble *subgridLength = grid->GetSubgridLength();

        for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
        {
            lengthOfLES[iCell] = CoeffOfDES * lowReynoldsNumberCorrection[iCell] * subgridLength[iCell];
        }
    }
    else
    {
        RDouble *largestLocalGridLength = grid->GetLargestLocalGridLength();

        for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
        {
            lengthOfLES[iCell] = CoeffOfDES * lowReynoldsNumberCorrection[iCell] * largestLocalGridLength[iCell];
        }
    }
}

void ComputeLowReynoldsNumberCorrection(Grid *gridIn, RDouble *lowReynoldsNumberCorrection)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    //! Compute low Reynolds number correction.
    RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *viscousLaminar   = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));

    const RDouble karm   = 0.41;
    const RDouble cb1    = 0.1355;
    const RDouble cb2    = 0.622;
    const RDouble cv1    = 7.1;
    const RDouble sigma  = 2.0/3.0;
    const RDouble cv13   =  cv1 * cv1 * cv1;
    const RDouble kap2   = karm * karm;
    const RDouble rkap2  = one / kap2;
    const RDouble cw1    = cb1 * rkap2 + (one + cb2) / sigma;
    const RDouble turbulentFwStar          = 0.424;

    using namespace IDX;
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
    {
        RDouble nuet = qTurbulence[ISA][iCell];
        RDouble olam = q[IR][iCell] / (viscousLaminar[iCell] + SMALL);
        RDouble xsi  = nuet * olam + SMALL;
        RDouble xsi3 = POWER3(xsi);

        RDouble turbulentFt2 = 0.0; 

        RDouble turbulentFv1 = xsi3 / (xsi3 + cv13);
        RDouble turbulentFv2 = one - xsi / (one + xsi * turbulentFv1);

        RDouble numerator   =  1.0  -  (cb1 / (cw1 * kap2 * turbulentFwStar)) * (turbulentFt2 + (1.0 - turbulentFt2) * turbulentFv2);
        RDouble denominator =  turbulentFv1 * MAX(1.0e-10, (1.0 - turbulentFt2));

        lowReynoldsNumberCorrection[iCell] = sqrt(MIN(1.0e2, numerator/denominator));
    }
}

void TurbSolverUnstr::ComputeGradient(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    if(grid->IsFinestGrid())
    {
        GetGradientField(gridIn);
    }
}

void TurbSolverUnstr::InitMixingPlane(RDouble ***MixingPlaneVar, int Dim1, int Dim2, int Dim3, RDouble Value)
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
void TurbSolverUnstr::AverageMixingPlane(Grid *grid)
{
    UnstructGrid *gridUnstruct   = UnstructGridCast(grid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1   = refGama - 1;

    RDouble refDensity  = parameters->GetRefDimensionalDensity();
    RDouble refPressure = parameters->GetRefDimensionalPressure();

    RDouble **q           = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q_turb"));

    RDouble **r_ori    = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Face"));
    RDouble **r_target = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Span"));

    RDouble **SpanArea = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("SpanArea"));

    //! nSpan is the number of spanwise point to do data average
    //! also interpolate data on these points.
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

    //! r_ori should be an array that contains coordinate of each node.
    //! the range of this array is equal to the face number of this boundary.
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

    RDouble ***SpanFlux     = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanFlux"));
    RDouble ***SpanTurbFlux = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanTurbFlux"));
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

        RDouble rotationAngle = PeriodicRotationAngle[iTurboZone] / 180.0 * PI;

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
                //! also need to compute flux of each spanwise.
                //! spanwise can be defined based on mesh or manually.
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
                RDouble Flux[5]     = { 0.0, 0.0, 0.0, 0.0, 0.0 };
                RDouble TurbFlux[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
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
                            TurbFlux[m] += q[IR][le] * vni * qTurbulence[m][le] * ns[iFace];
                        }
                    }
                }

                SpanArea[MixingFlag][iSpan] = TotalArea;

                for (int m = 0; m < nEquation; m++)
                {
                    SpanTurbFlux[MixingFlag][m][iSpan] = TurbFlux[m] / SpanArea[MixingFlag][iSpan];
                }
            }
        }
    }
}

//! SpanFlux from AverageMixingPlane should be input data.
//! Here we will obtain averaged primitive variable of each span from SpanFlux.
//! The averaged primitive variable will be treated as outlet data of upstream and inlet data of downsteam.
void TurbSolverUnstr::MixingPlaneDataTransfer(Grid *grid, Grid *NeighborGrid)
{
    UnstructGrid *gridUnstruct   = UnstructGridCast(grid);
    UnstructGrid *gridSource     = UnstructGridCast(NeighborGrid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1   = refGama - 1;

    RDouble refDensity  = parameters->GetRefDimensionalDensity();
    RDouble refPressure = parameters->GetRefDimensionalPressure();

    RDouble **q           = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast< RDouble **> (gridUnstruct->GetDataPtr("q_turb"));

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
    RDouble ***qTurb_Span    = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTurb_Span"));
    RDouble ***SourceFlux    = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanFlux"));
    RDouble ***SpanTurbFlux  = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanTurbFlux"));

    RDouble ***qTurbSpanNeighbor = reinterpret_cast< RDouble *** > (gridSource->GetDataPtr("qTurbSpanNeighbor"));

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
                    qTurb_Span[MixingFlag][m][iSpan] = SpanTurbFlux[MixingFlag][m][iSpan] / SourceFlux[MixingFlag][IR][iSpan];
                }

                for (int m = 0; m < nEquation; m++)
                {
                    qTurbSpanNeighbor[TargetFlag][m][iSpan] = qTurb_Span[MixingFlag][m][iSpan];
                }
            }
        }
    }
}

void TurbSolverUnstr::NonReflective(Grid *grid, Grid *NeighborGrid)
{
    UnstructGrid *gridUnstruct   = UnstructGridCast(grid);
    UnstructGrid *gridTarget     = UnstructGridCast(NeighborGrid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1   = refGama - 1;

    RDouble refDensity  = parameters->GetRefDimensionalDensity();
    RDouble refPressure = parameters->GetRefDimensionalPressure();

    RDouble **q           = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast< RDouble **> (gridUnstruct->GetDataPtr("q_turb"));

    RDouble **r_ori    = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Face"));
    RDouble **r_target = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Span"));

    RDouble ***qTurb_Span = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTurb_Span"));

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
    RDouble ***q_SpanCurrent     = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("q_Span"));
    RDouble ***q_SpanTarget      = reinterpret_cast< RDouble *** > (gridTarget->GetDataPtr("q_Span"));

    RDouble ***qTurb_SpanCurrent = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTurb_Span"));
    RDouble ***qTurbSpanNeighbor = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTurbSpanNeighbor"));
    RDouble ***SpanTurbFlux      = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanTurbFlux"));

    RDouble ***qTurb_SpanTarget = reinterpret_cast< RDouble *** > (gridTarget->GetDataPtr("qTurb_Span"));

    RDouble ***dqTurbSpanIn = reinterpret_cast <RDouble ***> (gridUnstruct->GetDataPtr("dqTurbSpanIn"));
    RDouble ***dqTurbSpanEx = reinterpret_cast <RDouble ***> (gridUnstruct->GetDataPtr("dqTurbSpanEx"));
    RDouble ***dcTurbSpan   = reinterpret_cast <RDouble ***> (gridUnstruct->GetDataPtr("dcTurbSpan"));

    using namespace IDX;

    int FlowDirection = -1;

    int nEquation = GetNumberOfEquations();

    RDouble *ExtraDqTurb = new RDouble[nEquation];
    RDouble *deltaQTurb  = new RDouble[nEquation];
    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion++)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        int bcType = bcRegion->GetBCType();
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
                        ExtraDqTurb[m] = 0.0;
                    }

                    for (int m = 0; m < nEquation; m++)
                    {
                        dqTurbSpanIn[MixingFlag][m][iSpan] = qTurb_SpanTarget[TargetFlag][m][iSpan] - qTurb_SpanCurrent[MixingFlag][m][iSpan];
                    }

                    //! extrapolation
                    vector<int> *faceIndex = bcRegion->GetFaceIndex();
                    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                    {
                        int iFace = *iter;

                        int re = rightCellofFace[iFace];
                        int le = leftCellofFace[iFace];

                        RDouble r_Face = sqrt(yfc[iFace] * yfc[iFace] + zfc[iFace] * zfc[iFace]);

                        if (r_Face > r_target[MixingFlag][iSpan] && r_Face <= r_target[MixingFlag][iSpan + 1])
                        {
                            for (int m = 0; m < nEquation; m++)
                            {
                                ExtraDqTurb[m] = (qTurb_Span[MixingFlag][m][iSpan] - qTurbulence[m][le]) * ns[iFace];
                            }
                        }
                    }

                    for (int m = 0; m < nEquation; m++)
                    {
                        dqTurbSpanEx[MixingFlag][m][iSpan] = ExtraDqTurb[m] / SpanArea[MixingFlag][iSpan];
                    }
                }
                for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
                {
                    //! after data transfer in NSSlover, IU, IV, IW refers to vx, vtheta, vr.
                    RDouble vn = nxsSpan[MixingFlag][iSpan] * q_SpanCurrent[MixingFlag][IU][iSpan]
                        + nrsSpan[MixingFlag][iSpan] * q_SpanCurrent[MixingFlag][IW][iSpan];

                    RDouble avgRho = half * (q_SpanCurrent[MixingFlag][IR][iSpan] + q_SpanTarget[TargetFlag][IR][iSpan]);
                    RDouble avgP   = half * (q_SpanCurrent[MixingFlag][IP][iSpan] + q_SpanTarget[TargetFlag][IP][iSpan]);

                    RDouble c = sqrt(refGama * (avgP / avgRho));

                    RDouble rhoc = half * (q_SpanCurrent[MixingFlag][IR][iSpan] + q_SpanTarget[TargetFlag][IR][iSpan]) * c;

                    if (vn > 0)
                    {
                        for (int m = 0; m < nEquation; m++)
                        {
                            dcTurbSpan[MixingFlag][m][iSpan] = rhoc * dqTurbSpanEx[MixingFlag][m][iSpan];
                        }
                    }
                    else
                    {
                        for (int m = 0; m < nEquation; m++)
                        {
                            dcTurbSpan[MixingFlag][m][iSpan] = rhoc * dqTurbSpanIn[MixingFlag][m][iSpan];
                        }
                    }

                    for (int m = 0; m < nEquation; m++)
                    {
                        deltaQTurb[m] = dcTurbSpan[MixingFlag][m][iSpan] / rhoc;
                    }

                    for (int m = 0; m < nEquation; m++)
                    {
                        qTurbSpanNeighbor[MixingFlag][m][iSpan] += deltaQTurb[m];
                    }
                }
            }
        }
    }
    delete[]ExtraDqTurb;
    delete[]deltaQTurb;
}

//! function to set radial profile of mixingin and mixingout.
//! treat mixingin as a kind of inlet and treat mixing as a kind of outlet.
void TurbSolverUnstr::SetMixingPlaneData(Grid *grid)
{
    UnstructGrid *gridUnstruct   = UnstructGridCast(grid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble refPressure = parameters->GetRefDimensionalPressure();
    RDouble refGama     = parameters->GetRefGama();
    RDouble gama1       = refGama - 1.0;

    RDouble **q           = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (gridUnstruct->GetDataPtr("q_turb"));

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

    int *leftCellofFace  = gridUnstruct->GetLeftCellOfFace();
    int *rightCellofFace = gridUnstruct->GetRightCellOfFace();
    int nTotalCell = gridUnstruct->GetNTotalCell();

    int nBCRegion = unstructBCSet->GetnBCRegion();

    int nEquation = GetNumberOfEquations();

    using namespace IDX;

    RDouble ***qTurb_Span    = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTurb_Span"));
    RDouble ***qSpanNeighbor = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qTurbSpanNeighbor"));

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
                    qTurb_Span[MixingFlag][m][iSpan] = qSpanNeighbor[MixingFlag][m][iSpan];
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
                            qTurbulence[m][re] = qTurb_Span[MixingFlag][m][iSpan];
                        }
                    }
                }
            }
        }
    }
}

void TurbSolverUnstr::RotateVectorFromInterface(Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    RDouble **rotTurbgradValueX = reinterpret_cast <RDouble **> (grid->GetDataPtr("rotTurbgradValueX"));
    RDouble **rotTurbgradValueY = reinterpret_cast <RDouble **> (grid->GetDataPtr("rotTurbgradValueY"));
    RDouble **rotTurbgradValueZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("rotTurbgradValueZ"));

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
        RDouble **fieldRecvTurbY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceY"));
        RDouble **fieldRecvTurbZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceZ"));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int t1;
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
                            rotTurbgradValueY[m][t1] = fieldRecvTurbY[m][t1] * cos(2 * PI - rotationAngle) - fieldRecvTurbZ[m][t1] * sin(2 * PI - rotationAngle);
                            rotTurbgradValueZ[m][t1] = fieldRecvTurbY[m][t1] * sin(2 * PI - rotationAngle) + fieldRecvTurbZ[m][t1] * cos(2 * PI - rotationAngle);
                        }

                        for (int m = 0; m < nEquation; ++m)
                        {
                            fieldRecvTurbY[m][t1] = rotTurbgradValueY[m][t1];
                            fieldRecvTurbZ[m][t1] = rotTurbgradValueZ[m][t1];
                        }
                    }
                    else if (bcName == Periodic_Name[2 * iTurboZone + 1])
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {
                            rotTurbgradValueY[m][t1] = fieldRecvTurbY[m][t1] * cos(rotationAngle) - fieldRecvTurbZ[m][t1] * sin(rotationAngle);
                            rotTurbgradValueZ[m][t1] = fieldRecvTurbY[m][t1] * sin(rotationAngle) + fieldRecvTurbZ[m][t1] * cos(rotationAngle);
                        }

                        for (int m = 0; m < nEquation; ++m)
                        {
                            fieldRecvTurbY[m][t1] = rotTurbgradValueY[m][t1];
                            fieldRecvTurbZ[m][t1] = rotTurbgradValueZ[m][t1];
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
                            rotTurbgradValueY[m][t1] = fieldRecvTurbY[m][t1] * cos(2 * PI - rotationAngle) - fieldRecvTurbZ[m][t1] * sin(2 * PI - rotationAngle);
                            rotTurbgradValueZ[m][t1] = fieldRecvTurbY[m][t1] * sin(2 * PI - rotationAngle) + fieldRecvTurbZ[m][t1] * cos(2 * PI - rotationAngle);
                        }

                        for (int m = 0; m < nEquation; ++m)
                        {
                            fieldRecvTurbY[m][t1] = rotTurbgradValueY[m][t1];
                            fieldRecvTurbZ[m][t1] = rotTurbgradValueZ[m][t1];
                        }
                    }
                    else if (bcName == "Periodic_down")
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {
                            rotTurbgradValueY[m][t1] = fieldRecvTurbY[m][t1] * cos(rotationAngle) - fieldRecvTurbZ[m][t1] * sin(rotationAngle);
                            rotTurbgradValueZ[m][t1] = fieldRecvTurbY[m][t1] * sin(rotationAngle) + fieldRecvTurbZ[m][t1] * cos(rotationAngle);
                        }

                        for (int m = 0; m < nEquation; ++m)
                        {
                            fieldRecvTurbY[m][t1] = rotTurbgradValueY[m][t1];
                            fieldRecvTurbZ[m][t1] = rotTurbgradValueZ[m][t1];
                        }
                    }
                }
            }
        }
    }
}

void TurbSolverUnstr::GetGradientField(Grid *gridIn)
{
    UnstructGrid *gridUnstruct = UnstructGridCast(gridIn);
    int numberOfEquations = this->GetNumberOfEquations();

    RDouble **gradientTurbulenceX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTurbulenceX"));
    RDouble **gradientTurbulenceY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTurbulenceY"));
    RDouble **gradientTurbulenceZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTurbulenceZ"));
    RDouble **qTurbulenceNode = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("qTurbNode"));
    RDouble **qTurbulence = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("q_turb"));

    for (int m = 0; m < numberOfEquations; ++ m)
    {
        string gradientName = GlobalDataBase::GetStrParaFromDB("gradientName");
        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGNode_NEW(qTurbulence[m], qTurbulenceNode[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }    
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }
}
#ifdef USE_GMRESSOLVER
//! GMRES solver -- a linear system solver
void TurbSolverUnstr::GMRESSolver(Grid *gridIn, FieldProxy *dqProxy)
{
    //! GMRESCoupled
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if ( viscousType == ONE_EQU )
    {
        return;
    }


    UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    //! Memory allocating
    RDouble **dRdq = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdq_turb"));

    int nTotalCells = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();

    RDouble **res = reinterpret_cast<RDouble **>(grid->GetDataPtr("res_turb"));
    RDouble* vol = grid->GetCellVolume();
    RDouble* dt = reinterpret_cast<RDouble*> (grid->GetDataPtr("dt"));
    vector<int> AI = grid->GetJacobianAI4GMRES();  // GMRESSparse GMRESCSR
    vector<int> AJ = grid->GetJacobianAJ4GMRES();  // GMRESSparse GMRESCSR

    //! adding diagoal for turbulence equations
    //! GMRESCSR
    for (int icell = 0; icell < nTotalCells; icell++)
    {
        vector<int>::iterator result = find(AJ.begin() + AI[icell], AJ.begin() + AI[icell + 1], icell); // find qcell
        int index = distance(AJ.begin(), result);
        index *= nTurbulenceEquation;
        for (int iequ = 0; iequ < nTurbulenceEquation; iequ++)
        {
             dRdq[iequ][index + iequ] += 1.0 / dt[icell];
        }
    }

    //! entry to PETSc for turbulence equations

    Vec x, b;
    Mat A;
    KSP ksp;                        // define solver
    PC pc;                          // define matrix
    PetscInt nTotalSize;            // dimension size of matrix and vector
    PetscInt nLevels = 3;           // define number of levels used by preconditioner ILU
    PetscInt MaxStepsRestart = 500;  // define maximum steps required by [restart]
    PetscInt MaxSteps = 500;        // define maximum steps
    PetscReal rTolerance = 1.e-1;   // define tolerance of relative error // -1
    PetscReal dropTolerance = 1.e5; // define tolerance of drop error

    nTotalSize = nTurbulenceEquation * nTotalCells;

    // initialize PETSc
    PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);

    // create vectors;
    VecCreate(PETSC_COMM_SELF, &x);
    PetscObjectSetName((PetscObject)x, "solution");
    VecSetSizes(x, PETSC_DECIDE, nTotalSize);
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
    MatCreate(PETSC_COMM_SELF, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nTotalSize, nTotalSize);
    MatSetType(A, MATBAIJ);
    MatSetBlockSize(A, nTurbulenceEquation);
    MatSeqBAIJSetPreallocation(A, nTurbulenceEquation, PETSC_DEFAULT, nnz);
    MatSetUp(A);

    // assemble matrix
    // improved approach to assemble matrix CSR
    double beforematrix = clock();
    PetscInt petsccol[nTurbulenceEquation];
    PetscInt petscrow;
    int jmax = AI[nTotalCells];
    int rcell = 0;
    int j = AI[rcell];
    do
    {
        int rowidx = rcell * nTurbulenceEquation;
        for (j = AI[rcell]; j < AI[rcell + 1]; j++)
        {
            int qcell = AJ[j]; // get qcell (colume index of matrix)
            // start to insert values: limited to dimension of nTotalSize
            if(qcell < nTotalCells)
            {
                int colidx = qcell * nTurbulenceEquation;

                for (int irow = 0; irow < nTurbulenceEquation; irow++)
                {
                    petscrow = rowidx + irow;
                    for (int n = 0; n < nTurbulenceEquation; n++)
                    {
                        petsccol[n] = colidx + n;
                    }
                    MatSetValues(A, 1, &petscrow, nTurbulenceEquation, petsccol, &dRdq[irow][j * nTurbulenceEquation], INSERT_VALUES);
                }
            }
        }
        rcell++; // update rcell, go to the next row index of matrix
    } while (j < jmax);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    // assemble vector
    VecSet(b, zero);
    for (PetscInt icell = 0; icell < nTotalCells; icell++)
    {
        for (PetscInt iequ = 0; iequ < nTurbulenceEquation; iequ++)
        {
            PetscInt isize = icell * nTurbulenceEquation + iequ;
            VecSetValues(b, 1, &isize, &res[iequ][icell], ADD_VALUES);
        }
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    //! define solver
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetPCSide(ksp, PC_RIGHT);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCILU);
    PCFactorSetLevels(pc, nLevels);
    KSPGMRESSetRestart(ksp, MaxStepsRestart);
    KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPSetTolerances(ksp, rTolerance, 1.e-20, dropTolerance, MaxSteps);

    //! solution
    KSPSolve(ksp, b, x);

    //! GMRESResidual
    RDouble **dq = dqProxy->GetField_UNS();

    //! convert x to dqProxy
    for (PetscInt icell = 0; icell < nTotalCells; icell++)
    {
        for (PetscInt iequ = 0; iequ < nTurbulenceEquation; iequ++)
        {
            PetscInt isize = icell * nTurbulenceEquation + iequ;
            VecGetValues(x, 1, &isize, &dq[iequ][icell]);
        }
    }

    //! free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
    PetscFinalize();
    return;
}
#endif
LIB_EXPORT void TurbSolverUnstr::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_TurbSolverUnstruct();
    controlParameters->Init();

    //! Define turbulent equation formation: Original SA / CFL3D SA / NSMB SA, or SST.
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    if (viscousType == ONE_EQU)
    {
        Turb_Equation_Formation = FORMATION_SA_NSMB;
    }
    else
    {
        Turb_Equation_Formation = FORMATION_ORGSA_SST;
    }

    //! Compute machin zero. reference from cfl3d.
    machZeroExp = 1;
    RDouble add, x11;
    RDouble compare = 1.0;
    for(int iCount = 0; iCount < 20;++ iCount)
    {
        add = 1.0;
        for(int n = 0; n < iCount + 1; ++ n)
        {
            add = add * 0.1;
        }

        x11 = compare +add;
        if(x11 == compare)
        {
            machZeroExp = iCount + 1;
            break;
        }
    }
}

LIB_EXPORT Param_TurbSolverUnstruct * TurbSolverUnstr::GetControlParameters()
{
    return static_cast<Param_TurbSolverUnstruct *>(controlParameters);
}

}
