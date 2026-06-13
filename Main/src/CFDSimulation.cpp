#include "GlobalDataBase.h"
#include "CFDSimulation.h"
#include "Simulation.h"
#include "TK_Exit.h"
#include "Constants.h"
#include "PHIO.h"
#include "TK_Log.h"
#include "Residual.h"
#include "Post_Visual.h"

#ifdef USE_CUDA
#include "GPUDeviceControl.h"
#endif

using namespace std;

namespace PHSPACE
{
void RunPHengLEI()
{
    //! Initialize global variables and MPI.

    PHMPI::Initialization();

#ifdef USE_CUDA
    //! Query and set the device
    int numberOfProcessors = PHMPI::GetNumberOfProcessor();
    GPUDeviceMultiGPU(numberOfProcessors, GPUNum, GPUProp);
#endif

    //! Read control parameters and do self-checking.
    ReadControlParameters();

    ShowWelcomeTitle();

    CompleteParamInfo();

    ControlParametersPreliminarySelfCheck();

    //! Start CFD simulation.
    int dimension = GlobalDataBase::GetIntParaFromDB("ndim");
    Simulation *simulation = new Simulation(dimension);

    simulation->Start();

    delete simulation;
}

void ShowWelcomeTitle()
{
    using namespace PHMPI;
    PH_Barrier();

    int myid = GetCurrentProcessorID();
    if (myid == PHMPI::GetServerProcessorID())
    {
        ostringstream oss;
        oss << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
        oss << "$#                                                     #$" << endl;
        oss << "$#  PPPPP H   H  EEEE N     N  GGGGG  L     EEEEE III  #$" << endl;
        oss << "$#  P   P H   H  E    N N   N  G      L     E      I   #$" << endl;
        oss << "$#  PPPPP HHHHH  EEEE N  N  N  G  GG  L     EEEEE  I   #$" << endl;
        oss << "$#  P     H   H  E    N   N N  G   G  L     E      I   #$" << endl;
        oss << "$#  P     H   H  EEEE N     N  GGGGG  LLLLL EEEEE III  #$" << endl;
        oss << "$#                                                     #$" << endl;
        oss << "$#=====================================================#$" << endl;
        oss << "$#                                                     #$" << endl;
        oss << "$# Platform for Hybrid Engineering Simulation of Flows #$" << endl;
        oss << "$#                                                     #$" << endl;
        oss << "$#       Computational Aerodynamics Institute          #$" << endl;
        oss << "$#  China Aerodynamics Research and Development Center #$" << endl;
        oss << "$#            (C) Copyright, Since 2010                #$" << endl;
        oss << "$#                  PHengLEI 2506                      #$" << endl;
        oss << "$#                                                     #$" << endl;
        oss << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

        cout << oss.str();
    }
}

void ParametersSelfCheckForCFD()
{
    //! Is in Multi-Grid initialization state.
    int isInMGIniting = 0;
    GlobalDataBase::UpdateData("isInMGIniting", &isInMGIniting, PHINT, 1);

    int isSolve = 0;
    GlobalDataBase::UpdateData("isSolve", &isSolve, PHINT, 1);

    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");
    if (isOversetSlip)
    {
        int oversetInterpolationMethod = 0;
        GlobalDataBase::UpdateData("oversetInterpolationMethod", &oversetInterpolationMethod, PHINT, 1);
    }

    //! Is in fine-Grid initialization state.
    int isInIniting = 1;
    GlobalDataBase::UpdateData("isInIniting", &isInIniting, PHINT, 1);

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

    int flowInitStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("flowInitStep");
    int maxSimuStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxSimuStep");
    if (flowInitStep > maxSimuStep || isUnsteady)
    {
        //! If the maximum simulation step is small than the initialization step,
        //! do not init then.
        flowInitStep = 0;
        GlobalDataBase::UpdateData("flowInitStep", &flowInitStep, PHINT, 1);
    }

    int visualVariablesSize = GlobalDataBase::GetSizeOfData("visualVariables");
    int nVisualVariables = visualVariablesSize;
    GlobalDataBase::UpdateData("nVisualVariables", &nVisualVariables, PHINT, 1);
    int visualWallVariablesSize = GlobalDataBase::GetSizeOfData("visualWallVariables");
    int nVisualWallVariables = visualWallVariablesSize;
    GlobalDataBase::UpdateData("nVisualWallVariables", &nVisualWallVariables, PHINT, 1);

    int *visualVariables = new int[nVisualVariables];
    GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);
    int *visualWallVariables = new int[nVisualWallVariables];
    GlobalDataBase::GetData("visualWallVariables", visualWallVariables, PHINT, nVisualWallVariables);
    RDouble wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    int plotFieldType = GlobalDataBase::GetIntParaFromDB("plotFieldType");
    string viscousName = GlobalDataBase::GetStrParaFromDB("viscousName");
    int isOverset = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");

    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
    int isUseNoneqCond = GlobalDataBase::GetIntParaFromDB("isUseNoneqCond");

    for (int i = 0; i < nVisualVariables; ++ i)
    {
        if (visualVariables[i] >= VISUAL_VISCOSITY_LAMINAR && visualVariables[i] < VISUAL_CP)
        {
            if (viscousType == INVISCID)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The euler solver cannot export the data related to viscosity." << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }

        if (visualVariables[i] >= VISUAL_VORTICITY_X && visualVariables[i] <= VISUAL_Q_CRITERIA)
        {
            if (plotFieldType == TEC_SPACE::BoundaryVisual)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The program can not visualization vorticity when plotFieldType = 0!" << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }

        if (visualVariables[i] == VISUAL_VISCOSITY_TURBULENT)
        {
            if (viscousType == LAMINAR && iLES == NOLES_SOLVER)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The laminar solver cannot export the data \"viscosityTurbulent\" related to turbulence." << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }

        if (visualVariables[i] >= VISUAL_MODELED_TKE && visualVariables[i] <= VISUAL_MODELED_DISSIPATION)
        {
            if (viscousType != TWO_EQU)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The solver cannot export the data \"modeledTKE,modeleddissipationrate\" related to 2eq-kw model." << endl;
                TK_Exit::ExceptionExit(oss);
            }
            else if (viscousType == TWO_EQU && viscousName.substr(0, 6) != "2eq-kw")
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The solver cannot export the data \"modeledTKE,modeleddissipationrate\" related to 2eq-kw model." << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }

        if (visualVariables[i] >= VISUAL_SST_F1 && visualVariables[i] <= VISUAL_SST_F2)
        {
            if (viscousType != TWO_EQU)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The solver cannot export the data \"SSTF1,SSTF2\" related to SST model." << endl;
                TK_Exit::ExceptionExit(oss);
            }
            else if (viscousType == TWO_EQU && viscousName.substr(0, 17) != "2eq-kw-menter-sst")
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The solver cannot export the data \"SSTF1,SSTF2\" related to SST model." << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }

        if (visualVariables[i] == VISUAL_INTERMITTENCY || visualVariables[i] == VISUAL_MOMENTUMTHICKREYNOLDS)
        {
            if (viscousType == INVISCID || viscousType == LAMINAR)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The solver cannot export the transition data \"intermittency,MomentumThicknessReynolds\" with current physical model." << endl;
                TK_Exit::ExceptionExit(oss);
            }
            else if (viscousType == ONE_EQU)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The solver cannot export the transition data \"intermittency,MomentumThicknessReynolds\" with SA model." << endl;
                TK_Exit::ExceptionExit(oss);
            }
            else if (viscousType == TWO_EQU)
            {
                int transitionType = GlobalDataBase::GetIntParaFromDB("transitionType");
                if (transitionType != IREGAMA)
                {
                    ostringstream oss;
                    oss << "transitionType = " << transitionType << endl;
                    oss << "The array \"visualVariables\" is error. The solver cannot export the transition data \"intermittency,MomentumThicknessReynolds\" with transition type Not Equal 2." << endl;
                    TK_Exit::ExceptionExit(oss);
                }
            }
        }

        if (visualVariables[i] >= VISUAL_TEMPERATURE_VIBRATION && visualVariables[i] <= VISUAL_ELECTRON_NUMBER)
        {
            if (nChemical == 0 || isUseNoneqCond == 0)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The solver cannot export the data \"Tv,Te,Ev,Ee,Ne\" related to Chemical Non-equilibrium" << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }

        if (visualVariables[i] == VISUAL_GAMA || visualVariables[i] == VISUAL_DAMKOHLER_NUMBER || visualVariables[i] == VISUAL_VIBNONEQ_NUMBER)
        {
            if (nChemical == 0 || isUseNoneqCond == 0)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The solver cannot export the data \"gama,Da,Vi\" related to Chemical Non-equilibrium" << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }

        if (visualVariables[i] == VISUAL_KNUDSEN_NUMBER)
        {
            ostringstream oss;
            oss << "The array \"visualVariables\" is error. The solver cannot export the data \"Kn\". Please dump out the data \"Kn_wall\" in wall_aircoef.dat" << endl;
            TK_Exit::ExceptionExit(oss);
        }

        if (visualVariables[i] == VISUAL_IBLANK)
        {
            if (!isOverset)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The solver cannot export the data related to OversetGrid." << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }
    }

    DelPointer(visualVariables);

    //! The parameters selfcheck of VisualWallVariables in wall_aircoef.dat.
    for (int i = 0; i < nVisualWallVariables; ++i)
    {
        if (visualWallVariables[i] >= VISUAL_WALL_QT && visualWallVariables[i] <= VISUAL_WALL_QE)
        {
            if (nChemical == 0)
            {
                ostringstream oss;
                oss << "The array \"visualWallVariables\" is error. The solver cannot export the data \"Qtr,Qs,Qv,Qe\" related to Chemical Non-equilibrium" << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }

        if (visualWallVariables[i] == VISUAL_SLIP_TS || visualWallVariables[i] == VISUAL_SLIP_TV 
            || visualWallVariables[i] == VISUAL_SLIP_TE || visualWallVariables[i] == VISUAL_SLIP_DTS)
        {
            if (nSlipBCModel == 0)
            {
                ostringstream oss;
                oss << "The array \"visualWallVariables\" is error. The solver cannot export the data \"Tts,Tvs,Tes,Vs\" related to slip temperature" << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }

        if (visualWallVariables[i] == VISUAL_WALL_ST)
        {
            if (wallTemperature < 0)
            {
                ostringstream oss;
                oss << "The array \"visualWallVariables\" is error. The solver cannot export the data \"St\" when the wall is adiabatic" << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }
    }

    DelPointer(visualWallVariables);

    if (isUnsteady)
    {
        int subTotalIterStep = 0;
        GlobalDataBase::UpdateData("subTotalIterStep", &subTotalIterStep, PHINT, 1);

        int minSubIter = GlobalDataBase::GetIntParaFromDB("min_sub_iter");
        if (minSubIter <= 0)
        {
            ostringstream oss;
            oss << "The variable \"min_sub_iter\" must be greater than 1." << endl;
            TK_Exit::ExceptionExit(oss);
        }

        int maxSubIter = GlobalDataBase::GetIntParaFromDB("max_sub_iter");
        if (minSubIter > maxSubIter)
        {
            ostringstream oss;
            oss << "Warning: The variable \"min_sub_iter\" cannot be greater than \"max_sub_iter\"! Now let them be equal: max_sub_iter = min_sub_iter!" << endl;
            PrintToWindow(oss.str());
            WriteLogFile(oss.str());

            maxSubIter = minSubIter;
            GlobalDataBase::UpdateData("max_sub_iter", &maxSubIter, PHINT, 1);
        }

        int ifLowSpeedPrecon =  GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");

        if (ifLowSpeedPrecon)
        {
            RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");
            RDouble referenceLength  = GlobalDataBase::GetDoubleParaFromDB("refLengthOfLowSpeed");
            RDouble gridScaleFactor  = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
            RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

            RDouble machUsteady = referenceLength / PI / physicalTimeStep * refMachNumber / gridScaleFactor;
            RDouble machUsteady2 = machUsteady * machUsteady;

            if (1.0 <= machUsteady2)
            {
                ifLowSpeedPrecon = 0;
                GlobalDataBase::UpdateData("ifLowSpeedPrecon", &ifLowSpeedPrecon, PHINT, 1);
                ostringstream oss;
                oss << "  Note: ifLowSpeedPrecon = 0, due to the physicalTimeStep being too small!" << endl;
                PrintToWindow(oss.str());
                WriteLogFile(oss.str());
    }
            GlobalDataBase::UpdateData("preconCoeffConst", &machUsteady2, PHDOUBLE, 1);
        }
    }

    //! If use Block-LUSGS, ROE scheme must be selected.
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    if (tscheme == Block_LU_SGS)
    {
        string uns_scheme_name = GlobalDataBase::GetStrParaFromDB("uns_scheme_name");
        if (uns_scheme_name.substr(0, 3) != "roe")
        {
            ostringstream oss;
            oss << "Error: if BLU is used, ROE scheme must be selected!" << endl;
            PrintToWindow(oss.str());
            WriteLogFile(oss.str());
            TK_Exit::ExceptionExit(oss);
        }
    }
    //! Judge if use wenn scheme, provided by Li Qin.
    int isWennScheme = 0;
    string str_limiter_name = GlobalDataBase::GetStrParaFromDB("str_limiter_name");
    if (str_limiter_name == "weno3_js"
        || str_limiter_name == "wenn3_zm"
        || str_limiter_name == "wenn3_zes2"
        || str_limiter_name == "wenn3_zes3"
        || str_limiter_name == "wenn3_prm211")
    {
        isWennScheme = 1;
    }
    GlobalDataBase::UpdateData("isWennScheme", &isWennScheme, PHINT, 1);
}

void ControlParametersPreliminarySelfCheck()
{
    int taskType = GetTaskCode();

    switch (taskType)
    {
        case SOLVE_FIELD:
            ParametersSelfCheckForCFD();
            break;
        case CREATE_GRID:
            break;
        case CAL_WALL_DIST:
            break;
        case PARTITION_GRID:
            break;
        case OVERSETGRID_VIEW:
            break;
        case OVERSET_CONFIG:
            break;
        case HO_SOLVER:
            break;
        case INTEGRATIVE_SOLVER:
            ParametersSelfCheckForCFD();
            break;
        case POST_PROCESSING:
            break;
        case SPECDIFFHYB_SOLVER:
            break;
        case SPEC_SOLVER:
            break;
 #ifdef USE_LBMSolverMPI
        case LBM_SOLVER_MPI:
            break;
        #endif
        #ifdef USE_LBMSolverOMP
        case LBM_SOLVER_OMP:
            break;
        #endif
        default:
            TK_Exit::UnexpectedVarValue("taskType = ", taskType);
            break;
    }
}

}



