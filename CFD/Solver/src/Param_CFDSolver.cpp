#include "Param_CFDSolver.h"
#include "GlobalDataBase.h"
#include "Constants.h"

namespace PHSPACE
{

LIB_EXPORT Param_CFDSolver::Param_CFDSolver()
{
    speedVaryStep = NULL;
    speedVaryCoef = NULL;
    CFLVaryMultiStep = NULL;
    CFLVaryCoef = NULL;
}

LIB_EXPORT Param_CFDSolver::~Param_CFDSolver()
{
    if (speedVaryStep != NULL) { delete speedVaryStep; speedVaryStep = NULL; }
    if (speedVaryCoef != NULL) { delete speedVaryCoef; speedVaryCoef = NULL; }
    if (CFLVaryMultiStep != NULL) { delete CFLVaryMultiStep; CFLVaryMultiStep = NULL; }
    if (CFLVaryCoef != NULL) { delete CFLVaryCoef; CFLVaryCoef = NULL; }
}

LIB_EXPORT void Param_CFDSolver::Init()
{
    restartNSFile            = GlobalDataBase::GetStrParaFromDB("restartNSFile");
    resSaveFile              = GlobalDataBase::GetStrParaFromDB("resSaveFile");
    viscousName              = GlobalDataBase::GetStrParaFromDB("viscousName");

    plotFieldType            = GlobalDataBase::GetIntParaFromDB("plotFieldType");

    viscousType              = GlobalDataBase::GetIntParaFromDB("viscousType");
    isUnsteady               = GlobalDataBase::GetIntParaFromDB("iunsteady");

    intervalStepRes          = GlobalDataBase::GetIntParaFromDB("intervalStepRes");
    intervalStepPlot         = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
    intervalStepFlow         = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");

    isAle                    = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    ifLocalTimeStep          = GlobalDataBase::GetIntParaFromDB("ifLocalTimeStep");
    ifLowSpeedPrecon         = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");

    CFLVaryStep              = GlobalDataBase::GetIntParaFromDB("CFLVaryStep");
    nLUSGSSweeps             = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");
    nMGLevel                 = GlobalDataBase::GetIntParaFromDB("nMGLevel");

    ifStartFromSteadyResults = GlobalDataBase::GetIntParaFromDB("ifStartFromSteadyResults");
    ifStaticsFlowField       = GlobalDataBase::GetIntParaFromDB("ifStaticsFlowField");
    ifStaticsReynoldsStress  = 0;
    if (GlobalDataBase::IsExist("ifStaticsReynoldsStress", PHINT, 1))
    {
        ifStaticsReynoldsStress = GlobalDataBase::GetIntParaFromDB("ifStaticsReynoldsStress");
    }
    statisticMethod = 0;
    if (GlobalDataBase::IsExist("statisticMethod", PHINT, 1))
    {
        statisticMethod = GlobalDataBase::GetIntParaFromDB("statisticMethod");
    }
    isRestartChangeInflow     = GlobalDataBase::GetIntParaFromDB("isRestartChangeInflow");
    CFLStart                  = GlobalDataBase::GetDoubleParaFromDB("CFLStart");
    CFLEnd                    = GlobalDataBase::GetDoubleParaFromDB("CFLEnd");

    LUSGSTolerance            = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");

    refMachNumber             = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    refReNumber               = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    oRefReNumber              = GlobalDataBase::GetDoubleParaFromDB("oreynolds");

    refGama                   = GlobalDataBase::GetDoubleParaFromDB("refGama");

    TorqueRefX                = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    TorqueRefY                = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    TorqueRefZ                = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    refDimensionalPressure    = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
    refDimensionalDensity     = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    refDimensionalVelocity    = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    refDimensionalSonicSpeed  = GlobalDataBase::GetDoubleParaFromDB("refDimensionalSonicSpeed");
    refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    GlobalDataBase::GetData("lowerPlotFieldBox", &lowerPlotFieldBox, PHDOUBLE, 3);
    GlobalDataBase::GetData("upperPlotFieldBox", &upperPlotFieldBox, PHDOUBLE, 3);

    isOverset = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    isWennScheme = GlobalDataBase::GetIntParaFromDB("isWennScheme");

    pMax = 0.5, pMin = 0.2, deltaMax = 0.5, magnifyFactor = 2.0;
    if (GlobalDataBase::IsExist("pMaxForCFL", PHDOUBLE, 1))
    {
        pMax = GlobalDataBase::GetDoubleParaFromDB("pMaxForCFL");
    }
    if (GlobalDataBase::IsExist("pMinForCFL", PHDOUBLE, 1))
    {
        pMin = GlobalDataBase::GetDoubleParaFromDB("pMinForCFL");
    }
    if (GlobalDataBase::IsExist("deltaMaxForCFL", PHDOUBLE, 1))
    {
        deltaMax = GlobalDataBase::GetDoubleParaFromDB("deltaMaxForCFL");
    }
    if (GlobalDataBase::IsExist("magnifyFactorForCFL", PHDOUBLE, 1))
    {
        magnifyFactor = GlobalDataBase::GetDoubleParaFromDB("magnifyFactorForCFL");
    }
    reduceFactor = 1.0 / magnifyFactor;
    if (GlobalDataBase::IsExist("reduceFactorForCFL", PHDOUBLE, 1))
    {
        reduceFactor = GlobalDataBase::GetDoubleParaFromDB("reduceFactorForCFL"); //reduceFactorForCFL
    }
    fineCFL = 0.5;
    if (GlobalDataBase::IsExist("fineCFL", PHDOUBLE, 1))
    {
        fineCFL = GlobalDataBase::GetDoubleParaFromDB("fineCFL"); //fineCFL
    }

    AusmpwPlusLimiter = 1.0;
    if (GlobalDataBase::IsExist("AusmpwPlusLimiter", PHDOUBLE, 1))
    {
        AusmpwPlusLimiter = GlobalDataBase::GetDoubleParaFromDB("AusmpwPlusLimiter"); //reduceFactorForCFL
    }

    skewnessAngle = GlobalDataBase::GetDoubleParaFromDB("skewnessAngle");
    turbSkewnessAngle = GlobalDataBase::GetDoubleParaFromDB("turbSkewnessAngle");

    nDiagonalModifiedTurb = 0;

    trTemperatureMinNonDim = 10.0 ;
    if (GlobalDataBase::IsExist("trTemperatureMin", PHDOUBLE, 1))
    {
        trTemperatureMinNonDim = GlobalDataBase::GetDoubleParaFromDB("trTemperatureMin") ;
    }
    trTemperatureMinNonDim /= refDimensionalTemperature;

//////////////////////////////////////////////////////////////////////
    nNumberOfSpeedStep = 0;
    if (GlobalDataBase::IsExist("nNumberOfSpeedStep", PHINT, 1))
    {
        nNumberOfSpeedStep = GlobalDataBase::GetIntParaFromDB("nNumberOfSpeedStep");
    }
    if (nNumberOfSpeedStep > 0)
    {
        speedVaryStep = new int[nNumberOfSpeedStep];
        speedVaryCoef = new RDouble[nNumberOfSpeedStep];

        if (GlobalDataBase::IsExist("speedVaryStep", PHINT, nNumberOfSpeedStep))
        {
            GlobalDataBase::GetData("speedVaryStep", &speedVaryStep[0], PHINT, nNumberOfSpeedStep);
        }
        else
        {
            nNumberOfSpeedStep = 0;
        }

        if (GlobalDataBase::IsExist("speedVaryCoef", PHDOUBLE, nNumberOfSpeedStep))
        {
            GlobalDataBase::GetData("speedVaryCoef", &speedVaryCoef[0], PHDOUBLE, nNumberOfSpeedStep);
        }
        else
        {
            nNumberOfSpeedStep = 0;
        }
    }

//////////////////////////////////////////////////////////////////////
    nNumberOfCFLStep = 0;
    if (GlobalDataBase::IsExist("nNumberOfCFLStep", PHINT, 1))
    {
        nNumberOfCFLStep = GlobalDataBase::GetIntParaFromDB("nNumberOfCFLStep");
    }

    if (nNumberOfCFLStep > 0)
    {
        CFLVaryMultiStep = new int[nNumberOfCFLStep];
        CFLVaryCoef = new RDouble[nNumberOfCFLStep];

        if (GlobalDataBase::IsExist("CFLVaryMultiStep", PHINT, nNumberOfCFLStep))
        {
            GlobalDataBase::GetData("CFLVaryMultiStep", &CFLVaryMultiStep[0], PHINT, nNumberOfCFLStep);
        }
        else
        {
            nNumberOfCFLStep = 0;
        }

        if (GlobalDataBase::IsExist("CFLVaryCoef", PHDOUBLE, nNumberOfCFLStep))
        {
            GlobalDataBase::GetData("CFLVaryCoef", &CFLVaryCoef[0], PHDOUBLE, nNumberOfCFLStep);
        }
        else
        {
            nNumberOfCFLStep = 0;
        }
    }
}

}