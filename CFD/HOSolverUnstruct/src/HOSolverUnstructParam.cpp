#include "HOSolverUnstructParam.h"
#include "GlobalDataBase.h"

#include "HODefine.h"

namespace PHSPACE
{

LIB_EXPORT HOSolverUnstructParam::HOSolverUnstructParam()
{
}

LIB_EXPORT HOSolverUnstructParam::~HOSolverUnstructParam()
{
}

LIB_EXPORT void HOSolverUnstructParam::Init()
{
    Param_CFDSolver::Init();

    restartNSFile = GlobalDataBase::GetStrParaFromDB("restartNSFile");
    resSaveFile = GlobalDataBase::GetStrParaFromDB("resSaveFile");
    intervalStepForce = GlobalDataBase::GetIntParaFromDB("intervalStepForce");
    wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    ifLowSpeedPrecon = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");
    AoA  = GlobalDataBase::GetDoubleParaFromDB("attackd");
    GlobalDataBase::GetData("alf_l", &RoeEntropyFixCoef1, PHDOUBLE, 1);
    GlobalDataBase::GetData("alf_n", &RoeEntropyFixCoef2, PHDOUBLE, 1);
    angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    dgSolOrder = GlobalDataBase::GetIntParaFromDB("dg_sol_order");
    fvSolOrder = GlobalDataBase::GetIntParaFromDB("fv_sol_order");

    dgSolDof = (dgSolOrder+1)*(dgSolOrder+2)*(dgSolOrder+3)/6;
    fvSolDof = (fvSolOrder+1)*(fvSolOrder+2)*(fvSolOrder+3)/6;

    pMultiGrid = GlobalDataBase::GetIntParaFromDB("p_multigrid");
    if (pMultiGrid < 1) pMultiGrid = 1;

    uns_scheme_name = GlobalDataBase::GetStrParaFromDB("uns_scheme_name");

    RDouble Ma = GetRefMachNumber();
    R = 1 / (HO_GAMA * Ma * Ma);
}

}
