#include "Param_SpecSolver.h"
#include "GlobalDataBase.h"
#include "Constants.h"

namespace PHSPACE
{

LIB_EXPORT Param_SpecSolver::Param_SpecSolver()
{

}

LIB_EXPORT Param_SpecSolver::~Param_SpecSolver()
{

}

LIB_EXPORT void Param_SpecSolver::Init()
{
    maxSimuStep = GlobalDataBase::GetIntParaFromDB("maxSimuStep");

    statisticMethod = GlobalDataBase::GetIntParaFromDB("statisticMethod");

    startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");

    intervalStepRes = GlobalDataBase::GetIntParaFromDB("intervalStepRes");

    intervalStepStatistic = GlobalDataBase::GetIntParaFromDB("intervalStepStatistic");

    intervalStepFlow = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");

    startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");

    refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");

    viscousName = GlobalDataBase::GetStrParaFromDB("viscousName");

    viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    filterName = GlobalDataBase::GetStrParaFromDB("filterName");

    smagConstant = GlobalDataBase::GetDoubleParaFromDB("smagConstant");

    testFilterScale = GlobalDataBase::GetDoubleParaFromDB("testFilterScale");

    timeStep = GlobalDataBase::GetDoubleParaFromDB("timeStep");

    //pressureGradient = GlobalDataBase::GetDoubleParaFromDB("pressureGradient");

    spectrumType = GlobalDataBase::GetIntParaFromDB("spectrumType");
}

}