#include "Param_NSSolverUnstruct.h"
#include "GlobalDataBase.h"

namespace PHSPACE
{

LIB_EXPORT Param_NSSolverUnstruct::Param_NSSolverUnstruct()
{

}

LIB_EXPORT Param_NSSolverUnstruct::~Param_NSSolverUnstruct()
{

}

LIB_EXPORT void Param_NSSolverUnstruct::Init()
{
    Param_NSSolver::Init();
    uns_scheme_name = GlobalDataBase::GetStrParaFromDB("uns_scheme_name");

    gradientName = GlobalDataBase::GetStrParaFromDB("gradientName");

    limitVector = GlobalDataBase::GetIntParaFromDB("limitVector");

    limitVariables = GlobalDataBase::GetIntParaFromDB("limitVariables");

    sliceAxis = GlobalDataBase::GetIntParaFromDB("sliceAxis");

    slicePostion = GlobalDataBase::GetDoubleParaFromDB("slicePostion");

    mgCFLScale = GlobalDataBase::GetDoubleParaFromDB("mgCFLScale");

    mgCorrectionLimit = GlobalDataBase::GetDoubleParaFromDB("mgCorrectionLimit");

    catalyticCoef              = GlobalDataBase::GetDoubleParaFromDB("catalyticCoef");
    chemicalSpectrumRadiusCoef = GlobalDataBase::GetDoubleParaFromDB("chemicalSpectrumRadiusCoef");
    chemicalRelaxCorf          = GlobalDataBase::GetDoubleParaFromDB("chemicalRelaxCorf");
    staticPressureRelaxCorf    = GlobalDataBase::GetDoubleParaFromDB("staticPressureRelaxCorf");

}

}