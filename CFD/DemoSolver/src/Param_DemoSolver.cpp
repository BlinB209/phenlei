#include "Param_DemoSolver.h"
#include "GlobalDataBase.h"

namespace PHSPACE
{

LIB_EXPORT Param_DemoSolver::Param_DemoSolver()
{

}

LIB_EXPORT Param_DemoSolver::~Param_DemoSolver()
{

}

LIB_EXPORT void Param_DemoSolver::Init()
{
    Param_CFDSolver::Init();
    wallTemperature     = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    AoA                 = GlobalDataBase::GetDoubleParaFromDB("attackd");
    angleSlide          = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
}

}
