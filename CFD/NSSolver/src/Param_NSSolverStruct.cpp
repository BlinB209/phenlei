#include "Param_NSSolverStruct.h"
#include "GlobalDataBase.h"
#include "Constants.h"

namespace PHSPACE
{

LIB_EXPORT Param_NSSolverStruct::Param_NSSolverStruct()
{

}

LIB_EXPORT Param_NSSolverStruct::~Param_NSSolverStruct()
{

}

LIB_EXPORT void Param_NSSolverStruct::Init()
{
    Param_NSSolver::Init();
    MUSCLCoefXk                = GlobalDataBase::GetDoubleParaFromDB("MUSCLCoefXk");
    MUSCLCoefXb                = GlobalDataBase::GetDoubleParaFromDB("MUSCLCoefXb");
    catalyticCoef              = GlobalDataBase::GetDoubleParaFromDB("catalyticCoef");
    chemicalSpectrumRadiusCoef = GlobalDataBase::GetDoubleParaFromDB("chemicalSpectrumRadiusCoef");
    viscousSpectrumRadiusCoef  = GlobalDataBase::GetDoubleParaFromDB("viscousSpectrumRadiusCoef");
    inviscidSpectrumRadiusCoef = GlobalDataBase::GetDoubleParaFromDB("inviscidSpectrumRadiusCoef");
    chemicalRelaxCorf          = GlobalDataBase::GetDoubleParaFromDB("chemicalRelaxCorf");
    staticPressureRelaxCorf    = GlobalDataBase::GetDoubleParaFromDB("staticPressureRelaxCorf");
    inviscidSchemeName         = GlobalDataBase::GetStrParaFromDB("inviscidSchemeName");

    if(GetChemicalFlag()== 0)
    {
        chemicalSpectrumRadiusCoef = 1.0;
        viscousSpectrumRadiusCoef  = 1.0;
        inviscidSpectrumRadiusCoef = 1.0;
        staticPressureRelaxCorf    = 1.0;
    }

    iniSpeedCoef = 1.0;
    if (GlobalDataBase::IsExist("iniSpeedCoef", PHDOUBLE, 1))
    {
        iniSpeedCoef = GlobalDataBase::GetDoubleParaFromDB("iniSpeedCoef");
    }
    iniSpeedMode = 0;
    if(GlobalDataBase::GetIntParaFromDB("nIsComputeWallDist")== 0 && GlobalDataBase::IsExist("iniSpeedMode", PHINT, 1))
    {
        iniSpeedMode = GlobalDataBase::GetIntParaFromDB("iniSpeedMode");
    }
}

}