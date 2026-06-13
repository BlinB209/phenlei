#include"Param_INCompSolverUnstruct.h"
#include "GlobalDataBase.h"
#include "Glb_Dimension.h"

namespace PHSPACE
{

LIB_EXPORT Param_INCompSolverUnstruct::Param_INCompSolverUnstruct()
{
    setInComEquaMap();
}

LIB_EXPORT Param_INCompSolverUnstruct::~Param_INCompSolverUnstruct()
{

}

LIB_EXPORT void Param_INCompSolverUnstruct::Init()
{

}


void Param_INCompSolverUnstruct::setInComEquaMap()
{
    using namespace IDX;

    int *varIndex = new int[50];
    for (int iIndex = 0; iIndex < 50; iIndex++)
    {
        varIndex[iIndex] = 0;
    }
    GlobalDataBase::UpdateDataPtr("varIndex", varIndex);


    string* varNameIncom = new string[50];
    for (int iIndex = 0; iIndex < 50; iIndex++)
    {
        varNameIncom[iIndex] = "";
    }
    varNameIncom[S_IU] = "U";
    varNameIncom[S_IV] = "V";
    varNameIncom[S_IW] = "W";
    varNameIncom[S_IP] = "P";
    varNameIncom[S_ITEMP] = "Enthalpy";
    varNameIncom[S_ITURBK] = "Kinetic";
    varNameIncom[S_IEPSILON] = "Epsilon";
    varNameIncom[S_ISA] = "Kinetic";
    varNameIncom[S_IAIR  ] = "AIR";
    varNameIncom[S_IH2   ] = "H2";
    varNameIncom[S_IH2_l ] = "H2_l";
    varNameIncom[S_INH3_l] = "NH3_l";
    varNameIncom[S_ICH4  ] = "CH4";
    varNameIncom[S_IC2H4 ] = "C2H4";
    varNameIncom[S_IC3H8 ] = "C3H8";
    varNameIncom[S_IC2H6 ] = "C2H6";
    varNameIncom[S_ING   ] = "NG";
    varNameIncom[S_ILNG  ] = "LNG";
    varNameIncom[S_IH2S  ] = "H2S";
    varNameIncom[S_ICl2  ] = "Cl2";
    varNameIncom[S_INH3  ] = "NH3";

    GlobalDataBase::UpdateDataPtr("varNameIncom", varNameIncom);

    RDouble *soluRelaxCoeff = new RDouble[50];
    for (int iIndex = 0; iIndex < 50; iIndex++)
    {
        soluRelaxCoeff[iIndex] = 0.0;
    }
    GlobalDataBase::UpdateDataPtr("soluRelaxCoeff", soluRelaxCoeff);

    RDouble *solutionRes = new RDouble[50];
    for (int iIndex = 0; iIndex < 50; iIndex++)
    {
        solutionRes[iIndex] = 0.0;
    }
    GlobalDataBase::UpdateDataPtr("solutionRes", solutionRes);

    RDouble **DiffusionCoeff = new RDouble* [50];
    for (int iIndex = 0; iIndex < 50; iIndex++)
    {
        DiffusionCoeff[iIndex] = nullptr;
    }
    GlobalDataBase::UpdateDataPtr("DiffusionCoeff", DiffusionCoeff);

    int *localSolverIndex = new int[50];
    for (int iSolverIndex = 0; iSolverIndex < 50; iSolverIndex++)
    {
        localSolverIndex[iSolverIndex] = -1;
    }
    GlobalDataBase::UpdateDataPtr("localSolverIndex", localSolverIndex);

    int *maxIterOfEachEq = new int[50];
    for (int iEquation = 0; iEquation < 50; iEquation++)
    {
        maxIterOfEachEq[iEquation] = -1;
    }
    GlobalDataBase::UpdateDataPtr("maxIterOfEachEq", maxIterOfEachEq);

    RDouble *iterTolOfEachEq = new RDouble [50];
    for (int iEquation = 0; iEquation < 50; iEquation++)
    {
        iterTolOfEachEq[iEquation] = 1.0e-2;
    }
    GlobalDataBase::UpdateDataPtr("iterTolOfEachEq", iterTolOfEachEq);

    int *solMethodOfLinEqSystem = new int[50];
    for (int iEquation = 0; iEquation < 50; iEquation++)
    {
        solMethodOfLinEqSystem[iEquation] = -1;
    }
    GlobalDataBase::UpdateDataPtr("solMethodOfLinEqSystem", solMethodOfLinEqSystem);

    int *precondMethodOfLinEqSystem = new int[50];
    for (int iEquation = 0; iEquation < 50; iEquation++)
    {
        precondMethodOfLinEqSystem[iEquation] = -1;
    }
    GlobalDataBase::UpdateDataPtr("precondMethodOfLinEqSystem", precondMethodOfLinEqSystem);
}

}
