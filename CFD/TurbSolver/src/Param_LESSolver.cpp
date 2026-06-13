#include "Param_LESSolver.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "TK_Exit.h"

namespace PHSPACE
{

LIB_EXPORT Param_LESSolver::Param_LESSolver()
{

}

LIB_EXPORT Param_LESSolver::~Param_LESSolver()
{
    if (filterDirection != NULL)
    {
        delete [] filterDirection;
        filterDirection = NULL;
    }

    if (averageDirection != NULL)
    {
        delete [] averageDirection;
        averageDirection = NULL;
    }

    if (testFilterCoef != NULL)
    {
        delete [] testFilterCoef;
        testFilterCoef = NULL;
    }
}

LIB_EXPORT void Param_LESSolver::Init()
{
    Param_CFDSolver::Init();

    testFilterCoef = new RDouble [3]();

    filterDirection = new int [3]();

    averageDirection = new int [3]();

    eddyViscosityLimit = GlobalDataBase::GetDoubleParaFromDB("eddyViscosityLimit");

    wallDampingFunctionType = GlobalDataBase::GetIntParaFromDB("wallDampingFunctionType");

    deltaFunctionType = GlobalDataBase::GetIntParaFromDB("deltaFunctionType");

    turbViscousCutType = GlobalDataBase::GetIntParaFromDB("turbViscousCutType");

    smagConstant = GlobalDataBase::GetDoubleParaFromDB("smagConstant");

    isotropicConstant = GlobalDataBase::GetDoubleParaFromDB("isotropicConstant");

    string sgsmodel;
    GlobalDataBase::GetData("sgsmodel", &sgsmodel, PHSTRING, 1);

    if (sgsmodel.substr(0, 11) == "smagorinsky")
    {
        subgridScaleModel = 1;
    }
    else if (sgsmodel.substr(0, 6) == "dsmCom")
    {
        subgridScaleModel = 2;

        testFilterCoef[0] = 0.166666666666;
        testFilterCoef[1] = 0.666666666667;
        testFilterCoef[2] = 0.166666666667;

        testFilterWidth = 1;

        GlobalDataBase::GetData("filterDirection", filterDirection, PHINT, 3);

        GlobalDataBase::GetData("averageDirection", averageDirection, PHINT, 3);

        averageWidth = GlobalDataBase::GetIntParaFromDB("averageWidth");

        testFilterScale = GlobalDataBase::GetDoubleParaFromDB("testFilterScale");
    }
    else if (sgsmodel.substr(0, 4) == "wale")
    {
        subgridScaleModel = 3;
        waleConstant = GlobalDataBase::GetDoubleParaFromDB("waleConstant");
        isotropicConstant = GlobalDataBase::GetDoubleParaFromDB("isotropicConstant");
    }
    else if (sgsmodel.substr(0, 6) == "dsmInc")
    {
        subgridScaleModel = 4;

        testFilterCoef[0] = 0.166666666666;
        testFilterCoef[1] = 0.666666666667;
        testFilterCoef[2] = 0.166666666667;

        testFilterWidth = 1;

        GlobalDataBase::GetData("filterDirection", filterDirection, PHINT, 3);

        GlobalDataBase::GetData("averageDirection", averageDirection, PHINT, 3);

        averageWidth = GlobalDataBase::GetIntParaFromDB("averageWidth");

        testFilterScale = GlobalDataBase::GetDoubleParaFromDB("testFilterScale");
    }
    else if (sgsmodel.substr(0, 6) == "dsmOld")
    {
        subgridScaleModel = 5;

        testFilterCoef[0] = 0.166666666666;
        testFilterCoef[1] = 0.666666666667;
        testFilterCoef[2] = 0.166666666667;

        testFilterWidth = 1;

        GlobalDataBase::GetData("filterDirection", filterDirection, PHINT, 3);

        GlobalDataBase::GetData("averageDirection", averageDirection, PHINT, 3);

        averageWidth = GlobalDataBase::GetIntParaFromDB("averageWidth");

        testFilterScale = GlobalDataBase::GetDoubleParaFromDB("testFilterScale");
    }
    else if (sgsmodel.substr(0, 5) == "sigma")
    {
        subgridScaleModel = 6;
        sigmaConstant = GlobalDataBase::GetDoubleParaFromDB("sigmaConstant");
        isotropicConstant = GlobalDataBase::GetDoubleParaFromDB("isotropicConstant");
    }
    else
    {
        TK_Exit::UnexpectedVarValue("sgsmodel", sgsmodel);
    }

    if (wallDampingFunctionType == 2)
    {
        yploge = 1000.0;
        rgam = 0.41;
    }

    prandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("prt");

    oPrandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("oprt");

    freeStreamViscosity = GlobalDataBase::GetDoubleParaFromDB("freeStreamViscosity");
}

}