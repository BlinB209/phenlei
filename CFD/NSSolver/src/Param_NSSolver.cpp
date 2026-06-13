#include "Param_NSSolver.h"
#include "GlobalDataBase.h"
#include "TK_Log.h"
#include "Constants.h"

namespace PHSPACE
{

LIB_EXPORT Param_NSSolver::Param_NSSolver()
{
    
}

LIB_EXPORT Param_NSSolver::~Param_NSSolver()
{

}

LIB_EXPORT void Param_NSSolver::Init()
{
    Param_CFDSolver::Init();
    intervalStepForce   = GlobalDataBase::GetIntParaFromDB("intervalStepForce");
    wallTemperature     = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    AoA                 = GlobalDataBase::GetDoubleParaFromDB("attackd");
    angleSlide          = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    transitionType      = GlobalDataBase::GetIntParaFromDB("transitionType");

    kPreconCoeff    = GlobalDataBase::GetDoubleParaFromDB("Kprec");
    nIdealState = GlobalDataBase::GetIntParaFromDB("nIdealState");
    preconFarfieldBCMethod    = GlobalDataBase::GetIntParaFromDB("preconFarfieldBCMethod");
    refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");

    roeEntropyFixMethod = GlobalDataBase::GetIntParaFromDB("roeEntropyFixMethod");
    roeEntropyScale     = GlobalDataBase::GetDoubleParaFromDB("roeEntropyScale");

    nNSEquation         = GlobalDataBase::GetIntParaFromDB("nm");
    nLaminar            = GlobalDataBase::GetIntParaFromDB("nl");
    nChemical           = GlobalDataBase::GetIntParaFromDB("nchem");
    nTemperatureModel   = GlobalDataBase::GetIntParaFromDB("ntmodel");
    numberOfSpecies     = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");

    nEquation = nLaminar + nChemical + nTemperatureModel - 1;

    nDensityModify      = GlobalDataBase::GetIntParaFromDB("nDensityModify");
    nDebug              = GlobalDataBase::GetIntParaFromDB("nDebug");
    nTEnergyModel       = GlobalDataBase::GetIntParaFromDB("nTEnergyModel");

    nEnergyRecycle      = GlobalDataBase::GetIntParaFromDB("nEnergyRecycle") * nChemical;

    nTurblenceForChemical = GlobalDataBase::GetIntParaFromDB("nTurblenceForChemical");

    nViscosityFluxSublevelModified = GlobalDataBase::GetIntParaFromDB("nViscosityFluxSublevelModified");
    isPorousZone = GlobalDataBase::GetIntParaFromDB("isPorousZone");

    if (nChemical == 0)
    {
        if (nTemperatureModel > 1)
        {
            nTemperatureModel = 1;
            GlobalDataBase::UpdateData("ntmodel", &nTemperatureModel, PHINT, 1);
        }
        if (nLaminar > nNSEquation)
        {
            nLaminar = nNSEquation;
            GlobalDataBase::UpdateData("nl", &nLaminar, PHINT, 1);
        }
        if (numberOfSpecies > 0)
        {
            numberOfSpecies = 0;
            GlobalDataBase::UpdateData("numberOfSpecies", &numberOfSpecies, PHINT, 1);
        }
    }

    prandtlLaminar    = GlobalDataBase::GetDoubleParaFromDB("prl");
    prandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("prt");

    oPrandtlLaminar    = GlobalDataBase::GetDoubleParaFromDB("oprl");
    oPrandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("oprt");

    nChemicalRadius = GlobalDataBase::GetIntParaFromDB("nchemrad");
    nChemicalSource = GlobalDataBase::GetIntParaFromDB("nchemsrc");

    timeIntegration = GlobalDataBase::GetIntParaFromDB("tscheme");

    nonDimensionalSutherlandTemperature = GlobalDataBase::GetDoubleParaFromDB("tsuth");

    mgProlongationType = GlobalDataBase::GetIntParaFromDB("mprol");

    if (nChemical > 0)
    {
        indexOfNitrogen = GlobalDataBase::GetIntParaFromDB("nNitrogenIndex");
        indexOfElectron = GlobalDataBase::GetIntParaFromDB("nElectronIndex");
        catalyticMassFraction = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("catalyticMassFraction"));
    }
    else
    {
        indexOfNitrogen = -1;
        indexOfElectron = -1;
        catalyticMassFraction = 0;
    }
    nEquilibriumGas = GlobalDataBase::GetIntParaFromDB("nEquilibriumGas");

    primitiveVarFarfield = reinterpret_cast <RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));

    mTT = 0;
    if(nTemperatureModel == 3)
    {
        mTV = 1;
        mTE = 2;
    }
    else if(nTemperatureModel ==2)
    {
        mTV = 1;
        mTE = 1;
    }
    else
    {
        mTV = 0;
        mTE = 0;
    }

    densityMin = 1.0e-8;
    if (GlobalDataBase::IsExist("densityMin", PHDOUBLE, 1))
    {
        densityMin = GlobalDataBase::GetDoubleParaFromDB("densityMin");
    }
    densityMin /= GetRefDimensionalDensity();

    densityMinFactor = SMALL;
    if (GlobalDataBase::IsExist("densityMinFactor", PHDOUBLE, 1) && nChemical == 1)
    {
        densityMinFactor = GlobalDataBase::GetDoubleParaFromDB("densityMinFactor");
    }

    isSelfAdaptionSolve = GlobalDataBase::GetIntParaFromDB("isAdaptiveSolver");

    nGradPrimtiveMethod = 0;
    if (GlobalDataBase::IsExist("nGradPrimtiveMethod", PHINT, 1) && nChemical == 1)
    {
        nGradPrimtiveMethod = GlobalDataBase::GetIntParaFromDB("nGradPrimtiveMethod");
    }

    nInviscidFluxModify = 0;
    if (GlobalDataBase::IsExist("nInviscidFluxModify", PHINT, 1))
    {
        nInviscidFluxModify = nChemical * GlobalDataBase::GetIntParaFromDB("nInviscidFluxModify");  //Chemical only
    }

    nQlLimitMethod = 0;
    if (GlobalDataBase::IsExist("nQlLimitMethod", PHINT, 1))
    {
        nQlLimitMethod = GlobalDataBase::GetIntParaFromDB("nQlLimitMethod");  //ql¡¢qr
    }

    //! The method to compute the gradient of variable on surface.
    nSurfGradMethod = 0;
    if (GlobalDataBase::IsExist("nSurfGradMethod", PHINT, 1))
    {
        nSurfGradMethod = GlobalDataBase::GetIntParaFromDB("nSurfGradMethod");
    }

    //! To utilize the rapid method that can directly obtain the initial flowfield in boundary layer.
    nRapidFlowfield = 0;
    if (GlobalDataBase::IsExist("nRapidFlowfield", PHINT, 1))
    {
        nRapidFlowfield = GlobalDataBase::GetIntParaFromDB("nRapidFlowfield");
    }

    //! To utilize the monitor which exam the surface heating change during the iteration.
    nSurfHeatMonitor = 0;
    if (GlobalDataBase::IsExist("nSurfHeatMonitor", PHINT, 1))
    {
        nSurfHeatMonitor = GlobalDataBase::GetIntParaFromDB("nSurfHeatMonitor");
    }

    wallMultiTemperature = 0;
    if (GlobalDataBase::IsExist("wallMultiTemperature", PHINT, 1))
    {
        wallMultiTemperature = GlobalDataBase::GetIntParaFromDB("wallMultiTemperature");
    }

    //! The steps to initialize the boundary variables with the rapid flowfield values.
    nInitPressureStep = 100;
    if (GlobalDataBase::IsExist("nInitPressureStep", PHINT, 1))
    {
        nInitPressureStep = GlobalDataBase::GetIntParaFromDB("nInitPressureStep");
    }

    //! Use local CFL number or not.
    isUseLocalCFL = 0;
    if (GlobalDataBase::IsExist("isUseLocalCFL", PHINT, 1))
    {
        isUseLocalCFL = GlobalDataBase::GetIntParaFromDB("isUseLocalCFL");
    }

    //! the reference length for computation of Knudsen number.
    knLength = 1.0;
    if (GlobalDataBase::IsExist("knudsenLength", PHDOUBLE, 1))
    {
        knLength = GlobalDataBase::GetDoubleParaFromDB("knudsenLength");
    }

    //! the flag of using the non-equilibrium condition.
    isUseNoneqCond = 0;
    if (GlobalDataBase::IsExist("isUseNoneqCond", PHINT, 1))
    {
        isUseNoneqCond = GlobalDataBase::GetIntParaFromDB("isUseNoneqCond");
    }

    if (isPorousZone != NON_POROUS)
    {
        GlobalDataBase::GetData("viscousResistanceCoeff", &viscousResistanceCoeff, PHDOUBLE, 3);
        GlobalDataBase::GetData("inertialResistanceCoeff", &inertialResistanceCoeff, PHDOUBLE, 3);
        GlobalDataBase::GetData("porosity", &porosity, PHDOUBLE, 1);
        GlobalDataBase::GetData("densitySolid", &densitySolid, PHDOUBLE, 1);
        GlobalDataBase::GetData("cpSolid", &cpSolid, PHDOUBLE, 1);
        GlobalDataBase::GetData("kSolid", &kSolid, PHDOUBLE, 1);

    }

}

void Param_NSSolver::SetEntropyFixCoefficients()
{
    RDouble entropyAcoustic = 0.01, entropyConvective = 0.01;
    RDouble referenceMachNumber = this->GetRefMachNumber();

    int ifLowSpeedPrecon = this->GetIfLowSpeedPrecon();
    if (ifLowSpeedPrecon)
    {
        //! Incompressible, precondition.
        entropyAcoustic   = 1.0e-8;
    }
    else
    {
        //! Compressible.
        RDouble x0, x1, y0, y1;
        if(referenceMachNumber <= 0.6)
        {
            entropyAcoustic = 0.0001;
        }
        else if(referenceMachNumber <= 2.0)
        {
            x0 = 0.6;    x1 = 2.0;
            y0 = 0.0001; y1 = 0.01;
            entropyAcoustic = y0 * (referenceMachNumber - x1) / (x0 - x1) + y1 * (referenceMachNumber - x0) / (x1 - x0);
        }
        else if(referenceMachNumber <= 5.0)
        {
            x0 = 2.0;  x1 = 5.0;
            y0 = 0.01; y1 = 0.5;
            entropyAcoustic = y0 * (referenceMachNumber - x1) / (x0 - x1) + y1 * (referenceMachNumber - x0) / (x1 - x0);
        }
        else
        {
            entropyAcoustic = 1.0;
        }  

        //! Do not scale the entropy for incompressible precondition.
        entropyAcoustic *= this->roeEntropyScale;
    }

    entropyConvective = entropyAcoustic;

    this->RoeEntropyFixCoef1 = entropyAcoustic;
    this->RoeEntropyFixCoef2 = entropyConvective;
}

}
