//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Param_NSSolver.h
//! @brief     Record paramters of NS Solver.
//! @author    Bell, Zhang Jian, Wan Yunbo, Meng Liyuan.

#pragma once
#include "Param_CFDSolver.h"

namespace PHSPACE
{

class Param_NSSolver : public Param_CFDSolver
{
public:
    LIB_EXPORT Param_NSSolver();

    LIB_EXPORT ~Param_NSSolver();

public:
    //! Init all parameters.
    LIB_EXPORT void Init();

    //! Get interval steps of dumping force.
    int GetIntervalStepForce() const;

    //! Get wall temperature.
    RDouble GetWallTemperature() const; 

    //! Get the flag of the gas model.
    int GetChemicalFlag() const;

    //! Get the set of WallMultiTemperature.
    int GetWallMultiTemperature() const;

    //! Get the number of N-S equations without the species transport equations and Tv-Te.
    int GetNSEquationNumber() const;

    //! Get the number of the equations including the species transport equations,but without Tv-Te.
    int GetLaminarNumber() const;

    //! Get the number of N-S equations for all.
    int GetnEquation() const;

    //! Get the flag of thermodynamic temperature model.
    int GetTemperatureModel() const;

    //! Get coefficient K for preconditioned reference Mach number.
    RDouble GetPreconCoefficient() const;

    //! Get precondition method for Farfield Boundary Condition.
    int GetPreconFarFieldBCMethod() const;

    //! Get if ideal gas state .
    int GetIfIdealGasState() const;

    //! Get Reference dimensional velocity .
    RDouble GetReferenceDimensionalVelocity() const;

    //! Get attacked of angle.
    RDouble GetAoA() const;

    //! Get the Entropy fix (correction) method.
    int GetRoeEntropyFixMethod() const;

    //! Get the Entropy coefficient scale, default is 1.0.
    RDouble GetRoeEntropyScale() const;

    //! Get the entropy fix coefficient 1 of Roe method.
    RDouble GetRoeEntropyFixCoef1() const; 
    
    //! Get the entropy fix coefficient 2 of Roe method.
    RDouble GetRoeEntropyFixCoef2() const; 

    //! Get angle of slide.
    RDouble GetAngleOfSlide () const;

    //! Set the entropy fix coefficients by the inflow mach number.
    void SetEntropyFixCoefficients();

    //! Get prandtlLaminar.
    RDouble GetPrandtlLaminar() const;

    //! Get prandtlTurbulence.
    RDouble GetPrandtlTurbulence() const;

    //! Get oPrandtlLaminar.
    RDouble GetoPrandtlLaminar() const;

    //! Get oPrandtlTurbulence.
    RDouble GetoPrandtlTurbulence() const;

    int GetNChemicalRadius() const;

    int GetNChemicalSource() const;

    //! Get timeIntegration.
    int GetTimeIntegration() const;

    //! Get nonDimensionalSutherlandTemperature.
    RDouble  GetNonDimensionalSutherlandTemperature() const;
    
    int GetMgProlongationType() const;

    RDouble * GetPrimitiveVarFarfield() const;

    RDouble * GetFullyCatalyticMassFraction() const;

    int GetNumberOfSpecies() const;

    int GetIndexOfNitrogen() const; //! The index of the species with the maximum mass fraction.
    //for Earth air, Nitrogen is the species with the maximum mass fraction, so it is used as the variable name.

    int GetIndexOfElectron() const; //! The index of Electron in the array of species.

    int GetFlagOfEquilibriumGas() const;

    int GetnDensityModify() const;

    int GetThermoEnergyModelFlag() const;

    int GetnDebug() const;

    int GetmTT() const;

    int GetmTV() const;

    int GetmTE() const;

    int GetnEnergyRecycle() const;

    int GetTransitionType() const;

    int GetnTurblenceForChemical() const;

    int GetnViscosityFluxSublevelModified() const;

    RDouble  GetDensityMin() const;

    RDouble  GetdensityMinFactor() const;

    int GetSelfAdaptionSolveFlag() const;

    int GetnGradPrimtiveMethod() const;

    int GetnInviscidFluxModify() const;    

    int GetnQlLimitMethod() const;

    int GetSurfaceGradientMethod() const;

    int GetRapidFlowfieldMethod() const;

    int GetSurfaceHeatingMonitor() const;

    int GetInitialPressureSteps() const;

    int GetLocalCFLFlag() const;

    RDouble GetCharacteristicLength() const;

    int GetNonequilibriumConditionFlag() const;

    //! Get if need plot volume field.
    int GetIsPorousZone() const;

    const RDouble *GetViscousResistanceCoeff() const;
    const RDouble *GetInertialResistanceCoeff() const;

    RDouble GetPorosity() const;
    RDouble GetDensitySolid() const;
    RDouble GetSpecificHeatSolid() const;
    RDouble GetHeatConductivitySolid() const;

    void SetNonequilibriumConditionFlag(int nUseNoneqCond) {isUseNoneqCond = nUseNoneqCond;}

private:
    int intervalStepForce;

    int transitionType;

    RDouble wallTemperature;

    //! The number of N-S equations without the species transport equations.
    int nNSEquation;
    //! The number of the equations including the species transport equations.
    int nLaminar;
    //! The flag of the gas model, 1 indicates the chemical non-equilibrium gas, 0 indicates the perfect gas.
    int nChemical;

    int nEquation;

    int wallMultiTemperature;

    //! The flag of thermodynamic temperature model, 1, 2 and 3 denote one/two/three temperature model for chemical non-equilibrium gas.
    int nTemperatureModel;

    //! Coefficient K for preconditioned reference Mach number.
    RDouble kPreconCoeff;

    //! For incompressible flow (mach < 0.3), choose precondition method for Farfield Boundary Condition.
    int preconFarfieldBCMethod;

    //! judge if ideal gas state or not.
    int nIdealState;

    //! Reference dimensional velocity.
    RDouble refDimensionalVelocity;

    //! Angle of attacked.
    double AoA;

    //! Angle of slide.
    RDouble angleSlide;

    //! Entropy fix (correction) method.\n
    //!   -# 1: direct fix, which limits the minimum eigenvalue directly.\n
    //!   -# 2: multi-dimensional fix, which is derived from structured solver and now is only valid for struct solver.\n
    //!   -# 3: Harten type, which is default used.
    int roeEntropyFixMethod;

    //! Entropy coefficient scale, default is 1.0.
    //! It is used to scale the default RoeEntropyFixCoef1 and RoeEntropyFixCoef2.
    RDouble roeEntropyScale;

    //! Entropy fix for Acoustic eigenvalue.
    RDouble RoeEntropyFixCoef1;

    //! Entropy fix for convective eigenvalue.
    RDouble RoeEntropyFixCoef2;

    RDouble prandtlLaminar;

    RDouble prandtlTurbulence;

    RDouble oPrandtlLaminar;

    RDouble oPrandtlTurbulence;

    int nChemicalRadius;

    int nChemicalSource;

    int timeIntegration;

    RDouble nonDimensionalSutherlandTemperature;

    int mgProlongationType;

    RDouble *primitiveVarFarfield;
    RDouble *catalyticMassFraction;

    int numberOfSpecies; //! The number of species.

    int indexOfNitrogen; //! The index of Nitrogen in the array of species, it indicates the species with the maximum mass fraction.
    int indexOfElectron; //! The index of Electron in the array of species.

    int nEquilibriumGas; //! The flag of equilibrium gas, 0 is for perfect gas, while > 0 is for equilibrium gas.
    //! Attention: this variable is valid when the condition nChemical=0 is satisfied, meanwhile its value denotes the number of gas components.

    int nDensityModify;

    int nDebug;

    //! nTEnergyModel=0 indicates the energy terms are computed using the conventional method.
    //! nTEnergyModel=1 indicates the energy terms are computed using the curve-fitting method.
    int nTEnergyModel;

    int mTT, mTV, mTE;
    int nEnergyRecycle;

    int nTurblenceForChemical;

    int nViscosityFluxSublevelModified;

    RDouble densityMin;

    RDouble densityMinFactor;

    int nGradPrimtiveMethod;

    int isSelfAdaptionSolve;

    int nInviscidFluxModify;

    int nQlLimitMethod;

    //! The method to compute the gradient of variable on surface.
    int nSurfGradMethod;

    //! To utilize the rapid method that can directly obtain the initial flowfield in boundary layer.
    int nRapidFlowfield;

    //! To utilize the monitor which exam the surface heating change during the iteration.
    int nSurfHeatMonitor;

    //! The steps to initialize the boundary variables with the rapid flowfield values.
    int nInitPressureStep;

    //! Use local CFL number or not.
    //!   -# 0: Global unified CFL number.\n
    //!   -# 1: Local CFL number.
    int isUseLocalCFL;

    RDouble knLength; //the characteristic length for computation of Knudsen number.
    int isUseNoneqCond; //Use the non-equilibrium condition or not.

    int isPorousZone;
    RDouble viscousResistanceCoeff[3];
    RDouble inertialResistanceCoeff[3];
    RDouble porosity;
    RDouble densitySolid;
    RDouble cpSolid;
    RDouble kSolid;
};

#include "Param_NSSolver.hxx"
}