inline int Param_NSSolver::GetIntervalStepForce() const
{
    return intervalStepForce;
}

inline RDouble Param_NSSolver::GetWallTemperature() const
{
    return wallTemperature;
}

inline int Param_NSSolver::GetChemicalFlag() const
{
    return nChemical;
}

inline int Param_NSSolver::GetWallMultiTemperature() const
{
    return wallMultiTemperature;
}

inline int Param_NSSolver::GetNSEquationNumber() const
{
    return nNSEquation;
}

inline int Param_NSSolver::GetLaminarNumber() const
{
    return nLaminar;
}

inline int Param_NSSolver::GetnEquation() const
{
    return nEquation;
}

inline int Param_NSSolver::GetTemperatureModel() const
{
    return nTemperatureModel;
}

inline RDouble Param_NSSolver::GetPreconCoefficient() const
{
    return kPreconCoeff;
}

inline int Param_NSSolver::GetPreconFarFieldBCMethod() const
{
    return preconFarfieldBCMethod;
}

inline int Param_NSSolver::GetIfIdealGasState() const
{
    return nIdealState;
}

inline RDouble Param_NSSolver::GetReferenceDimensionalVelocity() const
{
    return refDimensionalVelocity;
}

inline RDouble Param_NSSolver::GetAoA() const
{
	return AoA;
}

inline int Param_NSSolver::GetRoeEntropyFixMethod() const
{
    return roeEntropyFixMethod;
}

inline RDouble Param_NSSolver::GetRoeEntropyScale() const
{
    return roeEntropyScale;
}

inline RDouble Param_NSSolver::GetRoeEntropyFixCoef1() const
{
	return RoeEntropyFixCoef1;
}

inline RDouble Param_NSSolver::GetRoeEntropyFixCoef2() const
{
	return RoeEntropyFixCoef2;
}

inline RDouble Param_NSSolver::GetAngleOfSlide () const
{
    return angleSlide;
}

inline RDouble Param_NSSolver::GetPrandtlLaminar() const
{
    return prandtlLaminar;
}

inline RDouble Param_NSSolver::GetPrandtlTurbulence() const
{
    return prandtlTurbulence;
}

inline RDouble Param_NSSolver::GetoPrandtlLaminar() const
{
    return oPrandtlLaminar;
}

inline RDouble Param_NSSolver::GetoPrandtlTurbulence() const
{
    return oPrandtlTurbulence;
}

inline int Param_NSSolver::GetNChemicalRadius() const
{
    return nChemicalRadius;
}

inline int Param_NSSolver::GetNChemicalSource() const
{
    return nChemicalSource;
}

inline int Param_NSSolver::GetTimeIntegration() const
{
    return timeIntegration;
}

inline RDouble Param_NSSolver::GetNonDimensionalSutherlandTemperature() const
{
    return nonDimensionalSutherlandTemperature;
}

inline int Param_NSSolver::GetMgProlongationType() const
{
    return mgProlongationType;
}

inline RDouble* Param_NSSolver::GetPrimitiveVarFarfield() const
{
    return primitiveVarFarfield;
}

inline RDouble* Param_NSSolver::GetFullyCatalyticMassFraction() const
{
    return catalyticMassFraction;
}

inline int Param_NSSolver::GetNumberOfSpecies() const
{
    return numberOfSpecies;
}

inline int Param_NSSolver::GetIndexOfNitrogen() const
{
    //! The index of Nitrogen in the array of species, it indicates the species with the maximum mass fraction.
    return indexOfNitrogen;
}

inline int Param_NSSolver::GetIndexOfElectron() const
{
    //! The index of Electron in the array of species.
    return indexOfElectron;
}


inline int Param_NSSolver::GetnDensityModify() const
{
    return nDensityModify;
}

inline int Param_NSSolver::GetnDebug() const
{
    return nDebug;
}

inline int Param_NSSolver::GetmTT() const
{
    return mTT;
}

inline int Param_NSSolver::GetmTV() const
{
    return mTV;
}

inline int Param_NSSolver::GetmTE() const
{
    return mTE;
}

inline int Param_NSSolver::GetnEnergyRecycle() const
{
    return nEnergyRecycle;
}

inline int Param_NSSolver::GetTransitionType() const
{
    return transitionType;
}

inline int Param_NSSolver::GetnViscosityFluxSublevelModified() const
{
    return nViscosityFluxSublevelModified;
}

inline int Param_NSSolver::GetFlagOfEquilibriumGas() const
{
    return nEquilibriumGas;
}

inline int Param_NSSolver::GetThermoEnergyModelFlag() const
{
    return nTEnergyModel;
}

inline int Param_NSSolver::GetnTurblenceForChemical() const
{
    return nTurblenceForChemical;
}

inline RDouble Param_NSSolver::GetDensityMin() const
{
    return densityMin;
}

inline RDouble Param_NSSolver::GetdensityMinFactor() const
{
    return densityMinFactor;
}

inline int Param_NSSolver::GetnGradPrimtiveMethod() const
{
    return nGradPrimtiveMethod;
}

inline int Param_NSSolver::GetnInviscidFluxModify() const
{
    return nInviscidFluxModify;
}

inline int Param_NSSolver::GetSelfAdaptionSolveFlag() const
{
    return isSelfAdaptionSolve;
}

inline int Param_NSSolver::GetnQlLimitMethod() const
{
    return nQlLimitMethod;
}

inline int Param_NSSolver::GetSurfaceGradientMethod() const
{
    return nSurfGradMethod;
}

inline int Param_NSSolver:: GetRapidFlowfieldMethod() const
{
    return nRapidFlowfield;
}

inline int Param_NSSolver::GetSurfaceHeatingMonitor() const
{
    return nSurfHeatMonitor;
}

inline int Param_NSSolver::GetInitialPressureSteps() const
{
    return nInitPressureStep;
}

inline int Param_NSSolver::GetLocalCFLFlag() const
{
    return isUseLocalCFL;
}

inline RDouble Param_NSSolver::GetCharacteristicLength() const
{
    return knLength;
}

inline int Param_NSSolver::GetNonequilibriumConditionFlag() const
{
    return isUseNoneqCond;
}

inline int Param_NSSolver::GetIsPorousZone() const
{
    return isPorousZone;
}

inline const RDouble *Param_NSSolver::GetViscousResistanceCoeff() const
{
    return viscousResistanceCoeff;
}

inline const RDouble *Param_NSSolver::GetInertialResistanceCoeff() const
{
    return inertialResistanceCoeff;
}

inline RDouble Param_NSSolver::GetPorosity() const
{
    return porosity;
}

inline RDouble Param_NSSolver::GetDensitySolid() const
{
    return densitySolid;
}

inline RDouble Param_NSSolver::GetSpecificHeatSolid() const
{
    return cpSolid;
}

inline RDouble Param_NSSolver::GetHeatConductivitySolid() const
{
    return kSolid;
}
