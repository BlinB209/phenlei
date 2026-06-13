inline RDouble Param_CFDSolver::GetRefMachNumber() const
{
    return refMachNumber;
}

inline RDouble Param_CFDSolver::GetRefReNumber() const
{
    return refReNumber;
}

inline RDouble Param_CFDSolver::GetoRefReNumber() const
{
    return oRefReNumber;
}

inline int Param_CFDSolver::GetPlotFieldType() const
{
    return plotFieldType;
}

inline const RDouble * Param_CFDSolver::GetLowerPlotFieldBox() const
{
    return lowerPlotFieldBox;
}

inline const RDouble * Param_CFDSolver::GetUpperPlotFieldBox() const
{
    return upperPlotFieldBox;
}

inline int Param_CFDSolver::GetViscousType() const
{
    return viscousType;
}

inline const std::string & Param_CFDSolver::GetViscousName() const
{
    return viscousName;
}

inline bool Param_CFDSolver::IsViscous() const
{
    return (viscousType != 0);
}

inline int Param_CFDSolver::GetIsUnsteady() const
{
    return isUnsteady;
}

inline int Param_CFDSolver::GetIntervalStepRes() const
{
    return intervalStepRes;
}

inline int Param_CFDSolver::GetIntervalStepPlot() const
{
    return intervalStepPlot;
}

inline int Param_CFDSolver::GetIsCodeOfAleModel() const
{
    return isAle;
}

inline RDouble Param_CFDSolver::GetCFLStart() const
{
    return CFLStart;
}

inline RDouble Param_CFDSolver::GetPMax() const
{
    return pMax;
}

inline RDouble Param_CFDSolver::GetPMin() const
{
    return pMin;
}

inline RDouble Param_CFDSolver::GetDeltaMax() const
{
    return deltaMax;
}

inline RDouble Param_CFDSolver::GetMagnifyFactor() const
{
    return magnifyFactor;
}

inline RDouble Param_CFDSolver::GetReduceFactor() const
{
    return reduceFactor;
}

inline RDouble Param_CFDSolver::GetCFLEnd() const
{
    return CFLEnd;
}

inline RDouble Param_CFDSolver::GetFineCFL() const
{
    return fineCFL;
}

inline int Param_CFDSolver::GetCFLVaryStep() const
{
    return CFLVaryStep;
}

inline int Param_CFDSolver::GetNMGLevel() const
{
    return nMGLevel;
}

inline int Param_CFDSolver::GetNLUSGSSweeps() const
{
    return nLUSGSSweeps;
}

inline RDouble Param_CFDSolver::GetLUSGSTolerance() const
{
    return LUSGSTolerance;
}

inline int Param_CFDSolver::GetIntervalStepFlow() const
{
    return intervalStepFlow;
}

inline int Param_CFDSolver::GetIfLocalTimeStep() const
{
    return ifLocalTimeStep;
}

inline RDouble Param_CFDSolver::GetTorqueRefX() const
{
    return TorqueRefX;
}

inline RDouble Param_CFDSolver::GetTorqueRefY() const
{
    return TorqueRefY;
}

inline RDouble Param_CFDSolver::GetTorqueRefZ() const
{
    return TorqueRefZ;
}

//! Get reference gama vlaue
inline RDouble Param_CFDSolver::GetRefGama() const
{
    return refGama;
}

inline RDouble Param_CFDSolver::GetRefDimensionalPressure() const
{
    return refDimensionalPressure;
}

inline RDouble Param_CFDSolver::GetRefDimensionalDensity() const
{
    return refDimensionalDensity;
}

inline RDouble Param_CFDSolver::GetRefDimensionalVelocity() const
{
    return refDimensionalVelocity;
}

inline RDouble Param_CFDSolver::GetRefDimensionalSonicSpeed() const
{
    return refDimensionalSonicSpeed;
}

inline RDouble Param_CFDSolver::GetRefDimensionalTemperature() const
{
    return refDimensionalTemperature;
}

inline RDouble Param_CFDSolver::GetAusmpwPlusLimiter() const
{
    return AusmpwPlusLimiter;
}

inline RDouble Param_CFDSolver::GetTrTemperatureMinNonDim() const
{
    return trTemperatureMinNonDim;
}

inline int Param_CFDSolver::GetIfStartFromSteadyResults() const
{
    return ifStartFromSteadyResults;
}

inline int Param_CFDSolver::GetIfStaticsFlowField() const
{
    return ifStaticsFlowField;
}

inline int Param_CFDSolver::GetIsRestartChangeInflow() const
{
    return isRestartChangeInflow;
}

inline int Param_CFDSolver::GetIfStaticsReynoldsStress() const
{
    return ifStaticsReynoldsStress;
}

inline int Param_CFDSolver::GetStatisticMethod() const
{
    return statisticMethod;
}

inline int Param_CFDSolver::GetIsOverLapping() const
{
    return isOverset;
}

inline RDouble Param_CFDSolver::GetSkewnessAngle() const
{
    return skewnessAngle;
}
inline RDouble Param_CFDSolver::GetTurbSkewnessAngle() const
{
    return turbSkewnessAngle;
}
inline int Param_CFDSolver::GetnDiagonalModifiedTurb() const
{
    return nDiagonalModifiedTurb;
}

inline int Param_CFDSolver::GetWennSchemeFlag() const
{
    return isWennScheme;
}
inline int Param_CFDSolver::GetIfLowSpeedPrecon() const
{
    return ifLowSpeedPrecon;
}

inline int Param_CFDSolver::GetnNumberOfSpeedStep() const
{
    return nNumberOfSpeedStep;
}
inline int *Param_CFDSolver::GetSpeedVaryStep() const
{
    return speedVaryStep;
}
inline RDouble *Param_CFDSolver::GetSpeedVaryCoef() const
{
    return speedVaryCoef;
}
inline int Param_CFDSolver::GetnNumberOfCFLStep() const
{
    return nNumberOfCFLStep;
}
inline int *Param_CFDSolver::GetCFLVaryMultiStep() const
{
    return CFLVaryMultiStep;
}
inline RDouble *Param_CFDSolver::GetCFLVaryCoef() const
{
    return CFLVaryCoef;
}





