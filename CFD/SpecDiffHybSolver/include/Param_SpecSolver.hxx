inline int Param_SpecSolver::GetStatisticMethod() const
{
    return statisticMethod;
}

inline int Param_SpecSolver::GetStartStatisticStep() const
{
    return statisticMethod;
}

inline int Param_SpecSolver::GetMaxSimuStep() const
{
    return maxSimuStep;
}

inline int Param_SpecSolver::GetIntervalStepRes() const
{
    return intervalStepRes;
}

inline int Param_SpecSolver::GetIntervalStepStatistic() const
{
    return intervalStepStatistic;
}

inline int Param_SpecSolver::GetIntervalStepFlow() const
{
    return intervalStepFlow;
}

inline RDouble Param_SpecSolver::GetRefReNumber() const
{
    return refReNumber;
}

inline const std::string & Param_SpecSolver::GetViscousName() const
{
	return viscousName;
}

inline int Param_SpecSolver::GetViscousType() const
{
	return viscousType;
}

inline const std::string & Param_SpecSolver::GetFilterName() const
{
	return filterName;
}

inline RDouble Param_SpecSolver::GetSmagConstant() const
{
    return smagConstant;
}

inline RDouble Param_SpecSolver::GetTestFilterScale() const
{
    return testFilterScale;
}

inline RDouble Param_SpecSolver::GetTimeStep() const
{
    return timeStep;
}

/*
inline RDouble Param_SpecSolver::GetPressureGradient() const
{
    return pressureGradient;
}
*/

inline int Param_SpecSolver::GetSpectrumType() const
{
    return spectrumType;
}