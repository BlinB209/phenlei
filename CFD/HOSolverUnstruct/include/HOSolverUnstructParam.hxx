#pragma once

inline const string & HOSolverUnstructParam::GetRestartNSFile() const
{
	return restartNSFile;
}

inline const string & HOSolverUnstructParam::GetResSaveFile() const
{
	return resSaveFile;
}

inline const string & HOSolverUnstructParam::GetFlowFieldTecFile() const
{
	return flowFieldTecFile;
}

inline int HOSolverUnstructParam::GetIntervalStepForce() const
{
	return intervalStepForce;
}

inline RDouble HOSolverUnstructParam::GetWallTemperature() const
{
	return wallTemperature;
}

inline int HOSolverUnstructParam::GetNTemperatureModel() const
{
	return nTemperatureModel;
}

inline int HOSolverUnstructParam::GetIfLowSpeedPrecon() const
{
	return ifLowSpeedPrecon;
}

inline double HOSolverUnstructParam::GetAoA() const
{
	return AoA;
}
inline RDouble HOSolverUnstructParam::GetRoeEntropyFixCoef1() const
{
	return RoeEntropyFixCoef1;
}

inline RDouble HOSolverUnstructParam::GetRoeEntropyFixCoef2() const
{
	return RoeEntropyFixCoef2;
}

inline RDouble HOSolverUnstructParam::GetAngleOfSlide () const
{
    return angleSlide;
}

inline RDouble HOSolverUnstructParam::GetR() const
{
    return R;
}