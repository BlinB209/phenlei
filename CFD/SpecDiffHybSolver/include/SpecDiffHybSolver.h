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
//! @file      SpecDiffHybSolver.h
//! @brief     Solver of Spectral-Difference hybrid method.
//! @author    Zhang Zipei.

#pragma once

#include "TypeDefine.h"
#include "LIB_Macro.h"
#include "Complex.h"
#include "ActionKey.h"
#include "SpecSolver.h"

using namespace std;

namespace PHSPACE
{
class Param_SpecDiffHybSolver;

class SpecDiffHybSolver : public SpecSolver
{
public:
    SpecDiffHybSolver();
    ~SpecDiffHybSolver();

private:

    RDouble1D *realnut, *realnut_00;

    string screenOutFile;

public:
    void Run();

    void CreateH5VelocityRestartFile(ActionKey *actkey);

    void DumpVelocityRestartH5(ActionKey *actkey);

    void ReadVelocityRestartH5(ActionKey *actkey);

    LIB_EXPORT Param_SpecDiffHybSolver * GetControlParameters() const;

    LIB_EXPORT void InitControlParameters();

private:
    void Init();

    void InitGridData();

    void InitTimeDependentData();

    void InitBoundaryData();

    void InitCompactDifferenceFirstDerivativeData();

    void InitCompactDifferenceSecondDerivativeData();

    void AllocTimeDependentData();

    void InitGlobal();

    void InitRescaleData();

    void InitPoissonSolverData();

    void InitCorrectWallDivergenceData();

    void ReadConvectiveTerm();

    void ReadTimeDependentInformation();

    void InitStatisticsData();

    void InitLESData();

    void GetRHS();

    void Solve();

    void RescaleVelocityField(const int);

    void GetVelocityGradient();

    void GetNonlinearTerm();

    void GetViscousTerm();

    void AllocFieldData();

    void GetInitField();

    void InitProfile();

    void InitFromOldField();

    void OutputFieldRestartFile();

    void GetAssumedSpectrum( const double, const double, double & );

    void GetAssumedSpectrumNearWall( const double, const double, double & );

    void ReadinVelocityRestartFile();

    void ReadinPressureRestartFile();

    void WriteVelocityRestartFile(const string);

    void WriteVelocityRestartFile_HDF5(const string);

    void ReadVelocityRestartFile_HDF5(const string);

    void WritePressureRestartFile(const string);

    void VanishingDIVproj();

    void GetStatistics();

    void PrepareFormatedResidual(string &formatRes);

    void DumpRes(ActionKey *actkey);

    void OutputRestartData();

    void OutputStatisticVisual();

    void OutputOldTimeData();

    void WriteConvectiveTerm(const string);

    void WriteTimeDependentInformation();

    void mode00fix();

    SpecDiffHybGrid *GetGrid();

    void SetGrid(SpecDiffHybGrid *GridData);

    RescaleField * GetRescaleFieldData();

    void SetRescaleFieldData(RescaleField *RescaleFieldData);

    LES * GetLESData();

    void SetLESData(LES *LESData);

    Statistics *GetStatisticData();

    void SetStatisticData(Statistics *StatisticData);

    void DeAllocateGlobalVariables();

private:
    class SpecDiffHybGrid * GridData;

    class ExplicitDifferenceBoundary * DiffBoundaryData;

    class CompactDifferenceFirstDerivative * CompactDifferenceFirstDerivativeData;

    class CompactDifferenceSecondDerivative * CompactDifferenceSecondDerivativeData;

    class PoissonSolver * PoissonSolverData;

    class CorrectWallDivergence * CorrectWallDivergenceData;

    class RescaleField * RescaleFieldData;

    class Statistics *StatisticData;

    class LES *LESData;

protected:
    Param_SpecDiffHybSolver *controlParameters;
};

#include "SpecDiffHybSolver.hxx"

}