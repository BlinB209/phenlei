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
//! @file      SpecSolver.h
//! @brief     Solver of Fourier-Spectral method.
//! @author    Zhang Zipei.

#pragma once

#include "TypeDefine.h"
#include "Complex.h"
#include "LIB_Macro.h"
#include "Param_SpecSolver.h"
#include "SpecDiffHybGrid.h"
#include "Spectrum.h"
#include "SpecGrid.h"
#include "RescaleField.h"
#include "Statistics.h"
#include "LES.h"

using namespace std;

namespace PHSPACE
{
class Param_SpecSolver;
//class SpecDiffHybGrid;

class SpecSolver
{
public:
    SpecSolver();
    ~SpecSolver();

protected:
    RDouble1D *nonlinearTermCoef, *pressureGradientTermCoef;
    Complex4D *convectiveTerm, *complexLaplacianVelocityTerm;
    Complex5D *complexPressureGradientTerm;

    Complex4D *complexViscousTerm, *complexNonlinearTerm, *complexVelocityGradient, *complexVelocitySquare;
    RDouble4D *realVelocityGradient, *realVelocitySquare;

    Complex3D *complexU, *complexV, *complexW, *complexP, *complexDivergenceU;

    RDouble3D *realnut3d;

    RDouble1D *ek3DInit, *ek3D;

    string nonlinearTermFile, nonlinearTermBackupFile;
    string pressureGradientTermFile, pressureGradientTermBackupFile;
    string laplacianVelocityTermFile, laplacianVelocityTermBackupFile;
    string timeInformationRestartFile;

    string velocityRestartFile, velocityBackupFile;
    string pressureRestartFile, pressureBackupFile;
    string velocityNameFile, fieldNameFile;

    bool IsConvectiveTermFailed;

private:
    RDouble realnut;

public:
    LIB_EXPORT Param_SpecSolver * GetControlParameters() const;

    void Run();

    virtual const string GetResidualFileName();

    SpecGrid *GetGrid();

private:
    void Init();

    void AllocTimeDependentData();

    void InitGlobal();

    //void ReadConvectiveTerm();
    //void ReadTimeDependentInformation();

    void Solve();
    //void RescaleVelocityField(const int);
    //void GetVelocityGradient();
    //void GetNonlinearTerm();
    //void GetViscousTerm();
    //void InitFieldData();

    //void AllocFieldData();
    //void GetInitField();
    //void InitFieldFileName();
    //void WriteVelocityName();
    //void InitProfile();
    //void InitFromOldField();
    //void OutputFieldRestartFile(const int);
    //void GetAssumedSpectrum( const double, const double, double & );
    //void GetAssumedSpectrumNearWall( const double, const double, double & );
    //void ReadinVelocityRestartFile();
    //void ReadinPressureRestartFile();
    //void WriteVelocityRestartFile(const string);
    //void WritePressureRestartFile(const string);
    //void SetInitTimeDataCoef();
    //void SetTimeDataCoef();
    //void VanishingDIVproj();
    //void OutputFieldName();
    //void GetStatistics();
    //void OutputResidual(const int);
    //void OutputRestartData(const int);
    //void OutputOldTimeData(const int);
    //void WriteConvectiveTerm(const string);
    //void WriteTimeDependentInformation();
    //void mode00fix();
    RescaleFieldHIT * GetRescaleFieldData();

    void SetRescaleFieldData(RescaleFieldHIT *RescaleFieldData);

    Spectrum * GetSpectrumData();

    void SetSpectrumData(Spectrum *SpectrumData);

    LESHIT * GetLESData();

    void SetLESData(LESHIT *LESData);

    StatisticsHIT *GetStatisticData();

    void SetStatisticData(StatisticsHIT *StatisticData);

    void NormVector(RDouble *vec1, RDouble *vec2, int ndim);

    void GetAssumedSpectrum(RDouble modk, RDouble &ek);

    void GetCBCSpectrum(RDouble modk, RDouble &ek);

    static void Splint(RDouble *xa, RDouble *ya, RDouble *y2a, int n, RDouble x, RDouble &y);

    void InitPressure();

    void GetEkAll(RDouble1D *ek3D, int nMax);

    void GetInitSpectrum();

    void GetGradientVector(PHComplex *var, PHComplex *dVardX, PHComplex *dVardY, PHComplex *dVardZ);

    void GetDivergence(PHComplex *vecX, PHComplex *vecY, PHComplex *vecZ, Complex3D *divVar);

protected:
    virtual void FreeControlParameters();

    virtual void InitControlParameters();

    virtual bool JudgeIfRestart();

    void SetGrid(SpecGrid *GridData);

    void InitSpectrumData();

    virtual void InitGridData();

    virtual void InitTimeDependentData();

    virtual void InitTimeDependentFileName();

    virtual void SetInitTimeDataCoef();

    virtual void SetTimeDataCoef();

    virtual void InitFieldData();

    virtual void InitFieldFileName();

    virtual void AllocFieldData();

    virtual void GetInitField();

    virtual void WriteVelocityName();

    virtual void InitProfile();

    virtual void InitLESData();

    virtual void InitStatisticsData();

    virtual void InitRescaleData();

    virtual void InitFromOldField();

    virtual void ReadinVelocityRestartFile();

    virtual void ReadinPressureRestartFile();

    virtual void OutputFieldRestartFile();

    virtual void WriteVelocityRestartFile(const string velocityFileName);

    virtual void WritePressureRestartFile(const string pressureFileName);

    virtual void InputOldTimeData();

    virtual void ReadConvectiveTerm();

    virtual void ReadTimeDependentInformation();

    virtual void OutputFieldName();

    virtual void GetRHS();

    virtual void GetVelocityGradient();

    virtual void GetNonlinearTerm();

    virtual void RescaleVelocityField(const int nstep);

    virtual void GetStatistics();

    virtual void OutputResidual();

    void FillActionKey(ActionKey *actkey, int action, const int level);

    virtual void DumpRes(ActionKey *actkey);

    virtual void PrepareFormatedResidual(string &formatRes);

    virtual void OutputRestartData();

    virtual void OutputOldTimeData();

    virtual void WriteConvectiveTerm(const string);

    virtual void WriteTimeDependentInformation();

    virtual void OutputStatisticVisual();

    virtual void DeAllocateGlobalVariables();

private:
    class Spectrum * SpectrumData;

    class SpecGrid * GridData;

    class RescaleFieldHIT * RescaleFieldData;

    class StatisticsHIT * StatisticData;

    class LESHIT * LESData;

protected:
    Param_SpecSolver *controlParameters;

};

#include "SpecSolver.hxx"

}