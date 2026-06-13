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
//! @file      Statistics.h
//! @brief     statistics that can be computed from fourier planes.
//! @author    ZhangZipei.

#pragma once
#include "TypeDefine.h"
#include "SpecDiffHybGrid.h"
//#include "SpecDiffHybSolver.h"
#include "Param_SpecDiffHybSolver.h"
#include "LIB_Macro.h"
#include "ActionKey.h"

using namespace std;

namespace PHSPACE
{
class StatisticsHIT
{
public:
    StatisticsHIT();
    ~StatisticsHIT();

public:
    void InitStatisticsData(SpecGrid *GridData);

    void GetStatisticsLocal(SpecGrid *GridData, Complex3D *complexU, Complex3D *complexV, Complex3D *complexW);

    void OutputStatisticRestart(SpecGrid *GridData);

    void OutputStatisticVisual(SpecGrid *GridData);

protected:
    virtual void InitStatisticsFileName();

    virtual bool IfStatisticExist();

    virtual void ReadStatisticsInformation();

    virtual void WriteStatisticsName();

    virtual void WriteStatisticInformation();

    const string GetStatisticVisualFileName();

    void FillActionKey(ActionKey *actkey, int action, const int level);

private:
    void AllocStatisticsdata(SpecGrid *GridData);

    void InitStatisticsDataValue(SpecGrid *GridData);

    void InputStatistics(SpecGrid *GridData);

    void ReadStatisticRestart(SpecGrid *GridData);

    void WriteStatisticRestart(SpecGrid *GridData);

    void WriteStatisticVisual(SpecGrid *GridData);

    void PrepareFormatedStatistic(SpecGrid *GridData, string &formatStatistic);

    void WriteSpectraVisual(SpecGrid *GridData);

    void WriteSpectraKxVisual(SpecGrid *GridData);

    void WriteSpectraKyVisual(SpecGrid *GridData);

    void WriteSpectraKzVisual(SpecGrid *GridData);

    void PrepareFormatedSpectraKx(SpecGrid *GridData, string &formatSpectra);

    void PrepareFormatedSpectraKy(SpecGrid *GridData, string &formatSpectra);

    void PrepareFormatedSpectraKz(SpecGrid *GridData, string &formatSpectra);

protected:
    string zGridPlot3dFile, zHalfGridPlot3dFile, statisticsNameFile, statisticsInformationFile, statisticsRestartFile, statisticsRestartTimeFile;

    int istats;

    bool ifInitStatistic;

public:
    int iz_half;

    RDouble1D *statisticAll, *statisticLocal, *statistic;

    RDouble2D *energyOfUUKx, *energyOfUUKy, *energyOfUUKz, *energyOfUUKxAll, *energyOfUUKyAll, *energyOfUUKzAll, *energyOfUUKxLocal, *energyOfUUKyLocal, *energyOfUUKzLocal;

    RDouble nutAverage;
};

class Statistics : public StatisticsHIT
{
public:
    Statistics();
    ~Statistics();

public:
    //int iz_half;//, naddnum;
    RDouble2D *statisticAll, *statisticLocal, *statistic;
    RDouble3D *energyOfUUKx, *energyOfUUKy, *energyOfUUKxAll, *energyOfUUKyAll, *energyOfUUKxLocal, *energyOfUUKyLocal;

    RDouble1D *nutAverage;

private: 

    RDouble tauWAverage;

public:
    void InitStatisticsData(SpecDiffHybGrid *GridData);

    void GetStatisticsLocal(SpecDiffHybGrid *GridData, Complex3D *, Complex3D *, Complex3D *, Complex4D *);

    void OutputStatisticRestart(SpecDiffHybGrid *GridData);

    void OutputStatisticVisual(SpecDiffHybGrid *GridData);

    RDouble GetTauWAverage();

    void SetTauWAverage(const RDouble &tauWAverageIn);

    bool GetIfInitStatistic();

    void SetIfInitStatistic(const bool &ifInitStatisticIn);

private:

    void OutputPlot3dGrid(SpecDiffHybGrid *GridData);

    void AllocStatisticsdata(SpecDiffHybGrid *GridData);

    void InitStatisticsDataValue(SpecDiffHybGrid *GridData);

    void InputStatistics(SpecDiffHybGrid *GridData);

    void ReadStatisticRestart(SpecDiffHybGrid *GridData);

    void ReadStatisticsSelected(SpecDiffHybGrid *GridData);

    void WriteStatisticRestart(SpecDiffHybGrid *GridData);

    void WriteStatisticVisual(SpecDiffHybGrid *GridData);

    void PrepareFormatedStatistic(SpecDiffHybGrid *GridData, string &formatStatistic);

    void WriteSpectraVisual(SpecDiffHybGrid *GridData);

    void WriteSpectraKxVisual(SpecDiffHybGrid *GridData);

    void WriteSpectraKyVisual(SpecDiffHybGrid *GridData);

    void PrepareFormatedSpectraKx(SpecDiffHybGrid *GridData, string &formatSpectra);

    void PrepareFormatedSpectraKy(SpecDiffHybGrid *GridData, string &formatSpectra);
};

}
