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
//! @file      Residual.h
//! @brief     It defines the residual computing, collection, dumping method.
//! @author    Bell, He Xin.

#pragma once
#include "Geo_StructGrid.h"
#include "Geo_UnstructGrid.h"

using namespace std;
namespace PHSPACE
{
//! Residual Initialize, called when task action.
void InitResidual();

//! Residual collection, called when task action.
void CollectResidual(ActionKey *actkey);

//! Residual dump, called when task action.
void PostDumpResidual(ActionKey *actkey);

class ActionKey;
class Residual;

//! @brief ResidualGlobal defines the residual operation globally,
//!        including residual initialization, collection(MPI), dumping.
class ResidualGlobal
{
private:
    //! Residual on each local zone.
    Residual **residualOnEachZone;

    //! Global residual after collection, which is ONLY defined on server processor.
    RDouble *residualNormal;

    //! Number of residual variables.
    int numberOfResidualVariables;

    //! Control if print residual to window.
    bool IsFirstStep;

public:
    ResidualGlobal();
    ~ResidualGlobal();

public:
    //! Residual Initialize globally.
    void InitResidual();

    //! Residual collection globally.
    void CollectResidual(ActionKey *actkey);

    //! Residual dump globally, ONLY execute it on server processor.
    void DumpResidual(ActionKey *actkey, int outnstep);

private:
    void ComputeResidualNorm(RDouble *localNorm);
    void PrepareFormatedResidual(int outnstep, string &formatRes);
    void PrintResidualtoWindow(int outnstep, string &formatRes);
    void OutResResidualNorm(fstream &file, string &formatRes, std::ostringstream &oss);
    void FreeMemory();
};

//! @brief Residual defines the residual on each local zone in the current processor.
class Residual
{
private:
    //! Zone index.
    int zoneID;
    int npoint;
    int numberOfResidualVariables;
    RDouble *x, *y, *z, *vol;
    int    *i, *j, *k;
    RDouble *maximumResidual, *averageResidual;
public:
    Residual();
    ~Residual();
public:
    //! Set the number of residual variables.
    void SetNvar(int nvar);

    //! Set the zone index.
    void SetZoneID(int zoneID) {this->zoneID = zoneID;}

    //! Get the number of residual variables.
    int  GetNvar() const {return this->numberOfResidualVariables;}

    //! Get the zone index.
    int  GetZoneID() const {return this->zoneID;}

    //! Get the residual on this zone.
    RDouble *GetAverageResidual() {return this->averageResidual;}
};

class MaxResidual
{
private:
    int processorOfMaxResidual;
    RDouble x, y, z;
    int maxResVariableIndex;
    RDouble maximumResidual;
    DataContainer *data;

public:
    MaxResidual();
    ~MaxResidual();

public:
    void Init();

    int GetProcessorOfMaxResidual();
    RDouble GetMaxRes();
    RDouble GetMaxResCoorX();
    RDouble GetMaxResCoorY();
    RDouble GetMaxResCoorZ();
    int GetMaxResVariableIndex();

    void SetProcessorOfMaxResidual(const int &processorOfMaxResidualIn);
    void SetMaxRes(const RDouble &maxResIn);
    void SetMaxResCoorX(const RDouble &xIn);
    void SetMaxResCoorY(const RDouble &yIn);
    void SetMaxResCoorZ(const RDouble &zIn);
    void setMaxResVariableIndex(const int &maxResVariableIndexIn);

    void ServerCollection();

private:
    void InitMamory();
    void DeleteMamory();
};

void InitMaxResidual();
void EndMaxResidual();

void ComputeResidualonGrid(UnstructGrid *gridIn, ActionKey *actkey, RDouble **dq, int nEquation);
void ComputeResidualonGrid(StructGrid   *gridIn, ActionKey *actkey, RDouble4D &dq, int nEquation);

}
