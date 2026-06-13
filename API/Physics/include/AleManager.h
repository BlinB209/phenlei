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
//! @file      AleManager.h
//! @brief     Explain this file briefly.
//! @author    He Kun.

#pragma once
#include "Precision.h"
#include "Force.h"
#include "Pointer.h"
#include "AMRDef.h"
#include "Data_ParamFieldSuite.h"

using namespace std;

namespace PHSPACE
{
class SixDofSolverManager;
class MovingMeshManager;
class DeformingSolverManager;
class AleForceManager;
class Grid;
class Data_ParamFieldSuite;

class AleManager
{
public:
    AleManager();
    ~AleManager();

private:
    AleForceManager *aleForceManager;
    SixDofSolverManager *sixDofSolverManager;
    DeformingSolverManager *deformingManager;
    MovingMeshManager *movingMeshManager;

public:
    void Initialize();
    void Allocate();
    void Restart();
    void Run();
    void UpdateAleForce();
    void Post();
    void Dump();
    void DumpRestart();

protected:
    void SychonizeAleTime();
    void InitializeDimensionalReferenceTime();
    fstream *currentAleFile;

public:
    fstream & GetCurrentAleFile() { return *currentAleFile; };
    AleForceManager        * GetAleForceManager()        { return aleForceManager;        };
    SixDofSolverManager    * GetSixDofSolverManager()    { return sixDofSolverManager;    };
    DeformingSolverManager* GetDeformingSolverManager()  { return deformingManager; };
    void ComputeMetrics();
    void CommunicateCellCenterData();
    void ComputeGridFaceVelocity();
};

//! The global call of Ale.
void CreateAleManager();
void FreeAleManager();
void SolveSingleInnerIterationStepForAleSolvers();
void PostprocessAfterInnerIterationStepForAleSolvers();

AleManager * GetAleManager();

void BackupAllZonesOldGrid();

//! The backup file of Ale.
bool IsReadAleFile();
void SetIsReadAleFile(bool isReadAleFileIn);
fstream & GetCurrentAleFile();

int GetIntegerParameterFromDataBase(int bodyIndex, const string& parameterName);
int GetIntegerParameterFromDataBaseIfExist(int bodyIndex, const string& parameterName, int defaultValve);

RDouble GetRDoubleParameterFromDataBase(int bodyIndex, const string& parameterName);

string GetStringParameterFromDataBase(int bodyIndex, const string& parameterName);

RDouble GetRDoubleParameterFromDataBaseIfExist(int bodyIndex, const string &parameterName, RDouble defaultValue);

int * GetIntegerArrayFromDataBase(const string& nameOfIntegerArray, int numberOfElements);
void GetRDoubleVectorFromDataBase(vector< RDouble >& parameter, int bodyIndex, const string& parameterName, int numberOfElements);

RDouble GetAleTimeStep();

bool IsImplicitAleMethod();
bool IsExplicitAleMethod();

void DeAllocateOversetStorage();

void AllocateOversetStorage();

}