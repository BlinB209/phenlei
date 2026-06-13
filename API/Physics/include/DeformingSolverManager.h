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
//! @file      DeformingSolverManager.h
//! @brief     Explain this file briefly.
//! @author    He Kun.

#pragma once
#include "AleManager.h"
#include "Force.h"
#include "SixDofManager.h"

namespace PHSPACE
{

namespace MORPHINGMANNER
{
    const int NOMPRPHING  = 1;
    const int FISH        = 2;
    const int INSECT      = 3;
    const int BIRD        = 4;
    const int MISSILE     = 5;
    const int OSCILLATE   = 6;
    const int MORPHINGAIRCRAFTFOLD = 7;
    const int MORPHINGAIRCRAFTTAILEDGEUP = 8;
}

class MassCharacter;

//! This parameter is used to determine the motion of the grid.
class DeformingParameter
{
public:
    DeformingParameter();
    ~DeformingParameter();
protected:
    int numberOfParameters;
    vector< RDouble > parameters;
public:
    DeformingParameter & operator = (const DeformingParameter & rhs);
    bool operator == (const DeformingParameter & rhs);
public:
    void SetNumberOfParameters(int numberOfParameters) { this->numberOfParameters = numberOfParameters; };
    int  GetNumberOfParameters() { return numberOfParameters; };

    vector< RDouble > & GetParameters()   { return parameters; };
    void SetParameters(PHVectorRDouble1D& parameterVector);
};

DeformingParameter * GetDeformingParameter();

class DeformingSolverManager
{
public:
    DeformingSolverManager();
    ~DeformingSolverManager();

public:
    void Initialize();
    void Restart();
    void Run();
    void Post();
    void Dump();
    void DumpRestart();
};

}