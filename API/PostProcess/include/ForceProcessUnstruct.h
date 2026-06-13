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
//! @file      ForceProcessUnstruct.h
//! @brief     UnstructGrid Flow result post-processing
//! @author    Baka.

#pragma once
#include "Force.h"
#include "Geo_Grid.h"
#include "PostProcessData.h"
#include "Constants.h"

namespace PHSPACE
{

class ForceProcessUnstruct
{
public:
    Force globalAerodynamicForce;

public:
    ForceProcessUnstruct();
    ~ForceProcessUnstruct();

public:
    void Run();

private:
    void InitForce();
    void CollectionForce();
    void AirForceCoef(int iZone);
    void AirForceCoefParts(int iZone, int partID, int Coordinate = GlobalCoordinate);
    void CompForce();
    void DumpForce();

    void InitValueField();
    void CompValueField();
    void ComputeGamaAndTemperature(Grid *gridIn);
    void ComputeViscousLaminarCoefficient(Grid *gridIn);
    void GetGradientField(Grid *gridIn);

private:
    string filename;
    DataContainer *CollectData;

};

}