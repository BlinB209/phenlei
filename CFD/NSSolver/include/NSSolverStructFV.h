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
//! @file      NSSolverStructFV.h
//! @brief     NS solver for struct grid with finite volume method.
//! @author    xxx

#pragma once
#include "NSSolverStruct.h"
#include "Geo_Grid.h"
#include "TypeDefine.h"

namespace PHSPACE
{
class NSSolverStructFV : public NSSolverStruct
{
public:
    NSSolverStructFV();
    ~NSSolverStructFV();
    void InviscidFlux(Grid* gridIn) override;

private:
    void GKSFluxStruct(Grid* gridIn);
    void ComputeGKSFluxCore(
        RDouble rhoL, RDouble uL, RDouble vL, RDouble wL, RDouble eL,
        RDouble rhoR, RDouble uR, RDouble vR, RDouble wR, RDouble eR,
        RDouble nx, RDouble ny, RDouble nz,
        RDouble* flux
    );
};
}
