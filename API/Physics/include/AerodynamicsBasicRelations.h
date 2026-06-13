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
//! @file      Aerodynamics_BasicRelations.h
//! @brief     Define some basic aerodynamic relations, e.g isentropic relations, shock relations.
//! @author    Zhang Jian, He Xin.

#pragma once
#include "Precision.h"

namespace PHSPACE
{

//! isentropic relationships
struct IsentropicRelations{
    //! T0/T
    static RDouble T0OverT(RDouble gama, RDouble mach);
    
    //! P0/P
    static RDouble P0OverP(RDouble gama, RDouble mach);
    
    //! Rho0/Rho
    static RDouble Rho0OverRho(RDouble gama, RDouble mach);
};

struct ShockRelations{
    //! Behind/Before shock temperature relation
    //! @param gama gama
    //! @param m1   mach before the shock
    static RDouble T2OverT1(RDouble gama, RDouble m1);

    //! Behind/Before shock pressure relation
    //! @param gama gama
    //! @param m1   mach before the shock
    static RDouble P2OverP1(RDouble gama, RDouble m1);

    //! Behind/Before shock density relation
    //! @param gama gama
    //! @param m1   mach before the shock
    static RDouble Rho2OverRho1(RDouble gama, RDouble m1);

    //! Compute mach behind the shock by mach before the shock
    //! @param gama gama
    //! @param m1   mach before the shock
    static RDouble MachBehindShock(RDouble gama, RDouble m1);
};

}
