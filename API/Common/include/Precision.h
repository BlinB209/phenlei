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
//! @file      Precision.h
//! @brief     Explain this file briefly.
//! @author    He Xin, Bell.

#pragma once
#include<stdio.h>
#include<mpi.h>

namespace PHSPACE
{
#ifdef PH_DOUBLE_PREC
    typedef double RDouble;
    typedef double RFloat;

    const RDouble LARGE = 1.0e40;
    const RDouble SMALL = 1.0e-40;
    const RDouble TINY  = 1.0e-40;
    const RDouble EPSILON = 1.0e-12;

    #define PH_MPI_DOUBLE MPI_DOUBLE
#else
    typedef float RDouble;
    typedef float RFloat;

    const RDouble LARGE = 1.0e36f;
    const RDouble SMALL = 1.0e-36f;
    const RDouble TINY  = 1.0e-36f;
    const RDouble EPSILON = 1.0e-8f;

    #define PH_MPI_DOUBLE MPI_FLOAT
#endif

typedef size_t uint_t;
typedef unsigned long long uint_long;

#ifdef WIN32
    #ifdef _WIN64
        typedef long long int_t;
    #else
        typedef int int_t;
    #endif
#else
    typedef long long int_t;
#endif

}