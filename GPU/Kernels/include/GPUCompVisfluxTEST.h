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
//! @file      GPUCompVisfluxTEST.h
//! @brief     Compute viscous flux.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
#include "Param_NSSolverUnstruct.h"
namespace GPUKernels
{
    void CallGPUComputeVisflux(Grid *gridIn, Param_NSSolverUnstruct *parameters, const int neqn, const int nst,
                               const int ned);

    __global__ void GPUCompVisfluxTEST_S1(
        const int nst, const int ned, const int nTotalCell, const int nBoundFace, int neqn, const int nchem,
        const int nm, const int nl, const int numberOfSpecies, const int ITT, const int nrokplus,
        const double skew_angle, const double oreynolds, const double SMALL, const double TINY,
        const int *left_cell_of_face, const int *right_cell_of_face, const RDouble *xfn, const RDouble *yfn,
        const RDouble *zfn, const RDouble *xfc, const RDouble *yfc, const RDouble *zfc, const RDouble *xcc,
        const RDouble *ycc, const RDouble *zcc, const RDouble *area, const RFloat *qpmv, const RFloat *t,
        const RFloat *dqdx, const RFloat *dqdy, const RFloat *dqdz, const RFloat *dtdx, const RFloat *dtdy,
        const RFloat *dtdz, const RFloat *prim, const RFloat *tm, const RDouble *deltl, const RDouble *deltr,
        const RFloat *kcp, const RFloat *mul, const RFloat *mut, const int len, const int *boundaryType,
        const int bctype_interface, const int bctype_solidface, const double coefseq, RFloat *flux);
}//！ namespace GPUKernels
