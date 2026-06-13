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
//! @file      GPUGetVisFaceValue.h
//! @brief     Get NS viscous face value.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
#include "Param_NSSolverUnstruct.h"
namespace GPUKernels
{
    void CallGPUGetVisFaceValue(Grid *gridIn, Param_NSSolverUnstruct *parameters, const int neqn, const int nst,
                                const int ned);
    __global__ void GPUGetVisFaceValueMulMut(const int nst, const int ned, const int nBoundFace, const int nTotalCell,
                                             const int *left_cell_of_face, const int *right_cell_of_face,
                                             const RFloat *visl, const RFloat *vist, RFloat *mul, RFloat *mut);
    __global__ void GPUGetVisFaceValuePrimTm(const int nst, const int ned, const int nBoundFace, const int nTotalCell,
                                             const int neqn, const int nchem, const int len, const int ITT,
                                             const int *left_cell_of_face, const int *right_cell_of_face,
                                             const RFloat *t, const RFloat *q, const RFloat *mul, const RFloat *mut,
                                             const double oprl, const double oprt, double gama0, RFloat m2,
                                             RFloat *prim, RFloat *tmid, RFloat *kcp);
    __global__ void GPUGetVisFaceValuePrimTm_S1(const int nst, const int ned, const int nBoundFace,
                                                const int nTotalCell, const int neqn, const int nchem, const int len,
                                                const int ITT, const int *left_cell_of_face,
                                                const int *right_cell_of_face, const RFloat *t, const RFloat *q,
                                                const RFloat *mul, const RFloat *mut, const double oprl,
                                                const double oprt, double gama0, RFloat m2, RFloat *prim, RFloat *tmid,
                                                RFloat *kcp);
    __global__ void GPUGetVisFaceValuePrimTm_S2(const int nst, const int nMid, const int ned, const int nBoundFace,
                                                const int nTotalCell, const int neqn, const int nchem, const int len,
                                                const int ITT, const int *left_cell_of_face,
                                                const int *right_cell_of_face, const RFloat *t, const RFloat *q,
                                                const RFloat *mul, const RFloat *mut, const double oprl,
                                                const double oprt, double gama0, RFloat m2, RFloat *prim, RFloat *tmid,
                                                RFloat *kcp);
    __global__ void GPUGetVisFaceValueOnBound(const int nst, const int nMid, const int nBoundFace, const int nTotalCell,
                                              const int neqn, const int len, const int *boundaryType,
                                              const int *left_cell_of_face, const int *right_cell_of_face,
                                              const RFloat *q, const RDouble *xtn, const RDouble *ytn,
                                              const RDouble *ztn, RFloat *prim, RFloat *tmid, RFloat tw,
                                              const int bctype_interface, const int bctype_solid);
} //! namespace GPUKernels
