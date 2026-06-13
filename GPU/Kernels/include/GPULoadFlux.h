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
//! @file      GPULoadFlux.h
//! @brief     Functions for load NS flux.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPULoadFlux(Grid *grid_in, int nst, int ned);

    __global__ void GPULoadFluxOnBound(const int nst, const int nMid, const int nTotalCell, const int nBoundFace,
                                       const int nl, const int len, const int *left_cell_of_face, const RFloat *flux,
                                       RFloat *res);
    __global__ void GPULoadFluxInterior(const int nst, const int nMid, const int ned, const int nTotalCell,
                                        const int nBoundFace, const int nl, const int len, const int *left_cell_of_face,
                                        const int *right_cell_of_face, const RFloat *flux, RFloat *res);
    __global__ void GPULoadFluxInteriorOptSep(const int equationID, const int nst, const int nMid, const int ned,
                                              const int nTotalCell, const int nBoundFace, const int nl, const int len,
                                              const int *left_cell_of_face, const int *right_cell_of_face,
                                              const RFloat *flux, RFloat *res);
} //！ namespace GPUKernels
