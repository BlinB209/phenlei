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
//! @file      GPUKernelTestAPI.h
//! @brief     Kernel test function for GPU results.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#ifndef __GPUKernelTestAPI_H_
#define __GPUKernelTestAPI_H_
#include "Precision.h"
using namespace PHSPACE;
//! test the GPUVencatLimiter7Loop1
void TestGPUVencatLimiter7Eps(RFloat *epsCell, RFloat *d_epsCell, const int nTotalCell);
void TestGPUVencatLimiter7Boun(RFloat *limit, int *left_cell_of_face, const int iVariable, const int nBoundFace,
                               const int nTotal);
#endif
