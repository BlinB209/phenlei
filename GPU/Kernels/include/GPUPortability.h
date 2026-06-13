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
//! @file      GPUPortability.h
//! @brief     Atomic kernel functions.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include <cuda.h>
#include "Precision.h"
using namespace PHSPACE;
namespace GPUKernels
{
    __device__ void atomicMax(double *addr, double val);
    //!__device__ void atomicMax(RFloat *addr, RFloat val);
    __device__ void atomicMin(double *addr, double val);
//!__device__ void atomicMin(RFloat *addr, RFloat val);
#if __CUDA_ARCH__ < 600
    __device__ double atomicAddTest(double *address, double val);
#endif
    //!__device__ double kernelMAXDOUBLE(double a, double b);
    __device__ double kernelMAXDOUBLE(RFloat a, RFloat b);
    //!__device__ double kernelMINDOUBLE(double a, double b);
    __device__ double kernelMINDOUBLE(RFloat a, RFloat b);
} //! namespace GPUKernels
