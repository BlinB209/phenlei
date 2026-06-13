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
//! @file      GPUSpectrumRadius.h
//! @brief     Spectrum radius for NS solver.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUInitSpec(const int nTotal);
    void CallGPUAddInvSpec(const int nTotalCell);
    void CallGPUAddVisSpec(const int nTotalCell);

    __global__ void GPUInitSpec(const int nTotal, RFloat *spec, const RFloat *dt);
    __global__ void GPUAddInvSpec(const int nTotalCell, RFloat *spec, const RFloat *invSpectrumRadius);
    __global__ void GPUAddVisSpec(const int nTotalCell, RFloat *spec, const RFloat *visSpectrumRadius);
} //! namespace GPUKernels
