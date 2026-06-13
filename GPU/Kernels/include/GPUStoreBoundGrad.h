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
//! @file      GPUStoreBoundGrad.h
//! @brief     Compute gradient for boundary.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUStoreBoundGrad(const int nTotalCell, const int nBoundFace, const int nTVar);

    __global__ void GPUStoreBoundGrad(const int nTotalCell, const int nBoundFace, const int nTVar,
                                      const int *left_cell_of_face, RFloat *d_dqdx_proxy, RFloat *d_dqdy_proxy,
                                      RFloat *d_dqdz_proxy, RFloat *d_bdqx, RFloat *d_bdqy, RFloat *d_bdqz);
} // namespace GPUKernels
