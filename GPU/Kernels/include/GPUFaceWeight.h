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
//! @file      GPUFaceWeight.h
//! @brief     Compute face weight.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUFaceWeight(const int nst, const int ned, const int nBoundFace);

    __global__ void GPUFaceWeight(const int nst, const int ned, const int *d_left_cell_of_face,
                                  const int *d_right_cell_of_face, const RDouble *d_xfc, const RDouble *d_yfc,
                                  const RDouble *d_zfc, const RDouble *d_xcc, const RDouble *d_ycc,
                                  const RDouble *d_zcc, RDouble *d_deltl, RDouble *d_deltr);

    __global__ void GPUFaceWeightOnBound(const int nst, const int nMid, RDouble *deltl, RDouble *deltr,
                                         const int *boundaryType, int bctype_in);
} //! namespace GPUKernels
