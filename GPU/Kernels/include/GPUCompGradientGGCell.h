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
//! @file      GPUCompGradientGGCell.h
//! @brief     Compute gradient by GGCell.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUCompGradientGGCell(const string d_gradient_field_proxy, const int d_index, const int nTotalCell,
                                   const int nBoundFace, const int nTotalFace);

    __global__ void GPUCompGradientGGCellInit(const int index, const int nTotal, RFloat *dqdx, RFloat *dqdy,
                                              RFloat *dqdz);
    __global__ void GPUCompGradientGGCellBoundaryFaceCal(const int index, const int nTotalCell, const int nBoundFace,
                                                         const int *left_cell_of_face, const RDouble *nxs,
                                                         const RDouble *nys, const RDouble *nzs, const RDouble *ns,
                                                         RFloat *q, RFloat *dqdx, RFloat *dqdy, RFloat *dqdz);
    __global__ void GPUCompGradientGGCellInteriorFaceCal(const int index, const int nTotalCell, const int nTotalFace,
                                                         const int nBoundFace, const int *left_cell_of_face,
                                                         const int *right_cell_of_face, const RDouble *nxs,
                                                         const RDouble *nys, const RDouble *nzs, const RDouble *ns,
                                                         RFloat *q, RFloat *dqdx, RFloat *dqdy, RFloat *dqdz);
    __global__ void GPUCompGradientGGCellDiVol(const int index, const int nTotalCell, const int nTotal, RFloat *dqdx,
                                               RFloat *dqdy, RFloat *dqdz, RDouble *vol);
    __global__ void GPUCompGradientGGCellGhostSet(const int index, const int nTotalCell, const int nTotal,
                                                  const int nBoundFace, const int *left_cell_of_face, RFloat *dqdx,
                                                  RFloat *dqdy, RFloat *dqdz);
} //! namespace GPUKernels
