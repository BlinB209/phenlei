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
//! @file      GPUCompGradientGGNodeNew.h
//! @brief     Compute gradient by GGNodeNew.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUCompGradientGGNode_NEW(const string d_gradient_field_proxy, const int d_gradient_var_index,
                                       const int nTotalNode, const int nTotalCell, const int nBoundFace,
                                       const int nTotalFace);
    void CallGPUSetVelocityNode(const int nTotalNode);

    __global__ void GPUSetVelocityNode(const int nTotalNode, const RFloat *qNode, RFloat *qVelNode);

    __global__ void GPUCompGradientGGNodeNewInit(const int index, const int nTotal, RFloat *dqdx, RFloat *dqdy,
                                                 RFloat *dqdz);

    __global__ void GPUCompGradientGGNodeNewInteriorFaceCal(
        const int nTotalNode, const int index, const int nTotalCell, const int nTotalFace, const int nBoundFace,
        const int *left_cell_of_face, const int *right_cell_of_face, const int *node_number_of_each_face,
        const int *face2node, const RDouble *nxs, const RDouble *nys, const RDouble *nzs, const RDouble *ns,
        const long long int *nodepos, const double *q_n, RFloat *dqdx, RFloat *dqdy, RFloat *dqdz);

    __global__ void GPUCompGradientGGNodeNewBoundaryFaceCal(const int nTotalNode, const int index, const int nTotalCell,
                                                            const int nBoundFace, const int *left_cell_of_face,
                                                            const int *node_number_of_each_face, const int *face2node,
                                                            const RDouble *nxs, const RDouble *nys, const RDouble *nzs,
                                                            const RDouble *ns, const long long int *nodepos,
                                                            const double *q_n, RFloat *dqdx, RFloat *dqdy, RFloat *dqdz,
                                                            const RFloat *q, const int *boundaryType,
                                                            const int d_SYMMETRY, const int d_INTERFACE);

    __global__ void GPUCompGradientGGNodeNewDiVol(const int index, const int nTotalCell, const int nTotal, RFloat *dqdx,
                                                  RFloat *dqdy, RFloat *dqdz, RDouble *vol);

    __global__ void GPUCompGradientGGNodeNewGhostSet(const int index, const int nTotalCell, const int nTotal,
                                                     const int nBoundFace, const int *left_cell_of_face, RFloat *dqdx,
                                                     RFloat *dqdy, RFloat *dqdz);
} //! namespace GPUKernels
