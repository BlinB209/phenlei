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
//! @file      GPUFixBCNodeVar.h
//! @brief     Compute node valur for boundary face.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUFixBCNodeVarByCompNodeVar(const int index, const int nTotalNode, const int nBoundFace,
                                          const int nTotalCell, const int *left_cell_of_face,
                                          const int *right_cell_of_face, const int *face2node,
                                          const int *node_number_of_each_face, const int *boundaryType,
                                          const long long int *nodepos, RFloat *q, RFloat *q_n, int *n_count,
                                          const int bctype_in, const bool twoside);

    __global__ void GPUFixBCNodeVarInit(const int nBoundFace, const int nTotalCell, const int *left_cell_of_face,
                                        const int *face2node, const int *node_number_of_each_face,
                                        const int *boundaryType, const long long int *nodepos, RFloat *q_n,
                                        int *n_count, const int bctype_in);

    __global__ void GPUFixBCNodeVarCal(const int index, const int nBoundFace, const int nTotalCell,
                                       const int *left_cell_of_face, const int *face2node,
                                       const int *node_number_of_each_face, const int *boundaryType,
                                       const long long int *nodepos, RFloat *q, RFloat *q_n, int *n_count,
                                       const int bctype_in, const bool twoside);

} //! namespace GPUKernels
