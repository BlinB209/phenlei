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
//! @file      GPUCompNodeVar.h
//! @brief     Compute nodevalue.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUCompNodeVarByGradient(const string d_gradient_field_proxy, const int d_gradient_var_index,
                                      const int nTotalNode, const int nBoundFace, const int nTotalFace,
                                      const int nTotalCell);

    __global__ void GPUNcountQnInit(const int nTotalNode, double *q_n, int *n_count);
    __global__ void GPUBoundaryFaceNodeCal(const int index, const int nBoundFace, const int nTotalCell,
                                           const int *left_cell_of_face, const int *node_number_of_each_face,
                                           const int *face2node, const long long int *nodepos, RDouble *x, RDouble *y,
                                           RDouble *z, RDouble *xcc, RDouble *ycc, RDouble *zcc, double *q_n,
                                           int *n_count, RFloat *q);

    __global__ void GPUInteriorFaceNodeCal(const int index, const int nBoundFace, const int nTotalCell,
                                           const int nTotalFace, const int *left_cell_of_face,
                                           const int *right_cell_of_face, const int *node_number_of_each_face,
                                           const int *face2node, const long long int *nodepos, double *q_n,
                                           int *n_count, RFloat *q);

    __global__ void GPUNodeVarAve(const int nTotalNode, double *q_n, int *n_count, const double SMALL);
    __global__ void GPUNodeBoundaryCalc(const int nBoundFace, const int *leftCellOfFace, const int *rightCellOfFace,
                                        const RFloat *q, const int *d_boundaryType, const int SOLID_SURFACE,
                                        const int SYMMETRY, const int *face2node, RFloat *qNode,
                                        const int *nodeNumberOfEachFace, const long long int *nodepos,
                                        const int nTotalCell, int index);
} //! namespace GPUKernels
