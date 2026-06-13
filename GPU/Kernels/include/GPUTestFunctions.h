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
//! @file      GPUTestFunctions.h
//! @brief     Test functions for GPU results.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include "FaceProxy.h"
#include "GeomProxy.h"
#include "Limiter.h"
#include "Precision.h"
#include "cudaErrorHandle.h"

using namespace PHSPACE;
namespace GPUTestSpace
{
    //! Test mesh nodes data allocation, d_x, d_y, d_z is mesh nodes on device
    void TestGPUNodesDataAllocCopy(RDouble *org_x, RDouble *org_y, RDouble *org_z, RDouble *d_x, RDouble *d_y,
                                   RDouble *d_z, int nTotalNode);
    //! Test cell center data allocation, d_xcc, d_ycc, d_zcc is cell center on
    //! device
    void TestGPUCellCentAllocCopy(RDouble *org_xcc, RDouble *org_ycc, RDouble *org_zcc, RDouble *d_xcc, RDouble *d_ycc,
                                  RDouble *d_zcc, int nTotalCell);
    //! Test face cell data allocation, d_left_cell_of_face, d_right_cell_of_face is
    //! face cell relationship on device
    void TestGPUFaceCellRelAllocCopy(const int *org_left_cell_of_face, const int *org_right_cell_of_face,
                                     int *d_left_cell_of_face, int *d_right_cell_of_face, int nTotalFace);
    //! Test face node data allocation, d_face2node, d_node_number_of_each_face is
    //! face node relationship on device. noting that the total number of elements in
    //! d_face2node should be calculated in advance.
    void TestGPUFaceNodeRelAllocCopy(const int *org_face2node, const int *org_node_number_of_each_face,
                                     int *d_face2node, int *d_node_number_of_each_face, int nTotalFace2Node,
                                     int nTotalFace);
    //! Test q_proxy copy from host to device
    void TestGPUQFullCopy(RFloat **q2d, RFloat *d_q_proxy, int nTotal, int nTVar);
    //! Test face normal data allocation on gpu
    void TestGPUFaceNormAllocCopy(RDouble *xfn, RDouble *yfn, RDouble *zfn, RDouble *d_xfn, RDouble *d_yfn,
                                  RDouble *d_zfn, const int n);
    //! Test face center data allocation on gpu
    void TestGPUFaceCentAllocCopy(RDouble *xfc, RDouble *yfc, RDouble *zfc, RDouble *d_xfc, RDouble *d_yfc,
                                  RDouble *d_zfc, const int n);
    void TestGPUAreaVolmAllocCopy(RDouble *area, RDouble *vol, RDouble *d_area, RDouble *d_vol, const int nTotalFace,
                                  const int nTotalCell);
    void TestGPUPointer(RDouble *d_ptr);
    void TestGPUQCopy(RFloat **q, RFloat *d_q_ns, const int eqn, const int nTotal);

    //! test the copy of gama
    void TestGPUGamaCopy(RFloat *gama, RFloat *gamal, RFloat *gamar, RFloat *d_gama_ns, RFloat *d_gamaL_ns,
                         RFloat *d_gamaR_ns, const int d_nTotal, const int len_gamaLR);
    //! test the copy of limiter
    void TestGPULIMITCopy(RFloat **LIMIT, RFloat *d_LIMIT, const int neqn_ns, const int nTotal);
    void TestGPULimitCopy(RFloat *limit, RFloat *d_limit, const int nTotal);
    //! test the copy of dqdx ...
    void TestGPUDqdxDqDyDqDzCopy(RDouble **dqdx, RDouble **dqdy, RDouble **dqdz, RFloat *d_dqdx_proxy,
                                 RFloat *d_dqdy_proxy, RFloat *d_dqdz_proxy, const int neqn_ns, const int nTotal);
    //! test the copy of xtn ytn ztn
    void TestGPUXtnYtnZtnCopy(RDouble *xtn, RDouble *ytn, RDouble *ztn, RDouble *d_xtn, RDouble *d_ytn, RDouble *d_ztn,
                              const int nTotalFace);
    //! test the copy of vgn
    void TestGPUVgnCopy(RDouble *vgn, RDouble *d_vgn, const int nTotalFace);
    //! test the copy of res_ns
    void TestGPUResCopy(RFloat **res, RFloat *d_res_ns, const int neqn, const int nTotal);
    //! test the copy of epsCell
    void TestGPUepsCellCopy(RFloat *epsCell, RFloat *d_epsCell, const int nTotalCell);
    void TestGPUDminDmaxAllocCopy(RFloat *dmax, RFloat *dmin, RFloat *d_dmax, RFloat *d_dmin, const int nTotal);
    void TestFaceCellRel(int nTotalCell, int *face_number_of_each_cell, int **cell2face, int *acc_face_number_of_cell);
} //! namespace GPUTestSpace
