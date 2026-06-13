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
//! @file      GPUKernelTest.h
//! @brief     Kernel test function for GPU results.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#ifndef GPUKERNELTEST
#define GPUKERNELTEST
#include <stdio.h>
#include <iostream>
#include "Geo_Grid.h"
#include "NSSolver.h"
#include "Precision.h"
#include "FaceProxy.h"
#include "GeomProxy.h"
#include "Flux_Inviscid.h"
const double abCompare = 1.0e-8;
const double rtCompare  = 0.0;
const double smallLimit = 1.0e-40;

using namespace PHSPACE;
void TestGPUGetQlQr(FaceProxy *face_proxy, const int nst, const int ned, const int SEG_LEN);
void TestGPUGetGamaLR(FaceProxy *face_proxy, const int nst, const int ned, const int SEG_LEN);
void TestGPUCompNodeVarByGradientLoop1(const string d_gradient_field_proxy, const int nTotalNode, RFloat *org_q_n,
                                       const int *org_n_count);
// test the GPUReConstructFaceValueLoop2
void TestGPUReConstructFaceValueLoop(RFloat **ql, RFloat **qr, RFloat *d_ql_ns, RFloat *d_qr_ns, const int d_neqn,
                                     const int d_SEG_LEN, const int nst, const int ned);
void TestGPUReConstructFaceValue(FaceProxy *face_proxy, const int nst, const int ned, const int SEG_LEN);
// test the GPUBoundaryQlQrFix
void TestGPUBoundaryQlQrFix(FaceProxy *face_proxy, const int nst, const int ned, const int SEG_LEN);
void TestGPUInv_scheme(InviscidSchemeParameter *invSchemePara, RFloat *d_flux, const int nst, const int ned,
                       const int SEG_LEN);
//! test the GPUInviscidFluxWrapLoop1
void TestGPUInviscidFluxWrapLoop1(RFloat **flux, RFloat *d_flux, const int nl, const int nst, const int ned,
                                  const int SEG_LEN);
//! test the GPULoadFlux
void TestGPULoadFlux(Grid *grid_in, const int nst, const int ned);
void TestGPUGetQlQrTurb(FaceProxy *face_proxy, const int neqn, const int nst, const int ned, const int SEG_LEN);
void TestGPUFlux(RFloat **flux, RFloat *d_flux, const int nl, const int nst, const int ned, const int SEG_LEN);
//! test the GPULoadFluxTurb
//! test the GPULoadFluxTurb
void TestGPULoadFluxTurb(RFloat **res_turb, RFloat *d_res_turb, const int *left_cell_of_face, const int nTotal,
                         const int neqn, const int nst, const int nMid);
void TestGPUTurbComputeVisflux(RFloat **res_turb, RFloat *d_res_turb, const int *left_cell_of_face, const int nTotal,
    const int neqn, const int nst, const int nMid);
// test CallGPUCompGamaAndTField
void TestGPUCompGamaAndTField(Grid *grid_in);
void TestGPUBoundary(RFloat **q);

double Max(double a, double b);

void TestGPULoadQ(RFloat **q);
void TestGPUFillField(RFloat **q);
void TestGPULoadResiduals(RFloat **res);
void TestGPULhs(RFloat **dq);
void TestGPUUpdateFlowFieldM1(RFloat **q);
void TestGPUZeroResiduals(RFloat **res, const int nl, const int nTotal);
void TestGPUDtCFL(RFloat *dt, const int nTotalCell);
void TestGPUReduceMaxTimeStep(RFloat *dt, const int nTotalCell);
void TestGPUSetGhostCell(RFloat *dt, const int nTotalCell);
void TestGPUComputeMinTimeStep(RFloat minDt, RFloat maxDt);
void TestGPUCompViscousCoef(RFloat *visl, const int nTotal);
void TestGPUGPUTurbViscosity(double vistmax, double vistmin);
void TestGPUTurbViscosity(const RFloat *org_vist, const int nTotal);
void TestGPUStoreRhsByResidual(RFloat **rhs, const int nl, const int nTotal);
#endif
