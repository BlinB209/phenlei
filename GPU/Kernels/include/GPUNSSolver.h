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
//! @file      GPUNSSolver.h
//! @brief     NSSolver functions.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include <cuda.h>
#include <cuda_runtime_api.h>
#include "Geo_Grid.h"
#include "NSSolver.h"
//! #include "FYUtil.h"
#include "FieldProxy.h"
#include "GeomProxy.h"
#include "FieldProxy.h"
#include "GeneralFieldProxy.h"
#include "Param_NSSolverUnstruct.h"

using namespace PHSPACE;

namespace GPUNSSolverUnstruct
{
    void            CallGPUBoundary(Grid *gridIn, Param_NSSolverUnstruct *parameters);
    __global__ void GPUBoundary(RFloat *q, RFloat *prim_inf, const int *d_boundaryType, const int *left_cell_of_face,
                                const int *right_cell_of_face, RDouble *xfn, RDouble *yfn, RDouble *zfn, RDouble *xtn,
                                RDouble *ytn, RDouble *ztn, RDouble *vgn, RFloat *gama, RDouble refMachNumber,
                                double refGama, double tw, const int nl, const int nchem, const int nBoundFace,
                                const int nTotal, const int iviscous, const RFloat twall, const int nm,
                                const int EXTRAPOLATION, const int SYMMETRY, const int FARFIELD, const int OUTFLOW,
                                const int ifLowSpeedPrecon, const int SOLID_SURFACE, const int INFLOW);

    __device__ void Outflow_BC(RFloat *prims, RFloat *primt, int nl, int nchem);
    __device__ void Inflow_BC(RFloat *prims, RFloat *primt, int nl, int nchem);
    __device__ void Symmetry_BC(RFloat *prims, RFloat *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn,
                                int nl, int nchem);
    //!!viscous adiabatic wall
    __device__ void Vis_Adi_Wall_BC(RFloat *prims, RFloat *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble xtn,
                                    RDouble ytn, RDouble ztn, RFloat reference_mach_number, int nl, int nchem);
    //!3D
    __device__ void Farfield_BC(RFloat *prims, RFloat *prim_inf, RFloat *primt, RFloat gama0, RFloat gama, RDouble nxs,
                                RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem);
    __device__ void Farfield_BC_check(RFloat *prims, RFloat *prim_inf, RFloat *primt, RFloat gama0, RFloat gama,
                                      RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem, int iFace);
    __device__ void Farfield_BC_checkMG(RFloat *prims, RFloat *prim_inf, RFloat *primt, RFloat gama0, RFloat gama,
                                        RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem,
                                        int iFace);
    //!void CallGPULoadQ(const int nTotal, const int nl);
    void CallGPULoadQ(Grid *gridIn, Param_NSSolverUnstruct *parameters);

    __global__ void GPULoadQ(RFloat *q, RFloat *qold, const int nTotal, const int nl);
    //!void CallGPUFillField(const int nTotal, const int neqn);
    void CallGPUFillField(Grid *gridIn, const int neqn);

    __global__ void GPUFillField(RFloat *field1, RFloat *field2, const int nTotal, const int neqn);
    //!void CallGPULoadResiduals(const int neqn, const int nTotal);
    void CallGPULoadResiduals(Grid *gridIn, const int neqn);

    __global__ void GPULoadResiduals(RFloat *res, RFloat *dq, RFloat *rhs, const int neqn, const int nTotal);

    void CallGPURecoverResidual(Grid *gridIn, const int neqn);

    __global__ void GPURecoverResidualNS(RFloat *res, RFloat *rhs, const int neqn, const int nTotal);
    //!void CallGPULhs(double coef, const int neqn, const int nTotalCell);
    void CallGPURungeKuttaResidual(Grid *gridIn, const int neqn, double coef);

    __global__ void GPULhs(RFloat *dq, RFloat *res, RFloat *dt, double coef, const int neqn, const int nTotal,
                           const int nTotalCell);

    void CallGPUUpdateFlowField(Grid *gridIn, Param_NSSolverUnstruct *parameters, const int neqn);
    void CallGPUUpdateFlowFieldLUSGSM1(RFloat density_limit, RFloat pressure_limit, RFloat mostNegativePressure,
                                       const int nTotalCell, const int neqn);
    __global__ void GPUUpdateFlowFieldM1_S1(const int nTotalCell, const int nTotal, const int neqn, RFloat *q,
                                            RFloat *dq, RFloat *qnew, RFloat *gamma, RFloat *t, RFloat density_limit,
                                            RFloat pressure_limit, RFloat mostNegativePressure,
                                            int *d_face_number_of_each_cell, int *d_cell2face,
                                            int *d_acc_face_number_of_cell, int *d_left_cell_of_face,
                                            int *d_right_cell_of_face);
    __global__ void GPUUpdateFlowFieldM1(const int nTotalCell, const int nTotal, const int neqn, RFloat *q, RFloat *dq,
                                         RFloat *qnew, RFloat *gamma, RFloat *t, RFloat density_limit,
                                         RFloat pressure_limit, RFloat mostNegativePressure,
                                         int *d_face_number_of_each_cell, int *d_cell2face,
                                         int *d_acc_face_number_of_cell, int *d_left_cell_of_face,
                                         int *d_right_cell_of_face);
    __device__ void Primitive2Conservative(RFloat *prim, RFloat gama, RFloat *q);
    __device__ void ComputeInternalEnergy(RFloat *prim, RFloat gama_in, RFloat &em);
    __device__ void Conservative2Primitive(RFloat *q, RFloat gama, RFloat *prim, RFloat &temperature);
    __device__ void GPUSolutionFix(RFloat *prim, int *face_number_of_each_cell, int *cell2face, int *left_cell_of_face,
                                   int *right_cell_of_face, int neqn, RFloat *q, int *d_acc_face_number_of_cell,
                                   int icell, int nTotal);
    __device__ void GPUSolutionFix_Z(RFloat *prim, int *face_number_of_each_cell, int *cell2Face,
                                     int *left_cell_of_face, int *right_cell_of_face, int neqn, RFloat *q,
                                     int *posiCell2Face, int icell, int nTotal);
    //!void CallGPUZeroResiduals(const int nl, const int nTotal);
    void CallGPUZeroResiduals(Grid *gridIn, const int nl);
    void CallGPUDtIni(const int nTotalCell);
    void CallGPUDtCFL(RFloat cfl, const int nTotalCell);
    void CallGPUVisSpectrumRadiusIni(const int nTotalCell);
    void CallGPUDtvIni(RFloat cfl, const int nTotalCell);
    void CallGPUDtFinal(const int nTotalCell);

    __global__ void GPUZeroResiduals(const int nl, RFloat *res, const int nTotal);
    __global__ void GPUDtIni(RFloat *dt, const int nTotalCell);
    __global__ void GPUVisSpectrumRadiusIni(RFloat *d_visSpectrumRadius, const int nTotalCell);
    __global__ void GPUDtvIni(RFloat cfl, RFloat *dtv, RFloat *vol, RFloat *visSpectrumRadius, const int nTotalCell);
    __global__ void GPUDtFinal(RFloat *dt, RFloat *dtv, const int nTotalCell);

    void CallGPULocalTimeStep(Grid *gridIn, const Param_NSSolverUnstruct *parameters,
                              const Param_NSSolver *parametersns);
    void CallGPULimitCFL(Grid *gridIn, const Param_NSSolverUnstruct *parameters, const Param_NSSolver *parametersns);
    __global__ void GPUSpectralRadiusInvis(RFloat *d_spectralRadius, RFloat *d_invSpectrumRadius, const int nTotalCell);
    __global__ void GPUSpectralRadiusVis(RFloat *d_spectralRadius, RFloat *d_visSpectrumRadius, const int nTotalCell);
    __global__ void GPUCFLCellIni(RDouble *d_CFLCell, RDouble CFL, const int nTotalCell);
    __global__ void GPUComputeCFLCell(RDouble *d_CFLCell, RFloat *d_q_ns, RDouble CFL, RDouble CFLStart,
                                      RDouble pLowerLimit, RDouble pUpperLimit, const int nTotalCell, const int nTotal);
    __global__ void GPUCorCFLCell(RDouble *d_CFLCell, RFloat *d_cellSkewness, const int nTotalCell);
    __global__ void GPUDtCFL(RFloat *d_dt, RFloat *d_spectralRadius, RFloat *vol, RDouble *d_CFLCell,
                             const int nTotalCell);

    void CallGPUReduceMaxTimeStepLoop1(RDouble ktmax, RDouble globalMinDt, const int nTotalCell);
    void CallGPUReduceMaxTimeStepLoop2(const int nTotalCell);

    void CallGPUReduceMaxTimeStep(double ktmax, RFloat globalMinDt, const int nTotalCell);
    void CallGPUReduceMaxTimeStepOneProcess(double ktmax, const int nTotalCell);
    void CallGPUReduceMaxTimeStep(Grid *gridIn);

    __global__ void GPUReduceMaxTimeStep1OneProcess(RFloat *dt, double ktmax, RFloat *minDt, const int nTotalCell);
    __global__ void GPUReduceMaxTimeStep1AWARE(RFloat *dt, double ktmax, double *globalMinDt, const int nTotalCell);

    __global__ void GPUReduceMaxTimeStep1(RFloat *dt, const double ktmax, const double globalMinDt,
                                          const int nTotalCell);
    __global__ void GPUReduceMaxTimeStep2(RFloat *dt, RFloat *vol, const int nTotalCell);
    void            CallGPUSetGhostCell(Grid *gridIn);
    __global__ void GPUSetGhostCell(RFloat *f, int *d_boundaryType, int *left_cell_of_face, int *right_cell_of_face,
                                    const int nBoundFace);

    void CallGPUCopyMinMaxTime(const double minDt, const double maxDt);
    void CallGPUComputeMinTimeStep(const int nTotalCell);
    void CallGPUComputeViscousCoefficientWithoutChemical(Grid *gridIn, Param_NSSolverUnstruct *parameters);
    void CallGPUStoreRhsByResidual(Grid *gridIn, const int neqn);

    __global__ void GPUComputeMinMaxStep1(RFloat *dt, RFloat *bmin_tmp, RFloat *bmax_tmp, const int nTotalCell);
    __global__ void GPUComputeMinMaxStep2(RFloat *bmin_tmp, RFloat *bmax_tmp, RFloat *minDt, RFloat *maxDt);
    __global__ void GPUCompViscousCoef(RFloat *visl, RFloat *t, double tsuth, double visl_min, const int nTotal,
                                       const int ITT);
    __global__ void GPUStoreRhsByResidual(RFloat *rhs, RFloat *res, const int neqn, const int nTotal);
} //! namespace GPUNSSolverUnstruct
