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
//! @file      GPULUSGS.h
//! @brief     Functions for LUSGS method.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUTurbBackwardSweepOneColor(const int nTotalFace, const int nTotal, const int groupID,
                                          const int inGroupNum, const int inGroupPosi, const int nTurbulenceEquation);

    __global__ void GPUTurbBackwardSweepOneColor(const int nTotalFace, const int nTotal, const int groupID,
                                                 const int inGroupNum, const int inGroupPosi,
                                                 const int nTurbulenceEquation, const int *InteriorCellGroup,
                                                 const int *faceNumberOfEachCell, const int *cell2Face,
                                                 const int *posiCell2Face, const int *leftCellofFace,
                                                 const int *rightCellofFace, const RFloat *matTurbl,
                                                 const RFloat *matTurbr, const RFloat *specTurb, RFloat *dq);

    void CallGPUTurbForwardSweepOneColor(const int nTotalFace, const int nTotal, const int groupID,
                                         const int inGroupNum, const int inGroupPosi, const int nTurbulenceEquation);

    __global__ void GPUTurbForwardSweepOneColor(const int nTotalFace, const int nTotal, const int groupID,
                                                const int inGroupNum, const int inGroupPosi,
                                                const int nTurbulenceEquation, const int *InteriorCellGroup,
                                                const int *faceNumberOfEachCell, const int *cell2Face,
                                                const int *posiCell2Face, const int *leftCellofFace,
                                                const int *rightCellofFace, const RFloat *matTurbl,
                                                const RFloat *matTurbr, const RFloat *specTurb, RFloat *dq);
    void CallGPUNSForwardSweepOneColor(const int nTotal, const int groupID, const int inGroupNum, const int inGroupPosi,
                                       const int nl, const int nm, const int nchem, const int nEquation,
                                       const int ifLowSpeedPrecondition, const int vis_run, const RDouble refReNumber);
    void CallGPUSetGhostDQLUSGS(const int nBoundFace, const int nTotal, const int nm, const int nl, const int nchem,
                                const int ntmodel, const int nEquation, const int isViscous);
    void CallGPUNSBackwardSweepOneColor(const int nTotal, const int groupID, const int inGroupNum,
                                        const int inGroupPosi, const int nl, const int nm, const int nchem,
                                        const int nEquation, const int ifLowSpeedPrecondition, const int vis_run,
                                        const RDouble refReNumber);
    __global__ void GPUNSBackwardSweepOneColor(
        const int nTotal, const int groupID, const int inGroupNum, const int inGroupPosi, const int nl, const int nm,
        const int nchem, const int nEquation, const int ifLowSpeedPrecondition, const int vis_run,
        const RDouble refReNumber, const int *InteriorCellGroup, const int *iblank, const int *faceNumberOfEachCell,
        const int *cell2Face, const int *posiCell2Face, const int *leftCellofFace, const int *rightCellofFace,
        const RDouble *xfn, const RDouble *yfn, const RDouble *zfn, const RDouble *area, const RDouble *xcc,
        const RDouble *ycc, const RDouble *zcc, const RDouble *vgn, const RFloat *q, const RFloat *gamma,
        const RFloat *t, const RFloat *viscousLaminar, const RFloat *viscousTurbulent, const RFloat *spec, RFloat *dq);
    __global__ void GPUSetGhostDQLUSGS(const int nBoundFace, const int nTotal, const int nm, const int nl,
                                       const int nchem, const int ntmodel, const int nEquation, const int isViscous,
                                       const int *leftCellofFace, const int *rightCellofFace, const int *boundaryType,
                                       const RDouble *xfn, const RDouble *yfn, const RDouble *zfn, const RFloat *q,
                                       const RFloat *t, const RFloat *gama, RFloat *dq);
    __device__ void Primitive2Conservative(const int nchem, const int nm, const int nl, const int ntmodel, RFloat *prim,
                                           RFloat gama, RFloat Tv, RFloat Te, RFloat *q);
    __device__ void ComputeInternalEnergy(const int nchem, RFloat *prim, RFloat gama_in, RFloat Tv, RFloat Te,
                                          RFloat &em);

    __global__ void GPUNSForwardSweepOneColor(
        const int nTotal, const int groupID, const int inGroupNum, const int inGroupPosi, const int nl, const int nm,
        const int nchem, const int nEquation, const int ifLowSpeedPrecondition, const int vis_run,
        const RDouble refReNumber, const int *InteriorCellGroup, const int *iblank, const int *faceNumberOfEachCell,
        const int *cell2Face, const int *posiCell2Face, const int *leftCellofFace, const int *rightCellofFace,
        const RDouble *xfn, const RDouble *yfn, const RDouble *zfn, const RDouble *area, const RDouble *xcc,
        const RDouble *ycc, const RDouble *zcc, const RDouble *vgn, const RFloat *q, const RFloat *gamma,
        const RFloat *t, const RFloat *viscousLaminar, const RFloat *viscousTurbulent, const RFloat *spec, RFloat *dq);
    __device__ RFloat GetRad(RFloat *prim, RDouble *nxyz, RFloat gama);

    __device__ void MXDQ_STD(const int nchem, RDouble *prim, RDouble nxs, RDouble nys, RDouble nzs, RDouble area,
                             RDouble vgn, RDouble gama, RDouble *dq, RDouble *f, int nm, int nLaminar, RDouble radius,
                             int ipn);

    __device__ void ComputeTotalEnthalpyAndDH(const int nchem, RFloat *prim, RFloat &gama, RFloat *dq,
                                              RFloat &totalEnthalpy, RFloat &dh);

}    //! namespace GPUKernels
