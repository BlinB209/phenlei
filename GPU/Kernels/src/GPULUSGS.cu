#include "GPULUSGS.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
#include "GPUBasicFunctionsPart2.h" //!for GPUSQR
#include "Constants.h"
#include "Geo_SimpleBC.h"

using namespace PHENGLEI; //!for constant in Geo_SimpleBC.h
namespace GPUKernels
{
    void CallGPUTurbBackwardSweepOneColor(const int nTotalFace, const int nTotal, const int groupID,
                                          const int inGroupNum, const int inGroupPosi, const int nTurbulenceEquation)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = inGroupNum;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUTurbBackwardSweepOneColor, loopLen, 0, gridSize, blockSize, blockSizeLimit);

#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUTurbBackwardSweepOneColor, 0, gridSize, blockSize);
#endif

        GPUTurbBackwardSweepOneColor<<<gridSize, blockSize>>>(
            nTotalFace, nTotal, groupID, inGroupNum, inGroupPosi, nTurbulenceEquation, d_InteriorCellGroup,
            d_face_number_of_each_cell, d_cell2Face, d_posiCell2Face, d_left_cell_of_face, d_right_cell_of_face,
            d_mat_turbl, d_mat_turbr, d_spec_turb, d_dq_turb);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUTurbBackwardSweepOneColor(const int nTotalFace, const int nTotal, const int groupID,
                                                 const int inGroupNum, const int inGroupPosi,
                                                 const int nTurbulenceEquation, const int *InteriorCellGroup,
                                                 const int *faceNumberOfEachCell, const int *cell2face,
                                                 const int *posiCell2Face, const int *leftCellOfFace,
                                                 const int *rightCellOfFace, const RFloat *matTurbl,
                                                 const RFloat *matTurbr, const RFloat *specTurb, RFloat *dq)
    {
        int    threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int    offset   = 0;
        int    face, le, re;
        double df[1];
        double rhs0[1];
        for (offset = threadID; offset < inGroupNum; offset += gridDim.x * blockDim.x)
        {
            int iCell = InteriorCellGroup[inGroupPosi + offset];
            for (int m = 0; m < nTurbulenceEquation; ++m)
            {
                rhs0[m] = 0.0;
            }

            int posiCell = posiCell2Face[iCell];

            for (int j = 0; j < faceNumberOfEachCell[iCell]; ++j)
            {
                face = cell2face[posiCell + j];
                le   = leftCellOfFace[face];
                re   = rightCellOfFace[face];
                //! One of le and re must be cell itself.
                if (le < iCell || re < iCell) continue;
                //! Now its neighboring cell belongs to lower triangular

                if (re == iCell)
                {
                    //!SWAP(le, re);
                    int tmp = le;
                    le      = re;
                    re      = tmp;
                    for (int m = 0; m < nTurbulenceEquation; ++m)
                    {
                        df[m] = matTurbl[m * nTotalFace + face] * dq[m * nTotal + re];
                    }
                }
                else
                {
                    for (int m = 0; m < nTurbulenceEquation; ++m)
                    {
                        df[m] = matTurbr[m * nTotalFace + face] * dq[m * nTotal + re];
                    }
                }

                //! Add Flux together
                for (int m = 0; m < nTurbulenceEquation; ++m)
                {
                    rhs0[m] += df[m];
                }
            }
            for (int m = 0; m < nTurbulenceEquation; ++m)
            {
                dq[m * nTotal + iCell] = dq[m * nTotal + iCell] - rhs0[m] / specTurb[m * nTotal + iCell];
            }
        }
    }

    void CallGPUTurbForwardSweepOneColor(const int nTotalFace, const int nTotal, const int groupID,
                                         const int inGroupNum, const int inGroupPosi, const int nTurbulenceEquation)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = inGroupNum;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUTurbForwardSweepOneColor, loopLen, 0, gridSize, blockSize, blockSizeLimit);

#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUTurbForwardSweepOneColor, 0, gridSize, blockSize);
#endif

        GPUTurbForwardSweepOneColor<<<gridSize, blockSize>>>(
            nTotalFace, nTotal, groupID, inGroupNum, inGroupPosi, nTurbulenceEquation, d_InteriorCellGroup,
            d_face_number_of_each_cell, d_cell2Face, d_posiCell2Face, d_left_cell_of_face, d_right_cell_of_face,
            d_mat_turbl, d_mat_turbr, d_spec_turb, d_dq_turb);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUTurbForwardSweepOneColor(const int nTotalFace, const int nTotal, const int groupID,
                                                const int inGroupNum, const int inGroupPosi,
                                                const int nTurbulenceEquation, const int *InteriorCellGroup,
                                                const int *faceNumberOfEachCell, const int *cell2face,
                                                const int *posiCell2Face, const int *leftCellOfFace,
                                                const int *rightCellOfFace, const RFloat *matTurbl,
                                                const RFloat *matTurbr, const RFloat *specTurb, RFloat *dq)
    {
        int    threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int    offset   = 0;
        int    face, le, re;
        double df[1];
        double rhs0[1];
        for (offset = threadID; offset < inGroupNum; offset += gridDim.x * blockDim.x)
        {
            int iCell = InteriorCellGroup[inGroupPosi + offset];
            for (int m = 0; m < nTurbulenceEquation; ++m)
            {
                rhs0[m] = 0.0;
            }

            int posiCell = posiCell2Face[iCell];

            for (int j = 0; j < faceNumberOfEachCell[iCell]; ++j)
            {
                face = cell2face[posiCell + j];
                le   = leftCellOfFace[face];
                re   = rightCellOfFace[face];
                //! One of le and re must be cell itself.
                if (le > iCell || re > iCell) continue;
                //! Now its neighboring cell belongs to lower triangular

                if (re == iCell)
                {
                    //!SWAP(le, re);
                    int tmp = le;
                    le      = re;
                    re      = tmp;
                    for (int m = 0; m < nTurbulenceEquation; ++m)
                    {
                        df[m] = matTurbl[m * nTotalFace + face] * dq[m * nTotal + re];
                    }
                }
                else
                {
                    for (int m = 0; m < nTurbulenceEquation; ++m)
                    {
                        df[m] = matTurbr[m * nTotalFace + face] * dq[m * nTotal + re];
                    }
                }

                //! Add Flux together
                for (int m = 0; m < nTurbulenceEquation; ++m)
                {
                    rhs0[m] += df[m];
                }
            }
            for (int m = 0; m < nTurbulenceEquation; ++m)
            {
                dq[m * nTotal + iCell] = (dq[m * nTotal + iCell] - rhs0[m]) / specTurb[m * nTotal + iCell];
            }
        }
    }

    void CallGPUNSForwardSweepOneColor(const int nTotal, const int groupID, const int inGroupNum, const int inGroupPosi,
                                       const int nl, const int nm, const int nchem, const int nEquation,
                                       const int ifLowSpeedPrecondition, const int vis_run, const RDouble refReNumber)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = inGroupNum;
        int blockSizeLimit = 0;
        //!int regsPerThread= 224;
        //!int regsPerThread= 0;
        //!int residentBlockNo= 1;
        //!int residentBlockNo= 0;
        //!size_t dsMemPerThread= 8 * nEquation * sizeof(RFloat);  //! calculate by the intended shared memory use in the kernel
        //!size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUNSForwardSweepOneColor, loopLen, 0, gridSize, blockSize, blockSizeLimit);
//!KernelLaunchPara((void*)GPUNSForwardSweepOneColor, loopLen, regsPerThread, residentBlockNo, dsMemPerThread, dsMemPerBlock, gridSize, blockSize, GPUProp);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUNSForwardSweepOneColor, 0, gridSize, blockSize);
#endif
        //!GPUNSForwardSweepOneColor<<<gridSize, blockSize, dsMemPerBlock>>>(nTotal, groupID, inGroupNum, inGroupPosi, nl, nm, nchem, nEquation, ifLowSpeedPrecondition, vis_run, refReNumber, d_InteriorCellGroup, d_blank, d_face_number_of_each_cell, d_cell2Face, d_posiCell2Face, d_left_cell_of_face, d_right_cell_of_face, d_xfn, d_yfn, d_zfn, d_area, d_xcc, d_ycc, d_zcc, d_vgn, d_q_ns, d_gama_ns, d_t_proxy, d_visl, d_vist, d_spec_ns, d_dq_ns);
        //!printf("gridSize = %d, blockSize = %d\n", gridSize, blockSize);
        GPUNSForwardSweepOneColor<<<gridSize, blockSize>>>(
            nTotal, groupID, inGroupNum, inGroupPosi, nl, nm, nchem, nEquation, ifLowSpeedPrecondition, vis_run,
            refReNumber, d_InteriorCellGroup, d_blank, d_face_number_of_each_cell, d_cell2Face, d_posiCell2Face,
            d_left_cell_of_face, d_right_cell_of_face, d_xfn, d_yfn, d_zfn, d_area, d_xcc, d_ycc, d_zcc, d_vgn, d_q_ns,
            d_gama_ns, d_t_proxy, d_visl, d_vist, d_spec_ns, d_dq_ns);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        //!GPUNSForwardSweepOneColor<<<1, 1, dsMemPerBlock>>>(nTotal, groupID, inGroupNum, inGroupPosi, nl, nm, nchem, nEquation, ifLowSpeedPrecondition, vis_run, refReNumber, d_InteriorCellGroup, d_blank, d_face_number_of_each_cell, d_cell2Face, d_posiCell2Face, d_left_cell_of_face, d_right_cell_of_face, d_xfn, d_yfn, d_zfn, d_area, d_xcc, d_ycc, d_zcc, d_vgn, d_q_ns, d_gama_ns, d_t_proxy, d_visl, d_vist, d_spec_ns, d_dq_ns);
    }

    void CallGPUSetGhostDQLUSGS(const int nBoundFace, const int nTotal, const int nm, const int nl, const int nchem,
                                const int ntmodel, const int nEquation, const int isViscous)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nBoundFace;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUSetGhostDQLUSGS, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSetGhostDQLUSGS, 0, gridSize, blockSize);
#endif

        GPUSetGhostDQLUSGS<<<gridSize, blockSize>>>(nBoundFace, nTotal, nm, nl, nchem, ntmodel, nEquation, isViscous,
                                                    d_left_cell_of_face, d_right_cell_of_face, d_boundaryType, d_xfn,
                                                    d_yfn, d_zfn, d_q_ns, d_t_proxy, d_gama_ns, d_dq_ns);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUNSBackwardSweepOneColor(const int nTotal, const int groupID, const int inGroupNum,
                                        const int inGroupPosi, const int nl, const int nm, const int nchem,
                                        const int nEquation, const int ifLowSpeedPrecondition, const int vis_run,
                                        const RDouble refReNumber)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = inGroupNum;
        int blockSizeLimit = 0;
        //!int regsPerThread= 224;
        //!int regsPerThread= 0;
        //!int residentBlockNo= 1;
        //!int residentBlockNo= 0;
        //!size_t dsMemPerThread= 8 * nEquation * sizeof(RFloat);  //! calculate by the intended shared memory use in the kernel
        //!size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUNSBackwardSweepOneColor, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUNSBackwardSweepOneColor, 0, gridSize, blockSize);
#endif
        //!GPUNSForwardSweepOneColor<<<gridSize, blockSize, dsMemPerBlock>>>(nTotal, groupID, inGroupNum, inGroupPosi, nl, nm, nchem, nEquation, ifLowSpeedPrecondition, vis_run, refReNumber, d_InteriorCellGroup, d_blank, d_face_number_of_each_cell, d_cell2Face, d_posiCell2Face, d_left_cell_of_face, d_right_cell_of_face, d_xfn, d_yfn, d_zfn, d_area, d_xcc, d_ycc, d_zcc, d_vgn, d_q_ns, d_gama_ns, d_t_proxy, d_visl, d_vist, d_spec_ns, d_dq_ns);
        //!printf("gridSize = %d, blockSize = %d\n", gridSize, blockSize);
        GPUNSBackwardSweepOneColor<<<gridSize, blockSize>>>(
            nTotal, groupID, inGroupNum, inGroupPosi, nl, nm, nchem, nEquation, ifLowSpeedPrecondition, vis_run,
            refReNumber, d_InteriorCellGroup, d_blank, d_face_number_of_each_cell, d_cell2Face, d_posiCell2Face,
            d_left_cell_of_face, d_right_cell_of_face, d_xfn, d_yfn, d_zfn, d_area, d_xcc, d_ycc, d_zcc, d_vgn, d_q_ns,
            d_gama_ns, d_t_proxy, d_visl, d_vist, d_spec_ns, d_dq_ns);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUNSBackwardSweepOneColor(
        const int nTotal, const int groupID, const int inGroupNum, const int inGroupPosi, const int nl, const int nm,
        const int nchem, const int nEquation, const int ifLowSpeedPrecondition, const int vis_run,
        const RDouble refReNumber, const int *InteriorCellGroup, const int *iblank, const int *faceNumberOfEachCell,
        const int *cell2face, const int *posiCell2Face, const int *leftCellofFace, const int *rightCellofFace,
        const RDouble *xfn, const RDouble *yfn, const RDouble *zfn, const RDouble *area, const RDouble *xcc,
        const RDouble *ycc, const RDouble *zcc, const RDouble *vgn, const RFloat *q, const RFloat *gamma,
        const RFloat *t, const RFloat *viscousLaminar, const RFloat *viscousTurbulent, const RFloat *spec, RFloat *dq)
    {
        using namespace IDX;
        //!extern __shared__ RFloat array[];
        /*
            RFloat * q_loc  = array + blockDim.x * nEquation * 0;
            RFloat * dq_loc = array + blockDim.x * nEquation * 1;
        RFloat * dqtry  = array + blockDim.x * nEquation * 2;
        RFloat * rhs0   = array + blockDim.x * nEquation * 3;
        RFloat * df     = array + blockDim.x * nEquation * 4;
        RFloat * q_pre  = array + blockDim.x * nEquation * 5;
        RFloat *primitiveVariable = array + blockDim.x * nEquation * 6;
        RFloat *facePrimitiveVariable = array + blockDim.x * nEquation * 7;
        */
        RFloat  q_loc[5];
        RFloat  dq_loc[5];
        RFloat  dqtry[5];
        RFloat  rhs0[5];
        RFloat  df[5];
        RFloat  q_pre[5];
        RFloat  primitiveVariable[5];
        RFloat  facePrimitiveVariable[5];
        int     cellGroup = blockDim.x * blockIdx.x + threadIdx.x;
        int     le, re, face, posiCellFace;
        RDouble nxyz[5];

        for (int offset = cellGroup; offset < inGroupNum; offset += gridDim.x * blockDim.x)
        {
            int iCell = InteriorCellGroup[offset + inGroupPosi];
            //!Initialize temp variables
            for (int eqnID = 0; eqnID < nEquation; eqnID++)
            {
                df[eqnID] = 0.0;
            }
            //!Operations for one cell
            if (iblank[iCell] <= 0)
            {
                for (int m = 0; m < nl; ++m)
                {
                    //!dq[m][iCell] = 0.0;
                    dq[m * nTotal + iCell] = 0.0;
                }
                continue;
            }

            for (int m = 0; m < nl; ++m)
            {
                rhs0[m] = 0.0;
            }

            RFloat tm = t[ITT * nTotal + iCell];

            posiCellFace = posiCell2Face[iCell];
            for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++jFace)
            {
                face = cell2face[posiCellFace + jFace];
                le   = leftCellofFace[face];
                re   = rightCellofFace[face];
                //! One of le and re must be cell itself.
                if (le < iCell || re < iCell) continue;
                //! Now its neighboring cell belongs to lower triangular

                nxyz[0] = xfn[face];
                nxyz[1] = yfn[face];
                nxyz[2] = zfn[face];
                nxyz[3] = area[face];
                nxyz[4] = vgn[face];

                if (re == iCell)
                {
                    //!SWAP(le, re);
                    int tmp = re;
                    re      = le;
                    le      = tmp;
                    nxyz[0] = -nxyz[0];
                    nxyz[1] = -nxyz[1];
                    nxyz[2] = -nxyz[2];
                    nxyz[4] = -nxyz[4];
                }

                for (int m = 0; m < nEquation; ++m)
                {
                    primitiveVariable[m] = q[m * nTotal + re];
                }

                for (int m = 0; m < nEquation; ++m)
                {
                    q_loc[m] = q[m * nTotal + re];
                }

                for (int m = 0; m < nl; ++m)
                {
                    dq_loc[m] = dq[m * nTotal + re];
                }

                //!* add
                {
                    RFloat rl = q[IDX::IR * nTotal + le];
                    RFloat ul = q[IDX::IU * nTotal + le];
                    RFloat vl = q[IDX::IV * nTotal + le];
                    RFloat wl = q[IDX::IW * nTotal + le];
                    RFloat pl = q[IDX::IP * nTotal + le];

                    RFloat rr = q[IDX::IR * nTotal + re];
                    RFloat ur = q[IDX::IU * nTotal + re];
                    RFloat vr = q[IDX::IV * nTotal + re];
                    RFloat wr = q[IDX::IW * nTotal + re];
                    RFloat pr = q[IDX::IP * nTotal + re];

                    RFloat specificHeatRatioLeft  = gamma[le];
                    RFloat specificHeatRatioRight = gamma[re];

                    RFloat hl = specificHeatRatioLeft / (specificHeatRatioLeft - 1.0) * pl / rl
                                + half * (ul * ul + vl * vl + wl * wl);
                    RFloat hr = specificHeatRatioRight / (specificHeatRatioRight - 1.0) * pr / rr
                                + half * (ur * ur + vr * vr + wr * wr);

                    RFloat tmp0 = sqrt(rr / rl);
                    RFloat tmp1 = 1.0 / (1.0 + tmp0);

                    RFloat rm = sqrt(rr * rl);
                    RFloat um = (ul + ur * tmp0) * tmp1;
                    RFloat vm = (vl + vr * tmp0) * tmp1;
                    RFloat wm = (wl + wr * tmp0) * tmp1;
                    RFloat hm = (hl + hr * tmp0) * tmp1;
                    RFloat pm =
                        rm * (hm - half * GPUSQR(um, vm, wm)) * (specificHeatRatioLeft - 1.0) / specificHeatRatioLeft;

                    facePrimitiveVariable[IDX::IR] = rm;
                    facePrimitiveVariable[IDX::IU] = um;
                    facePrimitiveVariable[IDX::IV] = vm;
                    facePrimitiveVariable[IDX::IW] = wm;
                    facePrimitiveVariable[IDX::IP] = pm;
                }
                //!*/

                if (ifLowSpeedPrecondition == 0)
                {
                    RFloat gama = gamma[re];

                    RFloat rad = GetRad(facePrimitiveVariable, nxyz, gama);
                    MXDQ_STD(nchem, primitiveVariable, nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], gama, dq_loc, df,
                             nm, nl, rad, -1);

                    //!MatrixFreeDF(q_loc, dq_loc, nxyz, refGama, df, nl);
                }
                else
                {
                    printf("Error: in GPUNSBackwardSweepOneColor, the other branch\n");
                    //!exit(1);
                    //!Compute the inviscous difference flux (df).
                    //!PreconMatrixFreeDF(q_loc, dq_loc, nxyz, MG, tm, refGama, df, nl);
                }

                if (vis_run)
                {
                    //!*
                    double distance;
                    distance = fabs(xfn[face] * (xcc[re] - xcc[le]) + yfn[face] * (ycc[re] - ycc[le])
                                    + zfn[face] * (zcc[re] - zcc[le]));

                    RFloat viscousLaminarFace   = half * (viscousLaminar[le] + viscousLaminar[re]);
                    RFloat viscousTurbulentFace = half * (viscousTurbulent[le] + viscousTurbulent[re]);
                    RFloat density              = half * (q[IDX::IR * nTotal + le] + q[IDX::IR * nTotal + re]);

                    RFloat viscosity             = viscousLaminarFace + viscousTurbulentFace;
                    RFloat faceArea              = area[face];
                    RFloat viscousSpectrumRadius = 2.0 * viscosity / (density * distance * refReNumber + SMALL);
                    viscousSpectrumRadius *= half * faceArea;

                    for (int m = 0; m < nl; ++m)
                    {
                        RFloat deltaVisFlux = viscousSpectrumRadius * dq_loc[m];
                        df[m] -= deltaVisFlux;
                    }
                }

                //! Add Flux together
                for (int m = 0; m < nl; ++m)
                {
                    rhs0[m] += df[m];
                }
            }

            if (ifLowSpeedPrecondition == 1)
            {
                printf("Error: in GPUNSBackwardSweepOneColor, ifLowSpeedPrecondition = 1\n");
                //!exit(1);
                //!ComputePreconMatrix(1, MG, q_pre, tm, refGama);
                //!PreconMatrixDQ(dq, MG, nl, iCell);
            }

            for (int m = 0; m < nl; ++m)
            {
                dqtry[m] = dq[m * nTotal + iCell] - rhs0[m] / spec[iCell];
            }

            //!RestrictDQ(dqtry, nl, dqLimit);
            for (int m = 0; m < nl; ++m)
            {
                dq[m * nTotal + iCell] = dqtry[m];
            }
        }
    }

    __global__ void GPUSetGhostDQLUSGS(const int nBoundFace, const int nTotal, const int nm, const int nl,
                                       const int nchem, const int ntmodel, const int nEquation, const int isViscous,
                                       const int *leftCellofFace, const int *rightCellofFace, const int *boundaryType,
                                       const RDouble *xfn, const RDouble *yfn, const RDouble *zfn, const RFloat *q,
                                       const RFloat *t, const RFloat *gama, RFloat *dq)
    {
        using namespace IDX;
        RFloat tv, te;
        int    bcType;
        int    le, re;
        int    threadID = blockDim.x * blockIdx.x + threadIdx.x;
        int    iFace    = 0;
        RFloat rm, um, vm, wm, pm, ub, vb, wb;
        RFloat rvn, ru, rv, rw, p_new;
        RFloat gam1;
        RFloat primitiveVariable[5];
        RFloat Q[5];
        for (iFace = threadID; iFace < nBoundFace; iFace += blockDim.x * gridDim.x)
        {
            le = leftCellofFace[iFace];
            re = rightCellofFace[iFace];

            //! Conservative Variables on level n+1

            for (int m = 0; m < nEquation; ++m)
            {
                primitiveVariable[m] = q[m * nTotal + le];
            }

            tv = t[ITT * nTotal + le];
            te = tv;
            if (ntmodel > 1) //!Consider Two-Temperature model and Three-Temperature model
            {
                printf("Error in GPUSetGhostDQLUSGS, ntmodel > 1\n");
                /*
                    tv = t[ITV][le];
                        te = tv;
                        if (ntmodel == 3){
                            te = t[ITE][le];
                        }
                */
            }

            Primitive2Conservative(nchem, nm, nl, ntmodel, primitiveVariable, gama[le], tv, te, Q);

            for (int m = 0; m < nl; ++m)
            {
                Q[m] += dq[m * nTotal + le];
            }
            for (int m = 0; m < ntmodel - 1; m++)
            {
                Q[nl + nchem + m] += dq[(nl + nchem + m) * nTotal + le];
            }

            //! non-conservative variables in re in level n
            rm = q[IR * nTotal + re];
            um = q[IU * nTotal + re];
            vm = q[IV * nTotal + re];
            wm = q[IW * nTotal + re];
            pm = q[IP * nTotal + re];

            bcType = boundaryType[iFace];

            for (int m = 0; m < nl; ++m)
            {
                dq[m * nTotal + re] = dq[m * nTotal + le];
            }
            for (int m = 0; m < ntmodel - 1; m++)
            {
                dq[(nl + nchem + m) * nTotal + re] = dq[(nl + nchem + m) * nTotal + le];
            }

            if (bcType == PHENGLEI::SYMMETRY || ((bcType == PHENGLEI::SOLID_SURFACE) && !isViscous))
            {
                //! Inviscid wall or symmetry boundary
                dq[IR * nTotal + re] = dq[IR * nTotal + le];
                rvn                  = 2.0
                      * (xfn[iFace] * dq[IRU * nTotal + le] + yfn[iFace] * dq[IRV * nTotal + le]
                         + zfn[iFace] * dq[IRW * nTotal + le]);
                dq[IRU * nTotal + re] = dq[IRU * nTotal + le] - xfn[iFace] * rvn;
                dq[IRV * nTotal + re] = dq[IRV * nTotal + le] - yfn[iFace] * rvn;
                dq[IRW * nTotal + re] = dq[IRW * nTotal + le] - zfn[iFace] * rvn;
                rvn                   = 2.0 * (xfn[iFace] * Q[IRU] + yfn[iFace] * Q[IRV] + zfn[iFace] * Q[IRW]);
                ru                    = Q[IRU] - xfn[iFace] * rvn;
                rv                    = Q[IRV] - yfn[iFace] * rvn;
                rw                    = Q[IRW] - zfn[iFace] * rvn;
                p_new                 = Q[IRE] - 0.5 * (Q[IRU] * Q[IRU] + Q[IRV] * Q[IRV] + Q[IRW] * Q[IRW]) / Q[IR];

                gam1 = gama[le] - 1.0;
                dq[IRE * nTotal + re] =
                    p_new - pm / gam1
                    + half * ((ru * ru + rv * rv + rw * rw) / Q[IR] - rm * (um * um + vm * vm + wm * wm));
            }
            else if ((bcType == PHENGLEI::SOLID_SURFACE) && isViscous)
            {
                ub                    = 0.0;
                vb                    = 0.0;
                wb                    = 0.0;
                dq[IR * nTotal + re]  = dq[IR * nTotal + le];
                dq[IRU * nTotal + re] = -dq[IRU * nTotal + le] + two * ub * dq[IR * nTotal + le];
                dq[IRV * nTotal + re] = -dq[IRV * nTotal + le] + two * vb * dq[IR * nTotal + le];
                dq[IRW * nTotal + re] = -dq[IRW * nTotal + le] + two * wb * dq[IR * nTotal + le];
                dq[IRE * nTotal + re] = dq[IRE * nTotal + le] + ub * (dq[IRU * nTotal + re] - dq[IRU * nTotal + le])
                                        + vb * (dq[IRV * nTotal + re] - dq[IRV * nTotal + le])
                                        + wb * (dq[IRW * nTotal + re] - dq[IRW * nTotal + le]);
            }
            else if (bcType == PHENGLEI::EXTRAPOLATION)
            {
                for (int m = 0; m < nl; ++m)
                {
                    dq[m * nTotal + re] = dq[m * nTotal + le];
                }
                for (int m = 0; m < ntmodel - 1; m++)
                {
                    dq[(nl + nchem + m) * nTotal + re] = dq[(nl + nchem + m) * nTotal + le];
                }
            }
            else
            {
                //! FIX_ALL set DQ = 0
                //! Far_FIELD, PROFILE and FIX_V have not been implemented yet
                for (int m = 0; m < nl; ++m)
                {
                    dq[m * nTotal + re] = 0.0;
                }

                for (int m = 0; m < ntmodel - 1; m++)
                {
                    dq[(nl + nchem + m) * nTotal + re] = 0.0;
                }
            }
        }
    }

    __global__ void GPUNSForwardSweepOneColor(
        const int nTotal, const int groupID, const int inGroupNum, const int inGroupPosi, const int nl, const int nm,
        const int nchem, const int nEquation, const int ifLowSpeedPrecondition, const int vis_run,
        const RDouble refReNumber, const int *InteriorCellGroup, const int *iblank, const int *faceNumberOfEachCell,
        const int *cell2face, const int *posiCell2Face, const int *leftCellofFace, const int *rightCellofFace,
        const RDouble *xfn, const RDouble *yfn, const RDouble *zfn, const RDouble *area, const RDouble *xcc,
        const RDouble *ycc, const RDouble *zcc, const RDouble *vgn, const RFloat *q, const RFloat *gamma,
        const RFloat *t, const RFloat *viscousLaminar, const RFloat *viscousTurbulent, const RFloat *spec, RFloat *dq)
    {
        using namespace IDX;
        //!extern __shared__ RFloat array[];
        /*
            RFloat * q_loc  = array + blockDim.x * nEquation * 0;
            RFloat * dq_loc = array + blockDim.x * nEquation * 1;
        RFloat * dqtry  = array + blockDim.x * nEquation * 2;
        RFloat * rhs0   = array + blockDim.x * nEquation * 3;
        RFloat * df     = array + blockDim.x * nEquation * 4;
        RFloat * q_pre  = array + blockDim.x * nEquation * 5;
        RFloat *primitiveVariable = array + blockDim.x * nEquation * 6;
        RFloat *facePrimitiveVariable = array + blockDim.x * nEquation * 7;
        */
        RFloat  q_loc[5];
        RFloat  dq_loc[5];
        RFloat  dqtry[5];
        RFloat  rhs0[5];
        RFloat  df[5];
        RFloat  q_pre[5];
        RFloat  primitiveVariable[5];
        RFloat  facePrimitiveVariable[5];
        int     cellGroup = blockDim.x * blockIdx.x + threadIdx.x;
        int     le, re, face, posiCellFace;
        RDouble nxyz[5];
        for (int offset = cellGroup; offset < inGroupNum; offset += gridDim.x * blockDim.x)
        {
            int iCell = InteriorCellGroup[offset + inGroupPosi];
            //!Initialize temp variables
            for (int eqnID = 0; eqnID < nEquation; eqnID++)
            {
                //!q_loc[eqnID] = 0.0;
                //!dq_loc[eqnID] = 0.0;
                //!rhs0[eqnID] = 0.0;
                df[eqnID] = 0.0;
                //!dqtry[eqnID] = 0.0;
                //!q_pre[eqnID] = 0.0;
                //!primitiveVariable[eqnID] = 0.0;
                //!facePrimitiveVariable[eqnID] = 0.0;
            }
            //!Operations for one cell
            if (iblank[iCell] <= 0)
            {
                for (int m = 0; m < nl; ++m)
                {
                    //!dq[m][iCell] = 0.0;
                    dq[m * nTotal + iCell] = 0.0;
                }
                continue;
            }

            for (int m = 0; m < nl; ++m)
            {
                rhs0[m]  = 0.0;
                q_pre[m] = q[m * nTotal + iCell];
            }

            RFloat tm = t[ITT * nTotal + iCell];

            posiCellFace = posiCell2Face[iCell];
            for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++jFace)
            {
                face = cell2face[posiCellFace + jFace];
                le   = leftCellofFace[face];
                re   = rightCellofFace[face];
                //! One of le and re must be cell itself.
                if (le > iCell || re > iCell) continue;
                //! Now its neighboring cell belongs to lower triangular
                nxyz[0] = xfn[face];
                nxyz[1] = yfn[face];
                nxyz[2] = zfn[face];
                nxyz[3] = area[face];
                nxyz[4] = vgn[face];

                if (re == iCell)
                {
                    //!SWAP(le, re);
                    int tmp = re;
                    re      = le;
                    le      = tmp;
                    nxyz[0] = -nxyz[0];
                    nxyz[1] = -nxyz[1];
                    nxyz[2] = -nxyz[2];
                    nxyz[4] = -nxyz[4];
                }

                for (int m = 0; m < nEquation; ++m)
                {
                    primitiveVariable[m] = q[m * nTotal + re];
                }

                for (int m = 0; m < nEquation; ++m)
                {
                    q_loc[m] = q[m * nTotal + re];
                }

                for (int m = 0; m < nl; ++m)
                {
                    dq_loc[m] = dq[m * nTotal + re];
                }

                //!* add
                {
                    RFloat rl = q[IDX::IR * nTotal + le];
                    RFloat ul = q[IDX::IU * nTotal + le];
                    RFloat vl = q[IDX::IV * nTotal + le];
                    RFloat wl = q[IDX::IW * nTotal + le];
                    RFloat pl = q[IDX::IP * nTotal + le];

                    RFloat rr = q[IDX::IR * nTotal + re];
                    RFloat ur = q[IDX::IU * nTotal + re];
                    RFloat vr = q[IDX::IV * nTotal + re];
                    RFloat wr = q[IDX::IW * nTotal + re];
                    RFloat pr = q[IDX::IP * nTotal + re];

                    RFloat specificHeatRatioLeft  = gamma[le];
                    RFloat specificHeatRatioRight = gamma[re];

                    RFloat hl = specificHeatRatioLeft / (specificHeatRatioLeft - 1.0) * pl / rl
                                + half * (ul * ul + vl * vl + wl * wl);
                    RFloat hr = specificHeatRatioRight / (specificHeatRatioRight - 1.0) * pr / rr
                                + half * (ur * ur + vr * vr + wr * wr);

                    RFloat tmp0 = sqrt(rr / rl);
                    RFloat tmp1 = 1.0 / (1.0 + tmp0);

                    RFloat rm = sqrt(rr * rl);
                    RFloat um = (ul + ur * tmp0) * tmp1;
                    RFloat vm = (vl + vr * tmp0) * tmp1;
                    RFloat wm = (wl + wr * tmp0) * tmp1;
                    RFloat hm = (hl + hr * tmp0) * tmp1;
                    RFloat pm =
                        rm * (hm - half * GPUSQR(um, vm, wm)) * (specificHeatRatioLeft - 1.0) / specificHeatRatioLeft;

                    facePrimitiveVariable[IDX::IR] = rm;
                    facePrimitiveVariable[IDX::IU] = um;
                    facePrimitiveVariable[IDX::IV] = vm;
                    facePrimitiveVariable[IDX::IW] = wm;
                    facePrimitiveVariable[IDX::IP] = pm;
                }
                //!*/
                if (ifLowSpeedPrecondition == 0)
                {
                    RFloat gama = gamma[re];

                    RFloat rad = GetRad(facePrimitiveVariable, nxyz, gama);
                    MXDQ_STD(nchem, primitiveVariable, nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], gama, dq_loc, df,
                             nm, nl, rad, -1);

                    //!MatrixFreeDF(q_loc, dq_loc, nxyz, refGama, df, nl);
                }
                else
                {
                    printf("Error: in GPUNSForwardSweepOneColor, the other branch\n");
                    //!exit(1);
                    //!Compute the inviscous difference flux (df).
                    //!PreconMatrixFreeDF(q_loc, dq_loc, nxyz, MG, tm, refGama, df, nl);
                }

                if (vis_run)
                {
                    double distance;
                    distance = fabs(xfn[face] * (xcc[re] - xcc[le]) + yfn[face] * (ycc[re] - ycc[le])
                                    + zfn[face] * (zcc[re] - zcc[le]));

                    RFloat viscousLaminarFace   = half * (viscousLaminar[le] + viscousLaminar[re]);
                    RFloat viscousTurbulentFace = half * (viscousTurbulent[le] + viscousTurbulent[re]);
                    RFloat density              = half * (q[IDX::IR * nTotal + le] + q[IDX::IR * nTotal + re]);

                    RFloat viscosity             = viscousLaminarFace + viscousTurbulentFace;
                    RFloat faceArea              = area[face];
                    RFloat viscousSpectrumRadius = 2.0 * viscosity / (density * distance * refReNumber + SMALL);
                    viscousSpectrumRadius *= half * faceArea;

                    for (int m = 0; m < nl; ++m)
                    {
                        RFloat deltaVisFlux = viscousSpectrumRadius * dq_loc[m];
                        df[m] -= deltaVisFlux;
                    }
                }

                //! Add Flux together
                for (int m = 0; m < nl; ++m)
                {
                    rhs0[m] += df[m];
                }
            }

            if (ifLowSpeedPrecondition == 1)
            {
                printf("Error: in GPUNSForwardSweepOneColor, ifLowSpeedPrecondition = 1\n");
                //!exit(1);
                //!ComputePreconMatrix(1, MG, q_pre, tm, refGama);
                //!PreconMatrixDQ(dq, MG, nl, iCell);
            }

            for (int m = 0; m < nl; ++m)
            {
                dqtry[m] = (dq[m * nTotal + iCell] - rhs0[m]) / spec[iCell];
            }

            //!RestrictDQ(dqtry, nl, dqLimit);
            for (int m = 0; m < nl; ++m)
            {
                dq[m * nTotal + iCell] = dqtry[m];
            }
        }
    }

    __device__ RFloat GetRad(RFloat *prim, RDouble *nxyz, RFloat gama)
    {
        using namespace IDX;

        RDouble &xfn  = nxyz[0];
        RDouble &yfn  = nxyz[1];
        RDouble &zfn  = nxyz[2];
        RDouble &area = nxyz[3];
        RDouble &vgn  = nxyz[4];

        RFloat &rm = prim[IR];
        RFloat &um = prim[IU];
        RFloat &vm = prim[IV];
        RFloat &wm = prim[IW];
        RFloat &pm = prim[IP];

        RFloat vn    = xfn * um + yfn * vm + zfn * wm - vgn;
        RFloat eigen = fabs(vn) + sqrt(gama * pm / rm);

        return eigen * area;
    }

    __device__ void MXDQ_STD(const int nchem, RDouble *prim, RDouble nxs, RDouble nys, RDouble nzs, RDouble area,
                             RDouble vgn, RDouble gama, RDouble *dq, RDouble *f, int nm, int nLaminar, RDouble radius,
                             int ipn)
    {
        //! The relevant variables related to vgn has modified by LiPeng on Feb 19, 2019
        using namespace IDX;
        //!using namespace CHEMICAL_SPACE;

        RDouble rm = prim[IR];
        RDouble um = prim[IU];
        RDouble vm = prim[IV];
        RDouble wm = prim[IW];
        RDouble pm = prim[IP];

        RDouble v2 = um * um + vm * vm + wm * wm;
        RDouble c2 = fabs(gama * pm / rm);
        RDouble cm = sqrt(c2);

        RDouble lmd1, lmd2, lmd3, x1, x2, cits;

        //! The variable cits equals to Ur=U-Ub, which can be referred to the formula (A.5) in appendix A of the PHengLEI Theory manual.
        cits = nxs * um + nys * vm + nzs * wm - vgn;

        RDouble ctgm = cits;
        RDouble cmgm = cm;

        lmd1 = ctgm;
        lmd2 = ctgm + cmgm;
        lmd3 = ctgm - cmgm;

        RDouble ipnrad = ipn * radius / area;

        lmd1 = half * (lmd1 + ipnrad);
        lmd2 = half * (lmd2 + ipnrad);
        lmd3 = half * (lmd3 + ipnrad);

        //!This term is actually equal to 0.0
        x1 = (lmd1 + lmd1 - lmd2 - lmd3) / (c2 + c2);
        //!This term is actually equal to 0.5
        x2 = (lmd2 - lmd3) / (cm + cm);

        RDouble dc, dh, c2dc;
        //! dc=-b1,b1 is referred to the formula (A.7) in appendix A of the PHengLEI Theory manual.
        //! dc   = cits * dq[IR] - nxs * dq[IRU] - nys * dq[IRV] - nzs * dq[IRW];//!the previous form
        dc = (cits + vgn) * dq[IR] - nxs * dq[IRU] - nys * dq[IRV]
             - nzs * dq[IRW]; //!Modified by LiPeng on Feb 21, 2019.
        c2dc = c2 * dc;

        //! Compute total enthalpy hm and coefficient dh,dh=b2,see formula (A.7) and (A.8) in appendix A of the PHengLEI Theory manual.
        RDouble hm;
        ComputeTotalEnthalpyAndDH(nchem, prim, gama, dq, hm, dh);

        //!This term is actually equal to (0.5 * dc)
        RDouble dhx1_dcx2 = dh * x1 + dc * x2;
        //!This term is actually equal to (0.5 * dh)
        RDouble c2dcx1_dhx2 = c2dc * x1 + dh * x2;

        //! The following expressions can be referred to formula (A.6) in appendix A of the PHengLEI Theory manual.
        f[IR]  = lmd1 * dq[IR] - dhx1_dcx2;
        f[IRU] = lmd1 * dq[IRU] - um * dhx1_dcx2 + nxs * c2dcx1_dhx2;
        f[IRV] = lmd1 * dq[IRV] - vm * dhx1_dcx2 + nys * c2dcx1_dhx2;
        f[IRW] = lmd1 * dq[IRW] - wm * dhx1_dcx2 + nzs * c2dcx1_dhx2;

        //! The previous form is not suitable for dynamic mesh.
        //! f[IRE] = lmd1 * dq[IRE] - hm * dhx1_dcx2 + cits * c2dcx1_dhx2;
        //! Modified by LiPeng on Feb 21, 2019.
        f[IRE] = lmd1 * dq[IRE] - hm * dhx1_dcx2 + (cits + vgn) * c2dcx1_dhx2;
        for (int iLaminar = nm; iLaminar < nLaminar; ++iLaminar)
        {
            f[iLaminar] = lmd1 * dq[iLaminar] - prim[iLaminar] * dhx1_dcx2;
        }

        for (int iLaminar = 0; iLaminar < nLaminar; ++iLaminar) //!Modified by LiPeng on Feb. 25, 2019
        {
            f[iLaminar] *= area;
        }
    }
    __device__ void ComputeTotalEnthalpyAndDH(const int nchem, RFloat *prim, RFloat &gama, RFloat *dq,
                                              RFloat &totalEnthalpy, RFloat &dh)
    {
        //!@param[in ]: prim is an array of saving the primary variables, gama is the specific heat ratio of the mixture gas.
        //!@param[in ]: dq is an array of saving the difference of the conservative variables.
        //!@param[out]: totalEnthalpy is the total enthalpy of the mixture, dh denotes the coefficient of the vector.
        using namespace IDX;

        //!Obtain the primary variables
        RFloat &density  = prim[IR];
        RFloat &um       = prim[IU];
        RFloat &vm       = prim[IV];
        RFloat &wm       = prim[IW];
        RFloat &pressure = prim[IP];

        //!V2 is squared velocity£¬c1 is product coefficient of the variable dq[IR]£¬c2 = R*T/Mns - (gama-1)*ens£¬gama_1 = gama-1¡£
        RFloat V2, c1, c2, gama_1;
        V2     = um * um + vm * vm + wm * wm;
        gama_1 = gama - 1.0;
        c1     = half * gama_1 * V2;
        dh     = -gama_1 * (um * dq[IRU] + vm * dq[IRV] + wm * dq[IRW] - dq[IRE]);

        //!enthalpy is static enthalpy.
        RFloat enthalpy;
        if (nchem == 0) //!Compute the static enthalpy
        {
            //!h = cp*T = [gama/(gama-1)]*(R/M)*T = [gama/(gama-1)]*(p/rho)
            enthalpy = (gama / gama_1) * (pressure / density);
        }
        else //!Compute the static enthalpy of the multi-species model.
        {
            printf("Error: in ComputeTotalEnthalpyAndDH, nchem > 0\n");
            //!exit(1);
        }

        dh += c1 * dq[IR];
        //!h0 = h+v*v/2, total enthalpy equals to static enthalpy plus half of the squared velocity.
        totalEnthalpy = enthalpy + half * V2;
    }

    __device__ void Primitive2Conservative(const int nchem, const int nm, const int nl, const int ntmodel, RFloat *prim,
                                           RFloat gama, RFloat Tv, RFloat Te, RFloat *q)
    {
        using namespace IDX;
        RFloat  em;
        RFloat &density  = prim[IR];
        RFloat &um       = prim[IU];
        RFloat &vm       = prim[IV];
        RFloat &wm       = prim[IW];
        RFloat &pressure = prim[IP];

        //!Obtain the mass fraction of the last species whose label in the collection is ns-1.
        if (nchem == 1)
        {
            printf("Error: in Primitive2Conservative nchem=1\n");
            //!NormalizePrimitive(prim);
        }

        //!obtain the total internal energy Em.
        ComputeInternalEnergy(nchem, prim, gama, Tv, Te, em);

        q[IR]  = density;
        q[IRU] = density * um;
        q[IRV] = density * vm;
        q[IRW] = density * wm;
        q[IRE] = density * em;

        for (int m = nm; m < nl; ++m)
        {
            q[m] = density * prim[m];
        }
        if (nchem == 1)
        {
            printf("Error: in Primitive2Conservative nchem=1\n");
            //!q[nl] = density * prim[nl];
        }
        for (int i = 0; i < ntmodel - 1; i++) //!multi-temperature model
        {
            q[nl + nchem + i] = density * prim[nl + nchem + i];
        }
    }

    __device__ void ComputeInternalEnergy(const int nchem, RFloat *prim, RFloat gama_in, RFloat Tv, RFloat Te,
                                          RFloat &em)
    {
        using namespace IDX;
        RFloat tm, v2, omav, hm;

        RFloat &rm = prim[IR];
        RFloat &um = prim[IU];
        RFloat &vm = prim[IV];
        RFloat &wm = prim[IW];
        RFloat &pm = prim[IP];
        v2         = um * um + vm * vm + wm * wm;

        if (nchem == 0) //!Perfect gas
        {
            em = (pm / rm) / (gama_in - one) + half * v2;
        }
        else
        {
            printf("Error in ComputeInternalEnergy, nchem !=0\n");
            /*
            if (ntmodel == 1) //!Single-Temperature model
            {
                ComputeMolecularWeightReciprocal(prim, omav);
                tm = pm / ( coefficientOfStateEquation * rm * omav );
                ComputeSpeciesEnthalpy(tm, hintSpecies);
                ComputeMixtureByPrimitive(prim, hintSpecies, hm);
                em = hm - pm / rm + half * v2;
            }
                else //!multi-temperature model
                {
                RFloat *massFraction = new RFloat[this->numberOfSpecies];
                for (int i = 0; i < this->numberOfSpecies; i++)
                    massFraction[i] = prim[nm+i];

                RFloat R = this->GetUniversalGasConstant();
                //!To obtain the static enthalpy of the mixed gas.
                RFloat hinf = 0.0, Pe = 0.0, Ttr = 0.0, ce_div_Me;
                RFloat Mr = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFraction, ce_div_Me);

                if (ntmodel == 2) //!Two-Temperature Model
                {
                    //!RFloat Eve = prim[nl+nchem];
                    //!RFloat Tve = ComputeVibrationAndElectronTemperature(massFraction, Eve); //!The vibration-electron temperature.
                    RFloat Tve = Tv;
                    Pe = rm * Tve * R * ce_div_Me;
                    //!To obtain the translation-rotation temperature where the pressure is known.
                    Ttr = (pm - Pe) / (rm * R * Mr);
                    //!Compute the static enthalpy.
                    hinf = ComputeMixedGasEnthalpy(Ttr, Tve, Tve, massFraction);
                }
                        else if (ntmodel == 3) //!Three-Temperature Model
                {
                    //!RFloat Ev = prim[nl+nchem];
                    //!RFloat Ee = prim[nl+nchem+1];
                    //!RFloat Tv = ComputeVibrationTemperature(massFraction, Ev); //!The vibrationtemperature.
                    //!RFloat Te = ComputeElectronTemperature(massFraction, Ee); //!The electron temperature.
                    Pe = rm * Te * R * ce_div_Me;
                    //!To obtain the translation-rotation temperature where the pressure is known.
                    Ttr = (pm - Pe) / (rm * R * Mr);
                    //!Compute the static enthalpy.
                    hinf = ComputeMixedGasEnthalpy(Ttr, Tv, Te, massFraction);
                }
                delete [] massFraction;

                //!Compute the total internal energy.
                em = hinf + 0.5 * v2 - pm / rm;
            }
            */
        }
    }
} //! namespace GPUKernels
