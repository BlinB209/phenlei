#include "GPUSpectrumRadius1eq.h"
#include "GPUBasicFunctions.h"
#include "GPUDeviceControl.h"
#if __CUDA_ARCH__ < 600
#include "GPUPortability.h"
#endif

namespace GPUKernels
{
    void CallGPUSpecTurbCell(const int nTotalCell, const int nBoundFace, const RFloat dualTimeSpectrumC1,
                             const RFloat dualTimeSpectrumC2, const RFloat cflturb, const double SMALL)
    {
        using namespace GPUMemory;
        using namespace GPUGeomVariables;
        using namespace GPUFlowVariables;
        int nTotal = nTotalCell + nBoundFace;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUSpecTurbCell, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSpecTurbCell, 0, gridSize, blockSize);
#endif
        GPUSpecTurbCell<<<gridSize, blockSize>>>(nTotalCell, nTotal, dualTimeSpectrumC1, dualTimeSpectrumC2, cflturb,
                                                 SMALL, d_vol, d_dt, d_spec_turb);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUSpecTurbCell(const int nTotalCell, const int nTotal, const RFloat dualTimeSpectrumC1,
                                    const RFloat dualTimeSpectrumC2, const RFloat cflturb, const double SMALL,
                                    const RDouble *vol, const RFloat *dt, RFloat *spec_turb)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        for (icell = i; icell < nTotalCell; icell += blockDim.x * gridDim.x)
        {
            //!spec_turb[ ISA ][ icell ] += dualTimeSpectrumC1 * vol[ icell ] + dualTimeSpectrumC2 / ( cflturb * dt[ icell ]+ SMALL );
            spec_turb[0 * nTotal + icell] +=
                dualTimeSpectrumC1 * vol[icell] + dualTimeSpectrumC2 / (cflturb * dt[icell] + SMALL);
        }
    }

    void CallGPUSpecTurbMatTurbFaces(const int nTotalFace, const int nBoundFace, const int nTotalCell,
                                     const double SMALL, const RFloat oreynolds, const double osigma, const double cb2,
                                     const bool IfModOrgSaMode)
    {
        using namespace GPUMemory;
        using namespace GPUGeomVariables;
        using namespace GPUFlowVariables;

        int gridSize        = 1;
        int blockSize       = 1;
        int loopLen         = nTotalFace;
        int regsPerThread   = 66;
        int residentBlockNo = 4;

        size_t dsMemPerThread = 0; //! calculate by the intended shared memory use in the kernel
        size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUSpecTurbMatTurbFaces, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                         dsMemPerBlock, gridSize, blockSize, GPUProp);

//!KernelLaunchPara((void*)GPUReConstructFaceValueLoop2, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSpecTurbMatTurbFaces, dsMemPerBlock, gridSize, blockSize);
#endif
#ifndef CUDADEBUGATOMIC
        GPUSpecTurbMatTurbFaces<<<gridSize, blockSize, dsMemPerBlock>>>(
            nTotalFace, nBoundFace, nTotalCell, SMALL, oreynolds, osigma, IfModOrgSaMode, d_left_cell_of_face,
            d_right_cell_of_face, d_xfn, d_yfn, d_zfn, d_area, d_q_ns, d_gama_ns, d_q_turb_proxy, d_visl, d_vol, cb2,
            d_spec_turb, d_mat_turbl, d_mat_turbr);
#else
        GPUSpecTurbMatTurbFaces<<<1, 1, dsMemPerBlock>>>(nTotalFace, nBoundFace, nTotalCell, SMALL, oreynolds, osigma,
                                                         IfModOrgSaMode, d_left_cell_of_face, d_right_cell_of_face,
                                                         d_xfn, d_yfn, d_zfn, d_area, d_q_ns, d_gama_ns, d_q_turb_proxy,
                                                         d_visl, d_vol, cb2, d_spec_turb, d_mat_turbl, d_mat_turbr);

#endif

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUSpecTurbMatTurbFaces(const int nTotalFace, const int nBoundFace, const int nTotalCell,
                                            const double SMALL, const RFloat oreynolds, const double osigma,
                                            const bool IfModOrgSaMode, const int *left_cell_of_face,
                                            const int *right_cell_of_face, const RDouble *xfn, const RDouble *yfn,
                                            const RDouble *zfn, const RDouble *area, const RFloat *q,
                                            const RFloat *gama, const RFloat *q_turb, const RFloat *visl,
                                            const RDouble *vol, const RFloat cb2, RFloat *spec_turb, RFloat *mat_turbl,
                                            RFloat *mat_turbr)
    {
        int    i      = blockIdx.x * blockDim.x + threadIdx.x;
        int    iface  = 0;
        int    nTotal = nTotalCell + nBoundFace;
        int    le, re;
        RFloat rl, ul, vl, wl, pl, cl, vnl, orl, abs_vnl;
        RFloat rr, ur, vr, wr, pr, cr, vnr, orr, abs_vnr;
        RFloat nx, ny, nz, vn, abs_vn, volume, ovolume, ns2;
        RFloat nu_l, nu_r, nueff_l, nueff_r, nueff;
        double coeff_le, coeff_re;
        for (iface = i; iface < nTotalFace; iface += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[iface];
            re = right_cell_of_face[iface];

            if (iface < nBoundFace) re = iface + nTotalCell;

            nx = xfn[iface] * area[iface];
            ny = yfn[iface] * area[iface];
            nz = zfn[iface] * area[iface];

            //!calculate velocity at the cell interface
            rl = q[0 * nTotal + le];
            ul = q[1 * nTotal + le];
            vl = q[2 * nTotal + le];
            wl = q[3 * nTotal + le];
            pl = q[4 * nTotal + le];
            cl = sqrt(gama[le] * pl / rl);

            vnl     = xfn[iface] * ul + yfn[iface] * vl + zfn[iface] * wl;
            abs_vnl = GPUABS(vnl);

            rr = q[0 * nTotal + re];
            ur = q[1 * nTotal + re];
            vr = q[2 * nTotal + re];
            wr = q[3 * nTotal + re];
            pr = q[4 * nTotal + re];
            cr = sqrt(gama[re] * pr / rr);

            vnr     = xfn[iface] * ur + yfn[iface] * vr + zfn[iface] * wr;
            abs_vnr = GPUABS(vnr);

            vn = 0.5 * (vnl + vnr);
            //!abs_vn = ABS(vn);
            abs_vn = 0.5 * (abs_vnl + abs_vnr);
//!spec_turb[ISA][le] += area[iface] * half * abs_vn;
//!spec_turb[ISA][re] += area[iface] * half * abs_vn;
#if __CUDA_ARCH__ < 600
            atomicAddTest(spec_turb + 0 * nTotal + le, area[iface] * 0.5 * abs_vn);
            atomicAddTest(spec_turb + 0 * nTotal + re, area[iface] * 0.5 * abs_vn);
#else
            atomicAdd(spec_turb + 0 * nTotal + le, area[iface] * 0.5 * abs_vn);
            atomicAdd(spec_turb + 0 * nTotal + re, area[iface] * 0.5 * abs_vn);
#endif
            //!¼´Ê¹abs_vn»òÕßABS(vnl)ÕâÑùµÄÑ¡ÔñÒ²»áÔì³É²î±ð£¬ËäÈ»×îÖÕ¶Ô½á¹ûµÄÓ°Ïì»¹ÄÑÒÔÔ¤ÁÏ
            //!mat_turbl[ISA][iface]   = area[iface] * half * ( - vnl - abs_vn );
            //!mat_turbr[ISA][iface]   = area[iface] * half * (   vnr - abs_vn );
            mat_turbl[0 * nTotal + iface] = area[iface] * 0.5 * (-vnl - abs_vn);
            mat_turbr[0 * nTotal + iface] = area[iface] * 0.5 * (vnr - abs_vn);

            //!spec_turb[ISA][le] += area[iface] * half * ( abs_vn + cl );
            //!spec_turb[ISA][re] += area[iface] * half * ( abs_vn + cr );
            //!¼´Ê¹abs_vn»òÕßABS(vnl)ÕâÑùµÄÑ¡ÔñÒ²»áÔì³É²î±ð£¬ËäÈ»×îÖÕ¶Ô½á¹ûµÄÓ°Ïì»¹ÄÑÒÔÔ¤ÁÏ
            //!mat_turbl[ISA][iface]   = area[iface] * half * ( - vnl - ( ABS(abs_vn) + cl ) );
            //!mat_turbr[ISA][iface]   = area[iface] * half * (   vnr - ( ABS(abs_vn) + cr ) );

            //!calculate jacobians of the diffusion operator at the cell face i+1/2
            rl = q[0 * nTotal + le];
            rr = q[0 * nTotal + re];

            orl = 1.0 / (rl + SMALL);
            orr = 1.0 / (rr + SMALL);

            nu_l = q_turb[0 * nTotal + le];
            nu_r = q_turb[0 * nTotal + re];

            nueff_l = visl[le] * orl + GPUABS(nu_l);
            nueff_r = visl[re] * orr + GPUABS(nu_r);

            nueff = 0.5 * (nueff_l + nueff_r);

            volume  = 0.5 * (vol[le] + vol[re]);
            ovolume = 1.0 / volume;

            ns2 = area[iface] * area[iface];

            if (IfModOrgSaMode)
            {
//!spec_turb[ISA][le] +=   oreynolds * ovolume * ( 1.0 + cb2 ) * ( osigma * nueff * ns2 );
//!spec_turb[ISA][re] +=   oreynolds * ovolume * ( 1.0 + cb2 ) * ( osigma * nueff * ns2 );
#if __CUDA_ARCH__ < 600
                atomicAddTest(spec_turb + 0 * nTotal + le, oreynolds * ovolume * (1.0 + cb2) * (osigma * nueff * ns2));
                atomicAddTest(spec_turb + 0 * nTotal + re, oreynolds * ovolume * (1.0 + cb2) * (osigma * nueff * ns2));
#else
                atomicAdd(spec_turb + 0 * nTotal + le, oreynolds * ovolume * (1.0 + cb2) * (osigma * nueff * ns2));
                atomicAdd(spec_turb + 0 * nTotal + re, oreynolds * ovolume * (1.0 + cb2) * (osigma * nueff * ns2));
#endif

                mat_turbl[0 * nTotal + iface] += -oreynolds * ovolume * (1.0 + cb2) * (osigma * nueff * ns2);
                mat_turbr[0 * nTotal + iface] += -oreynolds * ovolume * (1.0 + cb2) * (osigma * nueff * ns2);

                coeff_le = cb2 * osigma * nueff_l;
                coeff_re = cb2 * osigma * nueff_r;
                mat_turbl[0 * nTotal + iface] += oreynolds * ovolume * ns2 * coeff_le;
                mat_turbr[0 * nTotal + iface] += oreynolds * ovolume * ns2 * coeff_re;
            }
            else
            {
//!spec_turb[ISA][le] +=   oreynolds * ovolume * ( osigma * nueff * ns2 );
//!spec_turb[ISA][re] +=   oreynolds * ovolume * ( osigma * nueff * ns2 );
#if __CUDA_ARCH__ < 600
                atomicAddTest(spec_turb + 0 * nTotal + le, oreynolds * ovolume * (osigma * nueff * ns2));
                atomicAddTest(spec_turb + 0 * nTotal + re, oreynolds * ovolume * (osigma * nueff * ns2));
#else
                atomicAdd(spec_turb + 0 * nTotal + le, oreynolds * ovolume * (osigma * nueff * ns2));
                atomicAdd(spec_turb + 0 * nTotal + re, oreynolds * ovolume * (osigma * nueff * ns2));
#endif

                mat_turbl[0 * nTotal + iface] += -oreynolds * ovolume * (osigma * nueff * ns2);
                mat_turbr[0 * nTotal + iface] += -oreynolds * ovolume * (osigma * nueff * ns2);
            }
        } //!end of for
    }
} //! namespace GPUKernels
