#include "GPUSourceFlux_1eq_Original.h"
#include "GPUBasicFunctionsPart2.h"
#include "GPUBasicFunctions.h"
#include "GPUDeviceControl.h"

namespace GPUKernels
{
    void CallGPUSpecTurbInit(const int nTotal, const int n_turb)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUSpecTurbInit, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSpecTurbInit, 0, gridSize, blockSize);
#endif
        GPUSpecTurbInit<<<gridSize, blockSize>>>(nTotal, n_turb, d_spec_turb);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUSpecTurbInit(const int nTotal, const int n_turb, RFloat *spec_turb)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        for (icell = i; icell < nTotal; icell += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < n_turb; ++m)
            {
                spec_turb[m * nTotal + icell] = 0.0;
            }
        }
    }

    void CallGPUSourceFlux_1eq_Original(const int nTotalCell, const int nBoundFace, const int isplt, const double SMALL,
                                        const double oreynolds, const RFloat cv13, const RFloat sac2, const RFloat sac3,
                                        const RFloat cw2, const RFloat cw36, const RFloat or6, const RFloat cb1,
                                        const RFloat cb2s, const RFloat cw1k, const RFloat sdilim,
                                        const bool IfModOrgSaMode, const RFloat rkap2, const RFloat cw1)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize        = 1;
        int blockSize       = 1;
        int loopLen         = nTotalCell;
        int regsPerThread   = 106;
        int residentBlockNo = 2;

        size_t dsMemPerThread = 0; //! calculate by the intended shared memory use in the kernel
        size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUSourceFlux_1eq_Original, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                         dsMemPerBlock, gridSize, blockSize, GPUProp);

//!KernelLaunchPara((void*)GPUReConstructFaceValueLoop2, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSourceFlux_1eq_Original, dsMemPerBlock, gridSize, blockSize);
#endif
        GPUSourceFlux_1eq_Original<<<gridSize, blockSize, dsMemPerBlock>>>(
            nTotalCell, nBoundFace, isplt, SMALL, oreynolds, d_dveldx_proxy, d_dveldy_proxy, d_dveldz_proxy,
            d_q_turb_proxy, d_q_ns, d_dqdx_proxy, d_dqdy_proxy, d_dqdz_proxy, d_visl, d_vol, d_walldist, cv13, sac2,
            sac3, cw2, cw36, or6, cb1, cb2s, cw1k, sdilim, IfModOrgSaMode, rkap2, cw1, d_spec_turb, d_res_turb);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUSourceFlux_1eq_Original(const int nTotalCell, const int nBoundFace, const int isplt,
                                               const double SMALL, const double oreynolds, const RFloat *dVdx,
                                               const RFloat *dVdy, const RFloat *dVdz, const RFloat *q_turb,
                                               const RFloat *q, const RFloat *dqdx, const RFloat *dqdy,
                                               const RFloat *dqdz, const RFloat *visl, const RDouble *vol,
                                               const RDouble *lengthScale, const RFloat cv13, const RFloat sac2,
                                               const RFloat sac3, const RFloat cw2, const RFloat cw36, const RFloat or6,
                                               const RFloat cb1, const RFloat cb2s, const RFloat cw1k,
                                               const RFloat sdilim, const bool IfModOrgSaMode, const RFloat rkap2,
                                               const RFloat cw1, RFloat *spec_turb, RFloat *res_turb)
    {
        int    nTotal = nTotalCell + nBoundFace;
        RFloat rs, str, std, ostd, nuet_rs;
        RFloat r, r5, g, g6, fw, grd2;
        //!RFloat nuet,volume, olam,cmachfix;
        RFloat nuet, volume, olam;
        //!RFloat work1,work2,work3,work4;
        RFloat work4;
        RFloat dsdnu, dfv2dk, dfwdg, dgdr, drdnu, dfwdnu;
        //!RFloat d2,od2,xsi,xsi2,xsi3,fv1,fv2,ft2;
        RFloat d2, od2, xsi, xsi3, fv1, fv2, ft2;
        RFloat sbar, omega;
        //!RFloat production, diffsource, wall_dest, wall_prod, source;
        RFloat production, diffsource, wall_dest, source;
        RFloat prod, dest, prodp, destp, prde, prdep, part1, part2;
        //!RFloat sac2,sac3;
        int i     = blockIdx.x * blockDim.x + threadIdx.x;
        int icell = 0;
        for (icell = i; icell < nTotalCell; icell += blockDim.x * gridDim.x)
        {
            //!work1 = dudx[icell] + dvdy[icell] + dwdz[icell];
            //!    work1 = dVdx[0 * nTotal + icell] + dVdy[1 * nTotal + icell] + dVdz[2 * nTotal + icell];
            //!work2 = SQR( dudx[icell] ) + SQR( dvdy[icell] ) + SQR( dwdz[icell] );
            //!    work2 = GPUSQR( dVdx[0 * nTotal + icell] ) + GPUSQR( dVdy[1 * nTotal + icell] ) + GPUSQR( dVdz[2 * nTotal + icell] );
            /*
                work3 = SQR( dudy[icell] + dvdx[icell] ) +
                        SQR( dvdz[icell] + dwdy[icell] ) +
                        SQR( dwdx[icell] + dudz[icell] );
            */
            //!    work3 = GPUSQR( dVdy[0 * nTotal + icell] + dVdx[1 * nTotal + icell] ) +
            //!        GPUSQR( dVdz[1 * nTotal + icell] + dVdy[2 * nTotal + icell] ) +
            //!        GPUSQR( dVdx[2* nTotal + icell] + dVdz[0 * nTotal + icell] );
            /*
                   work4 = sqrt( SQR( dwdy[icell] - dvdz[icell] ) +
                              SQR( dudz[icell] - dwdx[icell] ) +
                              SQR( dvdx[icell] - dudy[icell] ) );
            */
            work4 = sqrt(GPUSQR(dVdy[2 * nTotal + icell] - dVdz[1 * nTotal + icell])
                         + GPUSQR(dVdz[0 * nTotal + icell] - dVdx[2 * nTotal + icell])
                         + GPUSQR(dVdx[1 * nTotal + icell] - dVdy[0 * nTotal + icell]));
            //!noting that in fact, isplt = 2!!
            if (isplt == 1)
            {
            }
            else if (isplt == 4)
            {
            }
            else if (isplt == 5)
            {
            }
            else if (isplt == 6 || isplt == 7)
            {
            }

            str = work4;
            //!nuet     = q_turb[ISA][icell];
            nuet   = q_turb[0 * nTotal + icell];
            volume = vol[icell];
            //!cmachfix = nuet * ( dudx[icell] + dvdy[icell] + dwdz[icell] );
            //!cmachfix = 0.0;
            //!olam     = q[IR][icell] / ( visl[icell] + SMALL );
            olam = q[0 * nTotal + icell] / (visl[icell] + SMALL);
            d2   = lengthScale[icell] * lengthScale[icell];
            od2  = 1.0 / d2;

            xsi  = nuet * olam + SMALL;
            xsi3 = xsi * xsi * xsi;
            fv1  = xsi3 / (xsi3 + cv13);
            fv2  = 1.0 - xsi / (1.0 + xsi * fv1);

            //!¶ÔstdÏÞÖÆÊÇ·Ç³£ÓÐ±ØÒªµÄ
            //!std   = MAX(str + nuet * fv2 * ky2 * oreynolds, SMALL);
            //!oostd = SIGN(one,std) / ( ABS(std) + SMALL );

            //!ÕâÀï½«oreynoldsÎüÊÕÈërs
            rs      = nuet * rkap2 * od2 * oreynolds;
            nuet_rs = nuet * rs;

            sbar  = rs * fv2;
            omega = str;

            if (sbar >= -sac2 * omega)
            {
                std = omega + sbar;
            }
            else
            {
                std = omega + omega * (sac2 * sac2 * omega + sac3 * sbar) / ((sac3 - 2 * sac2) * omega - sbar);
            }

            //!std   i  = MAX(str + rs * fv2, 0.3*str);
            ostd = 1.0 / (std + SMALL);
            r    = rs * ostd;

            r = GPUMIN(r, static_cast<RFloat>(10.0));
#ifndef DEBUGPOW
            r5 = pow(r, 5);
#else
            r5 = r * 1.0;
#endif
            g = r + cw2 * (r5 * r - r);

#ifndef DEBUGPOW
            g6 = pow(g, 6);
#else
            g6 = g * 1.0;
#endif

#ifndef DEBUGPOW
            fw = g * pow((static_cast<RFloat>(1.0) + cw36) / (g6 + cw36), or6);
#else
            fw = g * 1.0;
#endif

            //!calculate "negative" production term in order to assure a zero
            //!solution in the laminar region

            //!ftrans = 0.0;
            //!xsi2   = xsi * xsi;
            //!ft2    = ct3 * exp( - ct4 * xsi2 ) * ftrans;
            ft2 = 0.0;

            //!grd2 = SQR( dqdx[ISA][icell] ) + SQR( dqdy[ISA][icell] ) + SQR( dqdz[ISA][icell] );
            grd2 =
                GPUSQR(dqdx[0 * nTotal + icell]) + GPUSQR(dqdy[0 * nTotal + icell]) + GPUSQR(dqdz[0 * nTotal + icell]);

            production = cb1 * (1.0 - ft2) * std * nuet; //!cb1  * omega * muet + cb1  * muet_rs * fv2;
            diffsource = cb2s * grd2 * oreynolds;        //!cb2s * rho * gradnue2;
            wall_dest  = (cw1k * fw - cb1 * ft2)
                        * nuet_rs; //!cw1k * muet_rs * fw£¬»¯¼òºó£¬ºÍTauÊÖ²áÉÏµÄDÒ»ÖÂ: cw1 * fw * (miu/d)^2;

            diffsource = GPUMIN(sdilim * production, diffsource);

            source = production - wall_dest;

            if (IfModOrgSaMode)
            {
                //!source  += diffsource;            //! bell 20130513 delete for modsa
            }
            else
            {
                source += diffsource;
            }

#ifndef DEBUGPOW
            dfv2dk = (3.0 * cv13 * pow(fv1 / (xsi + SMALL), 2) - 1.0) / pow(1.0 + xsi * fv1, 2);
#else
            dfv2dk = (3.0 * cv13 * 1.0 - 1.0) / (1.0 + xsi * fv1);
#endif

            dsdnu = oreynolds * rkap2 * od2 * (fv2 + xsi * dfv2dk);

            dfwdg  = fw / (g + SMALL) * (1.0 - g6 / (g6 + cw36));
            dgdr   = 1.0 + cw2 * (6.0 * r5 - 1.0);
            drdnu  = oreynolds * rkap2 * od2 * ostd * (1.0 - nuet * ostd * dsdnu);
            dfwdnu = dfwdg * dgdr * drdnu;

            //!dft2dnu = - two * ct3 * ct4 * xsi * exp ( - ct4 * xsi2 ) * olam * ftrans;

            prod = cb1 * (1.0 - ft2) * std;
            dest = (cw1 * fw - cb1 * ft2 * rkap2) * nuet * od2 * oreynolds;

            prodp = cb1 * ((1.0 - ft2) * dsdnu);
            destp = (cw1 * (fw + dfwdnu * nuet) - cb1 * ft2 * rkap2) * od2 * oreynolds;

            prde  = prod - dest;
            prdep = prodp - destp;

            part1 = -0.5 * (prde - GPUABS(prde));
            part2 = -0.5 * (prdep - GPUABS(prdep)) * nuet;

            spec_turb[0 * nTotal + icell] += (part1 + part2) * volume;
            res_turb[0 * nTotal + icell] += (source)*volume;
        }
    }
} //! namespace GPUKernels
