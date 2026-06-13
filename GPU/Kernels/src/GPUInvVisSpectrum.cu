#include "GPUInvVisSpectrum.h"
#include "GPUBasicFunctions.h"
#include "GPUDeviceControl.h"
#include "cudaErrorHandle.h"
#if __CUDA_ARCH__ < 600
#include "GPUPortability.h"
#endif
namespace GPUKernels
{
    void CallGPUSetFieldInvSpectrum(int nTotalCell, const RFloat value)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUControlVariables;
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUSetFieldInvSpectrum, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSetFieldInvSpectrum, 0, gridSize, blockSize);
#endif
        GPUSetFieldInvSpectrum<<<gridSize, blockSize>>>(nTotalCell, d_invSpectrumRadius, value);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUSetFieldVisSpectrum(int nTotalCell, const RFloat value)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUControlVariables;
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUSetFieldVisSpectrum, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSetFieldVisSpectrum, 0, gridSize, blockSize);
#endif
        GPUSetFieldVisSpectrum<<<gridSize, blockSize>>>(nTotalCell, d_visSpectrumRadius, value);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUSpectrumRadiusInviscid(Grid *gridIn)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        using namespace GPUControlVariables;
        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nTotalCell = grid->GetNTotalCell();
        int nTotalFace = grid->GetNTotalFace();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        //! int isUnsteady = parameters->GetIsUnsteady();
        //! int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

        //! RDouble *preconCoefficient = NULL;
        //! RDouble *timeCoefficientInverse = NULL;

        //! if (ifLowSpeedPrecon != 0)
        //! {
        //!     preconCoefficient = reinterpret_cast<RDouble *> (grid->GetDataPtr("preconCoefficient"));
        //!     if (isUnsteady)
        //!     {
        //!         timeCoefficientInverse = reinterpret_cast<RDouble*> (grid->GetDataPtr("timeCoefficientInverse"));
        //!     }
        //! }

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;
        cudaMemset(d_invSpectrumRadius, 0, nTotalCell * sizeof(RFloat));

        loopLen = nTotalFace;
        KernelLaunchPara((void *)GPUInvSpectrumRadiusCalculate, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUInvSpectrumRadiusCalculate, 0, gridSize, blockSize);
#endif
        GPUInvSpectrumRadiusCalculate<<<gridSize, blockSize>>>(nTotalCell, nTotal, nTotalFace, d_left_cell_of_face,
                                                               d_right_cell_of_face, d_xfn, d_yfn, d_zfn, d_area,
                                                               d_q_ns, d_gama_ns, d_invSpectrumRadius);
        //! #ifndef CUDADEBUGATOMIC
        //!     GPUInvSpectrumRadiusCalculate<<<gridSize, blockSize>>>(nTotalCell, nTotal, nTotalFace, d_left_cell_of_face, d_right_cell_of_face, d_xfn, d_yfn, d_zfn, d_area, d_q_ns, d_gama_ns, d_invSpectrumRadius);
        //! #else
        //!     GPUInvSpectrumRadiusCalculate<<<1, 1>>>(nTotalCell, nTotal, nTotalFace, d_left_cell_of_face, d_right_cell_of_face, d_xfn, d_yfn, d_zfn, d_area, d_q_ns, d_gama_ns, d_invSpectrumRadius);
        //! #endif
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    //!void CallGPUVisSpectrumRadiusCalculate(const int nTotalCell, const int nTotal, const int nTotalFace, const double refReNumber, const double SMALL)
    void CallGPUSpectrumRadiusViscous(Grid *gridIn, Param_NSSolverUnstruct *parameters)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        using namespace GPUControlVariables;

        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nTotalCell = grid->GetNTotalCell();
        int nTotalFace = grid->GetNTotalFace();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        const RDouble small = 1.0e-36f;

        RDouble refReNumber = parameters->GetRefReNumber();

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalFace;
        int blockSizeLimit = 0;

        cudaMemset(d_visSpectrumRadius, 0, nTotalCell * sizeof(RFloat));

        KernelLaunchPara((void *)GPUVisSpectrumRadiusCalculate, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUVisSpectrumRadiusCalculate, 0, gridSize, blockSize);
#endif
        GPUVisSpectrumRadiusCalculate<<<gridSize, blockSize>>>(
            nTotalCell, nTotal, nTotalFace, refReNumber, small, d_left_cell_of_face, d_right_cell_of_face, d_xcc, d_ycc,
            d_zcc, d_xfn, d_yfn, d_zfn, d_area, d_q_ns, d_visl, d_vist, d_visSpectrumRadius);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUSetFieldInvSpectrum(const int nTotalCell, RFloat *invSpectrum, const RFloat value)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        for (iface = i; iface < nTotalCell; iface += blockDim.x * gridDim.x)
        {
            invSpectrum[iface] = value;
        }
    }

    __global__ void GPUSetFieldVisSpectrum(const int nTotalCell, RFloat *visSpectrum, const RFloat value)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        for (iface = i; iface < nTotalCell; iface += blockDim.x * gridDim.x)
        {
            visSpectrum[iface] = value;
        }
    }

    __global__ void GPUInvSpectrumRadiusCalculate(const int nTotalCell, const int nTotal, const int nTotalFace,
                                                  const int *left_cell_of_face, int *right_cell_of_face,
                                                  const RDouble *xfn, const RDouble *yfn, const RDouble *zfn,
                                                  const RDouble *area, const RFloat *q, const RFloat *gama,
                                                  RFloat *invSpectrum)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iFace = 0;
        for (iFace = i; iFace < nTotalFace; iFace += blockDim.x * gridDim.x)
        {
            int le = left_cell_of_face[iFace];
            int re = right_cell_of_face[iFace];

            RFloat rl  = q[nTotal * 0 + le];
            RFloat ul  = q[nTotal * 1 + le];
            RFloat vl  = q[nTotal * 2 + le];
            RFloat wl  = q[nTotal * 3 + le];
            RFloat pl  = q[nTotal * 4 + le];
            RFloat cl  = sqrt(gama[le] * pl / rl);
            RFloat vnl = xfn[iFace] * ul + yfn[iFace] * vl + zfn[iFace] * wl;

            RFloat rr  = q[nTotal * 0 + re];
            RFloat ur  = q[nTotal * 1 + re];
            RFloat vr  = q[nTotal * 2 + re];
            RFloat wr  = q[nTotal * 3 + re];
            RFloat pr  = q[nTotal * 4 + re];
            RFloat cr  = sqrt(gama[re] * pr / rr);
            RFloat vnr = xfn[iFace] * ur + yfn[iFace] * vr + zfn[iFace] * wr;

            RFloat specificHeatRatio = 0.5 * (gama[le] + gama[re]);

            RFloat pm = 0.5 * (pl + pr);
            RFloat rm = 0.5 * (rl + rr);
            RFloat vn = 0.5 * (vnl + vnr);
            RFloat cm = sqrt(specificHeatRatio * pm / rm);

            //! if(ifLowSpeedPrecon != 0)
            //! {
            //!     RDouble c2 = specificHeatRatio * pm / rm;
            //!     RDouble preconCoeff = half * (preconCoefficient[le] + preconCoefficient[re]);

            //!     vn = half * vn * (1.0 + preconCoeff);
            //!     cm = half * sqrt(((1.0 - preconCoeff) * vn) * ((1.0 - preconCoeff) * vn) + 4.0 * preconCoeff * c2);

            //!     if (isUnsteady)
            //!     {
            //!         RDouble timeCoeff = half * (timeCoefficientInverse[le] + timeCoefficientInverse[re]);
            //!         cm *= timeCoeff;
            //!         vn *= timeCoeff;
            //!     }
            //! }

            RFloat inviscidSpectrumRadius = 0.5 * area[iFace] * (GPUABS(vn) + cm);

//!invSpectrum[le] += inviscidSpectrumRadius;
#if __CUDA_ARCH__ < 600
            atomicAddTest(invSpectrum + le, inviscidSpectrumRadius);
#else
            atomicAdd(invSpectrum + le, inviscidSpectrumRadius);
#endif
            if (re < nTotalCell)
            {
//!invSpectrum[re] += inviscidSpectrumRadius;
#if __CUDA_ARCH__ < 600
                atomicAddTest(invSpectrum + re, inviscidSpectrumRadius);
#else
                atomicAdd(invSpectrum + re, inviscidSpectrumRadius);
#endif
            }
        }
    }

    __global__ void GPUVisSpectrumRadiusCalculate(const int nTotalCell, const int nTotal, const int nTotalFace,
                                                  const double refReNumber, const double SMALL,
                                                  const int *left_cell_of_face, const int *right_cell_of_face,
                                                  const RDouble *xcc, const RDouble *ycc, const RDouble *zcc,
                                                  const RDouble *xfn, const RDouble *yfn, const RDouble *zfn,
                                                  const RDouble *area, const RDouble *q, const RFloat *visl,
                                                  const RFloat *vist, RFloat *visSpectrum)
    {
        int    i = blockDim.x * blockIdx.x + threadIdx.x;
        int    le, re, iFace;
        double visLaminar, visTurb;
        double dx, dy, dz, ds, viscousSpectrumRadius, viscosity, density, faceArea;
        for (iFace = i; iFace < nTotalFace; iFace += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[iFace];
            re = right_cell_of_face[iFace];

            dx = xcc[re] - xcc[le];
            dy = ycc[re] - ycc[le];
            dz = zcc[re] - zcc[le];
            ds = GPUABS(xfn[iFace] * dx + yfn[iFace] * dy + zfn[iFace] * dz);

            visLaminar = 0.5 * (visl[le] + visl[re]);
            visTurb    = 0.5 * (vist[le] + vist[re]);
            //!density    = 0.5 * (rho[le] + rho[re]);
            density   = 0.5 * (q[0 * nTotal + le] + q[0 * nTotal + re]);
            viscosity = visLaminar + visTurb;

            faceArea              = area[iFace];
            viscousSpectrumRadius = 2.0 * viscosity / (density * ds * refReNumber + SMALL);
            viscousSpectrumRadius *= 0.5 * faceArea;

//!visSpectrum[le] += viscousSpectrumRadius;
#if __CUDA_ARCH__ < 600
            atomicAddTest(visSpectrum + le, viscousSpectrumRadius);
#else
            atomicAdd(visSpectrum + le, viscousSpectrumRadius);
#endif
            if (re < nTotalCell)
            {
//!visSpectrum[re] += viscousSpectrumRadius;
#if __CUDA_ARCH__ < 600
                atomicAddTest(visSpectrum + re, viscousSpectrumRadius);
#else
                atomicAdd(visSpectrum + re, viscousSpectrumRadius);
#endif
            }
        }
    }
} //! namespace GPUKernels
