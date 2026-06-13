#include "GPUGetVisFaceValue.h"
#include "GPUDeviceControl.h"
#include "cudaErrorHandle.h"
#include "Geo_SimpleBC.h"

using namespace std;
using namespace PHENGLEI;
namespace GPUKernels
{
    void CallGPUGetVisFaceValue(Grid *gridIn, Param_NSSolverUnstruct *parameters, const int neqn, const int nst,
                                const int ned)
    {
        using namespace GPUMemory;
        using namespace GPUGeomVariables;
        using namespace GPUFlowVariables;
        using namespace IDX;

        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nBoundFace = grid->GetNBoundFace();
        int nTotalCell = grid->GetNTotalCell();

        RDouble refGama = 1.4;
        GlobalDataBase::GetData("refGama", &refGama, PHDOUBLE, 1);
        const RDouble gama0 = refGama;
        //!RDouble refGama = parameters->GetRefGama();

        RDouble oprl = parameters->GetoPrandtlLaminar();
        RDouble oprt = parameters->GetoPrandtlTurbulence();

        RDouble referenceMachNumber = parameters->GetRefMachNumber();

        RDouble twall = parameters->GetWallTemperature();

        RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

        RDouble m2 = referenceMachNumber * referenceMachNumber;

        int numberOfSpecies = parameters->GetNumberOfSpecies();
        int nchem           = parameters->GetChemicalFlag();
        /* original kernel size
        gridSize=2092;
        blockSize=640;
        */
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = ned - nst;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUGetVisFaceValueMulMut, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUGetVisFaceValueMulMut, 0, gridSize, blockSize);
#endif

        GPUGetVisFaceValueMulMut<<<gridSize, blockSize>>>(nst, ned, nBoundFace, nTotalCell, d_left_cell_of_face,
                                                          d_right_cell_of_face, d_visl, d_vist, d_mul, d_mut);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        int nMid;
        //!modified by zhang xi for bug
        loopLen = ned - nst;
        KernelLaunchPara((void *)GPUGetVisFaceValuePrimTm, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUGetVisFaceValuePrimTm, 0, gridSize, blockSize);
#endif

        GPUGetVisFaceValuePrimTm<<<gridSize, blockSize>>>(nst, ned, nBoundFace, nTotalCell, neqn, nchem, d_SEG_LEN, ITT,
                                                          d_left_cell_of_face, d_right_cell_of_face, d_t_proxy, d_q_ns,
                                                          d_mul, d_mut, oprl, oprt, gama0, m2, d_prim, d_tm, d_kcp);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        //!modified end
        nMid = nBoundFace;
        if (ned < nBoundFace) nMid = ned;
        if (nst < nBoundFace)
        {
            RFloat tw = twall / refDimensionalTemperature;
            loopLen   = nMid - nst;
            KernelLaunchPara((void *)GPUGetVisFaceValueOnBound, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUGetVisFaceValueOnBound, 0, gridSize, blockSize);
#endif

            GPUGetVisFaceValueOnBound<<<gridSize, blockSize>>>(
                nst, nMid, nBoundFace, nTotalCell, neqn, d_SEG_LEN, d_boundaryType, d_left_cell_of_face,
                d_right_cell_of_face, d_q_ns, d_xtn, d_ytn, d_ztn, d_prim, d_tm, tw, PHENGLEI::INTERFACE,
                PHENGLEI::SOLID_SURFACE);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }
    }

    __global__ void GPUGetVisFaceValueMulMut(const int nst, const int ned, const int nBoundFace, const int nTotalCell,
                                             const int *left_cell_of_face, const int *right_cell_of_face,
                                             const RFloat *visl, const RFloat *vist, RFloat *mul, RFloat *mut)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        int le, re, j;
        for (iface = i + nst; iface < ned; iface += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[iface];
            re = right_cell_of_face[iface];
            if (iface < nBoundFace) re = iface + nTotalCell;
            j = iface - nst;

            mul[j] = 0.5 * (visl[le] + visl[re]);
            mut[j] = 0.5 * (vist[le] + vist[re]);
        }
    }

    __global__ void GPUGetVisFaceValuePrimTm(const int nst, const int ned, const int nBoundFace, const int nTotalCell,
                                             const int neqn, const int nchem, const int len, const int ITT,
                                             const int *left_cell_of_face, const int *right_cell_of_face,
                                             const RFloat *t, const RFloat *q, const RFloat *mul, const RFloat *mut,
                                             const double oprl, const double oprt, double gama0, RFloat m2,
                                             RFloat *prim, RFloat *tmid, RFloat *kcp)
    {
        int    i = blockDim.x * blockIdx.x + threadIdx.x;
        int    le, re, j, m;
        RFloat cp;
        //!int nMid = nBoundFace;
        int nTotal = nTotalCell + nBoundFace;
        //!if ( ned < nBoundFace ) nMid = ned;

        for (int iface = nst + i; iface < ned; iface += blockDim.x * gridDim.x)
        {
            re = right_cell_of_face[iface];
            le = left_cell_of_face[iface];
            if (iface < nBoundFace) re = iface + nTotalCell;
            j = iface - nst;
            for (m = 0; m < neqn; ++m)
            {
                //!prim[m][j] = 0.5 * (q[m][le] + q[m][re]);
                prim[m * len + j] = 0.5 * (q[m * nTotal + le] + q[m * nTotal + re]);
            }

            tmid[ITT * len + j] = 0.5 * (t[ITT * nTotal + le] + t[ITT * nTotal + re]);

            if (nchem == 1)
            {
                /*Becasue in the case, nchem is not equal to zero. Nothing happens.*/
            }
            else
            {
                cp     = 1.0 / ((gama0 - 1.0) * m2);
                kcp[j] = (mul[j] * oprl + mut[j] * oprt) * cp;
            }
        }
    }
    /* all faces on boundary */
    __global__ void GPUGetVisFaceValuePrimTm_S1(const int nst, const int ned, const int nBoundFace,
                                                const int nTotalCell, const int neqn, const int nchem, const int len,
                                                const int ITT, const int *left_cell_of_face,
                                                const int *right_cell_of_face, const RFloat *t, const RFloat *q,
                                                const RFloat *mul, const RFloat *mut, const double oprl,
                                                const double oprt, double gama0, RFloat m2, RFloat *prim, RFloat *tmid,
                                                RFloat *kcp)
    {
        int    i = blockDim.x * blockIdx.x + threadIdx.x;
        int    le, re, j, m;
        RFloat cp;
        //!int nMid = nBoundFace;
        int nTotal = nTotalCell + nBoundFace;
        //!if ( ned < nBoundFace ) nMid = ned;

        for (int iface = nst + i; iface < ned; iface += blockDim.x * gridDim.x)
        {
            //!re = right_cell_of_face[ iface ];
            le = left_cell_of_face[iface];
            re = iface + nTotalCell;
            j  = iface - nst;
            for (m = 0; m < neqn; ++m)
            {
                //!prim[m][j] = 0.5 * (q[m][le] + q[m][re]);
                prim[m * len + j] = 0.5 * (q[m * nTotal + le] + q[m * nTotal + re]);
            }

            tmid[ITT * len + j] = 0.5 * (t[ITT * nTotal + le] + t[ITT * nTotal + re]);

            if (nchem == 1)
            {
                /*Becasue in the case, nchem is not equal to zero. Nothing happens.*/
            }
            else
            {
                cp     = 1.0 / ((gama0 - 1.0) * m2);
                kcp[j] = (mul[j] * oprl + mut[j] * oprt) * cp;
            }
        }
    }
    /* all face interior */
    __global__ void GPUGetVisFaceValuePrimTm_S2(const int nst, const int nMid, const int ned, const int nBoundFace,
                                                const int nTotalCell, const int neqn, const int nchem, const int len,
                                                const int ITT, const int *left_cell_of_face,
                                                const int *right_cell_of_face, const RFloat *t, const RFloat *q,
                                                const RFloat *mul, const RFloat *mut, const double oprl,
                                                const double oprt, double gama0, RFloat m2, RFloat *prim, RFloat *tmid,
                                                RFloat *kcp)
    {
        int    i = blockDim.x * blockIdx.x + threadIdx.x;
        int    le, re, j, m;
        RFloat cp;
        //!int nMid = nBoundFace;
        int nTotal = nTotalCell + nBoundFace;
        //!if ( ned < nBoundFace ) nMid = ned;

        for (int iface = nMid + i; iface < ned; iface += blockDim.x * gridDim.x)
        {
            re = right_cell_of_face[iface];
            le = left_cell_of_face[iface];
            j  = iface - nst;
            for (m = 0; m < neqn; ++m)
            {
                //!prim[m][j] = 0.5 * (q[m][le] + q[m][re]);
                prim[m * len + j] = 0.5 * (q[m * nTotal + le] + q[m * nTotal + re]);
            }

            tmid[ITT * len + j] = 0.5 * (t[ITT * nTotal + le] + t[ITT * nTotal + re]);

            if (nchem == 1)
            {
                /*Becasue in the case, nchem is not equal to zero. Nothing happens.*/
            }
            else
            {
                cp     = 1.0 / ((gama0 - 1.0) * m2);
                kcp[j] = (mul[j] * oprl + mut[j] * oprt) * cp;
            }
        }
    }
    __global__ void GPUGetVisFaceValueOnBound(const int nst, const int nMid, const int nBoundFace, const int nTotalCell,
                                              const int neqn, const int len, const int *boundaryType,
                                              const int *left_cell_of_face, const int *right_cell_of_face,
                                              const RFloat *q, const RDouble *xtn, const RDouble *ytn,
                                              const RDouble *ztn, RFloat *prim, RFloat *tmid, RFloat tw,
                                              const int bctype_interface, const int bctype_solid)
    {
        int    i     = blockDim.x * blockIdx.x + threadIdx.x;
        int    iface = 0;
        int    le, re, bctype, j, m;
        RFloat tmp;
        int    nTotal = nBoundFace + nTotalCell;
        for (iface = nst + i; iface < nMid; iface += blockDim.x * gridDim.x)
        {
            bctype = boundaryType[iface];
            j      = iface - nst;

            le = left_cell_of_face[iface];
            re = right_cell_of_face[iface];

            if (bctype == bctype_interface) continue;

            for (m = 0; m < neqn; ++m)
            {
                //!tmp =  0.5 * ( q[m][le] + q[m][re] );
                tmp = 0.5 * (q[m * nTotal + le] + q[m * nTotal + re]);

                //!prim[m][j] = tmp;
                prim[m * len + j] = tmp;
            }

            if (bctype == bctype_solid)
            {
                prim[1 * len + j] = xtn[iface];
                prim[2 * len + j] = ytn[iface];
                prim[3 * len + j] = ztn[iface];

                if (tw > 0.0)
                {
                    //!tmid[0][j] = tw;
                    tmid[0 * len + j] = tw;
                }
            }
        }
    }

} //! namespace GPUKernels
