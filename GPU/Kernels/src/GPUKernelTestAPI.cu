#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cmath>
#include "BasicDeviceVariables.h"
#include "Precision.h"
#include "GPUKernelTest.h"
#include "GPUKernelTestAPI.h"
#include "cudaErrorHandle.h"

using namespace PHSPACE;
//! test the GPUVencatLimiter7Loop1
void TestGPUVencatLimiter7Eps(RFloat *epsCell, RFloat *d_epsCell, const int nTotalCell)
{
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    //! update absolute error
    //! double abCompare = 1.0e-14; //!20w
    //! double abCompare = 1.0e-13; //!20w
    //! double abCompare = 1.0e-12; //!20w
    //! double abCompare = 1.0e-11;

    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;

    double err         = 0.0;
    double err_rlt     = 0.0;
    double flagCompare = 1.0e-8;
    size_t epsSize     = nTotalCell * sizeof(RFloat);

    RFloat *g_epsCell = new RFloat[nTotalCell];
    HANDLE_API_ERR(cudaMemcpy(g_epsCell, d_epsCell, epsSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotalCell; i++)
    {
        //! err= fabs(g_epsCell[i]- epsCell[i]);
        abErr = fabs(g_epsCell[i] - epsCell[i]);
        rtErr = abErr / Max(fabs(epsCell[i]), smallLimit);
        //! err_rlt= err/(epsCell[i])
        //! if(err> flagCompare)
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("%d term error in TestGPUVencatLimiter7Eps, epsCell on host: %.30e\t "
                   " on device:%.30e, abErr = %.30e, rtErr = %.30e\n",
                   i, epsCell[i], g_epsCell[i], abErr, rtErr);
            printf("abErr = %.30e, rtErr= %.30e\n", abErr, rtErr);
            exit(EXIT_FAILURE);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("%d term error in TestGPUVencatLimiter7Eps, epsCell on host: %.30e\t "
                   " on device:%.30e, abErr = %.30e, rtErr = %.30e\n",
                   i, epsCell[i], g_epsCell[i], abErr, rtErr);
            printf("abErr = %.30e, rtErr= %.30e\n", abErr, rtErr);
            exit(EXIT_FAILURE);
        }
    }
#ifdef UNITTESTOUTPUT
    printf("GPUVencatLimiter7Loop1 is successful!\n");
#endif
    delete[] g_epsCell;
}

void TestGPUVencatLimiter7Boun(RFloat *limit, int *left_cell_of_face, const int iVariable, const int nBoundFace,
                               const int nTotal)
{
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;

    RFloat       err         = 0.0;
    const double flagCompare = 1.0e-40;

    int     le          = 0;
    RFloat *g_limit     = new RFloat[nTotal];
    RFloat *d_limit_tmp = d_LIMIT + iVariable * nTotal;
    size_t  limitSize   = nTotal * sizeof(RFloat);

    HANDLE_API_ERR(cudaMemcpy(g_limit, d_limit_tmp, limitSize, cudaMemcpyDeviceToHost));
    for (int j = 0; j < nBoundFace; j++)
    {
        le = left_cell_of_face[j];

        //! err= fabs(g_limit[le]- limit[le]);
        abErr           = fabs(g_limit[le] - limit[le]);
        rtErr           = abErr / Max(fabs(limit[le]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("%d term(face id:%d) error in GPUVencatLimiter7Boun, LIMIT on host: "
                   "%f\t  on device:%f\n",
                   j, le, limit[le], g_limit[le]);
            exit(EXIT_FAILURE);
        }
        if (err > flagCompare)
        {
            printf("%d term(face id:%d) error in GPUVencatLimiter7Boun, LIMIT on host: "
                   "%f\t  on device:%f\n",
                   j, le, limit[le], g_limit[le]);
            exit(EXIT_FAILURE);
        }
    }

#ifdef UNITTESTOUTPUT
    printf("GPUVencatLimiter7Boun is successful!\n");
#endif
    delete[] g_limit;
}
