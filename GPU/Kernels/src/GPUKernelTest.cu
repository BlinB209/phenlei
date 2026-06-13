#include <stdlib.h>
#include <cmath>
#include "BasicDeviceVariables.h"
#include "GPUKernelTest.h"
#include "OutputDebug.h"
#include "Pointer.h"
#include "TemporaryOperations.h"

using namespace PHSPACE;
using namespace TemporaryOperations;
using namespace std;
void TestGPUGetQlQr(FaceProxy *face_proxy, const int nst, const int ned, const int SEG_LEN)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;

    int nl;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);
    int nchem;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
    int neqn = nl + nchem;

    int len = SEG_LEN;
    //! ql qr from device
    RFloat **g_ql_ns = GPUQl_NSCopyBack(d_ql_ns, neqn, len);
    RFloat **g_qr_ns = GPUQr_NSCopyBack(d_qr_ns, neqn, len);
    //! ql qr from host
    RFloat **ql = face_proxy->GetQL();
    RFloat **qr = face_proxy->GetQR();

    //! const double flagCompare = 1.0e-8;
    //! RFloat err= 0.0;
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    //! absolute error update
    //! double abCompare = 1.0e-15;
#ifdef UNITTESTOUTPUT
    printf("Test GPUGetQlQr...");
#endif
    for (int i = 0; i < neqn; i++)
    {
        for (int j = 0; j < (ned - nst); j++)
        {
            abErr = fabs(g_ql_ns[i][j] - ql[i][j]);
            rtErr = abErr / Max(fabs(ql[i][j]), smallLimit);
            //! if(err>flagCompare)
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPUReConstructFaceValueLoop1, ql on "
                       "host: %.30e\t on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, ql[i][j], g_ql_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPUReConstructFaceValueLoop1, ql on "
                       "host: %.30e\t on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, ql[i][j], g_ql_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            abErr = fabs(g_qr_ns[i][j] - qr[i][j]);
            rtErr = abErr / Max(fabs(qr[i][j]), smallLimit);
            //! if(err>flagCompare)
            checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPUReConstructFaceValueLoop1, qr on "
                       "host: %.30e\t on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, qr[i][j], g_qr_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPUReConstructFaceValueLoop1, qr on "
                       "host: %.30e\t on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, qr[i][j], g_qr_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("Kernel GPUGetQlQr operation is successful!\n");
#endif
    DelPointer2<RFloat>(g_ql_ns);
    DelPointer2<RFloat>(g_qr_ns);
}

void TestGPUGetGamaLR(FaceProxy *face_proxy, const int nst, const int ned, const int SEG_LEN)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;
    //! gamaL gamaR on device
    RFloat *g_gamaL_ns = GPUGamaL_NSCopyBack(d_gamaL_ns, SEG_LEN);
    RFloat *g_gamaR_ns = GPUGamaR_NSCopyBack(d_gamaR_ns, SEG_LEN);
    //! gamaL gamaR on host
    RFloat *h_gamaL_ns = face_proxy->GetGamaL();
    RFloat *h_gamaR_ns = face_proxy->GetGamaR();

    //! RFloat err= 0.0;
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    for (int j = 0; j < (ned - nst); j++)
    {
        abErr           = fabs(g_gamaL_ns[j] - h_gamaL_ns[j]);
        rtErr           = abErr / Max(fabs(h_gamaL_ns[j]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("%d term error in GPUGetGamaLR, gamaL on host: %e\t on device:%e\t, "
                   "abErr= %e\t, rtErr= %e\n",
                   j, h_gamaL_ns[j], g_gamaL_ns[j], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("%d term error in GPUGetGamaLR, gamaL on host: %e\t on device:%e\t, "
                   "abErr= %e\t, rtErr= %e\n",
                   j, h_gamaL_ns[j], g_gamaL_ns[j], abErr, rtErr);
            exit(EXIT_FAILURE);
        }

        abErr    = fabs(g_gamaR_ns[j] - h_gamaR_ns[j]);
        rtErr    = abErr / Max(fabs(h_gamaR_ns[j]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("%d term error in GPUGetGamaLR, gamaR on host: %e\t on device:%e\t, "
                   "abErr= %e\t, rtErr= %e\n",
                   j, h_gamaR_ns[j], g_gamaR_ns[j], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("%d term error in GPUGetGamaLR, gamaR on host: %e\t on device:%e\t, "
                   "abErr= %e\t, rtErr= %e\n",
                   j, h_gamaR_ns[j], g_gamaR_ns[j], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
    }
#ifdef UNITTESTOUTPUT
    printf("Kernel GPUGetGamaLR operation is successful!\n");
#endif
    delete[] g_gamaL_ns;
    delete[] g_gamaR_ns;
}
//! test the GPUReConstructFaceValueLoop
void TestGPUReConstructFaceValueLoop(RFloat **ql, RFloat **qr, RFloat *d_ql_ns, RFloat *d_qr_ns, const int d_neqn,
                                     const int d_SEG_LEN, const int nst, const int ned)
{
    using namespace PHSPACE;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    RFloat **g_ql_ns = NewPointer2<RFloat>(d_neqn, d_SEG_LEN);
    RFloat **g_qr_ns = NewPointer2<RFloat>(d_neqn, d_SEG_LEN);

    size_t qlrSize = d_neqn * d_SEG_LEN * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(g_ql_ns[0], d_ql_ns, qlrSize, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(g_qr_ns[0], d_qr_ns, qlrSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < d_neqn; i++)
    {
        for (int j = 0; j < (ned - nst); j++)
        {
            /*
      abErr= fabs(g_ql_ns[i][j]- ql[i][j]);
      if(err>flagCompare)
      if(err>flagCompare)
      {
              printf("%d, %d term error in TestGPUReConstructFaceValueLoop2, ql
      on host: %f\t on device:%f\n", i, j, ql[i][j], g_ql_ns[i][j]);
              exit(EXIT_FAILURE);
      }
      err= fabs(g_qr_ns[i][j]- qr[i][j]);
      if(err>flagCompare)
      {
              printf("%d, %d term error in TestGPUReConstructFaceValueLoop2, qr
      on host: %f\t on device:%f\n", i, j, qr[i][j], g_qr_ns[i][j]);
              exit(EXIT_FAILURE);
      }
      */
            abErr           = fabs(g_ql_ns[i][j] - ql[i][j]);
            rtErr           = abErr / Max(fabs(ql[i][j]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPUReConstructFaceValueLoop2, ql on "
                       "host: %.30e\t on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, ql[i][j], g_ql_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            //! if(err>flagCompare)
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPUReConstructFaceValueLoop2, ql on "
                       "host: %.30e\t on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, ql[i][j], g_ql_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            abErr    = fabs(g_qr_ns[i][j] - qr[i][j]);
            rtErr    = abErr / Max(fabs(qr[i][j]), smallLimit);
            checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPUReConstructFaceValueLoop2, qr on "
                       "host: %.30e\t on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, qr[i][j], g_qr_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            //! if(err>flagCompare)
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPUReConstructFaceValueLoop2, qr on "
                       "host: %.30e\t on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, qr[i][j], g_qr_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test TestGPUReConstructFaceValueLoop2 successfully" << endl;
#endif
    DelPointer2(g_ql_ns);
    DelPointer2(g_qr_ns);
}

//! test the GPUBoundaryQlQrFix
void TestGPUBoundaryQlQrFix(FaceProxy *face_proxy, const int nst, const int ned, const int SEG_LEN)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;

    int nl;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);
    int nchem;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
    int neqn = nl + nchem;

    RFloat **ql = face_proxy->GetQL();
    RFloat **qr = face_proxy->GetQR();

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    RFloat **g_ql_ns = NewPointer2<RFloat>(neqn, SEG_LEN);
    RFloat **g_qr_ns = NewPointer2<RFloat>(neqn, SEG_LEN);

    size_t qlrSize = neqn * SEG_LEN * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(g_ql_ns[0], d_ql_ns, qlrSize, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(g_qr_ns[0], d_qr_ns, qlrSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < neqn; i++)
    {
        for (int j = 0; j < (ned - nst); j++)
        {
            abErr           = fabs(g_ql_ns[i][j] - ql[i][j]);
            rtErr           = abErr / Max(fabs(ql[i][j]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPUBoundaryQlQrFix, ql on host: %.30e\t "
                       "on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, ql[i][j], g_ql_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            //! if(err>flagCompare)
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPUBoundaryQlQrFix, ql on host: %.30e\t "
                       "on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, ql[i][j], g_ql_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            abErr    = fabs(g_qr_ns[i][j] - qr[i][j]);
            rtErr    = abErr / Max(fabs(qr[i][j]), smallLimit);
            checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPUBoundaryQlQrFix, qr on host: %.30e\t "
                       "on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, qr[i][j], g_qr_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            //! if(err>flagCompare)
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPUBoundaryQlQrFix, qr on host: %.30e\t "
                       "on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, qr[i][j], g_qr_ns[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test TestGPUBoundaryQlQrFix successfully" << endl;
    cout << "!!!All Kernels In GetInvFaceValue are successfully" << endl;
#endif
    DelPointer2(g_ql_ns);
    DelPointer2(g_qr_ns);
}

//! ! test the GPURoe_Scheme_Old kernel function
void TestGPUInv_scheme(InviscidSchemeParameter *invSchemePara, RFloat *d_flux, const int nst, const int ned,
                       const int SEG_LEN)
{
    using namespace PHSPACE;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    int    neqn  = invSchemePara->GetNumberOfTotalEquation();

    RFloat **flux   = invSchemePara->GetFlux();
    RFloat **g_flux = NewPointer2<RFloat>(neqn, SEG_LEN);

    size_t fluxSize = neqn * SEG_LEN * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(g_flux[0], d_flux, fluxSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < neqn; i++)
    {
        for (int j = 0; j < (ned - nst); j++)
        {
            abErr           = fabs(g_flux[i][j] - flux[i][j]);
            rtErr           = abErr / Max(fabs(flux[i][j]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPURoe_Scheme_Old, flux on host: %.30e\t "
                       "on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, flux[i][j], g_flux[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            //! if(err>flagCompare)
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPURoe_Scheme_Old, flux on host: %.30e\t "
                       "on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, flux[i][j], g_flux[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPURoe_Scheme_Old successfully" << endl;
#endif
    DelPointer2(g_flux);
}
//! test the GPUInviscidFluxWrapLoop1
void TestGPUInviscidFluxWrapLoop1(RFloat **flux, RFloat *d_flux, const int nl, const int nst, const int ned,
                                  const int SEG_LEN)
{
    //! RFloat err= 0.0;
    //! const double flagCompare = 1.0e-8;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    RFloat **g_flux = NewPointer2<RFloat>(nl, SEG_LEN);

    size_t fluxSize = nl * SEG_LEN * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(g_flux[0], d_flux, fluxSize, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nl; i++)
    {
        for (int j = 0; j < (ned - nst); j++)
        {
            abErr           = fabs(g_flux[i][j] - flux[i][j]);
            rtErr           = abErr / Max(fabs(flux[i][j]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in GPUInviscidFluxWrapLoop1, flux on host: %e\t "
                       "on device:%e\t, abErr= %e\t, rtErr= %e\n",
                       i, j, flux[i][j], g_flux[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            //! if(err>flagCompare)
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in GPUInviscidFluxWrapLoop1, flux on host: %e\t "
                       "on device:%e\t, abErr= %e\t, rtErr= %e\n",
                       i, j, flux[i][j], g_flux[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUInviscidFluxWrapLoop1 successfully" << endl;
#endif

    DelPointer2(g_flux);
}
//! test the GPULoadFlux
void TestGPULoadFlux(Grid *grid_in, const int nst, const int ned)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    UnstructGrid *grid               = UnstructGridCast(grid_in);
    int           nBoundFace         = grid->GetNBoundFace();
    int           nTotalCell         = grid->GetNTotalCell();
    int           nTotal             = nTotalCell + nBoundFace;
    int           nl                 = GlobalDataBase::GetIntParaFromDB("nl");
    int          *left_cell_of_face  = grid->GetLeftCellOfFace();
    int          *right_cell_of_face = grid->GetRightCellOfFace();

    size_t resSize = nl * nTotal * sizeof(RFloat);

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    RFloat **g_res_ns = NewPointer2<RFloat>(nl, nTotal);
    HANDLE_API_ERR(cudaMemcpy(g_res_ns[0], d_res_ns, resSize, cudaMemcpyDeviceToHost));

    RFloat **res = reinterpret_cast<RFloat **>(grid->GetDataPtr("res"));
    //! cout<<"in TestGPULoadFlux, nst= "<<nst<<", ned= "<<ned<<endl;
    for (int i = nst; i < ned; ++i)
    {
        int le, j, re;
        le = left_cell_of_face[i];
        re = right_cell_of_face[i];
        for (int m = 0; m < nl; ++m)
        {
            //! err= fabs(g_res_ns[m][le]- res[m][le]);
            abErr           = fabs(g_res_ns[m][le] - res[m][le]);
            rtErr           = abErr / Max(fabs(res[m][le]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPULoadFlux by le index, res on host: "
                       "%f\t  on device:%f\n",
                       m, le, res[m][le], g_res_ns[m][le]);
                exit(EXIT_FAILURE);
            }
            //! if(err> flagCompare)
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPULoadFlux by le index, res on host: "
                       "%f\t  on device:%f\n",
                       m, le, res[m][le], g_res_ns[m][le]);
                exit(EXIT_FAILURE);
            }
            //! err= fabs(g_res_ns[m][re]- res[m][re]);
            abErr    = fabs(g_res_ns[m][re] - res[m][re]);
            rtErr    = abErr / Max(fabs(res[m][re]), smallLimit);
            checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPULoadFlux by r index, res on host: "
                       "%f\t  on device:%f\n",
                       m, le, res[m][re], g_res_ns[m][re]);
                exit(EXIT_FAILURE);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPULoadFlux by r index, res on host: "
                       "%f\t  on device:%f\n",
                       m, le, res[m][re], g_res_ns[m][re]);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("GPULoadFlux is successful!\n");
#endif

    DelPointer2(g_res_ns);
}
void TestGPUGetQlQrTurb(FaceProxy *face_proxy, const int neqn, const int nst, const int ned, const int SEG_LEN)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    //! update absolute error
    int len = SEG_LEN;
    //! ql qr from device
    RFloat **g_ql_turb = GPUQl_NSCopyBack(d_ql_turb, neqn, len);
    RFloat **g_qr_turb = GPUQr_NSCopyBack(d_qr_turb, neqn, len);
    //! ql qr from host
    RFloat **h_ql_turb = face_proxy->GetQL();
    RFloat **h_qr_turb = face_proxy->GetQR();

    //! const double flagCompare = 1.0e-8;
    //! RFloat err= 0.0;
    for (int i = 0; i < neqn; i++)
    {
        for (int j = 0; j < (ned - nst); j++)
        {
            abErr           = fabs(g_ql_turb[i][j] - h_ql_turb[i][j]);
            rtErr           = abErr / Max(fabs(h_ql_turb[i][j]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in GPUGetQlQr in TestGPUGetQlQrTurb : ql, on "
                       "host: %.30e\t  device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, h_ql_turb[i][j], g_ql_turb[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            //! if(err>flagCompare)
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in GPUGetQlQr in TestGPUGetQlQrTurb : ql, on "
                       "host: %.30e\t  device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, h_ql_turb[i][j], g_ql_turb[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }

            abErr    = fabs(g_qr_turb[i][j] - h_qr_turb[i][j]);
            rtErr    = abErr / Max(fabs(h_qr_turb[i][j]), smallLimit);
            checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in GPUGetQlQr in TestGPUGetQlQrTurb : qr on "
                       "host: %.30e\t on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, h_qr_turb[i][j], g_qr_turb[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            //! if(err>flagCompare)
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in GPUGetQlQr in TestGPUGetQlQrTurb : qr on "
                       "host: %.30e\t on device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, h_qr_turb[i][j], g_qr_turb[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("GPUGetQlQrTurb successfully!\n");
#endif
    DelPointer2<RFloat>(g_ql_turb);
    DelPointer2<RFloat>(g_qr_turb);
}

void TestGPUFlux(RFloat **flux, RFloat *d_flux, const int nl, const int nst, const int ned, const int SEG_LEN)
{
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    RFloat **g_flux = NewPointer2<RFloat>(nl, SEG_LEN);

    size_t fluxSize = nl * SEG_LEN * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(g_flux[0], d_flux, fluxSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nl; i++)
    {
        for (int j = 0; j < (ned - nst); j++)
        {
            abErr           = fabs(g_flux[i][j] - flux[i][j]);
            rtErr           = abErr / Max(flux[i][j], smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in flux for NNDFlux,  on host: %.30e\t on "
                       "device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, flux[i][j], g_flux[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            //! if(err>flagCompare)
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in flux for NNDFlux,  on host: %.30e\t on "
                       "device:%.30e\t, abErr= %.30e\t, rtErr= %.30e\n",
                       i, j, flux[i][j], g_flux[i][j], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }

    DelPointer2(g_flux);
}

//! test the GPULoadFluxTurb
void TestGPULoadFluxTurb(RFloat **res_turb, RFloat *d_res_turb, const int *left_cell_of_face, const int nTotal,
                         const int neqn, const int nst, const int nMid)
{
    using namespace GPUMemory;
    using namespace GPUControlVariables;
    if (!IsInviscidFluxLoad) return;
    int loop = nMid - nst;
    if (0 == loop) return;
    size_t resSize = neqn * nTotal * sizeof(RFloat);

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    RFloat **g_res_turb = NewPointer2<RFloat>(neqn, nTotal);
    HANDLE_API_ERR(cudaMemcpy(g_res_turb[0], d_res_turb, resSize, cudaMemcpyDeviceToHost));

    //! cout<<"in TestGPULoadFlux, nst= "<<nst<<", ned= "<<ned<<endl;
    for (int i = nst; i < nMid; ++i)
    {
        int le, j, re;
        le = left_cell_of_face[i];
        for (int m = 0; m < neqn; ++m)
        {
            //! err= fabs(g_res_turb[m][le]- res_turb[m][le]);
            abErr = fabs(g_res_turb[m][le] - res_turb[m][le]);
            rtErr = abErr / Max(res_turb[m][le], smallLimit);
            //! if(err> flagCompare)
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPULoadFluxTurb by le index, res on "
                       "host: %f\t  on device:%f\n",
                       m, le, res_turb[m][le], g_res_turb[m][le]);
                exit(EXIT_FAILURE);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPULoadFluxTurb by le index, res on "
                       "host: %f\t  on device:%f\n",
                       m, le, res_turb[m][le], g_res_turb[m][le]);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("GPULoadFluxTurb is successful!\n");
#endif

    DelPointer2(g_res_turb);
}
//! test CallGPUCompGamaAndTField
void TestGPUCompGamaAndTField(Grid *grid_in)
{
    using namespace PHSPACE;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;

    UnstructGrid *grid = UnstructGridCast(grid_in);
    RFloat      **t    = reinterpret_cast<RFloat **>(grid->GetDataPtr("t"));
    RFloat       *gama = reinterpret_cast<RFloat *>(grid->GetDataPtr("gama"));

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal     = nTotalCell + nBoundFace;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    RFloat *g_gama = new RFloat[nTotal];

    size_t gamaSize = nTotal * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(g_gama, d_gama_ns, gamaSize, cudaMemcpyDeviceToHost));

#ifdef UNITTESTOUTPUT
    printf("Test GPUCompGamaAndTField ...\n");
#endif

    for (int i = 0; i < nTotal; i++)
    {
        abErr = fabs(g_gama[i] - gama[i]);
        rtErr = abErr / Max(fabs(gama[i]), smallLimit);
        //! if(err> flagCompare)
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("%d term error in Gama, gama on host: %.30e\t  on device:%.30e\t, "
                   "abErr= %.30e\t, rtErr= %.30e\n",
                   i, gama[i], g_gama[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("%d term error in Gama, gama on host: %.30e\t  on device:%.30e\t, "
                   "abErr= %.30e\t, rtErr= %.30e\n",
                   i, gama[i], g_gama[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
    }
    //! delete g_gama;
    delete[] g_gama;

    int    ntmodel = d_ntmodel; //! this can  be changed by different model
    int    ITT     = 0;
    size_t tSize   = ntmodel * nTotal * sizeof(RFloat);
    //! cout<<"neqn= "<<neqn<<" nTotal= "<<nTotal<<endl;
    //! RFloat **g_t_proxy= NewPointer2<RFloat>(ntmodel, nTotal);
    RFloat *g_t_proxy = new RFloat[nTotal];
    HANDLE_API_ERR(cudaMemcpy(g_t_proxy, d_t_proxy, tSize, cudaMemcpyDeviceToHost));

    for (int j = 0; j < nTotal; j++)
    {
        abErr = fabs(g_t_proxy[j] - t[ITT][j]);
        rtErr = abErr / Max(fabs(t[ITT][j]), smallLimit);
        //! if(err> flagCompare)
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("%d, %d term error in t, t on host: %.30e\t  on device:%.30e\t, , "
                   "abErr= %.30e\t, rtErr= %.30e\n",
                   ITT, j, t[ITT][j], g_t_proxy[j], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("%d, %d term error in t, t on host: %.30e\t  on device:%.30e\t, , "
                   "abErr= %.30e\t, rtErr= %.30e\n",
                   ITT, j, t[ITT][j], g_t_proxy[j], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
    }
    delete[] g_t_proxy;

    //!    DelPointer2(g_t_proxy);
#ifdef UNITTESTOUTPUT
    printf("TestGPUCompGamaAndTField is successful!\n");
#endif
}
void TestGPUBoundary(RFloat **q)
{
    using namespace PHSPACE;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;

#ifdef UNITTESTOUTPUT
    printf("Test GPUBoundary ...\n");
#endif

    //! const double flagCompare = 1.0e-8;
    RFloat abErr  = 0.0;
    RFloat rtErr  = 0.0;
    int    neqn   = d_neqn_ns;
    int    nTotal = d_nTotal;

    size_t   qSize  = neqn * nTotal * sizeof(RFloat);
    RFloat **g_q_ns = NewPointer2<RFloat>(neqn, nTotal);
    HANDLE_API_ERR(cudaMemcpy(g_q_ns[0], d_q_ns, qSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; ++i)
    {
        for (int m = 0; m < neqn; ++m)
        {
            abErr           = fabs(g_q_ns[m][i] - q[m][i]);
            rtErr           = abErr / Max(fabs(q[m][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in GPUBoundary , q on host: %.30e,  on "
                       "device:%.30e, abErr= %.30e, rtErr= %.30e\n",
                       m, i, q[m][i], g_q_ns[m][i], abErr, rtErr);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                if (processNum > 1)
                    printf("In globalRank %d, localRank %d, %d, %d term error in "
                           "GPUBoundary , q on host: %.30e,  on device:%.30e, abErr= %.30e, "
                           "rtErr= %.30e\n",
                           globalRank, localRank, m, i, q[m][i], g_q_ns[m][i], abErr, rtErr);
                else
                    printf("%d, %d term error in GPUBoundary , q on host: %.30e,  on "
                           "device:%.30e, abErr= %.30e, rtErr= %.30e\n",
                           m, i, q[m][i], g_q_ns[m][i], abErr, rtErr);

                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUBoundary is successful!\n");
#endif

    DelPointer2(g_q_ns);
}
//! why not RFloat??
double Max(double a, double b) { return a > b ? a : b; }

void TestGPULoadQ(RFloat **q)
{
    using namespace PHSPACE;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;

    RFloat abErr  = 0.0;
    RFloat rtErr  = 0.0;
    int    neqn   = d_neqn_ns;
    int    nTotal = d_nTotal;

    size_t   qSize     = neqn * nTotal * sizeof(RFloat);
    RFloat **g_q_proxy = NewPointer2<RFloat>(neqn, nTotal);
    HANDLE_API_ERR(cudaMemcpy(g_q_proxy[0], d_q_proxy, qSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; ++i)
    {
        for (int m = 0; m < neqn; ++m)
        {
            abErr           = fabs(g_q_proxy[m][i] - q[m][i]);
            rtErr           = abErr / Max(fabs(q[m][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in  TestGPULoadQ, q_proxy on host: %.30e\t  on "
                       "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                       m, i, q[m][i], g_q_proxy[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in  TestGPULoadQ, q_proxy on host: %.30e\t  on "
                       "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                       m, i, q[m][i], g_q_proxy[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPULoadQ is successful!\n");
#endif

    DelPointer2(g_q_proxy);
}
void TestGPUFillField(RFloat **q)
{
    using namespace PHSPACE;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;
    using namespace GPUControlVariables;
    RFloat  abErr  = 0.0;
    RFloat  rtErr  = 0.0;
    int     neqn   = d_neqn_ns;
    int     nTotal = d_nTotal;
    RFloat *d_q;
    /*
          if(IsFillFieldBeforeStage)
          {
                #ifdef UNITTESTOUTPUT
                  printf("In FillField Before Stage\n");
                #endif
                  d_q= d_q_proxy;
          }
          else
          {
                #ifdef UNITTESTOUTPUT
                  printf("In FillField during Stage\n");
                #endif
                  d_q= d_q_proxy_tmp;
          }
  */
    if (d_FillField_ns == "q_proxyToq_proxy_temp")
    {
#ifdef UNITTESTOUTPUT
        printf("In FillField Before Stage\n");
#endif
        d_q = d_q_proxy;
    }
    else if (d_FillField_ns == "q_proxy_tempToq_proxy")
    {
#ifdef UNITTESTOUTPUT
        printf("In FillField during Stage\n");
#endif
        d_q = d_q_proxy_tmp;
    }
    else if (d_FillField_ns == "resTodq")
    {
#ifdef UNITTESTOUTPUT
        printf("In NSSolverUnstruct::SolveLUSGSForward\n");
#endif
        d_q = d_dq_ns;
    }
    else
    {
        cout << "Warning: d_FillField_ns =" << d_FillField_ns << endl;
        return;
    }
    size_t   qSize     = neqn * nTotal * sizeof(RFloat);
    RFloat **g_q_proxy = NewPointer2<RFloat>(neqn, nTotal);
    HANDLE_API_ERR(cudaMemcpy(g_q_proxy[0], d_q, qSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; ++i)
    {
        for (int m = 0; m < neqn; ++m)
        {
            abErr           = fabs(g_q_proxy[m][i] - q[m][i]);
            rtErr           = abErr / Max(fabs(q[m][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPUFillField, q_proxy_tmp on host: "
                       "%.30e\t  on device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                       m, i, q[m][i], g_q_proxy[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPUFillField, q_proxy_tmp on host: "
                       "%.30e\t  on device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                       m, i, q[m][i], g_q_proxy[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUFillField is successful!\n");
#endif
    DelPointer2(g_q_proxy);
}

void TestGPULoadResiduals(RFloat **res)
{
    using namespace PHSPACE;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;
    using namespace GPUControlVariables;

    RFloat abErr  = 0.0;
    RFloat rtErr  = 0.0;
    int    neqn   = d_neqn_ns;
    int    nTotal = d_nTotal;

    size_t   qSize    = neqn * nTotal * sizeof(RFloat);
    RFloat **g_res_ns = NewPointer2<RFloat>(neqn, nTotal);
    HANDLE_API_ERR(cudaMemcpy(g_res_ns[0], d_res_ns, qSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; ++i)
    {
        for (int m = 0; m < neqn; ++m)
        {
            abErr           = fabs(g_res_ns[m][i] - res[m][i]);
            rtErr           = abErr / Max(fabs(res[m][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPULoadResiduals, res on host: %e\t  on "
                       "device:%e\t, abErr= %e, rtErr= %e\n",
                       m, i, res[m][i], g_res_ns[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPULoadResiduals, res on host: %e\t  on "
                       "device:%e\t, abErr= %e, rtErr= %e\n",
                       m, i, res[m][i], g_res_ns[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPULoadResiduals is successful!\n");
#endif
    DelPointer2(g_res_ns);
}
void TestGPULhs(RFloat **dq)
{
    using namespace PHSPACE;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;
    using namespace GPUControlVariables;

    RFloat abErr  = 0.0;
    RFloat rtErr  = 0.0;
    int    neqn   = d_neqn_ns;
    int    nTotal = d_nTotal;

    size_t   qSize   = neqn * nTotal * sizeof(RFloat);
    RFloat **g_dq_ns = NewPointer2<RFloat>(neqn, nTotal);
    HANDLE_API_ERR(cudaMemcpy(g_dq_ns[0], d_dq_ns, qSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; ++i)
    {
        for (int m = 0; m < neqn; ++m)
        {
            abErr           = fabs(g_dq_ns[m][i] - dq[m][i]);
            rtErr           = abErr / Max(fabs(dq[m][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPULhs, res on host: %.30e\t  on "
                       "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                       m, i, dq[m][i], g_dq_ns[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPULhs, res on host: %.30e\t  on "
                       "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                       m, i, dq[m][i], g_dq_ns[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPULhs is successful!\n");
#endif
    DelPointer2(g_dq_ns);
}
void TestGPUUpdateFlowFieldM1(RFloat **q)
{
    using namespace PHSPACE;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;
    using namespace GPUControlVariables;
    RFloat abErr  = 0.0;
    RFloat rtErr  = 0.0;
    int    neqn   = d_neqn_ns;
    int    nTotal = d_nTotal;

    size_t   qSize  = neqn * nTotal * sizeof(RFloat);
    RFloat **g_q_ns = NewPointer2<RFloat>(neqn, nTotal);
    HANDLE_API_ERR(cudaMemcpy(g_q_ns[0], d_q_ns, qSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; ++i)
    {
        for (int m = 0; m < neqn; ++m)
        {
            abErr           = fabs(g_q_ns[m][i] - q[m][i]);
            rtErr           = abErr / Max(fabs(q[m][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestUpdateFlowFieldM1, q_proxy_tmp on host: "
                       "%.30e\t  on device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                       m, i, q[m][i], g_q_ns[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestUpdateFlowFieldM1, q_proxy_tmp on host: "
                       "%.30e\t  on device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                       m, i, q[m][i], g_q_ns[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUUpdateFlowFieldM1 is successful!\n");
#endif
    DelPointer2(g_q_ns);
}
//! test the GPULoadFlux
void TestGPUZeroResiduals(RFloat **res, const int nl, const int nTotal)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t resSize = nl * nTotal * sizeof(RFloat);

    RFloat **g_res_ns = NewPointer2<RFloat>(nl, nTotal);
    HANDLE_API_ERR(cudaMemcpy(g_res_ns[0], d_res_ns, resSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; ++i)
    {
        for (int m = 0; m < nl; ++m)
        {
            abErr           = fabs(g_res_ns[m][i] - res[m][i]);
            rtErr           = abErr / Max(fabs(res[m][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in TestGPUZeroResiduals , res on host: %e\t  on "
                       "device:%e\t, abErr= %e, rtErr= %e\n",
                       m, i, res[m][i], g_res_ns[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in TestGPUZeroResiduals , res on host: %e\t  on "
                       "device:%e\t, abErr= %e, rtErr= %e\n",
                       m, i, res[m][i], g_res_ns[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUZeroResiduals is successful!\n");
#endif
    DelPointer2(g_res_ns);
}
void TestGPUDtCFL(RFloat *dt, const int nTotalCell)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t dSize = nTotalCell * sizeof(RFloat);

    RFloat *g_dt = new RFloat[nTotalCell];
    HANDLE_API_ERR(cudaMemcpy(g_dt, d_dt, dSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotalCell; ++i)
    {
        abErr = fabs(g_dt[i] - dt[i]);
        rtErr = abErr / Max(fabs(dt[i]), smallLimit);

        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("%d term error in TestGPUDtCFL, dt on host: %e\t  on device:%e\t, "
                   "abErr= %e, rtErr= %e\n",
                   i, dt[i], g_dt[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }

        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("%d term error in TestGPUDtCFL, dt on host: %e\t  on device:%e\t, "
                   "abErr= %e, rtErr= %e\n",
                   i, dt[i], g_dt[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUDtCFL is successful!\n");
#endif
    delete[] g_dt;
}
void TestGPUReduceMaxTimeStep(RFloat *dt, const int nTotalCell)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    size_t dSize = nTotalCell * sizeof(RFloat);

    RFloat *g_dt = new RFloat[nTotalCell];
    HANDLE_API_ERR(cudaMemcpy(g_dt, d_dt, dSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotalCell; ++i)
    {
        abErr           = fabs(g_dt[i] - dt[i]);
        rtErr           = abErr / Max(fabs(dt[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("%d term error in TestGPUReduceMaxTimeStep, dt on host: %.30e\t  on "
                   "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                   i, dt[i], g_dt[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }

        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("%d term error in TestGPUReduceMaxTimeStep, dt on host: %.30e\t  on "
                   "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                   i, dt[i], g_dt[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUReduceMaxTimeStep is successful!\n");
#endif
    delete[] g_dt;
}

void TestGPUSetGhostCell(RFloat *dt, const int nTotalCell)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    size_t dSize = nTotalCell * sizeof(RFloat);

    RFloat *g_dt = new RFloat[nTotalCell];
    HANDLE_API_ERR(cudaMemcpy(g_dt, d_dt, dSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotalCell; ++i)
    {
        abErr           = fabs(g_dt[i] - dt[i]);
        rtErr           = abErr / Max(fabs(dt[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("%d term error in TestGPUSetGhostCell, dt on host: %.30e\t  on "
                   "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                   i, dt[i], g_dt[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }

        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("%d term error in TestGPUSetGhostCell, dt on host: %.30e\t  on "
                   "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                   i, dt[i], g_dt[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUSetGhostCell is successful!\n");
#endif
    delete[] g_dt;
}
void TestGPUComputeMinTimeStep(RFloat minDt, RFloat maxDt)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    int    n     = 1;
    size_t dSize = n * sizeof(RFloat);

    RFloat *g_minDt = new RFloat[n];
    RFloat *g_maxDt = new RFloat[n];
    HANDLE_API_ERR(cudaMemcpy(g_minDt, d_minDt, dSize, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(g_maxDt, d_maxDt, dSize, cudaMemcpyDeviceToHost));

    abErr           = fabs(g_minDt[0] - minDt);
    rtErr           = abErr / Max(fabs(minDt), smallLimit);
    double checknan = abErr;
    if (std::isnan(checknan))
    {
        printf("Error: TestGPUComputeMinTimeStep, minDt on host: %e\t  on "
               "device:%e\t, abErr= %e, rtErr= %e\n",
               minDt, g_minDt[0], abErr, rtErr);
        exit(EXIT_FAILURE);
    }

    if (abErr > abCompare && rtErr > rtCompare)
    {
        printf("Error: TestGPUComputeMinTimeStep, minDt on host: %e\t  on "
               "device:%e\t, abErr= %e, rtErr= %e\n",
               minDt, g_minDt[0], abErr, rtErr);
        exit(EXIT_FAILURE);
    }

    abErr    = fabs(g_maxDt[0] - maxDt);
    rtErr    = abErr / Max(fabs(maxDt), smallLimit);
    checknan = abErr;
    if (std::isnan(checknan))
    {
        printf("Error: TestGPUComputeMinTimeStep, maxDt on host: %e\t  on "
               "device:%e\t, abErr= %e, rtErr= %e\n",
               maxDt, g_maxDt[0], abErr, rtErr);
        exit(EXIT_FAILURE);
    }

    if (abErr > abCompare && rtErr > rtCompare)
    {
        printf("Error: TestGPUComputeMinTimeStep, maxDt on host: %e\t  on "
               "device:%e\t, abErr= %e, rtErr= %e\n",
               maxDt, g_maxDt[0], abErr, rtErr);
        exit(EXIT_FAILURE);
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUComputeMinTimeStep is successful!\n");
#endif
    delete[] g_minDt;
    delete[] g_maxDt;
}

void TestGPUCompViscousCoef(RFloat *visl, const int nTotal)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    size_t dSize = nTotal * sizeof(RFloat);

    RFloat *g_visl = new RFloat[nTotal];
    HANDLE_API_ERR(cudaMemcpy(g_visl, d_visl, dSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; ++i)
    {
        abErr = fabs(g_visl[i] - visl[i]);
        rtErr = abErr / Max(fabs(visl[i]), smallLimit);

        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("%d term error in TestGPUCompViscousCoef, visl on host: %.30e\t  on "
                   "device:%30e\t, abErr= %.30e, rtErr= %.30e\n",
                   i, visl[i], g_visl[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }

        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("%d term error in TestGPUCompViscousCoef, visl on host: %.30e\t  on "
                   "device:%30e\t, abErr= %.30e, rtErr= %.30e\n",
                   i, visl[i], g_visl[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUCompViscousCoef is successful!\n");
#endif
    delete[] g_visl;
}
void TestGPUGPUTurbViscosity(double vistmax, double vistmin)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    int    n     = 1;
    size_t dSize = n * sizeof(double);

    RFloat *g_vistmin = new RFloat[n];
    RFloat *g_vistmax = new RFloat[n];
    HANDLE_API_ERR(cudaMemcpy(g_vistmin, d_vistmin, dSize, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(g_vistmax, d_vistmax, dSize, cudaMemcpyDeviceToHost));

    abErr = fabs(g_vistmin[0] - vistmin);
    rtErr = abErr / Max(fabs(vistmin), smallLimit);

    double checknan = abErr;
    if (std::isnan(checknan))
    {
        printf("Error: TestGPUGPUTurbViscosity, vistmin on host: %.30e\t  on "
               "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
               vistmin, g_vistmin[0], abErr, rtErr);
        exit(EXIT_FAILURE);
    }

    if (abErr > abCompare && rtErr > rtCompare)
    {
        printf("Error: TestGPUGPUTurbViscosity, vistmin on host: %.30e\t  on "
               "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
               vistmin, g_vistmin[0], abErr, rtErr);
        exit(EXIT_FAILURE);
    }

    abErr = fabs(g_vistmax[0] - vistmax);
    rtErr = abErr / Max(fabs(vistmax), smallLimit);

    checknan = abErr;
    if (std::isnan(checknan))
    {
        printf("Error: TestGPUGPUTurbViscosity, vistmax on host: %.30e\t  on "
               "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
               vistmax, g_vistmax[0], abErr, rtErr);
        exit(EXIT_FAILURE);
    }

    if (abErr > abCompare && rtErr > rtCompare)
    {
        printf("Error: TestGPUGPUTurbViscosity, vistmax on host: %.30e\t  on "
               "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
               vistmax, g_vistmax[0], abErr, rtErr);
        exit(EXIT_FAILURE);
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUGPUTurbViscosity  for vistmin vistmax is successful!\n");
#endif
    delete[] g_vistmin;
    delete[] g_vistmax;
}
void TestGPUTurbViscosity(const RFloat *org_vist, const int nTotal)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat  abErr       = 0.0;
    RFloat  rtErr       = 0.0; //! 20w
    size_t  sizeVisCoef = nTotal * sizeof(RFloat);
    RFloat *h_vist      = new RFloat[nTotal];
    HANDLE_API_ERR(cudaMemcpy(h_vist, d_vist, sizeVisCoef, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; i++)
    {
        abErr = fabs(h_vist[i] - org_vist[i]);
        rtErr = abErr / Max(fabs(org_vist[i]), smallLimit);

        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("%d term error in TestGPUTurbViscosity, vist on host: %.30e\t  on "
                   "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                   i, org_vist[i], h_vist[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("%d term error in TestGPUTurbViscosity, vist on host: %.30e\t  on "
                   "device:%.30e\t, abErr= %.30e, rtErr= %.30e\n",
                   i, org_vist[i], h_vist[i], abErr, rtErr);
            exit(EXIT_FAILURE);
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUTurbViscosity for vist is successful!\n");
#endif
    delete[] h_vist;
}

void TestGPUStoreRhsByResidual(RFloat **rhs, const int nl, const int nTotal)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t rhsSize = nl * nTotal * sizeof(RFloat);

    RFloat **g_rhs_ns = NewPointer2<RFloat>(nl, nTotal);
    HANDLE_API_ERR(cudaMemcpy(g_rhs_ns[0], d_rhs_ns, rhsSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; ++i)
    {
        for (int m = 0; m < nl; ++m)
        {
            abErr           = fabs(g_rhs_ns[m][i] - rhs[m][i]);
            rtErr           = abErr / Max(fabs(rhs[m][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("%d, %d term error in GPUStoreRhsByResidual , rhs on host: %e\t  "
                       "on device:%e\t, abErr= %e, rtErr= %e\n",
                       m, i, rhs[m][i], g_rhs_ns[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("%d, %d term error in GPUStoreRhsByResidual , rhs on host: %e\t  "
                       "on device:%e\t, abErr= %e, rtErr= %e\n",
                       m, i, rhs[m][i], g_rhs_ns[m][i], abErr, rtErr);
                exit(EXIT_FAILURE);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    printf("TestGPUStoreRhsByResidual is successful!\n");
#endif
    DelPointer2(g_rhs_ns);
}
