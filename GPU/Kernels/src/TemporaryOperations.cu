#include "BasicDeviceVariables.h"
#include "Limiter.h"
#include "Pointer.h"
#include "TemporaryOperations.h"
#include "TemporaryOperationsPart2.h"

using namespace PHSPACE;
using namespace PHMPI;
using namespace std;
using namespace GPUMemory;
using namespace GPUFlowVariables;
using namespace GPUTestSpace;
using namespace GPUControlVariables;
namespace TemporaryOperations
{
    void GPUQCopy(Grid *grid_in, const int neqn, const int nTotal)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        RFloat **q = reinterpret_cast<RFloat **>(grid->GetDataPtr("q"));

        //! q ql qr copy to devcie
        GPUQ_NSCopy(q, neqn, nTotal);
//! test the copy
#ifdef CUDAUNITTEST
        TestGPUQCopy(q, d_q_ns, neqn, nTotal);
#endif
    }

    void GPUTurbUnsteadyCopy(Grid *grid_in)
    {
        UnstructGrid *grid   = UnstructGridCast(grid_in);
        int           n_turb = 1;
        GlobalDataBase::GetData("n_turb", &n_turb, PHINT, 1);
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        size_t sizeUnsteady = nTotal * n_turb * sizeof(RFloat);

        RFloat **qn1    = reinterpret_cast<RFloat **>(grid->GetDataPtr("q_turb_unsteady_n1"));
        RFloat **qn2    = reinterpret_cast<RFloat **>(grid->GetDataPtr("q_turb_unsteady_n2"));
        RFloat **resn1  = reinterpret_cast<RFloat **>(grid->GetDataPtr("res_turb_unsteady_n1"));
        RFloat **resn2  = reinterpret_cast<RFloat **>(grid->GetDataPtr("res_turb_unsteady_n2"));
        RFloat **restmp = reinterpret_cast<RFloat **>(grid->GetDataPtr("res_turb_unsteady_tmp"));
        HANDLE_API_ERR(cudaMemcpy(d_q_turb_unsteady_n1, qn1[0], sizeUnsteady, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_q_turb_unsteady_n2, qn2[0], sizeUnsteady, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_res_turb_unsteady_n1, resn1[0], sizeUnsteady, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_res_turb_unsteady_n2, resn2[0], sizeUnsteady, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_res_turb_unsteady_tmp, restmp[0], sizeUnsteady, cudaMemcpyHostToDevice));
    }

    void GPUNSUnsteadyCopy(Grid *grid_in)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nl         = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem      = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn       = nl + nchem;
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        size_t sizeUnsteady = nTotal * neqn * sizeof(RFloat);

        RFloat **qn1    = reinterpret_cast<RFloat **>(grid->GetDataPtr("q_unsteady_n1"));
        RFloat **qn2    = reinterpret_cast<RFloat **>(grid->GetDataPtr("q_unsteady_n2"));
        RFloat **resn1  = reinterpret_cast<RFloat **>(grid->GetDataPtr("res_unsteady_n1"));
        RFloat **resn2  = reinterpret_cast<RFloat **>(grid->GetDataPtr("res_unsteady_n2"));
        RFloat **restmp = reinterpret_cast<RFloat **>(grid->GetDataPtr("res_unsteady_tmp"));
        HANDLE_API_ERR(cudaMemcpy(d_q_ns_unsteady_n1, qn1[0], sizeUnsteady, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_q_ns_unsteady_n2, qn2[0], sizeUnsteady, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_res_ns_unsteady_n1, resn1[0], sizeUnsteady, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_res_ns_unsteady_n2, resn2[0], sizeUnsteady, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_res_ns_unsteady_tmp, restmp[0], sizeUnsteady, cudaMemcpyHostToDevice));
    }

    void GPUQCopy(Grid *grid_in)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nl         = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem      = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn       = nl + nchem;
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        RFloat **q = reinterpret_cast<RFloat **>(grid->GetDataPtr("q"));

        //! q ql qr copy to devcie
        GPUQ_NSCopy(q, neqn, nTotal);
//! test the copy
#ifdef CUDAUNITTESTTEST
        TestGPUQCopy(q, d_q_ns, neqn, nTotal);
#endif
    }

    void GPUQ_NSCopy(RFloat **q, const int m, const int n)
    {
        size_t qSize = m * n * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_q_ns, q[0], qSize, cudaMemcpyHostToDevice));
    }

    void GPUQl_NSCopy(RFloat **ql, const int m, const int n)
    {
        size_t qlSize = m * n * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_ql_ns, ql[0], qlSize, cudaMemcpyHostToDevice));
    }
    void GPUQr_NSCopy(RFloat **qr, const int m, const int n)
    {
        size_t qrSize = m * n * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_qr_ns, qr[0], qrSize, cudaMemcpyHostToDevice));
    }
    //! q copy back to host
    void GPUQ_NSCopyBack(RFloat **q, const int m, const int n)
    {
        size_t qSize = m * n * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(q[0], d_q_ns, qSize, cudaMemcpyDeviceToHost));
    }
    void CallGPUFlowVariablesCopyBack(Grid *grid_in)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;
        int nl         = 5;
        GlobalDataBase::GetData("nl", &nl, PHINT, 1);

        int nchem = 0;
        GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
        int      neqn = nl + nchem;
        RFloat **q    = reinterpret_cast<RFloat **>(grid->GetDataPtr("q"));
        GPUQ_NSCopyBack(q, neqn, nTotal);
        RFloat *visl = reinterpret_cast<RFloat *>(grid->GetDataPtr("visl"));
        GPUVislCopyBack(visl, nTotal);
        //!copy back of bdqx, bdqy, and bdqz should be also added as well
        RDouble    **gradPrimtiveVarX = reinterpret_cast<RDouble **>(grid->GetDataPtr("gradPrimtiveVarX"));
        RDouble    **gradPrimtiveVarY = reinterpret_cast<RDouble **>(grid->GetDataPtr("gradPrimtiveVarY"));
        RDouble    **gradPrimtiveVarZ = reinterpret_cast<RDouble **>(grid->GetDataPtr("gradPrimtiveVarZ"));
        const string q_name           = "d_q_proxy";
        GPUGradientCopyBack(q_name, nTotal, neqn, gradPrimtiveVarX, gradPrimtiveVarY, gradPrimtiveVarZ);
    }
    //! copy ql qr back to host
    RFloat **GPUQl_NSCopyBack(RFloat *d_ql_ns, const int m, const int n)
    {
        RFloat **g_ql_ns = NewPointer2<RFloat>(m, n);
        size_t   qlSize  = m * n * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(g_ql_ns[0], d_ql_ns, qlSize, cudaMemcpyDeviceToHost));
        return g_ql_ns;
    }
    RFloat **GPUQr_NSCopyBack(RFloat *d_qr_ns, const int m, const int n)
    {
        RFloat **g_qr_ns = NewPointer2<RFloat>(m, n);
        size_t   qrSize  = m * n * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(g_qr_ns[0], d_qr_ns, qrSize, cudaMemcpyDeviceToHost));
        return g_qr_ns;
    }
    void GPUGamaCopy(Grid *grid_in, const int nTotal)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        RFloat *gama = reinterpret_cast<RFloat *>(grid->GetDataPtr("gama"));

        size_t gamaSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_gama_ns, gama, gamaSize, cudaMemcpyHostToDevice));
    }
    RFloat *GPUGamaL_NSCopyBack(RFloat *d_gamaL_ns, const int len_gamaL)
    {
        RFloat *g_gamaL_ns = new RFloat[len_gamaL];
        size_t  gamaLSize  = len_gamaL * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(g_gamaL_ns, d_gamaL_ns, gamaLSize, cudaMemcpyDeviceToHost));
        return g_gamaL_ns;
    }
    RFloat *GPUGamaR_NSCopyBack(RFloat *d_gamaR_ns, const int len_gamaR)
    {
        RFloat *g_gamaR_ns = new RFloat[len_gamaR];
        size_t  gamaRSize  = len_gamaR * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(g_gamaR_ns, d_gamaR_ns, gamaRSize, cudaMemcpyDeviceToHost));
        return g_gamaR_ns;
    }
    //! copy limiter onto device
    void GPULimitCopy(Limiter *limiter, const int neqn, const int nTotal)
    {
        //! limit on device
        RFloat *limit     = limiter->GetLimiter(0);
        size_t  limitSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_limit, limit, limitSize, cudaMemcpyHostToDevice));
#ifdef KERNELLAUNCHTEST
        TestGPULimitCopy(limit, d_limit, nTotal);
#endif
        /*
  RFloat **LIMIT = new RFloat * [neqn];
  for ( int m = 0; m < neqn; ++ m )
  {
          LIMIT[m]        =  limiter->GetLimiter(m);
  }
  size_t LIMITSize= neqn * nTotal *sizeof(RFloat);
  HANDLE_API_ERR(cudaMemcpy(d_LIMIT, LIMIT[0], LIMITSize,
cudaMemcpyHostToDevice)); #ifdef KERNELLAUNCHTEST
  //!TestGPULIMITCopy(LIMIT, d_LIMIT, neqn, nTotal);//!The bug will induce crack!!
#endif
  delete [] LIMIT;
  */
    }
    //! copy LIMIT onto device
    void GPULIMITCopy(RFloat **LIMIT, const int neqn, const int nTotal)
    {
        size_t LIMITSize = neqn * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_LIMIT, LIMIT[0], LIMITSize, cudaMemcpyHostToDevice));
#ifdef KERNELLAUNCHTEST
        TestGPULIMITCopy(LIMIT, d_LIMIT, neqn, nTotal);
#endif
    }
    //! copy dqdx ... onto device
    void GPUDqdxDqDyDqDzCopy(Gradient *gradient, const int neqn, const int nTotal)
    {
        RFloat **dqdx = gradient->GetGradX()->GetField_UNS();
        RFloat **dqdy = gradient->GetGradY()->GetField_UNS();
        RFloat **dqdz = gradient->GetGradZ()->GetField_UNS();

        size_t dSize = neqn * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_dqdx_proxy, dqdx[0], dSize, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_dqdy_proxy, dqdy[0], dSize, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_dqdz_proxy, dqdz[0], dSize, cudaMemcpyHostToDevice));

#ifdef KERNELLAUNCHTEST
        TestGPUDqdxDqDyDqDzCopy(dqdx, dqdy, dqdz, d_dqdx_proxy, d_dqdy_proxy, d_dqdz_proxy, neqn, nTotal);
#endif
    }
    //! copy xtn .. onto device
    void GPUXtnYtnZtnCopy(Grid *grid_in, const int nTotalFace)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        RDouble *xtn = grid->GetFaceVelocityX();
        RDouble *ytn = grid->GetFaceVelocityY();
        RDouble *ztn = grid->GetFaceVelocityZ();

        size_t tnSize = nTotalFace * sizeof(RDouble);
        HANDLE_API_ERR(cudaMemcpy(d_xtn, xtn, tnSize, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_ytn, ytn, tnSize, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_ztn, ztn, tnSize, cudaMemcpyHostToDevice));
#ifdef KERNELLAUNCHTEST
        //! test the copy of xtn ytn ztn
        TestGPUXtnYtnZtnCopy(xtn, ytn, ztn, d_xtn, d_ytn, d_ztn, nTotalFace);
#endif
    }
    //! copy vgn ont to device
    void GPUVgnCopy(Grid *grid_in, const int nTotalFace)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        RDouble *vgn     = grid->GetFaceNormalVelocity();
        size_t   vgnSize = nTotalFace * sizeof(RDouble);
        HANDLE_API_ERR(cudaMemcpy(d_vgn, vgn, vgnSize, cudaMemcpyHostToDevice));
#ifdef KERNELLAUNCHTEST
        //! test the copy of vgn
        TestGPUVgnCopy(vgn, d_vgn, nTotalFace);
#endif
    }
    //! copy res onto device
    void GPUResCopy(Grid *grid_in)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        //! Get Residual
        RFloat **res     = reinterpret_cast<RFloat **>(grid->GetDataPtr("res"));
        size_t   resSize = neqn * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_res_ns, res[0], resSize, cudaMemcpyHostToDevice));
#ifdef KERNELLAUNCHTEST
        //! test the copy of res
        TestGPUResCopy(res, d_res_ns, neqn, nTotal);
#endif
    }
    //! copy epsCell onto device
    void GPUepsCellCopy(RFloat *epsCell, const int nTotalCell)
    {
        size_t epsCellSize = nTotalCell * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_epsCell, epsCell, epsCellSize, cudaMemcpyHostToDevice));
#ifdef KERNELLAUNCHTEST
        //! test the copy of q_tmp
        TestGPUepsCellCopy(epsCell, d_epsCell, nTotalCell);
#endif
    }
    void GPUQlQrCopyBack(FaceProxy *face_proxy, const int SEG_LEN)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;

        RFloat **ql = face_proxy->GetQL();
        RFloat **qr = face_proxy->GetQR();

        size_t qSize = neqn * SEG_LEN * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(ql[0], d_ql_ns, qSize, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(qr[0], d_qr_ns, qSize, cudaMemcpyDeviceToHost));
    }
    void GPUGamaLRCopyBack(FaceProxy *face_proxy, const int SEG_LEN)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        RFloat *gamal = face_proxy->GetGamaL();
        RFloat *gamar = face_proxy->GetGamaR();

        size_t gamaSize = SEG_LEN * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(gamal, d_gamaL_ns, gamaSize, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(gamar, d_gamaR_ns, gamaSize, cudaMemcpyDeviceToHost));
    }
    void GPUFluxCopyBack(FaceProxy *face_proxy, const int SEG_LEN)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;

        RFloat **flux = face_proxy->GetFlux();

        size_t fluxSize = neqn * SEG_LEN * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(flux[0], d_flux, fluxSize, cudaMemcpyDeviceToHost));
    }
    //! copy res back to cpu
    void GPUResCopyBack(Grid *grid_in, const int nTotal)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;
        //! Get Residual
        RFloat **res     = reinterpret_cast<RFloat **>(grid->GetDataPtr("res"));
        size_t   resSize = neqn * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(res[0], d_res_ns, resSize, cudaMemcpyDeviceToHost));
    }
    //! copy d_epsCell back to device
    void GPUepsCellCopyBack(RFloat *epsCell, RFloat *d_epsCell, const int nTotalCell)
    {
        size_t epsSize = nTotalCell * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(epsCell, d_epsCell, epsSize, cudaMemcpyDeviceToHost));
    }
    //! copy dmax dmin onto device
    void GPUDmaxDminCopy(RFloat *dmax, RFloat *dmin, const int nTotal)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        size_t dSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_dmin, dmin, dSize, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_dmax, dmax, dSize, cudaMemcpyHostToDevice));
#ifdef KERNELLAUNCHTEST
        //! test the copy of dmax, dmin
        TestGPUDminDmaxAllocCopy(dmax, dmin, d_dmax, d_dmin, nTotal);
#endif
    }
    void GPUQTurbCopy(RFloat **q, const int n_turb)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int    nTotal = d_nTotal;
        size_t qSize  = n_turb * nTotal * sizeof(RFloat);

        HANDLE_API_ERR(cudaMemcpy(d_q_turb_proxy, q[0], qSize, cudaMemcpyHostToDevice));
#ifdef KERNELLAUNCHTEST
        cout << "Test GPUQTurbCopy" << endl;
        TestGPUQCopy(q, d_q_turb_proxy, n_turb, nTotal);
#endif
    }
    void GPUQlQrTurbCopyBack(FaceProxy *face_proxy, const int n_turb, const int SEG_LEN)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        RFloat **ql = face_proxy->GetQL();
        RFloat **qr = face_proxy->GetQR();

        size_t qSize = n_turb * SEG_LEN * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(ql[0], d_ql_turb, qSize, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(qr[0], d_qr_turb, qSize, cudaMemcpyDeviceToHost));
    }
    void GPUQCopy(RFloat **q, RFloat *d_q, const int neqn, const int nTotal)
    {
        size_t qSize = neqn * nTotal * sizeof(RFloat);

        HANDLE_API_ERR(cudaMemcpy(d_q, q[0], qSize, cudaMemcpyHostToDevice));
#ifdef KERNELLAUNCHTEST
        TestGPUQCopy(q, d_q, neqn, nTotal);
#endif
    }
    void GPUFluxTurbCopyBack(RFloat **flux, RFloat *d_flux, const int neqn, const int SEG_LEN)
    {
        size_t fluxSize = neqn * SEG_LEN * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(flux[0], d_flux, fluxSize, cudaMemcpyDeviceToHost));
    }

    void GPUFluxSubTurbCopy(RFloat **flux_sub, RFloat *d_flux_sub, const int neqn, const int SEG_LEN)
    {
        size_t fluxSize = neqn * SEG_LEN * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_flux_sub, flux_sub[0], fluxSize, cudaMemcpyHostToDevice));
    }
    void GPUVislCopy(RFloat *visl, RFloat *d_visl, const int nTotal)
    {
        size_t vSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_visl, visl, vSize, cudaMemcpyHostToDevice));
    }
    void GPUVislCopy(RFloat *visl, const int nTotal)
    {
        size_t vSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_visl, visl, vSize, cudaMemcpyHostToDevice));
    }
    void GPUVislCopy(RFloat **res_turb, RFloat *d_res_turb, const int neqn, const int nTotal)
    {
        size_t vSize = neqn * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_res_turb, res_turb[0], vSize, cudaMemcpyHostToDevice));
    }
    void GPUResTurbCopyBack(RFloat **res_turb, RFloat *d_res_turb, const int neqn, const int nTotal)
    {
        size_t vSize = neqn * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(res_turb[0], d_res_turb, vSize, cudaMemcpyDeviceToHost));
    }
    void GPUGamaCopyBack(Grid *grid_in)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        RFloat *gama = reinterpret_cast<RFloat *>(grid->GetDataPtr("gama"));

        size_t gamaSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(gama, d_gama_ns, gamaSize, cudaMemcpyDeviceToHost));
    }

    void GPUTCopyBack(Grid *grid_in)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        int      nTotalCell = grid->GetNTotalCell();
        int      nBoundFace = grid->GetNBoundFace();
        int      nTotal     = nTotalCell + nBoundFace;
        int      ntmodel    = 1;
        RFloat **t          = reinterpret_cast<RFloat **>(grid->GetDataPtr("t"));

        size_t tSize = ntmodel * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(t[0], d_t_proxy, tSize, cudaMemcpyDeviceToHost));
    }
    //!void GPUPrim_infCopy(RFloat *prim_inf, const int neqn)
    void GPUPrim_infCopy()
    {
        int nl;
        int nchem;

        GlobalDataBase::GetData("nl", &nl, PHINT, 1);
        GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
        RFloat *prim_inf = reinterpret_cast<RFloat *>(GlobalDataBase::GetDataPtr("prim_inf"));
        size_t  sSize    = (nl + nchem) * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_prim_inf, prim_inf, sSize, cudaMemcpyHostToDevice));
    }
    //! q_proxy copy back to host
    void GPUQ_proxyCopyBack(RFloat **q, const int m, const int n)
    {
        size_t qSize = m * n * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(q[0], d_q_proxy, qSize, cudaMemcpyDeviceToHost));
    }
    void GPUQ_proxy_tmpCopyBack(RFloat **q, const int m, const int n)
    {
        size_t qSize = m * n * sizeof(RFloat);
        if (IsFillFieldBeforeStage)
        {
            HANDLE_API_ERR(cudaMemcpy(q[0], d_q_proxy_tmp, qSize, cudaMemcpyDeviceToHost));
        }
        else
        {
            HANDLE_API_ERR(cudaMemcpy(q[0], d_q_proxy, qSize, cudaMemcpyDeviceToHost));
        }
    }
    void GPUrhsCopy(RFloat **rhs, const int neqn, const int nTotal)
    {
        size_t sSize = neqn * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_rhs_ns, rhs[0], sSize, cudaMemcpyHostToDevice));
    }
    void GPUrhsCopyBack(RFloat **rhs, const int neqn, const int nTotal)
    {
        size_t sSize = neqn * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(rhs[0], d_rhs_ns, sSize, cudaMemcpyDeviceToHost));
    }
    void GPUDtCopy(RFloat *dt, const int nTotal)
    {
        size_t sSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_dt, dt, sSize, cudaMemcpyHostToDevice));
    }
    void GPUResCopyBack(Grid *grid_in)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        //! Get Residual
        RFloat **res     = reinterpret_cast<RFloat **>(grid->GetDataPtr("res"));
        size_t   resSize = neqn * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(res[0], d_res_ns, resSize, cudaMemcpyDeviceToHost));
    }

    void GPUResTurbCopyBack(Grid *grid_in)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nTurbulenceEquation = GlobalDataBase::GetIntParaFromDB("n_turb");

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        //! Get Residual
        RFloat **resTurb = reinterpret_cast<RFloat **>(grid->GetDataPtr("res_turb"));
        size_t   resSize = nTurbulenceEquation * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(resTurb[0], d_res_turb, resSize, cudaMemcpyDeviceToHost));
    }

    void GPUInvSpectrumRadiusCopyBack(RFloat *invSpectrumRadius, const int nTotalCell)
    {
        size_t iSize = nTotalCell * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(invSpectrumRadius, d_invSpectrumRadius, iSize, cudaMemcpyDeviceToHost));
    }
    void GPUInvSpectrumRadiusCopy(RFloat *invSpectrumRadius, const int nTotalCell)
    {
        size_t iSize = nTotalCell * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_invSpectrumRadius, invSpectrumRadius, iSize, cudaMemcpyHostToDevice));
    }
    void GPUVisSpectrumRadiusCopyBack(RFloat *visSpectrumRadius, const int nTotalCell)
    {
        size_t iSize = nTotalCell * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(visSpectrumRadius, d_visSpectrumRadius, iSize, cudaMemcpyDeviceToHost));
    }
    void GPUVisSpectrumRadiusCopy(RFloat *visSpectrumRadius, const int nTotalCell)
    {
        size_t iSize = nTotalCell * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_visSpectrumRadius, visSpectrumRadius, iSize, cudaMemcpyHostToDevice));
    }
    void GPUDtCopyBack(RFloat *dt, const int nTotal)
    {
        size_t sSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(dt, d_dt, sSize, cudaMemcpyDeviceToHost));
    }
    void GPUDtvCopyBack(RFloat *dtv, const int nTotalCell)
    {
        size_t sSize = nTotalCell * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(dtv, d_dtv, sSize, cudaMemcpyDeviceToHost));
    }
    void GPUMDtCopyBack(RFloat &minDt, RFloat &maxDt)
    {
        size_t sSize = 1 * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(&minDt, d_minDt, sSize, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(&maxDt, d_maxDt, sSize, cudaMemcpyDeviceToHost));
    }
    void GPUVislCopyBack(RFloat *visl, const int nTotal)
    {
        size_t sSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(visl, d_visl, sSize, cudaMemcpyDeviceToHost));
    }

    void GPUQ_TurbCopy(RFloat **q, const int m, const int n)
    {
        size_t qSize = m * n * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_q_turb_proxy, q[0], qSize, cudaMemcpyHostToDevice));
    }

    void GPUVistCopy(RFloat *vist, const int nTotal)
    {
        size_t vSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_vist, vist, vSize, cudaMemcpyHostToDevice));
    }
    void GPUVistCopyBack(RFloat *vist, const int nTotal)
    {
        size_t vSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(vist, d_vist, vSize, cudaMemcpyDeviceToHost));
    }
    void GPUMVistCopyBack(RFloat &vistmin, RFloat &vistmax)
    {
        size_t sSize = 1 * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(&vistmin, d_vistmin, sSize, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(&vistmax, d_vistmax, sSize, cudaMemcpyDeviceToHost));
    }

    void GPUDQNSCopy(RFloat **dq, const int nl, const int nTotal)
    {
        size_t qSize = nl * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_dq_ns, dq[0], qSize, cudaMemcpyHostToDevice));
    }
    void GPUDQNSCopyBack(RFloat **dq, const int nl, const int nTotal)
    {
        size_t qSize = nl * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(dq[0], d_dq_ns, qSize, cudaMemcpyDeviceToHost));
    }
} //! namespace TemporaryOperations
