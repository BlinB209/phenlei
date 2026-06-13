#include "TemporaryOperationsPart2.h"
#include "BasicDeviceVariables.h"
#include "GPUKernelTest.h"

using namespace GPUMemory;
using namespace GPUFlowVariables;
namespace TemporaryOperations
{
    void GPUGradientVarsCopy(string q_name, int nTotal, int nTVar, RFloat **q2d, RFloat **dqdx2d, RFloat **dqdy2d,
                             RFloat **dqdz2d)
    {
        if (q_name == "Empty")
        {
            //!cout<<"no GPU field is appointed, Error"<<endl;
            //!cout<<"q_name=Empty"<<endl;
        }
        if (q_name == "d_q_proxy")
        {
            //!cout<<"q_name=d_q_proxy"<<endl;
            RFloat *q_proxy    = q2d[0];
            RFloat *dqdx_proxy = dqdx2d[0];
            RFloat *dqdy_proxy = dqdy2d[0];
            RFloat *dqdz_proxy = dqdz2d[0];
            size_t  sizeQFull  = nTotal * nTVar * sizeof(RFloat);

            //!comment the copy for q because the former copy in NSSolver::UpdateResiduals
            //!HANDLE_API_ERR(cudaMemcpy(d_q_ns, q_proxy, sizeQFull, cudaMemcpyHostToDevice));
            /*
            //!the copy of dqdx_proxy into d_dqdx_proxy is no use, becasue d_dqdx_proxy can be caculated from d_q_ns
                            HANDLE_API_ERR(cudaMemcpy(d_dqdx_proxy, dqdx_proxy, sizeQFull, cudaMemcpyHostToDevice));
                            HANDLE_API_ERR(cudaMemcpy(d_dqdy_proxy, dqdy_proxy, sizeQFull, cudaMemcpyHostToDevice));
                            HANDLE_API_ERR(cudaMemcpy(d_dqdz_proxy, dqdz_proxy, sizeQFull, cudaMemcpyHostToDevice));
            */
        }
        else if (q_name == "d_t_proxy")
        {
            RFloat *q_proxy    = q2d[0];
            RFloat *dqdx_proxy = dqdx2d[0];
            RFloat *dqdy_proxy = dqdy2d[0];
            RFloat *dqdz_proxy = dqdz2d[0];
            size_t  sizeQFull  = nTotal * nTVar * sizeof(RFloat);
            //! comment the copy for t because the initialization
            //!has been put onto device in CallGPUCompGamaAndTField_S kernel
            //!HANDLE_API_ERR(cudaMemcpy(d_t_proxy, q_proxy, sizeQFull, cudaMemcpyHostToDevice));
            /*
                HANDLE_API_ERR(cudaMemcpy(d_dtdx_proxy, dqdx_proxy, sizeQFull, cudaMemcpyHostToDevice));
                HANDLE_API_ERR(cudaMemcpy(d_dtdy_proxy, dqdy_proxy, sizeQFull, cudaMemcpyHostToDevice));
                HANDLE_API_ERR(cudaMemcpy(d_dtdz_proxy, dqdz_proxy, sizeQFull, cudaMemcpyHostToDevice));
*/
        }
        else if (q_name == "d_q_turb_proxy")
        {
            RFloat *q_proxy    = q2d[0];
            RFloat *dqdx_proxy = dqdx2d[0];
            RFloat *dqdy_proxy = dqdy2d[0];
            RFloat *dqdz_proxy = dqdz2d[0];
            size_t  sizeQFull  = nTotal * nTVar * sizeof(RFloat);
            //!It should be cancelled if all of operation about d_q_turb_proxy is on the device
            //!HANDLE_API_ERR(cudaMemcpy(d_q_turb_proxy, q_proxy, sizeQFull, cudaMemcpyHostToDevice));
            /*
                HANDLE_API_ERR(cudaMemcpy(d_dq_turbdx_proxy, dqdx_proxy, sizeQFull, cudaMemcpyHostToDevice));
                HANDLE_API_ERR(cudaMemcpy(d_dq_turbdy_proxy, dqdy_proxy, sizeQFull, cudaMemcpyHostToDevice));
                HANDLE_API_ERR(cudaMemcpy(d_dq_turbdz_proxy, dqdz_proxy, sizeQFull, cudaMemcpyHostToDevice));
*/
        }
        else if (q_name == "d_velocity_proxy")
        {
            RFloat *q_proxy    = q2d[0];
            RFloat *dqdx_proxy = dqdx2d[0];
            RFloat *dqdy_proxy = dqdy2d[0];
            RFloat *dqdz_proxy = dqdz2d[0];
            size_t  sizeQFull  = nTotal * nTVar * sizeof(RFloat);

            //!HANDLE_API_ERR(cudaMemcpy(d_vel_proxy, q_proxy, sizeQFull, cudaMemcpyHostToDevice));
            d_vel_proxy = d_q_ns + IDX::IU * nTotal; //! 1 here means
             /*
                HANDLE_API_ERR(cudaMemcpy(d_dveldx_proxy, dqdx_proxy, sizeQFull, cudaMemcpyHostToDevice));
                HANDLE_API_ERR(cudaMemcpy(d_dveldy_proxy, dqdy_proxy, sizeQFull, cudaMemcpyHostToDevice));
                HANDLE_API_ERR(cudaMemcpy(d_dveldz_proxy, dqdz_proxy, sizeQFull, cudaMemcpyHostToDevice));
             */
        }
    }

    void GPUViscousCoefCopy(const RFloat *h_visl, const RFloat *h_vist, const int nTotal)
    {
        size_t sizeVisCoef = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_visl, h_visl, sizeVisCoef, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_vist, h_vist, sizeVisCoef, cudaMemcpyHostToDevice));
    }

    void GPUResCopy(RFloat **res, const int nl, const int nTotal)
    {
        size_t  resSize = nl * nTotal * sizeof(RFloat);
        RFloat *h_res   = res[0];
        HANDLE_API_ERR(cudaMemcpy(d_res_ns, h_res, resSize, cudaMemcpyHostToDevice));
    }

    void GPUResCopyBack(RFloat **res, const int nl, const int nTotal)
    {
        size_t  resSize = nl * nTotal * sizeof(RFloat);
        RFloat *h_res   = res[0];
        HANDLE_API_ERR(cudaMemcpy(h_res, d_res_ns, resSize, cudaMemcpyDeviceToHost));
    }

    void GPUGradientCopyBack(const string q_name, const int nTotal, const int nTVar, RFloat **dqdx2d, RFloat **dqdy2d,
                             RFloat **dqdz2d)
    {
        RFloat *d_dqdx;
        RFloat *d_dqdy;
        RFloat *d_dqdz;
        if (q_name == "d_q_proxy")
        {
            d_dqdx = d_dqdx_proxy;
            d_dqdy = d_dqdy_proxy;
            d_dqdz = d_dqdz_proxy;
        }
        else if (q_name == "d_q_turb_proxy")
        {
            d_dqdx = d_dq_turbdx_proxy;
            d_dqdy = d_dq_turbdy_proxy;
            d_dqdz = d_dq_turbdz_proxy;
        }
        else if (q_name == "d_t_proxy")
        {
            d_dqdx = d_dtdx_proxy;
            d_dqdy = d_dtdy_proxy;
            d_dqdz = d_dtdz_proxy;
        }
        else if (q_name == "d_velocity_proxy")
        {
            d_dqdx = d_dveldx_proxy;
            d_dqdy = d_dveldy_proxy;
            d_dqdz = d_dveldz_proxy;
        }
        RFloat *h_dqdx = dqdx2d[0];
        RFloat *h_dqdy = dqdy2d[0];
        RFloat *h_dqdz = dqdz2d[0];
        size_t  dqSize = nTVar * nTotal * sizeof(RFloat);
        cudaMemcpy(h_dqdx, d_dqdx, dqSize, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_dqdy, d_dqdy, dqSize, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_dqdz, d_dqdz, dqSize, cudaMemcpyDeviceToHost);
    }

    void GPUStoreBoundGradCopyBack(const int nBoundFace, const int nTVar, RFloat **bdqx, RFloat **bdqy, RFloat **bdqz)
    {
        RFloat *h_bdqx = bdqx[0];
        RFloat *h_bdqy = bdqy[0];
        RFloat *h_bdqz = bdqz[0];
        size_t  bGSize = nBoundFace * nTVar * sizeof(RFloat);
        cudaMemcpy(h_bdqx, d_bdqx, bGSize, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_bdqy, d_bdqy, bGSize, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_bdqz, d_bdqz, bGSize, cudaMemcpyDeviceToHost);
    }

    void GPUQ_turbCopy(RFloat **q_turb, const int nTotal, const int n_turb)
    {
        size_t  sizeQturb = nTotal * n_turb * sizeof(RFloat);
        RFloat *q_turb0   = q_turb[0];
        HANDLE_API_ERR(cudaMemcpy(d_q_turb_proxy, q_turb0, sizeQturb, cudaMemcpyHostToDevice));
    }

    void GPUDQ_turbCopy(RFloat **dq_turb, const int nTotal, const int n_turb)
    {
        size_t  sizeQturb = nTotal * n_turb * sizeof(RFloat);
        RFloat *dq_turb0  = dq_turb[0];
        HANDLE_API_ERR(cudaMemcpy(d_dq_turb, dq_turb0, sizeQturb, cudaMemcpyHostToDevice));
    }

    void GPUDQ_turbCopyBack(RFloat **dq_turb, const int nTotal, const int n_turb)
    {
        size_t  sizeQturb = nTotal * n_turb * sizeof(RFloat);
        RFloat *dq_turb0  = dq_turb[0];
        //!HANDLE_API_ERR(cudaMemcpy(d_dq_turb, dq_turb0, sizeQturb, cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(dq_turb0, d_dq_turb, sizeQturb, cudaMemcpyDeviceToHost));
    }

    void GPUFluxTurbCopy(RFloat **flux, RFloat **flux_sub, const int n_turb)
    {
        size_t  sizeFluxTurb = d_SEG_LEN * n_turb * sizeof(RFloat);
        RFloat *flux0        = flux[0];
        RFloat *flux_sub0    = flux_sub[0];
        cudaMemcpy(d_flux_turb, flux0, sizeFluxTurb, cudaMemcpyHostToDevice);
        cudaMemcpy(d_flux_sub_turb, flux_sub0, sizeFluxTurb, cudaMemcpyHostToDevice);
    }

    void GPUResTurbCopy(RFloat **res_turb, const int nTotal, const int n_turb)
    {
        size_t  resTurbSize = n_turb * nTotal * sizeof(RFloat);
        RFloat *h_res       = res_turb[0];
        HANDLE_API_ERR(cudaMemcpy(d_res_turb, h_res, resTurbSize, cudaMemcpyHostToDevice));
    }

    void GPUResTurbCopyBack(RFloat **res_turb, const int n_turb, const int nTotal)
    {
        size_t  resTurbSize = n_turb * nTotal * sizeof(RFloat);
        RFloat *h_res       = res_turb[0];
        /*
        RFloat flagCompare = 1.0e-8;
        RFloat error = 0.0;
        
        for (int i = 0; i < nTotal; i++){
            err = fabs(h_res[i] - );
            if (h_res[i] -)
    i    }
        */

        HANDLE_API_ERR(cudaMemcpy(h_res, d_res_turb, resTurbSize, cudaMemcpyDeviceToHost));
    }

    void GPUSpecTurbCopy(const int nTotalCell, const int nBoundFace, const int n_turb, RFloat **spec_turb)
    {
        int     nTotal       = nTotalCell + nBoundFace;
        size_t  sizeSpecTurb = n_turb * nTotal * sizeof(RFloat);
        RFloat *h_spec_turb  = spec_turb[0];
        HANDLE_API_ERR(cudaMemcpy(d_spec_turb, h_spec_turb, sizeSpecTurb, cudaMemcpyHostToDevice));
    }

    void GPUDtCopy(const int nTotalCell, const int nBoundFace, const RFloat *dt)
    {
        int    nTotal = nTotalCell + nBoundFace;
        size_t sizeDt = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(d_dt, dt, sizeDt, cudaMemcpyHostToDevice));
    }

    void GPUSpecTurbCopyBack(RFloat **spec_turb, const int n_turb, const int nTotal)
    {
        size_t  sizeSpecTurb = nTotal * n_turb * sizeof(RFloat);
        RFloat *h_spec_turb  = spec_turb[0];
        HANDLE_API_ERR(cudaMemcpy(h_spec_turb, d_spec_turb, sizeSpecTurb, cudaMemcpyDeviceToHost));
    }

    void GPUMatTurbCopyBack(RFloat **mat_turbl, RFloat **mat_turbr, const int n_turb, const int nTotalFace)
    {
        size_t  sizeMat     = nTotalFace * n_turb * sizeof(RFloat);
        RFloat *h_mat_turbl = mat_turbl[0];
        RFloat *h_mat_turbr = mat_turbr[0];
        cudaMemcpy(h_mat_turbl, d_mat_turbl, sizeMat, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_mat_turbr, d_mat_turbr, sizeMat, cudaMemcpyDeviceToHost);
    }
    void GPUQ_turbCopyBack(RFloat **q_turb, const int nTotal, const int n_turb)
    {
        size_t  sizeQturb = nTotal * n_turb * sizeof(RFloat);
        RFloat *q_turb0   = q_turb[0];
        HANDLE_API_ERR(cudaMemcpy(q_turb0, d_q_turb_proxy, sizeQturb, cudaMemcpyDeviceToHost));
    }
    void GPUTurbRhsCopy(RFloat **rhs, const int nTotal, const int n_turb)
    {
        size_t  sizeRhs = nTotal * n_turb * sizeof(RFloat);
        RFloat *h_rhs   = rhs[0];
        HANDLE_API_ERR(cudaMemcpy(d_rhs_turb, h_rhs, sizeRhs, cudaMemcpyHostToDevice));
    }

    void GPUQProxyTempTurbCopyBack(RFloat **field1, const int nTotal, const int neqn)
    {
        using namespace GPUControlVariables;
        size_t  sizeQProxyTemp = nTotal * neqn * sizeof(RFloat);
        RFloat *field10        = field1[0];
        if (d_FillField_turb == "q_proxy_tempToq_proxy")
        {
            //!HANDLE_API_ERR(cudaMemcpy(field10, d_q_proxy_turb, sizeQProxyTemp, cudaMemcpyDeviceToHost));
            HANDLE_API_ERR(cudaMemcpy(field10, d_q_turb_proxy, sizeQProxyTemp, cudaMemcpyDeviceToHost));
        }
        else
            return;
    }

    void GPUSpecCopyBack(RFloat *spec, const int nTotal)
    {
        size_t sizeSpec = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(spec, d_spec_ns, sizeSpec, cudaMemcpyDeviceToHost));
    }
} //! namespace TemporaryOperations
