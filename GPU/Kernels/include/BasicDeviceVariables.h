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
//! @file      BasicDeviceVariables.h
//! @brief     All variables in device.
//! @author    Zhang Yong, Huang Can, Zhang Xi, Bell.

#ifndef BasicDeviceVariables_H
#define BasicDeviceVariables_H
#include <iostream>
#include "Constants.h"
#include "Geo_Grid.h"
#include "Geo_UnstructGrid.h"
#include "GlobalDataBase.h"
#include "Precision.h"
#include "Region.h"
#include "Zone.h"
#include "cudaErrorHandle.h"
#ifdef FACEPROXYOPT
#include "FaceProxy.h"
#endif
using namespace PHSPACE;
namespace GPUMemory
{
    namespace GPUControlVariables
    {
        extern string d_gradient_field_proxy;
        extern int    d_gradient_var_index;
        extern int    d_m_count;
        extern string d_visflux_solver;
        extern string d_FillField_turb;
        extern int    IsInviscidFluxLoad;
        extern int    IsFillFieldBeforeStage;
        extern string d_FillField_ns;
        extern string d_Limiter_var;

        void SetGPUGradientControlVars(string &d_gradient_control_var, const string value);
        void SetGPUGradientVarIdx(int *d_gradient_var_index, int value);
        void SetGPUVisFluxSolverControlVars(string &d_gradient_control_var, const string value);
        void SetGPUTurbFillField(string &d_FillField_control_var, const string value);
        void SetGPUNSFillField(string &d_FillField_control_var, const string value);
        void SetGPULimiterControlVarsByValue(string &d_Limiter_control_var, const int value);
    } //! namespace GPUControlVariables

    namespace GPUFlowVariables
    {
        //! equation numbers
        extern int d_nl;
        extern int d_nchem;
        extern int d_neqn_ns;
        extern int d_neqn_turb;
        //! field proxy variable d_q_proxy
        extern RFloat *d_q_proxy;
        extern RFloat *d_q_proxy_tmp;
        //! field proxy dqdx, dqdy, dqdz variables d_dqdx_proxy,d_dqdy_proxy,d_dqdz_proxy
        extern RFloat *d_dqdx_proxy;
        extern RFloat *d_dqdy_proxy;
        extern RFloat *d_dqdz_proxy;
        //! store boundary gradient value of q
        extern RFloat *d_bdqx;
        extern RFloat *d_bdqy;
        extern RFloat *d_bdqz;
        //! field proxy variables for temporature t
        extern RFloat *d_t_proxy;
        extern RFloat *d_dtdx_proxy;
        extern RFloat *d_dtdy_proxy;
        extern RFloat *d_dtdz_proxy;
        //! field proxy varialbes for turbulence q_turb
        extern RFloat *d_q_turb_proxy;
        extern RFloat *d_dq_turbdx_proxy;
        extern RFloat *d_dq_turbdy_proxy;
        extern RFloat *d_dq_turbdz_proxy;
        //! field proxy variables for velocity velocity
        extern RFloat *d_vel_proxy;
        extern RFloat *d_dveldx_proxy;
        extern RFloat *d_dveldy_proxy;
        extern RFloat *d_dveldz_proxy;
        //! field variable d_q_old: d_q_proxy in the last iteration step
        extern RFloat *d_q_old;
        //! node tmp variable d_q_n for CompNodeVar
        extern RFloat *d_q_n_tmp;
        //! extern double *d_q_n_double;
        //! node tmp variable d_n_count for CompNodeVar
        extern int *d_n_count_tmp;
        //! stores the order of field variabe: it often varies in loop
        extern int d_nTVar_m;
        //! SEG_LEN On device
        extern int d_SEG_LEN;
        //! flow variables q
        extern RFloat *d_q_ns;
        extern RFloat *d_ql_ns;
        extern RFloat *d_qr_ns;
        //! flow variables gamma
        extern RFloat *d_gama_ns;
        extern RFloat *d_gamaL_ns;
        extern RFloat *d_gamaR_ns;
        //! limiter variables
        extern RFloat *d_limit;
        extern RFloat *d_LIMIT;
        //! vgn FaceNormalVelocity
        extern RDouble *d_vgn;
        //! xtn ytn ztn NodeVelocity
        extern RDouble *d_xtn;
        extern RDouble *d_ytn;
        extern RDouble *d_ztn;

        //! viscous coeff
        extern RFloat *d_visl;
        extern RFloat *d_vist;
        extern RFloat *d_minLocalTurbVis; //! local minimum turb viscous in one GPU
        extern RFloat *d_maxLocalTurbVis;
        //! face proxy varialbes
        extern RFloat  *d_flux;
        extern RDouble *d_deltl;
        extern RDouble *d_deltr;
        //! ns face value
        extern RFloat *d_prim;
        extern RFloat *d_tm;
        extern RFloat *d_kcp;
        extern RFloat *d_mul;
        extern RFloat *d_mut;
        extern RFloat *d_rho_ds;
        extern RFloat *d_hint_s;

        //! residuals
        extern RFloat *d_res_ns;
        //! tmp variable: epsCell
        extern RFloat *d_epsCell;
        //! dmin dmax for limit tmp variabels
        extern RFloat *d_dmin;
        extern RFloat *d_dmax;

        //! ql,qr for turb part
        //! extern RFloat  *d_q_turb;
        extern RFloat *d_ql_turb;
        extern RFloat *d_qr_turb;

        //! flux, flux_sub for turb part
        extern RFloat *d_flux_turb;
        extern RFloat *d_flux_sub_turb;
        //! res for turb
        extern RFloat *d_res_turb;
        //! get the refGama for device
        extern double d_refGama;
        //! get the t modle number
        extern int d_ntmodel;
        //! get the coefficientOfStateEquation
        extern RFloat d_coefficientOfStateEquation;

        //! turb face value
        extern RFloat *d_mul_turb;
        extern RFloat *d_mut_turb;
        extern RFloat *d_prim_turb;
        extern RFloat *d_mlt_turb;
        //! extern RFloat * d_flux_sub_turb;
        extern RFloat *d_flux_turb;

        //! turbulence parameters
        extern double *d_turboo;

        //! turbulence residual
        extern RFloat *d_res_turb;

        //! turbulence vars
        extern RFloat *d_spec_turb;
        extern RFloat *d_mat_turbl;
        extern RFloat *d_mat_turbr;
        //! ns vars
        extern RFloat  *d_dt;
        extern RDouble *d_CFLCell;
        //! wall distance for turbulence
        extern RDouble *d_walldist;
        //! turbulence rhs
        extern RFloat *d_rhs_turb;
        //! turbulence temp and proxy variables
        //! extern RFloat * d_q_proxy_turb;
        extern RFloat *d_q_proxy_temp_turb;
        //! turbulence update flow field atomicAdd variable
        extern int *d_n_neg;
        extern int *d_n_pos;

        //! sengy
        extern RFloat *d_sengy;
        //! farfield boundary
        extern RFloat *d_prim_inf;
        //! rhs for ns
        extern RFloat *d_rhs_ns;
        //! for cfl calc
        extern RFloat *d_invSpectrumRadius;
        extern RFloat *d_visSpectrumRadius;
        extern RFloat *d_dtv;
        extern double *d_minDt;
        extern double *d_maxDt;
        extern RFloat *d_minLocalTime; //! minimum time in one GPU (local value, not global)
        extern RFloat *d_maxLocalTime;

        //! vistmin vistmax
        extern double *d_vistmin;
        extern double *d_vistmax;

        //! ComputeNodeValue
        extern RFloat *d_qNode;
        extern RFloat *d_tNode;
        extern int    *d_nCount;
        extern RFloat *d_qTurbNode;
        extern RFloat *d_qVelNode;
        extern RFloat *d_nodeWeight;
        extern int    *d_nodeBC;
        extern RFloat *d_qOnBCFace;

        //! field value for LU-SGS
        extern RFloat *d_dq_ns;
        extern RFloat *d_dq_turb;
        extern RFloat *d_spec_ns;
        //! field value for unsteady case
        extern RFloat *d_q_ns_unsteady_n1;
        extern RFloat *d_q_ns_unsteady_n2;
        extern RFloat *d_res_ns_unsteady_n1;
        extern RFloat *d_res_ns_unsteady_n2;
        extern RFloat *d_res_ns_unsteady_tmp;

        extern RFloat *d_q_turb_unsteady_n1;
        extern RFloat *d_q_turb_unsteady_n2;
        extern RFloat *d_res_turb_unsteady_n1;
        extern RFloat *d_res_turb_unsteady_n2;
        extern RFloat *d_res_turb_unsteady_tmp;

        extern RFloat *d_sum1_ns_unsteady;
        extern RFloat *d_sum2_ns_unsteady;
        extern RFloat *d_sum1_turb_unsteady;
        extern RFloat *d_sum2_turb_unsteady;
        //! interpolation point value
        extern RFloat *d_qInterPoint;
        extern RFloat *d_tInterPoint;
        extern RFloat *d_qTurbInterPoint;

        //! transfer storage for interface value
        extern bool    d_interfaceValueExt;
        extern RFloat *d_fg_send_q;
        extern RFloat *d_fg_send_t;
        extern RFloat *d_fg_recv_q;
        extern RFloat *d_fg_recv_t;
        extern RFloat *d_fg_send_qTurb;
        extern RFloat *d_fg_recv_qTurb;
        extern RFloat *d_fg_send_dqdx;
        extern RFloat *d_fg_recv_dqdx;
        extern RFloat *d_fg_send_dqdy;
        extern RFloat *d_fg_recv_dqdy;
        extern RFloat *d_fg_send_dqdz;
        extern RFloat *d_fg_recv_dqdz;
        extern RFloat *d_fg_send_limit;
        extern RFloat *d_fg_recv_limit;
        extern RFloat *d_fg_send_dtdx;
        extern RFloat *d_fg_recv_dtdx;
        extern RFloat *d_fg_send_dtdy;
        extern RFloat *d_fg_recv_dtdy;
        extern RFloat *d_fg_send_dtdz;
        extern RFloat *d_fg_recv_dtdz;
        extern RFloat *d_fg_send_dqTurbdx;
        extern RFloat *d_fg_recv_dqTurbdx;
        extern RFloat *d_fg_send_dqTurbdy;
        extern RFloat *d_fg_recv_dqTurbdy;
        extern RFloat *d_fg_send_dqTurbdz;
        extern RFloat *d_fg_recv_dqTurbdz;
        //! transfer storage for interpoint value
        extern bool    d_interpointInfoExt;
        extern RFloat *d_fg_send_qNode;
        extern RFloat *d_fg_send_tNode;
        extern RFloat *d_fg_recv_qNode;
        extern RFloat *d_fg_recv_tNode;
        extern RFloat *d_fg_send_qTurbNode;
        extern RFloat *d_fg_recv_qTurbNode;
        extern double *HostSendBufferNS;
        extern double *HostRecvBufferNS;
        extern double *GPUSendBufferNS;
        extern double *GPURecvBufferNS;
        extern int     lengthBufferNS;
        extern double *HostSendBufferTurb;
        extern double *HostRecvBufferTurb;
        extern double *GPUSendBufferTurb;
        extern double *GPURecvBufferTurb;
        extern int     lengthBufferTurb;
        extern int    *offsetBufferNS;
        extern int    *offsetBufferTurb;
        extern int    *dimBufferNS;
        extern int    *dimBufferTurb;

        extern vector<string> *nameNSBuffer;
        extern vector<string> *nameTurbBuffer;

        extern double *HostSendInterpointBufferNS;
        extern double *HostRecvInterpointBufferNS;
        extern double *GPUSendInterpointBufferNS;
        extern double *GPURecvInterpointBufferNS;

        extern double *HostSendInterpointBufferTurb;
        extern double *HostRecvInterpointBufferTurb;
        extern double *GPUSendInterpointBufferTurb;
        extern double *GPURecvInterpointBufferTurb;

        extern int  lengthInterpointBufferNS;
        extern int  lengthInterpointBufferTurb;
        extern int *offsetInterpointBufferNS;
        extern int *offsetInterpointBufferTurb;
        extern int *dimInterpointBufferNS;
        extern int *dimInterpointBufferTurb;

        extern vector<string> *nameNSInterpointBuffer;
        extern vector<string> *nameTurbInterpointBuffer;
        //! for cuda-aware-mpi
        extern double  h_localMinDt;  //! It is original mode of CUDA_AWARE_MPI
        extern double *d_globalMinDt; //! It is for CUDA_AWARE_MPI of NEWMPI

        extern int *sendGlobalZone;
        extern int *recvProcess;
        extern int *recvGlobalZone;
        extern int *numNgbZones;
        extern int *sendGlobalZoneForPoint;
        extern int *recvProcessForPoint;
        extern int *recvGlobalZoneForPoint;
        extern int *numNgbZonesForPoint;
        extern int *offsetMPIDataNS;
        extern int *offsetMPIDataTurb;
        extern int *offsetMPIDataInterpointNS;
        extern int *offsetMPIDataInterpointTurb;
        //! large send/recv buffer for cuda aware mpi
        extern double *GPUMPIDataSendNS;
        extern double *GPUMPIDataRecvNS;
        extern double *HostMPIDataSendNS;
        extern double *HostMPIDataRecvNS;
        extern double *GPUMPIDataSendTurb;
        extern double *GPUMPIDataRecvTurb;
        extern double *HostMPIDataSendTurb;
        extern double *HostMPIDataRecvTurb;
        extern double *GPUMPIDataSendInterpointNS;
        extern double *GPUMPIDataRecvInterpointNS;
        extern double *HostMPIDataSendInterpointNS;
        extern double *HostMPIDataRecvInterpointNS;
        extern double *GPUMPIDataSendInterpointTurb;
        extern double *GPUMPIDataRecvInterpointTurb;
        extern double *HostMPIDataSendInterpointTurb;
        extern double *HostMPIDataRecvInterpointTurb;

        //! for small device buffer
        extern RFloat **GPUDataSendNS_q;
        extern RFloat **GPUDataSendNS_t;
        extern RFloat **GPUDataSendNS_dqdx;
        extern RFloat **GPUDataSendNS_dqdy;
        extern RFloat **GPUDataSendNS_dqdz;
        extern RFloat **GPUDataSendNS_limit;
        extern RFloat **GPUDataSendNS_dtdx;
        extern RFloat **GPUDataSendNS_dtdy;
        extern RFloat **GPUDataSendNS_dtdz;
        extern RFloat **GPUDataSendTurb_q;
        extern RFloat **GPUDataSendTurb_dqdx;
        extern RFloat **GPUDataSendTurb_dqdy;
        extern RFloat **GPUDataSendTurb_dqdz;

        extern RFloat **GPUDataRecvNS_q;
        extern RFloat **GPUDataRecvNS_t;
        extern RFloat **GPUDataRecvNS_dqdx;
        extern RFloat **GPUDataRecvNS_dqdy;
        extern RFloat **GPUDataRecvNS_dqdz;
        extern RFloat **GPUDataRecvNS_limit;
        extern RFloat **GPUDataRecvNS_dtdx;
        extern RFloat **GPUDataRecvNS_dtdy;
        extern RFloat **GPUDataRecvNS_dtdz;
        extern RFloat **GPUDataRecvTurb_q;
        extern RFloat **GPUDataRecvTurb_dqdx;
        extern RFloat **GPUDataRecvTurb_dqdy;
        extern RFloat **GPUDataRecvTurb_dqdz;
        //! faceIndexForSend and faceIndexForRecv on device
        extern int *d_faceIndexForSend;
        extern int *d_faceIndexForRecv;
        extern int *offsetFaceIndex;
        extern int *d_pointIndexForSend;
        extern int *d_pointIndexForRecv;
        extern int *offsetPointIndex;
        //! container for MPI_Request
        //! extern vector <MPI_Request> requestContainer;
        extern MPI_Request *requestContainerNS;
        extern MPI_Request *requestContainerNSFirst;
        extern MPI_Request *requestContainerNSRemains;
        extern int          nTotalRequestNS;
        extern MPI_Request *requestContainerTurb;
        extern MPI_Request *requestContainerTurbFirst;
        extern MPI_Request *requestContainerTurbRemains;
        extern int          nTotalRequestTurb;
        extern MPI_Request *requestContainerInterpointNS;
        extern int          nTotalRequestInterpointNS;
        extern MPI_Request *requestContainerInterpointTurb;
        extern int          nTotalRequestInterpointTurb;
        //! for macro ORGISENDIRECVORDER
        extern int *isSendOrRecv;                  //! send or recv MPI operation
        extern int *sendRecvNgbID;                 //! local neighbor zones'ID of current zone
        extern int *isSendOrRecvForPoint;          //! send or recv MPI operation
        extern int *sendRecvNgbIDForPoint;         //! local neighbor zones'ID of current zone
        extern int *ngbBufferLengthNS;             //! length of neighbor buffer
        extern int *ngbBufferLengthNSFirst;        //! length of q and t of neighbor buff
        extern int *ngbBufferLengthNSRemains;      //! length of remains
        extern int *ngbInterpointBufferLengthNS;   //! length of neighbor buffer
        extern int *ngbInterpointBufferLengthTurb; //! length of neighbor buffer
        extern int *ngbBufferLengthTurb;           //! length of neighbor buffer
        extern int *ngbBufferLengthTurbFirst;      //! length of d_qTurb of neighbor buff
        extern int *ngbBufferLengthTurbRemains;    //! length of remains

#ifdef FACEPROXYOPT
        extern FaceProxy *faceProxyNS;
        extern FaceProxy *faceProxyTurb;
#endif

        //! streams for HDCOVERLAP
        extern cudaStream_t dataStreamNS;
        extern cudaStream_t dataStreamTurb;
        extern cudaEvent_t  transferDone;
        extern cudaEvent_t  compressDoneNS;
        extern cudaEvent_t  compressDoneTurb;

        extern int *offsetNgbLBNS;
        extern int *offsetNgbLBTurb;
        extern int  offsetTagInterfaceNS;
        extern int  offsetTagInterpointNS;
        extern int  offsetTagInterfaceTurb;
        extern int  offsetTagInterpointTurb;

        //! sreams for MPI-CUDA data transfer (overlap of H2D and D2H)
        extern cudaStream_t dataTransferDToH;
        extern cudaStream_t dataTransferDToHInterpoint;
        extern cudaStream_t dataTransferHToD;
        extern cudaStream_t dataTransferHToDInterpoint;

        //! events for MPI-CUDA data transfer (overlap of H2D and D2H)
        extern cudaEvent_t downloadDToDNS_done;
        extern cudaEvent_t downloadDToDTurb_done;
        extern cudaEvent_t downloadHToDNSq_done;
        extern cudaEvent_t downloadHToDNSt_done;
        extern cudaEvent_t downloadHToDNSdqdx_done;
        extern cudaEvent_t downloadHToDNSdqdy_done;
        extern cudaEvent_t downloadHToDNSdqdz_done;
        extern cudaEvent_t downloadHToDNSlimit_done;
        extern cudaEvent_t downloadHToDNSdtdx_done;
        extern cudaEvent_t downloadHToDNSdtdy_done;
        extern cudaEvent_t downloadHToDNSdtdz_done;
        extern cudaEvent_t downloadHToDTurbq_done;
        extern cudaEvent_t downloadHToDTurbdqdx_done;
        extern cudaEvent_t downloadHToDTurbdqdy_done;
        extern cudaEvent_t downloadHToDTurbdqdz_done;

        extern cudaEvent_t downloadHToDNSqNode_done;
        extern cudaEvent_t downloadHToDNStNode_done;
        extern cudaEvent_t downloadHToDTurbqNode_done;

        //! equation number copy
        void EqnNumberCopy(int nl, int nchem, int neqn);
        void TurbEqnNumberCopy(int nEqn_turb);
        //! malloc the variables on device.
        void GPUFlowVarTurbMemAlloc(const int iSolver);

        //! q_ns variable Alloc and Copy
        void    GPUQ_NSAlloc(const int neqn, const int nTotal);
        RFloat *GPUQ_NSAlloc1(const int neqn, const int nTotal);
        //! ql_ns qr_ns variables Alloc and Copy
        void GPUQl_NSAlloc(const int neqn, const int len_ql);
        void GPUQr_NSAlloc(const int neqn, const int len_qr);
        //! call GPU memory allocation functions, in RunCFD
        void GPUFlowVarMemAlloc(const int iSolver);
        //! allocate GPU memory for field variables defined on cells
        //! void GPUFieldVarsMemAlloc(int nEqn, int nTotal);
        void GPUFieldQProxyMemAlloc(int nEqn, int nTotal);
        void GPUFieldTProxyMemAlloc(int nEqn, int nTotal);
        void GPUFieldQTurbProxyMemAlloc(int nEqn, int nTotal);
        void GPUFieldVelocityProxyMemAlloc(int nEqn, int nTotal);
        void GPUFieldQBoundaryGradientMemAlloc(int nEqn, int nBoundFace);
        //! allocate GPU memory for field variables defined on nodes
        void GPUNodeVarsMemAlloc(int nTotalNode);
        //! Gama GamaL GamaR alloc
        RFloat *GPUGamaAlloc(const int nTotal);
        RFloat *GPUGamaLAlloc(const int len_gamaL);
        RFloat *GPUGamaRAlloc(const int len_gamaR);
        //! limiter variables alloc
        RFloat *GPULIMITAlloc(const int nEqn, const int nTotal);
        RFloat *GPULimitAlloc(const int nTotal);
        //! velocity alloc, not velocity field
        void     GPUXtnYtnZtnAlloc(const int nTotalFace);
        RDouble *GPUVgnAlloc(const int nTotalFace);
        //! flux alloc
        //! residuals
        //! d_epsCell alloc
        RFloat *GPUepsCellAlloc(const int nTotalCell);
        //! dmin dmax alloc and copy
        void GPUDminDmaxAlloc(const int nTotal);
        //! face proxy alloc
        void GPUFaceProxyMemAlloc(const int d_SEG_LEN, const int nEqn);
        void GPUNSFaceVarMemAlloc(const int nEqn, const int numberOfSpecies, const int d_SEG_LEN);
        void GPUViscousCoeffMemAlloc(const int nTotal);
        void GPUMallocMinMaxTurbViscous(const int nTotalCell);
        void GPUResMemAlloc(const int neqn, const int nTotal);
        //! memory alloc for turb parts: q, ql, qr
        void GPUQLRTurbMemAlloc(const int n_turb, const int SEG_LEN);
        //! memory alloc for turb flux
        RFloat *GPUFluxTurbAlloc(const int n_turb, const int SEG_LEN);
        RFloat *GPUFluxTurbAllocSub(const int n_turb, const int SEG_LEN);
        RFloat *GPUResTurbAlloc(const int neqn, const int nTotal);

        void GPUTurbFaceVarMemAlloc(const int neqn, const int nlen);
        //! void GPUfluxTurbSubMemAlloc(const int nsize, const int neqn);
        void GPUParametersMemAllocCopy(const int neqn);
        //! void CallGPUTurbParameterSet(const double reference_density_farfield, const
        //! double turb_cmu, const double SMALL); void GPUResTurbMemAlloc(const int
        //! n_turb, const int nTotal);
        void GPUSpecTurbMemAlloc(const int n_turb, const int nTotal);
        void GPUMatTurbLRMemAlloc(const int n_turb, const int nTotalFace);
        void GPUDtMemAlloc(const int nTotal);
        void GPUWallDistMemAlloc(const int nTotalCell);
        //! for farfield variables
        RFloat *GPUPrim_infMemAlloc(const int neqn);
        //! for rhs_ns variables
        RFloat *GPURhs_nsMemAlloc(const int neqn, const int nTotal);
        RFloat *GPUInvSpectrumMalloc(const int nTotalCell);
        double *GPUDtMalloc(const int n);

        void GPUMallocMinMaxLocalTime(const int nTotalCell);
        void GPUTurbRhsMemAlloc(const int nEqn_turb, const int nTotal);
        void GPUTurbTempProxyMemAlloc(const int nEqn_turb, const int nTotal);
        void GPUTurbUpdateFlowNegPosMemAlloc();
        void GPUVistMMalloc();
        void GPUWallDistMemCopy(const int nTotalCell, const RDouble *h_wallDist);
        void GPUComputeNodeValueMemAlloc(const int nTotalNode, const int nEqn);
        void GPUComputeNodeQTurbValueMemAlloc(const int nTotalNode, const int nEqn);
        void GPUInterPointMemAlloc(const int nIPoint, const int nEqn);
        void GPUInterPointMemAllocCopy(const int nIPoint, const int nEqn, RFloat **h_qInterPoint,
                                       RFloat **h_tInterPoint);
        void GPUTurbInterPointMemAlloc(const int nIPoint, const int nEqn);
        void GPUInterfaceTransferStorageNS(const int nIFace, const int nEqn);
        void GPUInterpointTransferStorageNS(const int nIPoint, const int nEqn);
        void GPUInterfaceTransferStorageTurb(const int nIFace, const int nEqn);
        void GPUInterpointTransferStorageTurb(const int nIPoint, const int nEqn);
        //! for whole buffer
        void GPUInterfaceTransferBufferNS(InterfaceFields *interfaceFields, const int solverNS, const int nIFace);
        void GPUInterfaceTransferBufferTurb(InterfaceFields *interfaceFields, const int solverTurb, const int nIFace);
        int  GetIndexBufferNS(const string &name);
        int  GetIndexBufferTurb(const string &name);
        int  GetIndexInterpointBufferNS(const string &name);
        int  GetIndexInterpointBufferTurb(const string &name);

        void GPUTransferInterpointBufferNS(InterpointFields *interpointFields, const int solverNS, const int nIPoint);
        void GPUTransferInterpointBufferTurb(InterpointFields *interpointFields, const int solverTurb,
                                             const int nIPoint);
        //! for LU-SGS
        void GPULUSGSDqAlloc(const int nTotal, const int nEqn);
        void GPULUSGSDqTurbAlloc(const int nTotal, const int nEqn);
        void GPUUnsteadyNSAlloc(const int nTotal, const int nEqn);
        void GPUUnsteadyTurbAlloc(const int nTotal, const int nEqn);
    } //! namespace GPUFlowVariables
    //!        __global__ void GPUTurbParameterSet(const double
    //! reference_density_farfield, const double turb_cmu, const double SMALL, double
    //!* turboo);

    namespace GPUGeomVariables
    {
        extern int d_nTotalFace;
        extern int d_nTotalCell;
        //! The number of boundary faces which include interaces
        extern int d_nBoundFace;
        //! The number of interface between different zones
        extern int d_nIFace;
        //! sum of TotalCell and BoundFace
        extern int d_nTotal;
        //! points information
        extern RDouble *d_x;
        extern RDouble *d_y;
        extern RDouble *d_z;
        //! face information
        extern RDouble *d_xfn;
        extern RDouble *d_yfn;
        extern RDouble *d_zfn;
        extern RDouble *d_area;
        extern RDouble *d_xfc;
        extern RDouble *d_yfc;
        extern RDouble *d_zfc;
        extern RDouble *d_nxs;
        extern RDouble *d_nys;
        extern RDouble *d_nzs;
        //! cell center information
        extern RDouble *d_xcc;
        extern RDouble *d_ycc;
        extern RDouble *d_zcc;
        //! Cell skewness data, excluding ghosts.
        extern RDouble *d_cellSkewness;
        //! cell volume
        extern RDouble *d_vol;
        //! face2cell face2node relationship
        extern int           *d_left_cell_of_face;
        extern int           *d_right_cell_of_face;
        extern int           *d_face2node;
        extern int           *d_node_number_of_each_face;
        extern long long int *d_nodePosiFace;
        extern int           *d_boundaryType;

        //! face Cell realationship
        //! extern int * d_face_number_of_each_cell;
        extern int *d_cell2Face;
        extern int *d_posiCell2Face;

        extern int *d_cell2face;
        extern int *d_face_number_of_each_cell;
        extern int *d_acc_face_number_of_cell;

        //! cell2node for cell and node relationship
        extern int *d_cell2Node;
        extern int *d_posiCell2Node;
        extern int *d_cell2NodeCount;
        extern int *d_node_number_of_each_cell;

        //! cell2Cell for cell and its neighbor cell relationship
        extern int *d_cell2Cell;
        extern int *d_posiCell2Cell;
        extern int *d_neighbour_cell_number_of_each_cell;
        //! interpolation points infomation
        extern int *d_interPoint2GlobalPoint;
        extern int *d_cellNumberOfInterPoint;
        extern int *d_labelOfInterPoint;
        //! interface information
        extern int *d_interFace2BoundaryFace;
        //! cell coloring information
        extern int *d_InteriorCellGroup;
        //! overset cell information
        extern int *d_blank;

        //! void GPUGeomInfoAllocCopy(Grid *grid);
        //!void TestGetRegion(Region *);
        //! copy cell2Cell, neighbour_cell_number_of_each_cell from host to device
        void GPUCell2CellNumberOfEachCellCopy(const int nTotalCell, vector<int> *cell2Cell,
                                              const int *ngb_cell_number_of_each_cell);
        //! copy cell2node, node_number_of_each_cell from host to device
        void GPUCell2NodeNumberOfEachCellCopy(const int nTotalCell, const int nTotalFace, const int nBoundFace,
                                              const int *cell2node, const int *node_number_of_each_cell,
                                              const int *face2node, const int *node_number_of_each_face,
                                              const int *leftCellofFace, const int *rightCellofFace);
        //! get UnstructGrid* grid in RunCFD and call many functions for transfer
        //! geometry data from host into device.
        void GPUGeomInfoAllocCopy();
        //! copy nTotalXXX variables to gpu device
        void GPUTotalNumAllocCopy(Grid *grid);
        //! copy raw data: points coordinates to gpu device
        void GPUNodesDataAllocCopy(RDouble *x, RDouble *y, RDouble *z, const int n);
        //! copy xfn xfc xcc ... to gpu device
        void GPUFaceNormAllocCopy(RDouble *xfn, RDouble *yfn, RDouble *zfn, const int n);
        void GPUFaceCentAllocCopy(RDouble *xfc, RDouble *yfc, RDouble *zfc, const int n);
        //! copy xcc, ycc, zcc into d_xcc, d_ycc, d_zcc on device
        void GPUCellCentAllocCopy(RDouble *xcc, RDouble *ycc, RDouble *zcc, const int n);
        void GPUAreaVolmAllocCopy(RDouble *area, RDouble *vol, const int nArea, const int nVol);
        void GPUCellSkewnessAllocCopy(RDouble *cellSkewness, const int nTotalCell);
        //! copy left_cell_of_face and right_cell_of_face into d_left_cell_of_face,
        //! d_right_cell_of_face on device
        void GPUFaceCellRelAllocCopy(const int *h_left_cell_of_face, const int *h_right_cell_of_face,
                                     const int nTotalFace);
        //! copy face2node, node_number_of_each_face int d_face2node,
        //! d_node_number_of_each_face on device
        void GPUFaceNodeRelAllocCopy(const int *h_face2node, const int *h_node_number_of_each_face, int nTotalFace2Node,
                                     int nTotalFace);
        //! copy into d_nodePosiFace
        void GPUNodePosiFaceAllocCopy(const long long int *h_nodePosiFace, const int nTotalNodePosiFace);
        //! copy into d_boundaryType
        void GPUBoundaryTypeAllocCopy(const int *h_boundaryType, const int nBoundFace);
        //! copy face2node, facetocell relationship to gpu
        //! void GPUFaceRelaAllocCopy(int *left_cell_of_face, int *right_cell_of_face,
        //! int *face2node, int *node_number_of_each_face, const int n);
        void GPUFaceCellRelMemcpy(int *face_number_of_each_cell, int **cell2face, int nTotalCell);

        void GPUFaceNormalAllocCopy(const RDouble *h_nxs, const RDouble *h_nys, const RDouble *h_nzs, int nTotalFace);
        void GPUCellFaceAllocCopy(int **cell2face, const int *face_number_of_each_cell, const int nTotalCell);
        void GPUInterPointsAllocCopy(int nIPoint, const int *interPoint2GlobalPoint, const int *cellNumberOfInterPoint,
                                     const int *labelOfInterPoint);
        void GPUInterFace2BoundaryFace(const int nIFace, const int *interFace2BoundaryFace);
        void GPUInterPoint2GlobalPoint(const int nIPoint, const int *interPoint2GlobalPoint);
        void GPUCellColorAlloc();
        void GPUBlankAlloc(const int nTotalCell, const int *h_blank);
    } //! namespace GPUGeomVariables
    //! free the memory of GPU variables
    void GPUMemoryFree();
}
#endif
