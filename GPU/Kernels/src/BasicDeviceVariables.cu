#include <unistd.h>
#include "BasicDeviceVariables.h"
#include "GPUDeviceControl.h" //! for cudaMalloc d_minLocalTime
#include "GPUNSSolver.h"      //! for cudaMalloc d_maxLocalTime
#include "GPUTestFunctions.h"
#include "TemporaryOperations.h"
#ifdef CUDAUNITTEST
#include "GPUKernelTestPart2.h"
#include "GPUTestFunctions.h"
using namespace GPUTestSpace;
#endif
#include "OutputDebug.h"
#include "GPUFaceColor.h" //! for deleting device variables for GPUFaceColor
#ifdef MATRIXOUTPUT
#endif
#ifdef MCLUSGS
#include "Geo_CellColor.h"
#endif

using namespace PHSPACE;
using namespace PHMPI;
using namespace GPUTestSpace;
using namespace TemporaryOperations;
namespace GPUMemory
{
    namespace GPUControlVariables
    {
        string d_gradient_field_proxy = "Empty";
        int    d_gradient_var_index;
        string d_visflux_solver = "Empty";
        string d_FillField_turb = "Empty";
        string d_FillField_ns   = "Empty";
        string d_Limiter_var    = "Empty";
        void   SetGPUGradientControlVars(string &d_gradient_control_var, const string value)
        {
            d_gradient_control_var = value;
            //! cout<<"gradient_control_var"<<gradient_control_var<<endl;
        }

        void SetGPUGradientVarIdx(int *d_gradient_var_index, int value)
        {
            //! cout<<"In SetGPUGradientVarIdx, value="<<value<<",
            //! d_gradient_control_var="<<d_gradient_field_proxy<<endl;
            *d_gradient_var_index = value;
        }
        int  d_m_count;
        void SetGPUVisFluxSolverControlVars(string &d_solver_control_var, const string value)
        {
            d_solver_control_var = value;
        }

        void SetGPUTurbFillField(string &d_FillField_control_var, const string value)
        {
            d_FillField_control_var = value;
        }

        void SetGPUNSFillField(string &d_FillField_control_var, const string value) { d_FillField_control_var = value; }

        int  IsInviscidFluxLoad;
        int  IsFillFieldBeforeStage = 0;
        void SetGPULimiterControlVarsByValue(string &d_Limiter_control_var, const int value)
        {
            switch (value)
            {
            case 0:
                d_Limiter_control_var = "rho";
                break;
            case 4:
                d_Limiter_control_var = "pressure";
                break;
            case -1:
                d_Limiter_control_var = "Empty";
                break;
            default:
                cout << "Unexpected limiter control value\n" << endl;
                break;
            }
        }
    } //! namespace GPUControlVariables
    namespace GPUFlowVariables
    {
        RFloat *d_q_proxy;
        RFloat *d_q_proxy_tmp;
        RFloat *d_dqdx_proxy;
        RFloat *d_dqdy_proxy;
        RFloat *d_dqdz_proxy;
        RFloat *d_bdqx;
        RFloat *d_bdqy;
        RFloat *d_bdqz;
        RFloat *d_t_proxy;
        RFloat *d_dtdx_proxy;
        RFloat *d_dtdy_proxy;
        RFloat *d_dtdz_proxy;
        RFloat *d_q_turb_proxy;
        RFloat *d_dq_turbdx_proxy;
        RFloat *d_dq_turbdy_proxy;
        RFloat *d_dq_turbdz_proxy;
        RFloat *d_vel_proxy;
        RFloat *d_dveldx_proxy;
        RFloat *d_dveldy_proxy;
        RFloat *d_dveldz_proxy;
        RFloat *d_q_old;
        RFloat *d_q_n_tmp;
        //! double *d_q_n_double;s
        int *d_n_count_tmp;
        int  d_nTVar_m;
        //! equation numbers
        int d_nl        = 0;
        int d_nchem     = 0;
        int d_neqn_ns   = 0;
        int d_neqn_turb = 0;

        //! SEG_LEN On device
        int d_SEG_LEN = 0;
        //! flow variables
        RFloat *d_q_ns;
        RFloat *d_ql_ns;
        RFloat *d_qr_ns;
        //! flow variables gamma
        RFloat *d_gama_ns;
        RFloat *d_gamaL_ns;
        RFloat *d_gamaR_ns;
        //! limiter variables
        RFloat *d_limit;
        RFloat *d_LIMIT;
        //! vgn FaceNormalVelocity
        RDouble *d_vgn;
        //! xtn ytn ztn NodeVelocity
        RDouble *d_xtn;
        RDouble *d_ytn;
        RDouble *d_ztn;
        //! residuals
        RFloat *d_res_ns;
        //! tmp variable: epsCell
        RFloat *d_epsCell;
        //! dmin dmax for limit tmp variabels
        RFloat *d_dmin;
        RFloat *d_dmax;
        //! viscous coeff
        RFloat *d_visl;
        RFloat *d_vist;
        RFloat *d_minLocalTurbVis; //! local minimum turb viscous in one GPU
        RFloat *d_maxLocalTurbVis;
        //! face proxy varialbes
        RFloat  *d_flux;
        RDouble *d_deltl;
        RDouble *d_deltr;
        //! ns face value variables
        RFloat *d_prim;
        RFloat *d_tm;
        RFloat *d_kcp;
        RFloat *d_mul;
        RFloat *d_mut;
        RFloat *d_rho_ds;
        RFloat *d_hint_s;
        //! res

        //! q for turb part
        //! RFloat  *d_q_turb;
        RFloat *d_ql_turb;
        RFloat *d_qr_turb;
        //! res for turb
        RFloat *d_res_turb;
        //! get the refGama for device
        double d_refGama;
        //! get the t modle number
        int d_ntmodel;
        //! get the coefficientOfStateEquation
        RFloat d_coefficientOfStateEquation;
        //! turb face value
        RFloat *d_mul_turb;
        RFloat *d_mut_turb;
        RFloat *d_prim_turb;
        RFloat *d_mlt_turb;
        RFloat *d_flux_sub_turb;
        RFloat *d_flux_turb;

        //! turbulence parameters
        double *d_turboo;

        //! turbulence vars
        RFloat *d_spec_turb;
        RFloat *d_mat_turbl;
        RFloat *d_mat_turbr;

        //! ns vars
        RFloat  *d_dt;
        RDouble *d_CFLCell;
        //! wall distance for turbulence
        RDouble *d_walldist;
        //! farfield boundary
        RFloat *d_prim_inf;
        RFloat *d_rhs_ns;
        RFloat *d_invSpectrumRadius;
        RFloat *d_visSpectrumRadius;
        RFloat *d_dtv;
        double *d_minDt;
        double *d_maxDt;
        RFloat *d_minLocalTime; //! minimum time in one GPU (local value, not global)
        RFloat *d_maxLocalTime;
        //! turbulence rhs
        RFloat *d_rhs_turb;
        //! turbulence temp and proxy variables
        //! RFloat * d_q_proxy_turb;
        RFloat *d_q_proxy_temp_turb;
        int    *d_n_neg;
        int    *d_n_pos;
        //! vistmin vistmax
        double *d_vistmin;
        double *d_vistmax;

        //! ComputeNodeValue
        RFloat *d_qNode;
        RFloat *d_tNode;
        int    *d_nCount;
        RFloat *d_qTurbNode;
        RFloat *d_qVelNode;
        RFloat *d_nodeWeight;
        int    *d_nodeBC;
        RFloat *d_qOnBCFace;
        //! field value for LU-SGS
        RFloat *d_dq_ns;
        RFloat *d_dq_turb;
        RFloat *d_spec_ns;
        //! field value for unsteady case
        RFloat *d_q_ns_unsteady_n1;
        RFloat *d_q_ns_unsteady_n2;
        RFloat *d_res_ns_unsteady_n1;
        RFloat *d_res_ns_unsteady_n2;
        RFloat *d_res_ns_unsteady_tmp;

        RFloat *d_q_turb_unsteady_n1;
        RFloat *d_q_turb_unsteady_n2;
        RFloat *d_res_turb_unsteady_n1;
        RFloat *d_res_turb_unsteady_n2;
        RFloat *d_res_turb_unsteady_tmp;

        RFloat *d_sum1_ns_unsteady;
        RFloat *d_sum2_ns_unsteady;
        RFloat *d_sum1_turb_unsteady;
        RFloat *d_sum2_turb_unsteady;
        //! interpolation points value
        RFloat *d_qInterPoint;
        RFloat *d_tInterPoint;
        RFloat *d_qTurbInterPoint;

        //! interface transfer storage
        bool    d_interfaceValueExt = false;
        RFloat *d_fg_send_q;
        RFloat *d_fg_send_t;
        RFloat *d_fg_recv_q;
        RFloat *d_fg_recv_t;
        RFloat *d_fg_send_qTurb;
        RFloat *d_fg_recv_qTurb;
        RFloat *d_fg_send_dqdx;
        RFloat *d_fg_recv_dqdx;
        RFloat *d_fg_send_dqdy;
        RFloat *d_fg_recv_dqdy;
        RFloat *d_fg_send_dqdz;
        RFloat *d_fg_recv_dqdz;
        RFloat *d_fg_send_limit;
        RFloat *d_fg_recv_limit;
        RFloat *d_fg_send_dtdx;
        RFloat *d_fg_recv_dtdx;
        RFloat *d_fg_send_dtdy;
        RFloat *d_fg_recv_dtdy;
        RFloat *d_fg_send_dtdz;
        RFloat *d_fg_recv_dtdz;
        RFloat *d_fg_send_dqTurbdx;
        RFloat *d_fg_recv_dqTurbdx;
        RFloat *d_fg_send_dqTurbdy;
        RFloat *d_fg_recv_dqTurbdy;
        RFloat *d_fg_send_dqTurbdz;
        RFloat *d_fg_recv_dqTurbdz;
        //! interface transfer storage for interpoint value
        bool    d_interpointValueExt = false;
        RFloat *d_fg_send_qNode;
        RFloat *d_fg_send_tNode;
        RFloat *d_fg_recv_qNode;
        RFloat *d_fg_recv_tNode;
        RFloat *d_fg_send_qTurbNode;
        RFloat *d_fg_recv_qTurbNode;

        vector<string> *nameNSBuffer;
        vector<string> *nameTurbBuffer;
        vector<string> *nameNSInterpointBuffer;
        vector<string> *nameTurbInterpointBuffer;

        double *HostSendBufferNS;
        double *HostRecvBufferNS;
        double *GPUSendBufferNS;
        double *GPURecvBufferNS;
        int     lengthBufferNS;
        double *HostSendBufferTurb;
        double *HostRecvBufferTurb;
        double *GPUSendBufferTurb;
        double *GPURecvBufferTurb;
        int     lengthBufferTurb;
        int    *offsetBufferNS;
        int    *offsetBufferTurb;
        int    *dimBufferNS;
        int    *dimBufferTurb;

        double *HostSendInterpointBufferNS;
        double *HostRecvInterpointBufferNS;
        double *GPUSendInterpointBufferNS;
        double *GPURecvInterpointBufferNS;
        int     lengthInterpointBufferNS;
        double *HostSendInterpointBufferTurb;
        double *HostRecvInterpointBufferTurb;
        double *GPUSendInterpointBufferTurb;
        double *GPURecvInterpointBufferTurb;
        int     lengthInterpointBufferTurb;
        int    *offsetInterpointBufferNS;
        int    *offsetInterpointBufferTurb;
        int    *dimInterpointBufferNS;
        int    *dimInterpointBufferTurb;
        //! for cuda-aware-mpi
        double  h_localMinDt;  //! It is for NEWMPI original mode
        double *d_globalMinDt; //! It is for CUDA_AWARE_MPI of NEWMPI

        int *sendGlobalZone;
        int *recvProcess;
        int *recvGlobalZone;
        int *numNgbZones;
        int *sendGlobalZoneForPoint;
        int *recvProcessForPoint;
        int *recvGlobalZoneForPoint;
        int *numNgbZonesForPoint;
        int *offsetMPIDataNS;
        int *offsetMPIDataTurb;
        int *offsetMPIDataInterpointNS;
        int *offsetMPIDataInterpointTurb;

        double *GPUMPIDataSendNS;
        double *GPUMPIDataRecvNS;
        double *HostMPIDataSendNS;
        double *HostMPIDataRecvNS;
        double *GPUMPIDataSendTurb;
        double *GPUMPIDataRecvTurb;
        double *HostMPIDataSendTurb;
        double *HostMPIDataRecvTurb;
        double *GPUMPIDataSendInterpointNS;
        double *GPUMPIDataRecvInterpointNS;
        double *HostMPIDataSendInterpointNS;
        double *HostMPIDataRecvInterpointNS;
        double *GPUMPIDataSendInterpointTurb;
        double *GPUMPIDataRecvInterpointTurb;
        double *HostMPIDataSendInterpointTurb;
        double *HostMPIDataRecvInterpointTurb;

        //! for small device buffer
        RFloat **GPUDataSendNS_q;
        RFloat **GPUDataSendNS_t;
        RFloat **GPUDataSendNS_dqdx;
        RFloat **GPUDataSendNS_dqdy;
        RFloat **GPUDataSendNS_dqdz;
        RFloat **GPUDataSendNS_limit;
        RFloat **GPUDataSendNS_dtdx;
        RFloat **GPUDataSendNS_dtdy;
        RFloat **GPUDataSendNS_dtdz;
        RFloat **GPUDataSendTurb_q;
        RFloat **GPUDataSendTurb_dqdx;
        RFloat **GPUDataSendTurb_dqdy;
        RFloat **GPUDataSendTurb_dqdz;

        RFloat **GPUDataRecvNS_q;
        RFloat **GPUDataRecvNS_t;
        RFloat **GPUDataRecvNS_dqdx;
        RFloat **GPUDataRecvNS_dqdy;
        RFloat **GPUDataRecvNS_dqdz;
        RFloat **GPUDataRecvNS_limit;
        RFloat **GPUDataRecvNS_dtdx;
        RFloat **GPUDataRecvNS_dtdy;
        RFloat **GPUDataRecvNS_dtdz;
        RFloat **GPUDataRecvTurb_q;
        RFloat **GPUDataRecvTurb_dqdx;
        RFloat **GPUDataRecvTurb_dqdy;
        RFloat **GPUDataRecvTurb_dqdz;
        //! faceIndexForSend and faceIndexForRecv on device
        int *d_faceIndexForSend;
        int *d_faceIndexForRecv;
        int *offsetFaceIndex;
        int *d_pointIndexForSend;
        int *d_pointIndexForRecv;
        int *offsetPointIndex;
        //! container for MPI_Request
        //! vector <MPI_Request> requestContainer;
        MPI_Request *requestContainerNS;
        //! divide requestContainerNS into two parts for MPIOVERLAP
        MPI_Request *requestContainerNSFirst;
        MPI_Request *requestContainerNSRemains;
        int          nTotalRequestNS;
        MPI_Request *requestContainerTurb;
        //! divide requestContainerTurb into two parts for MPIOVERLAP
        MPI_Request *requestContainerTurbFirst;
        MPI_Request *requestContainerTurbRemains;
        int          nTotalRequestTurb;
        MPI_Request *requestContainerInterpointNS;
        int          nTotalRequestInterpointNS;
        MPI_Request *requestContainerInterpointTurb;
        int          nTotalRequestInterpointTurb;
        //! for macro ORGISENDIRECVORDER
        int *isSendOrRecv;          //! send or recv MPI operation
        int *sendRecvNgbID;         //! local neighbor zones'ID of current zone
        int *isSendOrRecvForPoint;  //! send or recv MPI operation
        int *sendRecvNgbIDForPoint; //! local neighbor zones'ID of current zone
        int *ngbBufferLengthNS;
        int *ngbBufferLengthNSFirst;
        int *ngbBufferLengthNSRemains;
        int *ngbBufferLengthTurb;
        int *ngbBufferLengthTurbFirst;
        int *ngbBufferLengthTurbRemains;
        int *ngbInterpointBufferLengthNS;
        int *ngbInterpointBufferLengthTurb;

#ifdef FACEPROXYOPT
        FaceProxy *faceProxyNS;
        FaceProxy *faceProxyTurb;
#endif
        cudaStream_t dataStreamNS;
        cudaStream_t dataStreamTurb;
        cudaEvent_t  transferDone;
        cudaEvent_t  compressDoneNS;
        cudaEvent_t  compressDoneTurb;

        int *offsetNgbLBNS;
        int *offsetNgbLBTurb;
        int  offsetTagInterfaceNS;
        int  offsetTagInterpointNS;
        int  offsetTagInterfaceTurb;
        int  offsetTagInterpointTurb;

        cudaStream_t dataTransferDToH;
        cudaStream_t dataTransferDToHInterpoint;
        cudaStream_t dataTransferHToD;
        cudaStream_t dataTransferHToDInterpoint;

        cudaEvent_t downloadDToDNS_done;
        cudaEvent_t downloadDToDTurb_done;
        cudaEvent_t downloadHToDNSq_done;
        cudaEvent_t downloadHToDNSt_done;
        cudaEvent_t downloadHToDNSdqdx_done;
        cudaEvent_t downloadHToDNSdqdy_done;
        cudaEvent_t downloadHToDNSdqdz_done;
        cudaEvent_t downloadHToDNSlimit_done;
        cudaEvent_t downloadHToDNSdtdx_done;
        cudaEvent_t downloadHToDNSdtdy_done;
        cudaEvent_t downloadHToDNSdtdz_done;
        cudaEvent_t downloadHToDTurbq_done;
        cudaEvent_t downloadHToDTurbdqdx_done;
        cudaEvent_t downloadHToDTurbdqdy_done;
        cudaEvent_t downloadHToDTurbdqdz_done;

        cudaEvent_t downloadHToDNSqNode_done;
        cudaEvent_t downloadHToDNStNode_done;
        cudaEvent_t downloadHToDTurbqNode_done;

        void GPUFlowVarMemAlloc(const int solverNS)
        {
//! just for test
//! getCodeID(numberOfProcessors);
#ifdef CUDAUNITTEST
#ifdef UNITTESTOUTPUT
            printf("program is running in GPUFlowVarMemAlloc\n");
#endif
#endif
            int nLocalZones = GetNumberofLocalZones();
            for (int izone = 0; izone < nLocalZones; izone++)
            {
                int zoneID = GetLocalZoneIDToGlobalZoneID(izone);
                //! Zone *zone = GetZone(zoneID);
                int   level   = 0;
                Grid *grid_in = GetGrid(zoneID, level);
                //! UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, level));
                UnstructGrid *grid       = UnstructGridCast(grid_in);
                const int     nTotalCell = grid->GetNTotalCell();
                const int     nBoundFace = grid->GetNBoundFace();
                const int     nTotalFace = grid->GetNTotalFace();
#ifdef CUDAUNITTEST
#ifdef UNITTESTOUTPUT
                printf("nTotalFace = %d\n", nTotalFace);
#endif
#endif
                int nTotal = nTotalCell + nBoundFace;
                int nl;
                GlobalDataBase::GetData("nl", &nl, PHINT, 1);
                int nchem;
                GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
                int nEqn = nl + nchem;

                //! int ntmodel= d_ntmodel; //! the d_ntmodel is initialize in
                //! CompGamaAndTField function during initialize
                int ntmodel = 1; //! the d_ntmodel is initialize in CompGamaAndTField
                                 //! function during initialize

                GPUFieldQProxyMemAlloc(nEqn, nTotal);
                GPUFieldTProxyMemAlloc(ntmodel, nTotal);
                GPUFieldVelocityProxyMemAlloc(3, nTotal);
                GPUFieldQBoundaryGradientMemAlloc(nEqn, nBoundFace);
                int nTotalNode = grid->GetNTotalNode();
                //! Canceled due to new method GGNODENEW
                //! GPUNodeVarsMemAlloc(nTotalNode);
                //! Just for test

                //! copy nl nchem neqn to d_xx respectively
                EqnNumberCopy(nl, nchem, nEqn);
                //! set the SEG_LEN on device
                d_SEG_LEN = SEGCTION_LENGTH;
                //! cp q ql qr  variables onto GPU device
                GPUQ_NSAlloc(nEqn, nTotal);
                GPUQl_NSAlloc(nEqn, d_SEG_LEN);
                GPUQr_NSAlloc(nEqn, d_SEG_LEN);
                //! cp gama gamaL gamaR onto device
                d_gama_ns  = GPUGamaAlloc(nTotal);
                d_gamaL_ns = GPUGamaLAlloc(d_SEG_LEN);
                d_gamaR_ns = GPUGamaRAlloc(d_SEG_LEN);
                //! memory alloc for limiter
                d_LIMIT = GPULIMITAlloc(nEqn, nTotal);
                d_limit = GPULimitAlloc(nTotal);
                //! velocity alloc
                GPUXtnYtnZtnAlloc(nTotalFace);
                d_vgn = GPUVgnAlloc(nTotalFace);
                //! d_epsCell alloc
                d_epsCell = GPUepsCellAlloc(nTotalCell);

                //! malloc face proxy variables
                GPUFaceProxyMemAlloc(d_SEG_LEN, nEqn);
                int numberOfSpecies;
                GlobalDataBase::GetData("numberOfSpecies", &numberOfSpecies, PHINT, 1);
#ifdef CUDAUNITTEST
                //! cout<<"numberOfSpecies ="<<numberOfSpecies<<endl;
#endif
                GPUNSFaceVarMemAlloc(nEqn, numberOfSpecies, d_SEG_LEN);
                GPUViscousCoeffMemAlloc(nTotal);
                GPUMallocMinMaxTurbViscous(nTotalCell);
                GPUResMemAlloc(nEqn, nTotal);
                GPUDminDmaxAlloc(nTotal);

                /*copy xtn, ytn, ztn*/
                GPUXtnYtnZtnCopy(grid, nTotalFace);
                //! copy Vgn onto device
                GPUVgnCopy(grid, nTotalFace);
                GPUDtMemAlloc(nTotal);
                GPUWallDistMemAlloc(nTotalCell);
                //! for farfield variables
                d_prim_inf          = GPUPrim_infMemAlloc(nEqn);
                d_rhs_ns            = GPURhs_nsMemAlloc(nEqn, nTotal);
                d_invSpectrumRadius = GPUInvSpectrumMalloc(nTotalCell);
                d_visSpectrumRadius = GPUInvSpectrumMalloc(nTotalCell);
                d_dtv               = GPUInvSpectrumMalloc(nTotalCell);
                d_minDt             = GPUDtMalloc(1);
                d_maxDt             = GPUDtMalloc(1);
                GPUMallocMinMaxLocalTime(nTotalCell);
                //! wall dist copy
                //! RDouble * h_wallDist = grid_in->GetWallDist();
                RDouble *h_wallDist = grid->GetWallDist();
                GPUWallDistMemCopy(nTotalCell, h_wallDist);
                GPUComputeNodeValueMemAlloc(nTotalNode, nEqn);
                //! alloc device memory for interpolation points
                InterpointInformation *interPointInfor = grid->GetInterpointInfo();
                if (interPointInfor)
                {
                    int      nIPoint     = interPointInfor->GetNumberOfInterpoints();
                    RFloat **qInterPoint = reinterpret_cast<RFloat **>(grid->GetDataPtr("qInterPoint"));
                    RFloat **tInterPoint = reinterpret_cast<RFloat **>(grid->GetDataPtr("tInterPoint"));
                    GPUInterPointMemAllocCopy(nIPoint, nEqn, qInterPoint, tInterPoint);
                }
                //! allocate device memroy for interface transfer storage
                InterfaceInfo *iinfo = grid->GetInterfaceInfo();
                //! only for multi-process
                //! if (iinfo) {
                //! CUDA_AWARE_MPI may be used in one process condition
#ifdef CUDA_AWARE_MPI
                HANDLE_API_ERR(cudaMalloc((void **)&d_globalMinDt, 1 * sizeof(double)));
#endif
                //!}

                if (iinfo)
                {
                    d_interfaceValueExt = true;
                    int nIFace          = iinfo->GetNIFace();
#ifdef MPIBUFFERSMALL
                    GPUInterfaceTransferStorageNS(nIFace, nEqn);
#endif
#ifdef MPIBUFFERLARGE
                    InterfaceFields *interfaceFields = grid->GetInterfaceFields();
                    GPUInterfaceTransferBufferNS(interfaceFields, solverNS, nIFace);
#ifdef HOSTDEVOVERLAP
                    //! small device send buffer is required
                    GPUInterfaceTransferStorageNS(nIFace, nEqn);
                    //! It should be noted that steams for H2D and D2H overlap are used for
                    //! both NS branch and Turb branch
                    HANDLE_API_ERR(cudaStreamCreateWithFlags(&dataTransferDToH, cudaStreamNonBlocking));
                    HANDLE_API_ERR(cudaStreamCreateWithFlags(&dataTransferHToD, cudaStreamNonBlocking));
                    //! Events just used for recording data transfer condition in NS branch
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadDToDNS_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDNSq_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDNSt_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDNSdqdx_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDNSdqdy_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDNSdqdz_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDNSlimit_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDNSdtdx_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDNSdtdy_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDNSdtdz_done, cudaEventDisableTiming));
#endif
#endif
#ifdef NEWMPI
                    //! just for test
                    int numberOfProcessors = PHMPI::GetNumberOfProcessor();
                    //! set offsetTagInterfaceNS for tags used in Isend and Irecv
                    offsetTagInterfaceNS = 0;
                    //! Getting zone connectivity information
                    int nTotalGloablZones = GetNumberofGlobalZones();
                    //! printf("nTotalGloablZones = %d\n", nTotalGloablZones);
                    int currentProcessorID  = GetCurrentProcessorID();
                    sendGlobalZone          = new int[1]; //! only consider one zone in one process
                    numNgbZones             = new int[1];
                    int *h_faceIndexForSend = new int[nIFace];
                    int *h_faceIndexForRecv = new int[nIFace];
                    for (int iZone = 0; iZone < nTotalGloablZones; iZone++)
                    {
                        ZoneNeighbor *zoneNeighbor     = zoneConnectivity->GetZoneNeighbor(iZone);
                        std::size_t   numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();
                        int           sendProcessID    = GetZoneProcessorIDSepMode(iZone);
                        if (sendProcessID == currentProcessorID)
                        {
                            sendGlobalZone[0]              = iZone;
                            ZoneNeighbor *zoneNeighbor     = zoneConnectivity->GetZoneNeighbor(iZone);
                            std::size_t   numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();
                            recvGlobalZone                 = new int[numberOfNeighbor];
                            recvProcess                    = new int[numberOfNeighbor];
                            offsetMPIDataNS                = new int[numberOfNeighbor];
                            offsetFaceIndex                = new int[numberOfNeighbor];
                            int sumNgbLengthBuffer         = 0;
                            int sumNgbFaces                = 0;
                            int countFace                  = 0;
                            numNgbZones[0]                 = numberOfNeighbor;

                            ngbBufferLengthNS = new int[numberOfNeighbor];
                            //! MPI_Request * requestContainerNS;
                            nTotalRequestNS    = 2 * numberOfNeighbor; //! send and receive
                            requestContainerNS = new MPI_Request[nTotalRequestNS];
                            for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++iNeighbor)
                            {
                                int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(iNeighbor);
                                //! the ngbID'th neighbor zone of iZone
                                int ngbID = iinfo->FindIthNeighbor(neighborZoneIndex);
                                //! just check whether or not iNeighbor is equal to ngbID
                                if (ngbID != iNeighbor)
                                {
                                    printf("ngbID != iNeighbor, ngbID = %d, iNeighbor = %d, iZone = %d, "
                                           "neighborZoneIndex = %d\n",
                                           ngbID, iNeighbor, iZone, neighborZoneIndex);
                                    exit(1);
                                }
                                //! else printf("ngbID is equal to iNeighbor!\n");
                                int  nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(ngbID);
                                int *faceIndexForSend = iinfo->GetFaceIndexForSend(ngbID);
                                int *faceIndexForRecv = iinfo->GetFaceIndexForRecv(ngbID);

                                int receiveProcessID      = GetZoneProcessorIDSepMode(neighborZoneIndex);
                                recvGlobalZone[iNeighbor] = neighborZoneIndex;
                                recvProcess[iNeighbor]    = receiveProcessID;
                                //! calculate offsetMPIDataNS
                                offsetMPIDataNS[iNeighbor] = sumNgbLengthBuffer;
                                offsetFaceIndex[iNeighbor] = sumNgbFaces;
                                sumNgbFaces += nIFaceOfNeighbor;
                                ngbBufferLengthNS[iNeighbor] = 0;
                                for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                                {
                                    int dimData = dimBufferNS[nameID];

                                    sumNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                                    ngbBufferLengthNS[iNeighbor] += dimData * nIFaceOfNeighbor;
                                } //! end nameID of buffer
                                for (int jFaceOfNeighbor = 0; jFaceOfNeighbor < nIFaceOfNeighbor; ++jFaceOfNeighbor)
                                {
                                    int sendIFace = faceIndexForSend[jFaceOfNeighbor];
                                    int recvIFace = faceIndexForRecv[jFaceOfNeighbor];
                                    //! just for check
                                    if (sendIFace != recvIFace)
                                    {
                                        printf("iZone = %d, ngbID = %d, countFace = %d\n", iZone, ngbID, countFace);
                                    }
                                    h_faceIndexForSend[countFace] = sendIFace;
                                    h_faceIndexForRecv[countFace] = recvIFace;
                                    countFace++;
                                }
                            } //! end for iNeighbor
                            //! Test sumNgbLengthBuffer and sumNgbFaces. It should be noted that
                            //! sumNgbLengthBuffer is not equal to lengthBufferNS, becasue of the
                            //! varaible limit printf("sumNgbLengthBuffer = %d, lengthBufferNS =
                            //! %d\n", sumNgbLengthBuffer, lengthBufferNS); printf("sumNgbFaces =
                            //! %d, nIFace = %d\n", sumNgbFaces, nIFace); Test countFace and
                            //! nIFace, they should be equal. printf("countFace = %d, nIFace =
                            //! %d\n", countFace, nIFace);
                        } //! end for sendProcessID == currentProcessorID
                    } //! end of for iZone
#ifdef MPIOVERLAP

                    requestContainerNSFirst   = new MPI_Request[nTotalRequestNS];
                    requestContainerNSRemains = new MPI_Request[nTotalRequestNS];
                    //! compute ngbBufferLengthNSFirst and ngbBufferLengthNSRemains
                    int numberOfNeighbor       = numNgbZones[0];
                    ngbBufferLengthNSFirst     = new int[numberOfNeighbor];
                    ngbBufferLengthNSRemains   = new int[numberOfNeighbor];
                    int           iZone        = sendGlobalZone[0]; //! only consider one zone in one process
                    ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
                    for (int ngbID = 0; ngbID < numberOfNeighbor; ngbID++)
                    {
                        int ngbZoneID                 = recvGlobalZone[ngbID];
                        int neighborZoneIndex         = zoneNeighbor->GetZoneIndexOfNeighbor(ngbID);
                        int iNeighbor                 = iinfo->FindIthNeighbor(neighborZoneIndex);
                        int nIFaceOfNeighbor          = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                        ngbBufferLengthNSFirst[ngbID] = 0;

                        for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                        {
                            const string &name = (*nameNSBuffer)[nameID];
                            //! just for test
                            //! if ((ngbID == 0)&&(globalRank==0)) cout<<"nameID ="<<nameID<<",
                            //! name ="<<name<<endl;
                            if ((name == "q") || (name == "t"))
                            {
                                int dimData = dimBufferNS[nameID];
                                ngbBufferLengthNSFirst[ngbID] += dimData * nIFaceOfNeighbor;
                                //! if ((ngbID == 0)&&(globalRank==0)) printf("nameID = %d,
                                //! ngbBufferLengthNSFirst = %d\n", nameID,
                                //! ngbBufferLengthNSFirst[ngbID]);
                            }
                        }

                        ngbBufferLengthNSRemains[ngbID] = ngbBufferLengthNS[ngbID] - ngbBufferLengthNSFirst[ngbID];
                        //! just for test
                        //! if ((ngbID == 0)&&(globalRank==0)) printf("First = %d, remains = %d,
                        //! total = %d\n", ngbBufferLengthNSFirst[ngbID],
                        //! ngbBufferLengthNSRemains[ngbID], ngbBufferLengthNS[ngbID]);
                    }
#endif
//! create GPU small buffer for interface data transfer
#ifdef GPUBUFFERSMALL
                    //! Create
                    //! Get total number of neighbor zone of iZone
                    int nTotalNgbZones  = numNgbZones[0];
                    GPUDataSendNS_q     = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataSendNS_t     = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataSendNS_limit = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataSendNS_dqdx  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataSendNS_dqdy  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataSendNS_dqdz  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataSendNS_dtdx  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataSendNS_dtdy  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataSendNS_dtdz  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));

                    GPUDataRecvNS_q     = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataRecvNS_t     = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataRecvNS_limit = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataRecvNS_dqdx  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataRecvNS_dqdy  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataRecvNS_dqdz  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataRecvNS_dtdx  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataRecvNS_dtdy  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataRecvNS_dtdz  = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));

                    for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
                    {
                        int ngbZoneID = recvGlobalZone[ngbID];
                        //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
                        //! ngbID=iNeighbor
                        int iNeighbor        = iinfo->FindIthNeighbor(ngbZoneID);
                        int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                        for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                        {
                            const string &name = (*nameNSBuffer)[nameID];
                            if (name == "q")
                            {
                                int dimData = dimBufferNS[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendNS_q[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvNS_q[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                            else if (name == "t")
                            {
                                int dimData = dimBufferNS[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendNS_t[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvNS_t[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                            else if (name == "limit")
                            {
                                int dimData = dimBufferNS[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendNS_limit[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvNS_limit[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                            else if (name == "dqdx")
                            {
                                int dimData = dimBufferNS[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendNS_dqdx[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvNS_dqdx[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                            else if (name == "dqdy")
                            {
                                int dimData = dimBufferNS[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendNS_dqdy[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvNS_dqdy[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                            else if (name == "dqdz")
                            {
                                int dimData = dimBufferNS[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendNS_dqdz[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvNS_dqdz[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                            else if (name == "dtdx")
                            {
                                int dimData = dimBufferNS[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendNS_dtdx[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvNS_dtdx[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                            else if (name == "dtdy")
                            {
                                int dimData = dimBufferNS[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendNS_dtdy[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvNS_dtdy[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                            else if (name == "dtdz")
                            {
                                int dimData = dimBufferNS[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendNS_dtdz[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvNS_dtdz[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                        }
                    }
//! HtoD, DtoH, and computing overlap should be based on GPUSMALLBUFFER
#ifdef HDCOVERLAP
                    HANDLE_API_ERR(cudaStreamCreateWithFlags(&dataStreamNS, cudaStreamNonBlocking));

                    HANDLE_API_ERR(cudaEventCreateWithFlags(&transferDone, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&compressDoneNS, cudaEventDisableTiming));
                    //! compute offsetNgbLBNS that is global offsetNgbLengthBuffer
                    offsetNgbLBNS = new int[nTotalNgbZones * nameNSBuffer->size()];
                    for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
                    {
                        int ngbZoneID             = recvGlobalZone[ngbID];
                        int neighborZoneIndex     = zoneNeighbor->GetZoneIndexOfNeighbor(ngbID);
                        int iNeighbor             = iinfo->FindIthNeighbor(neighborZoneIndex);
                        int nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                        int offsetNgbLengthBuffer = 0;
                        //! offsetNgbLBNS[ngbID*nameNSBuffer->size()+0] = 0;

                        //! for(int nameID = 0; nameID < nameNSBuffer->size(); nameID++){
                        for (int nameID = 0; nameID < nameNSBuffer->size(); nameID++)
                        {
                            int dimData                                          = dimBufferNS[nameID];
                            offsetNgbLBNS[ngbID * nameNSBuffer->size() + nameID] = offsetNgbLengthBuffer;
                            //! offsetNgbLBNS[ngbID*nameNSBuffer->size()+nameID+1] =
                            //! offsetNgbLBNS[ngbID*nameNSBuffer->size()+nameID] + dimData *
                            //! nIFaceOfNeighbor;
                            offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                        }
                    }
                    /*
      //!Test for offsetNgbLBNS
      for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++){
              int ngbZoneID = recvGlobalZone[ngbID];
              int neighborZoneIndex =
      zoneNeighbor->GetZoneIndexOfNeighbor(ngbID); int iNeighbor =
      iinfo->FindIthNeighbor(neighborZoneIndex); int nIFaceOfNeighbor =
      iinfo->GetNIFaceOfNeighbor(iNeighbor); printf("globalRank = %d, ngbID =
      %d, nIFaceOfNeighbor =%d\n", globalRank, ngbID, nIFaceOfNeighbor); for(int
      nameID = 0; nameID < nameNSBuffer->size(); nameID++){ printf("globalRank =
      %d, nameID = %d, ngbID = %d, offsetNgbLBNS=%d\n", globalRank, nameID,
      ngbID, offsetNgbLBNS[ngbID*nameNSBuffer->size()+nameID]);
              }

      }
      */
#endif

#endif
                    /*
      //!test zone connectivity
      int iZone = sendGlobalZone[0];
      int nTotalNgbZones = numNgbZones[0];
      printf("Zone %d owns %d neighbor zones\n", iZone, nTotalNgbZones);
      for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++){
              printf("Zone %d owns %d neighbor zone in process %d\n", iZone,
      recvGlobalZone[ngbID], recvProcess[ngbID]);
      }
      */

                    /*
      //!Added for test nameID and zone property
      if (globalRank == 0){
              for(int nameID = 0; nameID < nameNSBuffer->size(); nameID++){
                      int dimData = dimBufferNS[nameID];
                      const string & name = (*nameNSBuffer)[nameID];
                      cout<<"nameID="<<nameID<<" name="<<name<<"
      dim="<<dimData<<endl; } //!end nameID of buffer
      }

      int iZone = sendGlobalZone[0]; //!only consider one zone in one process
      ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
      int numberOfNeighbor = numNgbZones[0];
      for (int ngbID = 0; ngbID < numberOfNeighbor; ngbID++){
              int ngbZoneID = recvGlobalZone[ngbID];
              int neighborZoneIndex =
      zoneNeighbor->GetZoneIndexOfNeighbor(ngbID); int iNeighbor =
      iinfo->FindIthNeighbor(neighborZoneIndex); int nIFaceOfNeighbor =
      iinfo->GetNIFaceOfNeighbor(iNeighbor); int ngbProcessID =
      recvProcess[ngbID]; int offsetFace = offsetFaceIndex[ngbID]; int offsetNgb
      = offsetMPIDataNS[ngbID];
              printf("globalRank=%d,localRank=%d,iZone=%d,ngbID=%d,ngbZoneID=%d,ngbProcessID=%d,offsetFace=%d,offsetNgb=%d,
      nIFaceOfNeighbor=%d,
      nIFace=%d\n",globalRank,localRank,iZone,ngbID,ngbZoneID,ngbProcessID,offsetFace,offsetNgb,
      nIFaceOfNeighbor, nIFace);
      }
      //!test end
      */
                    size_t sizeBufferNS = lengthBufferNS * sizeof(double);
                    //! printf("lengthBufferNS = %d\n", lengthBufferNS);
                    HANDLE_API_ERR(cudaMalloc((void **)&GPUMPIDataSendNS, sizeBufferNS));
                    HANDLE_API_ERR(cudaMalloc((void **)&GPUMPIDataRecvNS, sizeBufferNS));
#ifndef CUDA_AWARE_MPI //! traditional MPI used with HostMPIDataSendNS and \
                           //! HostMPIDataRecvNS
                    HANDLE_API_ERR(cudaHostAlloc((void **)&HostMPIDataSendNS, sizeBufferNS, cudaHostAllocDefault));
                    HANDLE_API_ERR(cudaHostAlloc((void **)&HostMPIDataRecvNS, sizeBufferNS, cudaHostAllocDefault));
#endif
                    size_t sizeFaceIndex = nIFace * sizeof(int);
                    HANDLE_API_ERR(cudaMalloc((void **)&d_faceIndexForSend, sizeFaceIndex));
                    HANDLE_API_ERR(
                        cudaMemcpy(d_faceIndexForSend, h_faceIndexForSend, sizeFaceIndex, cudaMemcpyHostToDevice));
                    HANDLE_API_ERR(cudaMalloc((void **)&d_faceIndexForRecv, sizeFaceIndex));
                    HANDLE_API_ERR(
                        cudaMemcpy(d_faceIndexForRecv, h_faceIndexForRecv, sizeFaceIndex, cudaMemcpyHostToDevice));
                    delete[] h_faceIndexForSend;
                    delete[] h_faceIndexForRecv;

//! just for macro ISENDIRECVMPI, the original MPI_Isend and MPI_Irecv
#ifdef ORGISENDIRECVORDER
                    isSendOrRecv  = new int[nTotalRequestNS];
                    sendRecvNgbID = new int[nTotalRequestNS];
                    int nTotalNgb = numNgbZones[0];
                    int offset    = 0;
                    for (int iZone = 0; iZone < nTotalGloablZones; iZone++)
                    {
                        ZoneNeighbor *zoneNeighbor     = zoneConnectivity->GetZoneNeighbor(iZone);
                        std::size_t   numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();
                        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++iNeighbor)
                        {
                            int neighborZoneIndex  = zoneNeighbor->GetZoneIndexOfNeighbor(iNeighbor);
                            int sendProcessorID    = GetZoneProcessorIDSepMode(iZone);
                            int receiveProcessorID = GetZoneProcessorIDSepMode(neighborZoneIndex);
                            if (currentProcessorID != sendProcessorID && currentProcessorID != receiveProcessorID)
                            {
                                continue;
                            }
                            if (currentProcessorID == sendProcessorID)
                            {
                                isSendOrRecv[offset]  = 1;
                                sendRecvNgbID[offset] = iNeighbor;
                                offset++;
                            }
                            else if (currentProcessorID == receiveProcessorID)
                            {
                                isSendOrRecv[offset] = 0;
                                int jNeighbor        = 999;
                                //! jNeighbor is the local neighbor zone ID of zone in
                                //! currentProcessor
                                for (int ngbZoneID = 0; ngbZoneID < nTotalNgb; ngbZoneID++)
                                {
                                    if (recvGlobalZone[ngbZoneID] == iZone) jNeighbor = ngbZoneID;
                                }
                                sendRecvNgbID[offset] = jNeighbor;
                                offset++;
                            }
                        }
                    }
                    /*
      //!Test sendRecvNgbID and isSendOrRecv
      printf("globalRank = %d, offset = %d\n", globalRank, offset);
      for (int mpiID = 0; mpiID < nTotalRequestNS; mpiID++){
              int ngbID = sendRecvNgbID[mpiID];
              printf("globalRank = %d, mpiID = %d, isSendOrRecv = %d,
      sendRecvNgbID = %d, recvProcess = %d, recvGlobalZone = %d\n", globalRank,
      mpiID, isSendOrRecv[mpiID], ngbID, recvProcess[ngbID],
      recvGlobalZone[ngbID]);
      }
      //!Test end
      */
#endif //! end of ORGISENDIRECVORDER
#endif //! end if NEWMPI                                                                        \
    //! For test partition
                    char hostname[32];

                    if (gethostname(hostname, sizeof(hostname)))
                    {
                        printf("Error: gethostname calling error\n");
                        exit(1);
                    }
#ifndef NEWMPI
                    printf("globalRank = %d, localRank = %d, nTotalCell = %d, nTotalNode = %d, "
                           "nBoundFace = %d, nIFace = %d, hostname=%s\n",
                           globalRank, localRank, nTotalCell, nTotalNode, nBoundFace, nIFace, hostname);
#else
                    printf("globalRank = %d, localRank = %d, nTotalCell = %d, nTotalNode = %d, "
                           "nBoundFace = %d, nIFace = %d, numNgbZones = %d, hostname=%s\n",
                           globalRank, localRank, nTotalCell, nTotalNode, nBoundFace, nIFace, numNgbZones[0], hostname);
#endif
                    //! Test end
                }
                //! allocate device memory for interpoint transfer storage
                InterpointInformation *ipointinfo = grid->GetInterpointInfo();
                if (ipointinfo)
                {
                    d_interpointValueExt = true;
                    int nIPoint          = ipointinfo->GetNumberOfInterpoints();
#ifdef MPIBUFFERSMALL
                    GPUInterpointTransferStorageNS(nIPoint, nEqn);
#endif

#ifdef MPIBUFFERLARGE
                    InterpointFields *interpointFields = grid->GetInterpointFields();
                    GPUTransferInterpointBufferNS(interpointFields, solverNS, nIPoint);
#ifdef HOSTDEVOVERLAP
                    //! It should be noted that steams for H2D and D2H overlap are used for
                    //! both NS branch and Turb branch
                    HANDLE_API_ERR(cudaStreamCreateWithFlags(&dataTransferDToHInterpoint, cudaStreamNonBlocking));
                    HANDLE_API_ERR(cudaStreamCreateWithFlags(&dataTransferHToDInterpoint, cudaStreamNonBlocking));

                    //! Events just used for recording data transfer condition in NS branch
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDNSqNode_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDNStNode_done, cudaEventDisableTiming));

                    GPUInterpointTransferStorageNS(nIPoint, nEqn);
#endif
#endif
#ifdef NEWMPI
                    //! It is similar to treatment in NS
                    //! Get total number of zones
                    int nTotalGloablZones = GetNumberofGlobalZones();
                    //! set offsetTagInterpointNS for tags used in Isend and Irecv
                    offsetTagInterpointNS = offsetTagInterfaceNS + nTotalGloablZones;
                    //! printf("nTotalGloablZones = %d\n", nTotalGloablZones);
                    //! Get current process ID
                    int currentProcessorID = GetCurrentProcessorID();
                    //! One process owns only one zone
                    sendGlobalZoneForPoint = new int[1]; //! only consider one zone in one process
                    //! The current zone owns its neighbor zones
                    numNgbZonesForPoint = new int[1];
                    //! malloc memory for pointIndexForSend/Recv
                    int *h_pointIndexForSend = new int[nIPoint];
                    int *h_pointIndexForRecv = new int[nIPoint];
                    //! loop for all of zones
                    for (int iZone = 0; iZone < nTotalGloablZones; iZone++)
                    {
                        //! ZoneNeighbor *zoneNeighborForPoint =
                        //! zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone); std::size_t
                        //! numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor(); Get
                        //! the process ID of a zone
                        int sendProcessID = GetZoneProcessorIDSepMode(iZone);
                        //! Only consider send process is equal to current process
                        if (sendProcessID == currentProcessorID)
                        {
                            //! set sendGlobalZoneForPoint, in fact, it is the zone label for the
                            //! current process
                            sendGlobalZoneForPoint[0] = iZone;
                            //! Get zone connectivity of the zone in the current process
                            ZoneNeighbor *zoneNeighborForPoint =
                                zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);
                            //! Get neighbor zone number of the zone in the current process
                            std::size_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();
                            //! malloc memory for recvGlobalZoneForPoint and recvProcessForPoint
                            recvGlobalZoneForPoint = new int[numberOfNeighbor];
                            recvProcessForPoint    = new int[numberOfNeighbor];
                            //! malloc memory for offet of GPUMPIDataInterpointSendNS and
                            //! GPUMPIDataInterpointRecvNS
                            offsetMPIDataInterpointNS = new int[numberOfNeighbor];
                            //! malloc memory for offset of pointIndexForSend and pointIndexForRecv
                            offsetPointIndex       = new int[numberOfNeighbor];
                            int sumNgbLengthBuffer = 0;
                            int sumNgbPoints       = 0;
                            int countPoint         = 0;
                            //! set neighbor zone number of the zone in the current process
                            numNgbZonesForPoint[0] = numberOfNeighbor;
                            //! malloc memory for length of buffer of point for neighbor zones
                            ngbInterpointBufferLengthNS = new int[numberOfNeighbor];
                            //! send and recv for one neighbor zone
                            nTotalRequestInterpointNS = 2 * numberOfNeighbor; //! send and receive
                            //! malloc memory for requestContainer for point
                            requestContainerInterpointNS = new MPI_Request[nTotalRequestInterpointNS];
                            //! loop for all of neighbor zones of the current zone
                            for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++iNeighbor)
                            {
                                //! neighborZoneIndexForPoint is global zone label of neighbor zone
                                int neighborZoneIndexForPoint = zoneNeighborForPoint->GetZoneIndexOfNeighbor(iNeighbor);
                                //! ngbIDForPoint is a local zone label from 0 to numberOfNeighbor
                                int ngbIDForPoint = ipointinfo->FindIthNeighbor(neighborZoneIndexForPoint);
                                //! just check whether or not iNeighbor is equal to ngbID
                                if (ngbIDForPoint != iNeighbor)
                                {
                                    printf("ngbIDForPoint != iNeighbor, ngbIDForPoint = %d, iNeighbor = "
                                           "%d, iZone = %d, neighborZoneIndexForPoint = %d\n",
                                           ngbIDForPoint, iNeighbor, iZone, neighborZoneIndexForPoint);
                                    exit(1);
                                }
                                //! else printf("ngbIDForPoint is equal to iNeighbor!\n");
                                //! Get number of points shared by zone in the current process
                                //! (iZone) and its neighbor zone (jZone)
                                int nIPointOfNeighbor = ipointinfo->GetNumberOfInterpointsForNeighbor(ngbIDForPoint);
                                //! Get pointIndexForSend/Recv of jZone, it is reflection between
                                //! local point label in jZone and global point label in iZone
                                int *pointIndexForSend = ipointinfo->GetPointIndexForSend(ngbIDForPoint);
                                int *pointIndexForRecv = ipointinfo->GetPointIndexForRecv(ngbIDForPoint);
                                //! Get process ID of jZone
                                int receiveProcessID = GetZoneProcessorIDSepMode(neighborZoneIndexForPoint);
                                //! Store zone ID and process ID of jZone
                                recvGlobalZoneForPoint[iNeighbor] = neighborZoneIndexForPoint;
                                recvProcessForPoint[iNeighbor]    = receiveProcessID;
                                //! store offset in GPUMPIDataInterpointSend/RecvNS and
                                //! pointIndexForSend/Recv
                                offsetMPIDataInterpointNS[iNeighbor] = sumNgbLengthBuffer;
                                offsetPointIndex[iNeighbor]          = sumNgbPoints;
                                //! update sumNgbPoints with nIPointOfNeighbor in the neighbor zone
                                sumNgbPoints += nIPointOfNeighbor;
                                //! Init ngbInterpointBufferLengthNS as zero
                                ngbInterpointBufferLengthNS[iNeighbor] = 0;
                                //! loop for all of variables in nameNSInterpointBuffer
                                for (int nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
                                {
                                    //! Get dimension of the variable
                                    int dimData = dimInterpointBufferNS[nameID];
                                    //! add length of the variable in the sumNgbLengthBuffer

                                    sumNgbLengthBuffer += dimData * nIPointOfNeighbor;
                                    //! add length of the varialbe in ngbInterpointBufferLengthNS,
                                    //! noting the difference between sumNgbLengthBuffer and
                                    //! ngbInterpointBufferLengthNS
                                    ngbInterpointBufferLengthNS[iNeighbor] += dimData * nIPointOfNeighbor;
                                } //! end nameID of buffer
                                //! loop all of points shared by iZone and jZone
                                for (int jPointOfNeighbor = 0; jPointOfNeighbor < nIPointOfNeighbor; ++jPointOfNeighbor)
                                {
                                    //! Get global point label in iZone by local point in jZone with
                                    //! the help of pointIndexForSend and pointIndexForRecv
                                    int sendIPoint = pointIndexForSend[jPointOfNeighbor];
                                    int recvIPoint = pointIndexForRecv[jPointOfNeighbor];
                                    //! just for check
                                    if (sendIPoint != recvIPoint)
                                    {
                                        printf("iZone = %d, ngbIDForPoint = %d, countPoint = %d\n", iZone,
                                               ngbIDForPoint, countPoint);
                                    }
                                    //! Set global point label in a long 1D variable
                                    //! h_pointIndexForSend/Recv
                                    h_pointIndexForSend[countPoint] = sendIPoint;
                                    h_pointIndexForRecv[countPoint] = recvIPoint;
                                    //! update offset in h_pointIndexForSend/Recv
                                    countPoint++;
                                }
                            } //! end for iNeighbor
                            //! Test sumNgbLengthBuffer and sumNgbFaces. It should be noted that
                            //! sumNgbLengthBuffer is not equal to lengthBufferNS, becasue of the
                            //! varaible limit printf("sumNgbLengthBuffer = %d, lengthBufferNS =
                            //! %d\n", sumNgbLengthBuffer, lengthBufferNS); printf("sumNgbFaces =
                            //! %d, nIFace = %d\n", sumNgbFaces, nIFace); Test countFace and
                            //! nIFace, they should be equal. printf("countFace = %d, nIFace =
                            //! %d\n", countFace, nIFace);
                        } //! end for sendProcessID == currentProcessorID
                    } //! end of for iZone
                    /*
      //!test zone connectivity
      int iZone = sendGlobalZone[0];
      int nTotalNgbZones = numNgbZones[0];
      printf("Zone %d owns %d neighbor zones\n", iZone, nTotalNgbZones);
      for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++){
              printf("Zone %d owns %d neighbor zone in process %d\n", iZone,
      recvGlobalZone[ngbID], recvProcess[ngbID]);
      }
      */

                    /*
      //!Added for test
      if (globalRank == 0){
              for(int nameID = 0; nameID < nameNSInterpointBuffer->size();
      nameID++){ int dimData = dimInterpointBufferNS[nameID]; const string &
      name = (*nameNSInterpointBuffer)[nameID]; cout<<"nameID="<<nameID<<"
      name="<<name<<" dim="<<dimData<<endl; } //!end nameID of buffer
      }

      int iZone = sendGlobalZoneForPoint[0]; //!only consider one zone in one
      process ZoneNeighbor *zoneNeighborForPoint =
      zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone); int
      numberOfNeighbor = numNgbZonesForPoint[0]; for (int ngbID = 0; ngbID <
      numberOfNeighbor; ngbID++){ int ngbZoneID = recvGlobalZoneForPoint[ngbID];
              int neighborZoneIndexForPoint =
      zoneNeighborForPoint->GetZoneIndexOfNeighbor(ngbID); int iNeighbor =
      ipointinfo->FindIthNeighbor(neighborZoneIndexForPoint); int
      nIPointOfNeighbor =
      ipointinfo->GetNumberOfInterpointsForNeighbor(iNeighbor); int ngbProcessID
      = recvProcessForPoint[ngbID]; int offsetFace = offsetPointIndex[ngbID];
      int offsetNgb = offsetMPIDataInterpointNS[ngbID];
              printf("globalRank=%d,localRank=%d,iZone=%d,ngbID=%d,ngbZoneID=%d,ngbProcessID=%d,offsetFace=%d,offsetNgb=%d,
      nIPointOfNeighbor=%d,
      nIPoint=%d\n",globalRank,localRank,iZone,ngbID,ngbZoneID,ngbProcessID,offsetFace,offsetNgb,
      nIPointOfNeighbor, nIPoint);
      }
      //!test end
      */
                    size_t sizeBufferNS = lengthInterpointBufferNS * sizeof(double);
                    //! printf("lengthBufferNS = %d\n", lengthBufferNS);
                    //! malloc GPUMPIDataSend/RecvInterpointNS
                    HANDLE_API_ERR(cudaMalloc((void **)&GPUMPIDataSendInterpointNS, sizeBufferNS));
                    HANDLE_API_ERR(cudaMalloc((void **)&GPUMPIDataRecvInterpointNS, sizeBufferNS));
#ifndef CUDA_AWARE_MPI //! traditional MPI used with HostMPIDataSendNS and \
                           //! HostMPIDataRecvNS
                    //! malloc HostMPIDataSend/RecvInterpointNS for get device varialbe
                    HANDLE_API_ERR(
                        cudaHostAlloc((void **)&HostMPIDataSendInterpointNS, sizeBufferNS, cudaHostAllocDefault));
                    HANDLE_API_ERR(
                        cudaHostAlloc((void **)&HostMPIDataRecvInterpointNS, sizeBufferNS, cudaHostAllocDefault));
#endif
                    size_t sizeFaceIndex = nIPoint * sizeof(int);
                    HANDLE_API_ERR(cudaMalloc((void **)&d_pointIndexForSend, sizeFaceIndex));
                    HANDLE_API_ERR(
                        cudaMemcpy(d_pointIndexForSend, h_pointIndexForSend, sizeFaceIndex, cudaMemcpyHostToDevice));
                    HANDLE_API_ERR(cudaMalloc((void **)&d_pointIndexForRecv, sizeFaceIndex));
                    HANDLE_API_ERR(
                        cudaMemcpy(d_pointIndexForRecv, h_pointIndexForRecv, sizeFaceIndex, cudaMemcpyHostToDevice));
                    delete[] h_pointIndexForSend;
                    delete[] h_pointIndexForRecv;

//! just for macro ISENDIRECVMPI, the original MPI_Isend and MPI_Irecv
#ifdef ORGISENDIRECVORDER
                    //! malloc isSendOrRecvForPoint for storing whether send or recv operation;
                    //! malloc sendRecvNgbIDForPoint for storing send to or recv from local
                    //! neighbor zone ID (from 0 to numNgbZonesForPoint[0])
                    isSendOrRecvForPoint  = new int[nTotalRequestInterpointNS];
                    sendRecvNgbIDForPoint = new int[nTotalRequestInterpointNS];
                    int nTotalNgb         = numNgbZonesForPoint[0];
                    int offset            = 0;
                    //! loop all of Zones
                    for (int iZone = 0; iZone < nTotalGloablZones; iZone++)
                    {
                        //! Get zone connectivity info of iZone
                        ZoneNeighbor *zoneNeighborForPoint = zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);
                        //! Get neighbor zone number of iZone
                        std::size_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();
                        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++iNeighbor)
                        {
                            //! Get global zone ID of iZone's neighbor due to point
                            int neighborZoneIndexForPoint = zoneNeighborForPoint->GetZoneIndexOfNeighbor(iNeighbor);
                            //! Get send process and recv process
                            int sendProcessorID    = GetZoneProcessorIDSepMode(iZone);
                            int receiveProcessorID = GetZoneProcessorIDSepMode(neighborZoneIndexForPoint);
                            if (currentProcessorID != sendProcessorID && currentProcessorID != receiveProcessorID)
                            {
                                continue;
                            }
                            if (currentProcessorID == sendProcessorID)
                            {
                                //! For send process, set isSendOrRecvForPoint 1 and
                                //! sendRecvNgbIDForPoint the local zone label iNeighbor
                                isSendOrRecvForPoint[offset]  = 1;
                                sendRecvNgbIDForPoint[offset] = iNeighbor;
                                offset++;
                            }
                            else if (currentProcessorID == receiveProcessorID)
                            {
                                //! For recv process, set isSendOrRecvForPoint 0
                                isSendOrRecvForPoint[offset] = 0;
                                int jNeighbor                = 999;
                                //! jNeighbor is the local neighbor zone ID of zone in
                                //! currentProcessor loop all of neighbor zones of iZone, find the
                                //! local zone label of jZone
                                for (int ngbZoneID = 0; ngbZoneID < nTotalNgb; ngbZoneID++)
                                {
                                    if (recvGlobalZoneForPoint[ngbZoneID] == iZone) jNeighbor = ngbZoneID;
                                }
                                sendRecvNgbIDForPoint[offset] = jNeighbor;
                                offset++;
                            }
                        }
                    }
                    /*
      //!Test sendRecvNgbID and isSendOrRecv
      printf("globalRank = %d, offset = %d\n", globalRank, offset);
      for (int mpiID = 0; mpiID < nTotalRequestInterpointNS; mpiID++){
              int ngbID = sendRecvNgbIDForPoint[mpiID];
              printf("globalRank = %d, mpiID = %d, isSendOrRecv = %d,
      sendRecvNgbID = %d, recvProcess = %d, recvGlobalZone = %d\n", globalRank,
      mpiID, isSendOrRecvForPoint[mpiID], ngbID, recvProcessForPoint[ngbID],
      recvGlobalZoneForPoint[ngbID]);
      }
      //!Test end
      */
#endif //! end of ORGISENDIRECVORDER
#endif //! end if NEWMPI
                }
                //! for tScheme=1, LU-SGS scheme
                int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
                if (tscheme == LU_SGS)
                {
                    printf("The scheme for GPU is LU-SGS in NS equations\n");
                    GPULUSGSDqAlloc(nTotal, nEqn);
                }

                //! for unsteady case
                int iunsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

                if (iunsteady != 0)
                {
                    GPUUnsteadyNSAlloc(nTotal, nEqn);
                }

            } //! end for (nLocalZones)
        } //! end function

        void GPUFlowVarTurbMemAlloc(const int solverTurb)
        {
            using namespace PHMPI;
            using namespace PHSPACE;
//! just for test
#ifdef CUDAUNITTEST
#ifdef UNITTESTOUTPUT
            printf("program is running in GPUFlowVarTurbMemAlloc\n");
#endif
#endif
            int nLocalZones = GetNumberofLocalZones();
            for (int izone = 0; izone < nLocalZones; izone++)
            {
                int zoneID = GetLocalZoneIDToGlobalZoneID(izone);
                //! Zone *zone = GetZone(zoneID);
                int level = 0;
                //! Grid *grid = GetGrid(zoneID, level);
                UnstructGrid *grid       = UnstructGridCast(GetGrid(zoneID, level));
                const int     nTotalCell = grid->GetNTotalCell();
                const int     nBoundFace = grid->GetNBoundFace();
                const int     nTotalFace = grid->GetNTotalFace();
                int           nTotal     = nTotalCell + nBoundFace;

                int nEqn_turb;
                GlobalDataBase::GetData("n_turb", &nEqn_turb, PHINT, 1);
                TurbEqnNumberCopy(nEqn_turb);

                //! memory alloc for turb parts: q, ql, qr
                GPUFieldQTurbProxyMemAlloc(nEqn_turb, nTotal);
                GPUQLRTurbMemAlloc(nEqn_turb, d_SEG_LEN);
                //! memory alloc for flux_turb flux_sub_turb
                d_flux_turb     = GPUFluxTurbAlloc(nEqn_turb, d_SEG_LEN);
                d_flux_sub_turb = GPUFluxTurbAllocSub(nEqn_turb, d_SEG_LEN);
                d_res_turb      = GPUResTurbAlloc(nEqn_turb, nTotal);
                GPUTurbFaceVarMemAlloc(nEqn_turb, d_SEG_LEN);
                GPUParametersMemAllocCopy(nEqn_turb);
#ifdef CUDAUNITTEST
                TestGPUParametersMemAllocCopy(nEqn_turb);
#endif
                GPUSpecTurbMemAlloc(nEqn_turb, nTotal);
                GPUMatTurbLRMemAlloc(nEqn_turb, nTotalFace);
                GPUTurbRhsMemAlloc(nEqn_turb, nTotal);
                GPUTurbTempProxyMemAlloc(nEqn_turb, nTotal);
                GPUTurbUpdateFlowNegPosMemAlloc();
                GPUVistMMalloc();
                int nTotalNode = grid->GetNTotalNode();
                GPUComputeNodeQTurbValueMemAlloc(nTotalNode, nEqn_turb);
                //! alloc device memory for interpolation points
                InterpointInformation *interPointInfor = grid->GetInterpointInfo();
                if (interPointInfor)
                {
                    int nIPoint = interPointInfor->GetNumberOfInterpoints();
                    GPUTurbInterPointMemAlloc(nIPoint, nEqn_turb);
                }
                //! allocate device memroy for interface transfer storage
                InterfaceInfo *iinfo = grid->GetInterfaceInfo();
                //! only for multi-process
                if (iinfo)
                {
                    int nIFace = iinfo->GetNIFace();
#ifdef MPIBUFFERSMALL
                    GPUInterfaceTransferStorageTurb(nIFace, nEqn_turb);
#endif
#ifdef MPIBUFFERLARGE
                    InterfaceFields *interfaceFields = grid->GetInterfaceFields();
                    GPUInterfaceTransferBufferTurb(interfaceFields, solverTurb, nIFace);
#ifdef HOSTDEVOVERLAP
                    //! small send buffer is required
                    GPUInterfaceTransferStorageTurb(nIFace, nEqn_turb);

                    //! Events just used for recording data transfer condition in NS branch
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadDToDTurb_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDTurbq_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDTurbdqdx_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDTurbdqdy_done, cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDTurbdqdz_done, cudaEventDisableTiming));

#endif
#endif
#ifdef NEWMPI
                    //! Getting zone connectivity information
                    int nTotalGloablZones = GetNumberofGlobalZones();
                    //! set offsetTagInterfaceTurb for tags used in Isend and Irecv
                    offsetTagInterfaceTurb = offsetTagInterfaceNS + 2 * nTotalGloablZones;
                    //! printf("nTotalGloablZones = %d\n", nTotalGloablZones);
                    int currentProcessorID = GetCurrentProcessorID();
                    for (int iZone = 0; iZone < nTotalGloablZones; iZone++)
                    {
                        ZoneNeighbor *zoneNeighbor     = zoneConnectivity->GetZoneNeighbor(iZone);
                        std::size_t   numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();
                        int           sendProcessID    = GetZoneProcessorIDSepMode(iZone);
                        if (sendProcessID == currentProcessorID)
                        {
                            ZoneNeighbor *zoneNeighbor     = zoneConnectivity->GetZoneNeighbor(iZone);
                            std::size_t   numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();
                            offsetMPIDataTurb              = new int[numberOfNeighbor];
                            int sumNgbLengthBuffer         = 0;
                            int countFace                  = 0;
                            ngbBufferLengthTurb            = new int[numberOfNeighbor];
                            //! MPI_Request * requestContainerNS;
                            nTotalRequestTurb    = 2 * numberOfNeighbor; //! send and receive
                            requestContainerTurb = new MPI_Request[nTotalRequestTurb];
                            for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++iNeighbor)
                            {
                                int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(iNeighbor);
                                //! the ngbID'th neighbor zone of iZone
                                int ngbID = iinfo->FindIthNeighbor(neighborZoneIndex);
                                //! just check whether or not iNeighbor is equal to ngbID
                                if (ngbID != iNeighbor)
                                {
                                    printf("ngbID != iNeighbor, ngbID = %d, iNeighbor = %d, iZone = %d, "
                                           "neighborZoneIndex = %d\n",
                                           ngbID, iNeighbor, iZone, neighborZoneIndex);
                                    exit(1);
                                }
                                //! else printf("ngbID is equal to iNeighbor!\n");
                                int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(ngbID);

                                //! calculate offsetMPIDataTurb
                                offsetMPIDataTurb[iNeighbor]   = sumNgbLengthBuffer;
                                ngbBufferLengthTurb[iNeighbor] = 0;
                                for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                                {
                                    int dimData = dimBufferTurb[nameID];

                                    sumNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                                    ngbBufferLengthTurb[iNeighbor] += dimData * nIFaceOfNeighbor;
                                } //! end nameID of buffer
                            } //! end for iNeighbor
                            //! Test sumNgbLengthBuffer and sumNgbFaces. It should be noted that
                            //! sumNgbLengthBuffer is not equal to lengthBufferTurb, becasue of the
                            //! varaible limit printf("For turbulence sumNgbLengthBuffer = %d,
                            //! lengthBufferTurb = %d\n", sumNgbLengthBuffer, lengthBufferTurb);
                        } //! end for sendProcessID == currentProcessorID
                    } //! end of for iZone
#ifdef MPIOVERLAP

                    requestContainerTurbFirst   = new MPI_Request[nTotalRequestTurb];
                    requestContainerTurbRemains = new MPI_Request[nTotalRequestTurb];
                    //! compute ngbBufferLengthNSFirst and ngbBufferLengthNSRemains
                    int numberOfNeighbor       = numNgbZones[0];
                    ngbBufferLengthTurbFirst   = new int[numberOfNeighbor];
                    ngbBufferLengthTurbRemains = new int[numberOfNeighbor];
                    int           iZone        = sendGlobalZone[0]; //! only consider one zone in one process
                    ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
                    for (int ngbID = 0; ngbID < numberOfNeighbor; ngbID++)
                    {
                        int ngbZoneID                   = recvGlobalZone[ngbID];
                        int neighborZoneIndex           = zoneNeighbor->GetZoneIndexOfNeighbor(ngbID);
                        int iNeighbor                   = iinfo->FindIthNeighbor(neighborZoneIndex);
                        int nIFaceOfNeighbor            = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                        ngbBufferLengthTurbFirst[ngbID] = 0;

                        for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                        {
                            const string &name = (*nameTurbBuffer)[nameID];
                            //! just for test
                            //! if ((ngbID == 0)&&(globalRank==0)) cout<<"nameID ="<<nameID<<",
                            //! name ="<<name<<endl;
                            if (name == "turb::q")
                            {
                                int dimData = dimBufferTurb[nameID];
                                ngbBufferLengthTurbFirst[ngbID] += dimData * nIFaceOfNeighbor;
                                //! if ((ngbID == 0)&&(globalRank==0)) printf("nameID = %d,
                                //! ngbBufferLengthTurbFirst = %d\n", nameID,
                                //! ngbBufferLengthTurbFirst[ngbID]);
                            }
                        }

                        ngbBufferLengthTurbRemains[ngbID] =
                            ngbBufferLengthTurb[ngbID] - ngbBufferLengthTurbFirst[ngbID];
                        //! just for test
                        //! if ((ngbID == 0)&&(globalRank==0)) printf("In turbulence MPI
                        //! Interface buffer: First = %d, remains = %d, total = %d\n",
                        //! ngbBufferLengthTurbFirst[ngbID], ngbBufferLengthTurbRemains[ngbID],
                        //! ngbBufferLengthTurb[ngbID]);
                    }
#endif

#ifdef GPUBUFFERSMALL
                    //! Create
                    //! Get total number of neighbor zone of iZone
                    int nTotalNgbZones   = numNgbZones[0];
                    GPUDataSendTurb_q    = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataSendTurb_dqdx = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataSendTurb_dqdy = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataSendTurb_dqdz = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));

                    GPUDataRecvTurb_q    = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataRecvTurb_dqdx = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataRecvTurb_dqdy = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    GPUDataRecvTurb_dqdz = (RFloat **)malloc(nTotalNgbZones * sizeof(RFloat *));
                    //! loop all of neighbor zones for Parcel send data
                    for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
                    {
                        //! get neighbor zone's global index (jZone)
                        int ngbZoneID = recvGlobalZone[ngbID];
                        //! get jZone 's local index of iZone (from 0 to nTotalNgbZones), in face
                        //! ngbID=iNeighbor
                        int iNeighbor = iinfo->FindIthNeighbor(ngbZoneID);
                        //! get total faces of jZone's interface
                        int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                        for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                        {
                            const string &name = (*nameTurbBuffer)[nameID];
                            if (name == "turb::q")
                            {
                                int dimData = dimBufferTurb[nameID];
                                HANDLE_API_ERR(
                                    cudaMalloc((void **)&GPUDataSendTurb_q[ngbID], nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(
                                    cudaMalloc((void **)&GPUDataRecvTurb_q[ngbID], nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                            else if (name == "turb::dqdx")
                            {
                                int dimData = dimBufferTurb[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendTurb_dqdx[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvTurb_dqdx[ngbID],
                                                          dimData * nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                            else if (name == "turb::dqdy")
                            {
                                int dimData = dimBufferTurb[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendTurb_dqdy[ngbID],
                                                          nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvTurb_dqdy[ngbID],
                                                          nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                            else if (name == "turb::dqdz")
                            {
                                int dimData = dimBufferTurb[nameID];
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataSendTurb_dqdz[ngbID],
                                                          nIFaceOfNeighbor * sizeof(RFloat)));
                                HANDLE_API_ERR(cudaMalloc((void **)&GPUDataRecvTurb_dqdz[ngbID],
                                                          nIFaceOfNeighbor * sizeof(RFloat)));
                            }
                        } //! end for nameID

                    } //! end for ngbID

#ifdef HDCOVERLAP
                    HANDLE_API_ERR(cudaStreamCreateWithFlags(&dataStreamTurb, cudaStreamNonBlocking));

                    //! HANDLE_API_ERR(cudaEventCreateWithFlags(&transferDone,
                    //! cudaEventDisableTiming));
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&compressDoneTurb, cudaEventDisableTiming));
                    //! compute offsetNgbLBTurb that is global offsetNgbLengthBuffer
                    offsetNgbLBTurb = new int[nTotalNgbZones * nameTurbBuffer->size()];
                    for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++)
                    {
                        int ngbZoneID             = recvGlobalZone[ngbID];
                        int neighborZoneIndex     = zoneNeighbor->GetZoneIndexOfNeighbor(ngbID);
                        int iNeighbor             = iinfo->FindIthNeighbor(neighborZoneIndex);
                        int nIFaceOfNeighbor      = iinfo->GetNIFaceOfNeighbor(iNeighbor);
                        int offsetNgbLengthBuffer = 0;

                        for (int nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
                        {
                            int dimData                                              = dimBufferTurb[nameID];
                            offsetNgbLBTurb[ngbID * nameTurbBuffer->size() + nameID] = offsetNgbLengthBuffer;
                            offsetNgbLengthBuffer += dimData * nIFaceOfNeighbor;
                        }
                    }
                    /*
      //!Test for offsetNgbLBTurb
      for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++){
              int ngbZoneID = recvGlobalZone[ngbID];
              int neighborZoneIndex =
      zoneNeighbor->GetZoneIndexOfNeighbor(ngbID); int iNeighbor =
      iinfo->FindIthNeighbor(neighborZoneIndex); int nIFaceOfNeighbor =
      iinfo->GetNIFaceOfNeighbor(iNeighbor); printf("globalRank = %d, ngbID =
      %d, nIFaceOfNeighbor =%d\n", globalRank, ngbID, nIFaceOfNeighbor); for(int
      nameID = 0; nameID < nameTurbBuffer->size(); nameID++){ printf("globalRank
      = %d, nameID = %d, ngbID = %d, offsetNgbLBTurb=%d\n", globalRank, nameID,
      ngbID, offsetNgbLBTurb[ngbID*nameTurbBuffer->size()+nameID]);
              }

      }
      //!Test end
      */

#endif //! end ifdef HDCOVERLAP
#endif //! end ifdef \
               //! Create    \
               //! Get total number of neighbor zone of iZone
                    /*
      //!Added for test
      if (globalRank == 0){
              for(int nameID = 0; nameID < nameTurbBuffer->size(); nameID++){
                      int dimData = dimBufferTurb[nameID];
                      const string & name = (*nameTurbBuffer)[nameID];
                      cout<<"nameID="<<nameID<<" name="<<name<<"
      dim="<<dimData<<endl; } //!end nameID of buffer
      }
      */

                    size_t sizeBufferTurb = lengthBufferTurb * sizeof(double);
                    //! printf("lengthBufferTurb = %d\n", lengthBufferTurb);
                    HANDLE_API_ERR(cudaMalloc((void **)&GPUMPIDataSendTurb, sizeBufferTurb));
                    HANDLE_API_ERR(cudaMalloc((void **)&GPUMPIDataRecvTurb, sizeBufferTurb));
#ifndef CUDA_AWARE_MPI //! traditional MPI used with HostMPIDataSendTurb and \
                           //! HostMPIDataRecvTurb
                    HANDLE_API_ERR(cudaHostAlloc((void **)&HostMPIDataSendTurb, sizeBufferTurb, cudaHostAllocDefault));
                    HANDLE_API_ERR(cudaHostAlloc((void **)&HostMPIDataRecvTurb, sizeBufferTurb, cudaHostAllocDefault));
#endif

#endif //! end if NEWMPI
                }
                InterpointInformation *ipointinfo = grid->GetInterpointInfo();
                if (ipointinfo)
                {
                    int nIPoint = ipointinfo->GetNumberOfInterpoints();
#ifdef MPIBUFFERSMALL
                    GPUInterpointTransferStorageTurb(nIPoint, nEqn_turb);
#endif
#ifdef MPIBUFFERLARGE
                    InterpointFields *interpointFields = grid->GetInterpointFields();
                    GPUTransferInterpointBufferTurb(interpointFields, solverTurb, nIPoint);
#ifdef HOSTDEVOVERLAP
                    //! Events just used for recording data transfer condition in NS branch
                    HANDLE_API_ERR(cudaEventCreateWithFlags(&downloadHToDTurbqNode_done, cudaEventDisableTiming));

                    GPUInterpointTransferStorageTurb(nIPoint, nEqn_turb);
#endif
#endif
#ifdef NEWMPI
                    //! It is similar to treatment in NS
                    //! Get total number of zones
                    int nTotalGloablZones = GetNumberofGlobalZones();
                    //! set offsetTagInterpointTurb for tags used in Isend and Irecv
                    offsetTagInterpointTurb = offsetTagInterfaceTurb + nTotalGloablZones;
                    //! Get current process ID
                    int currentProcessorID = GetCurrentProcessorID();
                    //! loop for all of zones
                    for (int iZone = 0; iZone < nTotalGloablZones; iZone++)
                    {
                        //! Get the process ID of a zone
                        int sendProcessID = GetZoneProcessorIDSepMode(iZone);
                        //! Only consider send process is equal to current process
                        if (sendProcessID == currentProcessorID)
                        {
                            //! Get zone connectivity of the zone in the current process
                            ZoneNeighbor *zoneNeighborForPoint =
                                zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);
                            //! Get neighbor zone number of the zone in the current process
                            std::size_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();
                            //! malloc memory for offet of GPUMPIDataInterpointSendNS and
                            //! GPUMPIDataInterpointRecvNS
                            offsetMPIDataInterpointTurb = new int[numberOfNeighbor];
                            int sumNgbLengthBuffer      = 0;
                            //! malloc memory for length of buffer of point for neighbor zones
                            ngbInterpointBufferLengthTurb = new int[numberOfNeighbor];
                            //! send and recv for one neighbor zone
                            nTotalRequestInterpointTurb = 2 * numberOfNeighbor; //! send and receive
                            //! malloc memory for requestContainer for point
                            requestContainerInterpointTurb = new MPI_Request[nTotalRequestInterpointTurb];
                            //! loop for all of neighbor zones of the current zone
                            for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++iNeighbor)
                            {
                                //! neighborZoneIndexForPoint is global zone label of neighbor zone
                                int neighborZoneIndexForPoint = zoneNeighborForPoint->GetZoneIndexOfNeighbor(iNeighbor);
                                //! ngbIDForPoint is a local zone label from 0 to numberOfNeighbor
                                int ngbIDForPoint = ipointinfo->FindIthNeighbor(neighborZoneIndexForPoint);
                                //! just check whether or not iNeighbor is equal to ngbID
                                if (ngbIDForPoint != iNeighbor)
                                {
                                    printf("ngbIDForPoint != iNeighbor, ngbIDForPoint = %d, iNeighbor = "
                                           "%d, iZone = %d, neighborZoneIndexForPoint = %d\n",
                                           ngbIDForPoint, iNeighbor, iZone, neighborZoneIndexForPoint);
                                    exit(1);
                                }
                                //! else printf("ngbIDForPoint is equal to iNeighbor!\n");
                                //! Get number of points shared by zone in the current process
                                //! (iZone) and its neighbor zone (jZone)
                                int nIPointOfNeighbor = ipointinfo->GetNumberOfInterpointsForNeighbor(ngbIDForPoint);
                                //! store offset in GPUMPIDataInterpointSend/RecvTurb and
                                //! pointIndexForSend/Recv
                                offsetMPIDataInterpointTurb[iNeighbor] = sumNgbLengthBuffer;
                                //! Init ngbInterpointBufferLengthTurb as zero
                                ngbInterpointBufferLengthTurb[iNeighbor] = 0;
                                //! loop for all of variables in nameTurbInterpointBuffer
                                for (int nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
                                {
                                    //! Get dimension of the variable
                                    int dimData = dimInterpointBufferTurb[nameID];
                                    //! add length of the variable in the sumNgbLengthBuffer

                                    sumNgbLengthBuffer += dimData * nIPointOfNeighbor;
                                    //! add length of the varialbe in ngbInterpointBufferLengthTurb,
                                    //! noting the difference between sumNgbLengthBuffer and
                                    //! ngbInterpointBufferLengthTurb
                                    ngbInterpointBufferLengthTurb[iNeighbor] += dimData * nIPointOfNeighbor;
                                } //! end nameID of buffer
                            } //! end for iNeighbor
                            //! Test sumNgbLengthBuffer and sumNgbFaces. It should be noted that
                            //! sumNgbLengthBuffer is not equal to lengthBufferNS, becasue of the
                            //! varaible limit printf("sumNgbLengthBuffer = %d, lengthBufferNS =
                            //! %d\n", sumNgbLengthBuffer, lengthBufferNS); printf("sumNgbFaces =
                            //! %d, nIFace = %d\n", sumNgbFaces, nIFace); Test countFace and
                            //! nIFace, they should be equal. printf("countFace = %d, nIFace =
                            //! %d\n", countFace, nIFace);
                        } //! end for sendProcessID == currentProcessorID
                    } //! end of for iZone
                    /*
      //!test zone connectivity
      int iZone = sendGlobalZone[0];
      int nTotalNgbZones = numNgbZones[0];
      printf("Zone %d owns %d neighbor zones\n", iZone, nTotalNgbZones);
      for (int ngbID = 0; ngbID < nTotalNgbZones; ngbID++){
              printf("Zone %d owns %d neighbor zone in process %d\n", iZone,
      recvGlobalZone[ngbID], recvProcess[ngbID]);
      }
      */
                    /*
      //!Added for test
      if (globalRank == 0){
              for(int nameID = 0; nameID < nameTurbInterpointBuffer->size();
      nameID++){ int dimData = dimInterpointBufferTurb[nameID]; const string &
      name = (*nameTurbInterpointBuffer)[nameID]; cout<<"nameID="<<nameID<<"
      name="<<name<<" dim="<<dimData<<endl; } //!end nameID of buffer
      }

      int iZone = sendGlobalZoneForPoint[0]; //!only consider one zone in one
      process ZoneNeighbor *zoneNeighborForPoint =
      zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone); int
      numberOfNeighbor = numNgbZonesForPoint[0]; for (int ngbID = 0; ngbID <
      numberOfNeighbor; ngbID++){ int ngbZoneID = recvGlobalZoneForPoint[ngbID];
              int neighborZoneIndexForPoint =
      zoneNeighborForPoint->GetZoneIndexOfNeighbor(ngbID); int iNeighbor =
      ipointinfo->FindIthNeighbor(neighborZoneIndexForPoint); int
      nIPointOfNeighbor =
      ipointinfo->GetNumberOfInterpointsForNeighbor(iNeighbor); int ngbProcessID
      = recvProcessForPoint[ngbID]; int offsetFace = offsetPointIndex[ngbID];
      int offsetNgb = offsetMPIDataInterpointTurb[ngbID]; printf("For Turb
      globalRank=%d,localRank=%d,iZone=%d,ngbID=%d,ngbZoneID=%d,ngbProcessID=%d,offsetFace=%d,offsetNgb=%d,
      nIPointOfNeighbor=%d,
      nIPoint=%d\n",globalRank,localRank,iZone,ngbID,ngbZoneID,ngbProcessID,offsetFace,offsetNgb,
      nIPointOfNeighbor, nIPoint);
      }
      //!test end
      */
                    size_t sizeBufferTurb = lengthInterpointBufferTurb * sizeof(double);
                    //! printf("lengthBufferTurb = %d\n", lengthBufferTurb);
                    //! malloc GPUMPIDataSend/RecvInterpointTurb
                    HANDLE_API_ERR(cudaMalloc((void **)&GPUMPIDataSendInterpointTurb, sizeBufferTurb));
                    HANDLE_API_ERR(cudaMalloc((void **)&GPUMPIDataRecvInterpointTurb, sizeBufferTurb));
#ifndef CUDA_AWARE_MPI //! traditional MPI used with HostMPIDataSendTurb and \
                           //! HostMPIDataRecvTurb
                    //! malloc HostMPIDataSend/RecvInterpointTurb for get device varialbe
                    HANDLE_API_ERR(
                        cudaHostAlloc((void **)&HostMPIDataSendInterpointTurb, sizeBufferTurb, cudaHostAllocDefault));
                    HANDLE_API_ERR(
                        cudaHostAlloc((void **)&HostMPIDataRecvInterpointTurb, sizeBufferTurb, cudaHostAllocDefault));
#endif
#endif //! end if NEWMPI
                }

                //! for tScheme=1, LU-SGS scheme
                int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
                if (tscheme == LU_SGS)
                {
                    printf("The scheme for GPU is LU-SGS in Turb equations\n");
                    GPULUSGSDqTurbAlloc(nTotal, nEqn_turb);
                }
                //! for unsteady case
                int iunsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

                if (iunsteady != 0)
                {
                    GPUUnsteadyTurbAlloc(nTotal, nEqn_turb);
                }
            }
        }

        void GPUUnsteadyTurbAlloc(const int nTotal, const int nEqn)
        {
            size_t sizeUnsteady = nTotal * nEqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_q_turb_unsteady_n1, sizeUnsteady));
            HANDLE_API_ERR(cudaMalloc((void **)&d_q_turb_unsteady_n2, sizeUnsteady));
            HANDLE_API_ERR(cudaMalloc((void **)&d_res_turb_unsteady_n1, sizeUnsteady));
            HANDLE_API_ERR(cudaMalloc((void **)&d_res_turb_unsteady_n2, sizeUnsteady));
            HANDLE_API_ERR(cudaMalloc((void **)&d_res_turb_unsteady_tmp, sizeUnsteady));
            size_t sizeSum = sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_sum1_turb_unsteady, sizeSum));
            HANDLE_API_ERR(cudaMalloc((void **)&d_sum2_turb_unsteady, sizeSum));
        }
        void GPUUnsteadyNSAlloc(const int nTotal, const int nEqn)
        {
            size_t sizeUnsteady = nTotal * nEqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_q_ns_unsteady_n1, sizeUnsteady));
            HANDLE_API_ERR(cudaMalloc((void **)&d_q_ns_unsteady_n2, sizeUnsteady));
            HANDLE_API_ERR(cudaMalloc((void **)&d_res_ns_unsteady_n1, sizeUnsteady));
            HANDLE_API_ERR(cudaMalloc((void **)&d_res_ns_unsteady_n2, sizeUnsteady));
            HANDLE_API_ERR(cudaMalloc((void **)&d_res_ns_unsteady_tmp, sizeUnsteady));
            size_t sizeSum = sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_sum1_ns_unsteady, sizeSum));
            HANDLE_API_ERR(cudaMalloc((void **)&d_sum2_ns_unsteady, sizeSum));
        }

        void GPULUSGSDqAlloc(const int nTotal, const int nEqn)
        {
            //!size_t sizeDq = nTotal * nEqn * sizeof(RFloat);
            //!HANDLE_API_ERR(cudaMalloc((void **)&d_dq_ns, sizeDq));
            size_t sizeSpec = nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_spec_ns, sizeSpec));
        }
        void GPULUSGSDqTurbAlloc(const int nTotal, const int nEqn)
        {
            size_t sizeDq = nTotal * nEqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_dq_turb, sizeDq));
        }
        void GPUMallocMinMaxTurbViscous(const int nTotalCell)
        {
            //! cudaMalloc for computing dmin and dmax in one GPU
            using namespace GPUNSSolverUnstruct;
            int    loopLen         = nTotalCell;
            int    gridSize        = 1;
            int    blockSize       = 1;
            int    regsPerThread   = 32;
            int    residentBlockNo = 4;
            size_t dsMemPerThread  = 2 * sizeof(RFloat); //! calculate by the intended shared memory use in the kernel
            size_t dsMemPerBlock;
            KernelLaunchPara((void *)GPUComputeMinMaxStep1, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                             dsMemPerBlock, gridSize, blockSize, GPUProp);

            //! calculate the gridSize with register and shared memory limitation,
            //! because it will be used as blockSize in the second reduction.
            regsPerThread = 26; //! register for GPUComputeMinMaxStep2
            int gridSize_tmp;
            KernelLaunchPara((void *)GPUComputeMinMaxStep2, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                             dsMemPerBlock, gridSize_tmp, gridSize, GPUProp);
            size_t bSize = gridSize * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_minLocalTurbVis, bSize));
            HANDLE_API_ERR(cudaMalloc((void **)&d_maxLocalTurbVis, bSize));
        }

        void GPUMallocMinMaxLocalTime(const int nTotalCell)
        {
            //! cudaMalloc for computing dmin and dmax in one GPU
            using namespace GPUNSSolverUnstruct;
            int    gridSize        = 1;
            int    blockSize       = 1;
            int    loopLen         = nTotalCell;
            int    regsPerThread   = 32;
            int    residentBlockNo = 4;
            size_t dsMemPerThread  = 2 * sizeof(RFloat); //! calculate by the intended shared memory use in the kernel
            size_t dsMemPerBlock;
            KernelLaunchPara((void *)GPUComputeMinMaxStep1, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                             dsMemPerBlock, gridSize, blockSize, GPUProp);

            //! calculate the gridSize with register and shared memory limitation,
            //! because it will be used as blockSize in the second reduction.
            regsPerThread = 26; //! register for GPUComputeMinMaxStep2
            int gridSize_tmp;
            KernelLaunchPara((void *)GPUComputeMinMaxStep2, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                             dsMemPerBlock, gridSize_tmp, gridSize, GPUProp);

            size_t bSize = gridSize * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_minLocalTime, bSize));
            HANDLE_API_ERR(cudaMalloc((void **)&d_maxLocalTime, bSize));
        }

        void EqnNumberCopy(int nl, int nchem, int neqn)
        {
            d_nl      = nl;
            d_nchem   = nchem;
            d_neqn_ns = neqn;
        }

        void TurbEqnNumberCopy(int nEqn_turb) { d_neqn_turb = nEqn_turb; }

        void GPUQ_NSAlloc(const int m, const int n)
        {
            size_t qSize = m * n * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_q_ns, qSize));
        }
        RFloat *GPUQ_NSAlloc1(const int m, const int n)
        {
            RFloat *q_ns;
            size_t  qSize = m * n * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&q_ns, qSize));
            return q_ns;
        }
        //! ql_ns qr_ns variables Alloc and Copy
        void GPUQl_NSAlloc(const int m, const int n)
        {
            size_t qlSize = m * n * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_ql_ns, qlSize));
        }
        void GPUQr_NSAlloc(const int m, const int n)
        {
            size_t qrSize = m * n * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_qr_ns, qrSize));
        }

        void GPUFieldQProxyMemAlloc(int nEqn, int nTotal)
        {
            size_t sizeFieldVars = nEqn * nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_q_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_q_proxy_tmp, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dqdx_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dqdy_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dqdz_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_q_old, sizeFieldVars));
        }

        void GPUFieldTProxyMemAlloc(int nEqn, int nTotal)
        {
            size_t sizeFieldVars = nEqn * nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_t_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dtdx_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dtdy_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dtdz_proxy, sizeFieldVars));
        }

        void GPUFieldQTurbProxyMemAlloc(int nEqn, int nTotal)
        {
            size_t sizeFieldVars = nEqn * nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_q_turb_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dq_turbdx_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dq_turbdy_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dq_turbdz_proxy, sizeFieldVars));
        }

        void GPUFieldVelocityProxyMemAlloc(int nEqn, int nTotal)
        {
            size_t sizeFieldVars = nEqn * nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_vel_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dveldx_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dveldy_proxy, sizeFieldVars));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dveldz_proxy, sizeFieldVars));
        }

        void GPUNodeVarsMemAlloc(int nTotalNode)
        {
            size_t sizeNodeInt = nTotalNode * sizeof(int);
            HANDLE_API_ERR(cudaMalloc((void **)&d_n_count_tmp, sizeNodeInt));
            size_t sizeNodeRFloat = nTotalNode * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_q_n_tmp, sizeNodeRFloat));
            //! size_t sizeNodeDouble = nTotalNode*sizeof(double);
            //! HANDLE_API_ERR(cudaMalloc((void**) &d_q_n_double, sizeNodeDouble));
        }

        void GPUFieldQBoundaryGradientMemAlloc(int nEqn, int nBoundFace)
        {
            size_t sizeFaceRFloat = nEqn * nBoundFace * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_bdqx, sizeFaceRFloat));
            HANDLE_API_ERR(cudaMalloc((void **)&d_bdqy, sizeFaceRFloat));
            HANDLE_API_ERR(cudaMalloc((void **)&d_bdqz, sizeFaceRFloat));
        }

        //! Gama GamaL GamaR alloc
        RFloat *GPUGamaAlloc(const int nTotal)
        {
            RFloat *d_gama;
            size_t  gamaSize = nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_gama, gamaSize));
            return d_gama;
        }
        RFloat *GPUGamaLAlloc(const int len_gamaL)
        {
            RFloat *d_gamaL;
            size_t  gamaSize = len_gamaL * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_gamaL, gamaSize));
            return d_gamaL;
        }
        RFloat *GPUGamaRAlloc(const int len_gamaR)
        {
            RFloat *d_gamaR;
            size_t  gamaSize = len_gamaR * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_gamaR, gamaSize));
            return d_gamaR;
        }
        //! limiter variables alloc
        RFloat *GPULIMITAlloc(const int nEqn, const int nTotal)
        {
            RFloat *d_limiter;
            size_t  limiterSize = nTotal * nEqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_limiter, limiterSize));
            return d_limiter;
        }
        RFloat *GPULimitAlloc(const int nTotal)
        {
            RFloat *d_limiter;
            size_t  limiterSize = nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_limiter, limiterSize));
            return d_limiter;
        }
        void GPUXtnYtnZtnAlloc(const int nTotalFace)
        {
            size_t tnSize = nTotalFace * sizeof(RDouble);
            HANDLE_API_ERR(cudaMalloc((void **)&d_xtn, tnSize));
            HANDLE_API_ERR(cudaMalloc((void **)&d_ytn, tnSize));
            HANDLE_API_ERR(cudaMalloc((void **)&d_ztn, tnSize));
        }
        RDouble *GPUVgnAlloc(const int nTotalFace)
        {
            RDouble *vgn;
            size_t   vgnSize = nTotalFace * sizeof(RDouble);
            HANDLE_API_ERR(cudaMalloc((void **)&vgn, vgnSize));
            return vgn;
        }
        //! d_epsCell alloc
        RFloat *GPUepsCellAlloc(const int nTotalCell)
        {
            RFloat *epsCell;
            size_t  epsCellSize = nTotalCell * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&epsCell, epsCellSize));
            return epsCell;
        }
        //! dmin dmax alloc and copy
        void GPUDminDmaxAlloc(const int nTotal)
        {
            size_t dSize = nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_dmin, dSize));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dmax, dSize));
            /*
          HANDLE_API_ERR(cudaMemcpy(d_dmin, dmin, dSize,
     cudaMemcpyHostToDevice)); HANDLE_API_ERR(cudaMemcpy(d_dmax, dmax, dSize,
     cudaMemcpyHostToDevice)); #ifdef CUDAUNITTEST1
          //! test the copy of dmax, dmin
                  TestGPUDminDmaxAllocCopy(dmax, d_dmax, nTotal);
          #endif
  */
        }

        void GPUFaceProxyMemAlloc(const int d_SEG_LEN, const int nEqn)
        {
            size_t sizeFaceRFloat = d_SEG_LEN * nEqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_flux, sizeFaceRFloat));
            size_t sizeFaceDouble = d_SEG_LEN * sizeof(RDouble);
            HANDLE_API_ERR(cudaMalloc((void **)&d_deltl, sizeFaceDouble));
            HANDLE_API_ERR(cudaMalloc((void **)&d_deltr, sizeFaceDouble));
        }

        void GPUNSFaceVarMemAlloc(const int neqn, const int ns, const int nlen)
        {
            d_rho_ds = 0;
            d_hint_s = 0;
            if (!(ns == 0 || nlen == 0))
            {
                size_t sizeGeneralField = ns * nlen * sizeof(RFloat);
                HANDLE_API_ERR(cudaMalloc((void **)&d_rho_ds, sizeGeneralField));
                HANDLE_API_ERR(cudaMalloc((void **)&d_hint_s, sizeGeneralField));
            }

            size_t sizeFaceLen = nlen * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_kcp, sizeFaceLen));
            HANDLE_API_ERR(cudaMalloc((void **)&d_mul, sizeFaceLen));
            HANDLE_API_ERR(cudaMalloc((void **)&d_mut, sizeFaceLen));

            size_t sizeFaceNeqnLen  = neqn * nlen * sizeof(RFloat);
            size_t sizeFaceNeqnLenT = 1 * nlen * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_prim, sizeFaceNeqnLen));
            HANDLE_API_ERR(cudaMalloc((void **)&d_tm, sizeFaceNeqnLenT));
        }

        void GPUViscousCoeffMemAlloc(const int nTotal)
        {
            size_t sizeViscCoeff = nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_visl, sizeViscCoeff));
            HANDLE_API_ERR(cudaMalloc((void **)&d_vist, sizeViscCoeff));
        }

        void GPUResMemAlloc(const int neqn, const int nTotal)
        {
            size_t sizeRes = nTotal * neqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_res_ns, sizeRes));
            HANDLE_API_ERR(cudaMalloc((void **)&d_dq_ns, sizeRes));
        }

        //! memory alloc for turb parts:  ql, qr
        void GPUQLRTurbMemAlloc(const int n_turb, const int SEG_LEN)
        {
            //! size_t qSize= n_turb * nTotal* sizeof(RFloat);
            size_t qlrSize = n_turb * SEG_LEN * sizeof(RFloat);
            //! HANDLE_API_ERR(cudaMalloc((void **) &q_q_turb, qSize));
            HANDLE_API_ERR(cudaMalloc((void **)&d_ql_turb, qlrSize));
            HANDLE_API_ERR(cudaMalloc((void **)&d_qr_turb, qlrSize));
        }
        RFloat *GPUFluxTurbAlloc(const int n_turb, const int SEG_LEN)
        {
            RFloat *flux_turb;
            size_t  fluxSize = n_turb * SEG_LEN * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&flux_turb, fluxSize));
            return flux_turb;
        }
        RFloat *GPUFluxTurbAllocSub(const int n_turb, const int SEG_LEN)
        {
            RFloat *flux_turb;
            size_t  fluxSize = n_turb * SEG_LEN * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&flux_turb, fluxSize));
            return flux_turb;
        }
        RFloat *GPUResTurbAlloc(const int neqn, const int nTotal)
        {
            RFloat *res_turb;
            size_t  vSize = neqn * nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&res_turb, vSize));
            return res_turb;
        }

        void GPUTurbFaceVarMemAlloc(const int neqn, const int nlen)
        {
            size_t sizeMulMut = nlen * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_mul_turb, sizeMulMut));
            HANDLE_API_ERR(cudaMalloc((void **)&d_mut_turb, sizeMulMut));
            size_t sizeMltPrim = neqn * nlen * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_mlt_turb, sizeMltPrim));
            HANDLE_API_ERR(cudaMalloc((void **)&d_prim_turb, sizeMltPrim));
        }
        void GPUParametersMemAllocCopy(const int neqn)
        {
            size_t sizeTurboo = neqn * sizeof(double);
            HANDLE_API_ERR(cudaMalloc((void **)&d_turboo, sizeTurboo));
            double *turboo = reinterpret_cast<double *>(GlobalDataBase::GetDataPtr("turboo"));
            HANDLE_API_ERR(cudaMemcpy(d_turboo, turboo, sizeTurboo, cudaMemcpyHostToDevice));
        }

        /*
                void CallGPUTurbParameterSet(const double
   reference_density_farfield, const double turb_cmu, const double SMALL){

                        GPUTurbParameterSet<<<1,
   1>>>(reference_density_farfield, turb_cmu, SMALL, d_turboo);

                }
*/

        void GPUSpecTurbMemAlloc(const int n_turb, const int nTotal)
        {
            size_t sizeSpecTurb = n_turb * nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_spec_turb, sizeSpecTurb));
        }

        void GPUMatTurbLRMemAlloc(const int n_turb, const int nTotalFace)
        {
            size_t sizeMatTurbLR = n_turb * nTotalFace * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_mat_turbl, sizeMatTurbLR));
            HANDLE_API_ERR(cudaMalloc((void **)&d_mat_turbr, sizeMatTurbLR));
        }

        void GPUDtMemAlloc(const int nTotal)
        {
            size_t sizeDt = nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_dt, sizeDt));
        }

        void GPUWallDistMemAlloc(const int nTotalCell)
        {
            size_t sizeWallDist = nTotalCell * sizeof(RDouble);
            HANDLE_API_ERR(cudaMalloc((void **)&d_walldist, sizeWallDist));
        }
        //! for farfield variables
        RFloat *GPUPrim_infMemAlloc(const int neqn)
        {
            RFloat *prim_inf;
            size_t  pSize = neqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&prim_inf, pSize));
            return prim_inf;
        }
        RFloat *GPURhs_nsMemAlloc(const int neqn, const int nTotal)
        {
            RFloat *rhs_ns;
            size_t  pSize = neqn * nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&rhs_ns, pSize));
            return rhs_ns;
        }
        RFloat *GPUInvSpectrumMalloc(const int nTotalCell)
        {
            RFloat *invSpectrum;
            size_t  iSize = nTotalCell * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&invSpectrum, iSize));
            return invSpectrum;
        }
        double *GPUDtMalloc(const int n)
        {
            RFloat *dt;
            size_t  iSize = n * sizeof(double);
            HANDLE_API_ERR(cudaMalloc((void **)&dt, iSize));
            return dt;
        }
        void GPUTurbRhsMemAlloc(const int nEqn_turb, const int nTotal)
        {
            size_t sizeRhs = nEqn_turb * nTotal * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_rhs_turb, sizeRhs));
        }

        void GPUTurbTempProxyMemAlloc(const int nEqn_turb, const int nTotal)
        {
            size_t sizeTmp = nEqn_turb * nTotal * sizeof(RFloat);
            //! HANDLE_API_ERR(cudaMalloc((void **) &d_q_proxy_turb, sizeTmp));
            //! d_q_proxy_turb = d_q_turb_proxy;
            HANDLE_API_ERR(cudaMalloc((void **)&d_q_proxy_temp_turb, sizeTmp));
        }

        void GPUTurbUpdateFlowNegPosMemAlloc()
        {
            size_t sizeNegPos = 1 * sizeof(int);
            HANDLE_API_ERR(cudaMalloc((void **)&d_n_neg, sizeNegPos));
            HANDLE_API_ERR(cudaMalloc((void **)&d_n_pos, sizeNegPos));
        }
        void GPUVistMMalloc()
        {
            size_t vSize = 1 * sizeof(double);
            HANDLE_API_ERR(cudaMalloc((void **)&d_vistmax, vSize));
            HANDLE_API_ERR(cudaMalloc((void **)&d_vistmin, vSize));
        }

        void GPUWallDistMemCopy(const int nTotalCell, const RDouble *h_wallDist)
        {
            size_t sizeWallDist = nTotalCell * sizeof(RDouble);
            HANDLE_API_ERR(cudaMemcpy(d_walldist, h_wallDist, sizeWallDist, cudaMemcpyHostToDevice));
        }

        void GPUComputeNodeValueMemAlloc(int nTotalNode, int nEqn)
        {
            size_t nodeSize    = nTotalNode * sizeof(RFloat);
            size_t nodeSizeInt = nTotalNode * sizeof(int);
            HANDLE_API_ERR(cudaMalloc((void **)&d_qNode, nodeSize * nEqn));
            HANDLE_API_ERR(cudaMalloc((void **)&d_tNode, nodeSize));
            HANDLE_API_ERR(cudaMalloc((void **)&d_nCount, nodeSizeInt));
            HANDLE_API_ERR(cudaMalloc((void **)&d_nodeWeight, nodeSize));
            //! For three components of velocity
            HANDLE_API_ERR(cudaMalloc((void **)&d_qVelNode, nodeSize * 3));
        }

        void GPUComputeNodeQTurbValueMemAlloc(const int nTotalNode, const int nEqn)
        {
            size_t nodeSize = nTotalNode * nEqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_qTurbNode, nodeSize));
        }

        void GPUInterPointMemAllocCopy(const int nIPoint, const int nEqn, RFloat **h_qInterPoint,
                                       RFloat **h_tInterPoint)
        {
            //! size_t pointSizeQ = nIPoint * nEqn * sizeof(int);
            //! size_t pointSizeT = nIPoint * 1 * sizeof(int); //!just one equation for t
            size_t pointSizeQ = nIPoint * nEqn * sizeof(RFloat);
            size_t pointSizeT = nIPoint * 1 * sizeof(RFloat); //! just one equation for t
            HANDLE_API_ERR(cudaMalloc((void **)&d_qInterPoint, pointSizeQ));
            HANDLE_API_ERR(cudaMalloc((void **)&d_tInterPoint, pointSizeT));
            HANDLE_API_ERR(cudaMemcpy(d_qInterPoint, h_qInterPoint, pointSizeQ, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMemcpy(d_tInterPoint, h_tInterPoint, pointSizeT, cudaMemcpyHostToDevice));
        }
        void GPUTurbInterPointMemAlloc(const int nIPoint, const int nEqn)
        {
            //! size_t pointSizeQ = nIPoint * nEqn * sizeof(int);
            size_t pointSizeQ = nIPoint * nEqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_qTurbInterPoint, pointSizeQ));
        }

        void GPUInterfaceTransferStorageNS(const int nIFace, const int nEqn)
        {
            size_t sizeQ = nIFace * nEqn * sizeof(RFloat);
            size_t sizeT = nIFace * 1 * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_q, sizeQ));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_t, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_q, sizeQ));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_t, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_dqdx, sizeQ));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_dqdx, sizeQ));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_dqdy, sizeQ));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_dqdy, sizeQ));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_dqdz, sizeQ));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_dqdz, sizeQ));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_limit, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_limit, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_dtdx, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_dtdx, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_dtdy, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_dtdy, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_dtdz, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_dtdz, sizeT));
        }

        void GPUInterpointTransferStorageNS(const int nIPoint, const int nEqn)
        {
            size_t sizeQ = nIPoint * nEqn * sizeof(RFloat);
            size_t sizeT = nIPoint * 1 * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_qNode, sizeQ));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_tNode, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_qNode, sizeQ));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_tNode, sizeT));
        }
        void GPUInterfaceTransferStorageTurb(const int nIFace, const int nEqn)
        {
            size_t sizeT = nIFace * nEqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_qTurb, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_qTurb, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_dqTurbdx, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_dqTurbdx, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_dqTurbdy, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_dqTurbdy, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_dqTurbdz, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_dqTurbdz, sizeT));
        }

        void GPUInterpointTransferStorageTurb(const int nIPoint, const int nEqn)
        {
            size_t sizeT = nIPoint * nEqn * sizeof(RFloat); //! In turb solver, nEqn_turb = 1
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_send_qTurbNode, sizeT));
            HANDLE_API_ERR(cudaMalloc((void **)&d_fg_recv_qTurbNode, sizeT));
        }

        void GPUTransferInterpointBufferNS(InterpointFields *interpointFields, const int solverNS, const int nIPoint)
        {
            if (!interpointFields)
            {
                return;
            }
            int iData;
            int solverID;
            int dataDimension;
            int sum_dim            = 0;
            nameNSInterpointBuffer = new vector<string>;
            for (iData = 0; iData < interpointFields->Size(); iData++)
            {
                solverID = interpointFields->GetSolverIndex(iData);
                if (solverID != solverNS) continue;
                string &dataName = interpointFields->GetName(iData);
                nameNSInterpointBuffer->push_back(dataName);
                dataDimension = interpointFields->GetDimension(iData);
                sum_dim += dataDimension;
            }
            lengthInterpointBufferNS = sum_dim * nIPoint;

            int sumBufferLength      = 0;
            int countBuffer          = 0;
            offsetInterpointBufferNS = new int[nameNSInterpointBuffer->size()];
            dimInterpointBufferNS    = new int[nameNSInterpointBuffer->size()];
            for (iData = 0; iData < interpointFields->Size(); iData++)
            {
                solverID = interpointFields->GetSolverIndex(iData);
                if (solverID != solverNS) continue;
                dataDimension                         = interpointFields->GetDimension(iData);
                offsetInterpointBufferNS[countBuffer] = sumBufferLength;
                dimInterpointBufferNS[countBuffer]    = dataDimension;
                sumBufferLength += dataDimension * nIPoint;
                countBuffer++;
            }
            size_t sizeBufferNS = lengthInterpointBufferNS * sizeof(double);
            HANDLE_API_ERR(cudaHostAlloc((void **)&HostSendInterpointBufferNS, sizeBufferNS, cudaHostAllocDefault));
            HANDLE_API_ERR(cudaHostAlloc((void **)&HostRecvInterpointBufferNS, sizeBufferNS, cudaHostAllocDefault));
            HANDLE_API_ERR(cudaMalloc((void **)&GPUSendInterpointBufferNS, sizeBufferNS));
            HANDLE_API_ERR(cudaMalloc((void **)&GPURecvInterpointBufferNS, sizeBufferNS));
        }

        void GPUTransferInterpointBufferTurb(InterpointFields *interpointFields, const int solverTurb,
                                             const int nIPoint)
        {
            if (!interpointFields)
            {
                return;
            }
            int iData;
            int solverID;
            int dataDimension;
            int sum_dim              = 0;
            nameTurbInterpointBuffer = new vector<string>;
            for (iData = 0; iData < interpointFields->Size(); iData++)
            {
                solverID = interpointFields->GetSolverIndex(iData);
                if (solverID != solverTurb) continue;
                string &dataName = interpointFields->GetName(iData);
                nameTurbInterpointBuffer->push_back(dataName);
                dataDimension = interpointFields->GetDimension(iData);
                sum_dim += dataDimension;
            }
            lengthInterpointBufferTurb = sum_dim * nIPoint;

            int sumBufferLength        = 0;
            int countBuffer            = 0;
            offsetInterpointBufferTurb = new int[nameTurbInterpointBuffer->size()];
            dimInterpointBufferTurb    = new int[nameTurbInterpointBuffer->size()];
            for (iData = 0; iData < interpointFields->Size(); iData++)
            {
                solverID = interpointFields->GetSolverIndex(iData);
                if (solverID != solverTurb) continue;
                dataDimension                           = interpointFields->GetDimension(iData);
                offsetInterpointBufferTurb[countBuffer] = sumBufferLength;
                dimInterpointBufferTurb[countBuffer]    = dataDimension;
                sumBufferLength += dataDimension * nIPoint;
                countBuffer++;
            }
            size_t sizeBufferTurb = lengthInterpointBufferTurb * sizeof(double);
            HANDLE_API_ERR(cudaHostAlloc((void **)&HostSendInterpointBufferTurb, sizeBufferTurb, cudaHostAllocDefault));
            HANDLE_API_ERR(cudaHostAlloc((void **)&HostRecvInterpointBufferTurb, sizeBufferTurb, cudaHostAllocDefault));
            HANDLE_API_ERR(cudaMalloc((void **)&GPUSendInterpointBufferTurb, sizeBufferTurb));
            HANDLE_API_ERR(cudaMalloc((void **)&GPURecvInterpointBufferTurb, sizeBufferTurb));
        }

        void GPUInterfaceTransferBufferNS(InterfaceFields *interfaceFields, const int solverNS, const int nIFace)
        {
            if (!interfaceFields)
            {
                return;
            }
            int iData;
            int solverID;
            int dataDimension;
            int sum_dim      = 0;
            int sum_dimFirst = 0;
            //! nameBufferNS = new vector<string>;
            vector<string> nameBufferNS;
            //! set name and total size of buffer NS
            for (iData = 0; iData < interfaceFields->Size(); iData++)
            {
                solverID = interfaceFields->GetSolverIndex(iData);
                if (solverID != solverNS) continue;
                string &dataName = interfaceFields->GetName(iData);
                nameBufferNS.push_back(dataName);
                dataDimension = interfaceFields->GetDim(iData);
                //! if (dataName == "limit") dataDimension = 1;
                //! add for test
                //! printf("iData = %d, dataName = %s\n", iData, dataName.c_str());
                //! add end
                sum_dim += dataDimension;
                if (dataName == "q") sum_dimFirst += dataDimension;
            }
            nameNSBuffer = new vector<string>;
            for (int nameID = 0; nameID < nameBufferNS.size(); nameID++)
            {
                nameNSBuffer->push_back(nameBufferNS[nameID]);
            }
            lengthBufferNS = sum_dim * nIFace;
            //! set offset and dim
            int sumBufferLength = 0;
            int countBuffer     = 0;
            offsetBufferNS      = new int[nameBufferNS.size()];
            dimBufferNS         = new int[nameBufferNS.size()];
            for (iData = 0; iData < interfaceFields->Size(); iData++)
            {
                solverID = interfaceFields->GetSolverIndex(iData);
                if (solverID != solverNS) continue;
                string &dataName = interfaceFields->GetName(iData);
                dataDimension    = interfaceFields->GetDim(iData);
                //! determined in NSSolverUnstruct::CreateLimiter, in this case, nVariables
                //! = 1. Thus, dataDimension should be 1. It can also be referenced by
                //! PHSPACE::CommunicateInterfaceValue(grid, limit2D, limiterVarName, 1) in
                //! Limiter::Calculation
                if (dataName == "limit") dataDimension = 1;
                offsetBufferNS[countBuffer] = sumBufferLength;
                dimBufferNS[countBuffer]    = dataDimension;
                sumBufferLength += dataDimension * nIFace;
                countBuffer++;
            }

            size_t sizeBufferNS = lengthBufferNS * sizeof(double);
            HANDLE_API_ERR(cudaHostAlloc((void **)&HostSendBufferNS, sizeBufferNS, cudaHostAllocDefault));
            HANDLE_API_ERR(cudaHostAlloc((void **)&HostRecvBufferNS, sizeBufferNS, cudaHostAllocDefault));
            HANDLE_API_ERR(cudaMalloc((void **)&GPUSendBufferNS, sizeBufferNS));
            HANDLE_API_ERR(cudaMalloc((void **)&GPURecvBufferNS, sizeBufferNS));
        }

        void GPUInterfaceTransferBufferTurb(InterfaceFields *interfaceFields, const int solverTurb, const int nIFace)
        {
            if (!interfaceFields)
            {
                return;
            }
            int iData;
            int solverID;
            int dataDimension;
            int sum_dim = 0;
            //! nameBufferNS = new vector<string>;
            vector<string> nameBufferTurb;
            for (iData = 0; iData < interfaceFields->Size(); iData++)
            {
                solverID = interfaceFields->GetSolverIndex(iData);

                if (solverID != solverTurb) continue;
                string &dataName = interfaceFields->GetName(iData);
                nameBufferTurb.push_back(dataName);
                dataDimension = interfaceFields->GetDim(iData);
                //! add for test
                //! printf("iData = %d, dataName = %s\n", iData, dataName.c_str());
                //! add end
                sum_dim += dataDimension;
            }
            //! add for test
            nameTurbBuffer = new vector<string>;
            for (int nameID = 0; nameID < nameBufferTurb.size(); nameID++)
            {
                nameTurbBuffer->push_back(nameBufferTurb[nameID]);
            }

            lengthBufferTurb = sum_dim * nIFace;
            //! set offset
            int sumBufferLength = 0;
            int countBuffer     = 0;
            offsetBufferTurb    = new int[nameBufferTurb.size()];
            dimBufferTurb       = new int[nameBufferTurb.size()];
            for (iData = 0; iData < interfaceFields->Size(); iData++)
            {
                solverID = interfaceFields->GetSolverIndex(iData);
                if (solverID != solverTurb) continue;
                string &dataName              = interfaceFields->GetName(iData);
                dataDimension                 = interfaceFields->GetDim(iData);
                offsetBufferTurb[countBuffer] = sumBufferLength;
                dimBufferTurb[countBuffer]    = dataDimension;
                sumBufferLength += dataDimension * nIFace;
                countBuffer++;
            }
            /*
  //!add for test
  for (int nameID = 0; nameID < nameBufferTurb.size(); nameID++){
          printf("nameID = %d, nameNSBuffer = %s, offset = %d\n", nameID,
  (*nameTurbBuffer)[nameID].c_str(), offsetBufferTurb[nameID]);
  }
  //!add end
  */
            size_t sizeBufferTurb = lengthBufferTurb * sizeof(double);
            HANDLE_API_ERR(cudaHostAlloc((void **)&HostSendBufferTurb, sizeBufferTurb, cudaHostAllocDefault));
            HANDLE_API_ERR(cudaHostAlloc((void **)&HostRecvBufferTurb, sizeBufferTurb, cudaHostAllocDefault));
            HANDLE_API_ERR(cudaMalloc((void **)&GPUSendBufferTurb, sizeBufferTurb));
            HANDLE_API_ERR(cudaMalloc((void **)&GPURecvBufferTurb, sizeBufferTurb));
        }

        int GetIndexBufferNS(const string &name)
        {
            int nameID;
            for (nameID = 0; nameID < nameNSBuffer->size(); nameID++)
            {
                if ((*nameNSBuffer)[nameID] == name) return nameID;
            }
            return 666;
        }
        int GetIndexBufferTurb(const string &name)
        {
            int nameID;
            for (nameID = 0; nameID < nameTurbBuffer->size(); nameID++)
            {
                if ((*nameTurbBuffer)[nameID] == name) return nameID;
            }
            return 666;
        }

        int GetIndexInterpointBufferNS(const string &name)
        {
            int nameID;
            for (nameID = 0; nameID < nameNSInterpointBuffer->size(); nameID++)
            {
                if ((*nameNSInterpointBuffer)[nameID] == name) return nameID;
            }
            return 666;
        }
        int GetIndexInterpointBufferTurb(const string &name)
        {
            int nameID;
            for (nameID = 0; nameID < nameTurbInterpointBuffer->size(); nameID++)
            {
                if ((*nameTurbInterpointBuffer)[nameID] == name) return nameID;
            }
            return 666;
        }
    } //! namespace GPUFlowVariables
    /*
           __global__ void GPUTurbParameterSet(const double
      reference_density_farfield, const double turb_cmu, const double SMALL, double
      * turboo){
   
                   turboo[0] =  3.0 / ( reference_density_farfield * turb_cmu +
      SMALL );
           }
   */
    namespace GPUGeomVariables
    {
        int d_nTotalFace = 0;
        int d_nTotalCell = 0;
        //! The number of boundary faces which include interaces
        int d_nBoundFace = 0;
        //! The number of interface between different zones
        int d_nIFace = 0;
        //! sum of TotalCell and BoundFace
        int d_nTotal = 0;
        //! points information
        RDouble *d_x;
        RDouble *d_y;
        RDouble *d_z;
        //! face information
        RDouble *d_xfn;
        RDouble *d_yfn;
        RDouble *d_zfn;
        RDouble *d_area;
        RDouble *d_xfc;
        RDouble *d_yfc;
        RDouble *d_zfc;
        RDouble *d_nxs;
        RDouble *d_nys;
        RDouble *d_nzs;
        //! cell center information
        RDouble *d_xcc;
        RDouble *d_ycc;
        RDouble *d_zcc;
        //! Cell skewness data, excluding ghosts.
        RDouble *d_cellSkewness;
        //! cell volume
        RDouble *d_vol;
        //! face2cell face2node relationship
        int *d_left_cell_of_face;
        int *d_right_cell_of_face;
        int *d_face2node;
        int *d_node_number_of_each_face;

        long long int *d_nodePosiFace;
        int           *d_boundaryType;
        //! face Cell realationship
        //! int * d_face_number_of_each_cell;
        int *d_cell2Face;
        int *d_posiCell2Face;

        int *d_cell2face;
        int *d_face_number_of_each_cell;
        int *d_acc_face_number_of_cell;

        //! cell2node for cell and node relationship, only used in loop mode optimization
        int *d_cell2Node;
        int *d_posiCell2Node;
        int *d_cell2NodeCount;
        int *d_node_number_of_each_cell;
        //! cell2Cell for cell and its neighbor cell relationship
        int *d_cell2Cell;
        int *d_posiCell2Cell;
        int *d_neighbour_cell_number_of_each_cell;

        //! interpolation points infomation
        int *d_interPoint2GlobalPoint;
        int *d_cellNumberOfInterPoint;
        int *d_labelOfInterPoint;
        //! interface information
        int *d_interFace2BoundaryFace;
        //! interpoint information

        //! cell coloring information
        int *d_InteriorCellGroup;
        //! overset cell information
        int *d_blank;

        //! void GPUGeomInfoAllocCopy(Grid *grid)
        //!void TestGetRegion(Region *) {}

        void GPUGeomInfoAllocCopy()
        {
            int nLocalZones = GetNumberofLocalZones();
            for (int izone = 0; izone < nLocalZones; izone++)
            {
                int zoneID = GetLocalZoneIDToGlobalZoneID(izone);
                //! Zone *zone = GetZone(zoneID);
                int level = 0;
                //! Grid *grid = GetGrid(zoneID, level);
                UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, level));
                //! set d_nTotalFace, d_nTotalCell, d_nBoundFace, d_nIFace
                GPUTotalNumAllocCopy(grid);
                //! set d_x, d_y, d_z and test it
                const int nTotalNode = grid->GetNTotalNode();
                RDouble  *h_x        = grid->GetX();
                RDouble  *h_y        = grid->GetY();
                RDouble  *h_z        = grid->GetZ();
                GPUNodesDataAllocCopy(h_x, h_y, h_z, nTotalNode);
#ifdef CUDAUNITTEST
                using namespace GPUTestSpace;
#endif
//! GPUTestSpace::TestGPUNodesDataAllocCopy(h_x, h_y, h_z, d_x, d_y, d_z,
//! nTotalNode);
#ifdef CUDAUNITTEST
                TestGPUNodesDataAllocCopy(h_x, h_y, h_z, d_x, d_y, d_z, nTotalNode);
#endif
                //! set d_xcc, d_ycc, d_zcc and test it
                RDouble  *h_xcc      = grid->GetCellCenterX();
                RDouble  *h_ycc      = grid->GetCellCenterY();
                RDouble  *h_zcc      = grid->GetCellCenterZ();
                const int nTotalCell = grid->GetNTotalCell();
                const int nBoundFace = grid->GetNBoundFace(); //! change the size of xcc ycc zcc
                int       nTotal     = nTotalCell + nBoundFace;
                GPUCellCentAllocCopy(h_xcc, h_ycc, h_zcc, nTotal);
#ifdef CUDAUNITTEST
                TestGPUCellCentAllocCopy(h_xcc, h_ycc, h_zcc, d_xcc, d_ycc, d_zcc, nTotal);
#endif
                //! set face normal d_xfn ... and test it
                RDouble  *h_xfn      = grid->GetFaceNormalX();
                RDouble  *h_yfn      = grid->GetFaceNormalY();
                RDouble  *h_zfn      = grid->GetFaceNormalZ();
                const int nTotalFace = grid->GetNTotalFace();

                GPUFaceNormAllocCopy(h_xfn, h_yfn, h_zfn, nTotalFace);
//! test the AllocCopy
#ifdef CUDAUNITTEST
                TestGPUFaceNormAllocCopy(h_xfn, h_yfn, h_zfn, d_xfn, d_yfn, d_zfn, nTotalFace);
#endif

                //! set face center d_xfc ... and test it
                RDouble *h_xfc = grid->GetFaceCenterX();
                RDouble *h_yfc = grid->GetFaceCenterY();
                RDouble *h_zfc = grid->GetFaceCenterZ();

                GPUFaceCentAllocCopy(h_xfc, h_yfc, h_zfc, nTotalFace);
#ifdef CUDAUNITTEST
                TestGPUFaceCentAllocCopy(h_xfc, h_yfc, h_zfc, d_xfc, d_yfc, d_zfc, nTotalFace);
#endif
                //! set faace area and cell vol ... and test it
                RDouble *h_area = grid->GetFaceArea();
                RDouble *h_vol  = grid->GetCellVolume();
                //! Noting that for vol, ghost cells should be considered. Thus, nTotal
                //! should be used to replace nTotalCell.
                GPUAreaVolmAllocCopy(h_area, h_vol, nTotalFace, nTotal);
                RDouble *h_cellSkewness = grid->GetCellSkewness();
                GPUCellSkewnessAllocCopy(h_cellSkewness, nTotalCell);
#ifdef CUDAUNITTEST
                TestGPUAreaVolmAllocCopy(h_area, h_vol, d_area, d_vol, nTotalFace, nTotalCell);
#endif

                //! set left_cell_of_face and right_cell_of_face and test
                const int *h_left_cell_of_face  = grid->GetLeftCellOfFace();
                const int *h_right_cell_of_face = grid->GetRightCellOfFace();
                GPUFaceCellRelAllocCopy(h_left_cell_of_face, h_right_cell_of_face, nTotalFace);
#ifdef CUDAUNITTEST
                TestGPUFaceCellRelAllocCopy(h_left_cell_of_face, h_right_cell_of_face, d_left_cell_of_face,
                                            d_right_cell_of_face, nTotalFace);
#endif
                //! set face2node and node_numb er_of_each_face
                const int *h_face2node                = grid->GetFace2Node();
                const int *h_node_number_of_each_face = grid->GetNodeNumberOfEachFace();
                int        nTotalFace2Node            = 0;
                for (int i = 0; i < nTotalFace; i++)
                {
                    for (int j = 0; j < h_node_number_of_each_face[i]; j++)
                    {
                        nTotalFace2Node++;
                    }
                }
                GPUFaceNodeRelAllocCopy(h_face2node, h_node_number_of_each_face, nTotalFace2Node, nTotalFace);

#ifdef CUDAUNITTEST
                TestGPUFaceNodeRelAllocCopy(h_face2node, h_node_number_of_each_face, d_face2node,
                                            d_node_number_of_each_face, nTotalFace2Node, nTotalFace);
#endif

                int **cell2face                = grid->GetCell2Face();
                int  *face_number_of_each_cell = grid->GetFaceNumberOfEachCell();
#ifdef REORDER
                /*
        int numCellFaces = face_number_of_each_cell[0];
        for (int iFace = 0; iFace < numCellFaces; iFace++){
                int faceID = cell2face[0][iFace];
                int leftCellID = h_left_cell_of_face[faceID];
                int rightCellID = h_right_cell_of_face[faceID];
                printf("cellID = %d, faceID=%d, leftCellID = %d, rightCellID =
        %d\n", 0, faceID, leftCellID, rightCellID);
        }
        */
#endif

                GPUFaceCellRelMemcpy(face_number_of_each_cell, cell2face, nTotalCell);

                //! set nodePosiFace variable,
                //!TODO:long long int for h_nodePosiFace will be error in CUDAUNITTEST
                const long long int *h_nodePosiFace     = grid->GetFace2NodeSubscript();
                const int            nTotalNodePosiFace = grid->GetNTotalFace() + 1;
                GPUNodePosiFaceAllocCopy(h_nodePosiFace, nTotalNodePosiFace);
#ifdef CUDAUNITTEST
                TestGPUNodePosiFaceAllocCopy(h_nodePosiFace, d_nodePosiFace, nTotalNodePosiFace);
#endif
                //! set boundary type variable
                UnstructBCSet **bcr = grid->GetBCRecord();
                //! const int nBoundFace = grid->GetNBoundFace();
                int *h_boundaryType = new int[nBoundFace];
                for (int iface = 0; iface < nBoundFace; iface++)
                {
                    h_boundaryType[iface] = bcr[iface]->GetKey();
                }
                GPUBoundaryTypeAllocCopy(h_boundaryType, nBoundFace);
#ifdef CUDAUNITTEST
                TestGPUBoundaryTypeAllocCopy(h_boundaryType, d_boundaryType, nBoundFace);
#endif
                delete[] h_boundaryType;
                //! set nxs, nys, nzs
                RDouble *h_nxs = grid->GetFaceNormalX();
                RDouble *h_nys = grid->GetFaceNormalY();
                RDouble *h_nzs = grid->GetFaceNormalZ();
                GPUFaceNormalAllocCopy(h_nxs, h_nys, h_nzs, nTotalFace);
                //! establish cell face relationship on device
                //! int **cell2face  = grid->GetCell2Face();
                //! int * face_number_of_each_cell = grid->GetFaceNumberOfEachCell();
                GPUCellFaceAllocCopy(cell2face, face_number_of_each_cell, nTotalCell);
                //! transfer interpolation points information from host to device
                InterpointInformation *interPointInfor = grid->GetInterpointInfo();
                //! interPointInfor only exists in multi-process case
                if (interPointInfor)
                {
                    int  nIPoint                = interPointInfor->GetNumberOfInterpoints();
                    int *interPoint2GlobalPoint = interPointInfor->GetInterPoint2GlobalPoint();
                    int *cellNumberOfInterPoint = interPointInfor->GetCellNumberOfInterPoint();
                    int *labelOfInterPoint      = interPointInfor->GetLabelOfInterPoint();
                    GPUInterPointsAllocCopy(nIPoint, interPoint2GlobalPoint, cellNumberOfInterPoint, labelOfInterPoint);
#ifdef CUDAUNITTEST
                    TestGPUInterPointsAllocCopy(nIPoint, interPoint2GlobalPoint, cellNumberOfInterPoint,
                                                labelOfInterPoint);
#endif
                    //! TestGPUInterPointsAllocCopy();
                }
                InterfaceInfo *iinfo = grid->GetInterfaceInfo();
                //! whether or not interface exists
                if (iinfo)
                {
                    int  nIFace                 = iinfo->GetNIFace();
                    int *interFace2BoundaryFace = iinfo->GetInterFace2BoundaryFace();
                    GPUInterFace2BoundaryFace(nIFace, interFace2BoundaryFace);
#ifdef CUDAUNITTEST
                    TestGPUInterFace2BoundaryFace(nIFace, interFace2BoundaryFace);
#endif
                }

#ifdef MCLUSGS
                const int *h_blank = grid->GetBlankIndex();
                GPUBlankAlloc(nTotalCell, h_blank);
#endif

                //!#ifdef LOOPMODEOPT
                //! allocate d_cell2Node and corredponding variables
                const int *h_cell2node                = grid->GetCell2Node();
                const int *h_node_number_of_each_cell = grid->GetNodeNumberOfEachCell();
                //! noting that node_number_of_each_cell does not contain ghost cells!!!
                //! printf("h_node_number_of_each_cell[%d] = %d\n", nTotalCell,
                //! h_node_number_of_each_cell[nTotalCell]); Set d_cell2Node,
                GPUCell2NodeNumberOfEachCellCopy(nTotalCell, nTotalFace, nBoundFace, h_cell2node,
                                                 h_node_number_of_each_cell, h_face2node, h_node_number_of_each_face,
                                                 h_left_cell_of_face, h_right_cell_of_face);
                //! GPUCell2NodeSet(nTotalFace, );
                //! set cell2NodeCount

                //! allocate d_cell2Cell and corresponding variables
                vector<int> *h_cell2Cell                  = grid->GetCell2Cell();
                const int   *ngb_cell_number_of_each_cell = grid->GetNumberOfNeighborCell();
#ifdef MATRIXOUTPUT
                //! test cell index reorder results
                using namespace Reorder;
                //! CallBandwidthFTS(nTotalCell, h_cell2Cell, ngb_cell_number_of_each_cell);
                //! CallCell2CellMatrixOutput(nTotalCell, h_cell2Cell,
                //! ngb_cell_number_of_each_cell);
                OutputMatrixProperties();
#endif
                GPUCell2CellNumberOfEachCellCopy(nTotalCell, h_cell2Cell, ngb_cell_number_of_each_cell);
                /*
    //!test face_number_of_each_cell and ngb_cell_number_of_each_cell
    for (int iCell = 0; iCell < nTotalCell; iCell++){
            int diff = face_number_of_each_cell[iCell] -
    ngb_cell_number_of_each_cell[iCell]; if (diff != 0) { printf("iCell = %d,
    face_number_of_each_cell = %d, ngb_cell_number_of_each_cell = %d\n", iCell,
    face_number_of_each_cell[iCell], ngb_cell_number_of_each_cell[iCell]); for
    (int iNgbCell = 0; iNgbCell < ngb_cell_number_of_each_cell[iCell];
    iNgbCell++){ printf("h_cell2Cell[%d][%d] = %d\n", iCell, iNgbCell,
    h_cell2Cell[iCell][iNgbCell]);
                    }
                    for (int iFace = 0; iFace < face_number_of_each_cell[iCell];
    iFace++){ int faceID = cell2face[iCell][iFace]; int ownCell =
    h_left_cell_of_face[iCell]; int ngbCell = h_right_cell_of_face[iCell];
                            printf("iFace = %d, ownCell = %d, ngbCell = %d\n",
    iFace, ownCell, ngbCell);
                    }
                    exit(1);
            }
    }
    */

                //!#endif
            }
        }

        void GPUBlankAlloc(const int nTotalCell, const int *h_blank)
        {
            size_t sizeBlank = nTotalCell * sizeof(int);
            HANDLE_API_ERR(cudaMalloc((void **)&d_blank, sizeBlank));
            HANDLE_API_ERR(cudaMemcpy(d_blank, h_blank, sizeBlank, cudaMemcpyHostToDevice));
        }

        void GPUCell2CellNumberOfEachCellCopy(const int nTotalCell, vector<int> *cell2Cell,
                                              const int *ngb_cell_number_of_each_cell)
        {
            //! printf("program is running in GPUCell2CellNumberOfEachCellCopy\n");
            int numCells;
            numCells             = 0;
            int *h_posiCell2Cell = new int[nTotalCell];
            for (int iCell = 0; iCell < nTotalCell; iCell++)
            {
                h_posiCell2Cell[iCell] = numCells;
                numCells += ngb_cell_number_of_each_cell[iCell];
            }
            size_t sizeCell2Cell = numCells * sizeof(int);
            int   *h_cell2Cell1D = new int[numCells];

            int sumID = 0;
            for (int iCell = 0; iCell < nTotalCell; iCell++)
            {
                for (int iNgbCell = 0; iNgbCell < ngb_cell_number_of_each_cell[iCell]; iNgbCell++)
                {
                    h_cell2Cell1D[sumID] = cell2Cell[iCell][iNgbCell];
                    sumID++;
                }
            }

            HANDLE_API_ERR(cudaMalloc((void **)&d_cell2Cell, sizeCell2Cell));
            HANDLE_API_ERR(cudaMemcpy(d_cell2Cell, h_cell2Cell1D, sizeCell2Cell, cudaMemcpyHostToDevice));
            size_t sizeCell = nTotalCell * sizeof(int);
            HANDLE_API_ERR(cudaMalloc((void **)&d_neighbour_cell_number_of_each_cell, sizeCell));
            HANDLE_API_ERR(cudaMemcpy(d_neighbour_cell_number_of_each_cell, ngb_cell_number_of_each_cell, sizeCell,
                                      cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMalloc((void **)&d_posiCell2Cell, sizeCell));
            HANDLE_API_ERR(cudaMemcpy(d_posiCell2Cell, h_posiCell2Cell, sizeCell, cudaMemcpyHostToDevice));

            delete[] h_cell2Cell1D;
            delete[] h_posiCell2Cell;
        }

        void GPUCell2NodeNumberOfEachCellCopy(const int nTotalCell, const int nTotalFace, const int nBoundFace,
                                              const int *cell2node, const int *node_number_of_each_cell,
                                              const int *face2node, const int *node_number_of_each_face,
                                              const int *leftCellofFace, const int *rightCellofFace)
        {
            //! printf("program is running in GPUCell2NodeNumberOfEachCellCopy\n");
            int numNodes = 0;
            for (int iCell = 0; iCell < nTotalCell; iCell++)
            {
                numNodes += node_number_of_each_cell[iCell];
            }
            size_t sizeCell2Node = numNodes * sizeof(int);
            HANDLE_API_ERR(cudaMalloc((void **)&d_cell2Node, sizeCell2Node));
            HANDLE_API_ERR(cudaMalloc((void **)&d_cell2NodeCount, sizeCell2Node));
            int *cell2NodeCount = new int[numNodes];

            for (int element = 0; element < numNodes; element++)
            {
                cell2NodeCount[element] = 0;
            }

            HANDLE_API_ERR(cudaMemcpy(d_cell2Node, cell2node, sizeCell2Node, cudaMemcpyHostToDevice));
            size_t sizeNodeNumber = nTotalCell * sizeof(int);
            HANDLE_API_ERR(cudaMalloc((void **)&d_node_number_of_each_cell, sizeNodeNumber));
            HANDLE_API_ERR(cudaMemcpy(d_node_number_of_each_cell, node_number_of_each_cell, sizeNodeNumber,
                                      cudaMemcpyHostToDevice));
            int *h_posiCell2Node = new int[nTotalCell];
            HANDLE_API_ERR(cudaMalloc((void **)&d_posiCell2Node, sizeNodeNumber));
            numNodes = 0;
            for (int iCell = 0; iCell < nTotalCell; iCell++)
            {
                h_posiCell2Node[iCell] = numNodes;
                numNodes += node_number_of_each_cell[iCell];
            }
            HANDLE_API_ERR(cudaMemcpy(d_posiCell2Node, h_posiCell2Node, sizeNodeNumber, cudaMemcpyHostToDevice));
            /*
  //!check value, the last or the first one.
  //!int iCell = nTotalCell - 1;
  int iCell = 0;
  numNodes = node_number_of_each_cell[iCell];
  int startNode = h_posiCell2Node[iCell];
  for (int iNode = 0; iNode < numNodes; iNode++){
          int node = cell2node[startNode + iNode];
          printf("iCell = %d, numNodes = %d, iNode = %d, node = %d\n", iCell,
  numNodes, iNode, node);
  }
  */
            int *h_posiFace2Node = new int[nTotalFace];
            int  numFaces        = 0;
            for (int iFace = 0; iFace < nTotalFace; iFace++)
            {
                h_posiFace2Node[iFace] = numFaces;
                numFaces += node_number_of_each_face[iFace];
            }

            for (int iFace = nBoundFace; iFace < nTotalFace; iFace++)
            {
                int le = leftCellofFace[iFace];
                int re = rightCellofFace[iFace];
                int jNode;
                int nodePosition = h_posiFace2Node[iFace];
                int offset, nodeID;
                int point, start, numNodesInCell;
                for (jNode = 0; jNode < node_number_of_each_face[iFace]; ++jNode)
                {
                    point = face2node[nodePosition + jNode];
                    //! for left cell
                    start          = h_posiCell2Node[le];
                    numNodesInCell = node_number_of_each_cell[le];
                    for (offset = 0; offset < numNodesInCell; offset++)
                    {
                        nodeID = cell2node[start + offset];
                        if (nodeID == point) cell2NodeCount[start + offset]++;
                    }
                    //! for right cell
                    start          = h_posiCell2Node[re];
                    numNodesInCell = node_number_of_each_cell[re];
                    for (offset = 0; offset < numNodesInCell; offset++)
                    {
                        nodeID = cell2node[start + offset];
                        if (nodeID == point) cell2NodeCount[start + offset]++;
                    }
                }
            }
            HANDLE_API_ERR(cudaMemcpy(d_cell2NodeCount, cell2NodeCount, sizeCell2Node, cudaMemcpyHostToDevice));
            /*
                          //!check cell2NodeCount
                          int cellRecord[6];
                          for (int element = 0; element < 6; element++){
                                  cellRecord[element] = 0;
                          }
                          for (int iCell = 0; iCell < nTotalCell; iCell++ ){
                                  int numNodes =
     node_number_of_each_cell[iCell]; int start = h_posiCell2Node[iCell]; int
     jNode; int nodeID; int cellCount;

                                  for (jNode = 0; jNode < numNodes; jNode++){
                                          nodeID = cell2node[start + jNode];
                                          cellCount = cell2NodeCount[start +
     jNode]; cellRecord[cellCount]++; if (cellCount == 0) printf("nodeID = %d,
     cellCount, iCell  = %d\n", nodeID, cellCount, iCell);
                                          //!if (cellCount == 6) printf("nodeID =
     %d, cellCount = %d\n", nodeID, cellCount);
                                  }
                          }

                          for (int element = 0; element < 6; element++){
                                  printf("The number of node owned by %d cell is
     %d\n", element, cellRecord[element]);
                          }
  */

            delete[] h_posiCell2Node;
            delete[] h_posiFace2Node;
            delete[] cell2NodeCount;
        }

        //! copy nTotalXXX variables to gpu device
        void GPUTotalNumAllocCopy(Grid *grid)
        {
            d_nTotalFace = grid->GetNTotalFace();
            d_nTotalCell = grid->GetNTotalCell();
            d_nBoundFace = grid->GetNBoundFace();
            d_nIFace     = grid->GetNIFace();
            d_nTotal     = d_nTotalCell + d_nBoundFace;
        }
        //! copy raw data: points coordinates to gpu device
        void GPUNodesDataAllocCopy(RDouble *x, RDouble *y, RDouble *z, const int n)
        {
            size_t sizeNodes = n * sizeof(RDouble);
            HANDLE_API_ERR(cudaMalloc((void **)&d_x, sizeNodes));
            //! handleAPIErr(cudaMalloc((RDouble**) &d_x, sizeNodes), __FILE__, __LINE__);

            HANDLE_API_ERR(cudaMemcpy(d_x, x, sizeNodes, cudaMemcpyHostToDevice));
            //! handleAPIErr(cudaMemcpy(d_x, x, sizeNodes, cudaMemcpyHostToDevice),
            //! __FILE__, __LINE__);

            HANDLE_API_ERR(cudaMalloc((void **)&d_y, sizeNodes));
            //! handleAPIErr(cudaMalloc((RDouble**) &d_y, sizeNodes), __FILE__, __LINE__);

            HANDLE_API_ERR(cudaMemcpy(d_y, y, sizeNodes, cudaMemcpyHostToDevice));
            //! handleAPIErr(cudaMemcpy(d_y, y, sizeNodes, cudaMemcpyHostToDevice),
            //! __FILE__, __LINE__);

            HANDLE_API_ERR(cudaMalloc((void **)&d_z, sizeNodes));
            //! handleAPIErr(cudaMalloc((RDouble**) &d_z, sizeNodes), __FILE__, __LINE__);

            HANDLE_API_ERR(cudaMemcpy(d_z, z, sizeNodes, cudaMemcpyHostToDevice));
            //! handleAPIErr(cudaMemcpy(d_z, z, sizeNodes, cudaMemcpyHostToDevice),
            //! __FILE__, __LINE__);
        }
        //! copy xfn xfc xcc ... to gpu device
        void GPUFaceNormAllocCopy(RDouble *xfn, RDouble *yfn, RDouble *zfn, const int n)
        {
            size_t sizeFaces = n * sizeof(RDouble);
            HANDLE_API_ERR(cudaMalloc((void **)&d_xfn, sizeFaces));
            HANDLE_API_ERR(cudaMemcpy(d_xfn, xfn, sizeFaces, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMalloc((void **)&d_yfn, sizeFaces));
            HANDLE_API_ERR(cudaMemcpy(d_yfn, yfn, sizeFaces, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMalloc((void **)&d_zfn, sizeFaces));
            HANDLE_API_ERR(cudaMemcpy(d_zfn, zfn, sizeFaces, cudaMemcpyHostToDevice));
        }
        void GPUFaceCentAllocCopy(RDouble *xfc, RDouble *yfc, RDouble *zfc, const int n)
        {
            size_t sizeFaces = n * sizeof(RDouble);
            HANDLE_API_ERR(cudaMalloc((void **)&d_xfc, sizeFaces));
            HANDLE_API_ERR(cudaMemcpy(d_xfc, xfc, sizeFaces, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMalloc((void **)&d_yfc, sizeFaces));
            HANDLE_API_ERR(cudaMemcpy(d_yfc, yfc, sizeFaces, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMalloc((void **)&d_zfc, sizeFaces));
            HANDLE_API_ERR(cudaMemcpy(d_zfc, zfc, sizeFaces, cudaMemcpyHostToDevice));
        }
        void GPUCellCentAllocCopy(RDouble *xcc, RDouble *ycc, RDouble *zcc, const int n)
        {
            size_t sizeCells = n * sizeof(RDouble);

            HANDLE_API_ERR(cudaMalloc((void **)&d_xcc, sizeCells));
            //! handleAPIErr(cudaMalloc((RDouble**) &d_xcc, sizeCells), __FILE__,
            //! __LINE__);

            HANDLE_API_ERR(cudaMemcpy(d_xcc, xcc, sizeCells, cudaMemcpyHostToDevice));
            //! handleAPIErr(cudaMemcpy(d_xcc, xcc, sizeCells, cudaMemcpyHostToDevice),
            //! __FILE__, __LINE__);

            HANDLE_API_ERR(cudaMalloc((void **)&d_ycc, sizeCells));
            //! handleAPIErr(cudaMalloc((RDouble**) &d_ycc, sizeCells), __FILE__,
            //! __LINE__);

            HANDLE_API_ERR(cudaMemcpy(d_ycc, ycc, sizeCells, cudaMemcpyHostToDevice));
            //! handleAPIErr(cudaMemcpy(d_ycc, ycc, sizeCells, cudaMemcpyHostToDevice),
            //! __FILE__, __LINE__);

            HANDLE_API_ERR(cudaMalloc((void **)&d_zcc, sizeCells));
            //! handleAPIErr(cudaMalloc((RDouble**) &d_zcc, sizeCells), __FILE__,
            //! __LINE__);

            HANDLE_API_ERR(cudaMemcpy(d_zcc, zcc, sizeCells, cudaMemcpyHostToDevice));
            //! handleAPIErr(cudaMemcpy(d_zcc, zcc, sizeCells, cudaMemcpyHostToDevice),
            //! __FILE__, __LINE__);
        }
        void GPUAreaVolmAllocCopy(RDouble *area, RDouble *vol, const int nArea, const int nVol)
        {
            size_t sizeArea = nArea * sizeof(RDouble);
            size_t sizeVol  = nVol * sizeof(RDouble);

            HANDLE_API_ERR(cudaMalloc((void **)&d_area, sizeArea));
            HANDLE_API_ERR(cudaMemcpy(d_area, area, sizeArea, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMalloc((void **)&d_vol, sizeVol));
            HANDLE_API_ERR(cudaMemcpy(d_vol, vol, sizeVol, cudaMemcpyHostToDevice));
        }
        void GPUCellSkewnessAllocCopy(RDouble *cellSkewness, const int nTotalCell)
        {
            size_t sizecellSkewness = nTotalCell * sizeof(RDouble);
            HANDLE_API_ERR(cudaMalloc((void **)&d_cellSkewness, sizecellSkewness));
            HANDLE_API_ERR(cudaMemcpy(d_cellSkewness, cellSkewness, sizecellSkewness, cudaMemcpyHostToDevice));
        }
        void GPUFaceCellRelAllocCopy(const int *h_left_cell_of_face, const int *h_right_cell_of_face,
                                     const int nTotalFace)
        {
            size_t sizeFaces = nTotalFace * sizeof(int);

            HANDLE_API_ERR(cudaMalloc((void **)&d_left_cell_of_face, sizeFaces));
            //! handleAPIErr(cudaMalloc((int**) &d_left_cell_of_face, sizeFaces), __FILE__,
            //! __LINE__);

            HANDLE_API_ERR(cudaMemcpy(d_left_cell_of_face, h_left_cell_of_face, sizeFaces, cudaMemcpyHostToDevice));
            //! handleAPIErr(cudaMemcpy(d_left_cell_of_face, h_left_cell_of_face,
            //! sizeFaces, cudaMemcpyHostToDevice), __FILE__, __LINE__);

            HANDLE_API_ERR(cudaMalloc((void **)&d_right_cell_of_face, sizeFaces));
            //! handleAPIErr(cudaMalloc((int**) &d_right_cell_of_face, sizeFaces),
            //! __FILE__, __LINE__);

            HANDLE_API_ERR(cudaMemcpy(d_right_cell_of_face, h_right_cell_of_face, sizeFaces, cudaMemcpyHostToDevice));
            //! handleAPIErr(cudaMemcpy(d_right_cell_of_face, h_right_cell_of_face,
            //! sizeFaces, cudaMemcpyHostToDevice), __FILE__, __LINE__);
        }
        //! copy face2node, facetocell relationship to gpu
        //! void GPUFaceRelaAllocCopy(int *left_cell_of_face, int *right_cell_of_face int
        //! *face2node, int *node_number_of_each_face, const int n)
        //!{}
        void GPUFaceNodeRelAllocCopy(const int *h_face2node, const int *h_node_number_of_each_face, int nTotalFace2Node,
                                     int nTotalFace)
        {
            size_t sizeFaces = nTotalFace * sizeof(int);

            HANDLE_API_ERR(cudaMalloc((void **)&d_node_number_of_each_face, sizeFaces));
            //! handleAPIErr(cudaMalloc((int**) &d_node_number_of_each_face, sizeFaces),
            //! __FILE__, __LINE__);

            HANDLE_API_ERR(
                cudaMemcpy(d_node_number_of_each_face, h_node_number_of_each_face, sizeFaces, cudaMemcpyHostToDevice));
            //! handleAPIErr(cudaMemcpy(d_node_number_of_each_face,
            //! h_node_number_of_each_face, sizeFaces, cudaMemcpyHostToDevice), __FILE__,
            //! __LINE__);

            size_t sizeface2node = nTotalFace2Node * sizeof(int);

            HANDLE_API_ERR(cudaMalloc((void **)&d_face2node, sizeface2node));
            //! handleAPIErr(cudaMalloc((int**) &d_face2node, sizeface2node), __FILE__,
            //! __LINE__);

            HANDLE_API_ERR(cudaMemcpy(d_face2node, h_face2node, sizeface2node, cudaMemcpyHostToDevice));
            //! handleAPIErr(cudaMemcpy(d_face2node, h_face2node, sizeface2node,
            //! cudaMemcpyHostToDevice), __FILE__, __LINE__);
        }

        void GPUNodePosiFaceAllocCopy(const long long int *h_nodePosiFace, const int nTotalNodePosiFace)
        {
            size_t sizeNodePosiFace = nTotalNodePosiFace * sizeof(long long int);
            HANDLE_API_ERR(cudaMalloc((void **)&d_nodePosiFace, sizeNodePosiFace));
            HANDLE_API_ERR(cudaMemcpy(d_nodePosiFace, h_nodePosiFace, sizeNodePosiFace, cudaMemcpyHostToDevice));
        }

        void GPUBoundaryTypeAllocCopy(const int *h_boundaryType, const int nBoundFace)
        {
            size_t sizeBoundaryType = nBoundFace * sizeof(int);
            HANDLE_API_ERR(cudaMalloc((void **)&d_boundaryType, sizeBoundaryType));
            HANDLE_API_ERR(cudaMemcpy(d_boundaryType, h_boundaryType, sizeBoundaryType, cudaMemcpyHostToDevice));
        }

        void GPUFaceNormalAllocCopy(const RDouble *h_nxs, const RDouble *h_nys, const RDouble *h_nzs, int nTotalFace)
        {
            size_t sizeFaceNormal = nTotalFace * sizeof(RDouble);
            HANDLE_API_ERR(cudaMalloc((void **)&d_nxs, sizeFaceNormal));
            HANDLE_API_ERR(cudaMemcpy(d_nxs, h_nxs, sizeFaceNormal, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMalloc((void **)&d_nys, sizeFaceNormal));
            HANDLE_API_ERR(cudaMemcpy(d_nys, h_nys, sizeFaceNormal, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMalloc((void **)&d_nzs, sizeFaceNormal));
            HANDLE_API_ERR(cudaMemcpy(d_nzs, h_nzs, sizeFaceNormal, cudaMemcpyHostToDevice));
        }

        void GPUFaceCellRelMemcpy(int *face_number_of_each_cell, int **cell2face, int nTotalCell)
        {
            //! calculate the acc_face_number_of_cell
            int *acc_face_number_of_cell = new int[nTotalCell];
            acc_face_number_of_cell[0]   = 0;
            for (int i = 1; i < nTotalCell; i++)
                acc_face_number_of_cell[i] = acc_face_number_of_cell[i - 1] + face_number_of_each_cell[i - 1];
            //! copy face_number_of_each_cell and acc_face_number_of_cell
            size_t cSize = sizeof(int) * nTotalCell;
            HANDLE_API_ERR(cudaMalloc((void **)&d_face_number_of_each_cell, cSize));
            HANDLE_API_ERR(
                cudaMemcpy(d_face_number_of_each_cell, face_number_of_each_cell, cSize, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMalloc((void **)&d_acc_face_number_of_cell, cSize));
            HANDLE_API_ERR(
                cudaMemcpy(d_acc_face_number_of_cell, acc_face_number_of_cell, cSize, cudaMemcpyHostToDevice));
            //! copy cell2face
            size_t c2fSize = sizeof(int) * acc_face_number_of_cell[nTotalCell - 1];
            HANDLE_API_ERR(cudaMalloc((void **)&d_cell2face, c2fSize));
            HANDLE_API_ERR(cudaMemcpy(d_cell2face, cell2face[0], c2fSize, cudaMemcpyHostToDevice));
#ifdef CUDAUNITTEST
            TestFaceCellRel(nTotalCell, face_number_of_each_cell, cell2face, acc_face_number_of_cell);
#endif

            delete[] acc_face_number_of_cell;
        }

        void GPUCellFaceAllocCopy(int **cell2face, const int *face_number_of_each_cell, const int nTotalCell)
        {
            int    sizeCell2Face     = 0;
            size_t sizePosiCell2Face = nTotalCell * sizeof(int);
            int   *h_posiCell2Face   = new int[nTotalCell];
            for (int i = 0; i < nTotalCell; i++)
            {
                //! get the total number of face label from face_number_of_each_cell
                //! set the position of each cell in Cell2Face arranged in 1D
                h_posiCell2Face[i] = sizeCell2Face;
                sizeCell2Face += face_number_of_each_cell[i];
            }
            //! sizeCell2Face--;
            size_t numCell2Face       = sizeCell2Face * sizeof(int);
            int   *h_cell2face        = new int[sizeCell2Face];
            int    iFlagPosiCell2Face = 0;
            for (int i = 0; i < nTotalCell; i++)
            {
                for (int j = 0; j < face_number_of_each_cell[i]; j++)
                {
                    h_cell2face[iFlagPosiCell2Face] = cell2face[i][j];
                    iFlagPosiCell2Face++;
                }
            }
            HANDLE_API_ERR(cudaMalloc((void **)&d_cell2Face, numCell2Face));
            HANDLE_API_ERR(cudaMemcpy(d_cell2Face, h_cell2face, numCell2Face, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMalloc((void **)&d_posiCell2Face, sizePosiCell2Face));
            HANDLE_API_ERR(cudaMemcpy(d_posiCell2Face, h_posiCell2Face, sizePosiCell2Face, cudaMemcpyHostToDevice));
            //! comment by sunxu, 20191218
            //! HANDLE_API_ERR(cudaMalloc((void **)&d_face_number_of_each_cell,
            //! sizePosiCell2Face)); HANDLE_API_ERR(cudaMemcpy(d_face_number_of_each_cell,
            //! face_number_of_each_cell, sizePosiCell2Face, cudaMemcpyHostToDevice));
            delete[] h_posiCell2Face;
            delete[] h_cell2face;
        }
        //! Alloc device variables and transfer data  for interpolation points
        void GPUInterPointsAllocCopy(int nIPoint, const int *interPoint2GlobalPoint, const int *cellNumberOfInterPoint,
                                     const int *labelOfInterPoint)
        {
            size_t intPointSize = nIPoint * sizeof(int);
            HANDLE_API_ERR(cudaMalloc((void **)&d_interPoint2GlobalPoint, intPointSize));
            HANDLE_API_ERR(cudaMalloc((void **)&d_cellNumberOfInterPoint, intPointSize));
            HANDLE_API_ERR(cudaMalloc((void **)&d_labelOfInterPoint, intPointSize));
            HANDLE_API_ERR(
                cudaMemcpy(d_interPoint2GlobalPoint, interPoint2GlobalPoint, intPointSize, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(
                cudaMemcpy(d_cellNumberOfInterPoint, cellNumberOfInterPoint, intPointSize, cudaMemcpyHostToDevice));
            HANDLE_API_ERR(cudaMemcpy(d_labelOfInterPoint, labelOfInterPoint, intPointSize, cudaMemcpyHostToDevice));
        }

        void GPUInterFace2BoundaryFace(const int nIFace, const int *interFace2BoundaryFace)
        {
            size_t sizeIF2BF = nIFace * sizeof(int);
            HANDLE_API_ERR(cudaMalloc((void **)&d_interFace2BoundaryFace, sizeIF2BF));
            HANDLE_API_ERR(
                cudaMemcpy(d_interFace2BoundaryFace, interFace2BoundaryFace, sizeIF2BF, cudaMemcpyHostToDevice));
        }
#ifdef MCLUSGS
        void GPUCellColorAlloc()
        {
            using namespace CellColor;
            using namespace PHMPI;
            int nLocalZones = GetNumberofLocalZones();
            for (int izone = 0; izone < nLocalZones; izone++)
            {
                int           zoneID     = GetLocalZoneIDToGlobalZoneID(izone);
                int           level      = 0;
                UnstructGrid *grid       = UnstructGridCast(GetGrid(zoneID, level));
                const int     nTotalCell = grid->GetNTotalCell();
                size_t        sizeColor  = nTotalCell * sizeof(int);

                HANDLE_API_ERR(cudaMalloc((void **)&d_InteriorCellGroup, sizeColor));

                HANDLE_API_ERR(cudaMemcpy(d_InteriorCellGroup, InteriorCellGroup, sizeColor, cudaMemcpyHostToDevice));
            }
        }
#endif
    } //! namespace GPUGeomVariables
    //! free the memory of GPU variables
    void GPUMemoryFree()
    {
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
#ifdef CUDAUNITTEST
        using namespace GPUTestSpace;
#endif
        //! check device pointer
        //! TestGPUPointer(d_x);

        HANDLE_API_ERR(cudaFree(d_x));
        HANDLE_API_ERR(cudaFree(d_y));
        HANDLE_API_ERR(cudaFree(d_z));
        //! check device pointer after cudaFree
        //! TestGPUPointer(d_x);

        HANDLE_API_ERR(cudaFree(d_xfn));
        HANDLE_API_ERR(cudaFree(d_yfn));
        HANDLE_API_ERR(cudaFree(d_zfn));

        HANDLE_API_ERR(cudaFree(d_xcc));
        HANDLE_API_ERR(cudaFree(d_ycc));
        HANDLE_API_ERR(cudaFree(d_zcc));

        HANDLE_API_ERR(cudaFree(d_xfc));
        HANDLE_API_ERR(cudaFree(d_yfc));
        HANDLE_API_ERR(cudaFree(d_zfc));
        HANDLE_API_ERR(cudaFree(d_area));
        HANDLE_API_ERR(cudaFree(d_vol));
        HANDLE_API_ERR(cudaFree(d_cellSkewness));

        HANDLE_API_ERR(cudaFree(d_left_cell_of_face));
        HANDLE_API_ERR(cudaFree(d_right_cell_of_face));
        HANDLE_API_ERR(cudaFree(d_face2node));
        HANDLE_API_ERR(cudaFree(d_node_number_of_each_face));
        //! flow variables
        HANDLE_API_ERR(cudaFree(d_q_proxy));
        HANDLE_API_ERR(cudaFree(d_dqdx_proxy));
        HANDLE_API_ERR(cudaFree(d_dqdy_proxy));
        HANDLE_API_ERR(cudaFree(d_dqdz_proxy));
        HANDLE_API_ERR(cudaFree(d_bdqx));
        HANDLE_API_ERR(cudaFree(d_bdqy));
        HANDLE_API_ERR(cudaFree(d_bdqz));
        HANDLE_API_ERR(cudaFree(d_t_proxy));
        HANDLE_API_ERR(cudaFree(d_dtdx_proxy));
        HANDLE_API_ERR(cudaFree(d_dtdy_proxy));
        HANDLE_API_ERR(cudaFree(d_dtdz_proxy));
        HANDLE_API_ERR(cudaFree(d_q_turb_proxy));
        HANDLE_API_ERR(cudaFree(d_dq_turbdx_proxy));
        HANDLE_API_ERR(cudaFree(d_dq_turbdy_proxy));
        HANDLE_API_ERR(cudaFree(d_dq_turbdz_proxy));
        //! HANDLE_API_ERR(cudaFree(d_vel_proxy));
        HANDLE_API_ERR(cudaFree(d_dveldx_proxy));
        HANDLE_API_ERR(cudaFree(d_dveldy_proxy));
        HANDLE_API_ERR(cudaFree(d_dveldz_proxy));
        //! cancelled due to new method GGNODENEW
        //! HANDLE_API_ERR(cudaFree(d_q_n_tmp));
        //! HANDLE_API_ERR(cudaFree(d_n_count_tmp));
        HANDLE_API_ERR(cudaFree(d_visl));
        HANDLE_API_ERR(cudaFree(d_vist));
        HANDLE_API_ERR(cudaFree(d_flux));
        HANDLE_API_ERR(cudaFree(d_deltl));
        HANDLE_API_ERR(cudaFree(d_deltr));
        HANDLE_API_ERR(cudaFree(d_prim));
        HANDLE_API_ERR(cudaFree(d_tm));
        HANDLE_API_ERR(cudaFree(d_kcp));
        HANDLE_API_ERR(cudaFree(d_mul));
        HANDLE_API_ERR(cudaFree(d_mut));
        HANDLE_API_ERR(cudaFree(d_rho_ds));
        HANDLE_API_ERR(cudaFree(d_hint_s));
        HANDLE_API_ERR(cudaFree(d_q_ns));
        HANDLE_API_ERR(cudaFree(d_ql_ns));
        HANDLE_API_ERR(cudaFree(d_qr_ns));
        //! flow variables gama
        HANDLE_API_ERR(cudaFree(d_gama_ns));
        HANDLE_API_ERR(cudaFree(d_gamaL_ns));
        HANDLE_API_ERR(cudaFree(d_gamaR_ns));
        //! limiter variables
        HANDLE_API_ERR(cudaFree(d_limit));
        HANDLE_API_ERR(cudaFree(d_LIMIT));
        //! vgn FaceNormalVelocity
        HANDLE_API_ERR(cudaFree(d_vgn));
        //! xtn ytn ztn NodeVelocity
        HANDLE_API_ERR(cudaFree(d_xtn));
        HANDLE_API_ERR(cudaFree(d_ytn));
        HANDLE_API_ERR(cudaFree(d_ztn));
        //! boundary type
        HANDLE_API_ERR(cudaFree(d_boundaryType));
        //! flux
        //! res
        HANDLE_API_ERR(cudaFree(d_res_ns));
        HANDLE_API_ERR(cudaFree(d_dq_ns));
        HANDLE_API_ERR(cudaFree(d_dt));
        HANDLE_API_ERR(cudaFree(d_walldist));
        HANDLE_API_ERR(cudaFree(d_prim_inf));
        HANDLE_API_ERR(cudaFree(d_rhs_ns));
        HANDLE_API_ERR(cudaFree(d_invSpectrumRadius));
        HANDLE_API_ERR(cudaFree(d_visSpectrumRadius));
        HANDLE_API_ERR(cudaFree(d_dtv));
        HANDLE_API_ERR(cudaFree(d_minLocalTurbVis));
        HANDLE_API_ERR(cudaFree(d_maxLocalTurbVis));
        HANDLE_API_ERR(cudaFree(d_minDt));
        HANDLE_API_ERR(cudaFree(d_maxDt));
        HANDLE_API_ERR(cudaFree(d_minLocalTime));
        HANDLE_API_ERR(cudaFree(d_maxLocalTime));

        HANDLE_API_ERR(cudaFree(d_nxs));
        HANDLE_API_ERR(cudaFree(d_nys));
        HANDLE_API_ERR(cudaFree(d_nzs));
        HANDLE_API_ERR(cudaFree(d_qNode));
        HANDLE_API_ERR(cudaFree(d_tNode));
        HANDLE_API_ERR(cudaFree(d_nCount));
        HANDLE_API_ERR(cudaFree(d_nodeWeight));
        HANDLE_API_ERR(cudaFree(d_qTurbNode));
        HANDLE_API_ERR(cudaFree(d_qVelNode));
        //! remaining device vairalbes free
        HANDLE_API_ERR(cudaFree(d_nodePosiFace));
        HANDLE_API_ERR(cudaFree(d_cell2Face));
        HANDLE_API_ERR(cudaFree(d_cell2face));
        HANDLE_API_ERR(cudaFree(d_posiCell2Face));
        HANDLE_API_ERR(cudaFree(d_face_number_of_each_cell));
        HANDLE_API_ERR(cudaFree(d_acc_face_number_of_cell));
        if (d_interfaceValueExt)
        {
            //! delete interface information
            HANDLE_API_ERR(cudaFree(d_interFace2BoundaryFace));

            HANDLE_API_ERR(cudaFree(d_fg_send_q));
            HANDLE_API_ERR(cudaFree(d_fg_send_t));
            HANDLE_API_ERR(cudaFree(d_fg_recv_q));
            HANDLE_API_ERR(cudaFree(d_fg_recv_t));
            HANDLE_API_ERR(cudaFree(d_fg_send_dqdx));
            HANDLE_API_ERR(cudaFree(d_fg_recv_dqdx));
            HANDLE_API_ERR(cudaFree(d_fg_send_dqdy));
            HANDLE_API_ERR(cudaFree(d_fg_recv_dqdy));
            HANDLE_API_ERR(cudaFree(d_fg_send_dqdz));
            HANDLE_API_ERR(cudaFree(d_fg_recv_dqdz));
            HANDLE_API_ERR(cudaFree(d_fg_send_limit));
            HANDLE_API_ERR(cudaFree(d_fg_recv_limit));
            HANDLE_API_ERR(cudaFree(d_fg_send_qTurb));
            HANDLE_API_ERR(cudaFree(d_fg_recv_qTurb));
            HANDLE_API_ERR(cudaFree(d_fg_send_dtdx));
            HANDLE_API_ERR(cudaFree(d_fg_recv_dtdx));
            HANDLE_API_ERR(cudaFree(d_fg_send_dtdy));
            HANDLE_API_ERR(cudaFree(d_fg_recv_dtdy));
            HANDLE_API_ERR(cudaFree(d_fg_send_dtdz));
            HANDLE_API_ERR(cudaFree(d_fg_recv_dtdz));
            HANDLE_API_ERR(cudaFree(d_fg_send_dqTurbdx));
            HANDLE_API_ERR(cudaFree(d_fg_recv_dqTurbdx));
            HANDLE_API_ERR(cudaFree(d_fg_send_dqTurbdy));
            HANDLE_API_ERR(cudaFree(d_fg_recv_dqTurbdy));
            HANDLE_API_ERR(cudaFree(d_fg_send_dqTurbdz));
            HANDLE_API_ERR(cudaFree(d_fg_recv_dqTurbdz));
#ifdef MPIBUFFERLARGE
            HANDLE_API_ERR(cudaFreeHost(HostSendBufferNS));
            HANDLE_API_ERR(cudaFreeHost(HostRecvBufferNS));
            HANDLE_API_ERR(cudaFree(GPUSendBufferNS));
            HANDLE_API_ERR(cudaFree(GPURecvBufferNS));
            HANDLE_API_ERR(cudaFreeHost(HostSendBufferTurb));
            HANDLE_API_ERR(cudaFreeHost(HostRecvBufferTurb));
            HANDLE_API_ERR(cudaFree(GPUSendBufferTurb));
            HANDLE_API_ERR(cudaFree(GPURecvBufferTurb));
#ifdef HOSTDEVOVERLAP
            HANDLE_API_ERR(cudaEventDestroy(downloadDToDNS_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadDToDTurb_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDNSq_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDNSt_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDNSdqdx_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDNSdqdy_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDNSdqdz_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDNSlimit_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDNSdtdx_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDNSdtdy_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDNSdtdz_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDTurbq_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDTurbdqdx_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDTurbdqdy_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDTurbdqdz_done));
            HANDLE_API_ERR(cudaStreamDestroy(dataTransferDToH));
            HANDLE_API_ERR(cudaStreamDestroy(dataTransferHToD));
#endif
#endif
        }
        if (d_interpointValueExt)
        {
            //! delete interpolation point information
            HANDLE_API_ERR(cudaFree(d_interPoint2GlobalPoint));
            HANDLE_API_ERR(cudaFree(d_cellNumberOfInterPoint));
            HANDLE_API_ERR(cudaFree(d_labelOfInterPoint));
            HANDLE_API_ERR(cudaFree(d_qInterPoint));
            HANDLE_API_ERR(cudaFree(d_tInterPoint));
            HANDLE_API_ERR(cudaFree(d_qTurbInterPoint));

            HANDLE_API_ERR(cudaFree(d_fg_send_qNode));
            HANDLE_API_ERR(cudaFree(d_fg_send_tNode));
            HANDLE_API_ERR(cudaFree(d_fg_recv_qNode));
            HANDLE_API_ERR(cudaFree(d_fg_recv_tNode));
#ifdef MPIBUFFERLARGE
            HANDLE_API_ERR(cudaFreeHost(HostSendInterpointBufferNS));
            HANDLE_API_ERR(cudaFreeHost(HostRecvInterpointBufferNS));
            HANDLE_API_ERR(cudaFree(GPUSendInterpointBufferNS));
            HANDLE_API_ERR(cudaFree(GPURecvInterpointBufferNS));
            HANDLE_API_ERR(cudaFreeHost(HostSendInterpointBufferTurb));
            HANDLE_API_ERR(cudaFreeHost(HostRecvInterpointBufferTurb));
            HANDLE_API_ERR(cudaFree(GPUSendInterpointBufferTurb));
            HANDLE_API_ERR(cudaFree(GPURecvInterpointBufferTurb));
#ifdef HOSTDEVOVERLAP
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDNSqNode_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDNStNode_done));
            HANDLE_API_ERR(cudaEventDestroy(downloadHToDTurbqNode_done));
            HANDLE_API_ERR(cudaStreamDestroy(dataTransferDToHInterpoint));
            HANDLE_API_ERR(cudaStreamDestroy(dataTransferHToDInterpoint));
#endif
#endif
        }
        //! deleting device variables for GPUFaceColor
        using namespace GPUFaceColor;
        if (d_BoundFaceGroup) HANDLE_API_ERR(cudaFree(d_BoundFaceGroup));
        if (d_InteriorFaceGroup) HANDLE_API_ERR(cudaFree(d_InteriorFaceGroup));
        //! reset the gpu device
        //! cudaDeviceReset();
    }

} //! namespace GPUMemory
