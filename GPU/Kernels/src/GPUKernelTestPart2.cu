#include "GPUKernelTestPart2.h"
#include "GPUKernelTest.h"
#include "BasicDeviceVariables.h"
#include "Geo_SimpleBC.h"
#include "OutputDebug.h"

using namespace std;
void TestGPUCompNodeVarByGradientLoop1(const string d_gradient_field_proxy, const int nTotalNode, const RFloat *org_q_n,
                                       const int *org_n_count)
{
#ifdef UNITTESTOUTPUT
    cout << "Test GPUNcountQnInit ..." << endl;
#endif
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUControlVariables;
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    RFloat *h_q_n = new RFloat[nTotalNode];
    //!double *h_q_n = new double[nTotalNode];
    int *h_n_count = new int[nTotalNode];

    size_t sizeNodeRFloat = nTotalNode * sizeof(RFloat);
    //!size_t sizeNodeRFloat = nTotalNode * sizeof(double);
    //!size_t sizeNodeDouble = nTotalNode * sizeof(double);
    HANDLE_API_ERR(cudaMemcpy(h_q_n, d_q_n_tmp, sizeNodeRFloat, cudaMemcpyDeviceToHost));
    //!HANDLE_API_ERR(cudaMemcpy(h_q_n, d_q_n_double, sizeNodeDouble, cudaMemcpyDeviceToHost));
    size_t sizeNodeInt = nTotalNode * sizeof(int);
    HANDLE_API_ERR(cudaMemcpy(h_n_count, d_n_count_tmp, sizeNodeInt, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotalNode; i++)
    {
        abErr           = fabs(h_q_n[i] - org_q_n[i]);
        rtErr           = abErr / Max(fabs(org_q_n[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error: In TestGPUCompNodeVarByGradientLoop1 difference exists in %d term: "
                   "d_q_n=%e\t org_q_n=%e\t abErr= %e\t rtErr= %e\n",
                   i, h_q_n[i], org_q_n[i], abErr, rtErr);
            exit(1);
        }
        //!if ( (fabs(h_n_count[i]-org_n_count[i]) > flagDiff) || (fabs(h_q_n[i]-org_q_n[i]) > flagDiff) ) {
        //!if ( (fabs(h_n_count[i]-org_n_count[i]) > flagDiff) || (abErr>abCompare && rtErr>rtCompare)) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Error: In TestGPUCompNodeVarByGradientLoop1 difference exists in %d term: "
                   "d_q_n=%e\t org_q_n=%e\t abErr= %e\t rtErr= %e\n",
                   i, h_q_n[i], org_q_n[i], abErr, rtErr);
            exit(1);
        }
        abErr    = fabs(h_n_count[i] - org_n_count[i]);
        rtErr    = abErr / Max(fabs(org_n_count[i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error: In TestGPUCompNodeVarByGradientLoop1 difference exists in %d term: "
                   "d_n_count=%d, org_n_count=%d, abErr= %e\t rtErr= %e\n",
                   i, h_n_count[i], org_n_count[i], abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Error: In TestGPUCompNodeVarByGradientLoop1 difference exists in %d term: "
                   "d_n_count=%d, org_n_count=%d, abErr= %e\t rtErr= %e\n",
                   i, h_n_count[i], org_n_count[i], abErr, rtErr);
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUKernelCompNodeVarByGradientLoop1 successfully" << endl;
#endif
    delete[] h_q_n;
    delete[] h_n_count;
}

void TestGPUNodePosiFaceAllocCopy(const long long int *org_nodePosiFace, const long long int *d_nodePosiFace,
                                  const int nTotalNodePosiFace)
{
    using namespace GPUMemory;
    using namespace GPUGeomVariables;
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    long long int *h_nodePosiFace;
    h_nodePosiFace          = new long long int[nTotalNodePosiFace];
    size_t sizeNodePosiFace = nTotalNodePosiFace * sizeof(long long int);
    HANDLE_API_ERR(cudaMemcpy(h_nodePosiFace, d_nodePosiFace, sizeNodePosiFace, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotalNodePosiFace; i++)
    {
        abErr           = fabs(h_nodePosiFace[i] - org_nodePosiFace[i]);
        rtErr           = abErr / Max(fabs(org_nodePosiFace[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            cout << "Something wrong exists in memory allocate and copy for d_nodePosiFace" << endl;
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            cout << "Something wrong exists in memory allocate and copy for d_nodePosiFace" << endl;
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Allocate and copy d_nodePosiFace successfully" << endl;
#endif
    delete[] h_nodePosiFace;
}

void TestGPUBoundaryTypeAllocCopy(const int *org_boundaryType, const int *d_boundaryType, const int nBoundFace)
{
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t sizeBoundaryType = nBoundFace * sizeof(int);
    int   *h_boundaryType   = new int[nBoundFace];
    HANDLE_API_ERR(cudaMemcpy(h_boundaryType, d_boundaryType, sizeBoundaryType, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nBoundFace; i++)
    {
        /* for muli process, it comes out!!
                //!By output, no FANTASY::INTERFACE exists
                if (org_boundaryType[i] == FANTASY::INTERFACE) {
                    cout<<"FANTASY::INTERFACE exists in boundary face"<<i<<endl;
                    exit(1);
                }
                */
        abErr           = fabs(org_boundaryType[i] - h_boundaryType[i]);
        rtErr           = abErr / Max(fabs(h_boundaryType[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            cout << "Something wrong exists in memory allocate and copy for d_boundaryType" << endl;
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            cout << "Something wrong exists in memory allocate and copy for d_boundaryType" << endl;
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Allocate and copy d_boundaryType successfully" << endl;
#endif
    delete[] h_boundaryType;
}

void TestGPUCompNodeVarByGradientLoop2(const string d_gradient_field_proxy, const int nTotalNode, const RFloat *org_q_n,
                                       const int *org_n_count)
{
    //!Just a temporary treatment, because the other condition is not considered here
#ifdef UNITTESTOUTPUT
    cout << "Test GPUBoundaryFaceNodeCal ... " << endl;
#endif
    if (d_gradient_field_proxy != "d_q_proxy") return;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUControlVariables;
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    //!double *h_q_n = new double[nTotalNode];
    RFloat *h_q_n     = new RFloat[nTotalNode];
    int    *h_n_count = new int[nTotalNode];

    //!size_t sizeNodeDouble = nTotalNode * sizeof(double);
    size_t sizeNodeRFloat = nTotalNode * sizeof(RFloat);
    //!HANDLE_API_ERR(cudaMemcpy(h_q_n, d_q_n_double, sizeNodeDouble, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_q_n, d_q_n_tmp, sizeNodeRFloat, cudaMemcpyDeviceToHost));
    size_t sizeNodeInt = nTotalNode * sizeof(int);
    HANDLE_API_ERR(cudaMemcpy(h_n_count, d_n_count_tmp, sizeNodeInt, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotalNode; i++)
    {
        abErr           = fabs(h_n_count[i] - org_n_count[i]);
        rtErr           = abErr / Max(fabs(org_n_count[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error:  In TestGPUCompNodeVarByGradientLoop2, difference exists in %d term: "
                   "d_n_count=%d, org_n_count=%d\n",
                   i, h_n_count[i], org_n_count[i]);
            exit(1);
        }
        //!if  (abs(h_n_count[i]-org_n_count[i]) > flagDiff) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Error:  In TestGPUCompNodeVarByGradientLoop2, difference exists in %d term: "
                   "d_n_count=%d, org_n_count=%d\n",
                   i, h_n_count[i], org_n_count[i]);
            exit(1);
        }
        abErr    = fabs(h_q_n[i] - org_q_n[i]);
        rtErr    = abErr / Max(fabs(org_q_n[i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error:  In TestGPUCompNodeVarByGradientLoop2, difference exists in %d term: "
                   "d_n_count=%d, org_n_count=%d\n",
                   i, h_n_count[i], org_n_count[i]);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Error: In TestGPUCompNodeVarByGradientLoop2, difference exists in %d term: "
                   "d_q_n=%e\t, org_q_n=%e\t, abErr= %e\t rtErr= %e\n",
                   i, h_q_n[i], org_q_n[i], abErr, rtErr);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUKernelCompNodeVarByGradientLoop2 successfully" << endl;
#endif
    delete[] h_q_n;
    delete[] h_n_count;
}

void TestGPUCompNodeVarByGradientFinal(const string d_gradient_field_proxy, const int nTotalNode, RFloat *org_q_n)
{
    //!Just a temporary treatment, because the other condition is not considered here
    if ((d_gradient_field_proxy != "d_q_proxy") && (d_gradient_field_proxy != "d_t_proxy")
        && (d_gradient_field_proxy != "d_q_turb_proxy") && (d_gradient_field_proxy != "d_velocity_proxy"))
        return;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUControlVariables;
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    //!double *h_q_n = new double[nTotalNode];
    RFloat *h_q_n = new RFloat[nTotalNode];
    //!size_t sizeNodeDouble = nTotalNode * sizeof(double);
    size_t sizeNodeRFloat = nTotalNode * sizeof(RFloat);
    //!HANDLE_API_ERR(cudaMemcpy(h_q_n, d_q_n_double, sizeNodeDouble, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_q_n, d_q_n_tmp, sizeNodeRFloat, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotalNode; i++)
    {
        abErr           = fabs(h_q_n[i] - org_q_n[i]);
        rtErr           = abErr / Max(fabs(org_q_n[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error: In TestGPUCompNodeVarByGradientFinal, difference exists in %d term: "
                   "d_q_n=%e\t, org_q_n=%e\t, abErr= %e\t rtErr= %e\n",
                   i, h_q_n[i], org_q_n[i], abErr, rtErr);
            exit(1);
        }
        //!if (abs(h_q_n[i]-org_q_n[i])>flagCompare){
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Error: In TestGPUCompNodeVarByGradientFinal, difference exists in %d term: "
                   "d_q_n=%e\t, org_q_n=%e\t, abErr= %e\t rtErr= %e\n",
                   i, h_q_n[i], org_q_n[i], abErr, rtErr);
            exit(1);
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUKernelCompNodeVarByGradientFinal with d_gradient_field_proxy=" << d_gradient_field_proxy
         << " successfully" << endl;
#endif
    delete[] h_q_n;
}

void TestGPUCompGradientGGNode(const string q_name, const int index, const int nTotalCell, const int nBoundFace,
                               RFloat *org_dqdx, RFloat *org_dqdy, RFloat *org_dqdz)
{
    //!    cout<<"q_name="<<q_name<<endl;
    //!    if (( q_name != "d_q_proxy" )&&( q_name != "d_t_proxy" )) return;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    //!double flagCompare = 1.0e-10; //!By test: dqdx[icell] /= vol[icell], dividing vol induce large error.
    RFloat  abErr  = 0.0;
    RFloat  rtErr  = 0.0;
    int     nTotal = nTotalCell + nBoundFace;
    size_t  sizedQ;
    RFloat *h_dqdx;
    RFloat *h_dqdy;
    RFloat *h_dqdz;
    if (q_name == "d_q_proxy")
    {
        sizedQ = nTotal * d_neqn_ns * sizeof(RFloat);
        h_dqdx = new RFloat[nTotal * d_neqn_ns];
        HANDLE_API_ERR(cudaMemcpy(h_dqdx, d_dqdx_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdy = new RFloat[nTotal * d_neqn_ns];
        HANDLE_API_ERR(cudaMemcpy(h_dqdy, d_dqdy_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdz = new RFloat[nTotal * d_neqn_ns];
        HANDLE_API_ERR(cudaMemcpy(h_dqdz, d_dqdz_proxy, sizedQ, cudaMemcpyDeviceToHost));
    }
    else if (q_name == "d_t_proxy")
    {
        sizedQ = nTotal * 1 * sizeof(RFloat);
        h_dqdx = new RFloat[nTotal * 1];
        HANDLE_API_ERR(cudaMemcpy(h_dqdx, d_dtdx_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdy = new RFloat[nTotal * 1];
        HANDLE_API_ERR(cudaMemcpy(h_dqdy, d_dtdy_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdz = new RFloat[nTotal * 1];
        HANDLE_API_ERR(cudaMemcpy(h_dqdz, d_dtdz_proxy, sizedQ, cudaMemcpyDeviceToHost));
    }
    else if (q_name == "d_q_turb_proxy")
    {
        sizedQ = nTotal * d_neqn_turb * sizeof(RFloat);
        h_dqdx = new RFloat[nTotal * d_neqn_turb];
        HANDLE_API_ERR(cudaMemcpy(h_dqdx, d_dq_turbdx_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdy = new RFloat[nTotal * 1];
        HANDLE_API_ERR(cudaMemcpy(h_dqdy, d_dq_turbdy_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdz = new RFloat[nTotal * 1];
        HANDLE_API_ERR(cudaMemcpy(h_dqdz, d_dq_turbdz_proxy, sizedQ, cudaMemcpyDeviceToHost));
    }
    else if (q_name == "d_velocity_proxy")
    {
        sizedQ = nTotal * 3 * sizeof(RFloat);
        h_dqdx = new RFloat[nTotal * 3];
        HANDLE_API_ERR(cudaMemcpy(h_dqdx, d_dveldx_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdy = new RFloat[nTotal * 3];
        HANDLE_API_ERR(cudaMemcpy(h_dqdy, d_dveldy_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdz = new RFloat[nTotal * 3];
        HANDLE_API_ERR(cudaMemcpy(h_dqdz, d_dveldz_proxy, sizedQ, cudaMemcpyDeviceToHost));
    }
    else
    {
        return;
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUCompGradientGGNode with q_name =" << q_name << endl;
#endif
    for (int i = 0; i < nTotal; i++)
    {
        abErr           = fabs(org_dqdx[i] - h_dqdx[index * nTotal + i]);
        rtErr           = abErr / Max(fabs(org_dqdx[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error in GPUCompGradientGGNode: %d term, h_dqdx=%e, d_dqdx_proxy=%e\t, abErr= "
                   "%e\t rtErr= %e\n",
                   i, org_dqdx[i], h_dqdx[index * nTotal + i], abErr, rtErr);

            exit(1);
        }
        //!if (fabs(org_dqdx[i]-h_dqdx[index*nTotal+i])>flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            if (processNum > 1)
                printf("Error in GPUCompGradientGGNode on globalRank %d, localRank %d: %d term, "
                       "h_dqdx=%e, d_dqdx_proxy=%e\t, abErr= %e\t rtErr= %e\n",
                       globalRank, localRank, i, org_dqdx[i], h_dqdx[index * nTotal + i], abErr, rtErr);
            else
                printf("Error in GPUCompGradientGGNode: %d term, h_dqdx=%e, d_dqdx_proxy=%e\t, "
                       "abErr= %e\t rtErr= %e\n",
                       i, org_dqdx[i], h_dqdx[index * nTotal + i], abErr, rtErr);

            exit(1);
        }
        abErr    = fabs(org_dqdy[i] - h_dqdy[index * nTotal + i]);
        rtErr    = abErr / Max(fabs(org_dqdy[i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error in GPUCompGradientGGNode: %d term, h_dqdy=%.30e, d_dqdy_proxy=%.30e\t, "
                   "abErr= %.30e\t rtErr= %.30e\n",
                   i, org_dqdy[i], h_dqdy[index * nTotal + i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_dqdy[i]-h_dqdy[index*nTotal+i]));
            exit(1);
        }
        //!if (fabs(org_dqdy[i]-h_dqdy[index*nTotal+i])>flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Error in GPUCompGradientGGNode: %d term, h_dqdy=%.30e, d_dqdy_proxy=%.30e\t, "
                   "abErr= %.30e\t rtErr= %.30e\n",
                   i, org_dqdy[i], h_dqdy[index * nTotal + i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_dqdy[i]-h_dqdy[index*nTotal+i]));
            exit(1);
        }

        abErr    = fabs(org_dqdz[i] - h_dqdz[index * nTotal + i]);
        rtErr    = abErr / Max(fabs(org_dqdz[i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error in GPUCompGradientGGNode: %d term, h_dqdz=%e, d_dqdz_proxy=%e\t, abErr= "
                   "%e\t rtErr= %e\n",
                   i, org_dqdz[i], h_dqdz[index * nTotal + i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_dqdz[i]-h_dqdz[index*nTotal+i]));
            exit(1);
        }
        //!if (fabs(org_dqdz[i]-h_dqdz[index*nTotal+i])>flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Error in GPUCompGradientGGNode: %d term, h_dqdz=%e, d_dqdz_proxy=%e\t, abErr= "
                   "%e\t rtErr= %e\n",
                   i, org_dqdz[i], h_dqdz[index * nTotal + i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_dqdz[i]-h_dqdz[index*nTotal+i]));
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUCompGradientGGNode with q_name =" << q_name << " successfully" << endl;
#endif
    delete[] h_dqdx;
    delete[] h_dqdy;
    delete[] h_dqdz;
}
void TestGPUCompGradientGGCell(const string q_name, const int index, const int nTotalCell, const int nBoundFace,
                               RFloat *org_dqdx, RFloat *org_dqdy, RFloat *org_dqdz)
{
    //!    cout<<"q_name="<<q_name<<endl;
    //!    if (( q_name != "d_q_proxy" )&&( q_name != "d_t_proxy" )) return;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat  abErr  = 0.0;
    RFloat  rtErr  = 0.0;
    int     nTotal = nTotalCell + nBoundFace;
    size_t  sizedQ;
    RFloat *h_dqdx;
    RFloat *h_dqdy;
    RFloat *h_dqdz;
    if (q_name == "d_q_proxy")
    {
        sizedQ = nTotal * d_neqn_ns * sizeof(RFloat);
        h_dqdx = new RFloat[nTotal * d_neqn_ns];
        HANDLE_API_ERR(cudaMemcpy(h_dqdx, d_dqdx_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdy = new RFloat[nTotal * d_neqn_ns];
        HANDLE_API_ERR(cudaMemcpy(h_dqdy, d_dqdy_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdz = new RFloat[nTotal * d_neqn_ns];
        HANDLE_API_ERR(cudaMemcpy(h_dqdz, d_dqdz_proxy, sizedQ, cudaMemcpyDeviceToHost));
    }
    else if (q_name == "d_t_proxy")
    {
        sizedQ = nTotal * 1 * sizeof(RFloat);
        h_dqdx = new RFloat[nTotal * 1];
        HANDLE_API_ERR(cudaMemcpy(h_dqdx, d_dtdx_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdy = new RFloat[nTotal * 1];
        HANDLE_API_ERR(cudaMemcpy(h_dqdy, d_dtdy_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdz = new RFloat[nTotal * 1];
        HANDLE_API_ERR(cudaMemcpy(h_dqdz, d_dtdz_proxy, sizedQ, cudaMemcpyDeviceToHost));
    }
    else if (q_name == "d_q_turb_proxy")
    {
        sizedQ = nTotal * d_neqn_turb * sizeof(RFloat);
        h_dqdx = new RFloat[nTotal * d_neqn_turb];
        HANDLE_API_ERR(cudaMemcpy(h_dqdx, d_dq_turbdx_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdy = new RFloat[nTotal * 1];
        HANDLE_API_ERR(cudaMemcpy(h_dqdy, d_dq_turbdy_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdz = new RFloat[nTotal * 1];
        HANDLE_API_ERR(cudaMemcpy(h_dqdz, d_dq_turbdz_proxy, sizedQ, cudaMemcpyDeviceToHost));
    }
    else if (q_name == "d_velocity_proxy")
    {
        sizedQ = nTotal * 3 * sizeof(RFloat);
        h_dqdx = new RFloat[nTotal * 3];
        HANDLE_API_ERR(cudaMemcpy(h_dqdx, d_dveldx_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdy = new RFloat[nTotal * 3];
        HANDLE_API_ERR(cudaMemcpy(h_dqdy, d_dveldy_proxy, sizedQ, cudaMemcpyDeviceToHost));
        h_dqdz = new RFloat[nTotal * 3];
        HANDLE_API_ERR(cudaMemcpy(h_dqdz, d_dveldz_proxy, sizedQ, cudaMemcpyDeviceToHost));
    }
    else
    {
        return;
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUCompGradientGGCell with q_name =" << q_name << endl;
#endif
    for (int i = 0; i < nTotal; i++)
    {
        abErr           = fabs(org_dqdx[i] - h_dqdx[index * nTotal + i]);
        rtErr           = abErr / Max(fabs(org_dqdx[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error in GPUCompGradientGGCell: %d term, h_dqdx=%e, d_dqdx_proxy=%e\t, abErr= "
                   "%e\t rtErr= %e\n",
                   i, org_dqdx[i], h_dqdx[index * nTotal + i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_dqdx[i]-h_dqdx[index*nTotal+i]));
        }
        //!if (fabs(org_dqdx[i]-h_dqdx[index*nTotal+i])>flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Error in GPUCompGradientGGCell: %d term, h_dqdx=%e, d_dqdx_proxy=%e\t, abErr= "
                   "%e\t rtErr= %e\n",
                   i, org_dqdx[i], h_dqdx[index * nTotal + i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_dqdx[i]-h_dqdx[index*nTotal+i]));

            exit(1);
        }
        abErr    = fabs(org_dqdy[i] - h_dqdy[index * nTotal + i]);
        rtErr    = abErr / Max(fabs(org_dqdy[i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error in GPUCompGradientGGCell: %d term, h_dqdy=%e, d_dqdy_proxy=%e\t, abErr= "
                   "%e\t rtErr= %e\n",
                   i, org_dqdy[i], h_dqdy[index * nTotal + i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_dqdy[i]-h_dqdy[index*nTotal+i]));
            exit(1);
        }
        //!if (fabs(org_dqdy[i]-h_dqdy[index*nTotal+i])>flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Error in GPUCompGradientGGCell: %d term, h_dqdy=%e, d_dqdy_proxy=%e\t, abErr= "
                   "%e\t rtErr= %e\n",
                   i, org_dqdy[i], h_dqdy[index * nTotal + i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_dqdy[i]-h_dqdy[index*nTotal+i]));
            exit(1);
        }

        abErr    = fabs(org_dqdz[i] - h_dqdz[index * nTotal + i]);
        rtErr    = abErr / Max(fabs(org_dqdz[i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error in GPUCompGradientGGCell: %d term, h_dqdz=%e, d_dqdz_proxy=%e\t, abErr= "
                   "%e\t rtErr= %e\n",
                   i, org_dqdz[i], h_dqdz[index * nTotal + i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_dqdz[i]-h_dqdz[index*nTotal+i]));
            exit(1);
        }
        //!if (fabs(org_dqdz[i]-h_dqdz[index*nTotal+i])>flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Error in GPUCompGradientGGCell: %d term, h_dqdz=%e, d_dqdz_proxy=%e\t, abErr= "
                   "%e\t rtErr= %e\n",
                   i, org_dqdz[i], h_dqdz[index * nTotal + i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_dqdz[i]-h_dqdz[index*nTotal+i]));
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUCompGradientGGCell with q_name =" << q_name << " successfully" << endl;
#endif
    delete[] h_dqdx;
    delete[] h_dqdy;
    delete[] h_dqdz;
}

void TestGPUStoreBoundGrad(const int nBoundFace, const int nTVar, RFloat **org_bdqx, RFloat **org_bdqy,
                           RFloat **org_bdqz)
{
    RFloat  abErr  = 0.0;
    RFloat  rtErr  = 0.0;
    RFloat *h_bdqx = new RFloat[nBoundFace * nTVar];
    RFloat *h_bdqy = new RFloat[nBoundFace * nTVar];
    RFloat *h_bdqz = new RFloat[nBoundFace * nTVar];
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    size_t sizeFaceRFloat = nBoundFace * nTVar * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(h_bdqx, d_bdqx, sizeFaceRFloat, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_bdqy, d_bdqy, sizeFaceRFloat, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_bdqz, d_bdqz, sizeFaceRFloat, cudaMemcpyDeviceToHost));

    /*        
        for (int i = 0; i < nBoundFace * nTVar; i++){
            if (abs(h_bdqx[i] - org_bdqx[0][i]) > flagCompare) {
                cout<<"something wrong in GPUStoreBoundGrad"<<endl;
                exit(1);
            }
        }
*/
    for (int i = 0; i < nBoundFace; i++)
    {
        for (int j = 0; j < nTVar; j++)
        {
            abErr           = fabs(h_bdqx[j * nBoundFace + i] - org_bdqx[j][i]);
            rtErr           = abErr / Max(fabs(org_bdqx[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error in GPUStoreBoundGrad, %d, %d term, bdqx on host: %.30e, on device: "
                       "%.30e, abErr= %.30e, rtErr= %.30e\n",
                       j, i, org_bdqx[j][i], h_bdqx[j * nBoundFace + i], abErr, rtErr);
                //!cout<<"something wrong in GPUStoreBoundGrad at bdqx"<<endl;
                exit(1);
            }
            //!if (abs(h_bdqx[j*nBoundFace+i] - org_bdqx[j][i]) > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error in GPUStoreBoundGrad, %d, %d term, bdqx on host: %.30e, on device: "
                       "%.30e, abErr= %.30e, rtErr= %.30e\n",
                       j, i, org_bdqx[j][i], h_bdqx[j * nBoundFace + i], abErr, rtErr);
                //!cout<<"something wrong in GPUStoreBoundGrad at bdqx"<<endl;
                exit(1);
            }
            abErr    = fabs(h_bdqy[j * nBoundFace + i] - org_bdqy[j][i]);
            rtErr    = abErr / Max(fabs(org_bdqy[j][i]), smallLimit);
            checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error in GPUStoreBoundGrad, %d, %d term, bdqy on host: %.30e, on device: "
                       "%.30e, abErr= %.30e, rtErr= %.30e\n",
                       j, i, org_bdqy[j][i], h_bdqy[j * nBoundFace + i], abErr, rtErr);
                exit(1);
            }
            //!if (abs(h_bdqy[j*nBoundFace+i] - org_bdqy[j][i]) > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                //!cout<<"something wrong in GPUStoreBoundGrad at bdqy"<<endl;
                printf("Error in GPUStoreBoundGrad, %d, %d term, bdqy on host: %.30e, on device: "
                       "%.30e, abErr= %.30e, rtErr= %.30e\n",
                       j, i, org_bdqy[j][i], h_bdqy[j * nBoundFace + i], abErr, rtErr);
                exit(1);
            }
            abErr    = fabs(h_bdqz[j * nBoundFace + i] - org_bdqz[j][i]);
            rtErr    = abErr / Max(fabs(org_bdqz[j][i]), smallLimit);
            checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error in GPUStoreBoundGrad, %d, %d term, bdqz on host: %e, on device: %e, "
                       "abErr= %e, rtErr= %e\n",
                       j, i, org_bdqz[j][i], h_bdqz[j * nBoundFace + i], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                //!if (abs(h_bdqz[j*nBoundFace+i] - org_bdqz[j][i]) > flagCompare) {
                //!cout<<"something wrong in GPUStoreBoundGrad at bdqz"<<endl;
                printf("Error in GPUStoreBoundGrad, %d, %d term, bdqz on host: %e, on device: %e, "
                       "abErr= %e, rtErr= %e\n",
                       j, i, org_bdqz[j][i], h_bdqz[j * nBoundFace + i], abErr, rtErr);
                exit(1);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUStoreBoundGrad successfully" << endl;
#endif

    delete[] h_bdqx;
    delete[] h_bdqy;
    delete[] h_bdqz;
}

void TestGPUFaceWeight(const int nst, const int ned, const RDouble *org_deltl, const RDouble *org_deltr)
{
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RDouble *h_deltl  = new RDouble[d_SEG_LEN];
    RDouble *h_deltr  = new RDouble[d_SEG_LEN];
    size_t   sizeDelt = d_SEG_LEN * sizeof(RDouble);
    HANDLE_API_ERR(cudaMemcpy(h_deltl, d_deltl, sizeDelt, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_deltr, d_deltr, sizeDelt, cudaMemcpyDeviceToHost));

    //!for(int i = nst; i < ned; i++){
    for (int i = 0; i < d_SEG_LEN; i++)
    {
        abErr           = fabs(org_deltl[i] - h_deltl[i]);
        rtErr           = abErr / Max(fabs(org_deltl[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in GPUFaceWeight in %d term, h_deltl=%.30e, d_deltl=%.30e, "
                   "abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_deltl[i], h_deltl[i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_deltl[i] - h_deltl[i]));
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            //!if (fabs(org_deltl[i] - h_deltl[i]) > flagCompare) {
            printf("Something wrong in GPUFaceWeight in %d term, h_deltl=%.30e, d_deltl=%.30e, "
                   "abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_deltl[i], h_deltl[i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_deltl[i] - h_deltl[i]));
            exit(1);
        }
        abErr    = fabs(org_deltr[i] - h_deltr[i]);
        rtErr    = abErr / Max(fabs(org_deltr[i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in GPUFaceWeight in %d term, h_deltr=%.30e, d_deltr=%.30e, "
                   "abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_deltr[i], h_deltr[i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_deltr[i] - h_deltr[i]));
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            //!if (fabs(org_deltr[i] - h_deltr[i]) > flagCompare) {
            printf("Something wrong in GPUFaceWeight in %d term, h_deltr=%.30e, d_deltr=%.30e, "
                   "abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_deltr[i], h_deltr[i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_deltr[i] - h_deltr[i]));
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUFaceWeight succefully" << endl;
#endif
    delete[] h_deltl;
    delete[] h_deltr;
}

void TestGPUViscousCoefCopy(const RFloat *org_visl, const RFloat *org_vist, const int nTotal)
{
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    size_t  sizeVisCoef = nTotal * sizeof(RFloat);
    RFloat *h_visl      = new RFloat[nTotal];
    RFloat *h_vist      = new RFloat[nTotal];
    HANDLE_API_ERR(cudaMemcpy(h_visl, d_visl, sizeVisCoef, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_vist, d_vist, sizeVisCoef, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; i++)
    {
        abErr           = fabs(h_visl[i] - org_visl[i]);
        rtErr           = abErr / Max(fabs(org_visl[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            cout << "Something wrong in  GPUViscousCoefCopy for visl" << endl;
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            cout << "Something wrong in  GPUViscousCoefCopy for visl" << endl;
            exit(1);
        }
        abErr    = fabs(h_vist[i] - org_vist[i]);
        rtErr    = abErr / Max(fabs(org_vist[i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            cout << "Something wrong in  GPUViscousCoefCopy for vist" << endl;
            exit(1);
        }
        //!if (fabs(h_vist[i] - org_vist[i]) > flagCompare){
        if (abErr > abCompare && rtErr > rtCompare)
        {
            cout << "Something wrong in  GPUViscousCoefCopy for vist" << endl;
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test copy visl to d_visl and vist to d_vist successfully" << endl;
#endif
    delete[] h_visl;
    delete[] h_vist;
}

void TestGPUGetVisFaceValueMult(const int nst, const int ned, const RFloat *org_mul, const RFloat *org_mut)
{
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat *h_mul = new RFloat[d_SEG_LEN];
    RFloat *h_mut = new RFloat[d_SEG_LEN];

    size_t sizeMult = d_SEG_LEN * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(h_mul, d_mul, sizeMult, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_mut, d_mut, sizeMult, cudaMemcpyDeviceToHost));

    for (int i = 0; i < d_SEG_LEN; i++)
    {
        abErr           = fabs(h_mul[i] - org_mul[i]);
        rtErr           = abErr / Max(fabs(org_mul[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in GPUGetFaceValue in %d term, h_mul = %.30e, d_mul = %.30e, "
                   "abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_mul[i], h_mul[i], abErr, rtErr);
            exit(1);
        }
        //!if ( fabs(h_mul[i] - org_mul[i]) > flagCompare ) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in GPUGetFaceValue in %d term, h_mul = %.30e, d_mul = %.30e, "
                   "abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_mul[i], h_mul[i], abErr, rtErr);
            exit(1);
        }
        abErr    = fabs(h_mut[i] - org_mut[i]);
        rtErr    = abErr / Max(fabs(org_mut[i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in GPUGetFaceValue in %d term, h_mut = %.30e, d_mut = %.30e "
                   "abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_mut[i], h_mut[i], abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            //!if ( fabs(h_mut[i] - org_mut[i]) > flagCompare ) {
            printf("Something wrong in GPUGetFaceValue in %d term, h_mut = %.30e, d_mut = %.30e "
                   "abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_mut[i], h_mut[i], abErr, rtErr);
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUGetFaceValue in computing mul and mut  successfully" << endl;
#endif

    delete[] h_mul;
    delete[] h_mut;
}

void TestGPUGetVisFaceValuePrimTmKcp(const int nst, const int ned, RFloat **org_prim, RFloat **org_tmid,
                                     const RFloat *org_kcp)
{
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat *h_prim = new RFloat[d_neqn_ns * d_SEG_LEN];
    RFloat *h_tm   = new RFloat[1 * d_SEG_LEN];
    RFloat *h_kcp  = new RFloat[d_SEG_LEN];

    size_t sizePrim = d_neqn_ns * d_SEG_LEN * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(h_prim, d_prim, sizePrim, cudaMemcpyDeviceToHost));
    size_t sizeTm = 1 * d_SEG_LEN * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(h_tm, d_tm, sizeTm, cudaMemcpyDeviceToHost));
    size_t sizeKcp = d_SEG_LEN * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(h_kcp, d_kcp, sizeKcp, cudaMemcpyDeviceToHost));
    for (int i = 0; i < d_SEG_LEN; i++)
    {
        for (int j = 0; j < d_neqn_ns; j++)
        {
            abErr           = fabs(org_prim[j][i] - h_prim[j * d_SEG_LEN + i]);
            rtErr           = abErr / Max(fabs(org_prim[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("something wrong in GPUGetVisFaceValuePrim, %d equation %d term, h_prim = "
                       "%.30e, d_prim = %.30e abErr= %.30e, rtErr= %.30e\n",
                       j, i, org_prim[j][i], h_prim[j * d_SEG_LEN + i], abErr, rtErr);
                //!printf("error = %.30le\n", fabs(org_prim[j][i] - h_prim[j * d_SEG_LEN + i]));
                exit(1);
            }
            //!if (fabs(org_prim[j][i] - h_prim[j * d_SEG_LEN + i]) > flagComparePrimTm) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("something wrong in GPUGetVisFaceValuePrim, %d equation %d term, h_prim = "
                       "%.30e, d_prim = %.30e abErr= %.30e, rtErr= %.30e\n",
                       j, i, org_prim[j][i], h_prim[j * d_SEG_LEN + i], abErr, rtErr);
                //!printf("error = %.30le\n", fabs(org_prim[j][i] - h_prim[j * d_SEG_LEN + i]));
                exit(1);
            }
        }
    }

    for (int i = 0; i < d_SEG_LEN; i++)
    {
        abErr           = fabs(org_tmid[0][i] - h_tm[i]);
        rtErr           = abErr / Max(fabs(org_tmid[0][i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            //!if (fabs(org_tmid[0][i] - h_tm[i]) > flagComparePrimTm) {
            printf("something wrong in GPUGetVisFaceValueTm, %d term, h_tm = %.30e, d_tm = %.30e "
                   "abErr= %.30e, rtErr= %.30e\n",
                   nst + i, org_tmid[0][i], h_tm[i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_tmid[0][i] - h_tm[i]));
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            //!if (fabs(org_tmid[0][i] - h_tm[i]) > flagComparePrimTm) {
            printf("something wrong in GPUGetVisFaceValueTm, %d term, h_tm = %.30e, d_tm = %.30e "
                   "abErr= %.30e, rtErr= %.30e\n",
                   nst + i, org_tmid[0][i], h_tm[i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_tmid[0][i] - h_tm[i]));
            exit(1);
        }
    }

    for (int i = 0; i < d_SEG_LEN; i++)
    {
        abErr           = fabs(org_kcp[i] - h_kcp[i]);
        rtErr           = abErr / Max(fabs(org_kcp[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("something wrong in GPUGetVisFaceValueKcp, %d term, h_kcp = %.30e, d_kcp = "
                   "%.30e abErr= %.30e, rtErr= %.30e\n",
                   i, org_kcp[i], h_kcp[i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_kcp[i] - h_kcp[i]));
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            //!if (fabs(org_kcp[i] - h_kcp[i])>flagCompareKcp){
            printf("something wrong in GPUGetVisFaceValueKcp, %d term, h_kcp = %.30e, d_kcp = "
                   "%.30e abErr= %.30e, rtErr= %.30e\n",
                   i, org_kcp[i], h_kcp[i], abErr, rtErr);
            //!printf("error = %.30le\n", fabs(org_kcp[i] - h_kcp[i]));
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUGetVisFaceValuePrimTmKcp successfully" << endl;
#endif
    delete[] h_prim;
    delete[] h_tm;
    delete[] h_kcp;
}

void TESTGPUCompVisfluxTEST(const int nst, RFloat **org_flux, const int nl)
{
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat *h_flux   = new RFloat[nl * d_SEG_LEN];
    size_t  sizeFlux = nl * d_SEG_LEN * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(h_flux, d_flux, sizeFlux, cudaMemcpyDeviceToHost));
    for (int i = 0; i < d_SEG_LEN; i++)
    {
        for (int j = 0; j < nl; j++)
        {
            abErr           = fabs(h_flux[j * d_SEG_LEN + i] - org_flux[j][i]);
            rtErr           = abErr / Max(fabs(org_flux[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong exists in GPUCompVisfluxTEST in term nst = %d, j = %d, i = "
                       "%d,  h_flux = %e, org_flux = %e abErr= %e, rtErr= %e\n",
                       nst, j, i, h_flux[j * d_SEG_LEN + i], org_flux[j][i], abErr, rtErr);
                //!printf("error =%.30le\n", fabs(h_flux[j * d_SEG_LEN + i] - org_flux[j][i]));
                exit(1);
            }
            //!if (fabs(h_flux[j * d_SEG_LEN + i] - org_flux[j][i]) > flagCompare){
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Something wrong exists in GPUCompVisfluxTEST in term nst = %d, j = %d, i = "
                       "%d,  h_flux = %e, org_flux = %e abErr= %e, rtErr= %e\n",
                       nst, j, i, h_flux[j * d_SEG_LEN + i], org_flux[j][i], abErr, rtErr);
                //!printf("error =%.30le\n", fabs(h_flux[j * d_SEG_LEN + i] - org_flux[j][i]));
                exit(1);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUCompVisfluxTEST successfully" << endl;
#endif
    delete[] h_flux;
}

void TestGPULoadFlux(const int nst, const int ned, const int nl, const int nTotalCell, const int nBoundFace,
                     RFloat **org_res)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUControlVariables;
    //!if (d_visflux_solver != "NSSolverUnstruct") return;
    //!RFloat flagCompare = 1.0e-8;
    double  abErr   = 0.0;
    double  rtErr   = 0.0;
    int     nTotal  = nTotalCell + nBoundFace;
    RFloat *h_res   = new RFloat[nl * nTotal];
    size_t  sizeRes = nl * nTotal * sizeof(RFloat);
    cudaMemcpy(h_res, d_res_ns, sizeRes, cudaMemcpyDeviceToHost);
    for (int i = 0; i < nTotal; i++)
    {
        for (int j = 0; j < nl; j++)
        {
            abErr           = abs(org_res[0][j * nTotal + i] - h_res[j * nTotal + i]);
            rtErr           = abErr / Max(fabs(org_res[0][j * nTotal + i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error in GPULoadFlux at %d equation, i=%d,  h_res=%.30e, d_res_ns=%.30e "
                       "abErr= %.30e, rtErr= %.30e\n",
                       j, i, org_res[j][i], h_res[j * nTotal + i], abErr, rtErr);
                //!printf("absolute error =%.30le \t relative error = %.30le\n", abErr, reErr);
                exit(1);
            }
            //!if ( abs(org_res[j][i] - h_res[j * nTotal + i]) > flagCompare ) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error in GPULoadFlux at %d equation, i=%d,  h_res=%.30e, d_res_ns=%.30e "
                       "abErr= %.30e, rtErr= %.30e\n",
                       j, i, org_res[j][i], h_res[j * nTotal + i], abErr, rtErr);
                //!printf("absolute error =%.30le \t relative error = %.30le\n", abErr, reErr);
                exit(1);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPULoadFlux successfully" << endl;
#endif
    delete[] h_res;
}

void TestGPUResCopy(RFloat **res, const int nl, const int nTotal)
{
    double  abErr   = 0.0;
    double  rtErr   = 0.0;
    RFloat *h_res   = new RFloat[nl * nTotal];
    size_t  sizeRes = nl * nTotal * sizeof(RFloat);
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    cudaMemcpy(h_res, d_res_ns, sizeRes, cudaMemcpyDeviceToHost);
    RFloat *res_proxy = res[0];
    for (int i = 0; i < nTotal * nl; i++)
    {
        abErr           = fabs(res_proxy[i] - h_res[i]);
        rtErr           = abErr / Max(fabs(res_proxy[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in GPUResCopy: i=%d, h_res=%e, d_res_ns=%e abErr= %e, rtErr= %e\n", i, res_proxy[i],
                   h_res[i], abErr, rtErr);
            //!printf("error =%.30le\n", abs(res_proxy[i] - h_res[i]));
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            //!if ( abs(res_proxy[i] - h_res[i]) > flagCompare ) {
            printf("Something wrong in GPUResCopy: i=%d, h_res=%e, d_res_ns=%e abErr= %e, rtErr= %e\n", i, res_proxy[i],
                   h_res[i], abErr, rtErr);
            //!printf("error =%.30le\n", abs(res_proxy[i] - h_res[i]));
            exit(1);
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUResCopy successfully" << endl;
#endif
    delete[] h_res;
}

void TestGPUGetTurbVisFaceValueMulMut(const int nst, const int ned, const RFloat *org_mul, const RFloat *org_mut)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat *h_mul = new RFloat[d_SEG_LEN];
    RFloat *h_mut = new RFloat[d_SEG_LEN];

    size_t sizeMult = d_SEG_LEN * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(h_mul, d_mul_turb, sizeMult, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_mut, d_mut_turb, sizeMult, cudaMemcpyDeviceToHost));

    for (int i = 0; i < d_SEG_LEN; i++)
    {
        abErr           = fabs(h_mul[i] - org_mul[i]);
        rtErr           = abErr / Max(fabs(org_mul[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in GPUGetTurbVisFaceValueMulMut in %d term, h_mul = %.30e, "
                   "d_mul = %.30e abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_mul[i], h_mul[i], abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in GPUGetTurbVisFaceValueMulMut in %d term, h_mul = %.30e, "
                   "d_mul = %.30e abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_mul[i], h_mul[i], abErr, rtErr);
            exit(1);
        }
        abErr    = fabs(h_mut[i] - org_mut[i]);
        rtErr    = abErr / Max(fabs(org_mut[i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in GPUGetTurbVisFaceValueMulMut in %d term, h_mut = %.30e, "
                   "d_mut = %.30e abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_mut[i], h_mut[i], abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in GPUGetTurbVisFaceValueMulMut in %d term, h_mut = %.30e, "
                   "d_mut = %.30e abErr= %.30e, rtErr= %.30e\n",
                   i + nst, org_mut[i], h_mut[i], abErr, rtErr);
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUGetTurbVisFaceValue in computing mul and mut  successfully" << endl;
#endif
    delete[] h_mul;
    delete[] h_mut;
}

void TestGPUGetTurbVisFaceValuePrim(const int nst, const int ned, const int n_turb, RFloat **org_prim)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    size_t  sizePrim = n_turb * d_SEG_LEN * sizeof(RFloat);
    int     j;
    RFloat *h_prim = new RFloat[n_turb * d_SEG_LEN];
    cudaMemcpy(h_prim, d_prim_turb, sizePrim, cudaMemcpyDeviceToHost);
    for (int iface = nst; iface < ned; iface++)
    {
        j = iface - nst;
        for (int m = 0; m < n_turb; m++)
        {
            abErr           = fabs(h_prim[m * d_SEG_LEN + j] - org_prim[m][j]);
            rtErr           = abErr / Max(fabs(org_prim[m][j]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: Something wrong in GPUGetTurbVisFaceValuePrim in %d term, "
                       "h_prim=%.30e, d_prim=%.30e abErr= %.30e, rtErr= %.30e\n",
                       m * d_SEG_LEN, org_prim[m][j], h_prim[m * d_SEG_LEN + j], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                //!if (fabs(h_prim[m * d_SEG_LEN + j]-org_prim[m][j]) > flagCompare){
                printf("Error: Something wrong in GPUGetTurbVisFaceValuePrim in %d term, "
                       "h_prim=%.30e, d_prim=%.30e abErr= %.30e, rtErr= %.30e\n",
                       m * d_SEG_LEN, org_prim[m][j], h_prim[m * d_SEG_LEN + j], abErr, rtErr);
                exit(1);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUGetTurbVisFaceValuePrim successfully" << endl;
#endif

    delete[] h_prim;
}

void TestGPUGetTurbVisFaceValueMlt(const int nst, const int ned, const int n_turb, RFloat **org_mlt)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    size_t  sizeMlt = n_turb * d_SEG_LEN * sizeof(RFloat);
    int     j;
    RFloat *h_mlt = new RFloat[n_turb * d_SEG_LEN];
    cudaMemcpy(h_mlt, d_mlt_turb, sizeMlt, cudaMemcpyDeviceToHost);
    for (int iface = nst; iface < ned; iface++)
    {
        j = iface - nst;
        for (int m = 0; m < n_turb; m++)
        {
            abErr           = fabs(h_mlt[m * d_SEG_LEN + j] - org_mlt[m][j]);
            rtErr           = abErr / Max(fabs(org_mlt[m][j]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: Something wrong in GPUGetTurbVisFaceValueMlt in %d term, "
                       "h_mlt=%.30e, d_mlt=%.30e abErr= %.30e, rtErr= %.30e\n",
                       m * d_SEG_LEN, org_mlt[m][j], h_mlt[m * d_SEG_LEN + j], abErr, rtErr);
                exit(1);
            }
            //!if (fabs(h_mlt[m * d_SEG_LEN + j]-org_mlt[m][j]) > flagCompare){
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: Something wrong in GPUGetTurbVisFaceValueMlt in %d term, "
                       "h_mlt=%.30e, d_mlt=%.30e abErr= %.30e, rtErr= %.30e\n",
                       m * d_SEG_LEN, org_mlt[m][j], h_mlt[m * d_SEG_LEN + j], abErr, rtErr);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUGetTurbVisFaceValueMlt successfully" << endl;
#endif
    delete[] h_mlt;
}

void TestGPUParametersMemAllocCopy(const int neqn)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    size_t  sizeTurboo = neqn * sizeof(double);
    double *h_turboo   = new double[neqn];
    cudaMemcpy(h_turboo, d_turboo, sizeTurboo, cudaMemcpyDeviceToHost);
    double *org_turboo = reinterpret_cast<double *>(GlobalDataBase::GetDataPtr("turboo"));
    for (int i = 0; i < neqn; i++)
    {
        abErr           = fabs(h_turboo[i] - org_turboo[i]);
        rtErr           = abErr / Max(fabs(org_turboo[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in copy turboo in %d term, h_turboo=%.30le, d_turboo=%.30le\n", i, org_turboo[i],
                   d_turboo[i]);
            printf("error is %.30le\n", fabs(h_turboo[i] - org_turboo[i]));
            exit(1);
        }
        //!if (fabs(h_turboo[i] - org_turboo[i]) > flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in copy turboo in %d term, h_turboo=%.30le, d_turboo=%.30le\n", i, org_turboo[i],
                   d_turboo[i]);
            printf("error is %.30le\n", fabs(h_turboo[i] - org_turboo[i]));
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test TestGPUParametersMemAllocCopy successfully" << endl;
#endif
    delete[] h_turboo;
}

void TestCallGPUCompTurbVisfluxTEST(const int nst, const int ned, const int n_turb, RFloat **org_flux_turb,
                                    RFloat **org_flux_sub_turb)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    size_t  sizeFlux        = d_SEG_LEN * n_turb * sizeof(RFloat);
    RFloat *h_flux_turb     = new RFloat[d_SEG_LEN * n_turb];
    RFloat *h_flux_sub_turb = new RFloat[d_SEG_LEN * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_flux_turb, d_flux_turb, sizeFlux, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_flux_sub_turb, d_flux_sub_turb, sizeFlux, cudaMemcpyDeviceToHost));
    for (int iface = nst; iface < ned; iface++)
    {
        int j = iface - nst;
        for (int m = 0; m < n_turb; m++)
        {
            abErr           = fabs(h_flux_turb[m * d_SEG_LEN + j] - org_flux_turb[m][j]);
            rtErr           = abErr / Max(fabs(org_flux_turb[m][j]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUCompTurbVisfluxTEST in %d term, h_flux_turb =  %e, "
                       "d_flux_turb = %e abErr= %e, rtErr= %e\n",
                       m * d_SEG_LEN + j, org_flux_turb[m][j], h_flux_turb[m * d_SEG_LEN + j], abErr, rtErr);
                //!printf("Error = %.30le\n", fabs(h_flux_turb[m * d_SEG_LEN + j] - org_flux_turb[m][j]));
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                //!if (fabs(h_flux_turb[m * d_SEG_LEN + j] - org_flux_turb[m][j]) > flagCompare){
                printf("Something wrong with GPUCompTurbVisfluxTEST in %d term, h_flux_turb =  %e, "
                       "d_flux_turb = %e abErr= %e, rtErr= %e\n",
                       m * d_SEG_LEN + j, org_flux_turb[m][j], h_flux_turb[m * d_SEG_LEN + j], abErr, rtErr);
                //!printf("Error = %.30le\n", fabs(h_flux_turb[m * d_SEG_LEN + j] - org_flux_turb[m][j]));
                exit(1);
            }
            abErr    = fabs(h_flux_sub_turb[m * d_SEG_LEN + j] - org_flux_sub_turb[m][j]);
            rtErr    = abErr / Max(fabs(org_flux_sub_turb[m][j]), smallLimit);
            checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUCompTurbVisfluxTEST in %d term, h_flux_sub_turb =  "
                       "%e, d_flux_sub_turb = %e abErr= %e, rtErr= %e\n",
                       m * d_SEG_LEN + j, org_flux_sub_turb[m][j], h_flux_sub_turb[m * d_SEG_LEN + j], abErr, rtErr);
                //!printf("Error = %.30le\n", fabs(h_flux_sub_turb[m * d_SEG_LEN + iface] - org_flux_sub_turb[m][iface]));
                exit(1);
            }
            //!if (fabs(h_flux_sub_turb[m * d_SEG_LEN + j] - org_flux_sub_turb[m][j]) > flagCompare){
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Something wrong with GPUCompTurbVisfluxTEST in %d term, h_flux_sub_turb =  "
                       "%e, d_flux_sub_turb = %e abErr= %e, rtErr= %e\n",
                       m * d_SEG_LEN + j, org_flux_sub_turb[m][j], h_flux_sub_turb[m * d_SEG_LEN + j], abErr, rtErr);
                //!printf("Error = %.30le\n", fabs(h_flux_sub_turb[m * d_SEG_LEN + iface] - org_flux_sub_turb[m][iface]));
                exit(1);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUCompTurbVisfluxTEST successfully" << endl;
#endif
    delete[] h_flux_turb;
    delete[] h_flux_sub_turb;
}

void TestCompInviscid(const int nst, const int ned, const int n_turb, RFloat **org_flux_turb,
                      RFloat **org_flux_sub_turb, const int flag)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    size_t  sizeFlux        = d_SEG_LEN * n_turb * sizeof(RFloat);
    RFloat *h_flux_turb     = new RFloat[d_SEG_LEN * n_turb];
    RFloat *h_flux_sub_turb = new RFloat[d_SEG_LEN * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_flux_turb, d_flux_turb, sizeFlux, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_flux_sub_turb, d_flux_sub_turb, sizeFlux, cudaMemcpyDeviceToHost));
    for (int iface = nst; iface < ned; iface++)
    {
        int j = iface - nst;
        if (flag == 1)
        {
            for (int m = 0; m < n_turb; m++)
            {
                abErr           = fabs(h_flux_turb[m * d_SEG_LEN + j] - org_flux_turb[m][j]);
                rtErr           = abErr / Max(fabs(org_flux_turb[m][j]), smallLimit);
                double checknan = abErr;
                if (std::isnan(checknan))
                {
                    printf("Something wrong with TestCompInviscid in %d term, h_flux_turb =  "
                           "%.30e, d_flux_turb = %.30e abErr= %.30e, rtErr= %.30e\n",
                           m * d_SEG_LEN + j, org_flux_turb[m][j], h_flux_turb[m * d_SEG_LEN + j], abErr, rtErr);
                    //!printf("Error = %.30le\n", fabs(h_flux_turb[m * d_SEG_LEN + j] - org_flux_turb[m][j]));
                    exit(1);
                }
                if (abErr > abCompare && rtErr > rtCompare)
                {
                    //!if (fabs(h_flux_turb[m * d_SEG_LEN + j] - org_flux_turb[m][j]) > flagCompare){
                    printf("Something wrong with TestCompInviscid in %d term, h_flux_turb =  "
                           "%.30e, d_flux_turb = %.30e abErr= %.30e, rtErr= %.30e\n",
                           m * d_SEG_LEN + j, org_flux_turb[m][j], h_flux_turb[m * d_SEG_LEN + j], abErr, rtErr);
                    //!printf("Error = %.30le\n", fabs(h_flux_turb[m * d_SEG_LEN + j] - org_flux_turb[m][j]));
                    exit(1);
                }
                abErr    = fabs(h_flux_sub_turb[m * d_SEG_LEN + j] - org_flux_sub_turb[m][j]);
                rtErr    = abErr / Max(fabs(org_flux_sub_turb[m][j]), smallLimit);
                checknan = abErr;
                if (std::isnan(checknan))
                {
                    printf("Something wrong with TestCompInviscid in %d term, h_flux_sub_turb =  "
                           "%.30e, d_flux_sub_turb = %.30e abErr= %.30e, rtErr= %.30e\n",
                           m * d_SEG_LEN + j, org_flux_sub_turb[m][j], h_flux_sub_turb[m * d_SEG_LEN + j], abErr,
                           rtErr);
                    //!printf("Error = %.30le\n", fabs(h_flux_sub_turb[m * d_SEG_LEN + iface] - org_flux_sub_turb[m][iface]));
                    exit(1);
                }
                //!if (fabs(h_flux_sub_turb[m * d_SEG_LEN + j] - org_flux_sub_turb[m][j]) > flagCompare){
                if (abErr > abCompare && rtErr > rtCompare)
                {
                    printf("Something wrong with TestCompInviscid in %d term, h_flux_sub_turb =  "
                           "%.30e, d_flux_sub_turb = %.30e abErr= %.30e, rtErr= %.30e\n",
                           m * d_SEG_LEN + j, org_flux_sub_turb[m][j], h_flux_sub_turb[m * d_SEG_LEN + j], abErr,
                           rtErr);
                    //!printf("Error = %.30le\n", fabs(h_flux_sub_turb[m * d_SEG_LEN + iface] - org_flux_sub_turb[m][iface]));
                    exit(1);
                }
                double err0 = fabs(h_flux_sub_turb[m * d_SEG_LEN + j] - 0.0);
                if (err0 > abCompare)
                {
                    printf("d_flux_sub is not equal to zero\n");
                    printf("d_flux_sub  = %.30e, h_flux_sub = %.30e\n", h_flux_sub_turb[m * d_SEG_LEN + j],
                           org_flux_sub_turb[m][j]);
                    exit(1);
                }
            }
        } //!end if (flag == 1)
    }

#ifdef UNITTESTOUTPUT
    cout << "Test CompInviscid successfully" << endl;
#endif
    delete[] h_flux_turb;
    delete[] h_flux_sub_turb;
}

void TestGPUFlux_turb()
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    int     n_turb          = 1;
    size_t  sizeFlux        = d_SEG_LEN * n_turb * sizeof(RFloat);
    RFloat *h_flux_sub_turb = new RFloat[d_SEG_LEN * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_flux_sub_turb, d_flux_sub_turb, sizeFlux, cudaMemcpyDeviceToHost));
    //!Output of flux_turb
    //!printf("d_flux_sub_turb[%d] = %e\n", 3, h_flux_sub_turb[3]);
    //!printf("d_flux_sub_turb[%d] = %e\n", 2325, h_flux_sub_turb[2325]);
    delete[] h_flux_sub_turb;
}

void TestGPUQ_turbCopy(RFloat **org_q_turb, const int nTotal, const int n_turb)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    size_t  sizeQturb = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_q_turb  = new RFloat[nTotal * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_q_turb, d_q_turb_proxy, sizeQturb, cudaMemcpyDeviceToHost));
    RFloat *org_q_turb0 = org_q_turb[0];
    for (int i = 0; i < nTotal * n_turb; i++)
    {
        abErr           = fabs(org_q_turb0[i] - h_q_turb[i]);
        rtErr           = abErr / (fabs(org_q_turb0[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong with %d term, h_q_turb = %.30le, d_q_turb = %.30le\n", i, org_q_turb0[i],
                   h_q_turb[i]);
            printf("Error = %.30le\n", fabs(org_q_turb0[i] - h_q_turb[i]));
            exit(1);
        }
        //!if (fabs(org_q_turb0[i] - h_q_turb[i]) > flagCompare){
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong with %d term, h_q_turb = %.30le, d_q_turb = %.30le\n", i, org_q_turb0[i],
                   h_q_turb[i]);
            printf("Error = %.30le\n", fabs(org_q_turb0[i] - h_q_turb[i]));
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUQ_turbCopy successfully" << endl;
#endif
    delete[] h_q_turb;
}

void TestGPULoadTurbFlux(RFloat **org_res_turb, const int nTotalCell, const int nBoundFace, const int n_turb)
{
    double abErr = 0.0;
    double rtErr = 0.0;

    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUControlVariables;
    //!if (d_visflux_solver != "TurbSolverUnstr") return;
    int     nTotal      = nTotalCell + nBoundFace;
    size_t  sizeResTurb = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_res_turb  = new RFloat[nTotal * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_res_turb, d_res_turb, sizeResTurb, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotal; i++)
    {
        for (int m = 0; m < n_turb; m++)
        {
            abErr           = fabs(h_res_turb[m * nTotal + i] - org_res_turb[m][i]);
            rtErr           = abErr / Max(fabs(org_res_turb[m][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("something wrong with GPULoadTurbFlux in %d term, h_res_turb = %e, "
                       "d_res_turb = %e  abErr= %e, rtErr= %e\n",
                       m * nTotal + i, h_res_turb[m * nTotal + i], org_res_turb[m][i], abErr, rtErr);
                //!printf("Error = %.30le\n", fabs(h_res_turb[m * nTotal + i]-org_res_turb[m][i]));
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                //!if (fabs(h_res_turb[m * nTotal + i]-org_res_turb[m][i]) > flagCompare) {
                printf("something wrong with GPULoadTurbFlux in %d term, h_res_turb = %e, "
                       "d_res_turb = %e  abErr= %e, rtErr= %e\n",
                       m * nTotal + i, h_res_turb[m * nTotal + i], org_res_turb[m][i], abErr, rtErr);
                //!printf("Error = %.30le\n", fabs(h_res_turb[m * nTotal + i]-org_res_turb[m][i]));
                exit(1);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPULoadTurbFlux successfully" << endl;
#endif
    delete[] h_res_turb;
}

void TestGPUResTurbCopy(RFloat **org_res_turb, const int nTotalCell, const int nBoundFace, const int n_turb)
{
    double abErr = 0.0;
    double rtErr = 0.0;

    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUControlVariables;
    //!if (d_visflux_solver != "TurbSolverUnstr") return;
    int     nTotal      = nTotalCell + nBoundFace;
    size_t  sizeResTurb = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_res_turb  = new RFloat[nTotal * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_res_turb, d_res_turb, sizeResTurb, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotal; i++)
    {
        for (int m = 0; m < n_turb; m++)
        {
            abErr           = fabs(h_res_turb[m * nTotal + i] - org_res_turb[m][i]);
            rtErr           = abErr / Max(fabs(org_res_turb[m][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("something wrong with GPUResTurbCopy in %d term, h_res_turb = %e, "
                       "d_res_turb = %e  abErr= %e, rtErr= %e\n",
                       m * nTotal + i, h_res_turb[m * nTotal + i], org_res_turb[m][i], abErr, rtErr);
                //!printf("Error = %.30le\n", fabs(h_res_turb[m * nTotal + i]-org_res_turb[m][i]));
                exit(1);
            }
            //!if (fabs(h_res_turb[m * nTotal + i]-org_res_turb[m][i]) > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("something wrong with GPUResTurbCopy in %d term, h_res_turb = %e, "
                       "d_res_turb = %e  abErr= %e, rtErr= %e\n",
                       m * nTotal + i, h_res_turb[m * nTotal + i], org_res_turb[m][i], abErr, rtErr);
                //!printf("Error = %.30le\n", fabs(h_res_turb[m * nTotal + i]-org_res_turb[m][i]));
                exit(1);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUResTurbCopy successfully" << endl;
#endif
    delete[] h_res_turb;
}

void TestGPUSpecTurbCell(const int n_turb, const int nTotalCell, const int nBoundFace, RFloat **org_spec_turb)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    int     nTotal       = nTotalCell + nBoundFace;
    size_t  sizeSpecTurb = n_turb * nTotal * sizeof(RFloat);
    RFloat *h_spec_turb  = new RFloat[n_turb * nTotal];
    HANDLE_API_ERR(cudaMemcpy(h_spec_turb, d_spec_turb, sizeSpecTurb, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotalCell; i++)
    {
        abErr           = fabs(org_spec_turb[0][i] - h_spec_turb[0 * nTotal + i]);
        rtErr           = abErr / Max(fabs(org_spec_turb[0][i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong with GPUSpecTurbCell in term %d, d_spec_turb = %.30e, "
                   "h_spec_turb = %.30e abErr= %.30e, rtErr= %.30e\n",
                   i, h_spec_turb[0 * nTotal + i], org_spec_turb[0][i], abErr, rtErr);
            //!printf("error = %.30el\n", abs(org_spec_turb[0][i] - h_spec_turb[0 * nTotal + i]));
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            //!if (fabs(org_spec_turb[0][i] - h_spec_turb[0 * nTotal + i]) > flagCompare) {
            printf("Something wrong with GPUSpecTurbCell in term %d, d_spec_turb = %.30e, "
                   "h_spec_turb = %.30e abErr= %.30e, rtErr= %.30e\n",
                   i, h_spec_turb[0 * nTotal + i], org_spec_turb[0][i], abErr, rtErr);
            //!printf("error = %.30el\n", abs(org_spec_turb[0][i] - h_spec_turb[0 * nTotal + i]));
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUSpecTurbCell successfully" << endl;
#endif
    delete[] h_spec_turb;
}

void TestGPUSpecTurbMatTurbFaces(const int n_turb, const int nTotalFace, const int nBoundFace, const int nTotalCell,
                                 RFloat **org_spec_turb, RFloat **org_mat_turbl, RFloat **org_mat_turbr)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    int     nTotal       = nBoundFace + nTotalCell;
    size_t  sizeSpecTurb = n_turb * nTotal * sizeof(RFloat);
    RFloat *h_spec_turb  = new RFloat[n_turb * nTotal];
    HANDLE_API_ERR(cudaMemcpy(h_spec_turb, d_spec_turb, sizeSpecTurb, cudaMemcpyDeviceToHost));
    size_t  sizeMatTurb = n_turb * nTotalFace * sizeof(RFloat);
    RFloat *h_mat_turbl = new RFloat[nTotalFace * n_turb];
    RFloat *h_mat_turbr = new RFloat[nTotalFace * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_mat_turbl, d_mat_turbl, sizeMatTurb, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_mat_turbr, d_mat_turbr, sizeMatTurb, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; i++)
    {
        for (int j = 0; j < n_turb; j++)
        {
            abErr           = fabs(h_spec_turb[j * nTotal + i] - org_spec_turb[j][i]);
            rtErr           = abErr / Max(fabs(org_spec_turb[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUSpecTurbMatTurbFaces in term %d, h_spec_turb = "
                       "%.30e, d_spec_turb = %.30e abErr= %.30e, rtErr= %.30e\n",
                       i, org_spec_turb[j][i], h_spec_turb[j * nTotal + i], abErr, rtErr);
                //!printf("error = %.30el\n", error);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                //!if (error > flagCompare) {
                printf("Something wrong with GPUSpecTurbMatTurbFaces in term %d, h_spec_turb = "
                       "%.30e, d_spec_turb = %.30e abErr= %.30e, rtErr= %.30e\n",
                       i, org_spec_turb[j][i], h_spec_turb[j * nTotal + i], abErr, rtErr);
                //!printf("error = %.30el\n", error);
                exit(1);
            }
        }
    }
    for (int i = 0; i < nTotalFace; i++)
    {
        for (int j = 0; j < n_turb; j++)
        {
            abErr           = fabs(h_mat_turbr[j * nTotalFace + i] - org_mat_turbr[j][i]);
            rtErr           = abErr / Max(fabs(org_mat_turbr[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUSpecTurbMatTurbFaces in term %d, h_mat_turbr = "
                       "%.30e, d_mat_turbr = %.30e\n",
                       i, org_mat_turbr[j][i], h_mat_turbr[j * nTotal + i], abErr, rtErr);
                //!printf("error = %.30el\n", error);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                //!if (error > flagCompare){
                printf("Something wrong with GPUSpecTurbMatTurbFaces in term %d, h_mat_turbr = "
                       "%.30e, d_mat_turbr = %.30e\n",
                       i, org_mat_turbr[j][i], h_mat_turbr[j * nTotal + i], abErr, rtErr);
                //!printf("error = %.30el\n", error);
                exit(1);
            }
            abErr    = fabs(h_mat_turbl[j * nTotalFace + i] - org_mat_turbl[j][i]);
            rtErr    = abErr / Max(fabs(org_mat_turbl[j][i]), smallLimit);
            checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUSpecTurbMatTurbFaces in term %d, h_mat_turbr = "
                       "%.30e, d_mat_turbr = %.30e abErr= %.30e, rtErr= %.30e\n",
                       i, org_mat_turbl[j][i], h_mat_turbl[j * nTotal + i], abErr, rtErr);
                //!printf("error = %.30el\n", error);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                //!if (error > flagCompare){
                printf("Something wrong with GPUSpecTurbMatTurbFaces in term %d, h_mat_turbr = "
                       "%.30e, d_mat_turbr = %.30e abErr= %.30e, rtErr= %.30e\n",
                       i, org_mat_turbl[j][i], h_mat_turbl[j * nTotal + i], abErr, rtErr);
                //!printf("error = %.30el\n", error);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUSpecTurbMatTurbFaces successfully" << endl;
#endif
    delete[] h_spec_turb;
    delete[] h_mat_turbl;
    delete[] h_mat_turbr;
}

void TestGPUSourceFlux_1eq_Original(RFloat **org_res_turb, RFloat **org_spec_turb, const int nTotalCell,
                                    const int nBoundFace, const int n_turb)
{
    double abErr = 0.0;
    double rtErr = 0.0;

    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    int     nTotal       = nBoundFace + nTotalCell;
    size_t  sizeSpecTurb = n_turb * nTotal * sizeof(RFloat);
    RFloat *h_spec_turb  = new RFloat[n_turb * nTotal];
    HANDLE_API_ERR(cudaMemcpy(h_spec_turb, d_spec_turb, sizeSpecTurb, cudaMemcpyDeviceToHost));
    size_t  sizeResTurb = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_res_turb  = new RFloat[nTotal * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_res_turb, d_res_turb, sizeResTurb, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotal; i++)
    {
        for (int j = 0; j < n_turb; j++)
        {
            abErr           = fabs(h_spec_turb[j * nTotal + i] - org_spec_turb[j][i]);
            rtErr           = abErr / Max(fabs(org_spec_turb[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUSourceFlux_1eq_Original in term %d, h_spec_turb = "
                       "%e, d_spec_turb = %e abErr= %e, rtErr= %e\n",
                       i, org_spec_turb[j][i], h_spec_turb[j * nTotal + i], abErr, rtErr);
                //!printf("error = %.30el\n", error);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                //!if (error > flagCompare) {
                printf("Something wrong with GPUSourceFlux_1eq_Original in term %d, h_spec_turb = "
                       "%e, d_spec_turb = %e abErr= %e, rtErr= %e\n",
                       i, org_spec_turb[j][i], h_spec_turb[j * nTotal + i], abErr, rtErr);
                //!printf("error = %.30el\n", error);
                exit(1);
            }
        }
    }
    for (int i = 0; i < nTotal; i++)
    {
        for (int m = 0; m < n_turb; m++)
        {
            abErr           = fabs(h_res_turb[m * nTotal + i] - org_res_turb[m][i]);
            rtErr           = abErr / Max(fabs(org_res_turb[m][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("something wrong with GPUSourceFlux_1eq_Original in %d term, h_res_turb = "
                       "%e, d_res_turb = %ebErr= %e, rtErr= %e\n",
                       m * nTotal + i, h_res_turb[m * nTotal + i], org_res_turb[m][i], abErr, rtErr);
                //!printf("Error = %.30le\n", fabs(h_res_turb[m * nTotal + i]-org_res_turb[m][i]));
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                //!if (fabs(h_res_turb[m * nTotal + i]-org_res_turb[m][i]) > flagCompare) {
                printf("something wrong with GPUSourceFlux_1eq_Original in %d term, h_res_turb = "
                       "%e, d_res_turb = %ebErr= %e, rtErr= %e\n",
                       m * nTotal + i, h_res_turb[m * nTotal + i], org_res_turb[m][i], abErr, rtErr);
                //!printf("Error = %.30le\n", fabs(h_res_turb[m * nTotal + i]-org_res_turb[m][i]));
                exit(1);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUSourceFlux_1eq_Original successfully" << endl;
#endif
    delete[] h_spec_turb;
    delete[] h_res_turb;
}
void TestGPUTurbBoundary(RFloat **org_q_turb, const int nTotal, const int n_turb)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    size_t  sizeQturb = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_q_turb  = new RFloat[nTotal * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_q_turb, d_q_turb_proxy, sizeQturb, cudaMemcpyDeviceToHost));
    RFloat *org_q_turb0 = org_q_turb[0];
    for (int i = 0; i < nTotal * n_turb; i++)
    {
        abErr           = fabs(org_q_turb0[i] - h_q_turb[i]);
        rtErr           = abErr / Max(fabs(org_q_turb0[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in GPUTurbBoundary with %d term, h_q_turb = %.30e, d_q_turb = "
                   "%.30e\n",
                   i, org_q_turb0[i], h_q_turb[i]);
            printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
            exit(1);
        }
        //!if (fabs(org_q_turb0[i] - h_q_turb[i]) > flagCompare){
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in GPUTurbBoundary with %d term, h_q_turb = %.30e, d_q_turb = "
                   "%.30e\n",
                   i, org_q_turb0[i], h_q_turb[i]);
            printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUTurbBoundary successfully" << endl;
#endif
    delete[] h_q_turb;
}
void TestGPUTurbLoadResiduals(RFloat **org_res, const int nTotal, const int n_turb)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    size_t  sizeRes = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_res   = new RFloat[nTotal * n_turb];
    RFloat  Error   = 0.0;
    HANDLE_API_ERR(cudaMemcpy(h_res, d_res_turb, sizeRes, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotal; i++)
    {
        for (int j = 0; j < n_turb; j++)
        {
            //!Error = fabs(h_res[j * nTotal +i] - org_res[j][i]);
            abErr           = fabs(h_res[j * nTotal + i] - org_res[j][i]);
            rtErr           = abErr / Max(fabs(org_res[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUTurbLoadResiduals in %d term, h_res = %.30el, "
                       "d_res = %.30el\n",
                       j * nTotal + i, org_res[j][i], h_res[j * nTotal + i]);
                printf("Error = %.30el\n", Error);
                exit(1);
            }
            //!if (Error > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Something wrong with GPUTurbLoadResiduals in %d term, h_res = %.30el, "
                       "d_res = %.30el\n",
                       j * nTotal + i, org_res[j][i], h_res[j * nTotal + i]);
                printf("Error = %.30el\n", Error);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUTurbLoadResiduals successfully" << endl;
#endif
    delete[] h_res;
}

void TestGPUTurbFillField(RFloat **org_field, const int nTotalCell, const int nTotal, const int neqn)
{
    using namespace GPUMemory;
    using namespace GPUControlVariables;
    using namespace GPUFlowVariables;
    double  abErr     = 0.0;
    double  rtErr     = 0.0;
    RFloat  error     = 0.0;
    size_t  sizeField = nTotal * neqn * sizeof(RFloat);
    RFloat *h_field   = new RFloat[nTotal * neqn];
    if (d_FillField_turb == "q_proxyToq_proxy_temp")
    {
        HANDLE_API_ERR(cudaMemcpy(h_field, d_q_proxy_temp_turb, sizeField, cudaMemcpyDeviceToHost));
    }
    else if (d_FillField_turb == "q_proxy_tempToq_proxy")
    {
        //!HANDLE_API_ERR(cudaMemcpy(h_field, d_q_proxy_turb, sizeField, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(h_field, d_q_turb_proxy, sizeField, cudaMemcpyDeviceToHost));
    }
    else if (d_FillField_turb == "resTodq")
    {
        HANDLE_API_ERR(cudaMemcpy(h_field, d_dq_turb, sizeField, cudaMemcpyDeviceToHost));
    }
    else
    {
        cout << "Error : d_FillField_turb owns a wrong value " << d_FillField_turb
             << ", which should be q_proxyToq_proxy_temp or q_proxy_tempToq_proxy" << endl;
        exit(1);
    }

    for (int i = 0; i < nTotalCell; i++)
    {
        for (int j = 0; j < neqn; j++)
        {
            //!error = fabs(h_field[j * nTotal + i]- org_field[j][i]);
            abErr           = fabs(h_field[j * nTotal + i] - org_field[j][i]);
            rtErr           = abErr / Max(fabs(org_field[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUTurbFillField in %d term, d_field = %.30e, h_field "
                       "= %.30e\n",
                       j * nTotal + i, h_field[j * nTotal + i], org_field[j][i]);
                cout << "d_FillField_turb = " << d_FillField_turb << endl;
                printf("abErr = %.30e, rtErr = %.30e\n", error);
                exit(1);
            }
            //!if (error > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Something wrong with GPUTurbFillField in %d term, d_field = %.30e, h_field "
                       "= %.30e\n",
                       j * nTotal + i, h_field[j * nTotal + i], org_field[j][i]);
                cout << "d_FillField_turb = " << d_FillField_turb << endl;
                printf("abErr = %.30e, rtErr = %.30e\n", error);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUTurbFillField successfully with d_FillField_turb =" << d_FillField_turb << endl;
#endif

    delete[] h_field;
}

void TestGPUTurbLoadQ(RFloat **org_qq_turb, const int nTotalCell, const int nTotal, const int n_turb)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat  error          = 0.0;
    size_t  sizeQProxyTurb = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_qq_turb      = new RFloat[nTotal * n_turb];
    //!HANDLE_API_ERR(cudaMemcpy(h_qq_turb, d_q_proxy_turb, sizeQProxyTurb, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_qq_turb, d_q_turb_proxy, sizeQProxyTurb, cudaMemcpyDeviceToHost));
    //!if i < nTotal, error exists due to ghost cell.
    for (int i = 0; i < nTotalCell; i++)
    {
        for (int j = 0; j < n_turb; j++)
        {
            //!error = fabs(h_qq_turb[j*nTotal + i] - org_qq_turb[j][i]);
            abErr           = fabs(h_qq_turb[j * nTotal + i] - org_qq_turb[j][i]);
            rtErr           = abErr / Max(fabs(org_qq_turb[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUTurbLoadQ in %d term, h_q_proxy_turb = %.30e, "
                       "d_q_turb_proxy = %.30e\n",
                       j * nTotal + i, org_qq_turb[j][i], h_qq_turb[j * nTotal + i]);
                printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
                exit(1);
            }
            //!if (error > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                //!printf("Something wrong with GPUTurbLoadQ in %d term, h_q_proxy_turb = %.30el, d_q_proxy_turb = %.30el\n", j*nTotal+i, org_qq_turb[j][i], h_qq_turb[j*nTotal + i]);
                printf("Something wrong with GPUTurbLoadQ in %d term, h_q_proxy_turb = %.30e, "
                       "d_q_turb_proxy = %.30e\n",
                       j * nTotal + i, org_qq_turb[j][i], h_qq_turb[j * nTotal + i]);
                printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUTurbLoadQ successfully" << endl;
#endif
    delete[] h_qq_turb;
}

void TestGPUTurbStoreRhsByResidual(RFloat **org_rhs, const int nTotal, const int n_turb)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    double  abErr   = 0.0;
    double  rtErr   = 0.0;
    RFloat  error   = 0.0;
    size_t  rhsSize = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_rhs   = new RFloat[nTotal * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_rhs, d_rhs_turb, rhsSize, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; i++)
    {
        for (int j = 0; j < n_turb; j++)
        {
            //!error = fabs(h_rhs[j * nTotal + i] - org_rhs[j][i]);
            abErr           = fabs(h_rhs[j * nTotal + i] - org_rhs[j][i]);
            rtErr           = abErr / Max(fabs(org_rhs[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUTurbStoreRhsByResidual in %d term, h_rhs = %.30el, "
                       "d_rhs = %.30el\n",
                       j * nTotal + i, org_rhs[j][i], h_rhs[j * nTotal + i]);
                printf("error = %.30el\n", error);
                exit(1);
            }
            //!if (error > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Something wrong with GPUTurbStoreRhsByResidual in %d term, h_rhs = %.30el, "
                       "d_rhs = %.30el\n",
                       j * nTotal + i, org_rhs[j][i], h_rhs[j * nTotal + i]);
                printf("error = %.30el\n", error);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUTurbStoreRhsByResidual successfully" << endl;
#endif
    delete[] h_rhs;
}

void TestGPUTurbZeroResiduals(RFloat **org_res, const int numberOfGeneralCells, const int numberOfEquations)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat  error   = 0.0;
    size_t  sizeRes = numberOfGeneralCells * numberOfEquations * sizeof(RFloat);
    RFloat *h_res   = new RFloat[numberOfGeneralCells * numberOfEquations];
    HANDLE_API_ERR(cudaMemcpy(h_res, d_res_turb, sizeRes, cudaMemcpyDeviceToHost));
    for (int i = 0; i < numberOfGeneralCells; i++)
    {
        for (int j = 0; j < numberOfEquations; j++)
        {
            //!error = fabs(h_res[j * numberOfGeneralCells + i] - org_res[j][i]);
            abErr           = fabs(h_res[j * numberOfGeneralCells + i] - org_res[j][i]);
            rtErr           = abErr / Max(fabs(org_res[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUTurbZeroResiduals in %d term, h_res = %.30el, "
                       "d_res = %.30el\n",
                       j * numberOfGeneralCells + i, h_res[j * numberOfGeneralCells + i], org_res[j][i]);
                printf("error = %.30el\n", error);
                exit(1);
            }
            //!if (error > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Something wrong with GPUTurbZeroResiduals in %d term, h_res = %.30el, "
                       "d_res = %.30el\n",
                       j * numberOfGeneralCells + i, h_res[j * numberOfGeneralCells + i], org_res[j][i]);
                printf("error = %.30el\n", error);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUTurbZeroResiduals successfully" << endl;
#endif
    delete[] h_res;
}

void TestGPUSpecTurbInit(RFloat **org_spec_turb, const int nTotal, const int n_turb)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat  error        = 0.0;
    size_t  sizeSpecTurb = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_specTurb   = new RFloat[nTotal * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_specTurb, d_spec_turb, sizeSpecTurb, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotal; i++)
    {
        for (int j = 0; j < n_turb; j++)
        {
            //!error = fabs(h_specTurb[j * nTotal + i] - org_spec_turb[j][i]);
            abErr           = fabs(h_specTurb[j * nTotal + i] - org_spec_turb[j][i]);
            rtErr           = abErr / Max(fabs(org_spec_turb[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUSpecTurbInit in %d term, h_specTurb = %.30el, "
                       "d_specTurb = %.30el\n",
                       j * nTotal + i, org_spec_turb[j][i], h_specTurb[j * nTotal + i]);
                printf("error = %.30el\n", error);
                exit(1);
            }
            //!if (error > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Something wrong with GPUSpecTurbInit in %d term, h_specTurb = %.30el, "
                       "d_specTurb = %.30el\n",
                       j * nTotal + i, org_spec_turb[j][i], h_specTurb[j * nTotal + i]);
                printf("error = %.30el\n", error);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUSpecTurbInit successfully" << endl;
#endif
    delete[] h_specTurb;
}

void TestGPUTurbLhs(const int nTotalCell, const int nTotal, const int n_turb, RFloat **org_dq_turb)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat  error      = 0.0;
    size_t  sizeDqTurb = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_dq_turb  = new RFloat[nTotal * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_dq_turb, d_res_turb, sizeDqTurb, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotalCell; i++)
    {
        for (int j = 0; j < n_turb; j++)
        {
            //!error = fabs(h_dq_turb[j * nTotal + i] - org_dq_turb[j][i]);
            abErr           = fabs(h_dq_turb[j * nTotal + i] - org_dq_turb[j][i]);
            rtErr           = abErr / Max(fabs(org_dq_turb[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUTurbLhs in %d term, h_dq_turb = %.30el, d_dq_turb "
                       "= %.30el",
                       j * nTotal + i, org_dq_turb[j][i], h_dq_turb[j * nTotal + i]);
                printf("error = %.30el\n", error);
                exit(1);
            }
            //!if (error > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Something wrong with GPUTurbLhs in %d term, h_dq_turb = %.30el, d_dq_turb "
                       "= %.30el",
                       j * nTotal + i, org_dq_turb[j][i], h_dq_turb[j * nTotal + i]);
                printf("error = %.30el\n", error);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUTurbLhs successfully" << endl;
#endif
    delete[] h_dq_turb;
}

void TestGPUTurbUpdateFlowField(const int nTotalCell, const int nTotal, const int n_turb, RFloat **q, const int n_neg,
                                const int n_pos)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat error = 0.0;

    size_t  sizeQ = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_q   = new RFloat[nTotal * n_turb];
    //!HANDLE_API_ERR(cudaMemcpy(h_q, d_q_proxy_turb, sizeQ, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_q, d_q_turb_proxy, sizeQ, cudaMemcpyDeviceToHost));
    int *h_n_neg = new int;
    int *h_n_pos = new int;

    HANDLE_API_ERR(cudaMemcpy(h_n_neg, d_n_neg, sizeof(int), cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_n_pos, d_n_pos, sizeof(int), cudaMemcpyDeviceToHost));
    if (fabs((n_neg) - (*h_n_neg)) != 0)
    {
        printf("Something wrong with GPUTurbUpdateFlowField in neg count, h_n_neg = %d, d_n_neg = "
               "%d\n",
               n_neg, *h_n_neg);
        exit(1);
    }
    if (fabs((n_pos) - (*h_n_pos)) != 0)
    {
        printf("Something wrong with GPUTurbUpdateFlowField in neg count, h_n_pos = %d, d_n_pos = "
               "%d\n",
               n_pos, *h_n_pos);
        exit(1);
    }
    for (int i = 0; i < nTotalCell; i++)
    {
        for (int j = 0; j < n_turb; j++)
        {
            //!error = fabs(h_q[j * nTotal + i] - q[j][i]);
            abErr           = fabs(h_q[j * nTotal + i] - q[j][i]);
            rtErr           = abErr / Max(fabs(q[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with GPUTurbUpdateFlowField in %d equation %d term, h_q = "
                       "%.30e, d_q = %.30e\n",
                       j, i, q[j][i], h_q[j * nTotal + i]);
                printf(" abErr=%.30el, rtErr=%.30e\n", abErr, rtErr);
                exit(1);
            }
            //!if (error > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Something wrong with GPUTurbUpdateFlowField in %d equation %d term, h_q = "
                       "%.30e, d_q = %.30e\n",
                       j, i, q[j][i], h_q[j * nTotal + i]);
                printf(" abErr=%.30el, rtErr=%.30e\n", abErr, rtErr);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUTurbUpdateFlowField successfully" << endl;
#endif
    delete[] h_q;
}

void TestCallGPURecoverResidual(const int nTotal, const int n_turb, RFloat **org_res)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat  error   = 0.0;
    size_t  sizeRes = nTotal * n_turb * sizeof(RFloat);
    RFloat *h_res   = new RFloat[nTotal * n_turb];
    HANDLE_API_ERR(cudaMemcpy(h_res, d_res_turb, sizeRes, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotal; i++)
    {
        for (int j = 0; j < n_turb; j++)
        {
            //!error = fabs(h_res[j * nTotal + i] - org_res[j][i]);
            abErr           = fabs(h_res[j * nTotal + i] - org_res[j][i]);
            rtErr           = abErr / Max(fabs(org_res[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error in GPURecoverResidual\n");
                exit(1);
            }
            //!if (error > flagCompare) {
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error in GPURecoverResidual\n");
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPURecoverResidual successfully" << endl;
#endif
    delete[] h_res;
}

void TestGPUMinMax(const int nTotal, const int number, const RFloat *org_dmin, const RFloat *org_dmax)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUControlVariables;
    RFloat  error      = 0.0;
    size_t  sizeMinMax = nTotal * sizeof(RFloat);
    RFloat *h_dmin     = new RFloat[nTotal];
    RFloat *h_dmax     = new RFloat[nTotal];
    HANDLE_API_ERR(cudaMemcpy(h_dmin, d_dmin, sizeMinMax, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_dmax, d_dmax, sizeMinMax, cudaMemcpyDeviceToHost));
    for (int i = 0; i < number; i++)
    {
        //!error = fabs(h_dmin[i] - org_dmin[i]);
        abErr           = fabs(h_dmin[i] - org_dmin[i]);
        rtErr           = abErr / Max(fabs(org_dmin[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in GPUMinMax in %d term, d_dmin = %.30e, h_dmin = %.30e\n", i, h_dmin[i],
                   org_dmin[i]);
            cout << "d_Limiter_var = " << d_Limiter_var << endl;
            printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
            exit(1);
        }
        //!if (error > flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in GPUMinMax in %d term, d_dmin = %.30e, h_dmin = %.30e\n", i, h_dmin[i],
                   org_dmin[i]);
            cout << "d_Limiter_var = " << d_Limiter_var << endl;
            printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
            exit(1);
        }
        //!error = fabs(h_dmax[i] - org_dmax[i]);
        abErr    = fabs(h_dmax[i] - org_dmax[i]);
        rtErr    = abErr / Max(fabs(org_dmax[i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in GPUMinMax in %d term, d_dmax = %.30e, h_dmax = %.30e\n", i, h_dmax[i],
                   org_dmax[i]);
            cout << "d_Limiter_var = " << d_Limiter_var << endl;
            printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
            exit(1);
        }
        //!if (error > flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in GPUMinMax in %d term, d_dmax = %.30e, h_dmax = %.30e\n", i, h_dmax[i],
                   org_dmax[i]);
            cout << "d_Limiter_var = " << d_Limiter_var << endl;
            printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUMinMax successfully" << endl;
#endif
    delete[] h_dmin;
    delete[] h_dmax;
}

void TestGPULimit(const int nTotal, const RFloat *org_limit)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat error     = 0.0;
    size_t sizeLimit = nTotal * sizeof(RFloat);

    RFloat *h_limit = new RFloat[nTotal];
    HANDLE_API_ERR(cudaMemcpy(h_limit, d_limit, sizeLimit, cudaMemcpyDeviceToHost));

    for (int i = 0; i < nTotal; i++)
    {
        //!error = h_limit[i] - org_limit[i];
        abErr           = fabs(h_limit[i] - org_limit[i]);
        rtErr           = abErr / Max(fabs(org_limit[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong with limit calculation in %d term, h_limit = %.30e, d_limit = "
                   "%.30e\n",
                   i, org_limit[i], h_limit[i]);
            printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
            exit(1);
        }
        //!if (error > flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong with limit calculation in %d term, h_limit = %.30e, d_limit = "
                   "%.30e\n",
                   i, org_limit[i], h_limit[i]);
            printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPULimit calculation successfully" << endl;
#endif
    delete[] h_limit;
}

void TestGPUVisSpectrum(const int nTotalCell, const RFloat *org_visSpectrum)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat error           = 0.0;
    size_t sizeVisSpectrum = nTotalCell * sizeof(RFloat);

    RFloat *h_visSpectrum = new RFloat[nTotalCell];
    HANDLE_API_ERR(cudaMemcpy(h_visSpectrum, d_visSpectrumRadius, sizeVisSpectrum, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotalCell; i++)
    {
        //!error = fabs(h_visSpectrum[i] - org_visSpectrum[i]);
        abErr           = fabs(h_visSpectrum[i] - org_visSpectrum[i]);
        rtErr           = abErr / Max(fabs(org_visSpectrum[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong with GPUVisSpectrum in %d term, h_visSpectrum = %.30e, "
                   "d_visSpectrum = %.30e\n",
                   i, org_visSpectrum[i], h_visSpectrum[i]);
            printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
            exit(1);
        }
        //!if (error > flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            if (processNum > 1)
                printf("In globalRank %d, localRank %d, Something wrong with GPUVisSpectrum in %d "
                       "term, h_visSpectrum = %.30e, d_visSpectrum = %.30e\n",
                       globalRank, localRank, i, org_visSpectrum[i], h_visSpectrum[i]);
            else
                printf("Something wrong with GPUVisSpectrum in %d term, h_visSpectrum = %.30e, "
                       "d_visSpectrum = %.30e\n",
                       i, org_visSpectrum[i], h_visSpectrum[i]);
            printf("abErr = %.30e, rtErr = %.30e\n", abErr, rtErr);
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUVisSpectrum successfully" << endl;
#endif
    delete[] h_visSpectrum;
}
void TestGPUInvSpectrum(const int nTotalCell, const RFloat *org_invSpectrum)
{
    double abErr = 0.0;
    double rtErr = 0.0;
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat error           = 0.0;
    size_t sizeInvSpectrum = nTotalCell * sizeof(RFloat);

    RFloat *h_invSpectrum = new RFloat[nTotalCell];
    HANDLE_API_ERR(cudaMemcpy(h_invSpectrum, d_invSpectrumRadius, sizeInvSpectrum, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotalCell; i++)
    {
        //!error = fabs(h_invSpectrum[i] - org_invSpectrum[i]);
        abErr           = fabs(h_invSpectrum[i] - org_invSpectrum[i]);
        rtErr           = abErr / Max(fabs(org_invSpectrum[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong with GPUInvSpectrum in %d term, h_invSpectrum = %.30e, "
                   "d_invSpectrum = %.30e, abErr = %.30e\n",
                   i, org_invSpectrum[i], h_invSpectrum[i], abErr);
            //!printf("error = %.30el\n", error);
            exit(1);
        }
        //!if (error > flagCompare) {
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong with GPUInvSpectrum in %d term, h_invSpectrum = %.30e, "
                   "d_invSpectrum = %.30e, abErr = %.30e\n",
                   i, org_invSpectrum[i], h_invSpectrum[i], abErr);
            //!printf("error = %.30el\n", error);
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUInvSpectrum successfully" << endl;
#endif
    delete[] h_invSpectrum;
}
void TestComputeNodeValue(const int nTotalNode, const int nBoundFace, const int nEquation, RFloat **org_qNode,
                          RFloat **org_tNode, const int *org_nodeBC, RFloat **org_qOnBCFace)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    size_t  sizeNodeRFloat  = nTotalNode * sizeof(RFloat);
    size_t  sizeQNodeRFloat = nEquation * nTotalNode * sizeof(RFloat);
    size_t  sizeNodeInt     = nTotalNode * sizeof(int);
    RFloat *h_qNode         = new RFloat[nEquation * nTotalNode];
    RFloat *h_tNode         = new RFloat[nTotalNode];
    int    *h_nodeBC        = new int[nTotalNode];
    RFloat *h_qOnBCFace     = new RFloat[(nEquation + 1) * nBoundFace];
    HANDLE_API_ERR(cudaMemcpy(h_qNode, d_qNode, sizeQNodeRFloat, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_tNode, d_tNode, sizeNodeRFloat, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_nodeBC, d_nodeBC, sizeNodeInt, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(
        cudaMemcpy(h_qOnBCFace, d_qOnBCFace, (nEquation + 1) * nBoundFace * sizeof(RFloat), cudaMemcpyDeviceToHost));

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    for (int i = 0; i < nTotalNode; i++)
    {
        //!check nCount
        abErr           = fabs(h_nodeBC[i] - org_nodeBC[i]);
        rtErr           = abErr / Max(fabs(org_nodeBC[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in ComputeNodeValue in %d term, h_nodeBC = %d, d_nodeBC=%d\n", i, org_nodeBC[i],
                   h_nodeBC[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in ComputeNodeValue in %d term, h_nodeBC = %d, d_nodeBC=%d\n", i, org_nodeBC[i],
                   h_nodeBC[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
        for (int j = 0; j < nEquation; j++)
        {
            //!check qNode
            abErr    = fabs(h_qNode[j * nTotalNode + i] - org_qNode[j][i]);
            rtErr    = abErr / Max(fabs(org_qNode[j][i]), smallLimit);
            checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong in ComputeNodeValue in %d equation, %d node, h_qNode = %e, "
                       "d_qNode=%e\n",
                       j, i, org_qNode[j][i], h_qNode[j * nTotalNode + i]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                if (processNum > 1)
                    printf("globalRank = %d, localRank = %d, globalZoneID = %d, Something wrong in "
                           "ComputeNodeValue in %d equation, %d node, h_qNode = %e, d_qNode=%e\n",
                           globalRank, localRank, globalZoneID, j, i, org_qNode[j][i], h_qNode[j * nTotalNode + i]);
                else
                    printf("Something wrong in ComputeNodeValue in %d equation, %d node, h_qNode = "
                           "%e, d_qNode=%e\n",
                           j, i, org_qNode[j][i], h_qNode[j * nTotalNode + i]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }

        //!check tNode
        abErr    = fabs(h_tNode[0 * nTotalNode + i] - org_tNode[0][i]);
        rtErr    = abErr / Max(fabs(org_tNode[0][i]), smallLimit);
        checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in ComputeNodeValue in %d term, h_tNode = %e, d_tNode=%e\n", 0 * nTotalNode + i,
                   org_tNode[0][i], h_qNode[0 * nTotalNode + i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            if (processNum > 1)
                printf("globalRank = %d, localRank = %d, globalZoneID = %d, Something wrong in "
                       "ComputeNodeValue in %d term, h_tNode = %e, d_tNode=%e\n",
                       globalRank, localRank, globalZoneID, 0 * nTotalNode + i, org_tNode[0][i],
                       h_qNode[0 * nTotalNode + i]);
            else
                printf("Something wrong in ComputeNodeValue in %d term, h_tNode = %e, d_tNode=%e\n", 0 * nTotalNode + i,
                       org_tNode[0][i], h_qNode[0 * nTotalNode + i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
    }

    //!check qOnBCFace
    for (int i = 0; i < nBoundFace; i++)
    {
        for (int j = 0; j < nEquation + 1; j++)
        {
            //!check qNode
            abErr           = fabs(h_qOnBCFace[j * nBoundFace + i] - org_qOnBCFace[j][i]);
            rtErr           = abErr / Max(fabs(org_qOnBCFace[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong in ComputeNodeValue in %d equation, %d node, h_qOnBCFace = "
                       "%e, d_qOnBCFace=%e\n",
                       j, i, org_qOnBCFace[j][i], h_qOnBCFace[j * nTotalNode + i]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                if (processNum > 1)
                    printf("globalRank = %d, localRank = %d, globalZoneID = %d, Something wrong in "
                           "ComputeNodeValue in %d equation, %d node, h_qOnBCFace = %e, "
                           "d_qOnBCFace=%e\n",
                           globalRank, localRank, globalZoneID, j, i, org_qOnBCFace[j][i],
                           h_qOnBCFace[j * nTotalNode + i]);
                else
                    printf("Something wrong in ComputeNodeValue in %d equation, %d node, "
                           "h_qOnBCFace = %e, d_qOnBCFace=%e\n",
                           j, i, org_qOnBCFace[j][i], h_qOnBCFace[j * nBoundFace + i]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test the kernel ComputeNodeValue successfully!" << endl;
#endif

    delete[] h_qNode;
    delete[] h_tNode;
    delete[] h_qOnBCFace;
    delete[] h_nodeBC;
}

void TestNodeValue(const int nTotalNode, const int nEquation, RFloat **org_qNode, RFloat **org_tNode)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    size_t  sizeNodeRFloat  = nTotalNode * sizeof(RFloat);
    size_t  sizeQNodeRFloat = nEquation * nTotalNode * sizeof(RFloat);
    size_t  sizeNodeInt     = nTotalNode * sizeof(int);
    RFloat *h_qNode         = new RFloat[nEquation * nTotalNode];
    RFloat *h_tNode         = new RFloat[nTotalNode];
    HANDLE_API_ERR(cudaMemcpy(h_qNode, d_qNode, sizeQNodeRFloat, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_tNode, d_tNode, sizeNodeRFloat, cudaMemcpyDeviceToHost));

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    for (int i = 0; i < nTotalNode; i++)
    {
        for (int j = 0; j < nEquation; j++)
        {
            //!check qNode
            abErr           = fabs(h_qNode[j * nTotalNode + i] - org_qNode[j][i]);
            rtErr           = abErr / Max(fabs(org_qNode[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong in TestNodeValue in %d equation, %d node, h_qNode = %e, "
                       "d_qNode=%e\n",
                       j, i, org_qNode[j][i], h_qNode[j * nTotalNode + i]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                if (processNum > 1)
                    printf("globalRank = %d, localRank = %d, globalZoneID = %d, Something wrong in "
                           "ComputeNodeValue in %d equation, %d node, h_qNode = %e, d_qNode=%e\n",
                           globalRank, localRank, globalZoneID, j, i, org_qNode[j][i], h_qNode[j * nTotalNode + i]);
                else
                    printf("Something wrong in TestNodeValue in %d equation, %d node, h_qNode = "
                           "%e, d_qNode=%e\n",
                           j, i, org_qNode[j][i], h_qNode[j * nTotalNode + i]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }

        //!check tNode
        abErr           = fabs(h_tNode[0 * nTotalNode + i] - org_tNode[0][i]);
        rtErr           = abErr / Max(fabs(org_tNode[0][i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in TestNodeValue in %d term, h_tNode = %e, d_tNode=%e\n", 0 * nTotalNode + i,
                   org_tNode[0][i], h_qNode[0 * nTotalNode + i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            if (processNum > 1)
                printf("globalRank = %d, localRank = %d, globalZoneID = %d, Something wrong in "
                       "ComputeNodeValue in %d term, h_tNode = %e, d_tNode=%e\n",
                       globalRank, localRank, globalZoneID, 0 * nTotalNode + i, org_tNode[0][i],
                       h_qNode[0 * nTotalNode + i]);
            else
                printf("Something wrong in TestNodeValue in %d term, h_tNode = %e, d_tNode=%e\n", 0 * nTotalNode + i,
                       org_tNode[0][i], h_qNode[0 * nTotalNode + i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "TestNodeValue is successfull!" << endl;
#endif

    delete[] h_qNode;
    delete[] h_tNode;
}
void TestNodeWeight(const int nTotalNode, RFloat *org_NodeWeight)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    size_t  sizeNodeRFloat = nTotalNode * sizeof(RFloat);
    RFloat *h_NodeWeight   = new RFloat[nTotalNode];
    HANDLE_API_ERR(cudaMemcpy(h_NodeWeight, d_nodeWeight, sizeNodeRFloat, cudaMemcpyDeviceToHost));

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    for (int i = 0; i < nTotalNode; i++)
    {
        //!check tNode
        abErr           = fabs(h_NodeWeight[i] - org_NodeWeight[i]);
        rtErr           = abErr / Max(fabs(org_NodeWeight[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in TestNodeWeight in %d term, h_NodeWeight = %e, "
                   "d_NodeWeight=%e\n",
                   0 * nTotalNode + i, org_NodeWeight[i], h_NodeWeight[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            if (processNum > 1)
                printf("globalRank = %d, localRank = %d, globalZoneID = %d, Something wrong in "
                       "ComputeNodeValue in %d term, h_NodeWeight = %e, d_NodeWeight=%e\n",
                       globalRank, localRank, globalZoneID, i, org_NodeWeight[i], h_NodeWeight[i]);
            else
                printf("Something wrong in TestNodeWeight in %d term, h_NodeWeight = %e, "
                       "d_NodeWeight=%e\n",
                       i, org_NodeWeight[i], h_NodeWeight[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "TestNodeWeight is successfull!" << endl;
#endif

    delete[] h_NodeWeight;
}

void TestCellValue(const int nTotal, const RDouble *h_data, const RDouble *d_data)
{
    size_t  sizeTotalRDouble = nTotal * sizeof(RDouble);
    RFloat *h_data_tmp       = new RDouble[nTotal];
    HANDLE_API_ERR(cudaMemcpy(h_data_tmp, d_data, sizeTotalRDouble, cudaMemcpyDeviceToHost));

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    for (int icell = 0; icell < nTotal; icell++)
    {
        //!check d_data
        abErr           = fabs(h_data_tmp[icell] - h_data[icell]);
        rtErr           = abErr / Max(fabs(h_data[icell]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in ComputeNodeValue in %d term, h_data = %e, d_data=%e\n", icell, h_data[icell],
                   h_data_tmp[icell]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            if (processNum > 1)
                printf("globalRank = %d, localRank = %d, globalZoneID = %d, Something wrong in "
                       "ComputeNodeValue in %d term, h_data = %e, d_data=%e\n",
                       globalRank, localRank, globalZoneID, icell, h_data[icell], h_data_tmp[icell]);
            else
                printf("Something wrong in ComputeNodeValue in %d term, h_data = %e, d_data=%e\n", icell, h_data[icell],
                       h_data_tmp[icell]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "TestCellValue is successfull!" << endl;
#endif

    delete[] h_data_tmp;
}

void TestComputeQTurbNodeValue(const int nTotalNode, const int nEquation, RFloat **org_qNode, const int *org_nCount)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    size_t  sizeQNodeRFloat = nEquation * nTotalNode * sizeof(RFloat);
    size_t  sizeNodeInt     = nTotalNode * sizeof(int);
    RFloat *h_qNode         = new RFloat[nEquation * nTotalNode];
    int    *h_nCount        = new int[nTotalNode];
    HANDLE_API_ERR(cudaMemcpy(h_qNode, d_qTurbNode, sizeQNodeRFloat, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_nCount, d_nCount, sizeNodeInt, cudaMemcpyDeviceToHost));

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    for (int i = 0; i < nTotalNode; i++)
    {
        for (int j = 0; j < nEquation; j++)
        {
            //!check qNode
            abErr           = fabs(h_qNode[j * nTotalNode + i] - org_qNode[j][i]);
            rtErr           = abErr / Max(fabs(org_qNode[j][i]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong in ComputeQTurbNodeValue in %d term, h_qNode = %.30e, "
                       "d_qNode=%.30e\n",
                       j * nTotalNode + i, org_qNode[j][i], h_qNode[j * nTotalNode + i]);
                printf("abErr = %.30e, rtErr=%.30e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Something wrong in ComputeQTurbNodeValue in %d term, h_qNode = %.30e, "
                       "d_qNode=%.30e\n",
                       j * nTotalNode + i, org_qNode[j][i], h_qNode[j * nTotalNode + i]);
                printf("abErr = %.30e, rtErr=%.30e\n", abErr, rtErr);
                exit(1);
            }
        }

        //!check nCount
        abErr           = fabs(h_nCount[i] - org_nCount[i]);
        rtErr           = abErr / Max(fabs(org_nCount[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in ComputeQTurbNodeValue in %d term, h_nCount = %d, d_nCount=%d\n", i,
                   org_nCount[i], h_nCount[i]);
            printf("abErr = %.30e, rtErr=%.30e\n", abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in ComputeQTurbNodeValue in %d term, h_nCount = %d, d_nCount=%d\n", i,
                   org_nCount[i], h_nCount[i]);
            printf("abErr = %.30e, rtErr=%.30e\n", abErr, rtErr);
            exit(1);
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test the kernel ComputeQTurbNodeValue successfully!" << endl;
#endif

    delete[] h_qNode;
    delete[] h_nCount;
}

void TestGPUInterPointsAllocCopy(int nIpoint, const int *org_interPoint2GlobalPoint,
                                 const int *org_cellNumberOfInterPoint, const int *org_labelOfInterPoint)
{
    using namespace GPUMemory;
    using namespace GPUGeomVariables;

    size_t pointSize                = nIpoint * sizeof(int);
    int   *h_interPoint2GlobalPoint = new int[nIpoint];
    int   *h_cellNumberOfInterPoint = new int[nIpoint];
    int   *h_labelOfInterPoint      = new int[nIpoint];
    HANDLE_API_ERR(cudaMemcpy(h_interPoint2GlobalPoint, d_interPoint2GlobalPoint, pointSize, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_cellNumberOfInterPoint, d_cellNumberOfInterPoint, pointSize, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_labelOfInterPoint, d_labelOfInterPoint, pointSize, cudaMemcpyDeviceToHost));
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    for (int i = 0; i < nIpoint; i++)
    {
        abErr           = fabs(h_interPoint2GlobalPoint[i] - org_interPoint2GlobalPoint[i]);
        rtErr           = abErr / Max(fabs(org_interPoint2GlobalPoint[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in interPoint2GlobalPoint in %d term, h_interPoint2GlobalPoint "
                   "= %d, d_interPoint2GlobalPoint=%d\n",
                   i, org_interPoint2GlobalPoint[i], h_interPoint2GlobalPoint[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in interPoint2GlobalPoint in %d term, h_interPoint2GlobalPoint "
                   "= %d, d_interPoint2GlobalPoint=%d\n",
                   i, org_interPoint2GlobalPoint[i], h_interPoint2GlobalPoint[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
    }
    for (int i = 0; i < nIpoint; i++)
    {
        abErr           = fabs(h_cellNumberOfInterPoint[i] - org_cellNumberOfInterPoint[i]);
        rtErr           = abErr / Max(fabs(org_cellNumberOfInterPoint[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in cellNumberOfInterPoint in %d term, h_cellNumberOfInterPoint "
                   "= %d, d_cellNumberOfInterPoint=%d\n",
                   i, org_cellNumberOfInterPoint[i], h_cellNumberOfInterPoint[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in cellNumberOfInterPoint in %d term, h_cellNumberOfInterPoint "
                   "= %d, d_cellNumberOfInterPoint=%d\n",
                   i, org_cellNumberOfInterPoint[i], h_cellNumberOfInterPoint[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
    }
    for (int i = 0; i < nIpoint; i++)
    {
        abErr           = fabs(h_labelOfInterPoint[i] - org_labelOfInterPoint[i]);
        rtErr           = abErr / Max(fabs(org_labelOfInterPoint[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in labelOfInterPoint in %d term, h_labelOfInterPoint = %d, "
                   "d_labelOfInterPoint=%d\n",
                   i, org_labelOfInterPoint[i], h_labelOfInterPoint[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in labelOfInterPoint in %d term, h_labelOfInterPoint = %d, "
                   "d_labelOfInterPoint=%d\n",
                   i, org_labelOfInterPoint[i], h_labelOfInterPoint[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUInterPointsAllocCopy successfully!!" << endl;
#endif
    delete[] h_interPoint2GlobalPoint;
    delete[] h_cellNumberOfInterPoint;
    delete[] h_labelOfInterPoint;
}
void TestGPUInterFace2BoundaryFace(const int nIFace, const int *org_interFace2BoundaryFace)
{
    using namespace GPUMemory;
    using namespace GPUGeomVariables;

    size_t sizeIF2BF                = nIFace * sizeof(int);
    int   *h_interFace2BoundaryFace = new int[nIFace];
    HANDLE_API_ERR(cudaMemcpy(h_interFace2BoundaryFace, d_interFace2BoundaryFace, sizeIF2BF, cudaMemcpyDeviceToHost));
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    for (int i = 0; i < nIFace; i++)
    {
        abErr           = fabs(h_interFace2BoundaryFace[i] - org_interFace2BoundaryFace[i]);
        rtErr           = abErr / Max(fabs(org_interFace2BoundaryFace[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong in GPUInterFace2BoundaryFace in %d term, "
                   "h_cellNumberOfInterPoint = %d, d_cellNumberOfInterPoint=%d\n",
                   i, org_interFace2BoundaryFace[i], h_interFace2BoundaryFace[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Something wrong in GPUInterFace2BoundaryFace in %d term, "
                   "h_cellNumberOfInterPoint = %d, d_cellNumberOfInterPoint=%d\n",
                   i, org_interFace2BoundaryFace[i], h_interFace2BoundaryFace[i]);
            printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUInterFace2BoundaryFace successfully" << endl;
#endif
    delete[] h_interFace2BoundaryFace;
}
void TestGPUUploadInterfaceValue(RFloat **org_fg, const int neqn, const int nIFace, const string &name)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;
    RFloat       *deviceIFVar;
    const RFloat *org_fg0 = org_fg[0];
    if ("q" == name)
    {
        deviceIFVar = d_fg_send_q;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("t" == name)
    {
        deviceIFVar = d_fg_send_t;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("turb::q" == name)
    {
        deviceIFVar = d_fg_send_qTurb;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("dqdx" == name)
    {
        deviceIFVar = d_fg_send_dqdx;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("dqdy" == name)
    {
        deviceIFVar = d_fg_send_dqdy;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("dqdz" == name)
    {
        deviceIFVar = d_fg_send_dqdz;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("limit" == name)
    {
        deviceIFVar = d_fg_send_limit;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("dtdx" == name)
    {
        deviceIFVar = d_fg_send_dtdx;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("dtdy" == name)
    {
        deviceIFVar = d_fg_send_dtdy;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("dtdz" == name)
    {
        deviceIFVar = d_fg_send_dtdz;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("turb::dqdx" == name)
    {
        deviceIFVar = d_fg_send_dqTurbdx;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("turb::dqdy" == name)
    {
        deviceIFVar = d_fg_send_dqTurbdy;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
    else if ("turb::dqdz" == name)
    {
        deviceIFVar = d_fg_send_dqTurbdz;
        TestFGValue(org_fg0, deviceIFVar, nIFace, neqn, name);
    }
}

void TestDeviceUploadInterfaceValueHostLargeBuffer(RFloat **org_fg, const int neqn, const int nIFace,
                                                   const string &name)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    int     nameNSID, nameTurbID;
    int     offsetBuffer;
    int     faceID, equationID;
    RFloat *HostSendLargeBuffer;
    RFloat  abErr = 0.0;
    RFloat  rtErr = 0.0;

    if ("q" == name)
    {
        nameNSID            = GetIndexBufferNS(name);
        offsetBuffer        = offsetBufferNS[nameNSID];
        HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
    }
    else if ("t" == name)
    {
        nameNSID            = GetIndexBufferNS(name);
        offsetBuffer        = offsetBufferNS[nameNSID];
        HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
    }
    else if ("turb::q" == name)
    {
        nameTurbID          = GetIndexBufferTurb(name);
        offsetBuffer        = offsetBufferTurb[nameTurbID];
        HostSendLargeBuffer = HostSendBufferTurb + offsetBuffer;
    }
    else if ("dqdx" == name)
    {
        nameNSID            = GetIndexBufferNS(name);
        offsetBuffer        = offsetBufferNS[nameNSID];
        HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
    }
    else if ("dqdy" == name)
    {
        nameNSID            = GetIndexBufferNS(name);
        offsetBuffer        = offsetBufferNS[nameNSID];
        HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
    }
    else if ("dqdz" == name)
    {
        nameNSID            = GetIndexBufferNS(name);
        offsetBuffer        = offsetBufferNS[nameNSID];
        HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
    }
    else if ("limit" == name)
    {
        nameNSID            = GetIndexBufferNS(name);
        offsetBuffer        = offsetBufferNS[nameNSID];
        HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
    }
    else if ("dtdx" == name)
    {
        nameNSID            = GetIndexBufferNS(name);
        offsetBuffer        = offsetBufferNS[nameNSID];
        HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
    }
    else if ("dtdy" == name)
    {
        nameNSID            = GetIndexBufferNS(name);
        offsetBuffer        = offsetBufferNS[nameNSID];
        HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
    }
    else if ("dtdz" == name)
    {
        nameNSID            = GetIndexBufferNS(name);
        offsetBuffer        = offsetBufferNS[nameNSID];
        HostSendLargeBuffer = HostSendBufferNS + offsetBuffer;
    }
    else if ("turb::dqdx" == name)
    {
        nameTurbID          = GetIndexBufferTurb(name);
        offsetBuffer        = offsetBufferTurb[nameTurbID];
        HostSendLargeBuffer = HostSendBufferTurb + offsetBuffer;
    }
    else if ("turb::dqdy" == name)
    {
        nameTurbID          = GetIndexBufferTurb(name);
        offsetBuffer        = offsetBufferTurb[nameTurbID];
        HostSendLargeBuffer = HostSendBufferTurb + offsetBuffer;
    }
    else if ("turb::dqdz" == name)
    {
        nameTurbID          = GetIndexBufferTurb(name);
        offsetBuffer        = offsetBufferTurb[nameTurbID];
        HostSendLargeBuffer = HostSendBufferTurb + offsetBuffer;
    }
    else
    {
        printf("Error: name does not exist\n");
        exit(1);
    }
    //!cudaMemcpyAsync(HostSendLargeBuffer, deviceSendBuffer, sizeArray, cudaMemcpyDeviceToHost, dataTransferDToH) may not finish yet. In fact, due to host calculation, it can be overlapped in CUDAUNITTEST mode.
    HANDLE_API_ERR(cudaStreamSynchronize(dataTransferDToH));
    for (faceID = 0; faceID < nIFace; faceID++)
    {
        for (equationID = 0; equationID < neqn; equationID++)
        {
            abErr           = fabs(org_fg[equationID][faceID] - HostSendLargeBuffer[equationID * nIFace + faceID]);
            rtErr           = abErr / Max(fabs(org_fg[equationID][faceID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in DeviceUploadInterfaceValueHostLargeBuffer in %d , "
                       "%dterm, org_fg0 = %.30e, HostSendLargeBuffer=%.30e\n",
                       equationID, faceID, org_fg[equationID][faceID],
                       HostSendLargeBuffer[equationID * nIFace + faceID]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in DeviceUploadInterfaceValueHostLargeBuffer in %d , "
                       "%dterm, org_fg0 = %.30e, HostSendLargeBuffer=%.30e\n",
                       equationID, faceID, org_fg[equationID][faceID],
                       HostSendLargeBuffer[equationID * nIFace + faceID]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test DeviceUploadInterfaceValueHostLargeBuffer " << name << " successfully" << endl;
#endif
}

void TestGPUUploadInterfaceValueLargeBuffer(RFloat **org_fg, const int neqn, const int nIFace, const string &name)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    int           nameNSID, nameTurbID;
    int           offsetBuffer;
    const RFloat *org_fg0 = org_fg[0];
    if ("q" == name)
    {
        nameNSID          = GetIndexBufferNS(name);
        offsetBuffer      = offsetBufferNS[nameNSID];
        size_t bufferSize = lengthBufferNS * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNS, GPUSendBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferNS, nIFace, neqn, name, offsetBuffer);
    }
    else if ("t" == name)
    {
        nameNSID          = GetIndexBufferNS(name);
        offsetBuffer      = offsetBufferNS[nameNSID];
        size_t bufferSize = lengthBufferNS * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNS, GPUSendBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferNS, nIFace, neqn, name, offsetBuffer);
    }
    else if ("turb::q" == name)
    {
        nameTurbID        = GetIndexBufferTurb(name);
        offsetBuffer      = offsetBufferTurb[nameTurbID];
        size_t bufferSize = lengthBufferTurb * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferTurb, GPUSendBufferTurb, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferTurb, nIFace, neqn, name, offsetBuffer);
    }
    else if ("dqdx" == name)
    {
        nameNSID          = GetIndexBufferNS(name);
        offsetBuffer      = offsetBufferNS[nameNSID];
        size_t bufferSize = lengthBufferNS * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNS, GPUSendBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferNS, nIFace, neqn, name, offsetBuffer);
    }
    else if ("dqdy" == name)
    {
        nameNSID          = GetIndexBufferNS(name);
        offsetBuffer      = offsetBufferNS[nameNSID];
        size_t bufferSize = lengthBufferNS * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNS, GPUSendBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferNS, nIFace, neqn, name, offsetBuffer);
    }
    else if ("dqdz" == name)
    {
        nameNSID          = GetIndexBufferNS(name);
        offsetBuffer      = offsetBufferNS[nameNSID];
        size_t bufferSize = lengthBufferNS * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNS, GPUSendBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferNS, nIFace, neqn, name, offsetBuffer);
    }
    else if ("limit" == name)
    {
        nameNSID          = GetIndexBufferNS(name);
        offsetBuffer      = offsetBufferNS[nameNSID];
        size_t bufferSize = lengthBufferNS * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNS, GPUSendBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferNS, nIFace, neqn, name, offsetBuffer);
    }
    else if ("dtdx" == name)
    {
        nameNSID          = GetIndexBufferNS(name);
        offsetBuffer      = offsetBufferNS[nameNSID];
        size_t bufferSize = lengthBufferNS * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNS, GPUSendBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferNS, nIFace, neqn, name, offsetBuffer);
    }
    else if ("dtdy" == name)
    {
        nameNSID          = GetIndexBufferNS(name);
        offsetBuffer      = offsetBufferNS[nameNSID];
        size_t bufferSize = lengthBufferNS * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNS, GPUSendBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferNS, nIFace, neqn, name, offsetBuffer);
    }
    else if ("dtdz" == name)
    {
        nameNSID          = GetIndexBufferNS(name);
        offsetBuffer      = offsetBufferNS[nameNSID];
        size_t bufferSize = lengthBufferNS * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferNS, GPUSendBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferNS, nIFace, neqn, name, offsetBuffer);
    }
    else if ("turb::dqdx" == name)
    {
        nameTurbID        = GetIndexBufferTurb(name);
        offsetBuffer      = offsetBufferTurb[nameTurbID];
        size_t bufferSize = lengthBufferTurb * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferTurb, GPUSendBufferTurb, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferTurb, nIFace, neqn, name, offsetBuffer);
    }
    else if ("turb::dqdy" == name)
    {
        nameTurbID        = GetIndexBufferTurb(name);
        offsetBuffer      = offsetBufferTurb[nameTurbID];
        size_t bufferSize = lengthBufferTurb * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferTurb, GPUSendBufferTurb, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferTurb, nIFace, neqn, name, offsetBuffer);
    }
    else if ("turb::dqdz" == name)
    {
        nameTurbID        = GetIndexBufferTurb(name);
        offsetBuffer      = offsetBufferTurb[nameTurbID];
        size_t bufferSize = lengthBufferTurb * sizeof(double);
        HANDLE_API_ERR(cudaMemcpy(HostSendBufferTurb, GPUSendBufferTurb, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueLargeBuffer(org_fg0, HostSendBufferTurb, nIFace, neqn, name, offsetBuffer);
    }
}

void TestFGValueLargeBuffer(const RFloat *org_fg0, const double *HostSendBuffer, const int nIFace, const int neqn,
                            const string &name, const int offset)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    int    iFace, m;
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    for (iFace = 0; iFace < nIFace; iFace++)
    {
        for (m = 0; m < neqn; m++)
        {
            abErr           = fabs(org_fg0[m * nIFace + iFace] - HostSendBuffer[offset + m * nIFace + iFace]);
            rtErr           = abErr / Max(fabs(org_fg0[m * nIFace + iFace]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in GPUUploadInterfaceValue in %d , %dterm, org_fg0 = "
                       "%.30e, HostSendBuffer=%.30e\n",
                       m, iFace, org_fg0[m * nIFace + iFace], HostSendBuffer[offset + m * nIFace + iFace]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in GPUUploadInterfaceValue in %d , %dterm, org_fg0 = "
                       "%.30e, HostSendBuffer=%.30e\n",
                       m, iFace, org_fg0[m * nIFace + iFace], HostSendBuffer[offset + m * nIFace + iFace]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUUploadInterfaceValue in " << name << " successfully" << endl;
#endif
}

void TestFGValue(const RFloat *org_fg0, const RFloat *deviceIFVar, const int nIFace, const int neqn, const string &name)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    int     iFace, m;
    RFloat  abErr     = 0.0;
    RFloat  rtErr     = 0.0;
    RFloat *h_fg_send = new RFloat[nIFace * neqn];
    size_t  sizeIF    = nIFace * neqn * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(h_fg_send, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
    for (iFace = 0; iFace < nIFace; iFace++)
    {
        for (m = 0; m < neqn; m++)
        {
            abErr           = fabs(org_fg0[m * nIFace + iFace] - h_fg_send[m * nIFace + iFace]);
            rtErr           = abErr / Max(fabs(org_fg0[m * nIFace + iFace]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in GPUUploadInterfaceValue in %d , %dterm, org_fg0 = "
                       "%.30e, d_fg_send=%.30e\n",
                       m, iFace, org_fg0[m * nIFace + iFace], h_fg_send[m * nIFace + iFace]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in GPUUploadInterfaceValue in %d , %dterm, org_fg0 = "
                       "%.30e, d_fg_send=%.30e\n",
                       m, iFace, org_fg0[m * nIFace + iFace], h_fg_send[m * nIFace + iFace]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUUploadInterfaceValue in " << name << " successfully" << endl;
#endif
    delete[] h_fg_send;
}
void TestUploadInterfaceValueDeviceToHost(const int *interFace2BoundaryFace, const int *leftCellOfFace, RFloat **fg,
                                          RFloat **f, const int nIFace, const int neqn, const string &name)
{
    RFloat abErr = 0;
    RFloat rtErr = 0;
    if (("q" != name) && ("t" != name) && ("turb::q" != name) && ("dqdx" != name) && ("dqdy" != name)
        && ("dqdz" != name) && ("limit" != name) && ("dtdx" != name) && ("dtdy" != name) && ("dtdz" != name)
        && ("turb::dqdx" != name) && ("turb::dqdy" != name) && ("turb::dqdz" != name))
        return;
    for (int iFace = 0; iFace < nIFace; ++iFace)
    {
        int sourceCell = leftCellOfFace[interFace2BoundaryFace[iFace]];
        for (int m = 0; m < neqn; ++m)
        {
            abErr           = fabs(fg[m][iFace] - f[m][sourceCell]);
            rtErr           = abErr / Max(fabs(f[m][sourceCell]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in UploadInterfaceValueDeviceToHost in %d , %dterm, f = "
                       "%.30e, fg =%.30e\n",
                       m, iFace, f[m][sourceCell], fg[m][iFace]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in UploadInterfaceValueDeviceToHost in %d , %dterm, f = "
                       "%.30e, fg =%.30e\n",
                       m, iFace, f[m][sourceCell], fg[m][iFace]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test UploadInterfaceValueDeviceToHost in " << name << " successfully" << endl;
#endif
}

void TestGPUDownloadInterfaceValue(const int *interFace2BoundaryFace, const int *rightCellOfFace, RFloat **f,
                                   const int nIFace, const int neqn, const int nTotal, const string &name)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat *deviceFieldVar;
    if ("q" == name)
    {
        deviceFieldVar = d_q_ns;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("t" == name)
    {
        deviceFieldVar = d_t_proxy;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("turb::q" == name)
    {
        deviceFieldVar = d_q_turb_proxy;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("dqdx" == name)
    {
        deviceFieldVar = d_dqdx_proxy;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("dqdy" == name)
    {
        deviceFieldVar = d_dqdy_proxy;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("dqdz" == name)
    {
        deviceFieldVar = d_dqdz_proxy;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("limit" == name)
    {
        deviceFieldVar = d_limit;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("dtdx" == name)
    {
        deviceFieldVar = d_dtdx_proxy;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("dtdy" == name)
    {
        deviceFieldVar = d_dtdy_proxy;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("dtdz" == name)
    {
        deviceFieldVar = d_dtdz_proxy;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("turb::dqdx" == name)
    {
        deviceFieldVar = d_dq_turbdx_proxy;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("turb::dqdy" == name)
    {
        deviceFieldVar = d_dq_turbdy_proxy;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
    else if ("turb::dqdz" == name)
    {
        deviceFieldVar = d_dq_turbdz_proxy;
        TestFieldValue(f, deviceFieldVar, interFace2BoundaryFace, rightCellOfFace, nIFace, nTotal, neqn, name);
    }
}

void TestFieldValue(RFloat **f, const RFloat *deviceFieldVar, const int *interFace2BoundaryFace,
                    const int *rightCellOfFace, const int nIFace, const int nTotal, const int neqn, const string &name)
{
    size_t  sizeFieldVar = neqn * nTotal * sizeof(RFloat);
    RFloat *fg           = new RFloat[neqn * nTotal];
    HANDLE_API_ERR(cudaMemcpy(fg, deviceFieldVar, sizeFieldVar, cudaMemcpyDeviceToHost));

    RFloat abErr = 0;
    RFloat rtErr = 0;
    for (int iFace = 0; iFace < nIFace; ++iFace)
    {
        int sourceCell = rightCellOfFace[interFace2BoundaryFace[iFace]];
        for (int m = 0; m < neqn; ++m)
        {
            abErr           = fabs(fg[m * nTotal + sourceCell] - f[m][sourceCell]);
            rtErr           = abErr / Max(fabs(f[m][sourceCell]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in GPUDownloadInterfaceValue in %d , %dterm, f = %.30e, fg "
                       "=%.30e\n",
                       m, iFace, f[m][sourceCell], fg[m * nTotal + sourceCell]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in GPUDownloadInterfaceValue in %d , %dterm, f = %.30e, fg "
                       "=%.30e\n",
                       m, iFace, f[m][sourceCell], fg[m * nTotal + sourceCell]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUDownloadInterfaceValue in " << name << " successfully" << endl;
#endif
    delete[] fg;
}

void TestDeviceUploadInterpointValueHostLargeBuffer(RFloat **org_fg, const int neqn, const int nIPoint,
                                                    const string &name)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    int     nameNSID, nameTurbID;
    int     offsetBuffer;
    int     pointID, equationID;
    RFloat *HostSendLargeBuffer;
    RFloat  abErr = 0.0;
    RFloat  rtErr = 0.0;

    if ("qnode" == name)
    {
        nameNSID            = GetIndexInterpointBufferNS(name);
        offsetBuffer        = offsetInterpointBufferNS[nameNSID];
        HostSendLargeBuffer = HostSendInterpointBufferNS + offsetBuffer;
    }
    else if ("tnode" == name)
    {
        nameNSID            = GetIndexInterpointBufferNS(name);
        offsetBuffer        = offsetInterpointBufferNS[nameNSID];
        HostSendLargeBuffer = HostSendInterpointBufferNS + offsetBuffer;
    }
    else if ("qTurbNode" == name)
    {
        nameTurbID          = GetIndexInterpointBufferTurb(name);
        offsetBuffer        = offsetInterpointBufferTurb[nameTurbID];
        HostSendLargeBuffer = HostSendInterpointBufferTurb + offsetBuffer;
    }

    HANDLE_API_ERR(cudaStreamSynchronize(dataTransferDToHInterpoint));

    for (pointID = 0; pointID < nIPoint; pointID++)
    {
        for (equationID = 0; equationID < neqn; equationID++)
        {
            abErr           = fabs(org_fg[equationID][pointID] - HostSendLargeBuffer[equationID * nIPoint + pointID]);
            rtErr           = abErr / Max(fabs(org_fg[equationID][pointID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in DeviceUploadInterpointValueHostLargeBuffer in %d , "
                       "%dterm, org_fg0 = %.30e, HostSendLargeBuffer=%.30e\n",
                       equationID, pointID, org_fg[equationID][pointID],
                       HostSendLargeBuffer[equationID * nIPoint + pointID]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in DeviceUploadInterpointValueHostLargeBuffer in %d , "
                       "%dterm, org_fg0 = %.30e, HostSendLargeBuffer=%.30e\n",
                       equationID, pointID, org_fg[equationID][pointID],
                       HostSendLargeBuffer[equationID * nIPoint + pointID]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test DeviceUploadInterpointValueHostLargeBuffer " << name << " successfully" << endl;
#endif
}

void TestGPUUploadInterpointValueLargeBuffer(RFloat **org_fg, const int neqn, const int nIPoint, const string &name)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;
    int           nameNSID, nameTurbID;
    int           offsetBuffer;
    const RFloat *org_fg0 = org_fg[0];

    if ("qnode" == name)
    {
        nameNSID          = GetIndexInterpointBufferNS(name);
        offsetBuffer      = offsetInterpointBufferNS[nameNSID];
        size_t bufferSize = lengthInterpointBufferNS * sizeof(double);
        HANDLE_API_ERR(
            cudaMemcpy(HostSendInterpointBufferNS, GPUSendInterpointBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueInterpointLargeBuffer(org_fg0, HostSendInterpointBufferNS, nIPoint, neqn, name, offsetBuffer);
    }
    else if ("tnode" == name)
    {
        nameNSID          = GetIndexInterpointBufferNS(name);
        offsetBuffer      = offsetInterpointBufferNS[nameNSID];
        size_t bufferSize = lengthInterpointBufferNS * sizeof(double);
        HANDLE_API_ERR(
            cudaMemcpy(HostSendInterpointBufferNS, GPUSendInterpointBufferNS, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueInterpointLargeBuffer(org_fg0, HostSendInterpointBufferNS, nIPoint, neqn, name, offsetBuffer);
    }
    else if ("qTurbNode" == name)
    {
        nameTurbID        = GetIndexInterpointBufferTurb(name);
        offsetBuffer      = offsetInterpointBufferTurb[nameTurbID];
        size_t bufferSize = lengthInterpointBufferTurb * sizeof(double);
        HANDLE_API_ERR(
            cudaMemcpy(HostSendInterpointBufferTurb, GPUSendInterpointBufferTurb, bufferSize, cudaMemcpyDeviceToHost));
        TestFGValueInterpointLargeBuffer(org_fg0, HostSendInterpointBufferTurb, nIPoint, neqn, name, offsetBuffer);
    }

#ifdef UNITTESTOUTPUT
    cout << "Test GPUUploadInterpointValueLargeBuffer in " << name << " successfully" << endl;
#endif
}

void TestFGValueInterpointLargeBuffer(const RFloat *org_fg0, const double *HostSendBuffer, const int nIPoint,
                                      const int neqn, const string &name, const int offset)
{
    int    iPoint, m;
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    for (iPoint = 0; iPoint < nIPoint; iPoint++)
    {
        for (m = 0; m < neqn; m++)
        {
            abErr           = fabs(org_fg0[m * nIPoint + iPoint] - HostSendBuffer[offset + m * nIPoint + iPoint]);
            rtErr           = abErr / Max(fabs(org_fg0[m * nIPoint + iPoint]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: TestFGValueInterpointLargeBuffer in m = %d, iPoint = %d, name = %s, "
                       "org_fg0 = %.30e, HostSendBuffer = %.30e, abErr = %.30e, rtErr = %.30e\n",
                       m, iPoint, name.c_str(), org_fg0[m * nIPoint + iPoint],
                       HostSendBuffer[offset + m * nIPoint + iPoint], abErr, rtErr);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: TestFGValueInterpointLargeBuffer in m = %d, iPoint = %d, name = %s, "
                       "org_fg0 = %.30e, HostSendBuffer = %.30e, abErr = %.30e, rtErr = %.30e\n",
                       m, iPoint, name.c_str(), org_fg0[m * nIPoint + iPoint],
                       HostSendBuffer[offset + m * nIPoint + iPoint], abErr, rtErr);
            }
        }
    }
}

void TestGPUUploadInterpointValue(RFloat **org_fg, const int neqn, const int nIPoint, const string &name)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    using namespace GPUGeomVariables;
    RFloat       *deviceIFVar;
    const RFloat *org_fg0 = org_fg[0];
    if ("qnode" == name)
    {
        deviceIFVar = d_fg_send_qNode;
        TestFGValueInterpoint(org_fg0, deviceIFVar, nIPoint, neqn, name);
    }
    else if ("tnode" == name)
    {
        deviceIFVar = d_fg_send_tNode;
        TestFGValueInterpoint(org_fg0, deviceIFVar, nIPoint, neqn, name);
    }
    else if ("qTurbNode" == name)
    {
        deviceIFVar = d_fg_send_qTurbNode;
        TestFGValueInterpoint(org_fg0, deviceIFVar, nIPoint, neqn, name);
    }
}
void TestFGValueInterpoint(const RFloat *org_fg0, const RFloat *deviceIFVar, const int nIFace, const int neqn,
                           const string &name)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    int     iFace, m;
    RFloat  abErr     = 0.0;
    RFloat  rtErr     = 0.0;
    RFloat *h_fg_send = new RFloat[nIFace * neqn];
    size_t  sizeIF    = nIFace * neqn * sizeof(RFloat);
    HANDLE_API_ERR(cudaMemcpy(h_fg_send, deviceIFVar, sizeIF, cudaMemcpyDeviceToHost));
    for (iFace = 0; iFace < nIFace; iFace++)
    {
        for (m = 0; m < neqn; m++)
        {
            abErr           = fabs(org_fg0[m * nIFace + iFace] - h_fg_send[m * nIFace + iFace]);
            rtErr           = abErr / Max(fabs(org_fg0[m * nIFace + iFace]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in GPUUploadInterfaceValue in %d , %dterm, org_fg0 = "
                       "%.30e, h_fg_send=%.30e\n",
                       m, iFace, org_fg0[m * nIFace + iFace], h_fg_send[m * nIFace + iFace]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in GPUUploadInterfaceValue in %d , %dterm, org_fg0 = "
                       "%.30e, h_fg_send=%.30e\n",
                       m, iFace, org_fg0[m * nIFace + iFace], h_fg_send[m * nIFace + iFace]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUUploadInterpointValue in " << name << " successfully" << endl;
#endif
    delete[] h_fg_send;
}
void TestUploadInterpointValueDeviceToHost(const int *interPoint2GlobalPoint, RFloat **fg, RFloat **f, const int neqn,
                                           const int nIPoint, const string &name)
{
    RFloat abErr = 0;
    RFloat rtErr = 0;
    if (("qnode" != name) && ("tnode" != name) && ("qTurbNode" != name)) return;
    for (int iPoint = 0; iPoint < nIPoint; ++iPoint)
    {
        int sourcePoint = interPoint2GlobalPoint[iPoint];
        for (int m = 0; m < neqn; ++m)
        {
            abErr           = fabs(fg[m][iPoint] - f[m][sourcePoint]);
            rtErr           = abErr / Max(fabs(f[m][sourcePoint]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in UploadInterpointValueDeviceToHost in %d , %dterm, f = "
                       "%.30e, fg =%.30e\n",
                       m, iPoint, f[m][sourcePoint], fg[m][iPoint]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Test " << name << "fails" << endl;
                printf("Something wrong in UploadInterpointValueDeviceToHost in %d , %dterm, f = "
                       "%.30e, fg =%.30e\n",
                       m, iPoint, f[m][sourcePoint], fg[m][iPoint]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }
    }

#ifdef UNITTESTOUTPUT
    cout << "Test UploadInterpointValueDeviceToHost in " << name << " successfully" << endl;
#endif
}

void TestGPUDownloadInterpointValue(const int *interPoint2GlobalPoint, RFloat **f, const int nIPoint, const int neqn,
                                    const int nTotalNode, const string &name)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat *deviceFieldVar;
    if ("qnode" == name)
    {
        //!deviceFieldVar = d_qNode;
        deviceFieldVar = d_qInterPoint;
        TestFieldValueInterpoint(f, deviceFieldVar, interPoint2GlobalPoint, nIPoint, nTotalNode, neqn, name);
    }
    else if ("tnode" == name)
    {
        //!deviceFieldVar = d_tNode;
        deviceFieldVar = d_tInterPoint;
        TestFieldValueInterpoint(f, deviceFieldVar, interPoint2GlobalPoint, nIPoint, nTotalNode, neqn, name);
    }
    else if ("qTurbNode" == name)
    {
        //!deviceFieldVar = d_qTurbNode;
        deviceFieldVar = d_qTurbInterPoint;
        TestFieldValueInterpoint(f, deviceFieldVar, interPoint2GlobalPoint, nIPoint, nTotalNode, neqn, name);
    }
}

void TestFieldValueInterpoint(RFloat **f, const RFloat *deviceFieldVar, const int *interPoint2GlobalPoint,
                              const int nIPoint, const int nTotalNode, const int neqn, const string &name)
{
    //!size_t sizeFieldVar = neqn * nTotalNode * sizeof(RFloat);
    size_t sizeFieldVar = neqn * nIPoint * sizeof(RFloat);
    //!RFloat * fg = new RFloat[neqn * nTotalNode];
    RFloat *fg = new RFloat[neqn * nIPoint];
    HANDLE_API_ERR(cudaMemcpy(fg, deviceFieldVar, sizeFieldVar, cudaMemcpyDeviceToHost));

    RFloat abErr = 0;
    RFloat rtErr = 0;
    for (int iPoint = 0; iPoint < nIPoint; ++iPoint)
    {
        for (int m = 0; m < neqn; ++m)
        {
            //!abErr = fabs(f[m][iPoint] - fg[m * nTotalNode + iPoint]);
            abErr           = fabs(f[m][iPoint] - fg[m * nIPoint + iPoint]);
            rtErr           = abErr / Max(fabs(f[m][iPoint]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Test " << name << "fails" << endl;
                //!printf("Something wrong in GPUDownloadInterpointValue in %d , %dterm, f = %.30e, fg =%.30e\n", m, iPoint, f[m][iPoint], fg[m * nTotalNode + iPoint]);
                printf("Something wrong in GPUDownloadInterpointValue in %d , %dterm, f = %.30e, "
                       "fg =%.30e\n",
                       m, iPoint, f[m][iPoint], fg[m * nIPoint + iPoint]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Test " << name << "fails" << endl;
                //!printf("Something wrong in GPUDownloadInterpointValue in %d , %dterm, f = %.30e, fg =%.30e\n", m, iPoint, f[m][iPoint], fg[m * nTotalNode + iPoint]);
                printf("Something wrong in GPUDownloadInterpointValue in %d , %dterm, f = %.30e, "
                       "fg =%.30e\n",
                       m, iPoint, f[m][iPoint], fg[m * nIPoint + iPoint]);
                printf("abErr = %e, rtErr=%e\n", abErr, rtErr);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUDownloadInterpointValue in " << name << " successfully" << endl;
#endif
    delete[] fg;
}

void TestGPUModifyNodeValue(RFloat **org_qNode, RFloat **org_tNode, const int *interPoint2GlobalPoint,
                            const int nIPoint, const int neqn, const int nTotalNode)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    size_t  sizeQNode = neqn * nTotalNode * sizeof(RFloat);
    size_t  sizeTNode = 1 * nTotalNode * sizeof(RFloat);
    int     iPoint, globalPoint, m;
    RFloat *h_qNode = new RFloat[neqn * nTotalNode];
    RFloat *h_tNode = new RFloat[1 * nTotalNode];
    HANDLE_API_ERR(cudaMemcpy(h_qNode, d_qNode, sizeQNode, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(h_tNode, d_tNode, sizeTNode, cudaMemcpyDeviceToHost));
    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;
    for (iPoint = 0; iPoint < nIPoint; ++iPoint)
    {
        globalPoint = interPoint2GlobalPoint[iPoint];
        for (m = 0; m < neqn; m++)
        {
            abErr           = fabs(org_qNode[m][globalPoint] - h_qNode[m * nTotalNode + globalPoint]);
            rtErr           = abErr / Max(fabs(org_qNode[m][globalPoint]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with TestGPUModifyNodeValue qNode in iPoint %d, equation "
                       "%d, globalPoint %d, h_qNode = %.30e, d_qNode= %.30e\n",
                       iPoint, m, globalPoint, org_qNode[m][globalPoint], h_qNode[m * nTotalNode + globalPoint]);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                if (processNum > 1)
                    printf("In GloablRank %d, something wrong with TestGPUModifyNodeValue qNode in "
                           "iPoint %d, equation %d, globalPoint %d, h_qNode = %.30e, d_qNode= %.30e\n",
                           globalRank, iPoint, m, globalPoint, org_qNode[m][globalPoint],
                           h_qNode[m * nTotalNode + globalPoint]);
                else
                    printf("Something wrong with TestGPUModifyNodeValue qNode in iPoint %d, "
                           "equation %d, globalPoint %d, h_qNode = %.30e, d_qNode= %.30e\n",
                           iPoint, m, globalPoint, org_qNode[m][globalPoint], h_qNode[m * nTotalNode + globalPoint]);
                exit(1);
            }
        }
        abErr           = fabs(org_tNode[0][globalPoint] - h_tNode[0 * nTotalNode + globalPoint]);
        rtErr           = abErr / Max(fabs(org_tNode[0][globalPoint]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Something wrong with TestGPUModifyNodeValue tNode in iPoint %d, equation %d, "
                   "globalPoint %d, h_tNode = %.30e, d_tNode= %.30e\n",
                   iPoint, m, globalPoint, org_tNode[0][globalPoint], h_tNode[0 * nTotalNode + globalPoint]);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            if (processNum > 1)
                printf("In GlobalRank %d, Something wrong with TestGPUModifyNodeValue tNode in "
                       "iPoint %d, equation %d, globalPoint %d, h_tNode = %.30e, d_tNode= %.30e\n",
                       globalRank, iPoint, m, globalPoint, org_tNode[0][globalPoint],
                       h_tNode[0 * nTotalNode + globalPoint]);
            else
                printf("Something wrong with TestGPUModifyNodeValue tNode in iPoint %d, equation "
                       "%d, globalPoint %d, h_tNode = %.30e, d_tNode= %.30e\n",
                       iPoint, m, globalPoint, org_tNode[0][globalPoint], h_tNode[0 * nTotalNode + globalPoint]);
            exit(1);
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUModifyNodeValue successfully" << endl;
#endif
    delete[] h_qNode;
    delete[] h_tNode;
}

void TestModifyQTurbNodeValue(const int *interPoint2GlobalPoint, RFloat **org_qTurbNode, const int nIPoint,
                              const int neqn, const int nTotalNode)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    size_t  sizeQNode = neqn * nTotalNode * sizeof(RFloat);
    int     iPoint, globalPoint, m;
    RFloat  abErr       = 0.0;
    RFloat  rtErr       = 0.0;
    RFloat *h_qTurbNode = new RFloat[neqn * nTotalNode];
    HANDLE_API_ERR(cudaMemcpy(h_qTurbNode, d_qTurbNode, sizeQNode, cudaMemcpyDeviceToHost));
    for (iPoint = 0; iPoint < nIPoint; ++iPoint)
    {
        globalPoint = interPoint2GlobalPoint[iPoint];
        for (m = 0; m < neqn; m++)
        {
            abErr           = fabs(org_qTurbNode[m][globalPoint] - h_qTurbNode[m * nTotalNode + globalPoint]);
            rtErr           = abErr / Max(fabs(org_qTurbNode[m][globalPoint]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Something wrong with TestGPUModifyQTurbNodeValue, iPoint %d, equation %d, "
                       "globalPoint %d term, h_qTurbNode = %.30e, d_qTurbNode= %.30e\n",
                       iPoint, m, globalPoint, org_qTurbNode[m][globalPoint],
                       h_qTurbNode[m * nTotalNode + globalPoint]);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                if (processNum > 1)
                    printf("In GlobalRank %d, Something wrong with TestGPUModifyQTurbNodeValue, "
                           "iPoint %d, equation %d, globalPoint %d term, h_qTurbNode = %.30e, "
                           "d_qTurbNode= %.30e\n",
                           globalRank, iPoint, m, globalPoint, org_qTurbNode[m][globalPoint],
                           h_qTurbNode[m * nTotalNode + globalPoint]);
                else
                    printf("Something wrong with TestGPUModifyQTurbNodeValue, iPoint %d, equation "
                           "%d, globalPoint %d term, h_qTurbNode = %.30e, d_qTurbNode= %.30e\n",
                           iPoint, m, globalPoint, org_qTurbNode[m][globalPoint],
                           h_qTurbNode[m * nTotalNode + globalPoint]);
                exit(1);
            }
        }
    }
#ifdef UNITTESTOUTPUT
    cout << "Test GPUModifyQTurbNodeValue successfully" << endl;
#endif
    delete[] h_qTurbNode;
}

void TestCompressAndTransferInterpointDataTurb(const int iZone, const int jZone, const int nameID, const int ngbID,
                                               const int nIPoint, const int nIPointOfNeighbor,
                                               const int offsetNgbLengthBuffer, const int *pointIndexForSend,
                                               const double *sendBufferNSAWARE, const double *sendBufferNS)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat abErr        = 0.0;
    RFloat rtErr        = 0.0;
    int    offsetPoints = offsetPointIndex[ngbID];
    int    offsetNgb    = offsetMPIDataInterpointTurb[ngbID];

    const string &name             = (*nameTurbInterpointBuffer)[nameID];
    int           offsetVar        = offsetInterpointBufferTurb[nameID];
    int           dimData          = dimInterpointBufferTurb[nameID];
    int           offsetVarInAWARE = offsetNgb + offsetNgbLengthBuffer;

    for (int iPoint = 0; iPoint < nIPointOfNeighbor; iPoint++)
    {
        //!a bug exists here! faceIndexForSend is just belong to ngbID
        //!int faceID = faceIndexForSend[offsetFaces + iFace];
        int pointID      = pointIndexForSend[iPoint];
        int pointIDAWARE = iPoint;
        for (int iDim = 0; iDim < dimData; iDim++)
        {
            double GPUAware = sendBufferNSAWARE[offsetVarInAWARE + iDim * nIPointOfNeighbor + pointIDAWARE];
            double GPUOrg   = sendBufferNS[offsetVar + iDim * nIPoint + pointID];
            abErr           = fabs(GPUAware - GPUOrg);
            rtErr           = abErr / Max(fabs(GPUOrg), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in CompressAndTransferInterpointDataTurb, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, iPoint = %d, pointID = "
                       "%d, GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, iPoint, pointID, GPUAware, GPUOrg);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in CompressAndTransferInterpointDataTurb, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, iPoint = %d, pointID = "
                       "%d, GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, iPoint, pointID, GPUAware, GPUOrg);
                exit(1);
            }
        }
    }
    //!cout<<"Test send data "<<name<<" successfully in iZone "<<iZone<<" and jZone "<<jZone<<" globalRank "<<globalRank<<" and localRank "<<localRank<<" !"<<endl;
}

void TestCompressAndTransferInterpointDataNS(const int iZone, const int jZone, const int nameID, const int ngbID,
                                             const int nIPoint, const int nIPointOfNeighbor,
                                             const int offsetNgbLengthBuffer, const int *pointIndexForSend,
                                             const double *sendBufferNSAWARE, const double *sendBufferNS)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat abErr        = 0.0;
    RFloat rtErr        = 0.0;
    int    offsetPoints = offsetPointIndex[ngbID];
    int    offsetNgb    = offsetMPIDataInterpointNS[ngbID];

    const string &name             = (*nameNSInterpointBuffer)[nameID];
    int           offsetVar        = offsetInterpointBufferNS[nameID];
    int           dimData          = dimInterpointBufferNS[nameID];
    int           offsetVarInAWARE = offsetNgb + offsetNgbLengthBuffer;

    for (int iPoint = 0; iPoint < nIPointOfNeighbor; iPoint++)
    {
        //!a bug exists here! faceIndexForSend is just belong to ngbID
        //!int faceID = faceIndexForSend[offsetFaces + iFace];
        int pointID      = pointIndexForSend[iPoint];
        int pointIDAWARE = iPoint;
        for (int iDim = 0; iDim < dimData; iDim++)
        {
            double GPUAware = sendBufferNSAWARE[offsetVarInAWARE + iDim * nIPointOfNeighbor + pointIDAWARE];
            double GPUOrg   = sendBufferNS[offsetVar + iDim * nIPoint + pointID];
            abErr           = fabs(GPUAware - GPUOrg);
            rtErr           = abErr / Max(fabs(GPUOrg), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in CompressAndTransferInterpointDataNS, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, iPoint = %d, pointID = "
                       "%d, GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, iPoint, pointID, GPUAware, GPUOrg);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in CompressAndTransferInterpointDataNS, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, iPoint = %d, pointID = "
                       "%d, GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, iPoint, pointID, GPUAware, GPUOrg);
                exit(1);
            }
        }
    }
    //!cout<<"Test send data "<<name<<" successfully in iZone "<<iZone<<" and jZone "<<jZone<<" globalRank "<<globalRank<<" and localRank "<<localRank<<" !"<<endl;
}

void TestCompressAndTransferInterpolateDataTurb(const int iZone, const int jZone, const int nameID, const int ngbID,
                                                const int nIFace, const int nIFaceOfNeighbor,
                                                const int offsetNgbLengthBuffer, const int *faceIndexForSend,
                                                const double *sendBufferTurbAWARE, const double *sendBufferTurb)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat abErr       = 0.0;
    RFloat rtErr       = 0.0;
    int    offsetFaces = offsetFaceIndex[ngbID];
    int    offsetNgb   = offsetMPIDataTurb[ngbID];

    const string &name             = (*nameTurbBuffer)[nameID];
    int           offsetVar        = offsetBufferTurb[nameID];
    int           dimData          = dimBufferTurb[nameID];
    int           offsetVarInAWARE = offsetNgb + offsetNgbLengthBuffer;

    for (int iFace = 0; iFace < nIFaceOfNeighbor; iFace++)
    {
        //!a bug exists here! faceIndexForSend is just belong to ngbID
        //!int faceID = faceIndexForSend[offsetFaces + iFace];
        int faceID      = faceIndexForSend[iFace];
        int faceIDAWARE = iFace;
        for (int iDim = 0; iDim < dimData; iDim++)
        {
            double GPUAware = sendBufferTurbAWARE[offsetVarInAWARE + iDim * nIFaceOfNeighbor + faceIDAWARE];
            double GPUOrg   = sendBufferTurb[offsetVar + iDim * nIFace + faceID];
            abErr           = fabs(GPUAware - GPUOrg);
            rtErr           = abErr / Max(fabs(GPUOrg), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in CompressAndTransferInterpolateDataTurb, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, iFace = %d, faceID = %d, "
                       "GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, iFace, faceID, GPUAware, GPUOrg);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in CompressAndTransferInterpolateDataTurb, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, iFace = %d, faceID = %d, "
                       "GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, iFace, faceID, GPUAware, GPUOrg);
                exit(1);
            }
        }
    }
    //!cout<<"Test send data "<<name<<" successfully in iZone "<<iZone<<" and jZone "<<jZone<<" globalRank "<<globalRank<<" and localRank "<<localRank<<" !"<<endl;
}

void TestCompressAndTransferInterpolateDataNS(const int iZone, const int jZone, const int nameID, const int ngbID,
                                              const int nIFace, const int nIFaceOfNeighbor,
                                              const int offsetNgbLengthBuffer, const int *faceIndexForSend,
                                              const double *sendBufferNSAWARE, const double *sendBufferNS)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;
    RFloat abErr       = 0.0;
    RFloat rtErr       = 0.0;
    int    offsetFaces = offsetFaceIndex[ngbID];
    int    offsetNgb   = offsetMPIDataNS[ngbID];

    const string &name             = (*nameNSBuffer)[nameID];
    int           offsetVar        = offsetBufferNS[nameID];
    int           dimData          = dimBufferNS[nameID];
    int           offsetVarInAWARE = offsetNgb + offsetNgbLengthBuffer;

    for (int iFace = 0; iFace < nIFaceOfNeighbor; iFace++)
    {
        //!a bug exists here! faceIndexForSend is just belong to ngbID
        //!int faceID = faceIndexForSend[offsetFaces + iFace];
        int faceID      = faceIndexForSend[iFace];
        int faceIDAWARE = iFace;
        for (int iDim = 0; iDim < dimData; iDim++)
        {
            double GPUAware = sendBufferNSAWARE[offsetVarInAWARE + iDim * nIFaceOfNeighbor + faceIDAWARE];
            double GPUOrg   = sendBufferNS[offsetVar + iDim * nIFace + faceID];
            abErr           = fabs(GPUAware - GPUOrg);
            rtErr           = abErr / Max(fabs(GPUOrg), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in CompressAndTransferInterpolateDataNS, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, iFace = %d, faceID = %d, "
                       "GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, iFace, faceID, GPUAware, GPUOrg);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in CompressAndTransferInterpolateDataNS, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, iFace = %d, faceID = %d, "
                       "GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, iFace, faceID, GPUAware, GPUOrg);
                exit(1);
            }
        }
    }
    //!cout<<"Test send data "<<name<<" successfully in iZone "<<iZone<<" and jZone "<<jZone<<" globalRank "<<globalRank<<" and localRank "<<localRank<<" !"<<endl;
}

void TestGPUMPIDataDecompressInterpointTurb(const int iZone, const int jZone, const int nameID, const int ngbID,
                                            const int nIPoint, const int nIPointOfNeighbor,
                                            const int offsetNgbLengthBuffer, const int *pointIndexForRecv,
                                            const double *recvBufferNSAWARE, const double *recvBufferNS)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr     = 0.0;
    RFloat rtErr     = 0.0;
    int    offsetNgb = offsetMPIDataInterpointTurb[ngbID];

    const string &name        = (*nameTurbInterpointBuffer)[nameID];
    int           offsetVar   = offsetInterpointBufferTurb[nameID];
    int           dimData     = dimInterpointBufferTurb[nameID];
    int           offsetAWARE = offsetNgb + offsetNgbLengthBuffer;
    for (int iPoint = 0; iPoint < nIPointOfNeighbor; iPoint++)
    {
        int pointIDAware = iPoint;
        int pointID      = pointIndexForRecv[iPoint];
        for (int iDim = 0; iDim < dimData; iDim++)
        {
            double GPUAware = recvBufferNSAWARE[offsetAWARE + iDim * nIPointOfNeighbor + pointIDAware];
            double GPUOrg   = recvBufferNS[offsetVar + iDim * nIPoint + pointID];
            abErr           = fabs(GPUAware - GPUOrg);
            rtErr           = abErr / Max(fabs(GPUOrg), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in GPUMPIDataDecompressInterpointTurb, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, GPUAware=%.30e, "
                       "GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, GPUAware, GPUOrg);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in GPUMPIDataDecompressInterpointTurb, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, GPUAware=%.30e, "
                       "GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, GPUAware, GPUOrg);
                exit(1);
            }
        }
    }

    //!cout<<"Test receive data "<<name<<" successfully in iZone "<<iZone<<" and jZone "<<jZone<<" globalRank "<<globalRank<<" and localRank "<<localRank<<" !"<<endl;
}

void TestGPUMPIDataDecompressInterpointNS(const int iZone, const int jZone, const int nameID, const int ngbID,
                                          const int nIPoint, const int nIPointOfNeighbor,
                                          const int offsetNgbLengthBuffer, const int *pointIndexForRecv,
                                          const double *recvBufferNSAWARE, const double *recvBufferNS)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr     = 0.0;
    RFloat rtErr     = 0.0;
    int    offsetNgb = offsetMPIDataInterpointNS[ngbID];

    const string &name        = (*nameNSInterpointBuffer)[nameID];
    int           offsetVar   = offsetInterpointBufferNS[nameID];
    int           dimData     = dimInterpointBufferNS[nameID];
    int           offsetAWARE = offsetNgb + offsetNgbLengthBuffer;
    for (int iPoint = 0; iPoint < nIPointOfNeighbor; iPoint++)
    {
        int pointIDAware = iPoint;
        int pointID      = pointIndexForRecv[iPoint];
        for (int iDim = 0; iDim < dimData; iDim++)
        {
            double GPUAware = recvBufferNSAWARE[offsetAWARE + iDim * nIPointOfNeighbor + pointIDAware];
            double GPUOrg   = recvBufferNS[offsetVar + iDim * nIPoint + pointID];
            abErr           = fabs(GPUAware - GPUOrg);
            rtErr           = abErr / Max(fabs(GPUOrg), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in GPUMPIDataDecompressInterpointNS, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, GPUAware=%.30e, "
                       "GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, GPUAware, GPUOrg);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in GPUMPIDataDecompressInterpointNS, iZone = %d, jZone = %d, "
                       "globalRank=%d, localRank=%d, ngbID=%d, nameID=%d, GPUAware=%.30e, "
                       "GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, GPUAware, GPUOrg);
                exit(1);
            }
        }
    }

    //!cout<<"Test receive data "<<name<<" successfully in iZone "<<iZone<<" and jZone "<<jZone<<" globalRank "<<globalRank<<" and localRank "<<localRank<<" !"<<endl;
}

void TestGPUMPIDataDecompressTurb(const int iZone, const int jZone, const int nameID, const int ngbID, const int nIFace,
                                  const int nIFaceOfNeighbor, const int offsetNgbLengthBuffer,
                                  const int *faceIndexForRecv, const double *recvBufferTurbAWARE,
                                  const double *recvBufferTurb)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr     = 0.0;
    RFloat rtErr     = 0.0;
    int    offsetNgb = offsetMPIDataTurb[ngbID];

    const string &name        = (*nameTurbBuffer)[nameID];
    int           offsetVar   = offsetBufferTurb[nameID];
    int           dimData     = dimBufferTurb[nameID];
    int           offsetAWARE = offsetNgb + offsetNgbLengthBuffer;
    for (int iFace = 0; iFace < nIFaceOfNeighbor; iFace++)
    {
        int faceIDAware = iFace;
        int faceID      = faceIndexForRecv[iFace];
        for (int iDim = 0; iDim < dimData; iDim++)
        {
            double GPUAware = recvBufferTurbAWARE[offsetAWARE + iDim * nIFaceOfNeighbor + faceIDAware];
            double GPUOrg   = recvBufferTurb[offsetVar + iDim * nIFace + faceID];
            abErr           = fabs(GPUAware - GPUOrg);
            rtErr           = abErr / Max(fabs(GPUOrg), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in GPUMPIDataDecompressTurb, iZone = %d, jZone = %d, globalRank=%d, "
                       "localRank=%d, ngbID=%d, nameID=%d, GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, GPUAware, GPUOrg);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in GPUMPIDataDecompressTurb, iZone = %d, jZone = %d, globalRank=%d, "
                       "localRank=%d, ngbID=%d, nameID=%d, GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, GPUAware, GPUOrg);
                exit(1);
            }
        }
    }

    //!cout<<"Test receive data "<<name<<" successfully in iZone "<<iZone<<" and jZone "<<jZone<<" globalRank "<<globalRank<<" and localRank "<<localRank<<" !"<<endl;
}

void TestGPUMPIDataDecompressNS(const int iZone, const int jZone, const int nameID, const int ngbID, const int nIFace,
                                const int nIFaceOfNeighbor, const int offsetNgbLengthBuffer,
                                const int *faceIndexForRecv, const double *recvBufferNSAWARE,
                                const double *recvBufferNS)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr     = 0.0;
    RFloat rtErr     = 0.0;
    int    offsetNgb = offsetMPIDataNS[ngbID];

    const string &name        = (*nameNSBuffer)[nameID];
    int           offsetVar   = offsetBufferNS[nameID];
    int           dimData     = dimBufferNS[nameID];
    int           offsetAWARE = offsetNgb + offsetNgbLengthBuffer;
    for (int iFace = 0; iFace < nIFaceOfNeighbor; iFace++)
    {
        int faceIDAware = iFace;
        int faceID      = faceIndexForRecv[iFace];
        for (int iDim = 0; iDim < dimData; iDim++)
        {
            double GPUAware = recvBufferNSAWARE[offsetAWARE + iDim * nIFaceOfNeighbor + faceIDAware];
            double GPUOrg   = recvBufferNS[offsetVar + iDim * nIFace + faceID];
            abErr           = fabs(GPUAware - GPUOrg);
            rtErr           = abErr / Max(fabs(GPUOrg), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in GPUMPIDataDecompressNS, iZone = %d, jZone = %d, globalRank=%d, "
                       "localRank=%d, ngbID=%d, nameID=%d, GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, GPUAware, GPUOrg);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                cout << "Error in " << name << " MPI transfer" << endl;
                printf("Error in GPUMPIDataDecompressNS, iZone = %d, jZone = %d, globalRank=%d, "
                       "localRank=%d, ngbID=%d, nameID=%d, GPUAware=%.30e, GPUOrg=%.30e\n",
                       iZone, jZone, globalRank, localRank, ngbID, nameID, GPUAware, GPUOrg);
                exit(1);
            }
        }
    }

    //!cout<<"Test receive data "<<name<<" successfully in iZone "<<iZone<<" and jZone "<<jZone<<" globalRank "<<globalRank<<" and localRank "<<localRank<<" !"<<endl;
}

void TestAllreduceMinDt(const double h_localMinDt, const double globalMinDt)
{
    double abErr    = 0.0;
    double rtErr    = 0.0;
    abErr           = fabs(globalMinDt - h_localMinDt);
    rtErr           = abErr / Max(fabs(globalMinDt), smallLimit);
    double checknan = abErr;
    if (std::isnan(checknan))
    {
        printf("Error in Allreduce, globalMinDt = %.30e, h_localMinDt = %.30e, not equal\n", globalMinDt, h_localMinDt);
        exit(1);
    }
    if (abErr > abCompare && rtErr > rtCompare)
    {
        printf("Error in Allreduce, globalMinDt = %.30e, h_localMinDt = %.30e, not equal\n", globalMinDt, h_localMinDt);
        exit(1);
    }
    //!cout<<"Test Allreduce MinDt successfully"<<endl;
}

void TestGPUInitSpec(const RFloat *org_spec, const int nTotal)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t  sizeSpec = nTotal * sizeof(RFloat);
    RFloat *gpu_spec = (RFloat *)malloc(sizeSpec);
    HANDLE_API_ERR(cudaMemcpy(gpu_spec, d_spec_ns, sizeSpec, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nTotal; i++)
    {
        abErr           = fabs(gpu_spec[i] - org_spec[i]);
        rtErr           = abErr / Max(fabs(org_spec[i]), smallLimit);
        double checknan = abErr;
        if (std::isnan(checknan))
        {
            printf("Error: In TestGPUInitSpec difference exists in %d term: d_res_ns=%e\t "
                   "org_spec=%e\t abErr= %e\t rtErr= %e\n",
                   i, gpu_spec[i], org_spec[i], abErr, rtErr);
            exit(1);
        }
        if (abErr > abCompare && rtErr > rtCompare)
        {
            printf("Error: In TestGPUInitSpec difference exists in %d term: d_res_ns=%e\t "
                   "org_spec=%e\t abErr= %e\t rtErr= %e\n",
                   i, gpu_spec[i], org_spec[i], abErr, rtErr);
            exit(1);
        }
    }

    free(gpu_spec);
}

void TestGPUDq(const int nTotal, const int nl, RFloat **org_dq)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t  sizeDq = nTotal * nl * sizeof(RFloat);
    RFloat *gpu_dq = (RFloat *)malloc(sizeDq);
    HANDLE_API_ERR(cudaMemcpy(gpu_dq, d_dq_ns, sizeDq, cudaMemcpyDeviceToHost));
    for (int eqnID = 0; eqnID < nl; eqnID++)
    {
        for (int cellID = 0; cellID < nTotal; cellID++)
        {
            abErr           = fabs(gpu_dq[eqnID * nTotal + cellID] - org_dq[eqnID][cellID]);
            rtErr           = abErr / Max(fabs(org_dq[eqnID][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestGPUDq difference exists in eqnID %d cellID %d: d_dq_ns=%e\t "
                       "org_dq=%e\t abErr= %e\t rtErr= %e\n",
                       eqnID, cellID, gpu_dq[eqnID * nTotal + cellID], org_dq[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestGPUDq difference exists in eqnID %d cellID %d: d_dq_ns=%e\t "
                       "org_dq=%e\t abErr= %e\t rtErr= %e\n",
                       eqnID, cellID, gpu_dq[eqnID * nTotal + cellID], org_dq[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }

    free(gpu_dq);
}

void TestGPUTurbDq(const int nTotal, const int nl, RFloat **org_dq)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t  sizeDq = nTotal * nl * sizeof(RFloat);
    RFloat *gpu_dq = (RFloat *)malloc(sizeDq);
    HANDLE_API_ERR(cudaMemcpy(gpu_dq, d_dq_turb, sizeDq, cudaMemcpyDeviceToHost));
    for (int eqnID = 0; eqnID < nl; eqnID++)
    {
        for (int cellID = 0; cellID < nTotal; cellID++)
        {
            abErr           = fabs(gpu_dq[eqnID * nTotal + cellID] - org_dq[eqnID][cellID]);
            rtErr           = abErr / Max(fabs(org_dq[eqnID][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestGPUTurbDq difference exists in eqnID %d cellID %d: "
                       "d_dq_turb=%e\t org_dq=%e\t abErr= %e\t rtErr= %e\n",
                       eqnID, cellID, gpu_dq[eqnID * nTotal + cellID], org_dq[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestGPUTurbDq difference exists in eqnID %d cellID %d: "
                       "d_dq_turb=%e\t org_dq=%e\t abErr= %e\t rtErr= %e\n",
                       eqnID, cellID, gpu_dq[eqnID * nTotal + cellID], org_dq[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }
    free(gpu_dq);
}

void TestSetNSResTmpbyRes(const int nTotal, const int nTotalCell, const int nEquation, RFloat **resTmp)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t sizeRes = nTotal * nEquation * sizeof(RFloat);

    RFloat *gpu_res_tmp = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_res_tmp, d_res_ns_unsteady_tmp, sizeRes, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotalCell; cellID++)
    {
        for (int m = 0; m < nEquation; m++)
        {
            abErr           = fabs(gpu_res_tmp[m * nTotal + cellID] - resTmp[m][cellID]);
            rtErr           = abErr / Max(fabs(resTmp[m][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestSetNSResTmpbyRes difference exists in eqnID %d cellID %d: "
                       "d_res_ns_unsteady_tmp=%.30e\t org_resTmp=%.30e\t abErr= %.30e\t rtErr= "
                       "%.30e\n",
                       m, cellID, gpu_res_tmp[m * nTotal + cellID], resTmp[m][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestSetNSResTmpbyRes difference exists in eqnID %d cellID %d: "
                       "d_res_ns_unsteady_tmp=%.30e\t org_resTmp=%.30e\t abErr= %.30e\t rtErr= "
                       "%.30e\n",
                       m, cellID, gpu_res_tmp[m * nTotal + cellID], resTmp[m][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }
    delete[] gpu_res_tmp;
}

void TestGPUNSDualTimeSourceRes(const int nTotal, const int nTotalCell, const int nEquation, RFloat **res)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t sizeRes = nTotal * nEquation * sizeof(RFloat);

    RFloat *gpu_res = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_res, d_res_ns, sizeRes, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotalCell; cellID++)
    {
        for (int m = 0; m < nEquation; m++)
        {
            abErr           = fabs(gpu_res[m * nTotal + cellID] - res[m][cellID]);
            rtErr           = abErr / Max(fabs(res[m][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestGPUNSDualTimeSourceRes difference exists in eqnID %d cellID "
                       "%d: d_res_ns=%.30e\t org_res=%.30e\t abErr= %.30e\t rtErr= %.30e\n",
                       m, cellID, gpu_res[m * nTotal + cellID], res[cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestGPUNSDualTimeSourceRes difference exists in eqnID %d cellID "
                       "%d: d_res_ns=%.30e\t org_res=%.30e\t abErr= %.30e\t rtErr= %.30e\n",
                       m, cellID, gpu_res[m * nTotal + cellID], res[cellID], abErr, rtErr);
                exit(1);
            }
        }
    }

    delete[] gpu_res;
}

void TestGPUTurbSetResTmpByRes(const int nTotal, const int nTotalCell, const int nEquation, RFloat **resTmp)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t sizeRes = nTotal * nEquation * sizeof(RFloat);

    RFloat *gpu_res_tmp = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_res_tmp, d_res_turb_unsteady_tmp, sizeRes, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotalCell; cellID++)
    {
        for (int m = 0; m < nEquation; m++)
        {
            abErr           = fabs(gpu_res_tmp[m * nTotal + cellID] - resTmp[m][cellID]);
            rtErr           = abErr / Max(fabs(resTmp[m][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestGPUTurbSetResTmpByRes difference exists in eqnID %d cellID "
                       "%d: d_res_ns_unsteady_tmp=%.30e\t org_resTmp=%.30e\t abErr= %.30e\t rtErr= "
                       "%.30e\n",
                       m, cellID, gpu_res_tmp[m * nTotal + cellID], resTmp[m][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestGPUTurbSetResTmpByRes difference exists in eqnID %d cellID "
                       "%d: d_res_ns_unsteady_tmp=%.30e\t org_resTmp=%.30e\t abErr= %.30e\t rtErr= "
                       "%.30e\n",
                       m, cellID, gpu_res_tmp[m * nTotal + cellID], resTmp[m][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }
    delete[] gpu_res_tmp;
}

void TestGPUTurbDualTimeSourceRes(const int nTotal, const int nTotalCell, const int nEquation, RFloat **res)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t sizeRes = nTotal * nEquation * sizeof(RFloat);

    RFloat *gpu_res = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_res, d_res_turb, sizeRes, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotalCell; cellID++)
    {
        for (int m = 0; m < nEquation; m++)
        {
            abErr           = fabs(gpu_res[m * nTotal + cellID] - res[m][cellID]);
            rtErr           = abErr / Max(fabs(res[m][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestGPUTurbDualTimeSourceRes difference exists in eqnID %d cellID "
                       "%d: d_res_turb=%.30e\t org_res=%.30e\t abErr= %.30e\t rtErr= %.30e\n",
                       m, cellID, gpu_res[m * nTotal + cellID], res[cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestGPUTurbDualTimeSourceRes difference exists in eqnID %d cellID "
                       "%d: d_res_turb=%.30e\t org_res=%.30e\t abErr= %.30e\t rtErr= %.30e\n",
                       m, cellID, gpu_res[m * nTotal + cellID], res[cellID], abErr, rtErr);
                exit(1);
            }
        }
    }

    delete[] gpu_res;
}

void TestGPUNSUnsteadyConvergence(const RFloat sum1, const RFloat sum2)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t  sizeSum  = 1 * sizeof(RFloat);
    RFloat *gpu_sum1 = new RFloat[1];
    RFloat *gpu_sum2 = new RFloat[1];
    HANDLE_API_ERR(cudaMemcpy(gpu_sum1, d_sum1_ns_unsteady, sizeSum, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(gpu_sum2, d_sum2_ns_unsteady, sizeSum, cudaMemcpyDeviceToHost));

    abErr           = fabs(*gpu_sum1 - sum1);
    rtErr           = abErr / Max(fabs(sum1), smallLimit);
    double checknan = abErr;
    if (std::isnan(checknan))
    {
        printf("Error: In TestGPUNSUnsteadyConvergence difference exists in "
               "d_sum1_ns_unsteady=%.30e\t org_sum1 = %.30e\t abErr= %.30e\t rtErr= %.30e\n",
               *gpu_sum1, sum1, abErr, rtErr);
        exit(1);
    }

    if (abErr > abCompare && rtErr > rtCompare)
    {
        printf("Error: In TestGPUNSUnsteadyConvergence difference exists in "
               "d_sum1_ns_unsteady=%.30e\t org_sum1 = %.30e\t abErr= %.30e\t rtErr= %.30e\n",
               *gpu_sum1, sum1, abErr, rtErr);
        exit(1);
    }

    abErr    = fabs(*gpu_sum2 - sum2);
    rtErr    = abErr / Max(fabs(sum2), smallLimit);
    checknan = abErr;
    if (std::isnan(checknan))
    {
        printf("Error: In TestGPUNSUnsteadyConvergence difference exists in "
               "d_sum2_ns_unsteady=%.30e\t org_sum2 = %.30e\t abErr= %.30e\t rtErr= %.30e\n",
               *gpu_sum2, sum2, abErr, rtErr);
        exit(1);
    }

    if (abErr > abCompare && rtErr > rtCompare)
    {
        printf("Error: In TestGPUNSUnsteadyConvergence difference exists in "
               "d_sum2_ns_unsteady=%.30e\t org_sum2 = %.30e\t abErr= %.30e\t rtErr= %.30e\n",
               *gpu_sum2, sum2, abErr, rtErr);
        exit(1);
    }

    delete[] gpu_sum1;
    delete[] gpu_sum2;
}

void TestGPUTurbUnsteadyConvergence(const RFloat sum1, const RFloat sum2)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t  sizeSum  = 1 * sizeof(RFloat);
    RFloat *gpu_sum1 = new RFloat[1];
    RFloat *gpu_sum2 = new RFloat[1];
    HANDLE_API_ERR(cudaMemcpy(gpu_sum1, d_sum1_turb_unsteady, sizeSum, cudaMemcpyDeviceToHost));
    HANDLE_API_ERR(cudaMemcpy(gpu_sum2, d_sum2_turb_unsteady, sizeSum, cudaMemcpyDeviceToHost));

    abErr           = fabs(*gpu_sum1 - sum1);
    rtErr           = abErr / Max(fabs(sum1), smallLimit);
    double checknan = abErr;
    if (std::isnan(checknan))
    {
        printf("Error: In TestGPUNSUnsteadyConvergence difference exists in "
               "d_sum1_ns_unsteady=%.30e\t org_sum1 = %.30e\t abErr= %.30e\t rtErr= %.30e\n",
               *gpu_sum1, sum1, abErr, rtErr);
        exit(1);
    }

    if (abErr > abCompare && rtErr > rtCompare)
    {
        printf("Error: In TestGPUNSUnsteadyConvergence difference exists in "
               "d_sum1_ns_unsteady=%.30e\t org_sum1 = %.30e\t abErr= %.30e\t rtErr= %.30e\n",
               *gpu_sum1, sum1, abErr, rtErr);
        exit(1);
    }

    abErr    = fabs(*gpu_sum2 - sum2);
    rtErr    = abErr / Max(fabs(sum2), smallLimit);
    checknan = abErr;
    if (std::isnan(checknan))
    {
        printf("Error: In TestGPUNSUnsteadyConvergence difference exists in "
               "d_sum2_ns_unsteady=%.30e\t org_sum2 = %.30e\t abErr= %.30e\t rtErr= %.30e\n",
               *gpu_sum2, sum2, abErr, rtErr);
        exit(1);
    }

    if (abErr > abCompare && rtErr > rtCompare)
    {
        printf("Error: In TestGPUNSUnsteadyConvergence difference exists in "
               "d_sum2_ns_unsteady=%.30e\t org_sum2 = %.30e\t abErr= %.30e\t rtErr= %.30e\n",
               *gpu_sum2, sum2, abErr, rtErr);
        exit(1);
    }

    delete[] gpu_sum1;
    delete[] gpu_sum2;
}

void TestCallGPUNSUpdateUnsteadyFlow(const int nTotal, const int nEquation, RFloat **qn1, RFloat **qn2, RFloat **resn1,
                                     RFloat **resn2)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t  sizeUnsteady = nTotal * nEquation * sizeof(RFloat);
    RFloat *gpu_qn1      = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_qn1, d_q_ns_unsteady_n1, sizeUnsteady, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotal; cellID++)
    {
        for (int eqnID = 0; eqnID < nEquation; eqnID++)
        {
            abErr           = fabs(gpu_qn1[eqnID * nTotal + cellID] - qn1[eqnID][cellID]);
            rtErr           = abErr / Max(fabs(qn1[eqnID][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_q_ns_unsteady_n1=%.30e, org_qn1 = %.30e, abErr= %.30e, rtErr= "
                       "%.30e\n",
                       cellID, eqnID, gpu_qn1[eqnID * nTotal + cellID], qn1[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_q_ns_unsteady_n1=%.30e, org_qn1 = %.30e, abErr= %.30e, rtErr= "
                       "%.30e\n",
                       cellID, eqnID, gpu_qn1[eqnID * nTotal + cellID], qn1[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }

    delete[] gpu_qn1;

    RFloat *gpu_qn2 = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_qn2, d_q_ns_unsteady_n2, sizeUnsteady, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotal; cellID++)
    {
        for (int eqnID = 0; eqnID < nEquation; eqnID++)
        {
            abErr           = fabs(gpu_qn2[eqnID * nTotal + cellID] - qn2[eqnID][cellID]);
            rtErr           = abErr / Max(fabs(qn2[eqnID][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_q_ns_unsteady_n2=%.30e, org_qn2 = %.30e, abErr= %.30e, rtErr= "
                       "%.30e\n",
                       cellID, eqnID, gpu_qn2[eqnID * nTotal + cellID], qn2[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_q_ns_unsteady_n2=%.30e, org_qn2 = %.30e, abErr= %.30e, rtErr= "
                       "%.30e\n",
                       cellID, eqnID, gpu_qn2[eqnID * nTotal + cellID], qn2[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }

    delete[] gpu_qn2;

    RFloat *gpu_resn1 = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_resn1, d_res_ns_unsteady_n1, sizeUnsteady, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotal; cellID++)
    {
        for (int eqnID = 0; eqnID < nEquation; eqnID++)
        {
            abErr           = fabs(gpu_resn1[eqnID * nTotal + cellID] - resn1[eqnID][cellID]);
            rtErr           = abErr / Max(fabs(resn1[eqnID][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_res_ns_unsteady_n1=%.30e, org_resn1 = %.30e, abErr= %.30e, "
                       "rtErr= %.30e\n",
                       cellID, eqnID, gpu_resn1[eqnID * nTotal + cellID], resn1[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_res_ns_unsteady_n1=%.30e, org_resn1 = %.30e, abErr= %.30e, "
                       "rtErr= %.30e\n",
                       cellID, eqnID, gpu_resn1[eqnID * nTotal + cellID], resn1[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }

    delete[] gpu_resn1;

    RFloat *gpu_resn2 = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_resn2, d_res_ns_unsteady_n2, sizeUnsteady, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotal; cellID++)
    {
        for (int eqnID = 0; eqnID < nEquation; eqnID++)
        {
            abErr           = fabs(gpu_resn2[eqnID * nTotal + cellID] - resn2[eqnID][cellID]);
            rtErr           = abErr / Max(fabs(resn2[eqnID][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_res_ns_unsteady_n2=%.30e, org_resn2 = %.30e, abErr= %.30e, "
                       "rtErr= %.30e\n",
                       cellID, eqnID, gpu_resn2[eqnID * nTotal + cellID], resn2[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_res_ns_unsteady_n2=%.30e, org_resn2 = %.30e, abErr= %.30e, "
                       "rtErr= %.30e\n",
                       cellID, eqnID, gpu_resn2[eqnID * nTotal + cellID], resn2[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }

    delete[] gpu_resn2;
}

void TestCallGPUTurbUpdateUnsteadyFlow(const int nTotal, const int nEquation, RFloat **qn1, RFloat **qn2,
                                       RFloat **resn1, RFloat **resn2)
{
    using namespace GPUMemory;
    using namespace GPUFlowVariables;

    RFloat abErr = 0.0;
    RFloat rtErr = 0.0;

    size_t  sizeUnsteady = nTotal * nEquation * sizeof(RFloat);
    RFloat *gpu_qn1      = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_qn1, d_q_turb_unsteady_n1, sizeUnsteady, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotal; cellID++)
    {
        for (int eqnID = 0; eqnID < nEquation; eqnID++)
        {
            abErr           = fabs(gpu_qn1[eqnID * nTotal + cellID] - qn1[eqnID][cellID]);
            rtErr           = abErr / Max(fabs(qn1[eqnID][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_q_ns_unsteady_n1=%.30e, org_qn1 = %.30e, abErr= %.30e, rtErr= "
                       "%.30e\n",
                       cellID, eqnID, gpu_qn1[eqnID * nTotal + cellID], qn1[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_q_ns_unsteady_n1=%.30e, org_qn1 = %.30e, abErr= %.30e, rtErr= "
                       "%.30e\n",
                       cellID, eqnID, gpu_qn1[eqnID * nTotal + cellID], qn1[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }

    delete[] gpu_qn1;

    RFloat *gpu_qn2 = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_qn2, d_q_turb_unsteady_n2, sizeUnsteady, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotal; cellID++)
    {
        for (int eqnID = 0; eqnID < nEquation; eqnID++)
        {
            abErr           = fabs(gpu_qn2[eqnID * nTotal + cellID] - qn2[eqnID][cellID]);
            rtErr           = abErr / Max(fabs(qn2[eqnID][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_q_ns_unsteady_n2=%.30e, org_qn2 = %.30e, abErr= %.30e, rtErr= "
                       "%.30e\n",
                       cellID, eqnID, gpu_qn2[eqnID * nTotal + cellID], qn2[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_q_ns_unsteady_n2=%.30e, org_qn2 = %.30e, abErr= %.30e, rtErr= "
                       "%.30e\n",
                       cellID, eqnID, gpu_qn2[eqnID * nTotal + cellID], qn2[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }

    delete[] gpu_qn2;

    RFloat *gpu_resn1 = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_resn1, d_res_turb_unsteady_n1, sizeUnsteady, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotal; cellID++)
    {
        for (int eqnID = 0; eqnID < nEquation; eqnID++)
        {
            abErr           = fabs(gpu_resn1[eqnID * nTotal + cellID] - resn1[eqnID][cellID]);
            rtErr           = abErr / Max(fabs(resn1[eqnID][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_res_ns_unsteady_n1=%.30e, org_resn1 = %.30e, abErr= %.30e, "
                       "rtErr= %.30e\n",
                       cellID, eqnID, gpu_resn1[eqnID * nTotal + cellID], resn1[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_res_ns_unsteady_n1=%.30e, org_resn1 = %.30e, abErr= %.30e, "
                       "rtErr= %.30e\n",
                       cellID, eqnID, gpu_resn1[eqnID * nTotal + cellID], resn1[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }

    delete[] gpu_resn1;

    RFloat *gpu_resn2 = new RFloat[nTotal * nEquation];
    HANDLE_API_ERR(cudaMemcpy(gpu_resn2, d_res_turb_unsteady_n2, sizeUnsteady, cudaMemcpyDeviceToHost));
    for (int cellID = 0; cellID < nTotal; cellID++)
    {
        for (int eqnID = 0; eqnID < nEquation; eqnID++)
        {
            abErr           = fabs(gpu_resn2[eqnID * nTotal + cellID] - resn2[eqnID][cellID]);
            rtErr           = abErr / Max(fabs(resn2[eqnID][cellID]), smallLimit);
            double checknan = abErr;
            if (std::isnan(checknan))
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_res_ns_unsteady_n2=%.30e, org_resn2 = %.30e, abErr= %.30e, "
                       "rtErr= %.30e\n",
                       cellID, eqnID, gpu_resn2[eqnID * nTotal + cellID], resn2[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
            if (abErr > abCompare && rtErr > rtCompare)
            {
                printf("Error: In TestCallGPUNSUpdateUnsteadyFlow difference exists in cellID=%d, "
                       "eqnID=%d, d_res_ns_unsteady_n2=%.30e, org_resn2 = %.30e, abErr= %.30e, "
                       "rtErr= %.30e\n",
                       cellID, eqnID, gpu_resn2[eqnID * nTotal + cellID], resn2[eqnID][cellID], abErr, rtErr);
                exit(1);
            }
        }
    }

    delete[] gpu_resn2;
}
