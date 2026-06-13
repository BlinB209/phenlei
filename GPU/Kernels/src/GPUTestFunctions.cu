#include "GPUTestFunctions.h"
#include "Pointer.h"

namespace GPUTestSpace
{
    void TestGPUNodesDataAllocCopy(RDouble *org_x, RDouble *org_y, RDouble *org_z, RDouble *d_x, RDouble *d_y,
                                   RDouble *d_z, int nTotalNode)
    {
        const double flagCompare = 1.0e-40;
        RDouble     *h_x;
        RDouble     *h_y;
        RDouble     *h_z;
#ifdef UNITTESTOUTPUT
        printf("Test GPUNodesDataAllocCopy ...\n");
#endif
        size_t sizeNodes = nTotalNode * sizeof(RDouble);
        h_x              = (RDouble *)malloc(sizeNodes);
        handleAPIErr(cudaMemcpy(h_x, d_x, sizeNodes, cudaMemcpyDeviceToHost), __FILE__, __LINE__);
        h_y = (RDouble *)malloc(sizeNodes);
        handleAPIErr(cudaMemcpy(h_y, d_y, sizeNodes, cudaMemcpyDeviceToHost), __FILE__, __LINE__);
        h_z = (RDouble *)malloc(sizeNodes);
        handleAPIErr(cudaMemcpy(h_z, d_z, sizeNodes, cudaMemcpyDeviceToHost), __FILE__, __LINE__);
        for (int i = 0; i < nTotalNode; i++)
        {
            if ((fabs(h_x[i] - org_x[i]) > flagCompare) || ((fabs(h_y[i] - org_y[i]) > flagCompare))
                || ((fabs(h_z[i] - org_z[i]) > flagCompare)))
            {
                printf("%d term error in GPUNodesDataAllocCopy, org_x=%d,d_x=%d\n", i, org_x[i], h_x[i]);
                exit(EXIT_FAILURE);
            }
        }

#ifdef UNITTESTOUTPUT
        printf("GPUNodesDataAllocCopy is successful!\n");
#endif
        free(h_x);
        free(h_y);
        free(h_z);
    }

    void TestGPUCellCentAllocCopy(RDouble *org_xcc, RDouble *org_ycc, RDouble *org_zcc, RDouble *d_xcc, RDouble *d_ycc,
                                  RDouble *d_zcc, int nTotalCell)
    {
        const double flagCompare = 1.0e-40;
        RDouble     *h_xcc;
        RDouble     *h_ycc;
        RDouble     *h_zcc;
#ifdef UNITTESTOUTPUT
        printf("Test GPUCellCentAllocCopy ...\n");
#endif
        size_t sizeCells = nTotalCell * sizeof(RDouble);
        h_xcc            = (RDouble *)malloc(sizeCells);
        handleAPIErr(cudaMemcpy(h_xcc, d_xcc, sizeCells, cudaMemcpyDeviceToHost), __FILE__, __LINE__);
        h_ycc = (RDouble *)malloc(sizeCells);
        handleAPIErr(cudaMemcpy(h_ycc, d_ycc, sizeCells, cudaMemcpyDeviceToHost), __FILE__, __LINE__);
        h_zcc = (RDouble *)malloc(sizeCells);
        handleAPIErr(cudaMemcpy(h_zcc, d_zcc, sizeCells, cudaMemcpyDeviceToHost), __FILE__, __LINE__);
        for (int i = 0; i < nTotalCell; i++)
        {
            if ((fabs(h_xcc[i] - org_xcc[i]) > flagCompare) || ((fabs(h_ycc[i] - org_ycc[i]) > flagCompare))
                || ((fabs(h_zcc[i] - org_zcc[i]) > flagCompare)))
            {
                printf("%d term error in GPUNodesDataAllocCopy, org_xcc=%d,d_xcc=%d\n", i, org_xcc[i], h_xcc[i]);
                exit(EXIT_FAILURE);
            }
        }
#ifdef UNITTESTOUTPUT
        printf("GPUCellCentAllocCopy is successful!\n");
#endif

        free(h_xcc);
        free(h_ycc);
        free(h_zcc);
    }

    void TestGPUFaceCellRelAllocCopy(const int *org_left_cell_of_face, const int *org_right_cell_of_face,
                                     int *d_left_cell_of_face, int *d_right_cell_of_face, int nTotalFace)
    {
        const double flagCompare = 1.0e-40;
        int         *h_left_cell_of_face;
        int         *h_right_cell_of_face;
        size_t       sizeFaces = nTotalFace * sizeof(int);
        h_left_cell_of_face    = (int *)malloc(sizeFaces);
        handleAPIErr(cudaMemcpy(h_left_cell_of_face, d_left_cell_of_face, sizeFaces, cudaMemcpyDeviceToHost), __FILE__,
                     __LINE__);
        h_right_cell_of_face = (int *)malloc(sizeFaces);
        handleAPIErr(cudaMemcpy(h_right_cell_of_face, d_right_cell_of_face, sizeFaces, cudaMemcpyDeviceToHost),
                     __FILE__, __LINE__);
        for (int i = 0; i < nTotalFace; i++)
        {
            if (abs(h_left_cell_of_face[i] - org_left_cell_of_face[i]) > flagCompare)
            {
                printf("%d term error in GPUFaceCellRelAllocCopy, org_left_cell_of_face=%d, "
                       "d_left_cell_of_face=%d\n",
                       i, org_left_cell_of_face[i], h_left_cell_of_face[i]);
                exit(EXIT_FAILURE);
            }
            if (abs(h_right_cell_of_face[i] - org_right_cell_of_face[i]) > flagCompare)
            {
                printf("%d term error in GPUFaceCellRelAllocCopy, org_right_cell_of_face=%d "
                       ",d_right_cell_of_face=%d\n",
                       i, org_right_cell_of_face[i], h_right_cell_of_face[i]);
                exit(EXIT_FAILURE);
            }
        }

#ifdef UNITTESTOUTPUT
        printf("GPUFaceCellRelAllocCopy is successful!\n");
#endif
        free(h_left_cell_of_face);
        free(h_right_cell_of_face);
    }

    void TestGPUFaceNodeRelAllocCopy(const int *org_face2node, const int *org_node_number_of_each_face,
                                     int *d_face2node, int *d_node_number_of_each_face, int nTotalFace2Node,
                                     int nTotalFace)
    {
        const double flagCompare = 1.0e-40;
        int         *h_face2node;
        int         *h_node_number_of_each_face;
        size_t       sizeFaces     = nTotalFace * sizeof(int);
        h_node_number_of_each_face = (int *)malloc(sizeFaces);
        handleAPIErr(
            cudaMemcpy(h_node_number_of_each_face, d_node_number_of_each_face, sizeFaces, cudaMemcpyDeviceToHost),
            __FILE__, __LINE__);

        for (int i = 0; i < nTotalFace; i++)
        {
            if (abs(h_node_number_of_each_face[i] - org_node_number_of_each_face[i]) > flagCompare)
            {
                printf("%d term error in GPUFaceNodeRelAllocCopy, "
                       "org_node_number_of_each_face=%d, d_node_number_of_each_face=%d\n",
                       i, org_node_number_of_each_face[i], h_node_number_of_each_face[i]);
                exit(EXIT_FAILURE);
            }
        }
        size_t sizeface2node = nTotalFace2Node * sizeof(int);
        h_face2node          = (int *)malloc(sizeface2node);
        handleAPIErr(cudaMemcpy(h_face2node, d_face2node, sizeface2node, cudaMemcpyDeviceToHost), __FILE__, __LINE__);

        for (int i = 0; i < nTotalFace2Node; i++)
        {
            if (abs(h_face2node[i] - org_face2node[i]) > flagCompare)
            {
                printf("%d term error in GPUFaceNodeRelAllocCopy, org_face2node=%d, "
                       "d_face2node=%d\n",
                       i, org_face2node[i], h_face2node[i]);
                exit(EXIT_FAILURE);
            }
        }

#ifdef UNITTESTOUTPUT
        printf("GPUFaceNodeRelAllocCopy is successful!\n");
#endif
        free(h_node_number_of_each_face);
        free(h_face2node);
    }

    void TestGPUQFullCopy(RFloat **q2d, RFloat *d_q_proxy, int nTotal, int nTVar)
    {
#ifdef UNITTESTOUTPUT
        printf("checking d_q_proxy copy ...\n");
#endif
        const double flagCompare = 1.0e-40;
        size_t       sizeQFull   = nTotal * nTVar * sizeof(RFloat);
        RFloat      *h_q_proxy   = (RFloat *)malloc(sizeQFull);
        //! handleAPIErr(cudaMemcpy(h_q_proxy, d_q_proxy, sizeQFull,
        //! cudaMemcpyDeviceToHost), __FILE__, __LINE__);
        HANDLE_API_ERR(cudaMemcpy(h_q_proxy, d_q_proxy, sizeQFull, cudaMemcpyDeviceToHost));
        //!    int nTVar_m = nTVar - 1;
        for (int i = 0; i < nTVar; i++)
        {
            for (int j = 0; j < nTotal; j++)
            {
                if (fabs(q2d[i][j] - h_q_proxy[i * nTotal + j]) > flagCompare)
                {
                    printf("(%d, %d )term error in GPUQFullCopy, h_q_proxy=%d, d_q_proxy=%d\n", i, j, q2d[i][j],
                           h_q_proxy[i * nTVar + j]);
                    exit(EXIT_FAILURE);
                }
            }
        }
        //! printf("checking d_q_proxy copy ...\n");
#ifdef UNITTESTOUTPUT
        printf("GPUQFullCopy is successful!\n");
#endif
        free(h_q_proxy);
    }

    //! Test face normal data allocation on gpu
    void TestGPUFaceNormAllocCopy(RDouble *xfn, RDouble *yfn, RDouble *zfn, RDouble *d_xfn, RDouble *d_yfn,
                                  RDouble *d_zfn, const int n)
    {
        const double flagCompare = 1.0e-40;
        const size_t sizeFaces   = n * sizeof(RDouble);
        double       err         = 0.0;

        RDouble *g_xfn = new RDouble[n];
        RDouble *g_yfn = new RDouble[n];
        RDouble *g_zfn = new RDouble[n];

        HANDLE_API_ERR(cudaMemcpy(g_xfn, d_xfn, sizeFaces, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(g_yfn, d_yfn, sizeFaces, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(g_zfn, d_zfn, sizeFaces, cudaMemcpyDeviceToHost));

#ifdef UNITTESTOUTPUT
        printf("Test GPUFaceNormAllocCopy...\n");
#endif

        for (int i = 0; i < n; i++)
        {
            err = fabs(g_xfn[i] - xfn[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUFaceNormAllocCopy, xfn on host: %f\txfn on "
                       "device:%f\n",
                       i, xfn[i], g_xfn[i]);
                exit(EXIT_FAILURE);
            }
            err = fabs(g_yfn[i] - yfn[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUFaceNormAllocCopy, yfn on host: %f\tyfn on "
                       "device:%f\n",
                       i, yfn[i], g_yfn[i]);
                exit(EXIT_FAILURE);
            }
            err = fabs(g_zfn[i] - zfn[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUFaceNormAllocCopy, zfn on host: %f\tzfn on "
                       "device:%f\n",
                       i, zfn[i], g_zfn[i]);
                exit(EXIT_FAILURE);
            }
        }
#ifdef UNITTESTOUTPUT
        printf("GPUFaceNormAllocCopy is successful!\n");
#endif
        delete[] g_xfn;
        delete[] g_yfn;
        delete[] g_zfn;
    }
    //! Test face center data allocation on gpu
    void TestGPUFaceCentAllocCopy(RDouble *xfc, RDouble *yfc, RDouble *zfc, RDouble *d_xfc, RDouble *d_yfc,
                                  RDouble *d_zfc, const int n)
    {
        const double flagCompare = 1.0e-40;
        const size_t sizeFaces   = n * sizeof(RDouble);
        double       err         = 0.0;

        RDouble *g_xfc = new RDouble[n];
        RDouble *g_yfc = new RDouble[n];
        RDouble *g_zfc = new RDouble[n];

        HANDLE_API_ERR(cudaMemcpy(g_xfc, d_xfc, sizeFaces, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(g_yfc, d_yfc, sizeFaces, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(g_zfc, d_zfc, sizeFaces, cudaMemcpyDeviceToHost));

#ifdef UNITTESTOUTPUT
        printf("Test GPUFaceCentAllocCopy...\n");
#endif
        for (int i = 0; i < n; i++)
        {
            err = fabs(g_xfc[i] - xfc[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUFaceCentAllocCopy, xfc on host: %f\txfc on "
                       "device:%f\n",
                       i, xfc[i], g_xfc[i]);
                exit(EXIT_FAILURE);
            }
            err = fabs(g_yfc[i] - yfc[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUFaceCentAllocCopy, yfc on host: %f\tyfc on "
                       "device:%f\n",
                       i, yfc[i], g_yfc[i]);
                exit(EXIT_FAILURE);
            }
            err = fabs(g_zfc[i] - zfc[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUFaceCentAllocCopy, zfc on host: %f\tzfc on "
                       "device:%f\n",
                       i, zfc[i], g_zfc[i]);
                exit(EXIT_FAILURE);
            }
        }
#ifdef UNITTESTOUTPUT
        printf("GPUFaceCentAllocCopy is successful!\n");
#endif
        delete[] g_xfc;
        delete[] g_yfc;
        delete[] g_zfc;
    }
    void TestGPUAreaVolmAllocCopy(RDouble *area, RDouble *vol, RDouble *d_area, RDouble *d_vol, const int nTotalFace,
                                  const int nTotalCell)
    {
        const double flagCompare = 1.0e-40;
        const size_t sizeFaces   = nTotalFace * sizeof(RDouble);
        const size_t sizeCells   = nTotalCell * sizeof(RDouble);
        double       err         = 0.0;

        RDouble *g_area = new RDouble[nTotalFace];
        RDouble *g_vol  = new RDouble[nTotalCell];

        HANDLE_API_ERR(cudaMemcpy(g_area, d_area, sizeFaces, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(g_vol, d_vol, sizeCells, cudaMemcpyDeviceToHost));

#ifdef UNITTESTOUTPUT
        printf("Test GPUAreaVolmAllocCopy....\n");
#endif
        for (int i = 0; i < nTotalFace; i++)
        {
            err = fabs(g_area[i] - area[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUAreaVolmAllocCopy, area on host: %f\tarea on "
                       "device:%f\n",
                       i, area[i], g_area[i]);
                exit(EXIT_FAILURE);
            }
        }
#ifdef UNITTESTOUTPUT
        printf("GPUAreaAllocCopy is successful!\n");
#endif
        for (int i = 0; i < nTotalCell; i++)
        {
            //! if (i==(nTotalCell-1)) printf("g_vol[%d]=%f, vol[%d]=%f\n", i, g_vol[i],
            //! i, vol[i]);
            err = fabs(g_vol[i] - vol[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUAreaVolmAllocCopy, vol on host: %f\tvol on "
                       "device:%f\n",
                       i, vol[i], g_vol[i]);
                exit(EXIT_FAILURE);
            }
        }
#ifdef UNITTESTOUTPUT
        printf("GPUVolAllocCopy is successful!\n");
#endif
        delete[] g_area;
        delete[] g_vol;
    }

    void TestGPUPointer(RDouble *ptr)
    {
#ifdef UNITTESTOUTPUT
        printf("check if pointer on device is null\n");
#endif
        if (ptr)
            printf("Memory allocated for pointer\n");
        else
            printf("No memory Allocated for pointer\n");
    }

    void TestGPUQCopy(RFloat **q, RFloat *d_q_ns, const int neqn, const int nTotal)
    {
        using namespace PHSPACE;

        const double flagCompare = 1.0e-40;
        RFloat       err         = 0.0;
        size_t       qSize       = neqn * nTotal * sizeof(RFloat);

        //! printf("neqn= %d, nTotal= %d\n", neqn, nTotal);

        RFloat **g_q_ns = NewPointer2<RFloat>(neqn, nTotal);

        HANDLE_API_ERR(cudaMemcpy(g_q_ns[0], d_q_ns, qSize, cudaMemcpyDeviceToHost));

        //! printf("Test GPUQ_NSCopy ...\n");

        //! printQ<<<1,1>>>(d_q_ns);
        //! cudaStreamSynchronize(0);

        for (int i = 0; i < neqn; i++)
        {
            for (int j = 0; j < nTotal; j++)
            {
                err = fabs(g_q_ns[i][j] - q[i][j]);
                if (err > flagCompare)
                {
                    printf("%d , %d term error in GPUQCopy, q on host: %f\t q on device:%f\n", i, j, q[i][j],
                           g_q_ns[i][j]);
                    exit(EXIT_FAILURE);
                }
            }
        }
#ifdef UNITTESTOUTPUT
        printf("GPUQCopy is successful!\n");
#endif
        //! free the temp variabels
        DelPointer2(g_q_ns);
    }
    __global__ void printQ(RFloat *d_q_ns)
    {
#ifdef UNITTESTOUTPUT
        printf("On device: d_q[0][0]= %f\n", d_q_ns[0]);
        printf("On device: d_q[0][221759]= %f\n", d_q_ns[221759]);
        printf("On device: d_q[0][221760]= %f\n", d_q_ns[221760]);
#endif
    }
    //! test the copy of gama
    void TestGPUGamaCopy(RFloat *gama, RFloat *gamal, RFloat *gamar, RFloat *d_gama_ns, RFloat *d_gamaL_ns,
                         RFloat *d_gamaR_ns, const int d_nTotal, const int len_gamaLR)
    {
        using namespace PHSPACE;

        const double flagCompare = 1.0e-40;
        RFloat       err         = 0.0;

        RFloat *g_gama     = new RFloat[d_nTotal];
        RFloat *g_gamaL_ns = new RFloat[d_nTotal];
        RFloat *g_gamaR_ns = new RFloat[d_nTotal];

        size_t gamaSize = d_nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(g_gama, d_gama_ns, gamaSize, cudaMemcpyDeviceToHost));

#ifdef UNITTESTOUTPUT
        printf("Test GPUGamaCopy ...\n");
#endif
        for (int i = 0; i < d_nTotal; i++)
        {
            err = fabs(g_gama[i] - gama[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUGamaCopy, gama on host: %f\t  on device:%f\n", i, gama[i], g_gama[i]);
                exit(EXIT_FAILURE);
            }
        }
#ifdef UNITTESTOUTPUT
        printf("GPUGamaLRCopy is successful!\n");
#endif
        size_t gamaLRSize = len_gamaLR * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(g_gamaL_ns, d_gamaL_ns, gamaLRSize, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(g_gamaR_ns, d_gamaR_ns, gamaLRSize, cudaMemcpyDeviceToHost));

#ifdef UNITTESTOUTPUT
        printf("Test GPUGamaLRCopy ...\n");
#endif
        for (int i = 0; i < len_gamaLR; i++)
        {
            err = fabs(g_gamaL_ns[i] - gamal[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUGamaLRCopy, gamal on host: %f\t on device:%f\n", i, gamal[i],
                       g_gamaL_ns[i]);
                exit(EXIT_FAILURE);
            }

            err = fabs(g_gamaR_ns[i] - gamar[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUGamaLRCopy, gamar on host: %f\t on device:%f\n", i, gamar[i],
                       g_gamaR_ns[i]);
                exit(EXIT_FAILURE);
            }
        }
#ifdef UNITTESTOUTPUT
        printf("GPUGamaLRCopy is successful!\n");
#endif
        delete[] g_gama;
        delete[] g_gamaL_ns;
        delete[] g_gamaR_ns;
    }
    //! test the copy of limiter
    void TestGPULIMITCopy(RFloat **LIMIT, RFloat *d_LIMIT, const int neqn_ns, const int nTotal)
    {
        using namespace PHSPACE;

        RFloat       err         = 0.0;
        const double flagCompare = 1.0e-40;

        RFloat **g_LIMIT = NewPointer2<RFloat>(neqn_ns, nTotal);

        size_t limiterSize = nTotal * neqn_ns * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(g_LIMIT[0], d_LIMIT, limiterSize, cudaMemcpyDeviceToHost));

#ifdef UNITTESTOUTPUT
        printf("Test GPULIMITCopy ...\n");
#endif
        for (int i = 0; i < neqn_ns; i++)
        {
            for (int j = 0; j < nTotal; j++)
            {
                err = fabs(g_LIMIT[i][j] - LIMIT[i][j]);
                if (err > flagCompare)
                {
                    printf("%d ,%d term error in GPULIMITCopy, LIMIT on host: %f\t  on "
                           "device:%f\n",
                           i, j, LIMIT[i][j], g_LIMIT[i][j]);
                    exit(EXIT_FAILURE);
                }
            }
        }
#ifdef UNITTESTOUTPUT
        printf("GPULIMITCopy is successful!\n");
#endif
        DelPointer2(g_LIMIT);
    }
    void TestGPULimitCopy(RFloat *limit, RFloat *d_limit, const int nTotal)
    {
        RFloat       err         = 0.0;
        const double flagCompare = 1.0e-40;

        RFloat *g_limit     = new RFloat[nTotal];
        size_t  limiterSize = nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(g_limit, d_limit, limiterSize, cudaMemcpyDeviceToHost));

#ifdef UNITTESTOUTPUT
        printf("Test GPULimitCopy ...\n");
#endif
        for (int j = 0; j < nTotal; j++)
        {
            err = fabs(g_limit[j] - limit[j]);
            if (err > flagCompare)
            {
                printf("%d term error in GPULimitCopy, limit on host: %f\t  on device:%f\n", j, limit[j], g_limit[j]);
                exit(EXIT_FAILURE);
            }
        }

#ifdef UNITTESTOUTPUT
        printf("GPULimitCopy is successful!\n");
#endif
        delete[] g_limit;
    }
    //! test the copy of dqdx ...
    void TestGPUDqdxDqDyDqDzCopy(RDouble **dqdx, RDouble **dqdy, RDouble **dqdz, RFloat *d_dqdx_proxy,
                                 RFloat *d_dqdy_proxy, RFloat *d_dqdz_proxy, const int neqn_ns, const int nTotal)
    {
        using namespace PHSPACE;

        RFloat       err         = 0.0;
        const double flagCompare = 1.0e-40;

        RFloat **g_dqdx = NewPointer2<RFloat>(neqn_ns, nTotal);
        RFloat **g_dqdy = NewPointer2<RFloat>(neqn_ns, nTotal);
        RFloat **g_dqdz = NewPointer2<RFloat>(neqn_ns, nTotal);

        size_t dSize = neqn_ns * nTotal * sizeof(RFloat);
        HANDLE_API_ERR(cudaMemcpy(g_dqdx[0], d_dqdx_proxy, dSize, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(g_dqdy[0], d_dqdy_proxy, dSize, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(g_dqdz[0], d_dqdz_proxy, dSize, cudaMemcpyDeviceToHost));

#ifdef UNITTESTOUTPUT
        printf("Test GPUDqdxDqDyDqDzCopy ...\n");
#endif
        for (int i = 0; i < neqn_ns; i++)
        {
            for (int j = 0; j < nTotal; j++)
            {
                err = fabs(g_dqdx[i][j] - dqdx[i][j]);
                if (err > flagCompare)
                {
                    printf("%d, %d term error in GPUDqdxCopy, dqdx on host: %f\t  on "
                           "device:%f\n",
                           i, j, dqdx[i][j], g_dqdx[i][j]);
                    exit(EXIT_FAILURE);
                }
            }
        }
        for (int i = 0; i < neqn_ns; i++)
        {
            for (int j = 0; j < nTotal; j++)
            {
                err = fabs(g_dqdy[i][j] - dqdy[i][j]);
                if (err > flagCompare)
                {
                    printf("%d, %d term error in GPUDqdyCopy, dqdy on host: %f\t  on "
                           "device:%f\n",
                           i, j, dqdy[i][j], g_dqdy[i][j]);
                    exit(EXIT_FAILURE);
                }
            }
        }
        for (int i = 0; i < neqn_ns; i++)
        {
            for (int j = 0; j < nTotal; j++)
            {
                err = fabs(g_dqdz[i][j] - dqdz[i][j]);
                if (err > flagCompare)
                {
                    printf("%d, %d term error in GPUDqdzCopy, dqdz on host: %f\t  on "
                           "device:%f\n",
                           i, j, dqdz[i][j], g_dqdz[i][j]);
                    exit(EXIT_FAILURE);
                }
            }
        }

        DelPointer2(g_dqdx);
        DelPointer2(g_dqdy);
        DelPointer2(g_dqdz);

#ifdef UNITTESTOUTPUT
        printf("GPUDqdxDqdxDqdzCopy is successful!\n");
#endif
    }
    //! test the copy of xtn ytn ztn
    void TestGPUXtnYtnZtnCopy(RDouble *xtn, RDouble *ytn, RDouble *ztn, RDouble *d_xtn, RDouble *d_ytn, RDouble *d_ztn,
                              const int nTotalFace)
    {
        double   err         = 0.0;
        double   flagCompare = 1.0e-40;
        size_t   tnSize      = nTotalFace * sizeof(RDouble);
        RDouble *g_xtn       = new RDouble[nTotalFace];
        RDouble *g_ytn       = new RDouble[nTotalFace];
        RDouble *g_ztn       = new RDouble[nTotalFace];

        HANDLE_API_ERR(cudaMemcpy(g_xtn, d_xtn, tnSize, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(g_ytn, d_ytn, tnSize, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(g_ztn, d_ztn, tnSize, cudaMemcpyDeviceToHost));

#ifdef UNITTESTOUTPUT
        printf("Test TestGPUXtnYtnZtnCopy ... \n");
#endif
        for (int i = 0; i < nTotalFace; i++)
        {
            err = fabs(xtn[i] - g_xtn[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUXtnYtnZtnCopy, xtn on host: %d\t  on "
                       "device:%d\n",
                       i, xtn[i], g_xtn[i]);
                exit(1);
            }
            err = fabs(ytn[i] - g_ytn[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUXtnYtnZtnCopy, ytn on host: %d\t  on "
                       "device:%d\n",
                       i, ytn[i], g_ytn[i]);
                exit(1);
            }
            err = fabs(ztn[i] - g_ztn[i]);
            if (err > flagCompare)
            {
                printf("%d term error in GPUXtnYtnZtnCopy, ztn on host: %d\t  on "
                       "device:%d\n",
                       i, ztn[i], g_ztn[i]);
                exit(1);
            }
        }
#ifdef UNITTESTOUTPUT
        printf("GPUXtnYtnZtnCopy is successful!\n");
#endif
        delete[] g_xtn;
        delete[] g_ytn;
        delete[] g_ztn;
    }
    //! test the copy of vgn
    void TestGPUVgnCopy(RDouble *vgn, RDouble *d_vgn, const int nTotalFace)
    {
        double err         = 0.0;
        double flagCompare = 1.0e-40;
        size_t vgnSize     = nTotalFace * sizeof(RDouble);

        RDouble *g_vgn = new RDouble[nTotalFace];

        HANDLE_API_ERR(cudaMemcpy(g_vgn, d_vgn, vgnSize, cudaMemcpyDeviceToHost));

#ifdef UNITTESTOUTPUT
        printf("Test TestGPUVgnCopy ... \n");
#endif
        for (int i = 0; i < nTotalFace; i++)
        {
            err = fabs(vgn[i] - g_vgn[i]);
            if (err > flagCompare)
            {
                printf("%d term error in TestGPUVgnCopy, vgn on host: %d\t  on device:%d\n", i, vgn[i], g_vgn[i]);
                exit(1);
            }
        }
#ifdef UNITTESTOUTPUT
        printf("TestGPUVgnCopy is successful!\n");
#endif
        delete[] g_vgn;
    }
    //! test the copy of res_ns
    void TestGPUResCopy(RFloat **res, RFloat *d_res_ns, const int neqn, const int nTotal)
    {
        double err = 0.0;
        //! why is it too large??
        //! double flagCompare = 1.0e-8;
        double flagCompare = 1.0e-40;
        size_t resSize     = neqn * nTotal * sizeof(RFloat);

        RFloat **g_res_ns = NewPointer2<RFloat>(neqn, nTotal);
        HANDLE_API_ERR(cudaMemcpy(g_res_ns[0], d_res_ns, resSize, cudaMemcpyDeviceToHost));

        for (int i = 0; i < neqn; i++)
        {
            for (int j = 0; j < nTotal; j++)
            {
                err = abs(g_res_ns[i][j] - res[i][j]);
                if (err > flagCompare)
                {
                    printf("%d, %d term error in TestGPUResCopy, res on host: %f\t  on "
                           "device:%f\n",
                           i, j, res[i][j], g_res_ns[i][j]);
                    exit(EXIT_FAILURE);
                }
            }
        }
#ifdef UNITTESTOUTPUT
        printf("TestGPUResCopy is successful!\n");
#endif
        DelPointer2(g_res_ns);
    }
    //! test the copy of epsCell
    void TestGPUepsCellCopy(RFloat *epsCell, RFloat *d_epsCell, const int nTotalCell)
    {
        double err = 0.0;
        //! why is it too large??
        //! double flagCompare = 1.0e-8;
        double flagCompare = 1.0e-40;
        size_t epsSize     = nTotalCell * sizeof(RFloat);

        RFloat *g_epsCell = new RFloat[nTotalCell];
        HANDLE_API_ERR(cudaMemcpy(g_epsCell, d_epsCell, epsSize, cudaMemcpyDeviceToHost));
        for (int i = 0; i < nTotalCell; i++)
        {
            err = abs(g_epsCell[i] - epsCell[i]);
            if (err > flagCompare)
            {
                printf("%d term error in TestGPUepsCellCopy, epsCell on host: %f\t  on "
                       "device:%f\n",
                       i, epsCell[i], g_epsCell[i]);
                exit(EXIT_FAILURE);
            }
        }
#ifdef UNITTESTOUTPUT
        printf("TestGPUepsCellCopy is successful!\n");
#endif
        delete[] g_epsCell;
    }
    void TestGPUDminDmaxAllocCopy(RFloat *dmax, RFloat *dmin, RFloat *d_dmax, RFloat *d_dmin, const int nTotal)
    {
        double err = 0.0;
        //! why is it too large??
        //! double flagCompare = 1.0e-8;
        double flagCompare = 1.0e-40;
        size_t dmSize      = nTotal * sizeof(RFloat);

        RFloat *g_dmax = new RFloat[nTotal];
        RFloat *g_dmin = new RFloat[nTotal];

        HANDLE_API_ERR(cudaMemcpy(g_dmax, d_dmax, dmSize, cudaMemcpyDeviceToHost));
        HANDLE_API_ERR(cudaMemcpy(g_dmin, d_dmin, dmSize, cudaMemcpyDeviceToHost));

        for (int i = 0; i < nTotal; i++)
        {
            err = abs(g_dmax[i] - dmax[i]);
            if (err > flagCompare)
            {
                printf("%d term error in TestGPUDminDmaxAllocCopy for dmax, dmax on host: "
                       "%f\t  on device:%f\n",
                       i, dmax[i], g_dmax[i]);
                exit(EXIT_FAILURE);
            }
            err = abs(g_dmin[i] - dmin[i]);
            if (err > flagCompare)
            {
                printf("%d term error in TestGPUDminDmaxAllocCopy for dmin, dmin on host: "
                       "%f\t  on device:%f\n",
                       i, dmin[i], g_dmin[i]);
                exit(EXIT_FAILURE);
            }
        }
#ifdef UNITTESTOUTPUT
        printf("TestGPUDminDmaxAllocCopy is successful!\n");
#endif
        delete[] g_dmax;
        delete[] g_dmin;
    }

    void TestFaceCellRel(int nTotalCell, int *face_number_of_each_cell, int **cell2face, int *acc_face_number_of_cell)
    {
        int fidh = 0;
        int fidd = 0;
        for (int i = 0; i < nTotalCell; i++)
            for (int j = 0; j < face_number_of_each_cell[i]; j++)
            {
                fidh = cell2face[i][j];
                fidd = *(cell2face[0] + acc_face_number_of_cell[i] + j);
                if (fidh != fidd)
                {
                    cout << "err in transfer on cpu= " << fidh << " on gpu: " << fidd << endl;
                    exit(1);
                }
            }
    }
} //! namespace GPUTestSpace
