#include "GPUCompVisfluxTEST.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
#include "GPUBasicFunctions.h"
#include "Geo_SimpleBC.h"
#include "Gas.h"
#include "Constants.h"
using namespace std;
namespace GPUKernels
{
    void CallGPUComputeVisflux(Grid *gridIn, Param_NSSolverUnstruct *parameters, const int neqn, const int nst,
                               const int ned)
    {
        using namespace GPUMemory;
        using namespace GPUGeomVariables;
        using namespace GPUFlowVariables;
        using namespace IDX;

        UnstructGrid *grid       = UnstructGridCast(gridIn);
        int           nBoundFace = grid->GetNBoundFace();
        int           nTotalCell = grid->GetNTotalCell();

        int numberOfSpecies = parameters->GetNumberOfSpecies();

        int nl    = parameters->GetLaminarNumber();
        int nchem = parameters->GetChemicalFlag();

        int nm = parameters->GetNSEquationNumber();

        RDouble refReNumber = parameters->GetRefReNumber();
        RDouble oreynolds   = parameters->GetoRefReNumber();

        int nolstress = GlobalDataBase::GetIntParaFromDB("nolstress");
        int nrokplus  = GlobalDataBase::GetIntParaFromDB("nrokplus");

        RDouble skew_angle = parameters->GetSkewnessAngle();

        using namespace GAS_SPACE;
        RDouble coeffstateequation = gas->GetCoefficientOfStateEquation();

        const int bctypein        = PHENGLEI::INTERFACE;
        const int bctypesolidface = PHENGLEI::SOLID_SURFACE;

        /* original kernel size
        gridSize= 2;
        blockSize= 40;
        */
        int    gridSize        = 1;
        int    blockSize       = 1;
        int    loopLen         = ned - nst;
        int    regsPerThread   = 154;
        int    residentBlockNo = 2;
        size_t dsMemPerThread = 12 * neqn * sizeof(RFloat); //! calculate by the intended shared memory use in the kernel
        size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUCompVisfluxTEST_S1, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                         dsMemPerBlock, gridSize, blockSize, GPUProp);

#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompVisfluxTEST_S1, dsMemPerBlock, gridSize, blockSize);
#endif
        //!cout<<"gridSize= "<<gridSize<<" blockSize= "<<blockSize<<endl;
        GPUCompVisfluxTEST_S1<<<gridSize, blockSize, dsMemPerBlock>>>(
            nst, ned, nTotalCell, nBoundFace, neqn, nchem, nm, nl, numberOfSpecies, ITT, nrokplus, skew_angle,
            oreynolds, SMALL, TINY, d_left_cell_of_face, d_right_cell_of_face, d_xfn, d_yfn, d_zfn, d_xfc, d_yfc, d_zfc,
            d_xcc, d_ycc, d_zcc, d_area, d_q_ns, d_t_proxy, d_dqdx_proxy, d_dqdy_proxy, d_dqdz_proxy, d_dtdx_proxy,
            d_dtdy_proxy, d_dtdz_proxy, d_prim, d_tm, d_deltl, d_deltr, d_kcp, d_mul, d_mut, d_SEG_LEN, d_boundaryType,
            bctypein, bctypesolidface, coeffstateequation, d_flux);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUCompVisfluxTEST(const int nst, const int ned, const int nTotalCell, const int nBoundFace,
                                       int neqn, const int nchem, const int nm, const int nl, const int numberOfSpecies,
                                       const int ITT, const int nrokplus, const double skew_angle,
                                       const double oreynolds, const double SMALL, const double TINY,
                                       const int *left_cell_of_face, const int *right_cell_of_face, const RDouble *xfn,
                                       const RDouble *yfn, const RDouble *zfn, const RDouble *xfc, const RDouble *yfc,
                                       const RDouble *zfc, const RDouble *xcc, const RDouble *ycc, const RDouble *zcc,
                                       const RDouble *area, const RFloat *qpmv, const RFloat *t, const RFloat *dqdx,
                                       const RFloat *dqdy, const RFloat *dqdz, const RFloat *dtdx, const RFloat *dtdy,
                                       const RFloat *dtdz, const RFloat *prim, const RFloat *tm, const RDouble *deltl,
                                       const RDouble *deltr, const RFloat *kcp, const RFloat *mul, const RFloat *mut,
                                       const int len, const int *boundaryType, const int bctype_interface, RFloat *flux)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        int le, re, j, m, bctype;
        int nTotal = nTotalCell + nBoundFace;
        //!printf("PI=%e\n", PI);
        RFloat  c1, c2, cl, cr;
        RDouble nxs, nys, nzs, t1x, t1y, t1z, t2x, t2y, t2z;
        RFloat  t1, t2, tmid, dtd1, dtd2, dtdn, qnorm;
        RFloat  d1, d2, dtmp, angle1, angle2;
        RDouble dxnl, dynl, dznl, dxnr, dynr, dznr;
        RDouble dxl, dyl, dzl, dxr, dyr, dzr;
        RFloat  dudx, dudy, dudz;
        RFloat  dvdx, dvdy, dvdz;
        RFloat  dwdx, dwdy, dwdz;
        RFloat  vis, divv2p3;
        RFloat  txx, tyy, tzz;
        RFloat  txy, txz, tyz;
        RFloat  um, vm, wm;
        //!RFloat f1, f2;
        size_t  sizeF    = neqn;
        RFloat *f1       = (RFloat *)malloc(sizeF);
        RFloat *f2       = (RFloat *)malloc(sizeF);
        RFloat *fmid     = (RFloat *)malloc(sizeF);
        RFloat *dfdn     = (RFloat *)malloc(sizeF);
        RFloat *dfd1     = (RFloat *)malloc(sizeF);
        RFloat *dfd2     = (RFloat *)malloc(sizeF);
        RFloat *dfdx     = (RFloat *)malloc(sizeF);
        RFloat *dfdy     = (RFloat *)malloc(sizeF);
        RFloat *dfdz     = (RFloat *)malloc(sizeF);
        RFloat *dfdt1    = (RFloat *)malloc(sizeF);
        RFloat *dfdt2    = (RFloat *)malloc(sizeF);
        size_t  sizeFvis = nl;
        RFloat *fvis     = (RFloat *)malloc(sizeFvis);
        for (iface = nst + i; iface < ned; iface += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[iface];
            re = right_cell_of_face[iface];
            j  = iface - nst;

            nxs = xfn[iface];
            nys = yfn[iface];
            nzs = zfn[iface];

            //! Get first tangential vector on the face
            if (GPUABS(nxs) > SMALL)
            {
                t1x = nys;
                t1y = -nxs;
                t1z = 0.0;
            }
            else if (GPUABS(nys) > SMALL)
            {
                t1x = -nys;
                t1y = nxs;
                t1z = 0.0;
            }
            else if (GPUABS(nzs) > SMALL)
            {
                t1x = 0.0;
                t1y = -nzs;
                t1z = nys;
            }
            else
            {
                for (m = 0; m < nl; ++m)
                {
                    //!flux[m][j] = 0.0;
                    flux[m * len + j] = 0.0;
                }
                continue;
            }
            //!normalize the tangential vector
            dtmp = 1.0 / sqrt(t1x * t1x + t1y * t1y + t1z * t1z);
            t1x *= dtmp;
            t1y *= dtmp;
            t1z *= dtmp;

            //!Get second tangential vector by cross dot t1 to normal
            t2x = nys * t1z - nzs * t1y;
            t2y = nzs * t1x - nxs * t1z;
            t2z = nxs * t1y - nys * t1x;

            dxl = xcc[le] - xfc[iface];
            dyl = ycc[le] - yfc[iface];
            dzl = zcc[le] - zfc[iface];

            dxr = xcc[re] - xfc[iface];
            dyr = ycc[re] - yfc[iface];
            dzr = zcc[re] - zfc[iface];

            d1 = nxs * dxl + nys * dyl + nzs * dzl;
            d2 = nxs * dxr + nys * dyr + nzs * dzr;

            dtmp   = -d1 / (sqrt(dxl * dxl + dyl * dyl + dzl * dzl) + SMALL);
            dtmp   = GPUMIN(dtmp, 1.0);
            dtmp   = GPUMAX(dtmp, -1.0);
            angle1 = asin(dtmp) * 180.0 / PI;

            dtmp   = d2 / (sqrt(dxr * dxr + dyr * dyr + dzr * dzr) + SMALL);
            dtmp   = GPUMIN(dtmp, 1.0);
            dtmp   = GPUMAX(dtmp, -1.0);
            angle2 = asin(dtmp) * 180.0 / PI;

            dxnl = nxs * d1 - dxl;
            dynl = nys * d1 - dyl;
            dznl = nzs * d1 - dzl;

            dxnr = nxs * d2 - dxr;
            dynr = nys * d2 - dyr;
            dznr = nzs * d2 - dzr;

            //!d1=(d1x,d1y) = d1 * ( nxs, nys )
            //!(r1x,r1y) - (rLx,rLy) = (dr1x,dr1y) - (drLx,drLy)
            //!(dr1x,dr1y) - (drLx,drLy) = d1 * ( nxs, nys ) - ( dxl, dyl )

            //! quentities at points 1 and 2

            for (m = 0; m < neqn; ++m)
            {
                //!f1[m] = qpmv[m][le];
                //!f2[m] = qpmv[m][re];
                f1[m] = qpmv[m * nTotal + le];
                f2[m] = qpmv[m * nTotal + re];
            }

            //!t1 = t[ITT][le];
            //!t2 = t[ITT][re];
            t1 = t[ITT * nTotal + le];
            t2 = t[ITT * nTotal + re];

            for (m = 0; m < neqn; ++m)
            {
                fmid[m] = 0.5 * (f1[m] + f2[m]);
            }
            tmid = 0.5 * (t1 + t2);

            //! Theroretically, more accurate to include the following terms
            if (angle1 > skew_angle && angle2 > skew_angle)
            {
                for (m = 0; m < neqn; ++m)
                {
                    //!f1[m] += dqdx[m][le] * dxnl + dqdy[m][le] * dynl + dqdz[m][le] * dznl;
                    f1[m] += dqdx[m * nTotal + le] * dxnl + dqdy[m * nTotal + le] * dynl + dqdz[m * nTotal + le] * dznl;
                    //!f2[m] += dqdx[m][re] * dxnr + dqdy[m][re] * dynr + dqdz[m][re] * dznr;
                    f2[m] += dqdx[m * nTotal + re] * dxnr + dqdy[m * nTotal + re] * dynr + dqdz[m * nTotal + re] * dznr;
                }

                //!t1 += dtdx[ITT][le] * dxnl + dtdy[ITT][le] * dynl + dtdz[ITT][le] * dznl;
                t1 += dtdx[ITT * nTotal + le] * dxnl + dtdy[ITT * nTotal + le] * dynl + dtdz[ITT * nTotal + le] * dznl;
                //!t2 += dtdx[ITT][re] * dxnr + dtdy[ITT][re] * dynr + dtdz[ITT][re] * dznr;
                t2 += dtdx[ITT * nTotal + re] * dxnr + dtdy[ITT * nTotal + re] * dynr + dtdz[ITT * nTotal + re] * dznr;

                if (t1 < SMALL) t1 = t[ITT * nTotal + le];
                if (t2 < SMALL) t2 = t[ITT * nTotal + re];

                //! quentities at the face
                for (m = 0; m < neqn; ++m)
                {
                    //!fmid[m] = prim[m][j];
                    //!fmid[m] = prim[m * nTotal + j];
                    fmid[m] = prim[m * len + j];
                }

                //!tmid = tm[ITT][j];
                tmid = tm[ITT * nTotal + j];
            }

            for (m = 0; m < neqn; ++m)
            {
                dfdn[m] = 0.0;
            }
            dtdn = 0.0;

            //!if ( angle1 > 0.0 && angle2 > 0.0 && ABS( d1 ) > TINY && ABS( d2 ) > TINY ){
            if (angle1 > 0.0 && angle2 > 0.0 && fabs(d1) > TINY && fabs(d2) > TINY)
            {
                for (m = 0; m < neqn; ++m)
                {
                    dfd1[m] = (f1[m] - fmid[m]) / d1;
                }

                for (m = 0; m < neqn; ++m)
                {
                    dfd2[m] = (f2[m] - fmid[m]) / d2;
                }

                dtd1 = (t1 - tmid) / d1;
                dtd2 = (t2 - tmid) / d2;

                dtmp = d1 * d1 + d2 * d2;
                c1   = d1 * d1 / dtmp;
                c2   = 1.0 - c1;

                if (iface < nBoundFace)
                {
                    bctype = boundaryType[iface];

                    //!if ( bctype != FANTASY::INTERFACE ){
                    if (bctype != bctype_interface)
                    {
                        c1 = 1.0;
                        c2 = 0.0;
                    }
                }

                for (m = 0; m < neqn; ++m)
                {
                    dfdn[m] = c1 * dfd1[m] + c2 * dfd2[m];
                }

                dtdn = dtd1 * c1 + dtd2 * c2;
            }

            cl = deltl[j];
            cr = deltr[j];

            for (m = 0; m < neqn; ++m)
            {
                //!dfdx[m] = cl * dqdx[m][le] + cr * dqdx[m][re];
                //!dfdy[m] = cl * dqdy[m][le] + cr * dqdy[m][re];
                //!dfdz[m] = cl * dqdz[m][le] + cr * dqdz[m][re];
                dfdx[m] = cl * dqdx[m * nTotal + le] + cr * dqdx[m * nTotal + re];
                dfdy[m] = cl * dqdy[m * nTotal + le] + cr * dqdy[m * nTotal + re];
                dfdz[m] = cl * dqdz[m * nTotal + le] + cr * dqdz[m * nTotal + re];
            }

            //!for ( int m = IU; m <= IW; ++ m ){
            for (m = 1; m <= 3; ++m)
            {
                dfdt1[m] = t1x * dfdx[m] + t1y * dfdy[m] + t1z * dfdz[m];
                dfdt2[m] = t2x * dfdx[m] + t2y * dfdy[m] + t2z * dfdz[m];
            }

            //! now true gradients
            for (m = 1; m <= 3; ++m)
            {
                dfdx[m] = nxs * dfdn[m] + t1x * dfdt1[m] + t2x * dfdt2[m];
                dfdy[m] = nys * dfdn[m] + t1y * dfdt1[m] + t2y * dfdt2[m];
                dfdz[m] = nzs * dfdn[m] + t1z * dfdt1[m] + t2z * dfdt2[m];
            }

            dudx = dfdx[1];
            dudy = dfdy[1];
            dudz = dfdz[1];

            dvdx = dfdx[2];
            dvdy = dfdy[2];
            dvdz = dfdz[2];

            dwdx = dfdx[3];
            dwdy = dfdy[3];
            dwdz = dfdz[3];

            //!qnÖÃÁã
            qnorm = 0.0;

            if (nchem == 1)
            {
                //!In the case, nchem is equal to zero, which is not considered
            }

            qnorm += kcp[j] * dtdn;
            divv2p3 = 2.0 / 3.0 * (dudx + dvdy + dwdz);

            vis = mul[j] + mut[j];
            um  = fmid[1];
            vm  = fmid[2];
            wm  = fmid[3];

            //! stress components
            txx = vis * (2.0 * dudx - divv2p3);
            tyy = vis * (2.0 * dvdy - divv2p3);
            tzz = vis * (2.0 * dwdz - divv2p3);
            txy = vis * (dudy + dvdx);
            txz = vis * (dudz + dwdx);
            tyz = vis * (dvdz + dwdy);

            if (nrokplus > 0)
            {
                //!In this case, nrokplus is always equal to -1, which is not considered
            }

            fvis[0] = 0.0;
            fvis[1] = nxs * txx + nys * txy + nzs * txz;
            fvis[2] = nxs * txy + nys * tyy + nzs * tyz;
            fvis[3] = nxs * txz + nys * tyz + nzs * tzz;
            fvis[4] = um * fvis[1] + vm * fvis[2] + wm * fvis[3] + qnorm;

            for (m = 0; m < nl; ++m)
            {
                //!flux[m][j] = - oreynolds * area[iface] * fvis[m];
                flux[m * len + j] = -oreynolds * area[iface] * fvis[m];
            }

        } //!End of loop
        free(f1);
        free(f2);
        free(fmid);
        free(dfdn);
        free(dfd1);
        free(dfd2);
        free(dfdx);
        free(dfdy);
        free(dfdz);
        free(dfdt1);
        free(dfdt2);
        free(fvis);
    } //!End of function GPUCompVisfluxTEST
    __global__ void GPUCompVisfluxTEST_S1(
        const int nst, const int ned, const int nTotalCell, const int nBoundFace, int neqn, const int nchem,
        const int nm, const int nl, const int numberOfSpecies, const int ITT, const int nrokplus,
        const double skew_angle, const double oreynolds, const double SMALL, const double TINY,
        const int *left_cell_of_face, const int *right_cell_of_face, const RDouble *xfn, const RDouble *yfn,
        const RDouble *zfn, const RDouble *xfc, const RDouble *yfc, const RDouble *zfc, const RDouble *xcc,
        const RDouble *ycc, const RDouble *zcc, const RDouble *area, const RFloat *qpmv, const RFloat *t,
        const RFloat *dqdx, const RFloat *dqdy, const RFloat *dqdz, const RFloat *dtdx, const RFloat *dtdy,
        const RFloat *dtdz, const RFloat *prim, const RFloat *tm, const RDouble *deltl, const RDouble *deltr,
        const RFloat *kcp, const RFloat *mul, const RFloat *mut, const int len, const int *boundaryType,
        const int bctype_interface, const int bctype_solidface, const double coefseq, RFloat *flux)
    {
        const int tidx  = threadIdx.x;
        int       i     = blockDim.x * blockIdx.x + threadIdx.x;
        int       iface = 0;
        int       le, re, j, m, bctype;
        int       nTotal = nTotalCell + nBoundFace;
        //!printf("PI=%e\n", PI);
        RFloat  c1, c2, cl, cr;
        RDouble nxs, nys, nzs, t1x, t1y, t1z, t2x, t2y, t2z;
        RFloat  t1, t2, tmid, dtd1, dtd2, dtdn, qnorm;
        RFloat  d1, d2, dtmp, angle1, angle2;
        RDouble dxnl, dynl, dznl, dxnr, dynr, dznr;
        RDouble dxl, dyl, dzl, dxr, dyr, dzr;
        RFloat  dudx, dudy, dudz;
        RFloat  dvdx, dvdy, dvdz;
        RFloat  dwdx, dwdy, dwdz;
        RFloat  vis, divv2p3;
        RFloat  txx, tyy, tzz;
        RFloat  txy, txz, tyz;
        RFloat  um, vm, wm;
        //!RFloat f1, f2;
        extern __shared__ RFloat array[];
        RFloat                  *f1    = array + blockDim.x * neqn * 0;
        RFloat                  *f2    = array + blockDim.x * neqn * 1;
        RFloat                  *fmid  = array + blockDim.x * neqn * 2;
        RFloat                  *dfdn  = array + blockDim.x * neqn * 3;
        RFloat                  *dfd1  = array + blockDim.x * neqn * 4;
        RFloat                  *dfd2  = array + blockDim.x * neqn * 5;
        RFloat                  *dfdx  = array + blockDim.x * neqn * 6;
        RFloat                  *dfdy  = array + blockDim.x * neqn * 7;
        RFloat                  *dfdz  = array + blockDim.x * neqn * 8;
        RFloat                  *dfdt1 = array + blockDim.x * neqn * 9;
        RFloat                  *dfdt2 = array + blockDim.x * neqn * 10;
        RFloat                  *fvis  = array + blockDim.x * neqn * 11;
        /*
        size_t sizeF = neqn;
        RFloat *f1 = (RFloat*)malloc(sizeF);
        RFloat *f2 = (RFloat*)malloc(sizeF);
        RFloat *fmid = (RFloat*)malloc(sizeF);
        RFloat *dfdn = (RFloat*)malloc(sizeF);
        RFloat *dfd1 = (RFloat*)malloc(sizeF);
        RFloat *dfd2 = (RFloat*)malloc(sizeF);
        RFloat *dfdx = (RFloat*)malloc(sizeF);
        RFloat *dfdy = (RFloat*)malloc(sizeF);
        RFloat *dfdz = (RFloat*)malloc(sizeF);
        RFloat *dfdt1 = (RFloat*)malloc(sizeF);
        RFloat *dfdt2 = (RFloat*)malloc(sizeF);
        size_t sizeFvis = nl;
            RFloat *fvis = (RFloat*)malloc(sizeFvis);
        */
        for (iface = nst + i; iface < ned; iface += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[iface];
            re = right_cell_of_face[iface];
            j  = iface - nst;

            nxs = xfn[iface];
            nys = yfn[iface];
            nzs = zfn[iface];
            //! Get first tangential vector on the face
            if (GPUABS(nxs) > SMALL)
            {
                t1x = nys;
                t1y = -nxs;
                t1z = 0.0;
            }
            else if (GPUABS(nys) > SMALL)
            {
                t1x = -nys;
                t1y = nxs;
                t1z = 0.0;
            }
            else if (GPUABS(nzs) > SMALL)
            {
                t1x = 0.0;
                t1y = -nzs;
                t1z = nys;
            }
            else
            {
                for (m = 0; m < nl; ++m)
                {
                    //!flux[m][j] = 0.0;
                    flux[m * len + j] = 0.0;
                }
                continue;
            }
            //!normalize the tangential vector
            dtmp = 1.0 / sqrt(t1x * t1x + t1y * t1y + t1z * t1z);
            t1x *= dtmp;
            t1y *= dtmp;
            t1z *= dtmp;

            //!Get second tangential vector by cross dot t1 to normal
            t2x = nys * t1z - nzs * t1y;
            t2y = nzs * t1x - nxs * t1z;
            t2z = nxs * t1y - nys * t1x;

            dxl = xcc[le] - xfc[iface];
            dyl = ycc[le] - yfc[iface];
            dzl = zcc[le] - zfc[iface];

            dxr = xcc[re] - xfc[iface];
            dyr = ycc[re] - yfc[iface];
            dzr = zcc[re] - zfc[iface];

            d1 = nxs * dxl + nys * dyl + nzs * dzl;
            d2 = nxs * dxr + nys * dyr + nzs * dzr;

            dtmp   = -d1 / (sqrt(dxl * dxl + dyl * dyl + dzl * dzl) + SMALL);
            dtmp   = GPUMIN(dtmp, 1.0);
            dtmp   = GPUMAX(dtmp, -1.0);
            angle1 = asin(dtmp) * 180.0 / PI;

            dtmp   = d2 / (sqrt(dxr * dxr + dyr * dyr + dzr * dzr) + SMALL);
            dtmp   = GPUMIN(dtmp, 1.0);
            dtmp   = GPUMAX(dtmp, -1.0);
            angle2 = asin(dtmp) * 180.0 / PI;

            dxnl = nxs * d1 - dxl;
            dynl = nys * d1 - dyl;
            dznl = nzs * d1 - dzl;

            dxnr = nxs * d2 - dxr;
            dynr = nys * d2 - dyr;
            dznr = nzs * d2 - dzr;

            //!d1=(d1x,d1y) = d1 * ( nxs, nys )
            //!(r1x,r1y) - (rLx,rLy) = (dr1x,dr1y) - (drLx,drLy)
            //!(dr1x,dr1y) - (drLx,drLy) = d1 * ( nxs, nys ) - ( dxl, dyl )

            //! quentities at points 1 and 2

            for (m = 0; m < neqn; ++m)
            {
                //!f1[m*blockDim.x+tidx] = qpmv[m][le];
                //!f2[m*blockDim.x+tidx] = qpmv[m][re];
                f1[m * blockDim.x + tidx] = qpmv[m * nTotal + le];
                f2[m * blockDim.x + tidx] = qpmv[m * nTotal + re];
                //!f1[m*blockDim.x+tidx] = prim[m * nTotal + le];
                //!f2[m*blockDim.x+tidx] = prim[m * nTotal + re];
            }

            //!t1 = t[ITT][le];
            //!t2 = t[ITT][re];
            t1 = t[ITT * nTotal + le];
            t2 = t[ITT * nTotal + re];

            for (m = 0; m < neqn; ++m)
            {
                fmid[m * blockDim.x + tidx] = 0.5 * (f1[m * blockDim.x + tidx] + f2[m * blockDim.x + tidx]);
            }
            tmid = 0.5 * (t1 + t2);

            if (iface < nBoundFace)
            {
                if (boundaryType[iface] == bctype_solidface)
                {
                    //!tMid = tm[ITT][jFace];
                    tmid                        = tm[ITT * nTotal + j];
                    fmid[0 * blockDim.x + tidx] = fmid[4 * blockDim.x + tidx] / (coefseq * tmid);
                }
            }

            //! Theroretically, more accurate to include the following terms
            if (angle1 > skew_angle && angle2 > skew_angle)
            {
                for (m = 0; m < neqn; ++m)
                {
                    //!f1[m*blockDim.x+tidx] += dqdx[m][le] * dxnl + dqdy[m][le] * dynl + dqdz[m][le] * dznl;
                    f1[m * blockDim.x + tidx] +=
                        dqdx[m * nTotal + le] * dxnl + dqdy[m * nTotal + le] * dynl + dqdz[m * nTotal + le] * dznl;
                    //!f2[m*blockDim.x+tidx] += dqdx[m][re] * dxnr + dqdy[m][re] * dynr + dqdz[m][re] * dznr;
                    f2[m * blockDim.x + tidx] +=
                        dqdx[m * nTotal + re] * dxnr + dqdy[m * nTotal + re] * dynr + dqdz[m * nTotal + re] * dznr;
                }

                //!t1 += dtdx[ITT][le] * dxnl + dtdy[ITT][le] * dynl + dtdz[ITT][le] * dznl;
                t1 += dtdx[ITT * nTotal + le] * dxnl + dtdy[ITT * nTotal + le] * dynl + dtdz[ITT * nTotal + le] * dznl;
                //!t2 += dtdx[ITT][re] * dxnr + dtdy[ITT][re] * dynr + dtdz[ITT][re] * dznr;
                t2 += dtdx[ITT * nTotal + re] * dxnr + dtdy[ITT * nTotal + re] * dynr + dtdz[ITT * nTotal + re] * dznr;

                if (t1 < SMALL) t1 = t[ITT * nTotal + le];
                if (t2 < SMALL) t2 = t[ITT * nTotal + re];

                //! quentities at the face
                for (m = 0; m < neqn; ++m)
                {
                    //!fmid[m*blockDim.x+tidx] = prim[m][j];
                    //!fmid[m*blockDim.x+tidx] = prim[m * nTotal + j];
                    fmid[m * blockDim.x + tidx] = prim[m * len + j];
                }

                //!tmid = tm[ITT][j];
                tmid = tm[ITT * nTotal + j];
            }

            for (m = 0; m < neqn; ++m)
            {
                dfdn[m * blockDim.x + tidx] = 0.0;
            }
            dtdn = 0.0;

            //!if ( angle1 > 0.0 && angle2 > 0.0 && ABS( d1 ) > TINY && ABS( d2 ) > TINY ){
            if (angle1 > 0.0 && angle2 > 0.0 && fabs(d1) > TINY && fabs(d2) > TINY)
            {
                for (m = 0; m < neqn; ++m)
                {
                    dfd1[m * blockDim.x + tidx] = (f1[m * blockDim.x + tidx] - fmid[m * blockDim.x + tidx]) / d1;
                }

                for (m = 0; m < neqn; ++m)
                {
                    dfd2[m * blockDim.x + tidx] = (f2[m * blockDim.x + tidx] - fmid[m * blockDim.x + tidx]) / d2;
                }

                dtd1 = (t1 - tmid) / d1;
                dtd2 = (t2 - tmid) / d2;

                dtmp = d1 * d1 + d2 * d2;
                c1   = d1 * d1 / dtmp;
                c2   = 1.0 - c1;

                if (iface < nBoundFace)
                {
                    bctype = boundaryType[iface];

                    //!if ( bctype != FANTASY::INTERFACE ){
                    if (bctype != bctype_interface)
                    {
                        c1 = 1.0;
                        c2 = 0.0;
                    }
                }

                for (m = 0; m < neqn; ++m)
                {
                    dfdn[m * blockDim.x + tidx] = c1 * dfd1[m * blockDim.x + tidx] + c2 * dfd2[m * blockDim.x + tidx];
                }

                dtdn = dtd1 * c1 + dtd2 * c2;
            }

            cl = deltl[j];
            cr = deltr[j];

            for (m = 0; m < neqn; ++m)
            {
                //!dfdx[m*blockDim.x+tidx] = cl * dqdx[m][le] + cr * dqdx[m][re];
                //!dfdy[m*blockDim.x+tidx] = cl * dqdy[m][le] + cr * dqdy[m][re];
                //!dfdz[m*blockDim.x+tidx] = cl * dqdz[m][le] + cr * dqdz[m][re];
                dfdx[m * blockDim.x + tidx] = cl * dqdx[m * nTotal + le] + cr * dqdx[m * nTotal + re];
                dfdy[m * blockDim.x + tidx] = cl * dqdy[m * nTotal + le] + cr * dqdy[m * nTotal + re];
                dfdz[m * blockDim.x + tidx] = cl * dqdz[m * nTotal + le] + cr * dqdz[m * nTotal + re];
            }

            //!for ( int m = IU; m <= IW; ++ m ){
            for (m = 1; m <= 3; ++m)
            {
                dfdt1[m * blockDim.x + tidx] = t1x * dfdx[m * blockDim.x + tidx] + t1y * dfdy[m * blockDim.x + tidx]
                                               + t1z * dfdz[m * blockDim.x + tidx];
                dfdt2[m * blockDim.x + tidx] = t2x * dfdx[m * blockDim.x + tidx] + t2y * dfdy[m * blockDim.x + tidx]
                                               + t2z * dfdz[m * blockDim.x + tidx];
            }

            //! now true gradients
            for (m = 1; m <= 3; ++m)
            {
                dfdx[m * blockDim.x + tidx] = nxs * dfdn[m * blockDim.x + tidx] + t1x * dfdt1[m * blockDim.x + tidx]
                                              + t2x * dfdt2[m * blockDim.x + tidx];
                dfdy[m * blockDim.x + tidx] = nys * dfdn[m * blockDim.x + tidx] + t1y * dfdt1[m * blockDim.x + tidx]
                                              + t2y * dfdt2[m * blockDim.x + tidx];
                dfdz[m * blockDim.x + tidx] = nzs * dfdn[m * blockDim.x + tidx] + t1z * dfdt1[m * blockDim.x + tidx]
                                              + t2z * dfdt2[m * blockDim.x + tidx];
            }

            dudx = dfdx[1 * blockDim.x + tidx];
            dudy = dfdy[1 * blockDim.x + tidx];
            dudz = dfdz[1 * blockDim.x + tidx];

            dvdx = dfdx[2 * blockDim.x + tidx];
            dvdy = dfdy[2 * blockDim.x + tidx];
            dvdz = dfdz[2 * blockDim.x + tidx];

            dwdx = dfdx[3 * blockDim.x + tidx];
            dwdy = dfdy[3 * blockDim.x + tidx];
            dwdz = dfdz[3 * blockDim.x + tidx];

            //!qnÖÃÁã
            qnorm = 0.0;

            if (nchem == 1)
            {
                //!In the case, nchem is equal to zero, which is not considered
            }

            qnorm += kcp[j] * dtdn;
            divv2p3 = 2.0 / 3.0 * (dudx + dvdy + dwdz);

            vis = mul[j] + mut[j];
            um  = fmid[1 * blockDim.x + tidx];
            vm  = fmid[2 * blockDim.x + tidx];
            wm  = fmid[3 * blockDim.x + tidx];

            //! stress components
            txx = vis * (2.0 * dudx - divv2p3);
            tyy = vis * (2.0 * dvdy - divv2p3);
            tzz = vis * (2.0 * dwdz - divv2p3);
            txy = vis * (dudy + dvdx);
            txz = vis * (dudz + dwdx);
            tyz = vis * (dvdz + dwdy);

            if (nrokplus > 0)
            {
                //!In this case, nrokplus is always equal to -1, which is not considered
            }

            fvis[0 * blockDim.x + tidx] = 0.0;
            fvis[1 * blockDim.x + tidx] = nxs * txx + nys * txy + nzs * txz;
            fvis[2 * blockDim.x + tidx] = nxs * txy + nys * tyy + nzs * tyz;
            fvis[3 * blockDim.x + tidx] = nxs * txz + nys * tyz + nzs * tzz;
            fvis[4 * blockDim.x + tidx] = um * fvis[1 * blockDim.x + tidx] + vm * fvis[2 * blockDim.x + tidx]
                                          + wm * fvis[3 * blockDim.x + tidx] + qnorm;

            for (m = 0; m < nl; ++m)
            {
                //!flux[m][j] = - oreynolds * area[iface] * fvis[m*blockDim.x+tidx];
                flux[m * len + j] = -oreynolds * area[iface] * fvis[m * blockDim.x + tidx];
            }

        } //!End of loop
        /*
        free(f1);
        free(f2);
        free(fmid);
        free(dfdn);
        free(dfd1);
        free(dfd2);
        free(dfdx);
        free(dfdy);
        free(dfdz);
        free(dfdt1);
        free(dfdt2);
        free(fvis);
        */
    } //!End of function GPUCompVisfluxTEST

} //! namespace GPUKernels
