#include <iostream>
#include "GPUBasicFunctions.h"
#include "GPUNSSolver.h"
#include "BasicDeviceVariables.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
#include "Constants.h"
#include "Geo_SimpleBC.h"
#include "Gas.h"

using namespace std;
namespace GPUNSSolverUnstruct
{
    void CallGPUBoundary(Grid *gridIn, Param_NSSolverUnstruct *parameters)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        using namespace GPUControlVariables;
        using namespace PHENGLEI;

        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nm         = GlobalDataBase::GetIntParaFromDB("nm");
        int nl         = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem      = GlobalDataBase::GetIntParaFromDB("nchem");
        int nBoundFace = grid->GetNBoundFace();

        double  refGama                   = parameters->GetRefGama();
        RDouble refMachNumber             = parameters->GetRefMachNumber();
        RDouble twall                     = parameters->GetWallTemperature();
        double  refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
        double  tw                        = twall / refDimensionalTemperature;
        int     ifLowSpeedPrecon          = parameters->GetIfLowSpeedPrecon();
        bool    iviscous                  = parameters->IsViscous();

        int nTotal = d_nTotal;

        int gridSize  = 1;
        int blockSize = 1;

        int    loopLen         = nBoundFace;
        int    regsPerThread   = 106;
        int    residentBlockNo = 4;
        size_t dsMemPerThread =
            2 * (nl + nchem) * sizeof(RFloat); //! calculate by the intended shared memory use in the kernel
        size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUBoundary, loopLen, regsPerThread, residentBlockNo, dsMemPerThread, dsMemPerBlock,
                         gridSize, blockSize, GPUProp);

#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUBoundary, dsMemPerBlock, gridSize, blockSize);
#endif
        //! cannot get the positive acc
        GPUBoundary<<<gridSize, blockSize, dsMemPerBlock>>>(
            d_q_ns, d_prim_inf, d_boundaryType, d_left_cell_of_face, d_right_cell_of_face, d_xfn, d_yfn, d_zfn, d_xtn,
            d_ytn, d_ztn, d_vgn, d_gama_ns, refMachNumber, refGama, tw, nl, nchem, nBoundFace, nTotal, iviscous, twall,
            nm, PHENGLEI::EXTRAPOLATION, PHENGLEI::SYMMETRY, PHENGLEI::FARFIELD, PHENGLEI::OUTFLOW, ifLowSpeedPrecon,
            PHENGLEI::SOLID_SURFACE, PHENGLEI::INFLOW);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    __global__ void GPUBoundary(RFloat *q, RFloat *prim_inf, const int *d_boundaryType, const int *left_cell_of_face,
                                const int *right_cell_of_face, RDouble *xfn, RDouble *yfn, RDouble *zfn, RDouble *xtn,
                                RDouble *ytn, RDouble *ztn, RDouble *vgn, RFloat *gama, RDouble refMachNumber,
                                double refGama, double tw, const int nl, const int nchem, const int nBoundFace,
                                const int nTotal, const int iviscous, const RFloat twall, const int nm,
                                const int EXTRAPOLATION, const int SYMMETRY, const int FARFIELD, const int OUTFLOW,
                                const int ifLowSpeedPrecon, const int SOLID_SURFACE, const int INFLOW)
    {
        extern __shared__ RFloat array[];
        RFloat                  *prims = array + 0 * (nl + nchem) * blockDim.x;
        RFloat                  *primt = array + 1 * (nl + nchem) * blockDim.x;
        /*    
        RFloat *prims, *primt;
            prims = new RFloat[nl+nchem];
           primt = new RFloat[nl+nchem];
        */
        //! Boundary condition.
        int       count = 0;
        int       neqn  = nl + nchem;
        const int tidx  = threadIdx.x;
        const int bidx  = blockIdx.x;
        for (int iface = bidx * blockDim.x + tidx; iface < nBoundFace; iface += gridDim.x * blockDim.x)
        {
            int le     = left_cell_of_face[iface];
            int re     = right_cell_of_face[iface];
            int bctype = d_boundaryType[iface];
            for (int m = 0; m < neqn; ++m)
            {
                prims[neqn * tidx + m] = q
                    [m * nTotal
                     + le]; //! it's not a good idea to access sharedMemory by neqn* tidx+ m, m*blockDim.x+tidx recommanded because of access continuity for threads, but it's convenient
            }

            if (bctype < 0)
            {
                continue;
            }
            else if (bctype == EXTRAPOLATION)
            {
                Outflow_BC(prims + neqn * tidx, primt + neqn * tidx, nl, nchem);
            }
            //!else if ( IsWall( bctype ) )
            else if (bctype == SOLID_SURFACE)
            {
                if (iviscous == 0)
                {
                    Symmetry_BC(prims + neqn * tidx, primt + neqn * tidx, xfn[iface], yfn[iface], zfn[iface],
                                vgn[iface], nl, nchem);
                }
                else if (twall <= 0.0)
                {
                    //!viscous adiabatic wall
                    Vis_Adi_Wall_BC(prims + neqn * tidx, primt + neqn * tidx, xfn[iface], yfn[iface], zfn[iface],
                                    xtn[iface], ytn[iface], ztn[iface], refMachNumber, nl, nchem);
                }
                else
                {
                    //! not used, comment by suxnu20191210
                    //!viscous iso-thermal wall
                    //!Vis_Iso_Wall_BC(prims, primt, xfn[iface], yfn[iface], zfn[iface], xtn[iface], ytn[iface], ztn[iface], refMachNumber, tw, nl, nchem);
                }
            }
            else if (bctype == SYMMETRY)
            {
                Symmetry_BC(prims + neqn * tidx, primt + neqn * tidx, xfn[iface], yfn[iface], zfn[iface], vgn[iface],
                            nl, nchem);
            }
            else if (bctype == FARFIELD)
            {
                if (ifLowSpeedPrecon == 0)
                {
                    Farfield_BC(prims + neqn * tidx, prim_inf, primt + neqn * tidx, refGama, gama[le], xfn[iface],
                                yfn[iface], zfn[iface], vgn[iface], nl, nchem);
                    //!Farfield_BC_checkMG(prims+neqn* tidx, prim_inf, primt+neqn* tidx, refGama, gama[le], xfn[iface], yfn[iface], zfn[iface], vgn[iface], nl, nchem, iface);
                    //!    Farfield_BC_check(prims+neqn* tidx, prim_inf, primt+neqn* tidx, refGama, gama[le], xfn[iface], yfn[iface], zfn[iface], vgn[iface], nl, nchem, iface);
                }
                else
                {
                    //! not used, comment by sunxu20191210
                    //!Precon_Farfield_BC(prims, prim_inf, primt, refGama, gama[le], xfn[iface], yfn[iface], zfn[iface], vgn[iface],refMachNumber, nl, nchem);
                }
            }

            else if (bctype == INFLOW)
            {
                for (int m = 0; m < nl + nchem; ++m)
                {
                    primt[neqn * tidx + m] = q[m * nTotal + re];
                }

                //!Isentropic_Inflow_BC(prims, primt, refGama, nl, nchem);
                //!Inflow_JET_BC(prim_inf, primt, refGama, nl, nchem);

                Inflow_BC(prim_inf, primt + neqn * tidx, nl, nchem);

                //!RFloat pst = P0Wall(poo, refGama, reference_mach_number);
                //!RFloat rst = R0Wall(roo, refGama, reference_mach_number);
                //!EngineInflow_BC(prims, primt, pst, rst, refGama, nl, nchem);
            }

            else if (bctype == OUTFLOW)
            {
                Outflow_BC(prims + neqn * tidx, primt + neqn * tidx, nl, nchem);
                //!MassOutflow_BC(prims, primt, mass_alpha, nl, nchem);
                //!RFloat pres = 0.6 * poo;
                //!Pressure_OutletBC(prims, primt, pres, nl, nchem);
            }

            //! for ( int m = 0; m < nl + nchem; ++ m )
            //! {
            //!     q[m*nTotal+re] = primt[neqn* tidx+m];
            //! }

            for (int m = 0; m < nm; ++m)
            {
                q[m * nTotal + re] = primt[m + neqn * tidx];
            }
            /* this part is not used in the m6 case, so optimized later, by sunxu 20191210
            else if ( bctype == PRESSURE_OUTLET)
            {
                Pressure_OutletBC(prims, primt, totalP_outlet, nl, nchem);
            }
            else if ( bctype == PRESSURE_INLET )
            {
                if ( directionMethod == 2 )
                {
                    //! vertical to the face normal.
                    direction_inlet[0] = - xfn[iface];
                    direction_inlet[1] = - yfn[iface];
                    direction_inlet[2] = - zfn[iface];
                }
                //!Pressure_InletBC(prims, primt, refGama, reference_mach_number, totalP, totalT, direction, nl, nchem);
                Pressure_InletBC_Riemann(prims, primt, refGama, refMachNumber, xfn[iface], yfn[iface], zfn[iface],
                                         totalP_inlet, totalT_inlet, direction_inlet, nl, nchem);
            }
            else if ( bctype == PRESSURE_OUTLET_PLATE)
            {
                Pressure_OutletBC_Plate(prims, primt, totalP_outlet, nl, nchem);
            }
            else if ( bctype == PRESSURE_INLET_PLATE )
            {
                Pressure_InletBC_Plate(prims, primt, refGama, refMachNumber, totalP_inlet, totalT_inlet, direction_inlet, nl, nchem);            
            }
            else if ( SetSubsonicBC( grid, bcr, bctype, prims, primt, -xfn[iface], -yfn[iface], -zfn[iface],
                              t[ITT][le] ) )
            {
            //!             cout << "Error : Illegal BCtype ID " << bctype << endl;
            //!             exit(0);
            }
            else if ( bctype == OVERSET )
            {
                Farfield_BC(prims, prim_inf, primt, refGama, gama[le], xfn[iface], yfn[iface], zfn[iface], vgn[iface], nl, nchem);
            }
            else
            {
                TK_Exit::ExceptionExit("Error: this boundary type does not exist!\n");
            }
            */
            /* not functional part, so comment it. by sunxu20191210
            if ( ( primt[IP] <= 0.0 || primt[IR] <= 0.0 ) && count < 1 )
            {
                    RFloat * xfc = grid->GetFaceCenterX();
                    RFloat * yfc = grid->GetFaceCenterY();
                    RFloat * zfc = grid->GetFaceCenterZ();
                cout << " zone = " << grid->GetZoneID() << " bctype = " << bctype << " primt[IP] = " << primt[IP] <<  "primt[IR] = " << primt[IR] << "\n";
                cout << " at surface position of: " << xfc[iface] << ", " << yfc[iface] << ", " << zfc[iface] << endl;
                    cout << " iface = " << iface << ", left cell ID = " << le << ", nBoundFace = " << nBoundFace << ", on level " << grid->GetLevel() << "\n";
                count ++;
            }
            */
        }
        /*
        delete []prims;
        delete []primt;
        */
    }

    __device__ void Outflow_BC(RFloat *prims, RFloat *primt, int nl, int nchem)
    {
        int neqn = nl + nchem;
        for (int m = 0; m < neqn; ++m)
        {
            primt[m] = prims[m];
        }
    }
    __device__ void Inflow_BC(RFloat *prims, RFloat *primt, int nl, int nchem)
    {
        int neqn = nl + nchem;
        for (int m = 0; m < neqn; ++m)
        {
            primt[m] = prims[m];
        }
    }

    __device__ void Symmetry_BC(RFloat *prims, RFloat *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn,
                                int nl, int nchem)
    {
        const int IR = 0;
        const int IU = 1;
        const int IV = 2;
        const int IW = 3;
        const int IP = 4;

        int neqn = nl + nchem;
        for (int m = 0; m < neqn; ++m)
        {
            primt[m] = prims[m];
        }

        double um, vm, wm, vn;
        um = prims[IU];
        vm = prims[IV];
        wm = prims[IW];

        vn = nxs * um + nys * vm + nzs * wm - vgn;

        primt[IU] = um - two * nxs * vn;
        primt[IV] = vm - two * nys * vn;
        primt[IW] = wm - two * nzs * vn;
    }
    //!viscous adiabatic wall
    __device__ void Vis_Adi_Wall_BC(RFloat *prims, RFloat *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble xtn,
                                    RDouble ytn, RDouble ztn, RFloat reference_mach_number, int nl, int nchem)
    {
        const int IR = 0;
        const int IU = 1;
        const int IV = 2;
        const int IW = 3;
        const int IP = 4;

        int neqn = nl + nchem;
        for (int m = 0; m < neqn; ++m)
        {
            primt[m] = prims[m];
        }

        RDouble uwall = 0.0;
        RDouble vwall = 0.0;
        RDouble wwall = 0.0;

        //! primt[IU] = - prims[IU] + two * xtn;
        //! primt[IV] = - prims[IV] + two * ytn;
        //! primt[IW] = - prims[IW] + two * ztn;

        primt[IU] = -prims[IU] + two * uwall;
        primt[IV] = -prims[IV] + two * vwall;
        primt[IW] = -prims[IW] + two * wwall;
    }

    //!3D
    __device__ void Farfield_BC(RFloat *prims, RFloat *prim_inf, RFloat *primt, RFloat gama0, RFloat gama, RDouble nxs,
                                RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem)
    {
        const int IR = 0;
        const int IU = 1;
        const int IV = 2;
        const int IW = 3;
        const int IP = 4;

        int neqn = nl + nchem;

        RFloat rin, uin, vin, win, pin, cin, vni;
        RFloat roo, uoo, voo, woo, poo, coo, vno;
        RFloat lam1, lam2, lam3, gamm1, vei;
        RFloat vtx, vty, vtz, entr;
        RFloat riemm, riemp;
        RFloat rb, ub, vb, wb, pb, vnb, cb;

        //!ÄÚµã
        rin = prims[IR];
        uin = prims[IU];
        vin = prims[IV];
        win = prims[IW];
        pin = prims[IP];

        //!Ô¶³¡
        roo = prim_inf[IR];
        uoo = prim_inf[IU];
        voo = prim_inf[IV];
        woo = prim_inf[IW];
        poo = prim_inf[IP];

        vno = nxs * uoo + nys * voo + nzs * woo - vgn;
        vni = nxs * uin + nys * vin + nzs * win - vgn;

        coo = sqrt(GPUABS(gama0 * poo / roo));
        cin = sqrt(GPUABS(gama * pin / rin));

        gamm1 = gama - one;

        lam1 = vni;
        lam2 = vni - cin;
        lam3 = vni + cin;

        vei = sqrt(uin * uin + vin * vin + win * win);

        //!³¬ÉùËÙ
        if (vei > cin)
        {
            if (vni >= 0.0)
            {
                for (int m = 0; m < neqn; ++m)
                {
                    primt[m] = prims[m];
                }
            }
            else
            {
                for (int m = 0; m < neqn; ++m)
                {
                    primt[m] = prim_inf[m];
                }
            }
            return;
        }
        //!changed because host is modified
        riemp = vni + 2.0 * cin / gamm1;
        riemm = vno - 2.0 * coo / gamm1;
        vnb   = 0.50 * (riemp + riemm);
        cb    = 0.25 * (riemp - riemm) * gamm1;
        //!Èô³ÌÐòÖ´ÐÐµ½ÕâÀï£¬±íÃ÷Ò»¶¨ÊÇÑÇÉùËÙ
        //!if ( vni >= 0.0 )
        //!host code is upadated!!
        if (vnb >= 0.0)
        {
            //! exit
            //! #ifndef DEBUGPOW
            //! entr = pin / pow(rin, gama);
            //! #else
            //! entr = pin * 1.0;
            //! #endif
            entr = pin / exp(gama * log(rin));
            vtx  = uin - nxs * vni;
            vty  = vin - nys * vni;
            vtz  = win - nzs * vni;
            //!Because the host code is changed here!
            /* 
            riemp = GPUABS( vni ) + 2.0 * cin / gamm1;
            riemm = GPUABS( vno ) - 2.0 * coo / gamm1;
            vnb =  half * ( riemp + riemm );
            cb  =  0.25 * ( riemp - riemm ) * gamm1;
            */
        }
        else
        {
            //!inlet
            //! #ifndef DEBUGPOW
            //! entr = poo / pow(roo, gama);
            //! #else
            //! entr = poo / 1.0;
            //! #endif
            entr = poo / exp(gama * log(roo));
            vtx  = uoo - nxs * vno;
            vty  = voo - nys * vno;
            vtz  = woo - nzs * vno;
            //!Because the host code is changed here
            /*
            riemp = GPUABS( vno ) + 2.0 * coo / gamm1;
            riemm = GPUABS( vni ) - 2.0 * cin / gamm1;
            vnb = - half * ( riemp + riemm );
            cb  =   0.25 * ( riemp - riemm ) * gamm1;
            */
        }

        //!     #ifndef DEBUGPOW
        //! rb  = pow( ( cb * cb / ( entr * gama ) ), static_cast<RFloat>(one) / gamm1 );
        //! #else
        //!     rb  = cb * cb / ( entr * gama );
        //! #endif
        rb = exp(static_cast<RFloat>(one) / gamm1 * log(cb * cb / (entr * gama)));
        ub = vtx + nxs * vnb;
        vb = vty + nys * vnb;
        wb = vtz + nzs * vnb;
        pb = cb * cb * rb / gama;

        primt[IR] = rb;
        primt[IU] = ub;
        primt[IV] = vb;
        primt[IW] = wb;
        primt[IP] = pb;
    }

    __device__ void Farfield_BC_checkMG(RFloat *prims, RFloat *prim_inf, RFloat *primt, RFloat gama0, RFloat gama,
                                        RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem,
                                        int iFace)
    {
        const int IR = 0;
        const int IU = 1;
        const int IV = 2;
        const int IW = 3;
        const int IP = 4;

        int neqn = nl + nchem;

        RFloat rin, uin, vin, win, pin, cin, vni;
        RFloat roo, uoo, voo, woo, poo, coo, vno;
        RFloat lam1, lam2, lam3, gamm1, vei;
        RFloat vtx, vty, vtz, entr;
        RFloat riemm, riemp;
        RFloat rb, ub, vb, wb, pb, vnb, cb;

        //!ÄÚµã
        rin = prims[IR];
        uin = prims[IU];
        vin = prims[IV];
        win = prims[IW];
        pin = prims[IP];

        //!Ô¶³¡
        roo = prim_inf[IR];
        uoo = prim_inf[IU];
        voo = prim_inf[IV];
        woo = prim_inf[IW];
        poo = prim_inf[IP];

        vno = nxs * uoo + nys * voo + nzs * woo - vgn;
        vni = nxs * uin + nys * vin + nzs * win - vgn;

        coo = sqrt(GPUABS(gama0 * poo / roo));
        cin = sqrt(GPUABS(gama * pin / rin));

        gamm1 = gama - one;

        lam1 = vni;
        lam2 = vni - cin;
        lam3 = vni + cin;

        vei = sqrt(uin * uin + vin * vin + win * win);

        //!³¬ÉùËÙ
        if (vei > cin)
        {
            if (vni >= 0.0)
            {
                for (int m = 0; m < neqn; ++m)
                {
                    primt[m] = prims[m];
                }
            }
            else
            {
                for (int m = 0; m < neqn; ++m)
                {
                    primt[m] = prim_inf[m];
                }
            }
            return;
        }
        //!changed because host is modified
        riemp = vni + 2.0 * cin / gamm1;
        riemm = vno - 2.0 * coo / gamm1;
        vnb   = half * (riemp + riemm);
        cb    = 0.25 * (riemp - riemm) * gamm1;
        //!Èô³ÌÐòÖ´ÐÐµ½ÕâÀï£¬±íÃ÷Ò»¶¨ÊÇÑÇÉùËÙ
        //!if ( vni >= 0.0 )
        //!host code is updated here.
        if (vnb >= 0.0)
        {
            //! exit
            //! #ifndef DEBUGPOW
            //! entr = pin / pow(rin, gama);
            //! #else
            //! entr = pin * 1.0;
            //! #endif
            entr = pin / exp(gama * log(rin));
            vtx  = uin - nxs * vni;
            vty  = vin - nys * vni;
            vtz  = win - nzs * vni;
            //!Because the host code is changed here!
            /* 
            riemp = GPUABS( vni ) + 2.0 * cin / gamm1;
            riemm = GPUABS( vno ) - 2.0 * coo / gamm1;
            vnb =  half * ( riemp + riemm );
            cb  =  0.25 * ( riemp - riemm ) * gamm1;
            */
        }
        else
        {
            //!inlet
            //! #ifndef DEBUGPOW
            //! entr = poo / pow(roo, gama);
            //! #else
            //! entr = poo / 1.0;
            //! #endif
            entr = poo / exp(gama * log(roo));
            vtx  = uoo - nxs * vno;
            vty  = voo - nys * vno;
            vtz  = woo - nzs * vno;
            //!Because the host code is changed here
            /*
            riemp = GPUABS( vno ) + 2.0 * coo / gamm1;
            riemm = GPUABS( vni ) - 2.0 * cin / gamm1;
            vnb = - half * ( riemp + riemm );
            cb  =   0.25 * ( riemp - riemm ) * gamm1;
            */
        }

        //!     #ifndef DEBUGPOW
        //! rb  = pow( ( cb * cb / ( entr * gama ) ), static_cast<RFloat>(one) / gamm1 );
        //! #else
        //!     rb  = cb * cb / ( entr * gama );
        //! #endif
        rb        = exp(static_cast<RFloat>(one) / gamm1 * log(cb * cb / (entr * gama)));
        ub        = vtx + nxs * vnb;
        vb        = vty + nys * vnb;
        wb        = vtz + nzs * vnb;
        pb        = cb * cb * rb / gama;
        primt[IR] = rb;
        primt[IU] = ub;
        primt[IV] = vb;
        primt[IW] = wb;
        primt[IP] = pb;
    }

    __device__ void Farfield_BC_check(RFloat *prims, RFloat *prim_inf, RFloat *primt, RFloat gama0, RFloat gama,
                                      RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem, int iFace)
    {
        const int IR = 0;
        const int IU = 1;
        const int IV = 2;
        const int IW = 3;
        const int IP = 4;

        int neqn = nl + nchem;

        RFloat rin, uin, vin, win, pin, cin, vni;
        RFloat roo, uoo, voo, woo, poo, coo, vno;
        RFloat lam1, lam2, lam3, gamm1, vei;
        RFloat vtx, vty, vtz, entr;
        RFloat riemm, riemp;
        RFloat rb, ub, vb, wb, pb, vnb, cb;

        //!ÄÚµã
        rin = prims[IR];
        uin = prims[IU];
        vin = prims[IV];
        win = prims[IW];
        pin = prims[IP];

        //!Ô¶³¡
        roo = prim_inf[IR];
        uoo = prim_inf[IU];
        voo = prim_inf[IV];
        woo = prim_inf[IW];
        poo = prim_inf[IP];

        vno = nxs * uoo + nys * voo + nzs * woo - vgn;
        vni = nxs * uin + nys * vin + nzs * win - vgn;

        coo = sqrt(GPUABS(gama0 * poo / roo));
        cin = sqrt(GPUABS(gama * pin / rin));

        gamm1 = gama - one;

        lam1 = vni;
        lam2 = vni - cin;
        lam3 = vni + cin;

        vei = sqrt(uin * uin + vin * vin + win * win);

        //!³¬ÉùËÙ
        if (vei > cin)
        {
            if (vni >= 0.0)
            {
                for (int m = 0; m < neqn; ++m)
                {
                    primt[m] = prims[m];
                }
            }
            else
            {
                for (int m = 0; m < neqn; ++m)
                {
                    primt[m] = prim_inf[m];
                }
            }
            return;
        }
        //!changed because host is modified
        riemp = vni + 2.0 * cin / gamm1;
        riemm = vno - 2.0 * coo / gamm1;
        vnb   = half * (riemp + riemm);
        cb    = 0.25 * (riemp - riemm) * gamm1;
        //!Èô³ÌÐòÖ´ÐÐµ½ÕâÀï£¬±íÃ÷Ò»¶¨ÊÇÑÇÉùËÙ
        //!if ( vni >= 0.0 )
        //!host code is updated here!!!
        if (vnb >= 0.0)
        {
            //! exit
#ifndef DEBUGPOW
            entr = pin / pow(rin, gama);
#else
            entr = pin * 1.0;
#endif

            vtx = uin - nxs * vni;
            vty = vin - nys * vni;
            vtz = win - nzs * vni;
            /*
        //! exit
        entr = pin / pow(rin, gama);

        vtx = uin - nxs * vni;
        vty = vin - nys * vni;
        vtz = win - nzs * vni;
*/
            //!Because the host code is changed here!
            /* 
        riemp = GPUABS( vni ) + 2.0 * cin / gamm1;
        riemm = GPUABS( vno ) - 2.0 * coo / gamm1;
        vnb =  half * ( riemp + riemm );
        cb  =  0.25 * ( riemp - riemm ) * gamm1;
        */
        }
        else
        {
            //!inlet
#ifndef DEBUGPOW
            entr = poo / pow(roo, gama);
#else
            entr = poo / 1.0;
#endif
            vtx = uoo - nxs * vno;
            vty = voo - nys * vno;
            vtz = woo - nzs * vno;
            /*
        //!inlet
        entr = poo / pow(roo, gama);
        vtx = uoo - nxs * vno;
        vty = voo - nys * vno;
        vtz = woo - nzs * vno;
*/
            //!Because the host code is changed here
            /*
        riemp = GPUABS( vno ) + 2.0 * coo / gamm1;
        riemm = GPUABS( vni ) - 2.0 * cin / gamm1;
        vnb = - half * ( riemp + riemm );
        cb  =   0.25 * ( riemp - riemm ) * gamm1;
*/
        }
        //!Add by zhang xi for testing MAD
        RFloat ub_tmp, vb_tmp, wb_tmp;
        /*
        ub_tmp  =  nxs * vnb;
        vb_tmp  =  nys * vnb;
        wb_tmp  =  nzs * vnb;
        ub = ub_tmp + vtx;
        vb = vb_tmp + vty;
            wb = wb_tmp + vtz;
*/
        //!Add end
        rb = pow((cb * cb / (entr * gama)), static_cast<RFloat>(one) / gamm1);
        ub = vtx + nxs * vnb;
        vb = vty + nys * vnb;
        wb = vtz + nzs * vnb;
        pb = cb * cb * rb / gama;

        primt[IR] = rb;
        primt[IU] = ub;
        primt[IV] = vb;
        primt[IW] = wb;
        primt[IP] = pb;
    }
    /*
    void CallGPULhs()
    {

    }
    
    __global__ void GPULhs()
    {
        
        const int tidx= threadIdx.x;
        const int bidx= blockIdx.x;    
        for ( int m = 0; m < neqn; ++ m )
        {
            for (int icell = bidx* blockDim.x+ tidx;
                 icell < nTotalCell;
                 icell+= blockDim.x* gridDim.x)
            {
                dq[m* nTotalCell+ icell] *= dt[icell] * coef;
            }
        }

    }
    */
    //!void CallGPULoadQ(const int nTotal, const int nl)
    void CallGPULoadQ(Grid *gridIn, Param_NSSolverUnstruct *parameters)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid       = UnstructGridCast(gridIn);
        int           nTotalCell = grid->GetNTotalCell();
        int           nBoundFace = grid->GetNBoundFace();
        int           nTotal     = nTotalCell + nBoundFace;

        int nl = parameters->GetLaminarNumber();

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotal;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPULoadQ, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPULoadQ, 0, gridSize, blockSize);
#endif

        GPULoadQ<<<gridSize, blockSize>>>(d_q_proxy, d_q_ns, nTotal, nl);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPULoadQ(RFloat *q, RFloat *qold, const int nTotal, const int nl)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int m = 0; m < nl; ++m)
        {
            for (int icell = bidx * blockDim.x + tidx; icell < nTotal; icell += gridDim.x * blockDim.x)
            {
                q[m * nTotal + icell] = qold[m * nTotal + icell];
            }
        }
    }
    //!void CallGPUFillField(const int nTotal, const int neqn)
    void CallGPUFillField(Grid *gridIn, const int neqn)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        using namespace GPUControlVariables;

        UnstructGrid *grid       = UnstructGridCast(gridIn);
        int           nTotalCell = grid->GetNTotalCell();
        int           nBoundFace = grid->GetNBoundFace();
        int           nTotal     = nTotalCell + nBoundFace;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotal;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUFillField, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUFillField, 0, gridSize, blockSize);
#endif
        /*
        if(IsFillFieldBeforeStage)
            GPUFillField<<<gridSize, blockSize>>>(d_q_proxy_tmp, d_q_proxy, nTotal, neqn);
        else
            GPUFillField<<<gridSize, blockSize>>>(d_q_proxy, d_q_proxy_tmp, nTotal, neqn);
        */
        if (d_FillField_ns == "q_proxyToq_proxy_temp")
        {
            GPUFillField<<<gridSize, blockSize>>>(d_q_proxy_tmp, d_q_proxy, nTotal, neqn);
        }
        else if (d_FillField_ns == "q_proxy_tempToq_proxy")
        {
            GPUFillField<<<gridSize, blockSize>>>(d_q_proxy, d_q_proxy_tmp, nTotal, neqn);
        }
        else if (d_FillField_ns == "resTodq")
        {
            size_t sizeDq = nTotal * neqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMemcpy(d_dq_ns, d_res_ns, sizeDq, cudaMemcpyDeviceToDevice));
        }
        else
        {
            //!cout<<"Warning: d_FillField_ns ="<<d_FillField_ns<<endl;
            //!exit(1);
        }
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUFillField(RFloat *field1, RFloat *field2, const int nTotal, const int neqn)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int m = 0; m < neqn; ++m)
        {
            for (int icell = bidx * blockDim.x + tidx; icell < nTotal; icell += gridDim.x * blockDim.x)
            {
                field1[m * nTotal + icell] = field2[m * nTotal + icell];
            }
        }
    }
    //!void CallGPULoadResiduals(const int neqn, const int nTotal)
    void CallGPULoadResiduals(Grid *gridIn, const int neqn)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotal;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPULoadResiduals, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPULoadResiduals, 0, gridSize, blockSize);
#endif
        GPULoadResiduals<<<gridSize, blockSize>>>(d_res_ns, d_dq_ns, d_rhs_ns, neqn, nTotal);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    __global__ void GPULoadResiduals(RFloat *res, RFloat *dq, RFloat *rhs, const int neqn, const int nTotal)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int m = 0; m < neqn; ++m)
        {
            for (int icell = bidx * blockDim.x + tidx; icell < nTotal; icell += gridDim.x * blockDim.x)
            {
                res[m * nTotal + icell] = -rhs[m * nTotal + icell];
                dq[m * nTotal + icell]  = 0.0;
            }
        }
    }
    void CallGPURecoverResidual(Grid *gridIn, const int neqn)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        if (grid->GetLevel() != 0)
        {
            KernelLaunchPara((void *)GPURecoverResidualNS, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPURecoverResidualNS, 0, gridSize, blockSize);
#endif
            GPURecoverResidualNS<<<gridSize, blockSize>>>(d_res_ns, d_rhs_ns, neqn, nTotal);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }
    }
    __global__ void GPURecoverResidualNS(RFloat *res, RFloat *rhs, const int neqn, const int nTotal)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int m = 0; m < neqn; ++m)
        {
            for (int icell = bidx * blockDim.x + tidx; icell < nTotal; icell += gridDim.x * blockDim.x)
            {
                res[m * nTotal + icell] = rhs[m * nTotal + icell];
            }
        }
    }
    //!void CallGPULhs(double coef, const int neqn, const int nTotalCell)
    void CallGPURungeKuttaResidual(Grid *gridIn, const int neqn, double coef)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nTotalCell = grid->GetNTotalCell();
        int nTotal     = d_nTotal;
        int gridSize   = 1;
        int blockSize  = 2;

        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPULhs, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPULhs, 0, gridSize, blockSize);
#endif
        GPULhs<<<gridSize, blockSize>>>(d_dq_ns, d_res_ns, d_dt, coef, neqn, nTotal, nTotalCell);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPULhs(RFloat *dq, RFloat *res, RFloat *dt, double coef, const int neqn, const int nTotal,
                           const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int m = 0; m < neqn; ++m)
        {
            for (int icell = bidx * blockDim.x + tidx; icell < nTotalCell; icell += gridDim.x * blockDim.x)
            {
                dq[m * nTotal + icell] = dt[icell] * coef * res[m * nTotal + icell];
            }
        }
    }

    void CallGPUUpdateFlowFieldLUSGSM1(RFloat density_limit, RFloat pressure_limit, RFloat mostNegativePressure,
                                       const int nTotalCell, const int neqn)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int nTotal = d_nTotal;

        int gridSize  = 1;
        int blockSize = 1;

        int    loopLen         = nTotalCell;
        int    regsPerThread   = 54;
        int    residentBlockNo = 4;
        size_t dsMemPerThread  = 2 * neqn * sizeof(RFloat); //! calculate by the intended shared memory use in the kernel
        size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUUpdateFlowFieldM1_S1, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                         dsMemPerBlock, gridSize, blockSize, GPUProp);

        //!KernelLaunchPara((void*)GPUReConstructFaceValueLoop2, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUUpdateFlowFieldM1_S1, dsMemPerBlock, gridSize, blockSize);
#endif
#ifndef DEBUGSolution
        //!GPUUpdateFlowFieldM1_S1<<<gridSize, blockSize, dsMemPerBlock>>>(nTotalCell, nTotal, neqn, d_q_proxy, d_dq_ns, d_q_ns, d_gama_ns,
        GPUUpdateFlowFieldM1_S1<<<gridSize, blockSize, dsMemPerBlock>>>(
            nTotalCell, nTotal, neqn, d_q_ns, d_dq_ns, d_q_ns, d_gama_ns,
            //!d_t_proxy, density_limit, pressure_limit, mostNegativePressure, d_face_number_of_each_cell, d_cell2face,
            //!d_acc_face_number_of_cell, d_left_cell_of_face, d_right_cell_of_face);
            d_t_proxy, density_limit, pressure_limit, mostNegativePressure, d_face_number_of_each_cell, d_cell2Face,
            d_posiCell2Face, d_left_cell_of_face, d_right_cell_of_face);
#else //!GPUSolutionFix makes it is like atomic operation.                                                              \
    //!GPUUpdateFlowFieldM1_S1<<<1, 1, dsMemPerBlock>>>(nTotalCell, nTotal, neqn, d_q_proxy, d_dq_ns, d_q_ns, d_gama_ns,
        GPUUpdateFlowFieldM1_S1<<<1, 1, dsMemPerBlock>>>(nTotalCell, nTotal, neqn, d_q_ns, d_dq_ns, d_q_ns, d_gama_ns,
                                                         d_t_proxy, density_limit, pressure_limit, mostNegativePressure,
                                                         d_face_number_of_each_cell, d_cell2Face, d_posiCell2Face,
                                                         d_left_cell_of_face, d_right_cell_of_face);
#endif
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    //!void CallGPUUpdateFlowFieldM1(RFloat density_limit, RFloat pressure_limit,RFloat mostNegativePressure, const int nTotalCell, const int neqn)
    void CallGPUUpdateFlowField(Grid *gridIn, Param_NSSolverUnstruct *parameters, const int neqn)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid       = UnstructGridCast(gridIn);
        int           nTotalCell = grid->GetNTotalCell();
        int           nBoundFace = grid->GetNBoundFace();
        int           nTotal     = nTotalCell + nBoundFace;

        using namespace GAS_SPACE;
        using namespace IDX;

        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

        RDouble roo            = primitiveVarFarfield[IR];
        RDouble poo            = primitiveVarFarfield[IP];
        RDouble density_limit  = 1.0e-6 * roo;
        RDouble pressure_limit = 1.0e-6 * poo;

        int     nNegativeCell        = 0;
        int     mostNegativeCell     = -1;
        RDouble mostNegativePressure = PHSPACE::LARGE;

        //!int nTotal= d_nTotal;

        int gridSize  = 1;
        int blockSize = 1;

        int    loopLen         = nTotalCell;
        int    regsPerThread   = 54;
        int    residentBlockNo = 4;
        size_t dsMemPerThread  = 2 * neqn * sizeof(RFloat); //! calculate by the intended shared memory use in the kernel
        size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUUpdateFlowFieldM1_S1, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                         dsMemPerBlock, gridSize, blockSize, GPUProp);

        //!KernelLaunchPara((void*)GPUReConstructFaceValueLoop2, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUUpdateFlowFieldM1_S1, dsMemPerBlock, gridSize, blockSize);
#endif

        GPUUpdateFlowFieldM1_S1<<<gridSize, blockSize, dsMemPerBlock>>>(
            nTotalCell, nTotal, neqn, d_q_proxy, d_dq_ns, d_q_ns, d_gama_ns,
            //!d_t_proxy, density_limit, pressure_limit, mostNegativePressure, d_face_number_of_each_cell, d_cell2face,
            //!d_acc_face_number_of_cell, d_left_cell_of_face, d_right_cell_of_face);
            d_t_proxy, density_limit, pressure_limit, mostNegativePressure, d_face_number_of_each_cell, d_cell2Face,
            d_posiCell2Face, d_left_cell_of_face, d_right_cell_of_face);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUUpdateFlowFieldM1(const int nTotalCell, const int nTotal, const int neqn, RFloat *q, RFloat *dq,
                                         RFloat *qnew, RFloat *gamma, RFloat *t, RFloat density_limit,
                                         RFloat pressure_limit, RFloat mostNegativePressure,
                                         int *d_face_number_of_each_cell, int *d_cell2face,
                                         int *d_acc_face_number_of_cell, int *d_left_cell_of_face,
                                         int *d_right_cell_of_face)
    {
        const int IR  = 0;
        const int IU  = 1;
        const int IV  = 2;
        const int IW  = 3;
        const int IP  = 4;
        const int ITT = 0;

        RFloat *prim = new RFloat[neqn];
        RFloat *qtry = new RFloat[neqn];

        RFloat rm, pm, tm, gama;

        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int icell = bidx * blockDim.x + tidx; icell < nTotalCell; icell += gridDim.x * blockDim.x)
        {
            for (int m = 0; m < neqn; ++m)
            {
                prim[m] = q[m * nTotal + icell];
            }

            gama = gamma[icell];

            Primitive2Conservative(prim, gama, qtry);

            for (int m = 0; m < neqn; ++m)
            {
                qtry[m] += dq[m * nTotal + icell];
            }

            tm = t[ITT * nTotal + icell];

            Conservative2Primitive(qtry, gama, prim, tm);

            rm = prim[IR];
            pm = prim[IP];

            if (rm < density_limit || pm < pressure_limit)
            {
                if (pm < mostNegativePressure)
                {
                    mostNegativePressure = pm;
                    //!mostNegativeCell     = icell;
                }

                //!cout << "Negative rho or pressure found, level = " << grid->GetLevel() << ", cell = " << icell << endl;
                //!TK_Exit::ExceptionExit("Negative!");

                //!nNegativeCell += 1;   //! no nNegativeCell count. by sunxu20191215
                //!SolutionFix(grid, prim, icell);
                GPUSolutionFix(prim, d_face_number_of_each_cell, d_cell2face, d_left_cell_of_face, d_right_cell_of_face,
                               neqn, qnew, d_acc_face_number_of_cell, icell, nTotal);
            }

            /*
                if ( nchem == 1 )
                {
                    NormalizePrimitive(prim);
                }
                */
            for (int m = 0; m < neqn; ++m)
            {
                //! qnew is q in grid
                qnew[m * nTotal + icell] = prim[m];
            }
        }

        /*
            if ( nNegativeCell > 0 )
            {
                cout.setf(ios::scientific);
                cout.precision(4);
                cout << "      Warning: negative pressure or density appears in " << nNegativeCell << " cells ... " << endl;
                cout << "               level = " << grid->GetLevel() << endl;
                cout << "               The minimum pressure appears at: ("
                     << xcc[mostNegativeCell] << ", " << ycc[mostNegativeCell] << ", " << zcc[mostNegativeCell] << ")" << endl;
            }
            */
        delete[] prim;
        delete[] qtry;
    }

    __global__ void GPUUpdateFlowFieldM1_S1(const int nTotalCell, const int nTotal, const int neqn, RFloat *q,
                                            RFloat *dq, RFloat *qnew, RFloat *gamma, RFloat *t, RFloat density_limit,
                                            RFloat pressure_limit, RFloat mostNegativePressure,
                                            int *d_face_number_of_each_cell, int *d_cell2Face, int *d_posiCell2Face,
                                            int *d_left_cell_of_face, int *d_right_cell_of_face)
    {
        const int IR  = 0;
        const int IU  = 1;
        const int IV  = 2;
        const int IW  = 3;
        const int IP  = 4;
        const int ITT = 0;

        extern __shared__ RFloat array[];

        RFloat *prim = array + 0 * neqn * blockDim.x;
        RFloat *qtry = array + 1 * neqn * blockDim.x;

        RFloat rm, pm, tm, gama;

        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int icell = bidx * blockDim.x + tidx; icell < nTotalCell; icell += gridDim.x * blockDim.x)
        {
            for (int m = 0; m < neqn; ++m)
            {
                prim[tidx * neqn + m] = q[m * nTotal + icell];
            }

            gama = gamma[icell];

            Primitive2Conservative(prim + tidx * neqn, gama, qtry + tidx * neqn);

            for (int m = 0; m < neqn; ++m)
            {
                qtry[tidx * neqn + m] += dq[m * nTotal + icell];
            }

            tm = t[ITT * nTotal + icell];

            Conservative2Primitive(qtry + tidx * neqn, gama, prim + tidx * neqn, tm);

            rm = prim[tidx * neqn + IR];
            pm = prim[tidx * neqn + IP];

            if (rm < density_limit || pm < pressure_limit)
            {
                //! if(pm < mostNegativePressure)
                //! {
                //!     mostNegativePressure = pm;
                //!     //!mostNegativeCell     = icell;
                //! }

                GPUSolutionFix_Z(prim + tidx * neqn, d_face_number_of_each_cell, d_cell2Face, d_left_cell_of_face,
                                 d_right_cell_of_face, neqn, qnew, d_posiCell2Face, icell, nTotal);
            }

            //! #ifndef DEBUGSolution //!it only exits in CUDAUNITTEST mode.
            //!     if ( rm < density_limit || pm < pressure_limit )
            //!     {
            //!         if(pm < mostNegativePressure)
            //!         {
            //!             mostNegativePressure = pm;
            //!             //!mostNegativeCell     = icell;
            //!         }
            //!         //!cout << "Negative rho or pressure found, level = " << grid->GetLevel() << ", cell = " << icell << endl;
            //!         //!TK_Exit::ExceptionExit("Negative!");

            //!         //!nNegativeCell += 1;   //! no nNegativeCell count. by sunxu20191215
            //!         //!SolutionFix(grid, prim, icell);
            //!         //!GPUSolutionFix(prim+tidx* neqn, d_face_number_of_each_cell, d_cell2face, d_left_cell_of_face,
            //!         //!        d_right_cell_of_face, neqn, qnew, d_acc_face_number_of_cell, icell, nTotal);
            //!         GPUSolutionFix_Z(prim+tidx* neqn, d_face_number_of_each_cell, d_cell2Face, d_left_cell_of_face, d_right_cell_of_face, neqn, qnew, d_posiCell2Face, icell, nTotal);
            //!     }
            //! #else
            //!     //!if the macro is defined, no branch exists, all of cells will get modifications.
            //!         GPUSolutionFix_Z(prim+tidx* neqn, d_face_number_of_each_cell, d_cell2Face, d_left_cell_of_face, d_right_cell_of_face, neqn, qnew, d_posiCell2Face, icell, nTotal);
            //! #endif
            /*
            if ( nchem == 1 )
            {
                NormalizePrimitive(prim);
            }
            */
            for (int m = 0; m < neqn; ++m)
            {
                //! qnew is q in grid
                qnew[m * nTotal + icell] = prim[tidx * neqn + m];
            }
        }
        /*
            if ( nNegativeCell > 0 )
            {
                cout.setf(ios::scientific);
                cout.precision(4);
                cout << "      Warning: negative pressure or density appears in " << nNegativeCell << " cells ... " << endl;
                cout << "               level = " << grid->GetLevel() << endl;
                cout << "               The minimum pressure appears at: ("
                     << xcc[mostNegativeCell] << ", " << ycc[mostNegativeCell] << ", " << zcc[mostNegativeCell] << ")" << endl;
            }
            */
        //!delete [] prim;
        //!delete [] qtry;
    }
    __device__ void Primitive2Conservative(RFloat *prim, RFloat gama, RFloat *q)
    {
        const int IR  = 0;
        const int IU  = 1;
        const int IV  = 2;
        const int IW  = 3;
        const int IP  = 4;
        const int IRU = 1;
        const int IRV = 2;
        const int IRW = 3;
        const int IRE = 4;

        RFloat  em;
        RFloat &density  = prim[IR];
        RFloat &um       = prim[IU];
        RFloat &vm       = prim[IV];
        RFloat &wm       = prim[IW];
        RFloat &pressure = prim[IP];
        /* chemical part is not supported by sunxu 20191215
                //!Obtain the mass fraction of the last species whose label in the collection is ns-1.
            if ( nchem == 1 )
            {
            NormalizePrimitive(prim);
            }
            */

        //!obtain the total internal energy Em.
        ComputeInternalEnergy(prim, gama, em);

        q[IR]  = density;
        q[IRU] = density * um;
        q[IRV] = density * vm;
        q[IRW] = density * wm;
        q[IRE] = density * em;
    }
    __device__ void ComputeInternalEnergy(RFloat *prim, RFloat gama_in, RFloat &em)
    {
        const int IR   = 0;
        const int IU   = 1;
        const int IV   = 2;
        const int IW   = 3;
        const int IP   = 4;
        double    one  = 1.0;
        double    half = 0.5;
        RFloat    tm, v2, omav, hm;

        RFloat &rm = prim[IR];
        RFloat &um = prim[IU];
        RFloat &vm = prim[IV];
        RFloat &wm = prim[IW];
        RFloat &pm = prim[IP];
        v2         = um * um + vm * vm + wm * wm;

        //!if ( nchem == 0 )   //! chemical part is not supported by sunxu 20191215
        {
            em = (pm / rm) / (gama_in - one) + half * v2;
        }
        /*
        else
        {
        ComputeMolecularWeightReciprocal(prim, omav);
        tm = pm / ( coefficientOfStateEquation * rm * omav );
        ComputeSpeciesEnthalpy(tm, hintSpecies);
        ComputeMixtureByPrimitive(prim, hintSpecies, hm);
        em = hm - pm / rm + half * v2;
        }
        */
    }
    __device__ void Conservative2Primitive(RFloat *q, RFloat gama, RFloat *prim, RFloat &temperature)
    {
        const int IR   = 0;
        const int IU   = 1;
        const int IV   = 2;
        const int IW   = 3;
        const int IP   = 4;
        double    one  = 1.0;
        double    zero = 0.0;
        double    half = 0.5;

        RFloat density, oDensity, um, vm, wm, pressure, rem, v2, reint;
        density = q[IR];
        if (density <= zero)
        {
            prim[IR] = density;
            return;
        }

        oDensity = 1.0 / density;
        um       = q[IU] * oDensity;
        vm       = q[IV] * oDensity;
        wm       = q[IW] * oDensity;
        rem      = q[IP];
        v2       = um * um + vm * vm + wm * wm;

        reint = rem - half * density * v2;

        prim[IR] = density;
        prim[IU] = um;
        prim[IV] = vm;
        prim[IW] = wm;
        /* chemical is not supported by sunxu 20191215
            for ( int m = nm; m < nl; ++ m )
            {
            prim[m] = MIN(static_cast<RFloat>(one), ABS( q[m] * oDensity ));
            }
            
            if ( nchem == 1 )
            {
            NormalizePrimitive(prim);
            }
            */
        //!GetPressure(prim, gama, reint, temperature, pressure);
        pressure = (gama - one) * reint;

        prim[IP] = pressure;
    }

    __device__ void GPUSolutionFix(RFloat *prim, int *face_number_of_each_cell, int *cell2face, int *left_cell_of_face,
                                   int *right_cell_of_face, int neqn, RFloat *q, int *d_acc_face_number_of_cell,
                                   int icell, int nTotal)
    {
        int          IR   = 0;
        int          IP   = 4;
        const double zero = 0.0;
        for (int m = 0; m < neqn; ++m)
        {
            prim[m] = zero;
        }

        int    iface, le, re, isur;
        RFloat vol_sum, ovol_sum, vsur;
        vol_sum = 0.0;
        for (int icell2face = 0; icell2face < face_number_of_each_cell[icell]; ++icell2face)
        {
            //!iface  = cell2face[icell][icell2face];
            iface = *(cell2face + d_acc_face_number_of_cell[icell] + icell2face);
            le    = left_cell_of_face[iface];
            re    = right_cell_of_face[iface];

            isur = le;
            if (icell == le) isur = re;

            //!vsur = vol[isur];
            vsur = one;

            vol_sum += vsur;

            for (int m = 0; m < neqn; ++m)
            {
                RFloat f = q[m * nTotal + icell];
                if (m == IR || m == IP)
                {
                    f = GPUABS(f);
                }
                prim[m] += f * vsur;
            }
        }

        ovol_sum = 1.0 / vol_sum;

        for (int m = 0; m < neqn; ++m)
        {
            prim[m] *= ovol_sum;
        }
    }

    __device__ void GPUSolutionFix_Z(RFloat *prim, int *face_number_of_each_cell, int *cell2Face,
                                     int *left_cell_of_face, int *right_cell_of_face, int neqn, RFloat *q,
                                     int *posiCell2Face, int icell, int nTotal)
    {
        int          IR   = 0;
        int          IP   = 4;
        const double zero = 0.0;
        const double one  = 1.0;
        for (int m = 0; m < neqn; ++m)
        {
            prim[m] = zero;
        }

        int    iface, le, re, isur;
        RFloat vol_sum, ovol_sum, vsur;
        vol_sum = 0.0;
        for (int icell2face = 0; icell2face < face_number_of_each_cell[icell]; ++icell2face)
        {
            //!iface  = cell2face[icell][icell2face];
            //!iface= *(cell2face+d_acc_face_number_of_cell[icell]+icell2face);
            iface = cell2Face[posiCell2Face[icell] + icell2face];
            le    = left_cell_of_face[iface];
            re    = right_cell_of_face[iface];
            isur  = le;
            if (icell == le) isur = re;

            //!vsur = vol[isur];
            vsur = one;

            vol_sum += vsur;

            for (int m = 0; m < neqn; ++m)
            {
                //!RFloat f = q[m*nTotal+icell];
                RFloat f = q[m * nTotal + isur];
                if (m == IR || m == IP)
                {
                    f = GPUABS(f);
                }
                prim[m] += f * vsur;
            }
        }

        ovol_sum = 1.0 / vol_sum;

        for (int m = 0; m < neqn; ++m)
        {
            prim[m] *= ovol_sum;
        }
    }

    void CallGPUZeroResiduals(Grid *gridIn, const int nl)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid       = UnstructGridCast(gridIn);
        int           nTotalCell = grid->GetNTotalCell();
        int           nBoundFace = grid->GetNBoundFace();
        int           nTotal     = nTotalCell + nBoundFace;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotal;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUZeroResiduals, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUZeroResiduals, 0, gridSize, blockSize);
#endif

        GPUZeroResiduals<<<gridSize, blockSize>>>(nl, d_res_ns, nTotal);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    __global__ void GPUZeroResiduals(const int nl, RFloat *res, const int nTotal)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int m = 0; m < nl; ++m)
        {
            for (int icell = bidx * blockDim.x + tidx; icell < nTotal; icell += gridDim.x * blockDim.x)
            {
                res[m * nTotal + icell] = 0.0;
            }
        }
    }
    void CallGPUDtIni(const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUDtIni, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUDtIni, 0, gridSize, blockSize);
#endif
        GPUDtIni<<<gridSize, blockSize>>>(d_dt, nTotalCell);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    __global__ void GPUDtIni(RFloat *dt, const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            dt[iCell] = 0.0;
        }
    }
    void CallGPUDtCFL(RFloat cfl, const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUDtCFL, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUDtCFL, 0, gridSize, blockSize);
#endif
//!GPUDtCFL<<<gridSize, blockSize>>>(d_dt, d_invSpectrumRadius,d_vol, cfl,  nTotalCell);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    __global__ void GPUDtCFL(RFloat *d_dt, RFloat *d_spectralRadius, RFloat *vol, RDouble *d_CFLCell,
                             const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            d_dt[iCell] = d_CFLCell[iCell] * vol[iCell] / d_spectralRadius[iCell];
        }
    }
    void CallGPULimitCFL(Grid *gridIn, const Param_NSSolverUnstruct *parameters, const Param_NSSolver *parametersns)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid       = UnstructGridCast(gridIn);
        int           nTotalCell = grid->GetNTotalCell();
        int           nBoundFace = grid->GetNBoundFace();
        int           nTotal     = nTotalCell + nBoundFace;
        int           level      = grid->GetLevel();

        //!NSSolverUnstruct mySolver;
        //!Param_NSSolverUnstruct *parameters = mySolver.GetControlParameters();
        //!Param_NSSolverUnstruct *parameters = GetControlParameters();
        //! ComputeCFL
        RDouble CFL        = 0.0;
        RDouble partialCFL = -1.0;
        {
            int iterationStep = 0;

            //!Param_NSSolver *parametersns = mySolver.GetControlParameters();
            //! Steady.
            int nMGLevel = parametersns->GetNMGLevel();
            if (nMGLevel <= 1)
            {
                //! Finest grid.
                GlobalDataBase::GetData("outnstep", &iterationStep, PHINT, 1);
            }
            else
            {
                //! Coarse grid.
                GlobalDataBase::GetData("newnstep", &iterationStep, PHINT, 1);
            }
            //! Unsteady.
            int isUnsteady = parametersns->GetIsUnsteady();
            if (1 == isUnsteady)
            {
                iterationStep = GlobalDataBase::GetIntParaFromDB("innstep");
            }

            int cflMethod = 0;
            if (GlobalDataBase::IsExist("CFLMethod", PHINT, 1))
            {
                GlobalDataBase::GetData("CFLMethod", &cflMethod, PHINT, 1);
            }
            else
            {
                GlobalDataBase::UpdateData("CFLMethod", &cflMethod, PHINT, 1);
            }

            int     cflNstep = parametersns->GetCFLVaryStep();
            RDouble cflStart = parametersns->GetCFLStart();
            RDouble cflEnd   = parametersns->GetCFLEnd();
            RDouble cfl      = 0.0;

            if (partialCFL > 0)
            {
                cflEnd = partialCFL;
            }

            if (iterationStep >= cflNstep)
            {
                cfl = cflEnd;
            }
            else
            {
                if (0 == cflMethod)
                {
                    cfl = cflStart + (cflEnd - cflStart) * iterationStep / cflNstep;
                }
                else
                {
                    RDouble CFLratio = cflEnd / (cflStart + TINY);
                    cfl              = cflStart * pow(CFLratio, (iterationStep - 1) * 1.0 / cflNstep);
                }
            }
            CFL = cfl;
        }

        int isInMGIniting = GlobalDataBase::GetIntParaFromDB("isInMGIniting");
        if (!grid->IsFinestGrid() && !isInMGIniting)
        {
            RDouble mgCFLScale = parameters->GetMgCFLScale();
            CFL *= mgCFLScale;
        }

        int gridSize       = 1;
        int blockSize      = 1;
        int blockSizeLimit = 0;

        int loopLen = nTotalCell;
        KernelLaunchPara((void *)GPUSpectralRadiusInvis, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSpectralRadiusInvis, 0, gridSize, blockSize);
#endif
        GPUCFLCellIni<<<gridSize, blockSize>>>(d_CFLCell, CFL, nTotalCell);

        int nChemical = parameters->GetChemicalFlag();
        if (nChemical != 0)
        {
            return;
        }

        RDouble CFLStart = GlobalDataBase::GetDoubleParaFromDB("CFLStart");
        RDouble pTotal   = GlobalDataBase::GetDoubleParaFromDB("pTotal");

        const RDouble LOWER       = 0.001;
        const RDouble UPPER       = 0.02;
        RDouble       pLowerLimit = LOWER * pTotal;
        RDouble       pUpperLimit = UPPER * pTotal;

        //!！adjust CFLCell according to p in the cell.
        GPUComputeCFLCell<<<gridSize, blockSize>>>(d_CFLCell, d_q_ns, CFL, CFLStart, pLowerLimit, pUpperLimit,
                                                   nTotalCell, nTotal);
        if (0 == level)
        {
            GPUCorCFLCell<<<gridSize, blockSize>>>(d_CFLCell, d_cellSkewness, nTotalCell);
        }
    }
    __global__ void GPUCFLCellIni(RDouble *d_CFLCell, RDouble CFL, const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            d_CFLCell[iCell] = CFL;
        }
    }
    __global__ void GPUComputeCFLCell(RDouble *d_CFLCell, RFloat *d_q_ns, RDouble CFL, RDouble CFLStart,
                                      RDouble pLowerLimit, RDouble pUpperLimit, const int nTotalCell, const int nTotal)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        const int ip   = 4;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            if (d_q_ns[ip * nTotal + iCell] < pLowerLimit)
            {
                d_CFLCell[iCell] = CFLStart;
            }
            else if (d_q_ns[ip * nTotal + iCell] < pUpperLimit)
            {
                d_CFLCell[iCell] =
                    (d_q_ns[ip * nTotal + iCell] - pUpperLimit) * (CFL - CFLStart) / (pUpperLimit - pLowerLimit) + CFL;
            }
        }
    }
    __global__ void GPUCorCFLCell(RDouble *d_CFLCell, RFloat *d_cellSkewness, const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            if (d_cellSkewness[iCell] < 0.1)
            {
                d_CFLCell[iCell] *= 0.1;
            }
            else if (d_cellSkewness[iCell] < 10)
            {
                d_CFLCell[iCell] *= (1 - exp(0 - d_cellSkewness[iCell]));
            }
        }
    }
    void CallGPULocalTimeStep(Grid *gridIn, const Param_NSSolverUnstruct *parameters,
                              const Param_NSSolver *parametersns)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid       = UnstructGridCast(gridIn);
        int           nTotalCell = grid->GetNTotalCell();

        //!NSSolverUnstruct mySolver;
        //!Param_NSSolverUnstruct *parameters = mySolver.GetControlParameters();
        //!Param_NSSolverUnstruct *parameters = GetControlParameters();
        bool isViscous = parameters->IsViscous();

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;

        RDouble *d_spectralRadius;
        HANDLE_API_ERR(cudaMalloc((void **)&d_CFLCell, nTotalCell * sizeof(RDouble)));
        HANDLE_API_ERR(cudaMalloc((void **)&d_spectralRadius, nTotalCell * sizeof(RDouble)));

        cudaMemset(d_CFLCell, 0, nTotalCell * sizeof(RDouble));
        cudaMemset(d_spectralRadius, 0, nTotalCell * sizeof(RDouble));

        KernelLaunchPara((void *)GPUSpectralRadiusInvis, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSpectralRadiusInvis, 0, gridSize, blockSize);
#endif
        GPUSpectralRadiusInvis<<<gridSize, blockSize>>>(d_spectralRadius, d_invSpectrumRadius, nTotalCell);
        if (isViscous)
        {
            GPUSpectralRadiusVis<<<gridSize, blockSize>>>(d_spectralRadius, d_visSpectrumRadius, nTotalCell);
        }

        CallGPULimitCFL(gridIn, parameters, parametersns);

        GPUDtIni<<<gridSize, blockSize>>>(d_dt, nTotalCell);

        GPUDtCFL<<<gridSize, blockSize>>>(d_dt, d_spectralRadius, d_vol, d_CFLCell, nTotalCell);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        HANDLE_API_ERR(cudaFree(d_CFLCell));
        HANDLE_API_ERR(cudaFree(d_spectralRadius));
    }
    __global__ void GPUSpectralRadiusInvis(RFloat *d_spectralRadius, RFloat *d_invSpectrumRadius, const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            d_spectralRadius[iCell] = d_invSpectrumRadius[iCell];
        }
    }
    __global__ void GPUSpectralRadiusVis(RFloat *d_spectralRadius, RFloat *d_visSpectrumRadius, const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            d_spectralRadius[iCell] += d_visSpectrumRadius[iCell];
        }
    }
    void CallGPUVisSpectrumRadiusIni(const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUVisSpectrumRadiusIni, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUVisSpectrumRadiusIni, 0, gridSize, blockSize);
#endif
        GPUVisSpectrumRadiusIni<<<gridSize, blockSize>>>(d_visSpectrumRadius, nTotalCell);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUVisSpectrumRadiusIni(RFloat *d_visSpectrumRadius, const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            d_visSpectrumRadius[iCell] = 0.0;
        }
    }

    void CallGPUDtvIni(RFloat cfl, const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUDtvIni, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUDtvIni, 0, gridSize, blockSize);
#endif
        GPUDtvIni<<<gridSize, blockSize>>>(cfl, d_dtv, d_vol, d_visSpectrumRadius, nTotalCell);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    __global__ void GPUDtvIni(RFloat cfl, RFloat *dtv, RFloat *vol, RFloat *visSpectrumRadius, const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            dtv[iCell] = cfl * vol[iCell] / visSpectrumRadius[iCell];
        }
    }

    void CallGPUDtFinal(const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUDtFinal, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUDtFinal, 0, gridSize, blockSize);
#endif
        GPUDtFinal<<<gridSize, blockSize>>>(d_dt, d_dtv, nTotalCell);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    __global__ void GPUDtFinal(RFloat *dt, RFloat *dtv, const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            dt[iCell] = dt[iCell] * dtv[iCell] / (dt[iCell] + dtv[iCell]);
        }
    }

    void CallGPUReduceMaxTimeStep(double ktmax, RFloat globalMinDt, const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;

        if (ktmax > 0)
        {
            KernelLaunchPara((void *)GPUReduceMaxTimeStep1, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUReduceMaxTimeStep1, 0, gridSize, blockSize);
#endif
            //! #ifndef CUDA_AWARE_MPI
            //!             GPUReduceMaxTimeStep1<<<gridSize, blockSize>>>(d_dt, ktmax, globalMinDt, nTotalCell);
            //! #else
            //!             GPUReduceMaxTimeStep1AWARE<<<gridSize, blockSize>>>(d_dt, ktmax, d_globalMinDt, nTotalCell);

            //! #endif

            GPUReduceMaxTimeStep1<<<gridSize, blockSize>>>(d_dt, ktmax, globalMinDt, nTotalCell);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }
        KernelLaunchPara((void *)GPUReduceMaxTimeStep2, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUReduceMaxTimeStep2, 0, gridSize, blockSize);
#endif
        GPUReduceMaxTimeStep2<<<gridSize, blockSize>>>(d_dt, d_vol, nTotalCell);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUReduceMaxTimeStepLoop1(RDouble ktmax, RDouble globalMinDt, const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUReduceMaxTimeStep1, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUReduceMaxTimeStep1, 0, gridSize, blockSize);
#endif
        GPUReduceMaxTimeStep1<<<gridSize, blockSize>>>(d_dt, ktmax, globalMinDt, nTotalCell);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUReduceMaxTimeStepLoop2(const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUReduceMaxTimeStep2, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUReduceMaxTimeStep2, 0, gridSize, blockSize);
#endif

        GPUReduceMaxTimeStep2<<<gridSize, blockSize>>>(d_dt, d_vol, nTotalCell);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    //!void CallGPUReduceMaxTimeStepOneProcess(double ktmax, const int nTotalCell)
    void CallGPUReduceMaxTimeStep(Grid *gridIn)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid       = UnstructGridCast(gridIn);
        int           nTotalCell = grid->GetNTotalCell();

        RDouble ktmax = GlobalDataBase::GetDoubleParaFromDB("ktmax");

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;

        if (ktmax > 0)
        {
            KernelLaunchPara((void *)GPUReduceMaxTimeStep1OneProcess, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUReduceMaxTimeStep1OneProcess, 0, gridSize, blockSize);
#endif
            //!It should be noted that d_minDt is an array of double even only with one element. RFloat * should be used here.
            GPUReduceMaxTimeStep1OneProcess<<<gridSize, blockSize>>>(d_dt, ktmax, d_minDt, nTotalCell);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }
        KernelLaunchPara((void *)GPUReduceMaxTimeStep2, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUReduceMaxTimeStep2, 0, gridSize, blockSize);
#endif
        GPUReduceMaxTimeStep2<<<gridSize, blockSize>>>(d_dt, d_vol, nTotalCell);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        CallGPUSetGhostCell(gridIn);
    }

    __global__ void GPUReduceMaxTimeStep1OneProcess(RFloat *dt, double ktmax, RFloat *minDt, const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            dt[iCell] = GPUMIN(dt[iCell], double(ktmax * minDt[0]));
        }
    }

    __global__ void GPUReduceMaxTimeStep1AWARE(RFloat *dt, double ktmax, double *globalMinDt, const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            dt[iCell] = GPUMIN(dt[iCell], double(ktmax * globalMinDt[0]));
        }
    }

    __global__ void GPUReduceMaxTimeStep1(RFloat *dt, const double ktmax, const double globalMinDt,
                                          const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            dt[iCell] = GPUMIN(dt[iCell], ktmax * globalMinDt);
        }
    }
    __global__ void GPUReduceMaxTimeStep2(RFloat *dt, RFloat *vol, const int nTotalCell)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iCell = bidx * blockDim.x + tidx; iCell < nTotalCell; iCell += gridDim.x * blockDim.x)
        {
            //!dt[iCell] = GPUMAX(dt[iCell] / vol[iCell], 1.0e-12);
            dt[iCell] = dt[iCell] / vol[iCell];
        }
    }

    void CallGPUSetGhostCell(Grid *gridIn)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nBoundFace = grid->GetNBoundFace();

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nBoundFace;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUSetGhostCell, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSetGhostCell, 0, gridSize, blockSize);
#endif
        GPUSetGhostCell<<<gridSize, blockSize>>>(d_dt, d_boundaryType, d_left_cell_of_face, d_right_cell_of_face,
                                                 nBoundFace);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUSetGhostCell(RFloat *f, int *d_boundaryType, int *left_cell_of_face, int *right_cell_of_face,
                                    const int nBoundFace)
    {
        int       le, re;
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int iface = bidx * blockDim.x + tidx; iface < nBoundFace; iface += gridDim.x * blockDim.x)
        {
            //!if ( bctype == PHENGLEI::INTERFACE)
            if (d_boundaryType[iface] == -1) continue;
            le    = left_cell_of_face[iface];
            re    = right_cell_of_face[iface];
            f[re] = f[le];
        }
    }

    void CallGPUCopyMinMaxTime(const double minDt, const double maxDt)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        RDouble *h_minDt = new RDouble[1];
        RDouble *h_maxDt = new RDouble[1];
        h_minDt[0]       = minDt;
        h_maxDt[0]       = maxDt;
        HANDLE_API_ERR(cudaMemcpy(d_minDt, h_minDt, sizeof(RDouble), cudaMemcpyHostToDevice));
        HANDLE_API_ERR(cudaMemcpy(d_maxDt, h_maxDt, sizeof(RDouble), cudaMemcpyHostToDevice));
        delete[] h_minDt;
        delete[] h_maxDt;
    }

    void CallGPUComputeMinTimeStep(const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        /*
        int gridSize= 128;
                int blockSize= 256;
                int loopLen= nTotalCell;

        RFloat *bmin_tmp, *bmax_tmp;
        size_t bSize= gridSize* sizeof(RFloat);
        HANDLE_API_ERR(cudaMalloc((void**)&bmin_tmp, bSize));
        HANDLE_API_ERR(cudaMalloc((void**)&bmax_tmp, bSize));
        size_t dsMemPerBlock= 2*blockSize* sizeof(RFloat);    
        GPUComputeMinMaxStep1<<<gridSize, blockSize, dsMemPerBlock>>>(d_dt, bmin_tmp, bmax_tmp, nTotalCell);
        #ifdef KERNELLAUNCHTEST 
        HANDLE_KERNEL_ERR();
        #endif
        dsMemPerBlock= 2*gridSize* sizeof(RFloat);
        GPUComputeMinMaxStep2<<<1, gridSize, dsMemPerBlock>>>(bmin_tmp, bmax_tmp, d_minDt, d_maxDt);
        #ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
        #endif
        */
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
        /*
        RFloat *bmin_tmp, *bmax_tmp;
        size_t bSize= gridSize* sizeof(RFloat);
        HANDLE_API_ERR(cudaMalloc((void**)&bmin_tmp, bSize));
        HANDLE_API_ERR(cudaMalloc((void**)&bmax_tmp, bSize));
        */

//!KernelLaunchPara((void*)GPUReConstructFaceValueLoop2, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUComputeMinMaxStep1, dsMemPerBlock, gridSize, blockSize);
#endif

        //!GPUComputeMinMaxStep1<<<gridSize, blockSize, dsMemPerBlock>>>(d_dt, bmin_tmp, bmax_tmp, nTotalCell);
        GPUComputeMinMaxStep1<<<gridSize, blockSize, dsMemPerBlock>>>(d_dt, d_minLocalTime, d_maxLocalTime, nTotalCell);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        //!GPUComputeMinMaxStep2<<<1, gridSize, dsMemPerBlock>>>(bmin_tmp, bmax_tmp, d_minDt, d_maxDt);
        GPUComputeMinMaxStep2<<<1, gridSize, dsMemPerBlock>>>(d_minLocalTime, d_maxLocalTime, d_minDt, d_maxDt);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        //!HANDLE_API_ERR(cudaFree(bmin_tmp));
        //!HANDLE_API_ERR(cudaFree(bmax_tmp));
    }
    __global__ void GPUComputeMinMaxStep1(RFloat *dt, RFloat *bmin_tmp, RFloat *bmax_tmp, const int nTotalCell)
    {
        //!_shared__ RFloat tmp_min[blockDim.x];
        //!__shared__ RFloat tmp_max[blockDim.x];
        extern __shared__ RFloat array[];
        RFloat                  *tmp_min = array;
        RFloat                  *tmp_max = array + blockDim.x;

        double    temp_min = 1.0e40;
        double    temp_max = -1.0e40;
        const int tidx     = threadIdx.x;
        const int bidx     = blockIdx.x;
        const int t_n      = blockDim.x * gridDim.x;
        //! get the min, max among the elements a thread deals with
        for (int icell = bidx * blockDim.x + tidx; icell < nTotalCell; icell += t_n)
        {
            temp_min = GPUMIN(temp_min, dt[icell]);
            temp_max = GPUMAX(temp_max, dt[icell]);
        }
        tmp_min[tidx] = temp_min;
        tmp_max[tidx] = temp_max;
        __syncthreads();
        //! get the min max among the threads of same block
        for (int i = blockDim.x / 2; i > 0; i /= 2)
        {
            if (tidx < i)
            {
                tmp_min[tidx] = GPUMIN(tmp_min[tidx], tmp_min[tidx + i]);
                tmp_max[tidx] = GPUMAX(tmp_max[tidx], tmp_max[tidx + i]);
            }
            __syncthreads();
        }
        if (tidx == 0)
        {
            bmin_tmp[bidx] = tmp_min[0];
            bmax_tmp[bidx] = tmp_max[0];
        }
    }
    __global__ void GPUComputeMinMaxStep2(RFloat *bmin_tmp, RFloat *bmax_tmp, RFloat *minDt, RFloat *maxDt)
    {
        //!__shared__ RFloat tmp_min[blockDim.x];
        //!__shared__ RFloat tmp_max[blockDim.x];
        extern __shared__ RFloat array[];
        RFloat                  *tmp_min = array;
        RFloat                  *tmp_max = array + blockDim.x;

        const int tidx = threadIdx.x;
        tmp_min[tidx]  = bmin_tmp[tidx];
        tmp_max[tidx]  = bmax_tmp[tidx];
        __syncthreads();

        for (int i = blockDim.x / 2; i > 0; i /= 2)
        {
            if (tidx < i)
            {
                tmp_min[tidx] = GPUMIN(tmp_min[tidx], tmp_min[tidx + i]);
                tmp_max[tidx] = GPUMAX(tmp_max[tidx], tmp_max[tidx + i]);
            }
            __syncthreads();
        }
        if (tidx == 0)
        {
            minDt[0] = tmp_min[0];
            maxDt[0] = tmp_max[0];
        }
    }
    void CallGPUComputeViscousCoefficientWithoutChemical(Grid *gridIn, Param_NSSolverUnstruct *parameters)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        using namespace IDX;

        UnstructGrid *grid       = UnstructGridCast(gridIn);
        int           nTotalCell = grid->GetNTotalCell();
        int           nBoundFace = grid->GetNBoundFace();
        int           nTotal     = nTotalCell + nBoundFace;

        RDouble tsuth = parameters->GetNonDimensionalSutherlandTemperature();

        RDouble viscousLaminarMin;
        GlobalDataBase::GetData("visl_min", &viscousLaminarMin, PHDOUBLE, 1);

        int gridSize  = 1;
        int blockSize = 1;

        int loopLen        = nTotal;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUCompViscousCoef, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompViscousCoef, 0, gridSize, blockSize);
#endif
        GPUCompViscousCoef<<<gridSize, blockSize>>>(d_visl, d_t_proxy, tsuth, viscousLaminarMin, nTotal, IDX::ITT);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUCompViscousCoef(RFloat *visl, RFloat *t, double tsuth, double visl_min, const int nTotal,
                                       const int ITT)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int icell = bidx * blockDim.x + tidx; icell < nTotal; icell += gridDim.x * blockDim.x)
        {
            visl[icell] = t[ITT * nTotal + icell] * sqrt(t[ITT * nTotal + icell]) * (1.0 + tsuth)
                          / (t[ITT * nTotal + icell] + tsuth);
            visl[icell] = GPUMAX(visl_min, visl[icell]);
        }
    }
    void CallGPUStoreRhsByResidual(Grid *gridIn, const int neqn)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        UnstructGrid *grid       = UnstructGridCast(gridIn);
        int           nTotalCell = grid->GetNTotalCell();
        int           nBoundFace = grid->GetNBoundFace();
        int           nTotal     = nTotalCell + nBoundFace;
        int           gridSize   = 1;
        int           blockSize  = 1;

        int loopLen        = nTotal;
        int blockSizeLimit = 0;

        KernelLaunchPara((void *)GPUStoreRhsByResidual, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUStoreRhsByResidual, 0, gridSize, blockSize);
#endif
        GPUStoreRhsByResidual<<<gridSize, blockSize>>>(d_rhs_ns, d_res_ns, neqn, nTotal);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUStoreRhsByResidual(RFloat *rhs, RFloat *res, const int neqn, const int nTotal)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int m = 0; m < neqn; ++m)
        {
            for (int icell = bidx * blockDim.x + tidx; icell < nTotal; icell += gridDim.x * blockDim.x)
            {
                rhs[m * nTotal + icell] = res[m * nTotal + icell];
            }
        }
    }
} //! namespace GPUNSSolverUnstruct
