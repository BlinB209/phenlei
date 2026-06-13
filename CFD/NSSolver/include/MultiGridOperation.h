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
//! @file      MultiGridOperation.h
//! @brief     multigrid operation.
//! @author    He Xin.

#pragma once
#include "Geo_StructGrid.h"
#include "Geo_UnstructGrid.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "Glb_Dimension.h"

namespace PHSPACE
{
template < typename T >
void RestrictQ(StructGrid *fgrid, PHArray< T, 2> &fq, StructGrid *cgrid, PHArray< T, 2> &cq)
{
    RDouble3D & fvol = * (fgrid->GetCellVolume());
    RDouble3D & cvol = * (cgrid->GetCellVolume());

    int fni = fgrid->GetNI();
    int fnj = fgrid->GetNJ();

    int cni = cgrid->GetNI();
    int cnj = cgrid->GetNJ();

    cq = 0.0;

    int istp = (fni - 1) / (cni - 1);
    int jstp = (fnj - 1) / (cnj - 1);

    int id = istp / 2;
    int jd = jstp / 2;

    int i0,j0,i,j,ip,jp;

    int k = 1;

    RDouble v1,v2,v3,v4,ovol;
    for (j0 = 1; j0 <= cnj - 1; ++ j0)
    {
        j  = (j0 - 1) * jstp + 1;
        jp = j + jd;
        for (i0 = 1; i0 <= cni - 1; ++ i0)
        {
            i  = (i0 - 1) * istp + 1;
            ip = i + id;
            v1 = fvol(i ,j ,k);
            v2 = fvol(ip,j ,k);
            v3 = fvol(ip,jp,k);
            v4 = fvol(i ,jp,k);

            ovol = one / (v1 + v2 + v3 + v4);
            cq(i0,j0) = (fq(i ,j) * v1 + fq(ip,j) * v2 + fq(ip,jp) * v3 + fq(i ,jp) * v4) * ovol;
        }
    }

    //! And secondly, restrict the first layer of ghost cells.
    for (j0 = 0; j0 <= cnj; ++ j0)
    {
        j  = (j0 - 1) * jstp + 1;
        jp = j + jd;
        i  = 0;
        i0 = 0;
        v1 = fvol(i ,j ,k);
        v2 = fvol(i ,jp,k);
        ovol = one / (v1 + v2);
        cq(i0,j0) = (fq(i ,j) * v1 + fq(i ,jp) * v2) * ovol;

        i  = fni;
        i0 = cni;
        v1 = fvol(i ,j ,k);
        v2 = fvol(i ,jp,k);
        ovol = one / (v1 + v2);
        cq(i0,j0) = (fq(i ,j) * v1 + fq(i ,jp) * v2) * ovol;
    }

    for (i0 = 0; i0 <= cni; ++ i0)
    {
        i  = (i0 - 1) * istp + 1;
        ip = i + id;
        j  = 0;
        j0 = 0;
        v1 = fvol(i ,j ,k);
        v2 = fvol(ip,j ,k);
        ovol = one / (v1 + v2);
        cq(i0,j0) = (fq(i ,j) * v1 + fq(ip,j) * v2) * ovol;

        j  = fnj;
        j0 = cnj;
        v1 = fvol(i ,j ,k);
        v2 = fvol(ip,j ,k);
        ovol = one / (v1 + v2);
        cq(i0,j0) = (fq(i ,j) * v1 + fq(ip,j) * v2) * ovol;
    }
}

template < typename T >
void RestrictQ(StructGrid *fgrid, PHArray< T, 3> &fq, StructGrid *cgrid, PHArray< T, 3> &cq, int nm)
{
    RDouble3D & fvol = * (fgrid->GetCellVolume());
    RDouble3D & cvol = * (cgrid->GetCellVolume());

    int fni = fgrid->GetNI();
    int fnj = fgrid->GetNJ();

    int cni = cgrid->GetNI();
    int cnj = cgrid->GetNJ();

    cq = 0.0;

    int istp = (fni - 1) / (cni - 1);
    int jstp = (fnj - 1) / (cnj - 1);

    int id = istp / 2;
    int jd = jstp / 2;

    int i0,j0,i,j,ip,jp;

    int k = 1;

    RDouble v1,v2,v3,v4,ovol;
    for (j0 = 1; j0 <= cnj - 1; ++ j0)
    {
        j  = (j0 - 1) * jstp + 1;
        jp = j + jd;
        for (i0 = 1; i0 <= cni - 1; ++ i0)
        {
            i  = (i0 - 1) * istp + 1;
            ip = i + id;
            v1 = fvol(i ,j ,k);
            v2 = fvol(ip,j ,k);
            v3 = fvol(ip,jp,k);
            v4 = fvol(i ,jp,k);

            ovol = one / (v1 + v2 + v3 + v4);

            for (int m = 0; m < nm; ++ m)
            {
                cq(i0,j0,m) = (fq(i ,j ,m) * v1 + fq(ip,j ,m) * v2 +
                                fq(ip,jp,m) * v3 + fq(i ,jp,m) * v4) * ovol;
            }
        }
    }
    //! And secondly, restrict the first layer of ghost cells.
    for (j0 = 0; j0 <= cnj; ++ j0)
    {
        j  = (j0 - 1) * jstp + 1;
        jp = j + jd;
        i  = 0;
        i0 = 0;
        v1 = fvol(i ,j ,k);
        v2 = fvol(i ,jp,k);
        ovol = one / (v1 + v2);
        for (int m = 0; m < nm; ++ m)
        {
            cq(i0,j0,m) = (fq(i ,j ,m) * v1 + fq(i ,jp,m) * v2) * ovol;
        }

        i  = fni;
        i0 = cni;
        v1 = fvol(i ,j ,k);
        v2 = fvol(i ,jp,k);
        ovol = one / (v1 + v2);
        for (int m = 0; m < nm; ++ m)
        {
            cq(i0,j0,m) = (fq(i ,j ,m) * v1 + fq(i ,jp,m) * v2) * ovol;
        }
    }

    for (i0 = 0; i0 <= cni; ++ i0)
    {
        i  = (i0 - 1) * istp + 1;
        ip = i + id;
        j  = 0;
        j0 = 0;
        v1 = fvol(i ,j ,k);
        v2 = fvol(ip,j ,k);
        ovol = one / (v1 + v2);
        for (int m = 0; m < nm; ++ m)
        {
            cq(i0,j0,m) = (fq(i ,j ,m) * v1 + fq(ip,j ,m) * v2) * ovol;
        }

        j  = fnj;
        j0 = cnj;
        v1 = fvol(i ,j ,k);
        v2 = fvol(ip,j ,k);
        ovol = one / (v1 + v2);
        for (int m = 0; m < nm; ++ m)
        {
            cq(i0,j0,m) = (fq(i ,j ,m) * v1 + fq(ip,j ,m) * v2) * ovol;
        }
    }
}

template < typename T >
void RestrictQ(StructGrid *fgrid, PHArray< T, 3> &fq, StructGrid *cgrid, PHArray< T, 3> &cq)
{
    RDouble3D & fvol = * (fgrid->GetCellVolume());

    int fni = fgrid->GetNI();
    int fnj = fgrid->GetNJ();
    int fnk = fgrid->GetNK();

    int cni = cgrid->GetNI();
    int cnj = cgrid->GetNJ();
    int cnk = cgrid->GetNK();

    cq = 0.0;

    int istp = fgrid->GetMultigridStepI();
    int jstp = fgrid->GetMultigridStepJ();
    int kstp = fgrid->GetMultigridStepK();

    int id = istp / 2;
    int jd = jstp / 2;
    int kd = kstp / 2;

    Range I0(1,MAX(cni-1,1));
    Range J0(1,MAX(cnj-1,1));
    Range K0(1,MAX(cnk-1,1));

    Range CI0(0,cni);
    Range CJ0(0,cnj);
    Range CK0(0,cnk);

    if (cnk == 1) CK0.setRange(1,1);

    int nDim = fgrid->GetDim();

    int i0,j0,k0,i,j,k,ip,jp,kp;

    RDouble v1,v2,v3,v4,v5,v6,v7,v8,ovol;
    for (k0 = K0.first(); k0 <= K0.last(); ++ k0)
    {
        k  = (k0 - 1) * kstp + 1;
        kp = k + kd;
        for (j0 = J0.first(); j0 <= J0.last(); ++ j0)
        {
            j  = (j0 - 1) * jstp + 1;
            jp = j + jd;
            for (i0 = I0.first(); i0 <= I0.last(); ++ i0)
            {
                i  = (i0 - 1) * istp + 1;
                ip = i + id;
                v1 = fvol(i ,j ,k);
                v2 = fvol(ip,j ,k);
                v3 = fvol(ip,jp,k);
                v4 = fvol(i ,jp,k);
                v5 = fvol(i ,j ,kp);
                v6 = fvol(ip,j ,kp);
                v7 = fvol(ip,jp,kp);
                v8 = fvol(i ,jp,kp);

                ovol = one / (v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8);
                cq(i0,j0,k0) = (fq(i ,j ,k) * v1
                               + fq(ip,j ,k) * v2
                               + fq(ip,jp,k) * v3
                               + fq(i ,jp,k) * v4
                               + fq(i ,j ,kp) * v5
                               + fq(ip,j ,kp) * v6
                               + fq(ip,jp,kp) * v7
                               + fq(i ,jp,kp) * v8
                              ) * ovol;
            }
        }
    }

    //! And secondly, restrict the first layer of ghost cells.
    for (k0 = CK0.first(); k0 <= CK0.last(); ++ k0)
    {
        k  = (k0 - 1) * kstp + 1;
        kp = k + kd;
        for (j0 = CJ0.first(); j0 <= CJ0.last(); ++ j0)
        {
            j  = (j0 - 1) * jstp + 1;
            jp = j + jd;
            i  = 0;
            i0 = 0;
            v1 = fvol(i ,j ,k);
            v2 = fvol(i ,j ,kp);
            v3 = fvol(i ,jp,k);
            v4 = fvol(i ,jp,kp);
            ovol = one / (v1 + v2 + v3 + v4);
            cq(i0,j0,k0) = (fq(i ,j ,k) * v1
                           + fq(i ,j ,kp) * v2
                           + fq(i ,jp,k) * v3
                           + fq(i ,jp,kp) * v4
                          ) * ovol;

            //! Second layer.
            cq(i0-1,j0,k0) = cq(i0,j0,k0);

            i  = fni;
            i0 = cni;
            v1 = fvol(i ,j ,k);
            v2 = fvol(i ,j ,kp);
            v3 = fvol(i ,jp,k);
            v4 = fvol(i ,jp,kp);
            ovol = one / (v1 + v2 + v3 + v4);
            cq(i0,j0,k0) = (fq(i ,j ,k) * v1
                           + fq(i ,j ,kp) * v2
                           + fq(i ,jp,k) * v3
                           + fq(i ,jp,kp) * v4
                          ) * ovol;

            //! Second layer.
            cq(i0+1,j0,k0) = cq(i0,j0,k0);
        }

        for (i0 = CI0.first(); i0 <= CI0.last(); ++ i0)
        {
            i  = (i0 - 1) * istp + 1;
            ip = i + id;
            j  = 0;
            j0 = 0;
            v1 = fvol(i ,j ,k);
            v2 = fvol(ip,j ,k);
            v3 = fvol(i ,j ,kp);
            v4 = fvol(ip,j ,kp);
            ovol = one / (v1 + v2 + v3 + v4);
            cq(i0,j0,k0) = (fq(i ,j ,k) * v1
                           + fq(ip,j ,k) * v2
                           + fq(i ,j ,kp) * v3
                           + fq(ip,j ,kp) * v4
                          ) * ovol;

            //! Second layer.
            cq(i0,j0-1,k0) = cq(i0,j0,k0);

            j  = fnj;
            j0 = cnj;
            v1 = fvol(i ,j ,k);
            v2 = fvol(ip,j ,k);
            v3 = fvol(i ,j ,kp);
            v4 = fvol(ip,j ,kp);
            ovol = one / (v1 + v2 + v3 + v4);
            cq(i0,j0,k0) = (fq(i ,j ,k) * v1
                           + fq(ip,j ,k) * v2
                           + fq(i ,j ,kp) * v3
                           + fq(ip,j ,kp) * v4
                          ) * ovol;

            //! Second layer.
            cq(i0,j0+1,k0) = cq(i0,j0,k0);
        }
    }

    if (nDim == THREE_D)
    {
        for (j0 = CJ0.first(); j0 <= CJ0.last(); ++ j0)
        {
            j  = (j0 - 1) * jstp + 1;
            jp = j + jd;

            for (i0 = CI0.first(); i0 <= CI0.last(); ++ i0)
            {
                i  = (i0 - 1) * istp + 1;
                ip = i + id;
                k  = 0;
                k0 = 0;
       
                v1 = fvol(i ,j ,k);
                v2 = fvol(ip,j ,k);
                v3 = fvol(i ,jp,k);
                v4 = fvol(ip,jp,k);
                ovol = one / (v1 + v2 + v3 + v4);
                cq(i0,j0,k0) = (fq(i ,j ,k) * v1
                               + fq(ip,j ,k) * v2
                               + fq(i ,jp,k) * v3
                               + fq(ip,jp,k) * v4
                              ) * ovol;

                //! Second layer.
                cq(i0,j0,k0-1) = cq(i0,j0,k0);


                k  = fnk;
                k0 = cnk;
       
                v1 = fvol(i ,j ,k);
                v2 = fvol(ip,j ,k);
                v3 = fvol(i ,jp,k);
                v4 = fvol(ip,jp,k);
                ovol = one / (v1 + v2 + v3 + v4);
                cq(i0,j0,k0) = (fq(i ,j ,k) * v1
                               + fq(ip,j ,k) * v2
                               + fq(i ,jp,k) * v3
                               + fq(ip,jp,k) * v4
                              ) * ovol;

                //! Second layer.
                cq(i0,j0,k0+1) = cq(i0,j0,k0);
            }
        }
    }
}

template < typename T >
void RestrictQ(StructGrid *fgrid, PHArray< T, 4> &fq, StructGrid *cgrid, PHArray< T, 4> &cq, int nm)
{
    RDouble3D & fvol = * (fgrid->GetCellVolume());

    int fni = fgrid->GetNI();
    int fnj = fgrid->GetNJ();
    int fnk = fgrid->GetNK();

    int cni = cgrid->GetNI();
    int cnj = cgrid->GetNJ();
    int cnk = cgrid->GetNK();

    Range I0(1,MAX(cni-1,1));
    Range J0(1,MAX(cnj-1,1));
    Range K0(1,MAX(cnk-1,1));

    Range CI0(0,cni);
    Range CJ0(0,cnj);
    Range CK0(0,cnk);

    if (cnk == 1) CK0.setRange(1,1);

    int istp = fgrid->GetMultigridStepI();
    int jstp = fgrid->GetMultigridStepJ();
    int kstp = fgrid->GetMultigridStepK();

    int id = istp / 2;
    int jd = jstp / 2;
    int kd = kstp / 2;

    int i0,j0,k0,i,j,k,ip,jp,kp;
    RDouble v1,v2,v3,v4,v5,v6,v7,v8,ovol;

    int nDim = fgrid->GetDim();

    cq = 0.0;
    for (int m = 0; m < nm; ++ m)
    {
        for (k0 = K0.first(); k0 <= K0.last(); ++ k0)
        {
            k  = (k0 - 1) * kstp + 1;
            kp = k + kd;
            for (j0 = J0.first(); j0 <= J0.last(); ++ j0)
            {
                j  = (j0 - 1) * jstp + 1;
                jp = j + jd;
                for (i0 = I0.first(); i0 <= I0.last(); ++ i0)
                {
                    i  = (i0 - 1) * istp + 1;
                    ip = i + id;
                    v1 = fvol(i ,j ,k);
                    v2 = fvol(ip,j ,k);
                    v3 = fvol(ip,jp,k);
                    v4 = fvol(i ,jp,k);
                    v5 = fvol(i ,j ,kp);
                    v6 = fvol(ip,j ,kp);
                    v7 = fvol(ip,jp,kp);
                    v8 = fvol(i ,jp,kp);

                    ovol = one / (v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8);
                    cq(i0,j0,k0,m) = (fq(i ,j ,k ,m) * v1
                                     + fq(ip,j ,k ,m) * v2
                                     + fq(ip,jp,k ,m) * v3
                                     + fq(i ,jp,k ,m) * v4
                                     + fq(i ,j ,kp,m) * v5
                                     + fq(ip,j ,kp,m) * v6
                                     + fq(ip,jp,kp,m) * v7
                                     + fq(i ,jp,kp,m) * v8
                                    ) * ovol;
                }
            }
        }

        //! And secondly, restrict the first layer of ghost cells.
        for (k0 = CK0.first(); k0 <= CK0.last(); ++ k0)
        {
            k  = (k0 - 1) * kstp + 1;
            kp = k + kd;
            for (j0 = CJ0.first(); j0 <= CJ0.last(); ++ j0)
            {
                j  = (j0 - 1) * jstp + 1;
                jp = j + jd;
                i  = 0;
                i0 = 0;
                v1 = fvol(i ,j ,k);
                v2 = fvol(i ,j ,kp);
                v3 = fvol(i ,jp,k);
                v4 = fvol(i ,jp,kp);
                ovol = one / (v1 + v2 + v3 + v4);
                cq(i0,j0,k0,m) = (fq(i ,j ,k ,m) * v1
                                 + fq(i ,j ,kp,m) * v2
                                 + fq(i ,jp,k ,m) * v3
                                 + fq(i ,jp,kp,m) * v4
                                ) * ovol;
                //! Second layer.
                cq(i0-1,j0,k0,m) = cq(i0,j0,k0,m);

                i  = fni;
                i0 = cni;
                v1 = fvol(i ,j ,k);
                v2 = fvol(i ,j ,kp);
                v3 = fvol(i ,jp,k);
                v4 = fvol(i ,jp,kp);
                ovol = one / (v1 + v2 + v3 + v4);
                cq(i0,j0,k0,m) = (fq(i ,j ,k ,m) * v1
                                 + fq(i ,j ,kp,m) * v2
                                 + fq(i ,jp,k ,m) * v3
                                 + fq(i ,jp,kp,m) * v4
                                ) * ovol;
                //! Second layer.
                cq(i0+1,j0,k0,m) = cq(i0,j0,k0,m);
            }

            for (i0 = CI0.first(); i0 <= CI0.last(); ++ i0)
            {
                i  = (i0 - 1) * istp + 1;
                ip = i + id;
                j  = 0;
                j0 = 0;
                v1 = fvol(i ,j ,k);
                v2 = fvol(ip,j ,k);
                v3 = fvol(i ,j ,kp);
                v4 = fvol(ip,j ,kp);
                ovol = one / (v1 + v2 + v3 + v4);
                cq(i0,j0,k0,m) = (fq(i ,j ,k ,m) * v1
                                 + fq(ip,j ,k ,m) * v2
                                 + fq(i ,j ,kp,m) * v3
                                 + fq(ip,j ,kp,m) * v4
                                ) * ovol;
                //! Second layer.
                cq(i0,j0-1,k0,m) = cq(i0,j0,k0,m);

                j  = fnj;
                j0 = cnj;
                v1 = fvol(i ,j ,k);
                v2 = fvol(ip,j ,k);
                v3 = fvol(i ,j ,kp);
                v4 = fvol(ip,j ,kp);
                ovol = one / (v1 + v2 + v3 + v4);
                cq(i0,j0,k0,m) = (fq(i ,j ,k ,m) * v1
                                 + fq(ip,j ,k ,m) * v2
                                 + fq(i ,j ,kp,m) * v3
                                 + fq(ip,j ,kp,m) * v4
                                ) * ovol;
                //! Second layer.
                cq(i0,j0+1,k0,m) = cq(i0,j0,k0,m);
            }
        }

        if (nDim == TWO_D) continue;

        for (j0 = CJ0.first(); j0 <= CJ0.last(); ++ j0)
        {
            j  = (j0 - 1) * jstp + 1;
            jp = j + jd;

            for (i0 = CI0.first(); i0 <= CI0.last(); ++ i0)
            {
                i  = (i0 - 1) * istp + 1;
                ip = i + id;
                k  = 0;
                k0 = 0;
       
                v1 = fvol(i ,j ,k);
                v2 = fvol(ip,j ,k);
                v3 = fvol(i ,jp,k);
                v4 = fvol(ip,jp,k);
                ovol = one / (v1 + v2 + v3 + v4);
                cq(i0,j0,k0,m) = (fq(i ,j ,k,m) * v1
                                 + fq(ip,j ,k,m) * v2
                                 + fq(i ,jp,k,m) * v3
                                 + fq(ip,jp,k,m) * v4
                                ) * ovol;
                //! Second layer.
                cq(i0,j0,k0-1,m) = cq(i0,j0,k0,m);

                k  = fnk;
                k0 = cnk;
       
                v1 = fvol(i ,j ,k);
                v2 = fvol(ip,j ,k);
                v3 = fvol(i ,jp,k);
                v4 = fvol(ip,jp,k);
                ovol = one / (v1 + v2 + v3 + v4);
                cq(i0,j0,k0,m) = (fq(i ,j ,k,m) * v1
                                 + fq(ip,j ,k,m) * v2
                                 + fq(i ,jp,k,m) * v3
                                 + fq(ip,jp,k,m) * v4
                                ) * ovol;
                //! Second layer.
                cq(i0,j0,k0+1,m) = cq(i0,j0,k0,m);
            }
        }
    }
}

template < typename T >
void RestrictQ(Grid *fgrid_in, T *fq, Grid *cgrid_in, T *cq)
{
    UnstructGrid *fgrid = UnstructGridCast(fgrid_in);
    UnstructGrid *cgrid = UnstructGridCast(cgrid_in);

    int      fnTotalCell  = fgrid->GetNTotalCell();
    int      cnTotalCell  = cgrid->GetNTotalCell();
    int      *cell2coarsegridcell = fgrid->GetCell2CoarseGridCell();
    RDouble *fvol    = fgrid->GetCellVolume();
    RDouble *cvol    = cgrid->GetCellVolume();

    for (int iCell = 0; iCell < cnTotalCell; ++ iCell)
    {
        cq[iCell] = 0.0;
    }

    for (int iCell = 0; iCell < fnTotalCell; ++ iCell)
    {
        cq[ cell2coarsegridcell[iCell] ] += fvol[iCell] * fq[iCell];
    }

    for (int iCell = 0; iCell < cnTotalCell; ++ iCell)
    {
        cq[iCell] /= cvol[iCell];
    }
}

template < typename T >
void InterpolatQ(Grid *fgrid_in, T *fq, Grid *cgrid_in, T *cq)
{
    UnstructGrid *fgrid = UnstructGridCast(fgrid_in);

    int fnTotalCell = fgrid->GetNTotalCell();

    int *cell2coarsegridcell   = fgrid->GetCell2CoarseGridCell();
    int cc;
    for (int iCell = 0; iCell < fnTotalCell; ++ iCell)
    {
        cc = cell2coarsegridcell[iCell];

        fq[iCell] = cq[cc];
    }
}

template < typename T >
void InterpolatQ(StructGrid *fgrid, PHArray< T, 4> &fq, StructGrid *cgrid, PHArray< T, 4> &cq, int nm)
{
    int cni = cgrid->GetNI();
    int cnj = cgrid->GetNJ();
    int cnk = cgrid->GetNK();

    int istp = fgrid->GetMultigridStepI();
    int jstp = fgrid->GetMultigridStepJ();
    int kstp = fgrid->GetMultigridStepK();

    int id = istp / 2;
    int jd = jstp / 2;
    int kd = kstp / 2;

    Range I(1,MAX(cni-1,1));
    Range J(1,MAX(cnj-1,1));
    Range K(1,MAX(cnk-1,1));

    int i0,j0,k0,i,j,k,ip,jp,kp;
    for (k0 = K.first(); k0 <= K.last(); ++ k0)
    {
        k  = (k0 - 1) * kstp + 1;
        kp = k + kd;
        for (j0 = J.first(); j0 <= J.last(); ++ j0)
        {
            j  = (j0 - 1) * jstp + 1;
            jp = j + jd;
            for (i0 = I.first(); i0 <= I.last(); ++ i0)
            {
                i  = (i0 - 1) * istp + 1;
                ip = i + id;

                for (int m = 0; m < nm; ++ m)
                {
                    RDouble temp = cq(i0,j0,k0,m);

                    fq(i ,j ,k ,m) = temp;
                    fq(ip,j ,k ,m) = temp;
                    fq(ip,jp,k ,m) = temp;
                    fq(i ,jp,k ,m) = temp;
                    fq(i ,j ,kp,m) = temp;
                    fq(ip,j ,kp,m) = temp;
                    fq(ip,jp,kp,m) = temp;
                    fq(i ,jp,kp,m) = temp;
                }
            }
        }
    }

    GhostCell3D(cq, cni, cnj, cnk, nm);
}

template < typename T >
void InterpolatQ(StructGrid *fgrid, PHArray< T, 3> &fq, StructGrid *cgrid, PHArray< T, 3> &cq)
{
    int cni = cgrid->GetNI();
    int cnj = cgrid->GetNJ();
    int cnk = cgrid->GetNK();

    int istp = fgrid->GetMultigridStepI();
    int jstp = fgrid->GetMultigridStepJ();
    int kstp = fgrid->GetMultigridStepK();

    int id = istp / 2;
    int jd = jstp / 2;
    int kd = kstp / 2;

    Range I(1,MAX(cni-1,1));
    Range J(1,MAX(cnj-1,1));
    Range K(1,MAX(cnk-1,1));

    int i0,j0,k0,i,j,k,ip,jp,kp;
    for (k0 = K.first(); k0 <= K.last(); ++ k0)
    {
        k  = (k0 - 1) * kstp + 1;
        kp = k + kd;
        for (j0 = J.first(); j0 <= J.last(); ++ j0)
        {
            j  = (j0 - 1) * jstp + 1;
            jp = j + jd;
            for (i0 = I.first(); i0 <= I.last(); ++ i0)
            {
                i  = (i0 - 1) * istp + 1;
                ip = i + id;

                RDouble temp = cq(i0, j0, k0);

                fq(i ,j ,k) = temp;
                fq(ip,j ,k) = temp;
                fq(ip,jp,k) = temp;
                fq(i ,jp,k) = temp;
                fq(i ,j ,kp) = temp;
                fq(ip,j ,kp) = temp;
                fq(ip,jp,kp) = temp;
                fq(i ,jp,kp) = temp;
            }
        }
    }

    GhostCell3D(cq, cni, cnj, cnk);
}

}
