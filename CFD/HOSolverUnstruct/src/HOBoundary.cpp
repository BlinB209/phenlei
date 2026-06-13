#include "HODefine.h"
#include "HOBoundary.h"
#include "Constants.h"
#include "Math_BasisFunction.h"

using PHSPACE::RDouble;
namespace HOUnstruct
{
void SymmetryBoundary(int iGausspoint, int numberGaussPoint, RDouble xNormal, RDouble yNormal, RDouble zNormal, vector < vector < RDouble > > & conservationQ)
{
    RDouble rhom,um, vm, wm, pm, rhoem, vn;
    RDouble rhoR,uR, vR, wR;

    rhom  = conservationQ[0][iGausspoint];
    um    = conservationQ[1][iGausspoint]/rhom;
    vm    = conservationQ[2][iGausspoint]/rhom;
    wm    = conservationQ[3][iGausspoint]/rhom;
    rhoem = conservationQ[4][iGausspoint];
    pm    = HO_GAMAM1 * (rhoem - 0.5 * rhom * (um * um + vm * vm + wm * wm));

    vn = xNormal * um + yNormal * vm + zNormal * wm;

    int iGaussPointR = iGausspoint + numberGaussPoint;

    rhoR  = rhom;
    uR    = um - two * xNormal * vn;
    vR    = vm - two * yNormal * vn;
    wR    = wm - two * zNormal * vn;

    conservationQ[0][iGaussPointR] = rhoR;
    conservationQ[1][iGaussPointR] = rhoR*uR;
    conservationQ[2][iGaussPointR] = rhoR*vR;
    conservationQ[3][iGaussPointR] = rhoR*wR;
    conservationQ[4][iGaussPointR] = pm/HO_GAMAM1 + 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);
}

void SymmetryBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
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

//! 3D.
//! Viscous adiabatic wall.
void ViscousAdiabaticWallBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble xfv, RDouble yfv, RDouble zfv, RDouble refMachNumber, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }

    primt[IU] = - prims[IU] + two * xfv;
    primt[IV] = - prims[IV] + two * yfv;
    primt[IW] = - prims[IW] + two * zfv;
}

void ViscousIsotropicWallBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble xfv, RDouble yfv, RDouble zfv, RDouble refMachNumber, RDouble tw, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }

    primt[IU] = - prims[IU] + two * xfv;
    primt[IV] = - prims[IV] + two * yfv;
    primt[IW] = - prims[IW] + two * zfv;

    RDouble omav = one;
    RDouble coefficientOfStateEquation = 1.0 / (HO_GAMA * refMachNumber * refMachNumber);

    RDouble tm = (1.0 / coefficientOfStateEquation) * prims[IP] / prims[IR];

    RDouble rg, pg, tg; //g ´ú±íghost

    pg = prims[IP];
    tg = two * tw - tm;
    if (tg <= 0.0) tg = tw;

    rg = pg / (coefficientOfStateEquation * tg * omav);

    primt[IR] = rg;
    primt[IP] = pg;
}

void FarfieldBC(RDouble *prims, RDouble *primInflow, RDouble *primt, RDouble gama0, RDouble gama, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;

    RDouble rin,uin,vin,win,pin,cIN,vni;
    RDouble roo,uoo,voo,woo,poo,coo,vno;
    RDouble lam1,lam2,lam3,gamm1,vei;
    RDouble vtx ,vty ,vtz ,entr ;
    RDouble riemm,riemp;
    RDouble rb,ub,vb,wb,pb,vnb,cb;

    //!inner cells
    rin = prims[IR];
    uin = prims[IU];
    vin = prims[IV];
    win = prims[IW];
    pin = prims[IP];

    //!far field
    roo = primInflow[IR];
    uoo = primInflow[IU];
    voo = primInflow[IV];
    woo = primInflow[IW];
    poo = primInflow[IP];

    vno = nxs * uoo + nys * voo + nzs * woo - vgn;
    vni = nxs * uin + nys * vin + nzs * win - vgn;

    coo = sqrt(ABS(gama0 * poo / roo));
    cIN = sqrt(ABS(gama  * pin / rin));

    gamm1 = gama - one;

    lam1  = vni;
    lam2  = vni - cIN;
    lam3  = vni + cIN;

    vei = sqrt(uin * uin + vin * vin + win * win);

    //supersonic flows.
    if (vei > cIN)
    {
        if (vni >= 0.0)
        {
            for (int m = 0; m < neqn; ++ m)
            {
                primt[m] = prims[m];
            }
        }
        else
        {
            for (int m = 0; m < neqn; ++ m)
            {
                primt[m] = primInflow[m];
            }
        }
        return;
    }

    //! if can reach here, it must be the subsonic flow.
    riemp = vni + 2.0 * cIN / gamm1;
    riemm = vno - 2.0 * coo / gamm1;
    vnb =  half * (riemp + riemm);
    cb  =  0.25 * (riemp - riemm) * gamm1;
    if (vnb >= 0.0)
    {
        // exit
        entr = pin / pow(rin, gama);
        vtx  = uin - nxs * vni;
        vty  = vin - nys * vni;
        vtz  = win - nzs * vni;
    }
    else
    {
        //inlet
        entr = poo / pow(roo, gama);
        vtx  = uoo - nxs * vno;
        vty  = voo - nys * vno;
        vtz  = woo - nzs * vno;
        //cb  =   0.25 * (riemp - riemm) * gamm1;
    }

    rb  = pow((cb * cb / (entr * gama)), static_cast<RDouble>(one) / gamm1);
    ub  = vtx + nxs * vnb;
    vb  = vty + nys * vnb;
    wb  = vtz + nzs * vnb;
    pb  = cb * cb * rb / gama;

    primt[IR] = rb;
    primt[IU] = ub;
    primt[IV] = vb;
    primt[IW] = wb;
    primt[IP] = pb;
}

void FarFieldRiemannInvariants(RDouble *prims, RDouble *primInflow, RDouble *primt, RDouble gama0, RDouble gama, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;

    RDouble rin, uin, vin, win, pin, vnin, cinner, vein;
    RDouble roo, uoo, voo, woo, poo, vnoo, coo;
    RDouble rb , ub , vb , wb , pb , vnb , cb;

    //! inner point variables.
    rin = prims[ IR ];
    uin = prims[ IU ];
    vin = prims[ IV ];
    win = prims[ IW ];
    pin = prims[ IP ];

    //! infinite variables.
    roo = primInflow[ IR ];
    uoo = primInflow[ IU ];
    voo = primInflow[ IV ];
    woo = primInflow[ IW ];
    poo = primInflow[ IP ];

    vnin = nxs * uin + nys * vin + nzs * win - vgn;
    vnoo = nxs * uoo + nys * voo + nzs * woo - vgn;

    cinner = sqrt(ABS(gama  * pin / rin));
    coo    = sqrt(ABS(gama0 * poo / roo));

    vein   = sqrt(uin * uin + vin * vin + win * win);

    if (vein > cinner)
    {
        //! supersonic.
        if (vnin >= 0.0)
        {
            //! supersonic outlet.
            for (int m = 0; m < neqn; ++ m)
            {
                primt[ m ] = prims[ m ];
            }
        }
        else
        {
            //! supersonic inlet.
            for (int m = 0; m < neqn; ++ m)
            {
                primt[ m ] = primInflow[ m ];
            }
        }
    }
    else
    {
        //! subsonic.
        RDouble gama1 = gama - one;
        RDouble riemp; //! Riemann invariant for inner.
        RDouble riemm; //! Riemann invariant for infinite.
        riemp = vnin + 2.0 * cinner / gama1;
        riemm = vnoo - 2.0 * coo    / gama1;

        vnb   = half   * (riemp + riemm);
        cb    = fourth * (riemp - riemm) * gama1;

        double sb;                      //! entropy at the boundary.
        double uref, vref, wref, vnref; //! Outlet refers to inner point, and inlet refers to infinite.
        
        if (vnb > 0.0)
        {
            //! subsonic outlet.
            sb    = pin / pow(rin, gama);
            uref  = uin;
            vref  = vin;
            wref  = win;
            vnref = vnin;
        }
        else
        {
            //! subsonic inlet.
            sb    = poo / pow(roo, gama0);
            uref  = uoo;
            vref  = voo;
            wref  = woo;
            vnref = vnoo;
        }

        rb  = pow((cb * cb / (sb * gama)), 1.0 / gama1);
        ub  = uref + nxs * (vnb - vnref);
        vb  = vref + nys * (vnb - vnref);
        wb  = wref + nzs * (vnb - vnref);
        pb  = cb * cb * rb / gama0; //! The same to "pb = sb * pow(rb, gama0);"

        primt[IR] = rb;
        primt[IU] = ub;
        primt[IV] = vb;
        primt[IW] = wb;
        primt[IP] = pb;
    }
}

void InflowBC(RDouble *prims, RDouble *primt, int nl, int nchem)
{
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }
}

void OutflowBC(RDouble *prims, RDouble *primt, int nl, int nchem)
{
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }
}
}

