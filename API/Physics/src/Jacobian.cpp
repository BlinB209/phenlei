#include "Jacobian.h"
#include "Constants.h"
#include "Math_BasisFunction.h"
#include "GlobalDataBase.h"

#include "Gas.h"

using namespace std;
namespace PHSPACE
{
//! The fuction computes the formula (M +- R) * dQ, M and R denote the Jacobian matrix and spectrum radius of inviscid flux,respectively,\n
//! dQ is the increment of conservative vector.
void MXDQ_STD(RDouble *prim, RDouble nxs, RDouble nys, RDouble nzs, RDouble area, RDouble vgn, RDouble gama, RDouble *dq, RDouble *f, 
              int nm, int nLaminar, RDouble radius, int ipn)
{
    //! The relevant variables related to vgn has modified by LiPeng on Feb 19, 2019
    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble rm = prim[IDX::IR];
    RDouble um = prim[IDX::IU];
    RDouble vm = prim[IDX::IV];
    RDouble wm = prim[IDX::IW];
    RDouble pm = prim[IDX::IP];

    RDouble c2 = ABS(gama * pm / rm);
    RDouble cm = sqrt(c2);

    RDouble lmd1, lmd2, lmd3, x1, x2, cits;

    //! The variable cits equals to Ur=U-Ub, which can be referred to the formula (A.5) in appendix A of the PHengLEI Theory manual.
    cits = nxs * um + nys * vm + nzs * wm - vgn;

    RDouble ctgm = cits;
    RDouble cmgm = cm;

    lmd1 = ctgm;
    lmd2 = ctgm + cmgm;
    lmd3 = ctgm - cmgm;

    RDouble ipnrad = ipn * radius / MAX(SMALL, area);

    lmd1 = half * (lmd1 + ipnrad);
    lmd2 = half * (lmd2 + ipnrad);
    lmd3 = half * (lmd3 + ipnrad);

    //! This term is actually equal to 0.0 
    x1 = (lmd1 + lmd1 - lmd2 - lmd3) / (c2 + c2);
    //! This term is actually equal to 0.5
    x2 = (lmd2 - lmd3) / (cm + cm);

    RDouble dc, dh, c2dc;
    //! dc=-b1,b1 is referred to the formula (A.7) in appendix A of the PHengLEI Theory manual.
    //! dc   = cits * dq[IR] - nxs * dq[IRU] - nys * dq[IRV] - nzs * dq[IRW];//!the previous form
    dc   = (cits + vgn) * dq[IDX::IR] - nxs * dq[IDX::IRU] - nys * dq[IDX::IRV] - nzs * dq[IDX::IRW];//!Modified by LiPeng on Feb 21, 2019.
    c2dc = c2 * dc;

    //! Compute total enthalpy hm and coefficient dh,dh=b2,see formula (A.7) and (A.8) in appendix A of the PHengLEI Theory manual.
    RDouble hm;
    gas->ComputeTotalEnthalpyAndDH(prim, gama, dq, hm, dh);

    //! This term is actually equal to (0.5 * dc)
    RDouble dhx1_dcx2 =  dh   * x1 + dc * x2;
    //! This term is actually equal to (0.5 * dh)
    RDouble c2dcx1_dhx2 = c2dc * x1 + dh * x2;

    //! The following expressions can be referred to formula (A.6) in appendix A of the PHengLEI Theory manual.
    f[IDX::IR ] = lmd1 * dq[IR ] - dhx1_dcx2;
    f[IDX::IRU] = lmd1 * dq[IRU] - um * dhx1_dcx2 + nxs  * c2dcx1_dhx2;
    f[IDX::IRV] = lmd1 * dq[IRV] - vm * dhx1_dcx2 + nys  * c2dcx1_dhx2;
    f[IDX::IRW] = lmd1 * dq[IRW] - wm * dhx1_dcx2 + nzs  * c2dcx1_dhx2;

    //! The previous form is not suitable for dynamic mesh.
    //! f[IRE] = lmd1 * dq[IRE] - hm * dhx1_dcx2 + cits * c2dcx1_dhx2;
    //! Modified by LiPeng on Feb 21, 2019.
    f[IDX::IRE] = lmd1 * dq[IDX::IRE] - hm * dhx1_dcx2 + (cits + vgn) * c2dcx1_dhx2;
    for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
    {
        f[iLaminar] = lmd1 * dq[iLaminar] - prim[iLaminar] * dhx1_dcx2;
    }

    for (int iLaminar = 0; iLaminar < nLaminar; ++ iLaminar) //!Modified by LiPeng on Feb. 25, 2019
    {
        f[iLaminar] *= area;
    }
}

void MXDQ_STD(RDouble *prim, RDouble *gykb, RDouble *dq, RDouble *f, RDouble gama, int nm, int nLaminar, RDouble radius, int ipn)
{
    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble rm  = prim[IR];
    RDouble um  = prim[IU];
    RDouble vm  = prim[IV];
    RDouble wm  = prim[IW];
    RDouble pm  = prim[IP];

    RDouble c2 = ABS(gama * pm / rm);
    RDouble cm = sqrt(c2);

    RDouble lmd1, lmd2, lmd3, x1, x2, cits;

    //! the variable cits equals to Ur=U-Ub, which can be referred to the formula (A.5) in appendix A of the PHengLEI Theory manual.
    cits = gykb[0] * um + gykb[1] * vm + gykb[2] * wm - gykb[4]; 

    RDouble ctgm = cits;
    RDouble cmgm = cm  ;

    lmd1 = ctgm;
    lmd2 = ctgm + cmgm;
    lmd3 = ctgm - cmgm;

    RDouble ipnrad = ipn * radius / MAX(SMALL, gykb[3]);
    //RDouble ipnrad = ipn * radius;  // org

    lmd1 = half * (lmd1 + ipnrad);
    lmd2 = half * (lmd2 + ipnrad);
    lmd3 = half * (lmd3 + ipnrad);

    x1 = (lmd1 + lmd1 - lmd2 - lmd3) / (c2 + c2);    //this term is actually equal to 0.0 
    x2 = (lmd2 - lmd3) / (cm + cm);                  //this term is actually equal to 0.5

    RDouble dc,dh,c2dc;
    //! dc=-b1, b1 is referred to the formula (A.7) in appendix A of the PHengLEI Theory manual.
    //! the previous form
    dc   = (cits + gykb[4]) * dq[IR] - gykb[0] * dq[IRU] - gykb[1] * dq[IRV] - gykb[2] * dq[IRW];
    c2dc = c2 * dc;

    //! compute total enthalpy hm and coefficient dh,dh=b2,see formula (A.7) and (A.8) in appendix A of the PHengLEI Theory manual.
    RDouble hm;
    gas->ComputeTotalEnthalpyAndDH(prim, gama, dq, hm, dh);

    RDouble dhx1_dcx2   =  dh   * x1 + dc * x2;//this term is actually equal to (0.5 * dc)
    RDouble c2dcx1_dhx2 =  c2dc * x1 + dh * x2;//this term is actually equal to (0.5 * dh)

    //the following expressions can be referred to formula (A.6) in appendix A of the PHengLEI Theory manual.
    f[IR ] = lmd1 * dq[IR ] - dhx1_dcx2;
    f[IRU] = lmd1 * dq[IRU] - um * dhx1_dcx2 + gykb[0]  * c2dcx1_dhx2;
    f[IRV] = lmd1 * dq[IRV] - vm * dhx1_dcx2 + gykb[1]  * c2dcx1_dhx2;
    f[IRW] = lmd1 * dq[IRW] - wm * dhx1_dcx2 + gykb[2]  * c2dcx1_dhx2;
    
    //! the previous form is not suitable for dynamic mesh.
    f[IRE] = lmd1 * dq[IRE] - hm * dhx1_dcx2 + (cits + gykb[4]) * c2dcx1_dhx2;
    for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
    {
        f[iLaminar] = lmd1 * dq[iLaminar] - prim[iLaminar] * dhx1_dcx2;
    }

    for (int iLaminar = 0; iLaminar < nLaminar; ++ iLaminar) //modified by LiPeng on Feb. 25, 2019
    {
        f[iLaminar] *= gykb[3];
    }
}

void ComputeConservativePreconMatrix(int sign, RDouble **MG, RDouble *prim, RDouble gama, RDouble preconCoeff, int nm, int neqn)
{
    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble rm, um, vm, wm, pm, Hm, v2, c2, psi = 1.0, gama1;

    for (int i = 0; i < neqn; ++ i)
    {
        for (int j = 0; j < neqn; ++ j)
        {
            MG[i][j] = 0.0;
        }
    }

    rm = prim[IR];
    um = prim[IU];
    vm = prim[IV];
    wm = prim[IW];
    pm = prim[IP];

    gama1 = gama - 1.0;
    c2 = gama * pm / rm;
    v2 = um * um + vm * vm + wm * wm;

    Hm = gama / gama1 * pm / rm + 0.5 * v2;


    switch(sign) 
    {

    case -1: // M * T-

        psi  = (preconCoeff - 1.0) * gama1 / c2;
        break;

    case 1: // T * M-

        psi  = (1.0 - preconCoeff) / preconCoeff * gama1 / c2 ;
        break;
    }

    MG[IR][0] = 1.0 + 0.5 * v2 * psi;
    MG[IR][1] = - um * psi;
    MG[IR][2] = - vm * psi;
    MG[IR][3] = - wm * psi;
    MG[IR][4] = psi;

    MG[IRU][0] = 0.5 * um * v2 * psi;
    MG[IRU][1] = 1.0 - um * um * psi;
    MG[IRU][2] = -um * vm * psi;
    MG[IRU][3] = -um * wm * psi;
    MG[IRU][4] = um * psi;

    MG[IRV][0] = 0.5 * vm * v2 * psi;
    MG[IRV][1] = -vm * um * psi;
    MG[IRV][2] = 1.0 - vm * vm * psi;
    MG[IRV][3] = -vm * wm * psi;
    MG[IRV][4] = vm * psi;

    MG[IRW][0] = 0.5 * wm * v2 * psi;
    MG[IRW][1] = -wm * um * psi;
    MG[IRW][2] = -wm * vm * psi;
    MG[IRW][3] = 1.0 - wm * wm * psi;
    MG[IRW][4] = wm * psi;

    MG[IRE][0] = 0.5 * Hm * v2 * psi;
    MG[IRE][1] = - um * Hm * psi;
    MG[IRE][2] = - vm * Hm * psi;
    MG[IRE][3] = - wm * Hm * psi;
    MG[IRE][4] = 1.0 + Hm * psi;

    for (int i = nm; i < neqn; ++ i)
    {
        for (int j = 0; j < neqn; ++ j)
        {
            MG[i][j] = prim[i] * MG[IR][j];
            if (j == i)
            {
                MG[i][j] = 1.0;
            }
        }
    }

}

void ConservativePreconMatrixVariables(RDouble *field, RDouble *preconMatrix, int neqn)
{
    using namespace IDX;

    RDouble *field_tmp = new RDouble [neqn];

    for (int m = 0; m < neqn; ++ m)
    {
        field_tmp[m]  = 0;
    }

    for (int i = 0; i < neqn; ++ i)
    {
        for (int j = 0; j < neqn; ++ j)
        {
            int index = i * neqn + j;

            field_tmp[i] += field[j] * preconMatrix[index];
        }
    }

    for (int m = 0; m < neqn; ++ m)
    {
        field[m] = field_tmp[m];
    }

    delete [] field_tmp;    field_tmp = nullptr;
}

void ConservativePreconMatrixVariables(RDouble *field, RDouble *preconMatrix, RDouble timeCoeff, int neqn)
{
    using namespace IDX;

    RDouble *field_tmp = new RDouble [neqn];

    for (int m = 0; m < neqn; ++ m)
    {
        field_tmp[m]  = 0;
    }

    for (int i = 0; i < neqn; ++ i)
    {
        for (int j = 0; j < neqn; ++ j)
        {
            int index = i * neqn + j;

            field_tmp[i] += field[j] * preconMatrix[index] * timeCoeff;
        }
    }

    for (int m = 0; m < neqn; ++ m)
    {
        field[m] = field_tmp[m];
    }

    delete [] field_tmp;    field_tmp = nullptr;
}




//! The fuction added by LiPeng computes the formula (Mv +- Rv) * dQ, Mv and Rv denote the Jacobian matrix and spectrum radius of viscous flux,respectively,\n
//! dQ is the increment of conservative vector.
void MVXDQ_STD_LP(RDouble *prim, RDouble nxs, RDouble nys, RDouble nzs, RDouble area, RDouble vol, RDouble visl, RDouble vist, RDouble gama, 
                  RDouble *dq, RDouble *f, int nm, int nLaminar, RDouble radius, int ipn)
{
    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble rm = prim[IR];//!the density
    RDouble um = prim[IU];//!the velocity u in x-direction
    RDouble vm = prim[IV];//!the velocity v in y-direction
    RDouble wm = prim[IW];//!the velocity w in z-direction
    RDouble UU = um * um + vm * vm + wm * wm;
    RDouble Un = um * nxs + vm * nys + wm * nzs;
    RDouble Uq = um * dq[1] + vm * dq[2] + wm * dq[3];

    RDouble prl = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prt = GlobalDataBase::GetDoubleParaFromDB("prt");

    RDouble scl = GlobalDataBase::GetDoubleParaFromDB("sc_l");
    RDouble sct = GlobalDataBase::GetDoubleParaFromDB("sc_t");

    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber"); //! Obtain the reynold number.
    RDouble coef = area * area / (vol * refReNumber);
    RDouble bb = (visl + vist)/rm;
    RDouble b1 = 1.0/3.0 * bb * (-Un * dq[0] + nxs * dq[1] + nys * dq[2] + nzs * dq[3]);
    RDouble phi = gama * (visl/prl + vist/prt) / rm;

    RDouble cp = 0.0, Tm = 0.0; //!the constant pressure specific heat and temperature.
    gas->GetTemperatureAndCP(prim, Tm, cp);

    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int ns = nLaminar + nChemical - nm;//!the number of species
    RDouble Psi = -cp * (visl/prl + vist/prt) * Tm/rm;
    RDouble kn = 0.0;//!the kinetic energy of turbulence
    Psi += -2.0*bb * UU + phi * 0.5 * UU - bb * kn;

    RDouble *Ds = new RDouble[ns]();
    RDouble *Phi_s = new RDouble[ns]();
    RDouble *theta_s = new RDouble[ns]();
    RDouble fsum = 0.0;

    if (nChemical == 1) //!the relevant terms of chemical reaction
    {
        RDouble K = visl/prl + vist/prt;
        RDouble rD = visl/scl + vist/sct;

        RDouble *hs = new RDouble[ns]();
        gas->GetEverySpeciesEnthalpy(Tm, Tm, Tm, hs);//!Obtain the enthalpy of each species.

        //!Compute the coefficient theta_s, Phi_s and Ds.
        gas->ComputeCoefficientInMVXDQ(prim, K, rD, theta_s, Phi_s, Ds);

        Psi += Ds[ns - 1] * prim[nm + ns - 1] * hs[ns - 1];
        for (int i = 0; i < ns - 1; i++)
        {
            fsum += (Phi_s[i]/rm - phi * theta_s[i]) * dq[nm + i];
            Psi += Ds[i] * prim[nm + i] * hs[i];
            Psi += phi * prim[nm + i] * theta_s[i];
            Psi += -prim[nm + i] * Phi_s[i]/rm;
        }
        delete [] hs;    hs= NULL;
    }

    //!To compute M*dQ
    f[IR ] = 0.0; //!Simplification
    f[IRU] = bb  * dq[IRU] + b1 * nxs - dq[0] * bb * um;
    f[IRV] = bb  * dq[IRV] + b1 * nys - dq[0] * bb * vm;
    f[IRW] = bb  * dq[IRW] + b1 * nzs - dq[0] * bb * wm;
    f[IRE] = phi * dq[IRE] + Uq * (2.0 * bb - phi) + b1 * (2.0*Un) + dq[0] * Psi + fsum; //!fsum denotes the computing result of species terms.
    for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
    {
        //!f[i] = Ds[i - nm] * dq[i] + b1 * 0.0 + dq[0] * 0.0;
        f[iLaminar] = Ds[iLaminar - nm] * dq[iLaminar]; //!Simplification
    }

    for (int iLaminar = 0; iLaminar < nLaminar; ++ iLaminar)
    {
        f[iLaminar] *= 0.5 * coef;
    }

    //!To add rM * dQ that is the product of the viscous spectrum radius rM and increment dQ.
    RDouble ipnrad = 0.5 * ipn * radius;
    for (int iLaminar = 0; iLaminar < nLaminar; ++ iLaminar)
    {
        f[iLaminar] += ipnrad * dq[iLaminar];
    }

    //!Remove the dynamic memory space.
    delete [] Ds;    Ds = NULL;
    delete [] Phi_s;    Ds = NULL;
    delete [] theta_s;    theta_s = NULL;
}

void MatrixFreeDF(RDouble *q, RDouble *dq, RDouble *nxyz, RDouble gama, RDouble *df, int nm, int neqn, int ipn)
{
    using namespace IDX;
    RDouble *Qtmp = new RDouble[neqn];

    RDouble gamm1 = gama - 1.0;
    RDouble nxs = 0.0, nys = 0.0, nzs = 0.0, vgn = 0.0, ns = 0.0;

    nxs = nxyz[0];
    nys = nxyz[1];
    nzs = nxyz[2];
    ns  = nxyz[3];
    vgn = nxyz[4];

    RDouble rm,um,vm,wm,pm,vn,eig;

    rm = q[IR];
    um = q[IU];
    vm = q[IV];
    wm = q[IW];
    pm = q[IP];

    Qtmp[IR]  = rm ;
    Qtmp[IRU] = rm * um;
    Qtmp[IRV] = rm * vm;
    Qtmp[IRW] = rm * wm;
    Qtmp[IRE] = pm / gamm1 + half * rm * (um * um + vm * vm + wm * wm);

    for (int m = nm; m < neqn; ++ m)
    {
        Qtmp[m] = q[IR] * q[m];
    }

    //! Normal velocity and Eigenvalues.
    vn = nxs * um + nys * vm + nzs * wm - vgn;
    eig   = ABS( vn ) + sqrt( gama * pm / rm );

    //! Need to find out the fluxes on level n.
    df[IR ] = - (Qtmp[IR] * vn);
    df[IRU] = - (Qtmp[IRU] * vn + nxs * pm);
    df[IRV] = - (Qtmp[IRV] * vn + nys * pm);
    df[IRW] = - (Qtmp[IRW] * vn + nzs * pm);
    df[IRE] = - ((Qtmp[IRE] + pm) * vn + vgn * pm);

    for (int m = nm; m < neqn; ++ m)
    {
        df[m] -= ( Qtmp[m] * vn );
    }

    //! Conservative variable on level n+1.
    Qtmp[IR]  += dq[IR ];
    Qtmp[IRU] += dq[IRU];
    Qtmp[IRV] += dq[IRV];
    Qtmp[IRW] += dq[IRW];
    Qtmp[IRE] += dq[IRE];

    for (int m = nm; m < neqn; ++ m)
    {
        Qtmp[m] += dq[m];
    }

    rm = Qtmp[IR];
    um = Qtmp[IRU] / rm;
    vm = Qtmp[IRV] / rm;
    wm = Qtmp[IRW] / rm;
    pm = gamm1 * (Qtmp[IRE] - half * rm * (um * um + vm * vm + wm * wm));

    rm = ABS(rm);
    pm = ABS(pm);

    vn  = nxs * um + nys * vm + nzs * wm - vgn;

    //! Now the flux difference due to dq.
    df[IR ] += (Qtmp[IR] * vn);
    df[IRU] += (Qtmp[IRU] * vn + nxs * pm);
    df[IRV] += (Qtmp[IRV] * vn + nys * pm);
    df[IRW] += (Qtmp[IRW] * vn + nzs * pm);
    df[IRE] += ((Qtmp[IRE] + pm) * vn + vgn * pm);

    for (int m = nm; m < neqn; ++ m)
    {
        df[m] += ( Qtmp[m] * vn );
    }

    //! Subtract eigenvalue terms from the flux difference.
    for (int m = 0; m < neqn; ++ m)
    {
        //df[m] -= eig * dq[m];
        df[m] += ipn * eig * dq[m];
    }

    //! Divide by 2 out of necessity.
    for (int m = 0; m < neqn; ++ m)
    {
        df[m] *= half * ns;
    }

    delete [] Qtmp;    Qtmp = nullptr;
}

void PreconMatrixFreeDF(RDouble *q, RDouble *dq, RDouble *nxyz, RDouble *preconMatrix, RDouble gama, RDouble preconCoeff, RDouble *df, int nm, int neqn, int ipn)
{
    using namespace IDX;
    RDouble *Qtmp = new RDouble[neqn];

    RDouble gamm1 = gama - 1.0;
    RDouble nxs = 0.0, nys = 0.0, nzs = 0.0, vgn = 0.0, ns = 0.0;

    nxs = nxyz[0];
    nys = nxyz[1];
    nzs = nxyz[2];
    ns  = nxyz[3];
    vgn = nxyz[4];

    RDouble rm,um,vm,wm,pm,vn,v2,c2;

    rm = q[IR];
    um = q[IU];
    vm = q[IV];
    wm = q[IW];
    pm = q[IP];

    Qtmp[IR]  = rm ;
    Qtmp[IRU] = rm * um;
    Qtmp[IRV] = rm * vm;
    Qtmp[IRW] = rm * wm;
    Qtmp[IRE] = pm / gamm1 + half * rm * (um * um + vm * vm + wm * wm);

    int isPorousZone = GlobalDataBase::GetIntParaFromDB("isPorousZone");
    if (isPorousZone != NON_POROUS)
    {
        RDouble porosity = GlobalDataBase::GetDoubleParaFromDB("porosity");
        RDouble densitySolid = GlobalDataBase::GetDoubleParaFromDB("densitySolid");
        RDouble cpSolid = GlobalDataBase::GetDoubleParaFromDB("cpSolid");
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
        RDouble refTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble referenceEnergy = refVelocity * refVelocity;
        RDouble tm = gama * refMachNumber * refMachNumber * pm / rm;
        Qtmp[IRE] = Qtmp[IRE] * porosity + densitySolid * cpSolid * tm * refTemperature / referenceEnergy * (1.0 - porosity);
    }

    for (int m = nm; m < neqn; ++ m)
    {
        Qtmp[m] = q[IR] * q[m];
    }

    //! Normal velocity and Eigenvalues.
    vn = nxs * um + nys * vm + nzs * wm - vgn;
    v2 = um * um + vm * vm + wm * wm;
    c2  = gama * pm / rm;

    RDouble vnew,cnew,eig;

    vnew  = half * vn * (1.0 + preconCoeff);
    cnew  = half * sqrt(((1.0 - preconCoeff) * vn) * ((1.0 - preconCoeff) * vn) + 4.0 * preconCoeff * c2);
    eig   = (ABS(vnew) +  cnew );

    //! Need to find out the fluxes on level n.
    df[IR ] = - (Qtmp[IR] * vn);
    df[IRU] = - (Qtmp[IRU] * vn + nxs * pm);
    df[IRV] = - (Qtmp[IRV] * vn + nys * pm);
    df[IRW] = - (Qtmp[IRW] * vn + nzs * pm);
    df[IRE] = - ((Qtmp[IRE] + pm) * vn + vgn * pm);

    for (int m = nm; m < neqn; ++ m)
    {
        df[m] -= ( Qtmp[m] * vn );
    }

    //! Conservative variable on level n+1.
    Qtmp[IR]  += dq[IR ];
    Qtmp[IRU] += dq[IRU];
    Qtmp[IRV] += dq[IRV];
    Qtmp[IRW] += dq[IRW];
    Qtmp[IRE] += dq[IRE];

    for (int m = nm; m < neqn; ++ m)
    {
        Qtmp[m] += dq[m];
    }

    rm = Qtmp[IR];
    um = Qtmp[IRU] / rm;
    vm = Qtmp[IRV] / rm;
    wm = Qtmp[IRW] / rm;
    pm = gamm1 * (Qtmp[IRE] - half * rm * (um * um + vm * vm + wm * wm));

    rm = ABS(rm);
    pm = ABS(pm);

    vn  = nxs * um + nys * vm + nzs * wm - vgn;

    //! Now the flux difference due to dq.
    df[IR ] += (Qtmp[IR] * vn);
    df[IRU] += (Qtmp[IRU] * vn + nxs * pm);
    df[IRV] += (Qtmp[IRV] * vn + nys * pm);
    df[IRW] += (Qtmp[IRW] * vn + nzs * pm);
    df[IRE] += ((Qtmp[IRE] + pm) * vn + vgn * pm);

    for (int m = nm; m < neqn; ++ m)
    {
        df[m] += ( Qtmp[m] * vn );
    }

    ConservativePreconMatrixVariables(df, preconMatrix, neqn);

    //! Subtract eigenvalue terms from the flux difference.

    for (int m = 0; m < neqn; ++ m)
    {
        //df[m] -= eig * dq[m];
        df[m] += ipn * eig * dq[m];
    }

    //! Divide by 2 out of necessity.
    for (int m = 0; m < neqn; ++ m)
    {
        df[m] *= half * ns;
    }

    delete [] Qtmp;    Qtmp = nullptr;
}

void PreconMatrixFreeDF(RDouble *q, RDouble *dq, RDouble *nxyz, RDouble *preconMatrix, RDouble gama, RDouble preconCoeff, RDouble timeCoeff, RDouble *df, int nm, int neqn, int ipn)
{
    using namespace IDX;
    RDouble *Qtmp = new RDouble[neqn];

    RDouble gamm1 = gama - 1.0;
    RDouble nxs = 0.0, nys = 0.0, nzs = 0.0, vgn = 0.0, ns = 0.0;

    nxs = nxyz[0];
    nys = nxyz[1];
    nzs = nxyz[2];
    ns  = nxyz[3];
    vgn = nxyz[4];

    RDouble rm,um,vm,wm,pm,vn,v2,c2;

    rm = q[IR];
    um = q[IU];
    vm = q[IV];
    wm = q[IW];
    pm = q[IP];

    Qtmp[IR]  = rm ;
    Qtmp[IRU] = rm * um;
    Qtmp[IRV] = rm * vm;
    Qtmp[IRW] = rm * wm;
    Qtmp[IRE] = pm / gamm1 + half * rm * (um * um + vm * vm + wm * wm);

    for (int m = nm; m < neqn; ++ m)
    {
        Qtmp[m] = q[IR] * q[m];
    }

    //! Normal velocity and Eigenvalues.
    vn = nxs * um + nys * vm + nzs * wm - vgn;
    v2 = um * um + vm * vm + wm * wm;
    c2 = gama * pm / rm;

    RDouble vnew,cnew,eig;

    vnew = half * timeCoeff * vn * (1.0 + preconCoeff);
    cnew = half * timeCoeff * sqrt(((1.0 - preconCoeff) * vn) * ((1.0 - preconCoeff) * vn) + 4.0 * preconCoeff * c2);
    eig  = (ABS(vnew) +  cnew );

    //! Need to find out the fluxes on level n.
    df[IR ] = - (Qtmp[IR] * vn);
    df[IRU] = - (Qtmp[IRU] * vn + nxs * pm);
    df[IRV] = - (Qtmp[IRV] * vn + nys * pm);
    df[IRW] = - (Qtmp[IRW] * vn + nzs * pm);
    df[IRE] = - ((Qtmp[IRE] + pm) * vn + vgn * pm);

    for (int m = nm; m < neqn; ++ m)
    {
        df[m] -= ( Qtmp[m] * vn );
    }

    //! Conservative variable on level n+1.
    Qtmp[IR]  += dq[IR ];
    Qtmp[IRU] += dq[IRU];
    Qtmp[IRV] += dq[IRV];
    Qtmp[IRW] += dq[IRW];
    Qtmp[IRE] += dq[IRE];

    for (int m = nm; m < neqn; ++ m)
    {
        Qtmp[m] += dq[m];
    }

    rm = Qtmp[IR];
    um = Qtmp[IRU] / rm;
    vm = Qtmp[IRV] / rm;
    wm = Qtmp[IRW] / rm;
    pm = gamm1 * (Qtmp[IRE] - half * rm * (um * um + vm * vm + wm * wm));

    rm = ABS(rm);
    pm = ABS(pm);

    vn = nxs * um + nys * vm + nzs * wm - vgn;

    //! Now the flux difference due to dq.
    df[IR ] += (Qtmp[IR] * vn);
    df[IRU] += (Qtmp[IRU] * vn + nxs * pm);
    df[IRV] += (Qtmp[IRV] * vn + nys * pm);
    df[IRW] += (Qtmp[IRW] * vn + nzs * pm);
    df[IRE] += ((Qtmp[IRE] + pm) * vn + vgn * pm);

    for (int m = nm; m < neqn; ++ m)
    {
        df[m] += ( Qtmp[m] * vn );
    }

    ConservativePreconMatrixVariables(df, preconMatrix, timeCoeff, neqn);

    //! Subtract eigenvalue terms from the flux difference.

    for (int m = 0; m < neqn; ++ m)
    {
        //df[m] -= eig * dq[m];
        df[m] += ipn * eig * dq[m];
    }

    //! Divide by 2 out of necessity.
    for (int m = 0; m < neqn; ++ m)
    {
        df[m] *= half * ns;
    }

    delete [] Qtmp;    Qtmp = nullptr;
}
}
