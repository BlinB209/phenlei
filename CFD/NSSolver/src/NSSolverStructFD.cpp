#include "NSSolverStructFD.h"
#include "Glb_Dimension.h"
#include "Geo_StructBC.h"
#include "Param_NSSolverStruct.h"
#include "TK_Exit.h"
#include "FieldProxy.h"
#include "Math_BasisFunction.h"

namespace PHSPACE
{

NSSolverStructFD::NSSolverStructFD()
{
}

NSSolverStructFD::~NSSolverStructFD()
{
    DeAllocateGlobalVariables();
}

void NSSolverStructFD::AllocateGlobalVar(Grid *gridIn)
{
    CFDSolver::AllocateGlobalVar(gridIn);

    InitLocalMemory(gridIn);

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int nDim    = GetDim();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K); 
    
    Range I4, J4, K4;
    GetRange(ni, nj, nk, -4, 3, I4, J4, K4);     

    Range M(0, 4);    
    Range T(0, 0);
    Range D(0, nDim - 1);

    RDouble4D *residual = new RDouble4D(I, J, K, M, fortranArray);
    RDouble3D *dt       = new RDouble3D(I, J, K,    fortranArray);
    RDouble3D *gama     = new RDouble3D(I, J, K,    fortranArray);
    RDouble4D *visSpectralRadius = new RDouble4D(I, J, K, D, fortranArray);
    RDouble4D *invSpectralRadius = new RDouble4D(I, J, K, D, fortranArray);
    RDouble4D *diagonal = new RDouble4D(I, J, K, M, fortranArray);

    RDouble4D *q    = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *t    = new RDouble4D(I, J, K, T, fortranArray);
    RDouble4D *q_FD = new RDouble4D(I4,J4,K4,M, fortranArray);
    RDouble4D *t_FD = new RDouble4D(I4,J4,K4,T, fortranArray);
    RDouble3D *rtem = new RDouble3D(I, J, K, fortranArray);

    grid->UpdateDataPtr("q"   , q);    //! Scalar flow field variable(rho/u/v/w/p).
    grid->UpdateDataPtr("q_FD", q_FD);    //! Scalar flow field variable(rho/u/v/w/p), only used for FD, 4 ghost layers.
    grid->UpdateDataPtr("res" , residual);    //! Residual or right-hand side.
    grid->UpdateDataPtr("dt"  , dt);    //! Time step.
    grid->UpdateDataPtr("gama", gama);    //! Ratio of specific heat coefficients at constant pressure and volume.
    grid->UpdateDataPtr("t"   , t);    //! Static temperature.
    grid->UpdateDataPtr("t_FD", t_FD);    //! Temperature field, only used for FD, 4 ghost layers.
    grid->UpdateDataPtr("visSpectralRadius", visSpectralRadius);    //! Viscous spectral radius.
    grid->UpdateDataPtr("invSpectralRadius", invSpectralRadius);    //! Inviscid spectral radius.
    grid->UpdateDataPtr("diagonal", diagonal);    //! Sacrificing space to improve efficiency.
    grid->UpdateDataPtr("rtem", rtem);    //! Pressure factor.

    RDouble4D *ql1 = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *qr1 = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *ql2 = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *qr2 = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *ql3 = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *qr3 = new RDouble4D(I, J, K, M, fortranArray);
    grid->UpdateDataPtr("ql1", ql1);
    grid->UpdateDataPtr("qr1", qr1);
    grid->UpdateDataPtr("ql2", ql2);
    grid->UpdateDataPtr("qr2", qr2);
    grid->UpdateDataPtr("ql3", ql3);
    grid->UpdateDataPtr("qr3", qr3);

    RDouble4D *faceInviscidfluxI = new RDouble4D(I4, J4, K4, M, fortranArray);
    RDouble4D *faceInviscidfluxJ = new RDouble4D(I4, J4, K4, M, fortranArray);
    RDouble4D *faceInviscidfluxK = new RDouble4D(I4, J4, K4, M, fortranArray);
    grid->UpdateDataPtr("faceInviscidfluxI", faceInviscidfluxI);    //! Inviscid flux at face, for I direction.
    grid->UpdateDataPtr("faceInviscidfluxJ", faceInviscidfluxJ);    //! Inviscid flux at face, for J direction.
    grid->UpdateDataPtr("faceInviscidfluxK", faceInviscidfluxK);    //! Inviscid flux at face, for K direction.

    Param_NSSolverStruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D *q_unsteady_n1   = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *q_unsteady_n2   = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *res_unsteady_n1 = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *res_unsteady_n2 = new RDouble4D(I, J, K, M, fortranArray);

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble4D *res_unsteady_tmp = new RDouble4D(I, J, K, M, fortranArray);

        grid->UpdateDataPtr("q_unsteady_n1"  , q_unsteady_n1);
        grid->UpdateDataPtr("q_unsteady_n2"  , q_unsteady_n2);
        grid->UpdateDataPtr("res_unsteady_n1", res_unsteady_n1);
        grid->UpdateDataPtr("res_unsteady_n2", res_unsteady_n2);
        grid->UpdateDataPtr("res_unsteady_tmp", res_unsteady_tmp);

        //! Statistical variables for unsteady simulation.
        RDouble4D *qAverage = new RDouble4D(I, J, K, M, fortranArray);
        *qAverage = 0.0;
        grid->UpdateDataPtr("qAverage", qAverage);

        //! Statistical Reynolds stress for unsteady simulation.
        RDouble4D *tauAverage = new RDouble4D(I, J, K, Range(0, 5), fortranArray);
        *tauAverage = 0.0;
        grid->UpdateDataPtr("tauAverage", tauAverage);

        RDouble4D *q2Average = new RDouble4D(I, J, K, Range(0, 5), fortranArray);
        *q2Average = 0.0;
        grid->UpdateDataPtr("q2Average", q2Average);
    }
    
    int viscousType = parameters->GetViscousType();
    if (viscousType > INVISCID)
    {
        RDouble3D *viscousLaminar = new RDouble3D(I, J, K, fortranArray);
        RDouble3D *viscousTurbulence = new RDouble3D(I, J, K, fortranArray);
        grid->UpdateDataPtr("visl", viscousLaminar);
        grid->UpdateDataPtr("vist", viscousTurbulence);
        *viscousTurbulence = 0.0;

        RDouble4D * faceViscousfluxI = new RDouble4D(I4, J4, K4, M, fortranArray);
        RDouble4D * faceViscousfluxJ = new RDouble4D(I4, J4, K4, M, fortranArray);
        RDouble4D * faceViscousfluxK = new RDouble4D(I4, J4, K4, M, fortranArray);
        grid->UpdateDataPtr("faceViscousfluxI", faceViscousfluxI);    //! Viscous flux at face, for I direction.
        grid->UpdateDataPtr("faceViscousfluxJ", faceViscousfluxJ);    //! Viscous flux at face, for J direction.
        grid->UpdateDataPtr("faceViscousfluxK", faceViscousfluxK);    //! Viscous flux at face, for K direction.
    }

    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    if (iLES == LES_SOLVER)
    {
        RDouble3D *subgridScaleEnergy = new RDouble3D(I, J, K, fortranArray);
        grid->UpdateDataPtr("subgridScaleEnergy", subgridScaleEnergy);
        *subgridScaleEnergy = 0.0;
    }

    Range N(0, 3);

    Range I3, J3, K3;
    GetRange(ni, nj, nk, -3, 2, I3, J3, K3);

    RDouble4D *dqdkxi = new RDouble4D(I3, J3, K3, N, fortranArray);
    RDouble4D *dqdeta = new RDouble4D(I3, J3, K3, N, fortranArray);
    RDouble4D *dqdcta = new RDouble4D(I3, J3, K3, N, fortranArray);
    grid->UpdateDataPtr("dqdkxi", dqdkxi);    //! Partial derivative of q with kxi, at cell-center.
    grid->UpdateDataPtr("dqdeta", dqdeta);    //! Partial derivative of q with eta, at cell-center.
    grid->UpdateDataPtr("dqdcta", dqdcta);    //! Partial derivative of q with cta, at cell-center.

    RDouble4D *dqdx = new RDouble4D(I3, J3, K3, N, fortranArray);
    RDouble4D *dqdy = new RDouble4D(I3, J3, K3, N, fortranArray);
    RDouble4D *dqdz = new RDouble4D(I3, J3, K3, N, fortranArray);
    grid->UpdateDataPtr("gradUVWTCellCenterX", dqdx);
    grid->UpdateDataPtr("gradUVWTCellCenterY", dqdy);
    grid->UpdateDataPtr("gradUVWTCellCenterZ", dqdz);

    if (gradient_method.substr(0,10) == "chain_rule")
    {
        RDouble4D *dqdxnokxi = new RDouble4D(I3, J3, K3, N, fortranArray);
        RDouble4D *dqdynokxi = new RDouble4D(I3, J3, K3, N, fortranArray);
        RDouble4D *dqdznokxi = new RDouble4D(I3, J3, K3, N, fortranArray);
        grid->UpdateDataPtr("dqdxnokxi", dqdxnokxi);    //! dqdxnokxi = dqdx - nx * dqdkxi.
        grid->UpdateDataPtr("dqdynokxi", dqdynokxi);    //! dqdynokxi = dqdy - ny * dqdkxi.
        grid->UpdateDataPtr("dqdznokxi", dqdznokxi);    //! dqdznokxi = dqdz - nz * dqdkxi.
        
        RDouble4D *dqdxnoeta = new RDouble4D(I3, J3, K3, N, fortranArray);
        RDouble4D *dqdynoeta = new RDouble4D(I3, J3, K3, N, fortranArray);
        RDouble4D *dqdznoeta = new RDouble4D(I3, J3, K3, N, fortranArray);
        grid->UpdateDataPtr("dqdxnoeta", dqdxnoeta);    //! dqdxnoeta = dqdx - nx * dqdeta.
        grid->UpdateDataPtr("dqdynoeta", dqdynoeta);    //! dqdynoeta = dqdy - ny * dqdeta.
        grid->UpdateDataPtr("dqdznoeta", dqdznoeta);    //! dqdznoeta = dqdz - nz * dqdeta.
        
        RDouble4D *dqdxnocta = new RDouble4D(I3, J3, K3, N, fortranArray);
        RDouble4D *dqdynocta = new RDouble4D(I3, J3, K3, N, fortranArray);
        RDouble4D *dqdznocta = new RDouble4D(I3, J3, K3, N, fortranArray);
        grid->UpdateDataPtr("dqdxnocta", dqdxnocta);    //! dqdxnocta = dqdx - nx * dqdcta.
        grid->UpdateDataPtr("dqdynocta", dqdynocta);    //! dqdynocta = dqdy - ny * dqdcta.
        grid->UpdateDataPtr("dqdznocta", dqdznocta);    //! dqdznocta = dqdz - nz * dqdcta.
    }
    //! For conservation of flux.
    int systemgridtype = GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (MIXGRID == systemgridtype)
    {
        RDouble4D *dqdx_cc_ruvwpt = new RDouble4D(I, J, K, Range(0,5), fortranArray);
        RDouble4D *dqdy_cc_ruvwpt = new RDouble4D(I, J, K, Range(0,5), fortranArray);
        RDouble4D *dqdz_cc_ruvwpt = new RDouble4D(I, J, K, Range(0,5), fortranArray);
        grid->UpdateDataPtr("dqdx_cc_ruvwpt", dqdx_cc_ruvwpt);    //! dqdx(q represents r/u/v/w/p/t) at cell center.
        grid->UpdateDataPtr("dqdy_cc_ruvwpt", dqdy_cc_ruvwpt);    //! dqdy(q represents r/u/v/w/p/t) at cell center.
        grid->UpdateDataPtr("dqdz_cc_ruvwpt", dqdz_cc_ruvwpt);    //! dqdz(q represents r/u/v/w/p/t) at cell center.
    }

    GetSurfaceCellNumber(gridIn, nWallBC, nMaxSurfaceCell);
    RDouble2D *firstLayerHeight = NULL;
    if (nWallBC > 0)
    {
        firstLayerHeight = new RDouble2D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), fortranArray);
    }
    grid->UpdateDataPtr("firstLayerHeight", firstLayerHeight);
    ComputeFirstLayerGridHeight(gridIn, firstLayerHeight);

    RDouble Twall = parameters->GetWallTemperature();
    if (nWallBC > 0)
    {
        RDouble refTemperature = parameters->GetRefDimensionalTemperature();
        RDouble refT = Twall / refTemperature;
        if (Twall <= EPSILON)
        {
            refT = 1.0; //! Radiation equilibrium temperature wall.
        }

        RDouble2D *temperatureWall = new RDouble2D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), fortranArray);

        //! Initialization with the freestream temperature.
        for (int iWall = 0; iWall < nWallBC; ++ iWall)
        {
            for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
            {
                (*temperatureWall)(iWall, iCell) = refT;
            }
        }
        grid->UpdateDataPtr("surfaceTemperature", temperatureWall);
    }

    int nTotalNumber = 16;
    Int1D *speciesOrder = new Int1D(Range(0, nTotalNumber - 1), fortranArray);
    for (int n = 0; n < nTotalNumber; ++ n)
    {
        (*speciesOrder)(n) = -1;
    }
    grid->UpdateDataPtr("speciesOrder", speciesOrder);

    Range IFACE, JFACE, KFACE;
    GetRange(ni, nj, nk, -2, 2, IFACE, JFACE, KFACE);
    RDouble4D *heatTransferCoeff = new RDouble4D(I, J, K, T, fortranArray);
    grid->UpdateDataPtr("heatTransferCoeff", heatTransferCoeff);    //!heat transfer coefficient.
}

void NSSolverStructFD::DeAllocateGlobalVar(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D *q = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D *q_FD = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));
    RDouble4D *residual = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
    RDouble3D *dt = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));
    RDouble3D *gama = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));
    RDouble4D *t = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
    RDouble4D *t_FD = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t_FD"));
    RDouble4D *visSpectralRadius = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble4D *invSpectralRadius = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble4D *diagonal = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
    RDouble3D *rtem = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("rtem"));

    delete q; q = 0;
    delete q_FD; q_FD = 0;
    delete residual; residual = 0;
    delete dt; dt = 0;
    delete gama; gama = 0;
    delete t; t = 0;
    delete t_FD; t_FD = 0;
    delete visSpectralRadius; visSpectralRadius = 0;
    delete invSpectralRadius; invSpectralRadius = 0;
    delete diagonal; diagonal = 0;
    delete rtem; rtem = 0;
    //Remove data pointer.
    grid->DeleteDataPtr("q");
    grid->DeleteDataPtr("q_FD");
    grid->DeleteDataPtr("res");
    grid->DeleteDataPtr("dt");
    grid->DeleteDataPtr("gama");
    grid->DeleteDataPtr("t");
    grid->DeleteDataPtr("t_FD");
    grid->DeleteDataPtr("visSpectralRadius");
    grid->DeleteDataPtr("invSpectralRadius");
    grid->DeleteDataPtr("diagonal");
    grid->DeleteDataPtr("rtem");

    RDouble4D *ql1 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1"));
    RDouble4D *qr1 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1"));
    RDouble4D *ql2 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2"));
    RDouble4D *qr2 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2"));
    RDouble4D *ql3 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3"));
    RDouble4D *qr3 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3"));

    delete ql1; ql1 = 0;
    delete qr1; qr1 = 0;
    delete ql2; ql2 = 0;
    delete qr2; qr2 = 0;
    delete ql3; ql3 = 0;
    delete qr3; qr3 = 0;
    grid->DeleteDataPtr("ql1");
    grid->DeleteDataPtr("qr1");
    grid->DeleteDataPtr("ql2");
    grid->DeleteDataPtr("qr2");
    grid->DeleteDataPtr("ql3");
    grid->DeleteDataPtr("qr3");

    RDouble4D *faceInviscidfluxI = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxI"));
    RDouble4D *faceInviscidfluxJ = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxJ"));
    RDouble4D *faceInviscidfluxK = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxK"));

    delete faceInviscidfluxI; faceInviscidfluxI = 0;
    delete faceInviscidfluxJ; faceInviscidfluxJ = 0;
    delete faceInviscidfluxK; faceInviscidfluxK = 0;
    grid->DeleteDataPtr("faceInviscidfluxI");
    grid->DeleteDataPtr("faceInviscidfluxJ");
    grid->DeleteDataPtr("faceInviscidfluxK");

    Param_NSSolverStruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D *q_unsteady_n1 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
        RDouble4D *q_unsteady_n2 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));
        RDouble4D *res_unsteady_n1 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n1"));
        RDouble4D *res_unsteady_n2 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n2"));
        RDouble4D *res_unsteady_tmp = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_tmp"));

        delete q_unsteady_n1; q_unsteady_n1 = 0;
        delete q_unsteady_n2; q_unsteady_n2 = 0;
        delete res_unsteady_n1; res_unsteady_n1 = 0;
        delete res_unsteady_n2; res_unsteady_n2 = 0;
        delete res_unsteady_tmp; res_unsteady_tmp = 0;
        grid->DeleteDataPtr("q_unsteady_n1");
        grid->DeleteDataPtr("q_unsteady_n2");
        grid->DeleteDataPtr("res_unsteady_n1");
        grid->DeleteDataPtr("res_unsteady_n2");
        grid->DeleteDataPtr("res_unsteady_tmp");

        RDouble4D *qAverage = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage"));
        RDouble4D *tauAverage = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("tauAverage"));

        delete qAverage; qAverage = 0;
        delete tauAverage; tauAverage = 0;
        delete res_unsteady_n1; res_unsteady_n1 = 0;
        grid->DeleteDataPtr("qAverage");
        grid->DeleteDataPtr("tauAverage");
        grid->DeleteDataPtr("q2Average");
    }

    int viscousType = parameters->GetViscousType();
    if (viscousType)
    {
        RDouble3D *viscousLaminar = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        RDouble3D *viscousTurbulence = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
        RDouble4D *faceViscousfluxI = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxI"));
        RDouble4D *faceViscousfluxJ = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxJ"));
        RDouble4D *faceViscousfluxK = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxK"));

        delete viscousLaminar; viscousLaminar = 0;
        delete viscousTurbulence; viscousTurbulence = 0;
        delete faceViscousfluxI; faceViscousfluxI = 0;
        delete faceViscousfluxJ; faceViscousfluxJ = 0;
        delete faceViscousfluxK; faceViscousfluxK = 0;
        grid->DeleteDataPtr("visl");
        grid->DeleteDataPtr("vist");
        grid->DeleteDataPtr("faceViscousfluxI");
        grid->DeleteDataPtr("faceViscousfluxJ");
        grid->DeleteDataPtr("faceViscousfluxK");
    }

    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    if (iLES)
    {
        RDouble3D *subgridScaleEnergy = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));

        delete subgridScaleEnergy; subgridScaleEnergy = 0;
        grid->DeleteDataPtr("subgridScaleEnergy");
    }

    RDouble4D *dqdkxi = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdkxi"));
    RDouble4D *dqdeta = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdeta"));
    RDouble4D *dqdcta = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdcta"));

    delete dqdkxi; dqdkxi = 0;
    delete dqdeta; dqdeta = 0;
    delete dqdcta; dqdcta = 0;
    grid->DeleteDataPtr("dqdkxi");
    grid->DeleteDataPtr("dqdeta");
    grid->DeleteDataPtr("dqdcta");

    RDouble4D *dqdx = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D *dqdy = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D *dqdz = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    delete dqdx; dqdx = 0;
    delete dqdy; dqdy = 0;
    delete dqdz; dqdz = 0;
    grid->DeleteDataPtr("gradUVWTCellCenterX");
    grid->DeleteDataPtr("gradUVWTCellCenterY");
    grid->DeleteDataPtr("gradUVWTCellCenterZ");

    if (gradient_method.substr(0, 10) == "chain_rule")
    {
        RDouble4D *dqdxnokxi = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdxnokxi"));
        RDouble4D *dqdynokxi = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdynokxi"));
        RDouble4D *dqdznokxi = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdznokxi"));

        delete dqdxnokxi; dqdxnokxi = 0;
        delete dqdynokxi; dqdynokxi = 0;
        delete dqdznokxi; dqdznokxi = 0;
        grid->DeleteDataPtr("dqdxnokxi");
        grid->DeleteDataPtr("dqdynokxi");
        grid->DeleteDataPtr("dqdznokxi");
        
        RDouble4D *dqdxnoeta = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdxnoeta"));
        RDouble4D *dqdynoeta = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdynoeta"));
        RDouble4D *dqdznoeta = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdznoeta"));

        delete dqdxnoeta; dqdxnoeta = 0;
        delete dqdynoeta; dqdynoeta = 0;
        delete dqdznoeta; dqdznoeta = 0;
        grid->DeleteDataPtr("dqdxnoeta");
        grid->DeleteDataPtr("dqdynoeta");
        grid->DeleteDataPtr("dqdznoeta");

        RDouble4D *dqdxnocta = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdxnocta"));
        RDouble4D *dqdynocta = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdynocta"));
        RDouble4D *dqdznocta = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdznocta"));

        delete dqdxnocta; dqdxnocta = 0;
        delete dqdynocta; dqdynocta = 0;
        delete dqdznocta; dqdznocta = 0;
        grid->DeleteDataPtr("dqdxnocta");
        grid->DeleteDataPtr("dqdynocta");
        grid->DeleteDataPtr("dqdznocta");
    }

    int systemgridtype = GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (MIXGRID == systemgridtype)
    {
        RDouble4D *dqdx_cc_ruvwpt = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdx_cc_ruvwpt"));
        RDouble4D *dqdy_cc_ruvwpt = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdy_cc_ruvwpt"));
        RDouble4D *dqdz_cc_ruvwpt = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdz_cc_ruvwpt"));

        delete dqdx_cc_ruvwpt; dqdx_cc_ruvwpt = 0;
        delete dqdy_cc_ruvwpt; dqdy_cc_ruvwpt = 0;
        delete dqdz_cc_ruvwpt; dqdz_cc_ruvwpt = 0;
        grid->DeleteDataPtr("dqdx_cc_ruvwpt");
        grid->DeleteDataPtr("dqdy_cc_ruvwpt");
        grid->DeleteDataPtr("dqdz_cc_ruvwpt");
    }

    GetSurfaceCellNumber(gridIn, nWallBC, nMaxSurfaceCell);
    if (nWallBC > 0)
    {
        RDouble2D *firstLayerHeight = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("firstLayerHeight"));
        RDouble2D *temperatureWall = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));

        delete firstLayerHeight; firstLayerHeight = 0;
        delete temperatureWall; temperatureWall = 0;
        grid->DeleteDataPtr("firstLayerHeight");
        grid->DeleteDataPtr("surfaceTemperature");
    }

    Int1D *speciesOrder = reinterpret_cast <Int1D *> (grid->GetDataPtr("speciesOrder"));
    RDouble4D *heatTransferCoeff = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("heatTransferCoeff"));

    delete speciesOrder; speciesOrder = 0;
    delete heatTransferCoeff; heatTransferCoeff = 0;
    grid->DeleteDataPtr("speciesOrder");
    grid->DeleteDataPtr("heatTransferCoeff");

    delete UVWTproxy; UVWTproxy = 0;
    delete UVWTfaceIproxy; UVWTfaceIproxy = 0;
    delete UVWTfaceJproxy; UVWTfaceJproxy = 0;
    delete UVWTfaceKproxy; UVWTfaceKproxy = 0;

    delete dqdxFaceProxy; dqdxFaceProxy = 0;
    delete dqdyFaceProxy; dqdyFaceProxy = 0;
    delete dqdzFaceProxy; dqdzFaceProxy = 0;

    if (solvername.substr(0, 4) == "HDCS")
    {
        delete CellfluxProxy; CellfluxProxy = 0;
    }
}

void NSSolverStructFD::GetSurfaceCellNumber(Grid *gridIn, int &nWall, int &nCell)
{
    StructGrid *grid = StructGridCast(gridIn);

    //! Obtain the number of boundary condition regions.
    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();
    nWall = 0;
    nCell = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
        int BCType = structBC->GetBCType();

        if (!IsWall(BCType))
        {
            continue;
        }
        ++ nWall; //! The number of surface regions.

        //! The following are surface boundary codes.
        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        int ni = ied - ist + 1;
        int nj = jed - jst + 1;
        int nk = ked - kst + 1;
        int nCurCell = MAX(MAX(ni * nk, nj * nk), ni * nj);
        nCell = MAX(nCurCell, nCell);
    }
}

void NSSolverStructFD::ComputeFirstLayerGridHeight(Grid *gridIn, RDouble2D *firstLayerHeight)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());

    using namespace IDX;

    StructBCSet *structBCSet = grid->GetStructBCSet();
    //! Obtain the number of boundary condition regions.
    int nBCRegion = structBCSet->GetnBCRegion();
    int indexOfWall = 0, indexOfCell = 0;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
        int BCType = structBC->GetBCType();

        if (!IsWall(BCType))
        {
            continue;
        }
        //! The following are surface boundary codes.
        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        indexOfCell = 0;
        int nSurface = structBC->GetFaceDirection() + 1;
        int leftOrRightIndex = structBC->GetFaceLeftOrRightIndex();
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int iWall, jWall, kWall;
                    structBC->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);

                    //! Obtain the centroid of the surface.
                    RDouble xCenterWall, yCenterWall, zCenterWall;
                    grid->FaceCoor(iWall, jWall, kWall, nSurface, xCenterWall, yCenterWall, zCenterWall);
                    //! Obtain the centroid of the cell on the positive first layer.
                    RDouble xCellCenter, yCellCenter, zCellCenter;
                    grid->CenterCoor(i, j, k, xCellCenter, yCellCenter, zCellCenter);

                    RDouble normalComponetX = leftOrRightIndex * faceNormalComponentX(iWall, jWall, kWall, nSurface);
                    RDouble normalComponetY = leftOrRightIndex * faceNormalComponentY(iWall, jWall, kWall, nSurface);
                    RDouble normalComponetZ = leftOrRightIndex * faceNormalComponentZ(iWall, jWall, kWall, nSurface);

                    #ifdef USE_ALTERNATIVE_CODE
                    RDouble cellVolume = gridCellVolume(iWall, jWall, kWall); //The volume of first cell.
                    RDouble deltaHeight = half * cellVolume / surfaceArea;
                    #else
                    //! The projection of the dy called the distance between the center of the cell and that of the surface on normal vector.
                    RDouble deltaHeight = fabs((xCellCenter - xCenterWall) * normalComponetX + (yCellCenter - yCenterWall) * normalComponetY + (zCellCenter - zCenterWall) * normalComponetZ);
                    #endif

                    (*firstLayerHeight)(indexOfWall, indexOfCell) = deltaHeight;
                    ++ indexOfCell; //! Next grid cell.
                }
            }
        }
        ++ indexOfWall; //! Next surface regions.
    }
}


void NSSolverStructFD::Primitive2Conservative(RDouble prim[5], RDouble q[5])
{
    using namespace IDX;

    const RDouble beta = 1.4 - 1.0;

    RDouble density = prim[IR];
    RDouble um = prim[IU];
    RDouble vm = prim[IV];
    RDouble wm = prim[IW];
    RDouble pressure = prim[IP];

    RDouble V2 = um * um + vm * vm + wm * wm;

    RDouble totalEnergy = pressure / (beta * density) + half * V2;

    q[IR ] = density;
    q[IRU] = density * um;
    q[IRV] = density * vm;
    q[IRW] = density * wm;
    q[IRE] = density * totalEnergy;    
}

void NSSolverStructFD::Conservative2Primitive(RDouble q[5], RDouble prim[5])
{
    using namespace IDX;

    const RDouble beta = 1.4 - 1.0;

    RDouble density = q[IR];
    if (density <= zero)
    {
        prim[IR] = density;
        return;
    }

    RDouble um   = q[IU] / density;
    RDouble vm   = q[IV] / density;
    RDouble wm   = q[IW] / density;
    RDouble V2   = um * um + vm * vm + wm * wm;

    prim[IR] = density;
    prim[IU] = um;
    prim[IV] = vm;
    prim[IW] = wm;        
    prim[IP] = beta * (q[IP] - half * density * V2);
}

void NSSolverStructFD::GetIndependentVariablesforStructHighOrder(Grid *gridIn)
{
    StructGrid * grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D & q    = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));
    RDouble4D & qFD = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q_FD"));
    
    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

    for (int m = 0; m <= 4; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    qFD(i, j, k, m) = q(i, j, k, m);
                }
            }
        }
    }

    ist = 1;
    jst = 1;
    kst = 1;
    ied = ni - 1;
    jed = nj - 1;
    ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    for (int m = 0; m <= 4; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                qFD(  - 2, j, k, m) = qFD(  - 1, j, k, m);
                qFD(  - 3, j, k, m) = qFD(  - 2, j, k, m);
                qFD(ni + 2, j, k, m) = qFD(ni + 1, j, k, m);
                qFD(ni + 3, j, k, m) = qFD(ni + 2, j, k, m);
            }
        }
    }

    for (int m = 0; m <= 4; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                qFD(i,    - 2, k, m) = qFD(i,    - 1, k, m);
                qFD(i,    - 3, k, m) = qFD(i,    - 2, k, m);
                qFD(i, nj + 2, k, m) = qFD(i, nj + 1, k, m);
                qFD(i, nj + 3, k, m) = qFD(i, nj + 2, k, m);
            }
        }
    }

    if (nk != 1)
    {
        for (int m = 0; m <= 4; ++ m)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    qFD(i, j,    - 2, m) = qFD(i, j,    - 1, m);
                    qFD(i, j,    - 3, m) = qFD(i, j,    - 2, m);
                    qFD(i, j, nk + 2, m) = qFD(i, j, nk + 1, m);
                    qFD(i, j, nk + 3, m) = qFD(i, j, nk + 2, m);
                }
            }
        }
    }
}

void NSSolverStructFD::GetDependentVariablesforStructHighOrder(Grid *gridIn)
{
    Param_NSSolverStruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    ObtainBoundaryValue(gridIn);
    ObtainGamaAndTemperature(gridIn);
    ObtainPressureFactor(gridIn);
    if (viscousType != INVISCID)
    {
        ObtainViscousCoef(gridIn);
        ObtainGradientCellCenter(gridIn);
    }
}

void NSSolverStructFD::ObtainGamaAndTemperature(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &qFD = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));
    RDouble4D &tFD = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t_FD"));
    RDouble3D &gama = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));
    gama = 1.4;

    Param_NSSolverStruct *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble refMachNumber = parameters->GetRefMachNumber();

    Range I, J, K;
    GetRange(ni, nj, nk, -4, 3, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);
    
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble rs = qFD(i, j, k, 0);
                RDouble ps = qFD(i, j, k, 4);
                tFD(i, j, k, 0) = refGama * refMachNumber * refMachNumber * ps / rs;
            }
        }
    }
}

void NSSolverStructFD::ObtainPressureFactor(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));
    RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("rtem"));
    
    RDouble pressureSensorCoefficient = third;  //! For 3-D

    int istp = 1;
    int jstp = 1;
    int kstp = 1;

    Range I, J, K;
    GetRange(ni, nj, nk, 0, -1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

    if (nk == 1)
    {
        kstp = 0;
        kst = 1;
        ked = 1;
        pressureSensorCoefficient = half;
    }

    using namespace IDX;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)    //! Corrected by clz, 2011-12-20.
            {
                int im = i - istp;
                int jm = j - jstp;
                int km = k - kstp;

                int ip = i + istp;
                int jp = j + jstp;
                int kp = k + kstp;

                RDouble pim = q(im, j , k , IP);
                RDouble pjm = q(i , jm, k , IP);
                RDouble pkm = q(i , j , km, IP);
                                            
                RDouble pip = q(ip, j , k , IP);
                RDouble pjp = q(i , jp, k , IP);
                RDouble pkp = q(i , j , kp, IP);
                                            
                RDouble ppp = q(i , j , k , IP);

                RDouble dpi = ABS((pip - two * ppp + pim) / (pip  + two * ppp + pim));
                RDouble dpj = ABS((pjp - two * ppp + pjm) / (pjp  + two * ppp + pjm));
                RDouble dpk = ABS((pkp - two * ppp + pkm) / (pkp  + two * ppp + pkm));

                rtem(i, j, k) =  pressureSensorCoefficient * (dpi + dpj + dpk);
            }
        }
    }

    GhostCell3D(rtem, ni, nj, nk);
}

void NSSolverStructFD::ObtainViscousCoef(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_NSSolverStruct *parameters = GetControlParameters();

    RDouble3D &visl = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble4D &temperature   = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t_FD"));

    RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble tm = temperature(i, j, k,0);
                visl(i, j, k) = tm * sqrt(tm) * (1.0 + nonDimensionalSutherlandTemperature) / (tm + nonDimensionalSutherlandTemperature);
            }
        }
    }
}

void NSSolverStructFD::ObtainGradientCellCenter(Grid *gridIn)
{
    GetUVWTproxy(gridIn);

    int nDim = GetDim();
    for (int nsurf = 1; nsurf <= nDim; ++nsurf)
    {
        GetUVWTfaceproxy(gridIn, nsurf);
        GetdqdkxietactaCellCenter(gridIn, nsurf);
    }

    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &xfv  = *(strgrid->GetFaceVectorX());
    RDouble4D &yfv  = *(strgrid->GetFaceVectorY());
    RDouble4D &zfv  = *(strgrid->GetFaceVectorZ());
    RDouble3D &vol  = *(strgrid->GetCellVolume());  

    RDouble4D &dqdx = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &dqdy = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &dqdz = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterZ"));  

    RDouble4D &dqdkxi = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdkxi"));
    RDouble4D &dqdeta = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdeta"));
    RDouble4D &dqdcta = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdcta"));

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1;
    int jed = nj - 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    } 

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble kxix = half * (xfv(i, j, k, 1) + xfv(i + 1, j, k, 1)); 
                RDouble kxiy = half * (yfv(i, j, k, 1) + yfv(i + 1, j, k, 1));
                RDouble kxiz = half * (zfv(i, j, k, 1) + zfv(i + 1, j, k, 1));

                RDouble etax = half * (xfv(i, j, k, 2) + xfv(i, j + 1, k, 2)); 
                RDouble etay = half * (yfv(i, j, k, 2) + yfv(i, j + 1, k, 2));
                RDouble etaz = half * (zfv(i, j, k, 2) + zfv(i, j + 1, k, 2));
                for (int m = 0; m <= 3; ++ m)
                {
                    dqdx(i, j, k, m) = kxix * dqdkxi(i, j, k, m) + etax * dqdeta(i, j, k, m);
                    dqdy(i, j, k, m) = kxiy * dqdkxi(i, j, k, m) + etay * dqdeta(i, j, k, m);
                    dqdz(i, j, k, m) = kxiz * dqdkxi(i, j, k, m) + etaz * dqdeta(i, j, k, m);
                }
            }
        }
    }

    if (nDim == 3)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble ctax = half * (xfv(i, j, k, 3) + xfv(i, j, k + 1, 3));
                    RDouble ctay = half * (yfv(i, j, k, 3) + yfv(i, j, k + 1, 3));
                    RDouble ctaz = half * (zfv(i, j, k, 3) + zfv(i, j, k + 1, 3));
                    for (int m = 0; m <= 3; ++ m)
                    {
                        dqdx(i, j, k, m) += ctax * dqdcta(i, j, k, m);
                        dqdy(i, j, k, m) += ctay * dqdcta(i, j, k, m);
                        dqdz(i, j, k, m) += ctaz * dqdcta(i, j, k, m);
                    }
                }
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m <= 3; ++ m)
                {
                    dqdx(i, j, k, m) /= vol(i, j, k);
                    dqdy(i, j, k, m) /= vol(i, j, k);
                    dqdz(i, j, k, m) /= vol(i, j, k);
                }
            }
        }
    }

    GhostCell3D(dqdx, ni, nj, nk, 4);
    GhostCell3D(dqdy, ni, nj, nk, 4);
    GhostCell3D(dqdz, ni, nj, nk, 4);

    if (gradient_method.substr(0,10) != "chain_rule")
    {
        return;
    }

    RDouble4D &dqdxnokxi = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdxnokxi"));
    RDouble4D &dqdynokxi = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdynokxi"));
    RDouble4D &dqdznokxi = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdznokxi"));

    RDouble4D &dqdxnoeta = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdxnoeta"));
    RDouble4D &dqdynoeta = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdynoeta"));
    RDouble4D &dqdznoeta = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdznoeta"));

    RDouble4D &dqdxnocta = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdxnocta"));
    RDouble4D &dqdynocta = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdynocta"));
    RDouble4D &dqdznocta = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdznocta"));

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble kxix = half * (xfv(i, j, k, 1) + xfv(i + 1, j, k, 1)); 
                RDouble kxiy = half * (yfv(i, j, k, 1) + yfv(i + 1, j, k, 1));
                RDouble kxiz = half * (zfv(i, j, k, 1) + zfv(i + 1, j, k, 1));

                RDouble etax = half * (xfv(i, j, k, 2) + xfv(i, j + 1, k, 2)); 
                RDouble etay = half * (yfv(i, j, k, 2) + yfv(i, j + 1, k, 2));
                RDouble etaz = half * (zfv(i, j, k, 2) + zfv(i, j + 1, k, 2));

                for (int m = 0; m <= 3; ++ m)
                {
                    dqdxnokxi(i, j, k, m) = dqdx(i, j, k, m) - kxix * dqdkxi(i, j, k, m) / vol(i, j, k);
                    dqdynokxi(i, j, k, m) = dqdy(i, j, k, m) - kxiy * dqdkxi(i, j, k, m) / vol(i, j, k);
                    dqdznokxi(i, j, k, m) = dqdz(i, j, k, m) - kxiz * dqdkxi(i, j, k, m) / vol(i, j, k);

                    dqdxnoeta(i, j, k, m) = dqdx(i, j, k, m) - etax * dqdeta(i, j, k, m) / vol(i, j, k);
                    dqdynoeta(i, j, k, m) = dqdy(i, j, k, m) - etay * dqdeta(i, j, k, m) / vol(i, j, k);
                    dqdznoeta(i, j, k, m) = dqdz(i, j, k, m) - etaz * dqdeta(i, j, k, m) / vol(i, j, k);
                }
            }
        }
    }

    if (nDim == 3)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble ctax = half * (xfv(i, j, k, 3) + xfv(i, j, k + 1, 3));
                    RDouble ctay = half * (yfv(i, j, k, 3) + yfv(i, j, k + 1, 3));
                    RDouble ctaz = half * (zfv(i, j, k, 3) + zfv(i, j, k + 1, 3));
                    for (int m = 0; m <= 3; ++ m)
                    {
                        dqdxnocta(i, j, k, m) = dqdx(i, j, k, m) - ctax * dqdcta(i, j, k, m) / vol(i, j, k);
                        dqdynocta(i, j, k, m) = dqdy(i, j, k, m) - ctay * dqdcta(i, j, k, m) / vol(i, j, k);
                        dqdznocta(i, j, k, m) = dqdz(i, j, k, m) - ctaz * dqdcta(i, j, k, m) / vol(i, j, k);
                    }
                }
            }
        }
    }
    GhostCell3D(dqdxnokxi, ni, nj, nk, 4);
    GhostCell3D(dqdynokxi, ni, nj, nk, 4);
    GhostCell3D(dqdznokxi, ni, nj, nk, 4);

    GhostCell3D(dqdxnoeta, ni, nj, nk, 4);
    GhostCell3D(dqdynoeta, ni, nj, nk, 4);
    GhostCell3D(dqdznoeta, ni, nj, nk, 4);

    GhostCell3D(dqdxnocta, ni, nj, nk, 4);
    GhostCell3D(dqdynocta, ni, nj, nk, 4);
    GhostCell3D(dqdznocta, ni, nj, nk, 4);
}

void NSSolverStructFD::Boundary(Grid *gridIn)
{
    //! Already solved in Controller::GetDependentVariablesforStructHighOrder().
}

void NSSolverStructFD::ComputeGamaAndTemperature(Grid *gridIn)
{
    //! Already solved in Controller::GetDependentVariablesforStructHighOrder().
}

void NSSolverStructFD::CalPressureFactor(Grid *gridIn)
{
    //! Already solved in Controller::GetDependentVariablesforStructHighOrder().
}

void NSSolverStructFD::ComputeViscousCoeff(Grid *gridIn)
{
    //! Already solved in Controller::GetDependentVariablesforStructHighOrder().
}

void NSSolverStructFD::ComputeViscousCoefficientWithoutChemical(Grid *gridIn)
{
    //! Already solved in Controller::GetDependentVariablesforStructHighOrder().
}

void NSSolverStructFD::ComputeGradientCellCenter(Grid *gridIn)
{
    //! Already solved in Controller::GetDependentVariablesforStructHighOrder().
}

void NSSolverStructFD::ObtainBoundaryValue(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    Param_NSSolverStruct *parameters = GetControlParameters();
    
    int viscousType         = parameters->GetViscousType();
    RDouble wallTemperature = parameters->GetWallTemperature();

    using namespace PHENGLEI;

    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
        int BCType = structBC->GetBCType();

        if (IsInterface(BCType))
        {
            continue;
        }        
        else if (IsWall(BCType))
        {
            if (viscousType == INVISCID)
            {
                SymmetryBC3D(grid, structBC);
            }
            else if (wallTemperature <= 0.0)
            {
                ViscousAdiabaticWall(grid, structBC);
            }
            else
            {
                ViscousIsotropicWall(grid, structBC);
            }
        }
        else if (BCType == SYMMETRY)
        {
            SymmetryBC3D(grid, structBC);
        }
        else if (BCType == FARFIELD)
        {
            FarFieldRiemannInvariants(grid, structBC);
        }
        else if (BCType == INFLOW)
        {
            InFlowBC3D(grid, structBC);
        }
        else if (BCType == OUTFLOW)
        {
            OutFlowBC3D(grid, structBC);
        }
        else if (BCType == POLE || BCType/10 == POLE)
        {
            OutFlowBC3D(grid, structBC);
        }
        else if (BCType == OVERSET)
        {
            OutFlowBC3D(grid, structBC);
        }
        else
        {
            //! In order to avoid illegal BCtype ID.
            ostringstream oss;
            oss << "Error : Illegal BCtype ID " << BCType << endl;
            TK_Exit::ExceptionExit(oss);
        }
    }
    CornerPoint(grid);
}

void NSSolverStructFD::InFlowBC3D(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &primitiveVars = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));

    Param_NSSolverStruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int *s_lr3d = structBC->GetFaceDirectionIndex();
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int ntt = 1; ntt <= 4; ++ ntt)
                {
                    int it = i + s_lr3d[0] * ntt;
                    int jt = j + s_lr3d[1] * ntt;
                    int kt = k + s_lr3d[2] * ntt;
                    for (int m = 0; m < 5; ++ m)
                    {
                        primitiveVars(it, jt, kt, m) = primitiveVarFarfield[m];
                    }
                }
            }
        }
    }
}

void NSSolverStructFD::OutFlowBC3D(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int *s_lr3d = structBC->GetFaceDirectionIndex();    
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int ntt = 1; ntt <= 4; ++ ntt)
                {
                    int it = i + s_lr3d[0] * ntt;
                    int jt = j + s_lr3d[1] * ntt;
                    int kt = k + s_lr3d[2] * ntt;
                    for (int m = 0; m < 5; ++ m)
                    {
                        q(it, jt, kt, m) = q(i, j, k, m);
                    }
                }
            }
        }
    }
}

void NSSolverStructFD::SymmetryBC3D(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &xfn  = *(grid->GetFaceNormalX());
    RDouble4D &yfn  = *(grid->GetFaceNormalY());
    RDouble4D &zfn  = *(grid->GetFaceNormalZ());
    
    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int *s_lr3d = structBC->GetFaceDirectionIndex();
    int iSurface = structBC->GetFaceDirection() + 1;

    int id,jd,kd;
    GetBCFaceIDX(s_lr3d, id, jd, kd);

    using namespace IDX;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int in = i + id;
                int jn = j + jd;
                int kn = k + kd;
                RDouble nxs = xfn(in,jn,kn,iSurface);
                RDouble nys = yfn(in,jn,kn,iSurface);
                RDouble nzs = zfn(in,jn,kn,iSurface);

                for (int ntt = 1; ntt <= 4; ++ ntt)
                {
                    int it = i + s_lr3d[0] * ntt;
                    int jt = j + s_lr3d[1] * ntt;
                    int kt = k + s_lr3d[2] * ntt;
                    
                    int nss = ntt - 1;
                    int is = i - s_lr3d[0] * nss;
                    int js = j - s_lr3d[1] * nss;
                    int ks = k - s_lr3d[2] * nss;
                    RestrictIndex(is, js, ks, ni - 1,nj - 1,nk - 1);

                    RDouble ro = q(is, js, ks, IR);
                    RDouble vx = q(is, js, ks, IU);
                    RDouble vy = q(is, js, ks, IV);
                    RDouble vz = q(is, js, ks, IW);
                    RDouble ps = q(is, js, ks, IP);
                    RDouble vn = nxs * vx + nys * vy + nzs * vz;

                    q(it, jt, kt, IR) = ro;
                    q(it, jt, kt, IU) = vx - 2.0 * nxs * vn;
                    q(it, jt, kt, IV) = vy - 2.0 * nys * vn;
                    q(it, jt, kt, IW) = vz - 2.0 * nzs * vn;
                    q(it, jt, kt, IP) = ps;
                }
            }
        }
    }
}

void NSSolverStructFD::ViscousAdiabaticWall(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int *s_lr3d = structBC->GetFaceDirectionIndex();

    using namespace IDX;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int ntt = 1; ntt <= 4; ++ ntt)
                {
                    int it = i + s_lr3d[0] * ntt;
                    int jt = j + s_lr3d[1] * ntt;
                    int kt = k + s_lr3d[2] * ntt;
                    
                    int nss = ntt - 1;
                    int is = i - s_lr3d[0] * nss;
                    int js = j - s_lr3d[1] * nss;
                    int ks = k - s_lr3d[2] * nss;
                    RestrictIndex(is, js, ks, ni - 1,nj - 1,nk - 1);

                    RDouble ro = q(is, js, ks, IR);
                    RDouble vx = q(is, js, ks, IU);
                    RDouble vy = q(is, js, ks, IV);
                    RDouble vz = q(is, js, ks, IW);
                    RDouble ps = q(is, js, ks, IP);

                    q(it, jt, kt, IR) =   ro;
                    q(it, jt, kt, IU) = - vx;
                    q(it, jt, kt, IV) = - vy;
                    q(it, jt, kt, IW) = - vz;
                    q(it, jt, kt, IP) =   ps;
                }
            }
        }
    }
}

void NSSolverStructFD::ViscousIsotropicWall(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_NSSolverStruct *parameters = GetControlParameters();
    
    RDouble refGama = parameters->GetRefGama();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble wallTemperature = parameters->GetWallTemperature();
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int *s_lr3d = structBC->GetFaceDirectionIndex();

    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));
    
    RDouble temperatureWallNonDimensional = wallTemperature / refDimensionalTemperature; //! non-dimensional temperature.
    RDouble temperaturelimitation = 0.1 * temperatureWallNonDimensional;

    using namespace IDX;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int layer = 1; layer <= 4; ++ layer)
                {
                    int nss = layer - 1;
                    int is = i - s_lr3d[0] * nss;
                    int js = j - s_lr3d[1] * nss;
                    int ks = k - s_lr3d[2] * nss;
                    RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);
                    
                    int ntt = layer;
                    int it = i + s_lr3d[0] * ntt;
                    int jt = j + s_lr3d[1] * ntt;
                    int kt = k + s_lr3d[2] * ntt;                    

                    RDouble rs = q(is, js, ks, 0);
                    RDouble us = q(is, js, ks, 1);
                    RDouble vs = q(is, js, ks, 2);
                    RDouble ws = q(is, js, ks, 3);
                    RDouble ps = q(is, js, ks, 4);
                    RDouble ts = refGama * refMachNumber * refMachNumber * ps / rs;

                    RDouble ut = - us;
                    RDouble vt = - vs;
                    RDouble wt = - ws;
                    RDouble pt =   ps;

                    RDouble tt = 2.0 * temperatureWallNonDimensional - ts;
                    if (tt < temperaturelimitation)
                    {
                        tt = temperaturelimitation;
                    }
                    RDouble rt = refGama * refMachNumber * refMachNumber * pt / tt;
                   
                    q(it, jt, kt, 0) = rt;
                    q(it, jt, kt, 1) = ut;
                    q(it, jt, kt, 2) = vt;
                    q(it, jt, kt, 3) = wt;
                    q(it, jt, kt, 4) = pt;
                }                
            }
        }
    }
}

void NSSolverStructFD::FarFieldRiemannInvariants(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    Param_NSSolverStruct *parameters = GetControlParameters();

    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble refGama = parameters->GetRefGama();

    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));

    RDouble4D &xfn = * (grid->GetFaceNormalX());
    RDouble4D &yfn = * (grid->GetFaceNormalY());
    RDouble4D &zfn = * (grid->GetFaceNormalZ());

    int iSurface = structBC->GetFaceDirection() + 1;
    int lr       = structBC->GetFaceLeftOrRightIndex();
    int *s_lr3d  = structBC->GetFaceDirectionIndex();

    int id, jd, kd;
    GetBCFaceIDX(s_lr3d, id, jd, kd);

    using namespace IDX;
    RDouble roo = primitiveVarFarfield[IR];
    RDouble uoo = primitiveVarFarfield[IU];
    RDouble voo = primitiveVarFarfield[IV];
    RDouble woo = primitiveVarFarfield[IW];
    RDouble poo = primitiveVarFarfield[IP];
         
    int nEquation = GetNumberOfEquations();
    RDouble *prims;  //! Inner point.
    RDouble *primt;  //! Ghost point.
    prims = new RDouble [nEquation];
    primt = new RDouble [nEquation];

    RDouble rin, uin, vin, win, pin, vnin, cinner, vein; //! inner    point variables.
    RDouble rb , ub , vb , wb , pb , vnb , cb;           //! boundary point variables.

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int is = i;
                int js = j;
                int ks = k; 

                //! The needed (nx,ny,nz) is stored in the cell of (in,jn,kn).
                int in = i + id;
                int jn = j + jd;
                int kn = k + kd;
                RDouble nxs = lr * xfn(in, jn, kn, iSurface);
                RDouble nys = lr * yfn(in, jn, kn, iSurface);
                RDouble nzs = lr * zfn(in, jn, kn, iSurface);

                //! The first inner point.
                for (int m = 0; m < nEquation; ++ m)
                {
                    prims[m] = q(is, js, ks, m);
                }
                rin    = prims[IR];
                uin    = prims[IU];
                vin    = prims[IV];
                win    = prims[IW];
                pin    = prims[IP];
                vnin   = nxs * uin + nys * vin + nzs * win; //! Not consider the "vgn".
                cinner = sqrt(ABS(refGama * pin / rin));
                vein   = sqrt(uin * uin + vin * vin + win * win);

                //! Infinite.
                RDouble vnoo = nxs * uoo + nys * voo + nzs * woo;
                RDouble coo  = sqrt(ABS(refGama * poo / roo));

                if (vein > cinner)
                {
                    //! Supersonic.
                    if (vnin >= 0.0)
                    {
                        //! Supersonic outlet.
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            primt[m] = prims[m];
                        }
                    }
                    else
                    {
                        //! Supersonic inlet.
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            primt[m] = primitiveVarFarfield[m];
                        }
                    }
                }
                else
                {
                    //! Subsonic.
                    RDouble gama1 = refGama - 1.0;
                    RDouble riemp;     //! Riemann invariant for inner.
                    RDouble riemm;     //! Riemann invariant for infinite.
                    riemp = vnin + 2.0 * cinner / gama1;
                    riemm = vnoo - 2.0 * coo    / gama1;

                    vnb   = half   * (riemp + riemm);
                    cb    = fourth * (riemp - riemm) * gama1;

                    RDouble sb;                         //! Entropy at the boundary.
                    RDouble uref, vref, wref, vnref;    //! Outlet refers to inner point, and inlet refers to infinite.

                    if (vnb > 0.0)
                    {
                        //! Subsonic outlet.
                        sb    = pin / pow(rin, refGama);
                        uref  = uin;
                        vref  = vin;
                        wref  = win;
                        vnref = vnin;
                    }
                    else
                    {
                        //! Subsonic inlet.
                        sb    = poo / pow(roo, refGama);
                        uref  = uoo;
                        vref  = voo;
                        wref  = woo;
                        vnref = vnoo;
                    }

                    rb = pow((cb * cb / (sb * refGama)), 1.0 / gama1);
                    ub = uref + nxs * (vnb - vnref);
                    vb = vref + nys * (vnb - vnref);
                    wb = wref + nzs * (vnb - vnref);
                    pb = cb * cb * rb / refGama;    //! The same to "pb = sb * pow(rb, refGama);".

                    primt[IR] = rb;
                    primt[IU] = ub;
                    primt[IV] = vb;
                    primt[IW] = wb;
                    primt[IP] = pb; 
                }
                for (int ntt = 1; ntt <= 4; ++ ntt)
                {
                    int it = i + s_lr3d[0] * ntt;
                    int jt = j + s_lr3d[1] * ntt;
                    int kt = k + s_lr3d[2] * ntt;

                    q(it, jt, kt, IR) = primt[IR];
                    q(it, jt, kt, IU) = primt[IU];
                    q(it, jt, kt, IV) = primt[IV];
                    q(it, jt, kt, IW) = primt[IW];
                    q(it, jt, kt, IP) = primt[IP];
                }
            }
        }
    }
    delete [] prims;
    delete [] primt;
}

void NSSolverStructFD::CornerPoint(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int nEquation = GetNumberOfEquations();

    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));

    FillCornerPoint3D(q, ni, nj, nk, nEquation);
}

void NSSolverStructFD::ZeroResiduals(Grid *grid)
{
    RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
    residual = 0.0;
}

void NSSolverStructFD::InitFlowAsRestart()
{
    StructGrid *grid = StructGridCast(GetGrid());
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

    Param_NSSolverStruct *parameters = GetControlParameters();

    int outputStepInterval = 0;
    GlobalDataBase::UpdateData("outnstep", &outputStepInterval, PHINT, 1);

    RDouble physicalTime = 0.0;
    GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);

    RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

    Range M(0, nEquation - 1);
    int mst = M.first();
    int med = M.last();

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    primitiveVariables(i, j, k, m) = primitiveVarFarfield[m];
                }
            }
        }
    }

    GetIndependentVariablesforStructHighOrder(grid);
    GetDependentVariablesforStructHighOrder(grid);
    
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        residual   = 0.0;
        return;
    }

    RDouble4D &primitiveVariablen1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
    RDouble4D &primitiveVariablen2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));

    primitiveVariablen1 = primitiveVariables;
    primitiveVariablen2 = primitiveVariablen1;

    RDouble4D &residual   = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
    RDouble4D &residualn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n1"));
    RDouble4D &residualn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n2"));

    residual   = 0.0;
    residualn1 = 0.0;
    residualn2 = 0.0;

    int nStatisticalStep = 0;
    GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

    //! Store the historical sum of InviscidFlux, ViscousFlux, et al., \n
    //! exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
    RDouble4D &resTmp = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_tmp"));

    resTmp = 0.0;
}

void NSSolverStructFD::InitLocalMemory(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    solvername        = GlobalDataBase::GetStrParaFromDB("str_highorder_solver");
    EPSILON           = GlobalDataBase::GetDoubleParaFromDB("str_highorder_interpolation_epsilon");
    fluxname          = GlobalDataBase::GetStrParaFromDB("str_highorder_flux_name");
    gradient_method   = GlobalDataBase::GetStrParaFromDB("structhighordergradient");
    
    Range I(-1, ni + 1);
    Range J(-1, nj + 1);
    Range K(-1, nk + 1);
    if (nk == 1) K.setRange(1, 1);

    UVWTfaceIproxy = new FieldProxy();
    UVWTfaceJproxy = new FieldProxy();
    UVWTfaceKproxy = new FieldProxy();
    UVWTfaceIproxy->SetField_STR(new RDouble4D(I, J, K, Range(0, 3), fortranArray), true); 
    UVWTfaceJproxy->SetField_STR(new RDouble4D(I, J, K, Range(0, 3), fortranArray), true);
    UVWTfaceKproxy->SetField_STR(new RDouble4D(I, J, K, Range(0, 3), fortranArray), true);

    dqdxFaceProxy = new FieldProxy();
    dqdyFaceProxy = new FieldProxy();
    dqdzFaceProxy = new FieldProxy();
    dqdxFaceProxy->SetField_STR(new RDouble4D(I, J, K, Range(0, 3), fortranArray), true); 
    dqdyFaceProxy->SetField_STR(new RDouble4D(I, J, K, Range(0, 3), fortranArray), true);
    dqdzFaceProxy->SetField_STR(new RDouble4D(I, J, K, Range(0, 3), fortranArray), true);

    Range I4(- 3, ni + 3);
    Range J4(- 3, nj + 3);
    Range K4(- 3, nk + 3);
    if (nk == 1) K4.setRange(1, 1);
    
    UVWTproxy = new FieldProxy();
    UVWTproxy->SetField_STR(new RDouble4D(I4, J4, K4, Range(0, 3), fortranArray), true);
    
    if (solvername.substr(0, 4) == "HDCS")
    {
        CellfluxProxy = new FieldProxy();
        CellfluxProxy->SetField_STR(new RDouble4D(I4, J4, K4, Range(0, 4), fortranArray), true);
    }

    WeightedInterpolationCoef = new RDouble6D(I, J, K, Range(1, 3), Range(0, 1), Range(0, 29), fortranArray);

    SetWeightedInterpolationCoef(gridIn);
}

void NSSolverStructFD::SetWeightedInterpolationCoef(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    string InterpolationType = GlobalDataBase::GetStrParaFromDB("str_highorder_interpolation_type");

    if (InterpolationType.substr(0, 4) == "test")
    {
        //! non-uniform interpolation
        RDouble3D &xcc = *grid->GetCellCenterX();
        RDouble3D &ycc = *grid->GetCellCenterY();
        RDouble3D &zcc = *grid->GetCellCenterZ();
        
        int nDim = GetDim();
        for (int nsurf = 1; nsurf <= nDim; ++ nsurf)
        {
            int il1, jl1, kl1;
            GetNsurfIndex(nsurf, il1, jl1, kl1);
            
            int ist = 1 - il1;
            int jst = 1 - jl1;
            int kst = 1 - kl1;
            int ied = ni - 1 + 2*il1;
            int jed = nj - 1 + 2*jl1;
            int ked = nk - 1 + 2*kl1;
            if (nk == 1)
            {
                kst = 1;
                ked = 1;
            }
            
            RDouble LengthPlate_L[4];
            for (int k = kst; k <= ked; ++k)
            {
                for (int j = jst; j <= jed; ++j)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        for (int ip = 0; ip <= 3; ++ ip)
                        {
                            int iL = i + (ip - 2) * il1;
                            int jL = j + (ip - 2) * jl1;
                            int kL = k + (ip - 2) * kl1;
                            
                            RDouble lx = xcc(iL, jL, kL) - xcc(iL - il1, jL - jl1, kL - kl1);
                            RDouble ly = ycc(iL, jL, kL) - ycc(iL - il1, jL - jl1, kL - kl1);
                            RDouble lz = zcc(iL, jL, kL) - zcc(iL - il1, jL - jl1, kL - kl1);
                            
                            LengthPlate_L[ip] = sqrt(lx * lx + ly * ly + lz * lz);
                        }
                        RDouble Lengthbase = LengthPlate_L[2];
                        RDouble LL = LengthPlate_L[0] / Lengthbase;
                        RDouble L  = LengthPlate_L[1] / Lengthbase;
                        RDouble RR = LengthPlate_L[3] / Lengthbase;

                        RDouble weight[3];
                        weight[0] = (1.0 + 2.0 * RR) / (4.0 * (1.0 + L + LL) * (1.0 + L + LL + RR));
                        weight[2] = (1.0 + 2.0 * L) * (1.0 + 2.0 * L + 2.0 * LL) / (4.0 * (1.0 + L + RR) * (1.0 + L + LL + RR));
                        weight[1] = 1.0 - weight[0] - weight[2];

                        RDouble InterCoef0[3];
                        InterCoef0[0] = (1.0 + 2.0 * L) / (4.0 * LL * (L + LL));
                        InterCoef0[1] = - (1.0 + 2.0 * L + 2.0 * LL) / (4.0 * L * LL);
                        InterCoef0[2] = 1.0 - InterCoef0[0] - InterCoef0[1];
                        
                        RDouble InterCoef1[3];
                        InterCoef1[0] = - 0.25 / (L * (1.0 + L));
                        InterCoef1[1] =   0.25 * (1.0 + 2.0 * L) / L;
                        InterCoef1[2] = 1.0 - InterCoef1[0] - InterCoef1[1];
                        
                        RDouble InterCoef2[3];
                        InterCoef2[0] = 0.25 * (1.0 + 2.0 * RR) / (1.0 + RR);
                        InterCoef2[1] = 0.25 * (1.0 + 2.0 * RR) / RR;
                        InterCoef2[2] = 1.0 - InterCoef2[0] - InterCoef2[1];
                        
                        RDouble Diff1stCoef0[3];
                        Diff1stCoef0[0] = L / (LL * (L + LL));
                        Diff1stCoef0[1] = - (L + LL) / (L * LL);
                        Diff1stCoef0[2] = - Diff1stCoef0[0] - Diff1stCoef0[1];
                        
                        RDouble Diff1stCoef1[3];
                        Diff1stCoef1[0] = - 1.0 / (L * (1.0 + L));
                        Diff1stCoef1[1] = - 1.0 + 1.0 / L;
                        Diff1stCoef1[2] = - Diff1stCoef1[0] - Diff1stCoef1[1];
                        
                        RDouble Diff1stCoef2[3];
                        Diff1stCoef2[0] = - (2.0 + RR) / (1.0 + RR);
                        Diff1stCoef2[1] = 1.0 + 1.0 / RR;
                        Diff1stCoef2[2] = - Diff1stCoef2[0] - Diff1stCoef2[1];
                        
                        RDouble Diff2ndCoef0[3];
                        Diff2ndCoef0[0] =   2.0 / (LL * (L + LL));
                        Diff2ndCoef0[1] = - 2.0 / (L * LL);
                        Diff2ndCoef0[2] = - Diff2ndCoef0[0] - Diff2ndCoef0[1];
                        
                        RDouble Diff2ndCoef1[3];
                        Diff2ndCoef1[0] =   2.0 / (L * (1.0 + L));
                        Diff2ndCoef1[1] = - 2.0 / L;
                        Diff2ndCoef1[2] = - Diff2ndCoef1[0] - Diff2ndCoef1[1];
                        
                        RDouble Diff2ndCoef2[3];
                        Diff2ndCoef2[0] =   2.0 / (1.0 + RR);
                        Diff2ndCoef2[1] = - 2.0 / RR;
                        Diff2ndCoef2[2] = - Diff2ndCoef2[0] - Diff2ndCoef2[1];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0,  0) = weight[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0,  1) = weight[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0,  2) = weight[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0,  3) = InterCoef0[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0,  4) = InterCoef0[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0,  5) = InterCoef0[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0,  6) = InterCoef1[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0,  7) = InterCoef1[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0,  8) = InterCoef1[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0,  9) = InterCoef2[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 10) = InterCoef2[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 11) = InterCoef2[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 12) = Diff1stCoef0[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 13) = Diff1stCoef0[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 14) = Diff1stCoef0[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 15) = Diff1stCoef1[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 16) = Diff1stCoef1[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 17) = Diff1stCoef1[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 18) = Diff1stCoef2[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 19) = Diff1stCoef2[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 20) = Diff1stCoef2[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 21) = Diff2ndCoef0[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 22) = Diff2ndCoef0[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 23) = Diff2ndCoef0[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 24) = Diff2ndCoef1[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 25) = Diff2ndCoef1[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 26) = Diff2ndCoef1[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 27) = Diff2ndCoef2[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 28) = Diff2ndCoef2[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, 29) = Diff2ndCoef2[2];
                    }
                }
            }
            
            RDouble LengthPlate_R[4];
            for (int k = kst; k <= ked; ++k)
            {
                for (int j = jst; j <= jed; ++j)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        for (int ip = 0; ip <= 3; ++ ip)
                        { 
                            int iR = i + (ip - 1) * il1;
                            int jR = j + (ip - 1) * jl1;
                            int kR = k + (ip - 1) * kl1;
                            
                            RDouble lx = xcc(iR, jR, kR) - xcc(iR - il1, jR - jl1, kR - kl1);
                            RDouble ly = ycc(iR, jR, kR) - ycc(iR - il1, jR - jl1, kR - kl1);
                            RDouble lz = zcc(iR, jR, kR) - zcc(iR - il1, jR - jl1, kR - kl1);
                            
                            LengthPlate_R[3 - ip] = sqrt(lx * lx + ly * ly + lz * lz);
                        }
                        RDouble Lengthbase = LengthPlate_R[2];
                        RDouble LL = LengthPlate_R[0] / Lengthbase;
                        RDouble L  = LengthPlate_R[1] / Lengthbase;
                        RDouble RR = LengthPlate_R[3] / Lengthbase;
                        
                        RDouble weight[3];
                        weight[0] = (1.0 + 2.0 * RR) / (4.0 * (1.0 + L + LL) * (1.0 + L + LL + RR));
                        weight[2] = (1.0 + 2.0 * L) * (1.0 + 2.0 * L + 2.0 * LL) / (4.0 * (1.0 + L + RR) * (1.0 + L + LL + RR));
                        weight[1] = 1.0 - weight[0] - weight[2];
                        
                        RDouble InterCoef0[3];
                        InterCoef0[0] = (1.0 + 2.0 * L) / (4.0 * LL * (L + LL));
                        InterCoef0[1] = - (1.0 + 2.0 * L + 2.0 * LL) / (4.0 * L * LL);
                        InterCoef0[2] = 1.0 - InterCoef0[0] - InterCoef0[1];
                        
                        RDouble InterCoef1[3];
                        InterCoef1[0] = - 0.25 / (L * (1.0 + L));
                        InterCoef1[1] =   0.25 * (1.0 + 2.0 * L) / L;
                        InterCoef1[2] = 1.0 - InterCoef1[0] - InterCoef1[1];
                        
                        RDouble InterCoef2[3];
                        InterCoef2[0] = 0.25 * (1.0 + 2.0 * RR) / (1.0 + RR);
                        InterCoef2[1] = 0.25 * (1.0 + 2.0 * RR) / RR;
                        InterCoef2[2] = 1.0 - InterCoef2[0] - InterCoef2[1];
                        
                        RDouble Diff1stCoef0[3];
                        Diff1stCoef0[0] = L / (LL * (L + LL));
                        Diff1stCoef0[1] = - (L + LL) / (L * LL);
                        Diff1stCoef0[2] = - Diff1stCoef0[0] - Diff1stCoef0[1];
                        
                        RDouble Diff1stCoef1[3];
                        Diff1stCoef1[0] = - 1.0 / (L * (1.0 + L));
                        Diff1stCoef1[1] = - 1.0 + 1.0 / L;
                        Diff1stCoef1[2] = - Diff1stCoef1[0] - Diff1stCoef1[1];
                        
                        RDouble Diff1stCoef2[3];
                        Diff1stCoef2[0] = - (2.0 + RR) / (1.0 + RR);
                        Diff1stCoef2[1] = 1.0 + 1.0 / RR;
                        Diff1stCoef2[2] = - Diff1stCoef2[0] - Diff1stCoef2[1];
                        
                        RDouble Diff2ndCoef0[3];
                        Diff2ndCoef0[0] =   2.0 / (LL * (L + LL));
                        Diff2ndCoef0[1] = - 2.0 / (L * LL);
                        Diff2ndCoef0[2] = - Diff2ndCoef0[0] - Diff2ndCoef0[1];
                        
                        RDouble Diff2ndCoef1[3];
                        Diff2ndCoef1[0] =   2.0 / (L * (1.0 + L));
                        Diff2ndCoef1[1] = - 2.0 / L;
                        Diff2ndCoef1[2] = - Diff2ndCoef1[0] - Diff2ndCoef1[1];
                        
                        RDouble Diff2ndCoef2[3];
                        Diff2ndCoef2[0] =   2.0 / (1.0 + RR);
                        Diff2ndCoef2[1] = - 2.0 / RR;
                        Diff2ndCoef2[2] = - Diff2ndCoef2[0] - Diff2ndCoef2[1];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1,  0) = weight[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1,  1) = weight[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1,  2) = weight[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1,  3) = InterCoef0[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1,  4) = InterCoef0[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1,  5) = InterCoef0[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1,  6) = InterCoef1[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1,  7) = InterCoef1[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1,  8) = InterCoef1[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1,  9) = InterCoef2[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 10) = InterCoef2[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 11) = InterCoef2[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 12) = Diff1stCoef0[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 13) = Diff1stCoef0[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 14) = Diff1stCoef0[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 15) = Diff1stCoef1[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 16) = Diff1stCoef1[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 17) = Diff1stCoef1[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 18) = Diff1stCoef2[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 19) = Diff1stCoef2[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 20) = Diff1stCoef2[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 21) = Diff2ndCoef0[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 22) = Diff2ndCoef0[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 23) = Diff2ndCoef0[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 24) = Diff2ndCoef1[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 25) = Diff2ndCoef1[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 26) = Diff2ndCoef1[2];
                        
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 27) = Diff2ndCoef2[0];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 28) = Diff2ndCoef2[1];
                        (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, 29) = Diff2ndCoef2[2];
                    }
                }
            }
        }
    }
    else
    {
        //! uniform interpolation
        RDouble Coef[30];
        
        Coef[0] =  1.0  / 16.0;
        Coef[1] = 10.0  / 16.0;
        Coef[2] =  5.0  / 16.0;
        
        Coef[3] =    3.0 / 8.0;
        Coef[4] = - 10.0 / 8.0;
        Coef[5] =   15.0 / 8.0;
        
        Coef[6] = - 1.0 / 8.0;
        Coef[7] =   6.0 / 8.0;
        Coef[8] =   3.0 / 8.0;
        
        Coef[ 9] =   3.0 / 8.0;
        Coef[10] =   6.0 / 8.0;
        Coef[11] = - 1.0 / 8.0;
        
        Coef[12] =   0.5;
        Coef[13] = - 2.0;
        Coef[14] =   1.5;
        
        Coef[15] = - 0.5;
        Coef[16] =   0.0;
        Coef[17] =   0.5;
        
        Coef[18] = - 1.5;
        Coef[19] =   2.0;
        Coef[20] = - 0.5;
        
        Coef[21] =   1.0;
        Coef[22] = - 2.0;
        Coef[23] =   1.0;
        
        Coef[24] =   1.0;
        Coef[25] = - 2.0;
        Coef[26] =   1.0;
        
        Coef[27] =   1.0;
        Coef[28] = - 2.0;
        Coef[29] =   1.0;
        
        int nDim = GetDim();
        for (int nsurf = 1; nsurf <= nDim; ++ nsurf)
        {
            int il1, jl1, kl1;
            GetNsurfIndex(nsurf, il1, jl1, kl1);
            
            int ist = 1 - il1;
            int jst = 1 - jl1;
            int kst = 1 - kl1;
            int ied = ni - 1 + 2*il1;
            int jed = nj - 1 + 2*jl1;
            int ked = nk - 1 + 2*kl1;
            
            if (nk == 1)
            {
                kst = 1;
                ked = 1;
            }
            for (int k = kst; k <= ked; ++k)
            {
                for (int j = jst; j <= jed; ++j)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        for (int coefid = 0; coefid <= 29; ++coefid)
                        {
                            (*WeightedInterpolationCoef)(i, j, k, nsurf, 0,  coefid) = Coef[coefid];
                            (*WeightedInterpolationCoef)(i, j, k, nsurf, 1,  coefid) = Coef[coefid];
                        }
                    }
                }
            }
        }
    }    
}

void NSSolverStructFD::InviscidFlux(Grid *gridIn)
{
    StructGrid * grid = StructGridCast(gridIn);
        
    total_negative = 0;

    FieldProxy *qlProxy = new FieldProxy();
    FieldProxy *qrProxy = new FieldProxy();
    FieldProxy * fluxProxy = new FieldProxy();

    int nDim = GetDim();
    for (int nsurf = 1; nsurf <= nDim; ++ nsurf)
    {
        if (nsurf == 1)
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1")));
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxI")));
        }
        else if (nsurf == 2)
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2")));
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxJ")));
        }
        else
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3")));
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxK")));
        }

        GetFaceVariable_WCNS(grid, qlProxy, qrProxy, nsurf);
        GetFaceVariableCorrectionAtPhysicalBoundary(grid, qlProxy, qrProxy, nsurf);
        GetFaceInviscidFlux(grid, qlProxy, qrProxy, fluxProxy, nsurf);
        FluxDifference(grid, fluxProxy, nsurf);
    }

    /*if (total_negative > 0)
    {
        std::cout << "warning:    " << total_negative << std::endl;
    }*/

    delete qlProxy;
    delete qrProxy;
    delete fluxProxy;
}

void NSSolverStructFD::GetFaceVariable_WCNS(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf)
{
    StructGrid * grid = StructGridCast(gridIn);
    RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q_FD"));
    RDouble4D & ql = qlProxy->GetField_STR();
    RDouble4D & qr = qrProxy->GetField_STR();

    Int1D *DiffOrder;
    if (nsurf == 1)
    {
        DiffOrder = grid->DiscretePrecisionKXI;
    }
    else if (nsurf == 2)
    {
        DiffOrder = grid->DiscretePrecisionETA;
    }
    else
    {
        DiffOrder = grid->DiscretePrecisionCTA;
    }

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);
    
    int ist = 1 - il1;
    int jst = 1 - jl1;
    int kst = 1 - kl1;
    int ied = ni - 1 + 2*il1;
    int jed = nj - 1 + 2*jl1;
    int ked = nk - 1 + 2*kl1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }   

    RDouble PrimPlate_L[5][5];
    RDouble PrimPlate_R[5][5];

    RDouble PrimHalf_L[5];
    RDouble PrimHalf_R[5];

    RDouble Coef_L[30];
    RDouble Coef_R[30];
    
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {            
            for (int i = ist; i <= ied; ++ i)
            {
                for (int ip = 0; ip <= 4; ++ ip)
                {
                    int iL = i + (ip - 3) * il1;
                    int jL = j + (ip - 3) * jl1;
                    int kL = k + (ip - 3) * kl1;

                    int iR = i + (ip - 2) * il1;
                    int jR = j + (ip - 2) * jl1;
                    int kR = k + (ip - 2) * kl1;

                    for (int m = 0; m <= 4; ++m)
                    {
                        PrimPlate_L[    ip][m] = q(iL, jL, kL, m);
                        PrimPlate_R[4 - ip][m] = q(iR, jR, kR, m);
                    }
                }

                int localcellID = i * il1 + j * jl1 + k * kl1;
                
                if ((*DiffOrder)(localcellID) == 2 || (*DiffOrder)(localcellID - 1) == 2)
                {
                    InterpolationRUVWP3rd(PrimPlate_L, PrimHalf_L);
                    InterpolationRUVWP3rd(PrimPlate_R, PrimHalf_R);
                }
                else
                {
                    for (int ic = 0; ic <= 29; ++ic)
                    {
                        Coef_L[ic] = (*WeightedInterpolationCoef)(i, j, k, nsurf, 0, ic);
                        Coef_R[ic] = (*WeightedInterpolationCoef)(i, j, k, nsurf, 1, ic);
                    }
                    InterpolationRUVWP5th(Coef_L, PrimPlate_L, PrimHalf_L);
                    InterpolationRUVWP5th(Coef_R, PrimPlate_R, PrimHalf_R);
                }                

                for (int m = 0; m <= 4; ++m)
                {
                    ql(i - il1, j - jl1, k - kl1, m) = PrimHalf_L[m];
                    qr(i,       j,       k,       m) = PrimHalf_R[m];
                }
            }
        }
    }
}

void NSSolverStructFD::GetFaceVariableCorrectionAtPhysicalBoundary(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q_FD"));
    RDouble4D & ql = qlProxy->GetField_STR();
    RDouble4D & qr = qrProxy->GetField_STR();
    
    RDouble4D &xfn  = * (grid->GetFaceNormalX());
    RDouble4D &yfn  = * (grid->GetFaceNormalY());
    RDouble4D &zfn  = * (grid->GetFaceNormalZ());

    Param_NSSolverStruct *parameters = GetControlParameters();
    int viscousType         = parameters->GetViscousType();
    RDouble wallTemperature = parameters->GetWallTemperature();
    RDouble refGama = parameters->GetRefGama();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
    RDouble temperatureWallNonDimensional = wallTemperature / refDimensionalTemperature;
    
    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int BCType = bcregion->GetBCType();
        int nsurf_bc = bcregion->GetFaceDirection() + 1;
        if (!IsWall(BCType) || nsurf_bc != nsurf)
        {
            continue;
        }
        int ist, ied, jst, jed, kst, ked;
        bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);
        
        int il1, jl1, kl1;
        GetNsurfIndex(nsurf, il1, jl1, kl1);

        int *s_lr3d = bcregion->GetFaceDirectionIndex();
        int id, jd, kd;
        GetBCFaceIDX(s_lr3d, id, jd, kd);
        
        RDouble q_1st[5];
        RDouble rtmp, utmp, vtmp, wtmp, ptmp;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int m = 0; m <= 4; ++ m)
                    {
                        q_1st[m] = q(i, j, k, m);
                    }
                    int ibc = i + id;
                    int jbc = j + jd;
                    int kbc = k + kd;
                    
                    if (viscousType == INVISCID)
                    {
                        RDouble vx  = q_1st[1];
                        RDouble vy  = q_1st[2];
                        RDouble vz  = q_1st[3];
                        
                        RDouble nxs = xfn(ibc, jbc, kbc, nsurf);
                        RDouble nys = yfn(ibc, jbc, kbc, nsurf);
                        RDouble nzs = zfn(ibc, jbc, kbc, nsurf);
                        
                        RDouble vn  = nxs * vx + nys * vy + nzs * vz;
                        
                        rtmp = q_1st[0];
                        utmp = vx - nxs * vn;
                        vtmp = vy - nys * vn;
                        wtmp = vz - nzs * vn;
                        ptmp = q_1st[4];
                    }
                    else if (wallTemperature <= 0.0)
                    {
                        rtmp = q_1st[0];
                        utmp = 0.0;
                        vtmp = 0.0;
                        wtmp = 0.0;
                        ptmp = q_1st[4];
                    }
                    else
                    {
                        utmp = 0.0;
                        vtmp = 0.0;
                        wtmp = 0.0;
                        ptmp = q_1st[4];
                        rtmp = refGama * refMachNumber * refMachNumber * ptmp / temperatureWallNonDimensional;
                    }
                    
                    qr(ibc, jbc, kbc, 0) = rtmp;
                    qr(ibc, jbc, kbc, 1) = utmp;
                    qr(ibc, jbc, kbc, 2) = vtmp;
                    qr(ibc, jbc, kbc, 3) = wtmp;
                    qr(ibc, jbc, kbc, 4) = ptmp;
                    
                    ql(ibc - il1, jbc - jl1, kbc - kl1, 0) = rtmp;
                    ql(ibc - il1, jbc - jl1, kbc - kl1, 1) = utmp;
                    ql(ibc - il1, jbc - jl1, kbc - kl1, 2) = vtmp;
                    ql(ibc - il1, jbc - jl1, kbc - kl1, 3) = wtmp;
                    ql(ibc - il1, jbc - jl1, kbc - kl1, 4) = ptmp;
                }
            }
        } 
    }    
}

void NSSolverStructFD::InterpolationRUVWP5th(RDouble Coef[30], RDouble PrimPlate[5][5], RDouble Primhalf[5])
{
    RDouble rBase = PrimPlate[2][0];
    RDouble pBase = PrimPlate[2][4];
    RDouble VBase = sqrt(pBase / rBase);

    RDouble rPlate[5];
    RDouble uPlate[5];
    RDouble vPlate[5];
    RDouble wPlate[5];
    RDouble pPlate[5];

    for (int j = 0; j <= 4; ++ j)
    {
        rPlate[j] = PrimPlate[j][0] / rBase;
        uPlate[j] = PrimPlate[j][1] / VBase;
        vPlate[j] = PrimPlate[j][2] / VBase;
        wPlate[j] = PrimPlate[j][3] / VBase;
        pPlate[j] = PrimPlate[j][4] / pBase;
    }

    Primhalf[0] = rBase * BasicWeightedInterpolation5th(Coef, rPlate);
    Primhalf[4] = pBase * BasicWeightedInterpolation5th(Coef, pPlate);

    if (Primhalf[0] < SMALL || Primhalf[4] < SMALL)
    {
        Primhalf[0] = PrimPlate[2][0];
        Primhalf[1] = PrimPlate[2][1];
        Primhalf[2] = PrimPlate[2][2];
        Primhalf[3] = PrimPlate[2][3];
        Primhalf[4] = PrimPlate[2][4];
        total_negative += 1;
    }
    else
    {
        Primhalf[1] = VBase * BasicWeightedInterpolation5th(Coef, uPlate);
        Primhalf[2] = VBase * BasicWeightedInterpolation5th(Coef, vPlate);
        Primhalf[3] = VBase * BasicWeightedInterpolation5th(Coef, wPlate);
    }
}

void NSSolverStructFD::InterpolationRUVWP3rd(RDouble Prim[5][5], RDouble Primhalf[5])
{
    RDouble Qk[5];

    for (int m = 0; m <= 4; ++ m)
    {
        for (int j = 0; j <= 4; ++ j)
        {
            Qk[j] = Prim[j][m];
        }
        Primhalf[m] = BasicInterpolation3rd(Qk, m);
    }

    if (Primhalf[0] < SMALL || Primhalf[4] < SMALL)
    {
        Primhalf[0] = Prim[2][0];
        Primhalf[1] = Prim[2][1];
        Primhalf[2] = Prim[2][2];
        Primhalf[3] = Prim[2][3];
        Primhalf[4] = Prim[2][4];
        total_negative += 1;
    }
}

RDouble NSSolverStructFD::BasicWeightedInterpolation5th(RDouble Coef[30], RDouble QmPlate[5])
{
    RDouble g0 = Coef[12] * QmPlate[0] + Coef[13] * QmPlate[1] + Coef[14] * QmPlate[2];
    RDouble g1 = Coef[15] * QmPlate[1] + Coef[16] * QmPlate[2] + Coef[17] * QmPlate[3];
    RDouble g2 = Coef[18] * QmPlate[2] + Coef[19] * QmPlate[3] + Coef[20] * QmPlate[4];

    RDouble s0 = Coef[21] * QmPlate[0] + Coef[22] * QmPlate[1] + Coef[23] * QmPlate[2];
    RDouble s1 = Coef[24] * QmPlate[1] + Coef[25] * QmPlate[2] + Coef[26] * QmPlate[3];
    RDouble s2 = Coef[27] * QmPlate[2] + Coef[28] * QmPlate[3] + Coef[29] * QmPlate[4];

    RDouble IS0 = s0 * s0 + g0 * g0;
    RDouble IS1 = s1 * s1 + g1 * g1;
    RDouble IS2 = s2 * s2 + g2 * g2;
    
    RDouble alfa0 = Coef[0] / (IS0 * IS0 + EPSILON) ;
    RDouble alfa1 = Coef[1] / (IS1 * IS1 + EPSILON) ;
    RDouble alfa2 = Coef[2] / (IS2 * IS2 + EPSILON) ;
    
    RDouble Qhalf0 = Coef[3] * QmPlate[0] + Coef[ 4] * QmPlate[1] + Coef[ 5] * QmPlate[2] ;
    RDouble Qhalf1 = Coef[6] * QmPlate[1] + Coef[ 7] * QmPlate[2] + Coef[ 8] * QmPlate[3] ;
    RDouble Qhalf2 = Coef[9] * QmPlate[2] + Coef[10] * QmPlate[3] + Coef[11] * QmPlate[4] ;

    RDouble Qhalf = (alfa0 * Qhalf0 + alfa1 * Qhalf1 + alfa2 * Qhalf2) / (alfa0 + alfa1 + alfa2);

    return Qhalf;
}

RDouble NSSolverStructFD::BasicInterpolation3rd(RDouble Qp[5], int m)
{
    RDouble dql = Qp[2] - Qp[1];
    RDouble dqr = Qp[3] - Qp[2];

    if (m == 0 || m == 4)
    {
        dql = Qp[2] * (dql / (Qp[2] - 0.5 * dql));
        dqr = Qp[2] * (dqr / (Qp[2] + 0.5 * dqr));
    }    
    
    RDouble ql = Qp[2] + half * Limiter_ThirdSmooth(dqr, dql);

    return ql;
}

RDouble NSSolverStructFD::Limiter_ThirdSmooth(const RDouble &x, const RDouble &y)
{
   const RDouble eps2 = 1.0e-6;
    RDouble Smooth1 = x * (y * y + 2.0 * eps2) +  y * (2.0 * x * x + eps2);
    RDouble Smooth2 = 2.0 * x * x - x * y + 2.0 * y * y + 3.0 * eps2;
    RDouble Smooth = Smooth1 / Smooth2;
    return  Smooth;
}

void NSSolverStructFD::GetFaceVariable_HDCS_firstlayerExplicit(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf)
{
    StructGrid * grid = StructGridCast(gridIn);
    RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q_FD"));
    RDouble4D & ql = qlProxy->GetField_STR();
    RDouble4D & qr = qrProxy->GetField_STR();
    
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble PrimPlateL[5];
    RDouble PrimPlateR[5];

    if (nsurf == 1)
    {
        int ist = 1;
        int jst = 1;
        int kst = 1;
        int ied = ni;
        int jed = nj - 1;
        int ked = nk - 1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }
        
        int ibc[2];
        ibc[0] = ist;
        ibc[1] = ied;
        
        for (int m = 0; m <= 4; ++m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int bcindex = 0; bcindex <= 1; ++ bcindex)
                    {
                        int i = ibc[bcindex];
                        
                        for (int ip = 0; ip <= 4; ++ ip)
                        {
                            int iL = i + (ip - 3);
                            int iR = i + (ip - 2);
                            
                            PrimPlateL[ip]     = q(iL, j, k, m);
                            PrimPlateR[4 - ip] = q(iR, j, k, m);
                        }
                        ql(i - 1, j, k, m) = LinearExplicitInterpolation5th(PrimPlateL);
                        qr(i,     j, k, m) = LinearExplicitInterpolation5th(PrimPlateR);
                    }
                }
            }
        }
        
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int bcindex = 0; bcindex <= 1; ++ bcindex)
                {
                    int i = ibc[bcindex];
                    
                    if (ql(i - 1, j, k, 0) < SMALL || ql(i - 1, j, k, 4) < SMALL)
                    {
                        total_negative += 1;
                        ql(i - 1, j, k, 0) = q(i - 1, j, k, 0);
                        ql(i - 1, j, k, 1) = q(i - 1, j, k, 1);
                        ql(i - 1, j, k, 2) = q(i - 1, j, k, 2);
                        ql(i - 1, j, k, 3) = q(i - 1, j, k, 3);
                        ql(i - 1, j, k, 4) = q(i - 1, j, k, 4);
                    }
                    
                    if (qr(i, j, k, 0) < SMALL || qr(i, j, k, 4) < SMALL)
                    {
                        total_negative += 1;
                        qr(i, j, k, 0) = q(i, j, k, 0);
                        qr(i, j, k, 1) = q(i, j, k, 1);
                        qr(i, j, k, 2) = q(i, j, k, 2);
                        qr(i, j, k, 3) = q(i, j, k, 3);
                        qr(i, j, k, 4) = q(i, j, k, 4);
                    }
                }
            }
        }
    }
    else if (nsurf == 2)
    {
        int ist = 1;
        int jst = 1;
        int kst = 1;
        int ied = ni - 1;
        int jed = nj;
        int ked = nk - 1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }
        
        int jbc[2];
        jbc[0] = jst;
        jbc[1] = jed;
        
        for (int m = 0; m <= 4; ++m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int bcindex = 0; bcindex <= 1; ++ bcindex)
                    {
                        int j = jbc[bcindex];
                        
                        for (int ip = 0; ip <= 4; ++ ip)
                        {
                            int jL = j + (ip - 3);
                            int jR = j + (ip - 2);
                            
                            PrimPlateL[ip]     = q(i, jL, k, m);
                            PrimPlateR[4 - ip] = q(i, jR, k, m);
                        }
                        ql(i, j - 1, k, m) = LinearExplicitInterpolation5th(PrimPlateL);
                        qr(i, j,     k, m) = LinearExplicitInterpolation5th(PrimPlateR);
                    }
                }
            }
        }
        
        for (int k = kst; k <= ked; ++ k)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int bcindex = 0; bcindex <= 1; ++ bcindex)
                {
                    int j = jbc[bcindex];
                    
                    if (ql(i, j - 1, k, 0) < SMALL || ql(i, j - 1, k, 4) < SMALL)
                    {
                        total_negative += 1;
                        ql(i, j - 1, k, 0) = q(i, j - 1, k, 0);
                        ql(i, j - 1, k, 1) = q(i, j - 1, k, 1);
                        ql(i, j - 1, k, 2) = q(i, j - 1, k, 2);
                        ql(i, j - 1, k, 3) = q(i, j - 1, k, 3);
                        ql(i, j - 1, k, 4) = q(i, j - 1, k, 4);
                    }
                    
                    if (qr(i, j, k, 0) < SMALL || qr(i, j, k, 4) < SMALL)
                    {
                        total_negative += 1;
                        qr(i, j, k, 0) = q(i, j, k, 0);
                        qr(i, j, k, 1) = q(i, j, k, 1);
                        qr(i, j, k, 2) = q(i, j, k, 2);
                        qr(i, j, k, 3) = q(i, j, k, 3);
                        qr(i, j, k, 4) = q(i, j, k, 4);
                    }
                }
            }
        }
    }
    else
    {
        int ist = 1;
        int jst = 1;
        int kst = 1;
        int ied = ni - 1;
        int jed = nj - 1;
        int ked = nk;
        
        int kbc[2];
        kbc[0] = kst;
        kbc[1] = ked;
        
        for (int m = 0; m <= 4; ++m)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int bcindex = 0; bcindex <= 1; ++ bcindex)
                    {
                        int k = kbc[bcindex];
                        
                        for (int ip = 0; ip <= 4; ++ ip)
                        {
                            int kL = k + (ip - 3);
                            int kR = k + (ip - 2);
                            
                            PrimPlateL[ip]     = q(i, j, kL, m);
                            PrimPlateR[4 - ip] = q(i, j, kR, m);
                        }
                        ql(i, j, k - 1, m) = LinearExplicitInterpolation5th(PrimPlateL);
                        qr(i, j, k,     m) = LinearExplicitInterpolation5th(PrimPlateR);
                    }
                }
            }
        }
        
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int bcindex = 0; bcindex <= 1; ++ bcindex)
                {
                    int k = kbc[bcindex];
                    
                    if (ql(i, j, k - 1, 0) < SMALL || ql(i, j, k - 1, 4) < SMALL)
                    {
                        total_negative += 1;
                        ql(i, j, k - 1, 0) = q(i, j, k - 1, 0);
                        ql(i, j, k - 1, 1) = q(i, j, k - 1, 1);
                        ql(i, j, k - 1, 2) = q(i, j, k - 1, 2);
                        ql(i, j, k - 1, 3) = q(i, j, k - 1, 3);
                        ql(i, j, k - 1, 4) = q(i, j, k - 1, 4);
                    }
                    
                    if (qr(i, j, k, 0) < SMALL || qr(i, j, k, 4) < SMALL)
                    {
                        total_negative += 1;
                        qr(i, j, k, 0) = q(i, j, k, 0);
                        qr(i, j, k, 1) = q(i, j, k, 1);
                        qr(i, j, k, 2) = q(i, j, k, 2);
                        qr(i, j, k, 3) = q(i, j, k, 3);
                        qr(i, j, k, 4) = q(i, j, k, 4);
                    }
                }
            }
        }
    }

    
}

void NSSolverStructFD::GetFaceVariable_HDCS_otherlayerImplicit(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf)
{
    if (nsurf == 1)
    {
        GetFaceVariable_HDCS_otherlayerImplicit_nsurf1(gridIn, qlProxy, qrProxy);
    }
    else if (nsurf == 2)
    {
        GetFaceVariable_HDCS_otherlayerImplicit_nsurf2(gridIn, qlProxy, qrProxy);
    }
    else
    {
        GetFaceVariable_HDCS_otherlayerImplicit_nsurf3(gridIn, qlProxy, qrProxy);
    }
}

void NSSolverStructFD::GetFaceVariable_HDCS_otherlayerImplicit_nsurf1(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy)
{
    StructGrid * grid = StructGridCast(gridIn);
    RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q_FD"));
    RDouble4D & ql = qlProxy->GetField_STR();
    RDouble4D & qr = qrProxy->GetField_STR();

    const RDouble BL = 0.5;
    const RDouble BR = 3.0 / 14.0;

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni;
    int jed = nj - 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    RDouble *EquivalentPrimhalfL = new RDouble[ni+1];
    RDouble *EquivalentPrimhalfR = new RDouble[ni+1];

    RDouble PrimPlateL[5];
    RDouble PrimPlateR[5];

    int ibc[2];
    ibc[0] = ist;
    ibc[1] = ied;

    for (int m = 0; m <= 4; ++m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int bcindex = 0; bcindex <= 1; ++ bcindex)
                {
                    int i = ibc[bcindex];
                    
                    EquivalentPrimhalfL[i] = ql(i - 1, j, k, m);
                    EquivalentPrimhalfR[i] = qr(i,     j, k, m);
                }
                for (int i = ist+1; i <= ied -1; ++ i)
                {
                    for (int ip = 0; ip <= 4; ++ ip)
                    {
                        int iL = i + (ip - 3);
                        int iR = i + (ip - 2);

                        PrimPlateL[ip]     = q(iL, j, k, m);
                        PrimPlateR[4 - ip] = q(iR, j, k, m);
                    }
                    EquivalentPrimhalfL[i] = LinearImplicitInterpolation7th(PrimPlateL);
                    EquivalentPrimhalfR[i] = LinearImplicitInterpolation7th(PrimPlateR);
                }

                SolveTridiagonalEquations(ni, BL, BR, EquivalentPrimhalfL);
                SolveTridiagonalEquations(ni, BR, BL, EquivalentPrimhalfR);

                for (int i = ist + 1; i <= ied - 1; ++ i)
                {
                    ql(i - 1, j, k, m) = EquivalentPrimhalfL[i];
                    qr(i,     j, k, m) = EquivalentPrimhalfR[i];
                }
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist + 1; i <= ied - 1; ++ i)
            {
                if (ql(i - 1, j, k, 0) < SMALL || ql(i - 1, j, k, 4) < SMALL)
                {
                    total_negative += 1;
                    ql(i - 1, j, k, 0) = q(i - 1, j, k, 0);
                    ql(i - 1, j, k, 1) = q(i - 1, j, k, 1);
                    ql(i - 1, j, k, 2) = q(i - 1, j, k, 2);
                    ql(i - 1, j, k, 3) = q(i - 1, j, k, 3);
                    ql(i - 1, j, k, 4) = q(i - 1, j, k, 4);
                }

                if (qr(i, j, k, 0) < SMALL || qr(i, j, k, 4) < SMALL)
                {
                    total_negative += 1;
                    qr(i, j, k, 0) = q(i, j, k, 0);
                    qr(i, j, k, 1) = q(i, j, k, 1);
                    qr(i, j, k, 2) = q(i, j, k, 2);
                    qr(i, j, k, 3) = q(i, j, k, 3);
                    qr(i, j, k, 4) = q(i, j, k, 4);
                }
            }
        }
    }
    delete [] EquivalentPrimhalfL;
    delete [] EquivalentPrimhalfR;
}

void NSSolverStructFD::GetFaceVariable_HDCS_otherlayerImplicit_nsurf2(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy)
{
    StructGrid * grid = StructGridCast(gridIn);
    RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q_FD"));
    RDouble4D & ql = qlProxy->GetField_STR();
    RDouble4D & qr = qrProxy->GetField_STR();

    const RDouble BL = 0.5;
    const RDouble BR = 3.0 / 14.0;

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1;
    int jed = nj;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    RDouble *EquivalentPrimhalfL = new RDouble[nj+1];
    RDouble *EquivalentPrimhalfR = new RDouble[nj+1];

    RDouble PrimPlateL[5];
    RDouble PrimPlateR[5];

    int jbc[2];
    jbc[0] = jst;
    jbc[1] = jed;

    for (int m = 0; m <= 4; ++m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int bcindex = 0; bcindex <= 1; ++ bcindex)
                {
                    int j = jbc[bcindex];

                    EquivalentPrimhalfL[j] = ql(i, j - 1, k, m);
                    EquivalentPrimhalfR[j] = qr(i, j,     k, m);
                }
                for (int j = jst+1; j <= jed -1; ++ j)
                {
                    for (int ip = 0; ip <= 4; ++ ip)
                    {
                        int jL = j + (ip - 3);
                        int jR = j + (ip - 2);

                        PrimPlateL[ip]     = q(i, jL, k, m);
                        PrimPlateR[4 - ip] = q(i, jR, k, m);
                    }
                    EquivalentPrimhalfL[j] = LinearImplicitInterpolation7th(PrimPlateL);
                    EquivalentPrimhalfR[j] = LinearImplicitInterpolation7th(PrimPlateR);
                }

                SolveTridiagonalEquations(nj, BL, BR, EquivalentPrimhalfL);
                SolveTridiagonalEquations(nj, BR, BL, EquivalentPrimhalfR);

                for (int j = jst + 1; j <= jed - 1; ++ j)
                {
                    ql(i, j - 1, k, m) = EquivalentPrimhalfL[j];
                    qr(i, j,     k, m) = EquivalentPrimhalfR[j];
                }
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst + 1; j <= jed - 1; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                if (ql(i, j - 1, k, 0) < SMALL || ql(i, j - 1, k, 4) < SMALL)
                {
                    total_negative += 1;
                    ql(i, j - 1, k, 0) = q(i, j - 1, k, 0);
                    ql(i, j - 1, k, 1) = q(i, j - 1, k, 1);
                    ql(i, j - 1, k, 2) = q(i, j - 1, k, 2);
                    ql(i, j - 1, k, 3) = q(i, j - 1, k, 3);
                    ql(i, j - 1, k, 4) = q(i, j - 1, k, 4);
                }

                if (qr(i, j, k, 0) < SMALL || qr(i, j, k, 4) < SMALL)
                {
                    total_negative += 1;
                    qr(i, j, k, 0) = q(i, j, k, 0);
                    qr(i, j, k, 1) = q(i, j, k, 1);
                    qr(i, j, k, 2) = q(i, j, k, 2);
                    qr(i, j, k, 3) = q(i, j, k, 3);
                    qr(i, j, k, 4) = q(i, j, k, 4);
                }
            }
        }
    }
    delete [] EquivalentPrimhalfL;
    delete [] EquivalentPrimhalfR;
}

void NSSolverStructFD::GetFaceVariable_HDCS_otherlayerImplicit_nsurf3(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy)
{
    StructGrid * grid = StructGridCast(gridIn);
    RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q_FD"));
    RDouble4D & ql = qlProxy->GetField_STR();
    RDouble4D & qr = qrProxy->GetField_STR();

    const RDouble BL = 0.5;
    const RDouble BR = 3.0 / 14.0;

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1;
    int jed = nj - 1;
    int ked = nk;

    RDouble *EquivalentPrimhalfL = new RDouble[nk+1];
    RDouble *EquivalentPrimhalfR = new RDouble[nk+1];

    RDouble PrimPlateL[5];
    RDouble PrimPlateR[5];

    int kbc[2];
    kbc[0] = kst;
    kbc[1] = ked;

    for (int m = 0; m <= 4; ++m)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int bcindex = 0; bcindex <= 1; ++ bcindex)
                {
                    int k = kbc[bcindex];

                    EquivalentPrimhalfL[k] = ql(i, j, k - 1, m);
                    EquivalentPrimhalfR[k] = qr(i, j, k,     m);
                }
                for (int k = kst+1; k <= ked -1; ++ k)
                {
                    for (int ip = 0; ip <= 4; ++ ip)
                    {
                        int kL = k + (ip - 3);
                        int kR = k + (ip - 2);

                        PrimPlateL[ip]     = q(i, j, kL, m);
                        PrimPlateR[4 - ip] = q(i, j, kR, m);
                    }
                    EquivalentPrimhalfL[k] = LinearImplicitInterpolation7th(PrimPlateL);
                    EquivalentPrimhalfR[k] = LinearImplicitInterpolation7th(PrimPlateR);
                }

                SolveTridiagonalEquations(nk, BL, BR, EquivalentPrimhalfL);
                SolveTridiagonalEquations(nk, BR, BL, EquivalentPrimhalfR);

                for (int k = kst + 1; k <= ked - 1; ++ k)
                {
                    ql(i, j, k - 1, m) = EquivalentPrimhalfL[k];
                    qr(i, j, k,     m) = EquivalentPrimhalfR[k];
                }
            }
        }
    }

    for (int k = kst + 1; k <= ked - 1; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                if (ql(i, j, k - 1, 0) < SMALL || ql(i, j, k - 1, 4) < SMALL)
                {
                    total_negative += 1;
                    ql(i, j, k - 1, 0) = q(i, j, k - 1, 0);
                    ql(i, j, k - 1, 1) = q(i, j, k - 1, 1);
                    ql(i, j, k - 1, 2) = q(i, j, k - 1, 2);
                    ql(i, j, k - 1, 3) = q(i, j, k - 1, 3);
                    ql(i, j, k - 1, 4) = q(i, j, k - 1, 4);
                }

                if (qr(i, j, k, 0) < SMALL || qr(i, j, k, 4) < SMALL)
                {
                    total_negative += 1;
                    qr(i, j, k, 0) = q(i, j, k, 0);
                    qr(i, j, k, 1) = q(i, j, k, 1);
                    qr(i, j, k, 2) = q(i, j, k, 2);
                    qr(i, j, k, 3) = q(i, j, k, 3);
                    qr(i, j, k, 4) = q(i, j, k, 4);
                }
            }
        }
    }
    delete [] EquivalentPrimhalfL;
    delete [] EquivalentPrimhalfR;
}

RDouble NSSolverStructFD::LinearExplicitInterpolation5th(RDouble prim[5])
{
    RDouble primhalf = (3.0 * prim[0] - 20.0 * prim[1] + 90.0 * prim[2] + 60.0 * prim[3] - 5.0 * prim[4]) / 128.0;

    return primhalf;
}

RDouble NSSolverStructFD::LinearImplicitInterpolation7th(RDouble prim[5])
{
    const RDouble a = - 1.0 / 224.0;
    const RDouble b =   1.0 /   8.0;
    const RDouble c =  15.0 /  16.0;
    const RDouble d =   5.0 /   8.0;
    const RDouble e =   1.0 /  32.0;

    RDouble primhalf = a * prim[0] + b * prim[1] + c * prim[2] + d * prim[3] + e * prim[4];

    return primhalf;
}

void NSSolverStructFD::SolveTridiagonalEquations(int Nod, RDouble BL, RDouble BR, RDouble *F)
{
    RDouble *C = new RDouble[Nod+1];

    C[2] = BR;
    F[2] = F[2] - BL * F[1];

    for (int k = 3; k <= Nod - 1; ++ k)
    {
        RDouble Inv = 1.0 / (1.0 - BL * C[k-1]);
        C[k] = Inv * BR;
        F[k] = Inv * (F[k] - BL * F[k-1]);
    }

    for (int k = Nod - 1; k >= 2; -- k)
    {
        F[k] = F[k] - C[k] * F[k+1];
    }
    delete [] C;
}

void NSSolverStructFD::GetFaceInviscidFlux(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, FieldProxy *fluxProxy, int nsurf)
{
    if (fluxname.substr(0, 3) == "roe")
    {
        FluxRoe(gridIn, qlProxy, qrProxy, fluxProxy, nsurf);
    }
    else
    {
        FluxSteger(gridIn, qlProxy, qrProxy, fluxProxy, nsurf);
    }
}

void NSSolverStructFD::FluxSteger(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, FieldProxy *fluxProxy, int nsurf)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &ql = qlProxy->GetField_STR();
    RDouble4D &qr = qrProxy->GetField_STR();
    RDouble4D &Flux_Steger = fluxProxy->GetField_STR();

    RDouble4D &xfn  = *(grid->GetFaceNormalX());
    RDouble4D &yfn  = *(grid->GetFaceNormalY());
    RDouble4D &zfn  = *(grid->GetFaceNormalZ());
    RDouble4D &area = *(grid->GetFaceArea());

    const RDouble gama = 1.4;
    const RDouble epsilon_steger = 0.12 * 0.12;

    using namespace IDX;
    
    int il1 = 0, jl1 = 0, kl1 = 0;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    int ist = 0, jst = 0, kst = 0;
    int ied = 0, jed = 0, ked = 0;

    if (solvername.substr(0, 4) == "WCNS")
    {
        ist = 1 - il1;
        jst = 1 - jl1;
        kst = 1 - kl1;
        ied = ni - 1 + 2*il1;
        jed = nj - 1 + 2*jl1;
        ked = nk - 1 + 2*kl1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }
    }
    else if (solvername.substr(0, 4) == "HDCS")
    {
        ist = 1;
        jst = 1;
        kst = 1;
        ied = ni - 1 + il1;
        jed = nj - 1 + jl1;
        ked = nk - 1 + kl1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }
    }
    

    RDouble fl[5] = {0};
    RDouble fr[5] = {0};
    RDouble priml[5] = {0};
    RDouble primr[5] = {0};

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m < 5; ++ m)
                {
                    priml[m] = ql(i - il1, j - jl1, k - kl1, m);
                    primr[m] = qr(i,       j,       k,       m);
                }

                RDouble ds = area(i, j, k, nsurf);
                RDouble nx = xfn (i, j, k, nsurf);
                RDouble ny = yfn (i, j, k, nsurf);
                RDouble nz = zfn (i, j, k, nsurf);

                RDouble rl = priml[IR];
                RDouble ul = priml[IU];
                RDouble vl = priml[IV];
                RDouble wl = priml[IW];
                RDouble pl = priml[IP];
                RDouble hint_l = (gama / (gama - one)) * (pl / rl);
                RDouble hl = hint_l + half * (ul * ul + vl * vl + wl * wl);
                RDouble el = hl - pl / (rl + SMALL);
                RDouble vnl = nx * ul + ny * vl + nz * wl;
                RDouble cl  = sqrt(ABS(gama * pl / rl));
                RDouble eigl1 = vnl      ;
                RDouble eigl2 = vnl + cl ;
                RDouble eigl3 = vnl - cl ;
                
                RDouble rr = primr[IR];
                RDouble ur = primr[IU];
                RDouble vr = primr[IV];
                RDouble wr = primr[IW];
                RDouble pr = primr[IP];
                RDouble hint_r = (gama / (gama - one)) * (pr / rr);
                RDouble hr = hint_r + half * (ur * ur + vr * vr + wr * wr);
                RDouble er = hr - pr / (rr + SMALL);
                RDouble vnr = nx * ur + ny * vr + nz * wr;
                RDouble cr  = sqrt(ABS(gama * pr / rr)); 
                RDouble eigr1 = vnr      ;
                RDouble eigr2 = vnr + cr ;
                RDouble eigr3 = vnr - cr ;
                
                eigl1 = half * (eigl1 +  sqrt (eigl1 * eigl1 + epsilon_steger));
                eigl2 = half * (eigl2 +  sqrt (eigl2 * eigl2 + epsilon_steger));
                eigl3 = half * (eigl3 +  sqrt (eigl3 * eigl3 + epsilon_steger)); 
                
                eigr1 = half * (eigr1 -  sqrt (eigr1 * eigr1 + epsilon_steger));
                eigr2 = half * (eigr2 -  sqrt (eigr2 * eigr2 + epsilon_steger));
                eigr3 = half * (eigr3 -  sqrt (eigr3 * eigr3 + epsilon_steger));
                
                RDouble xil1 = half * (eigl1 + eigl1 - eigl2 - eigl3) / gama;
                RDouble xil2 = half * (eigl2 - eigl3) * cl / gama;
                RDouble eig_xil = eigl1 - xil1;

                RDouble xir1 = half * (eigr1 + eigr1 - eigr2 - eigr3) / gama;
                RDouble xir2 = half * (eigr2 - eigr3) * cr / gama;
                RDouble eig_xir = eigr1 - xir1;
                
                fl[IR ] = rl * (eig_xil);
                fl[IRU] = rl * (eig_xil * ul + xil2 * nx);
                fl[IRV] = rl * (eig_xil * vl + xil2 * ny);
                fl[IRW] = rl * (eig_xil * wl + xil2 * nz);
                fl[IRE] = rl * (eigl1   * el - xil1 * hl + xil2 * vnl);  
                
                fr[IR ] = rr * (eig_xir);
                fr[IRU] = rr * (eig_xir * ur + xir2 * nx);
                fr[IRV] = rr * (eig_xir * vr + xir2 * ny);
                fr[IRW] = rr * (eig_xir * wr + xir2 * nz);
                fr[IRE] = rr * (eigr1   * er - xir1 * hr + xir2 * vnr);
                
                for (int m = 0; m < 5; ++ m)
                {
                    Flux_Steger(i, j, k, m) = (fl[m] + fr[m]) * ds;
                }
            }
        }
    }
}

void NSSolverStructFD::FluxRoe(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, FieldProxy *fluxProxy, int nsurf)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();
    int nl = 5;

    RDouble4D &ql = qlProxy->GetField_STR();
    RDouble4D &qr = qrProxy->GetField_STR();
    RDouble4D &Flux_Roe = fluxProxy->GetField_STR();

    RDouble4D &xfn  = *(grid->GetFaceNormalX());
    RDouble4D &yfn  = *(grid->GetFaceNormalY());
    RDouble4D &zfn  = *(grid->GetFaceNormalZ());
    RDouble4D &area = *(grid->GetFaceArea());
    RDouble KX = 0, KY = 0, KZ = 0, KA = 0;

    const RDouble gama = 1.4;
    const RDouble K1 = 0.25;
    const RDouble K2 =  5.0;
    const RDouble K3 = 15.0;
    RDouble3D &rtem = * reinterpret_cast< RDouble3D * > (grid->GetDataPtr("rtem"));

    using namespace IDX;
    
    int il1 = 0, jl1 = 0, kl1 = 0;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    int ist = 0, jst = 0, kst = 0;
    int ied = 0, jed = 0, ked = 0;

    if (solvername.substr(0, 4) == "WCNS")
    {
        ist = 1 - il1;
        jst = 1 - jl1;
        kst = 1 - kl1;
        ied = ni - 1 + 2*il1;
        jed = nj - 1 + 2*jl1;
        ked = nk - 1 + 2*kl1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }
    }
    else if (solvername.substr(0, 4) == "HDCS")
    {
        ist = 1;
        jst = 1;
        kst = 1;
        ied = ni - 1 + il1;
        jed = nj - 1 + jl1;
        ked = nk - 1 + kl1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }
    }

    RDouble *flux  = new RDouble[nl]();
    RDouble *priml = new RDouble[nl]();
    RDouble *primr = new RDouble[nl]();
    
    RDouble RRL = 0,UUL = 0,VVL = 0,WWL = 0,PPL = 0,HHL = 0,KPL = 0;
    RDouble RRR = 0,UUR = 0,VVR = 0,WWR = 0,PPR = 0,HHR = 0,KPR = 0;
    RDouble RAT = 0,RCF = 0,RF = 0,UF = 0,VF = 0,WF = 0,HF = 0;
    RDouble VN = 0,RVN = 0,DVN = 0;
    RDouble U2 = 0,A2 = 0,AF = 0;
    RDouble DR = 0,DU = 0,DV = 0,DW = 0,DP = 0;
    RDouble LA1 = 0,LA3 = 0,LA5 = 0;
    RDouble DLP = 0,DLR = 0;
    RDouble MF = 0,CP = 0,LAM = 0,ESP1 = 0,ESP3 = 0,ESP5 = 0;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m < nl; ++ m)
                {
                    priml[m] = ql(i - il1, j - jl1, k - kl1, m);
                    primr[m] = qr(i,       j,       k,       m);
                }

                KA = area(i, j, k, nsurf);
                KX = xfn (i, j, k, nsurf);
                KY = yfn (i, j, k, nsurf);
                KZ = zfn (i, j, k, nsurf);

                KPL = rtem(i - il1, j - jl1, k - kl1);
                KPR = rtem(i, j, k);

                RRL = priml[IR];
                UUL = priml[IU];
                VVL = priml[IV];
                WWL = priml[IW];
                PPL = priml[IP];
                
                RRR = primr[IR];
                UUR = primr[IU];
                VVR = primr[IV];
                WWR = primr[IW];
                PPR = primr[IP];

                RAT = sqrt(RRR / RRL);
                RCF = one / (one + RAT);
                
                RF = sqrt(RRL * RRR);
                UF = (UUL + UUR * RAT) * RCF;
                VF = (VVL + VVR * RAT) * RCF;
                WF = (WWL + WWR * RAT) * RCF;
                
                HHL = (gama / (gama - one)) * PPL / (RRL + SMALL) + half * (UUL * UUL + VVL * VVL + WWL * WWL);
                HHR = (gama / (gama - one)) * PPR / (RRR + SMALL) + half * (UUR * UUR + VVR * VVR + WWR * WWR);
                HF = (HHL + HHR * RAT) * RCF;
                
                //! Left full flux.
                VN  = KX * UUL + KY * VVL + KZ * WWL;
                RVN = VN * RRL;

                flux[IR ] = RVN                 ;
                flux[IRU] = RVN * UUL + KX * PPL;
                flux[IRV] = RVN * VVL + KY * PPL;
                flux[IRW] = RVN * WWL + KZ * PPL;
                flux[IRE] = RVN * HHL;
                
                //! Right full flux.
                VN  = KX * UUR + KY * VVR + KZ * WWR;
                RVN = VN * RRR;
                
                flux[IR ] += RVN                 ;
                flux[IRU] += RVN * UUR + KX * PPR;
                flux[IRV] += RVN * VVR + KY * PPR;
                flux[IRW] += RVN * WWR + KZ * PPR;
                flux[IRE] += RVN * HHR;

                //! Delta flux.
                U2 = half * (UF * UF + VF * VF + WF * WF);
                A2 = (gama - one) * (HF - U2);
                AF = sqrt(A2);
                
                DR = RRR - RRL;
                DU = UUR - UUL;
                DV = VVR - VVL;
                DW = WWR - WWL;
                DP = PPR - PPL;
                
                VN  = KX * UF + KY * VF + KZ * WF;
                DVN = KX * DU + KY * DV + KZ * DW;
                
                LA1 = abs(VN - AF);
                LA3 = abs(VN);
                LA5 = abs(VN + AF);
                
                //! Delta flux: entropy fix.
                MF   = sqrt (UF * UF + VF * VF + WF * WF) / (AF + SMALL);
                CP   = half * (KPL + KPR) / MAX(one,MF);
                LAM  = MAX(LA1, LA5);
                ESP1 = (K1 + K2 * CP) * LAM;
                ESP3 =        K3 * CP   * LAM;
                ESP5 =  ESP1;
                
                if (LA1 < ESP1)
                {
                    LA1 = half * (LA1 * LA1 + ESP1 * ESP1) / (ESP1 + SMALL);
                }
                if (LA3 < ESP3)
                {
                    LA3 = half * (LA3 * LA3 + ESP3 * ESP3) / (ESP3 + SMALL);
                }
                if (LA5 < ESP5)
                {
                    LA5 = half * (LA5 * LA5 + ESP5 * ESP5) / (ESP5 + SMALL);
                }
                
                //! Delta flux: LA1 & LA5.
                LA1 *= half * (DP - RF * AF * DVN) / (A2  + SMALL);
                LA5 *= half * (DP + RF * AF * DVN) / (A2  + SMALL);
                
                flux[IR ] -= LA1                    + LA5                   ;
                flux[IRU] -= LA1 * (UF - KX * AF) + LA5 * (UF + KX * AF);
                flux[IRV] -= LA1 * (VF - KY * AF) + LA5 * (VF + KY * AF);
                flux[IRW] -= LA1 * (WF - KZ * AF) + LA5 * (WF + KZ * AF);
                flux[IRE] -= LA1 * (HF - VN * AF) + LA5 * (HF + VN * AF);
                
                //! Delta flux: LA3.
                DLP = LA3 * (DR - DP / (A2  + SMALL));
                DLR = LA3 * RF;
                
                flux[IR ] -= DLP;
                flux[IRU] -= DLP * UF + DLR * (- KX * DVN + DU);
                flux[IRV] -= DLP * VF + DLR * (- KY * DVN + DV);
                flux[IRW] -= DLP * WF + DLR * (- KZ * DVN + DW);
                flux[IRE] -= DLP * U2 + DLR * (- VN * DVN + UF * DU + VF * DV + WF * DW);
                                
                for (int m = 0; m < nl; ++ m)
                {
                    Flux_Roe(i, j, k, m) = flux[m] * half * KA;
                }
            }
        }
    }
    delete [] flux;    flux = NULL;
    delete [] priml;    priml = NULL;
    delete [] primr;    primr = NULL;
}

void NSSolverStructFD::FluxDifference(Grid *gridIn, FieldProxy *fluxProxy, int nsurf)
{
    if (solvername.substr(0, 4) == "WCNS")
    {
        FluxDifference_WCNS(gridIn, fluxProxy, nsurf);
    }
    else if (solvername.substr(0, 4) == "HDCS")
    {
        FluxDifference_HDCS(gridIn, fluxProxy, nsurf);
    }
    else
    {
        FluxDifference_WCNS(gridIn, fluxProxy, nsurf);
    }
}

void NSSolverStructFD::FluxDifference_WCNS(Grid *gridIn, FieldProxy *fluxProxy, int nsurf)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &flux = fluxProxy->GetField_STR();
    RDouble4D &res  = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("res"));
    
    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    int mst = 0;
    int med = 4;

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1;
    int jed = nj - 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        ked = 1;
    } 

    Int1D *DiffOrder;
    if (nsurf == 1)
    {
        DiffOrder = grid->DiscretePrecisionKXI;
    }
    else if (nsurf == 2)
    {
        DiffOrder = grid->DiscretePrecisionETA;
    }
    else
    {
        DiffOrder = grid->DiscretePrecisionCTA;
    }

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble f0 = flux(i - il1,     j - jl1,     k - kl1,     m);
                    RDouble f1 = flux(i,           j,           k,           m);
                    RDouble f2 = flux(i + il1,     j + jl1,     k + kl1,     m);
                    RDouble f3 = flux(i + il1 * 2, j + jl1 * 2, k + kl1 * 2, m);

                    int localcellID = i * il1 + j * jl1 + k * kl1;

                    if ((*DiffOrder)(localcellID) == 2)
                    {
                        res(i, j, k, m) -= (f2 - f1);
                    }
                    else
                    {
                        res(i, j, k, m) -= (f2 - f1) * 9.0 / 8.0 - (f3 - f0) / 24.0;
                    }                    
                }
            }
        }
    }    
}

void NSSolverStructFD::FluxDifference_HDCS(Grid *gridIn, FieldProxy *fluxProxy, int nsurf)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();
    
    RDouble4D &faceflux = fluxProxy->GetField_STR();
    RDouble4D &Cellflux = CellfluxProxy->GetField_STR();
    RDouble4D &res  = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("res"));

    const RDouble alfa_sixth =  64.0 / 45.0;
    const RDouble    a_sixth = - 2.0 / 9.0;
    const RDouble    b_sixth =   1.0 / 180.0;

    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);
    
    int mst = 0;
    int med = 4;

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1;
    int jed = nj - 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        ked = 1;
    }

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble leftface  = faceflux(i,           j,           k,           m);
                    RDouble rightface = faceflux(i + il1,     j + jl1,     k + kl1,     m);

                    RDouble celll1 = Cellflux(i - il1, j - jl1, k - kl1, m);
                    RDouble cellr1 = Cellflux(i + il1, j + jl1, k + kl1, m);

                    RDouble celll2 = Cellflux(i - 2*il1, j - 2*jl1, k - 2*kl1, m);
                    RDouble cellr2 = Cellflux(i + 2*il1, j + 2*jl1, k + 2*kl1, m);

                    res(i, j, k, m) -= alfa_sixth * (rightface - leftface) + a_sixth * (cellr1 - celll1) + b_sixth * (cellr2 - celll2);
                }
            }
        }
    }
}


void NSSolverStructFD::GetCellInviscidFlux(Grid *gridIn, int nsurf)
{
    if (solvername.substr(0, 4) != "HDCS")
    {
        return;
    }
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();
    
    RDouble4D &q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q_FD"));
    RDouble4D &Cellflux = CellfluxProxy->GetField_STR();
    
    RDouble4D &xfv  = *(grid->GetFaceVectorX());
    RDouble4D &yfv  = *(grid->GetFaceVectorY());
    RDouble4D &zfv  = *(grid->GetFaceVectorZ());
    
    const RDouble gama = 1.4;
    
    using namespace IDX;
    
    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    int ist = 1 - 2*il1;
    int jst = 1 - 2*jl1;
    int kst = 1 - 2*kl1;
    int ied = ni - 1 + 2*il1;
    int jed = nj - 1 + 2*jl1;
    int ked = nk - 1 + 2*kl1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }
    
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble KX = half * (xfv(i, j, k, nsurf) + xfv(i + il1, j + jl1, k + kl1, nsurf));
                RDouble KY = half * (yfv(i, j, k, nsurf) + yfv(i + il1, j + jl1, k + kl1, nsurf));
                RDouble KZ = half * (zfv(i, j, k, nsurf) + zfv(i + il1, j + jl1, k + kl1, nsurf));
                
                RDouble R = q(i, j, k, IR);
                RDouble U = q(i, j, k, IU);
                RDouble V = q(i, j, k, IV);
                RDouble W = q(i, j, k, IW);
                RDouble P = q(i, j, k, IP);
                RDouble H = (gama / (gama - one)) * P / R + half * (U * U + V * V + W * W);
                
                RDouble RVN = R * (KX * U + KY * V + KZ * W);
                
                Cellflux(i, j, k, IR) = RVN;
                Cellflux(i, j, k, IRU) = RVN * U + KX * P;
                Cellflux(i, j, k, IRV) = RVN * V + KY * P;
                Cellflux(i, j, k, IRW) = RVN * W + KZ * P;
                Cellflux(i, j, k, IRE) = RVN * H;
            }
        }
    }
}

void NSSolverStructFD::GetCellViscousFlux(Grid *gridIn, int iSurface)
{
    if (solvername.substr(0, 4) != "HDCS")
    {
        return;
    }
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();
    
    RDouble4D &dqdx = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &dqdy = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &dqdz = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));
    RDouble4D &UVWT = UVWTproxy->GetField_STR();

    RDouble4D &Cellflux   = CellfluxProxy->GetField_STR();

    RDouble4D &xfv  = *(grid->GetFaceVectorX());
    RDouble4D &yfv  = *(grid->GetFaceVectorY());
    RDouble4D &zfv  = *(grid->GetFaceVectorZ());
    
    RDouble3D &visl = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &vist = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));

    Param_NSSolverStruct *parameters = GetControlParameters();
    RDouble oPrandtlLaminar    = parameters->GetoPrandtlLaminar();
    RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble specificHeatAtConstantPressure = 1.0 / ((1.4 - 1.0) * refMachNumber * refMachNumber);

    int il1, jl1, kl1;
    GetNsurfIndex(iSurface, il1, jl1, kl1);

    int ist = 1 - 2*il1;
    int jst = 1 - 2*jl1;
    int kst = 1 - 2*kl1;
    int ied = ni - 1 + 2*il1;
    int jed = nj - 1 + 2*jl1;
    int ked = nk - 1 + 2*kl1;
    
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    RDouble fvis[5] = {0};
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble nx = half * (xfv(i, j, k, iSurface) + xfv(i + il1, j + jl1, k + kl1, iSurface));
                RDouble ny = half * (yfv(i, j, k, iSurface) + yfv(i + il1, j + jl1, k + kl1, iSurface));
                RDouble nz = half * (zfv(i, j, k, iSurface) + zfv(i + il1, j + jl1, k + kl1, iSurface));

                RDouble um  = UVWT(i, j, k, 0);
                RDouble vm  = UVWT(i, j, k, 1);
                RDouble wm  = UVWT(i, j, k, 2);

                RDouble mul = visl(i, j, k);
                RDouble mut = vist(i, j, k);
                RDouble vis = mul + mut;

                RDouble dudx  = dqdx(i, j, k, 0);
                RDouble dudy  = dqdy(i, j, k, 0);
                RDouble dudz  = dqdz(i, j, k, 0);

                RDouble dvdx  = dqdx(i, j, k, 1);
                RDouble dvdy  = dqdy(i, j, k, 1);
                RDouble dvdz  = dqdz(i, j, k, 1);

                RDouble dwdx  = dqdx(i, j, k, 2);
                RDouble dwdy  = dqdy(i, j, k, 2);
                RDouble dwdz  = dqdz(i, j, k, 2);

                RDouble dtdx  = dqdx(i, j, k, 3);
                RDouble dtdy  = dqdy(i, j, k, 3);
                RDouble dtdz  = dqdz(i, j, k, 3);
                
                RDouble kcp = (mul * oPrandtlLaminar + mut * oPrandtlTurbulence) * specificHeatAtConstantPressure;

                RDouble qx = kcp * dtdx;
                RDouble qy = kcp * dtdy;
                RDouble qz = kcp * dtdz;

                RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

                RDouble txx = vis * (two * dudx - divv2p3);
                RDouble tyy = vis * (two * dvdy - divv2p3);
                RDouble tzz = vis * (two * dwdz - divv2p3);
                RDouble txy = vis * (dudy + dvdx);
                RDouble txz = vis * (dudz + dwdx);
                RDouble tyz = vis * (dvdz + dwdy);

                fvis[0] = 0.0;
                fvis[1] = nx * txx + ny * txy + nz * txz;
                fvis[2] = nx * txy + ny * tyy + nz * tyz;
                fvis[3] = nx * txz + ny * tyz + nz * tzz;
                fvis[4] = um * fvis[1] + vm * fvis[2] + wm * fvis[3] + (nx * qx + ny * qy + nz * qz);

                for (int m = 0; m < 5; ++ m)
                {
                    Cellflux(i, j, k, m) = - fvis[m] / refReNumber;
                }
            }
        }
    }
}

void NSSolverStructFD::GetUVWTproxy(Grid *gridIn)
{
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &q = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("q_FD"));
    RDouble4D &t = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("t_FD"));
    RDouble4D &UVWT = UVWTproxy->GetField_STR();

    int ist = - 3;
    int jst = - 3;
    int kst = - 3;
    int ied = ni + 3;
    int jed = nj + 3;
    int ked = nk + 3;
    
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                UVWT(i, j, k, 0) = q(i, j, k, 1);
                UVWT(i, j, k, 1) = q(i, j, k, 2);
                UVWT(i, j, k, 2) = q(i, j, k, 3);
                UVWT(i, j, k, 3) = t(i, j, k, 0);
            }
        }
    }
}

void NSSolverStructFD::GetUVWTfaceproxy(Grid *gridIn, int nsurf)
{
    if (solvername.substr(0, 4) == "HDCS")
    {
        GetUVWTfaceproxy_HDCS(gridIn, nsurf);
        return;
    }
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &UVWT = UVWTproxy->GetField_STR();

    RDouble4D *UVWTface;
    Int1D *DiffOrder;
    if (nsurf == 1)
    {
        UVWTface = &(UVWTfaceIproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionKXI;
    }
    else if (nsurf == 2)
    {
        UVWTface = &(UVWTfaceJproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionETA;
    }
    else
    {
        UVWTface = &(UVWTfaceKproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionCTA;
    }
    
    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    int ist = 1 - il1;
    int jst = 1 - jl1;
    int kst = 1 - kl1;
    int ied = ni - 1 + 2*il1;
    int jed = nj - 1 + 2*jl1;
    int ked = nk - 1 + 2*kl1;    
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int localcellID = i * il1 + j * jl1 + k * kl1;
                
                if ((*DiffOrder)(localcellID) == 2 || (*DiffOrder)(localcellID - 1) == 2)
                {
                    for (int m = 0; m <= 3; ++m)
                    {                    
                        RDouble b = UVWT(i - il1,   j - jl1,   k - kl1,   m);
                        RDouble c = UVWT(i,         j,         k,         m);
                        (*UVWTface)(i, j, k, m) = half * (b + c);
                    }
                }
                else
                {
                    for (int m = 0; m <= 3; ++m)
                    {                    
                        RDouble a = UVWT(i - il1*2, j - jl1*2, k - kl1*2, m);
                        RDouble b = UVWT(i - il1,   j - jl1,   k - kl1,   m);
                        RDouble c = UVWT(i,         j,         k,         m);
                        RDouble d = UVWT(i + il1,   j + jl1,   k + kl1,   m);
                        (*UVWTface)(i, j, k, m) = (- a + 9.0 * b + 9.0 * c - d) / 16.0;
                    }
                    if ((*UVWTface)(i, j, k, 3) < SMALL)
                    {
                        (*UVWTface)(i, j, k, 3) = half * (UVWT(i - il1, j - jl1, k - kl1,  3) + UVWT(i, j, k, 3));
                    }
                }                
            }
        }
    }   
}

void NSSolverStructFD::GetUVWTfaceproxy_HDCS(Grid *gridIn, int nsurf)
{
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &UVWT = UVWTproxy->GetField_STR();

    RDouble4D *UVWTface;
    if (nsurf == 1)
    {
        UVWTface = &(UVWTfaceIproxy->GetField_STR());
    }
    else if (nsurf == 2)
    {
        UVWTface = &(UVWTfaceJproxy->GetField_STR());
    }
    else
    {
        UVWTface = &(UVWTfaceKproxy->GetField_STR());
    }
    
    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1 + il1;
    int jed = nj - 1 + jl1;
    int ked = nk - 1 + kl1;    
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m <= 3; ++m)
                {
                    RDouble a = UVWT(i - il1*3, j - jl1*3, k - kl1*3, m);
                    RDouble b = UVWT(i - il1*2, j - jl1*2, k - kl1*2, m);
                    RDouble c = UVWT(i - il1,   j - jl1,   k - kl1,   m);
                    RDouble d = UVWT(i,         j,         k,         m);
                    RDouble e = UVWT(i + il1,   j + jl1,   k + kl1,   m);
                    RDouble f = UVWT(i + il1*2, j + jl1*2, k + kl1*2, m);
                    (*UVWTface)(i, j, k, m) = (150.0 * (c + d) - 25.0 * (b + e) + 3.0 * (a + f)) / 256.0;
                }
                if ((*UVWTface)(i, j, k, 3) < SMALL)
                {
                    (*UVWTface)(i, j, k, 3) = half * (UVWT(i - il1, j - jl1, k - kl1,  3) + UVWT(i, j, k, 3));
                }
            }
        }
    }   
}

void NSSolverStructFD::GetdqdkxietactaCellCenter(Grid *gridIn, int nsurf)
{
    if (solvername.substr(0, 4) == "HDCS")
    {
        GetdqdkxietactaCellCenter_HDCS(gridIn, nsurf);
        return;
    }
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D *UVWTface;
    Int1D *DiffOrder;
    RDouble4D *dqdkxietacta;
    if (nsurf == 1)
    {
        UVWTface = &(UVWTfaceIproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionKXI;
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdkxi")));
    }
    else if (nsurf == 2)
    {
        UVWTface = &(UVWTfaceJproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionETA;
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdeta")));
    }
    else
    {
        UVWTface = &(UVWTfaceKproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionCTA;
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdcta")));
    }

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1;
    int jed = nj - 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int localcellID = i * il1 + j * jl1 + k * kl1;
                
                if ((*DiffOrder)(localcellID) == 2)
                {
                    for (int m = 0; m <= 3; ++ m)
                    {
                        RDouble q1 = (*UVWTface)(i      , j      , k      , m);
                        RDouble q2 = (*UVWTface)(i + il1, j + jl1, k + kl1, m);

                        (*dqdkxietacta)(i, j, k, m) = q2 - q1;
                    }
                }
                else
                {
                    for (int m = 0; m <= 3; ++ m)
                    {
                        RDouble q0 = (*UVWTface)(i - il1    , j - jl1    , k - kl1    , m);
                        RDouble q1 = (*UVWTface)(i          , j          , k          , m);
                        RDouble q2 = (*UVWTface)(i + il1    , j + jl1    , k + kl1    , m);
                        RDouble q3 = (*UVWTface)(i + il1 * 2, j + jl1 * 2, k + kl1 * 2, m);

                        (*dqdkxietacta)(i, j, k, m) = (27.0 * (q2 - q1) - (q3 - q0)) / 24.0;
                    }
                }
            }
        }
    }   
}

void NSSolverStructFD::GetdqdkxietactaCellCenter_HDCS(Grid *gridIn, int nsurf)
{
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    const RDouble alfa_sixth =  64.0 / 45.0;
    const RDouble    a_sixth = - 2.0 / 9.0;
    const RDouble    b_sixth =   1.0 / 180.0;

    RDouble4D &UVWT = UVWTproxy->GetField_STR();
    RDouble4D *UVWTface;
    RDouble4D *dqdkxietacta;
    if (nsurf == 1)
    {
        UVWTface = &(UVWTfaceIproxy->GetField_STR());
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdkxi")));
    }
    else if (nsurf == 2)
    {
        UVWTface = &(UVWTfaceJproxy->GetField_STR());
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdeta")));
    }
    else
    {
        UVWTface = &(UVWTfaceKproxy->GetField_STR());
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdcta")));
    }

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1;
    int jed = nj - 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m <= 3; ++ m)
                {
                    RDouble fL = (*UVWTface)(i          , j          , k          , m);
                    RDouble fR = (*UVWTface)(i + il1    , j + jl1    , k + kl1    , m);

                    RDouble CL2 = UVWT(i - il1*2, j - jl1*2, k - kl1*2, m);
                    RDouble CL1 = UVWT(i - il1,   j - jl1,   k - kl1,   m);                    
                    RDouble CR1 = UVWT(i + il1,   j + jl1,   k + kl1,   m);
                    RDouble CR2 = UVWT(i + il1*2, j + jl1*2, k + kl1*2, m);
                    
                    (*dqdkxietacta)(i, j, k, m) = alfa_sixth * (fR - fL) + a_sixth * (CR1 - CL1) + b_sixth * (CR2 - CL2);
                }
            }
        }
    }   
}


void NSSolverStructFD::ViscousFlux(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    Param_NSSolverStruct *parameters = GetControlParameters();

    int viscousType = parameters->GetViscousType();
    if (viscousType <= INVISCID)
    {
        return;
    }

    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");

    FieldProxy * fluxProxy = new FieldProxy();

    int nDim = GetDim();
    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        if (iSurface == 1)
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxI")));
        }
        else if (iSurface == 2)
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxJ")));
        }
        else
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxK")));
        }
        GetGradientAtFace(gridIn, iSurface);
        GetGradientAtFaceCorrectionAtPhysicalBoundary(gridIn, iSurface);

        if (iLES == NOLES_SOLVER)
        {
            GetFaceViscousFlux(gridIn, fluxProxy, iSurface);
        }
        else
        {
            GetFaceViscousFluxLES(gridIn, fluxProxy, iSurface);
        }
        
        FluxDifference(gridIn, fluxProxy, iSurface);
    }
    delete fluxProxy;
}

void NSSolverStructFD::GetGradientAtFace(Grid *gridIn, int nSurface)
{
    if (gradient_method.substr(0,12) == "conservation")
    {
        GetGradientAtFaceconservation(gridIn, nSurface);
    }
    else if (gradient_method.substr(0,10) == "chain_rule")
    {
        GetGradientAtFacechainrule(gridIn, nSurface);
    }
    else
    {
        GetGradientAtFaceconservation(gridIn, nSurface);
    }
}

void NSSolverStructFD::GetGradientAtFaceconservation(Grid *gridIn, int nSurface)
{
    if (solvername.substr(0, 4) == "HDCS")
    {
        GetGradientAtFaceconservation_HDCS(gridIn, nSurface);
        return;
    }
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();    

    RDouble4D &dqdx = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &dqdy = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &dqdz = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    RDouble4D &UVWT = UVWTproxy->GetField_STR();

    Int1D *DiffOrder;
    if (nSurface == 1)
    {
        DiffOrder = grid->DiscretePrecisionKXI;
    }
    else if (nSurface == 2)
    {
        DiffOrder = grid->DiscretePrecisionETA;
    }
    else
    {
        DiffOrder = grid->DiscretePrecisionCTA;
    }

    RDouble4D &dqdxface = dqdxFaceProxy->GetField_STR();
    RDouble4D &dqdyface = dqdyFaceProxy->GetField_STR();
    RDouble4D &dqdzface = dqdzFaceProxy->GetField_STR();
    dqdxface = 0.0;
    dqdyface = 0.0;
    dqdzface = 0.0;   

    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    int il1, jl1, kl1;
    GetNsurfIndex(nSurface, il1, jl1, kl1);

    int ist = 1 - il1;
    int jst = 1 - jl1;
    int kst = 1 - kl1;
    int ied = ni - 1 + 2*il1;
    int jed = nj - 1 + 2*jl1;
    int ked = nk - 1 + 2*kl1;
    
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }
    
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int localcellID = i * il1 + j * jl1 + k * kl1;
                
                if ((*DiffOrder)(localcellID) == 2 || (*DiffOrder)(localcellID - 1) == 2)
                {
                    RDouble dx = xcc(i, j, k) - xcc(i - il1, j - jl1, k - kl1);
                    RDouble dy = ycc(i, j, k) - ycc(i - il1, j - jl1, k - kl1);
                    RDouble dz = zcc(i, j, k) - zcc(i - il1, j - jl1, k - kl1);
                    RDouble dr = sqrt(dx * dx + dy * dy + dz * dz);
                    
                    dx = dx / dr;
                    dy = dy / dr;
                    dz = dz / dr;
                    for (int m = 0; m <= 3; ++ m)
                    {
                        RDouble qxb = dqdx(i - il1,   j - jl1,   k - kl1,   m);
                        RDouble qxc = dqdx(i,         j,         k,         m);
                        
                        RDouble qyb = dqdy(i - il1,   j - jl1,   k - kl1,   m);
                        RDouble qyc = dqdy(i,         j,         k,         m);
                        
                        RDouble qzb = dqdz(i - il1,   j - jl1,   k - kl1,   m);
                        RDouble qzc = dqdz(i,         j,         k,         m);
                        
                        RDouble gradqx = half * (qxb + qxc);
                        RDouble gradqy = half * (qyb + qyc);
                        RDouble gradqz = half * (qzb + qzc);
                        
                        RDouble dq = UVWT(i, j, k, m) - UVWT(i - il1, j - jl1, k - kl1, m);
                        RDouble tmp = gradqx * dx + gradqy * dy + gradqz * dz - dq / dr;
                        
                        dqdxface(i, j, k, m) = gradqx - tmp * dx;
                        dqdyface(i, j, k, m) = gradqy - tmp * dy;
                        dqdzface(i, j, k, m) = gradqz - tmp * dz;
                    }
                }
                else
                {
                    RDouble dx = (27.0 * (xcc(i, j, k) - xcc(i - il1, j - jl1, k - kl1)) - (xcc(i + il1, j + jl1, k + kl1) - xcc(i - il1 * 2, j - jl1 * 2, k - kl1 * 2))) / 24.0;
                    RDouble dy = (27.0 * (ycc(i, j, k) - ycc(i - il1, j - jl1, k - kl1)) - (ycc(i + il1, j + jl1, k + kl1) - ycc(i - il1 * 2, j - jl1 * 2, k - kl1 * 2))) / 24.0;
                    RDouble dz = (27.0 * (zcc(i, j, k) - zcc(i - il1, j - jl1, k - kl1)) - (zcc(i + il1, j + jl1, k + kl1) - zcc(i - il1 * 2, j - jl1 * 2, k - kl1 * 2))) / 24.0;
                    RDouble dr = sqrt(dx * dx + dy * dy + dz * dz);
                    
                    dx = dx / dr;
                    dy = dy / dr;
                    dz = dz / dr;
                    for (int m = 0; m <= 3; ++ m)
                    {
                        RDouble qxa = dqdx(i - il1*2, j - jl1*2, k - kl1*2, m);
                        RDouble qxb = dqdx(i - il1,   j - jl1,   k - kl1,   m);
                        RDouble qxc = dqdx(i,         j,         k,         m);
                        RDouble qxd = dqdx(i + il1,   j + jl1,   k + kl1,   m);
                        
                        RDouble qya = dqdy(i - il1*2, j - jl1*2, k - kl1*2, m);
                        RDouble qyb = dqdy(i - il1,   j - jl1,   k - kl1,   m);
                        RDouble qyc = dqdy(i,         j,         k,         m);
                        RDouble qyd = dqdy(i + il1,   j + jl1,   k + kl1,   m);
                        
                        RDouble qza = dqdz(i - il1*2, j - jl1*2, k - kl1*2, m);
                        RDouble qzb = dqdz(i - il1,   j - jl1,   k - kl1,   m);
                        RDouble qzc = dqdz(i,         j,         k,         m);
                        RDouble qzd = dqdz(i + il1,   j + jl1,   k + kl1,   m);
                        
                        RDouble gradqx = (- qxa + 9.0 * qxb + 9.0 * qxc - qxd) / 16.0;
                        RDouble gradqy = (- qya + 9.0 * qyb + 9.0 * qyc - qyd) / 16.0;
                        RDouble gradqz = (- qza + 9.0 * qzb + 9.0 * qzc - qzd) / 16.0;
                        
                        RDouble dq = (27.0 * (UVWT(i, j, k, m) - UVWT(i - il1, j - jl1, k - kl1, m)) - (UVWT(i + il1, j + jl1, k + kl1,m) - UVWT(i - il1 * 2, j - jl1 * 2, k - kl1 * 2,m))) / 24.0;
                        RDouble tmp = gradqx * dx + gradqy * dy + gradqz * dz - dq / dr;
                        
                        dqdxface(i, j, k, m) = gradqx - tmp * dx;
                        dqdyface(i, j, k, m) = gradqy - tmp * dy;
                        dqdzface(i, j, k, m) = gradqz - tmp * dz;
                    }
                }
            }
        }
    }
}

void NSSolverStructFD::GetGradientAtFaceconservation_HDCS(Grid *gridIn, int nSurface)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();    

    RDouble4D &dqdx = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &dqdy = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &dqdz = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    RDouble4D &UVWT = UVWTproxy->GetField_STR();
    
    RDouble4D &dqdxface = dqdxFaceProxy->GetField_STR();
    RDouble4D &dqdyface = dqdyFaceProxy->GetField_STR();
    RDouble4D &dqdzface = dqdzFaceProxy->GetField_STR();
    dqdxface = 0.0;
    dqdyface = 0.0;
    dqdzface = 0.0;   

    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    int il1, jl1, kl1;
    GetNsurfIndex(nSurface, il1, jl1, kl1);

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1 + il1;
    int jed = nj - 1 + jl1;
    int ked = nk - 1 + kl1;
    
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }
    
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble dx = xcc(i, j, k) - xcc(i - il1, j - jl1, k - kl1);
                RDouble dy = ycc(i, j, k) - ycc(i - il1, j - jl1, k - kl1);
                RDouble dz = zcc(i, j, k) - zcc(i - il1, j - jl1, k - kl1);
                RDouble dr = sqrt(dx * dx + dy * dy + dz * dz);
                
                dx = dx / dr;
                dy = dy / dr;
                dz = dz / dr;
                for (int m = 0; m <= 3; ++ m)
                {
                    RDouble qxb = dqdx(i - il1,   j - jl1,   k - kl1,   m);
                    RDouble qxc = dqdx(i,         j,         k,         m);
                    
                    RDouble qyb = dqdy(i - il1,   j - jl1,   k - kl1,   m);
                    RDouble qyc = dqdy(i,         j,         k,         m);
                    
                    RDouble qzb = dqdz(i - il1,   j - jl1,   k - kl1,   m);
                    RDouble qzc = dqdz(i,         j,         k,         m);
                    
                    RDouble gradqx = half * (qxb + qxc);
                    RDouble gradqy = half * (qyb + qyc);
                    RDouble gradqz = half * (qzb + qzc);
                    
                    RDouble dq = UVWT(i, j, k, m) - UVWT(i - il1, j - jl1, k - kl1, m);
                    RDouble tmp = gradqx * dx + gradqy * dy + gradqz * dz - dq / dr;
                    
                    dqdxface(i, j, k, m) = gradqx - tmp * dx;
                    dqdyface(i, j, k, m) = gradqy - tmp * dy;
                    dqdzface(i, j, k, m) = gradqz - tmp * dz;
                }
            }
        }
    }
}

void NSSolverStructFD::GetGradientAtFacechainrule(Grid *gridIn, int nSurface)
{
    if (solvername.substr(0, 4) == "HDCS")
    {
        GetGradientAtFaceconservation_HDCS(gridIn, nSurface);
        return;
    }
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &xfv  = *(grid->GetFaceVectorX());
    RDouble4D &yfv  = *(grid->GetFaceVectorY());
    RDouble4D &zfv  = *(grid->GetFaceVectorZ());
    RDouble3D &vol  = *(grid->GetCellVolume());

    RDouble4D *dqdxnomaindirection;
    RDouble4D *dqdynomaindirection;
    RDouble4D *dqdznomaindirection;
    Int1D *DiffOrder;
    if (nSurface == 1)
    {
        dqdxnomaindirection = &(*reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdxnokxi")));
        dqdynomaindirection = &(*reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdynokxi")));
        dqdznomaindirection = &(*reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdznokxi")));

        DiffOrder = grid->DiscretePrecisionKXI;
    }
    else if (nSurface == 2)
    {
        dqdxnomaindirection = &(*reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdxnoeta")));
        dqdynomaindirection = &(*reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdynoeta")));
        dqdznomaindirection = &(*reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdznoeta")));

        DiffOrder = grid->DiscretePrecisionETA;
    }
    else
    {
        dqdxnomaindirection = &(*reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdxnocta")));
        dqdynomaindirection = &(*reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdynocta")));
        dqdznomaindirection = &(*reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdznocta")));

        DiffOrder = grid->DiscretePrecisionCTA;
    }

    RDouble4D &UVWT = UVWTproxy->GetField_STR();

    RDouble4D &dqdxface = dqdxFaceProxy->GetField_STR();
    RDouble4D &dqdyface = dqdyFaceProxy->GetField_STR();
    RDouble4D &dqdzface = dqdzFaceProxy->GetField_STR();
    dqdxface = 0.0;
    dqdyface = 0.0;
    dqdzface = 0.0;

    int il1, jl1, kl1;
    GetNsurfIndex(nSurface, il1, jl1, kl1);

    int ist = 1 - il1;
    int jst = 1 - jl1;
    int kst = 1 - kl1;
    int ied = ni - 1 + 2*il1;
    int jed = nj - 1 + 2*jl1;
    int ked = nk - 1 + 2*kl1;
    
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble nx = xfv(i, j, k, nSurface);
                RDouble ny = yfv(i, j, k, nSurface);
                RDouble nz = zfv(i, j, k, nSurface);

                RDouble Volume = half * (vol(i, j, k) + vol(i - il1, j - jl1, k - kl1));

                nx /= Volume;
                ny /= Volume;
                nz /= Volume;

                int localcellID = i * il1 + j * jl1 + k * kl1;
                if ((*DiffOrder)(localcellID) == 2 || (*DiffOrder)(localcellID - 1) == 2)
                {
                    for (int m = 0; m <= 3; ++ m)
                    {
                        RDouble Cellb_q = UVWT(i - il1, j - jl1, k - kl1, m);
                        RDouble Cellc_q = UVWT(i      , j      , k      , m);

                        RDouble Face_dqdmaindirection = Cellc_q - Cellb_q;

                        RDouble Cellb_dqdxnomaindirection = (*dqdxnomaindirection)(i - il1, j - jl1, k - kl1, m);
                        RDouble Cellc_dqdxnomaindirection = (*dqdxnomaindirection)(i      , j      , k      , m);

                        RDouble Cellb_dqdynomaindirection = (*dqdynomaindirection)(i - il1, j - jl1, k - kl1, m);
                        RDouble Cellc_dqdynomaindirection = (*dqdynomaindirection)(i      , j      , k      , m);

                        RDouble Cellb_dqdznomaindirection = (*dqdznomaindirection)(i - il1, j - jl1, k - kl1, m);
                        RDouble Cellc_dqdznomaindirection = (*dqdznomaindirection)(i      , j      , k      , m);

                        RDouble Face_dqdxnomaindirection = half * (Cellb_dqdxnomaindirection + Cellc_dqdxnomaindirection);
                        RDouble Face_dqdynomaindirection = half * (Cellb_dqdynomaindirection + Cellc_dqdynomaindirection);
                        RDouble Face_dqdznomaindirection = half * (Cellb_dqdznomaindirection + Cellc_dqdznomaindirection);

                        dqdxface(i, j, k, m) = nx * Face_dqdmaindirection + Face_dqdxnomaindirection;
                        dqdyface(i, j, k, m) = ny * Face_dqdmaindirection + Face_dqdynomaindirection;
                        dqdzface(i, j, k, m) = nz * Face_dqdmaindirection + Face_dqdznomaindirection;
                    }
                }
                else
                {
                    for (int m = 0; m <= 3; ++ m)
                    {
                        RDouble Cella_q = UVWT(i - il1 * 2, j - jl1 * 2, k - kl1 * 2, m);
                        RDouble Cellb_q = UVWT(i - il1    , j - jl1    , k - kl1    , m);
                        RDouble Cellc_q = UVWT(i          , j          , k          , m);
                        RDouble Celld_q = UVWT(i + il1    , j + jl1    , k + kl1    , m);

                        RDouble Face_dqdmaindirection = (27.0 * (Cellc_q - Cellb_q) - (Celld_q - Cella_q)) / 24.0;

                        RDouble Cella_dqdxnomaindirection = (*dqdxnomaindirection)(i - il1 * 2, j - jl1 * 2, k - kl1 * 2, m);
                        RDouble Cellb_dqdxnomaindirection = (*dqdxnomaindirection)(i - il1    , j - jl1    , k - kl1    , m);
                        RDouble Cellc_dqdxnomaindirection = (*dqdxnomaindirection)(i          , j          , k          , m);
                        RDouble Celld_dqdxnomaindirection = (*dqdxnomaindirection)(i + il1    , j + jl1    , k + kl1    , m);

                        RDouble Cella_dqdynomaindirection = (*dqdynomaindirection)(i - il1 * 2, j - jl1 * 2, k - kl1 * 2, m);
                        RDouble Cellb_dqdynomaindirection = (*dqdynomaindirection)(i - il1    , j - jl1    , k - kl1    , m);
                        RDouble Cellc_dqdynomaindirection = (*dqdynomaindirection)(i          , j          , k          , m);
                        RDouble Celld_dqdynomaindirection = (*dqdynomaindirection)(i + il1    , j + jl1    , k + kl1    , m);

                        RDouble Cella_dqdznomaindirection = (*dqdznomaindirection)(i - il1 * 2, j - jl1 * 2, k - kl1 * 2, m);
                        RDouble Cellb_dqdznomaindirection = (*dqdznomaindirection)(i - il1    , j - jl1    , k - kl1    , m);
                        RDouble Cellc_dqdznomaindirection = (*dqdznomaindirection)(i          , j          , k          , m);
                        RDouble Celld_dqdznomaindirection = (*dqdznomaindirection)(i + il1    , j + jl1    , k + kl1    , m);

                        RDouble Face_dqdxnomaindirection = (9.0 * (Cellb_dqdxnomaindirection + Cellc_dqdxnomaindirection) - (Cella_dqdxnomaindirection + Celld_dqdxnomaindirection)) / 16.0;
                        RDouble Face_dqdynomaindirection = (9.0 * (Cellb_dqdynomaindirection + Cellc_dqdynomaindirection) - (Cella_dqdynomaindirection + Celld_dqdynomaindirection)) / 16.0;
                        RDouble Face_dqdznomaindirection = (9.0 * (Cellb_dqdznomaindirection + Cellc_dqdznomaindirection) - (Cella_dqdznomaindirection + Celld_dqdznomaindirection)) / 16.0;

                        dqdxface(i, j, k, m) = nx * Face_dqdmaindirection + Face_dqdxnomaindirection;
                        dqdyface(i, j, k, m) = ny * Face_dqdmaindirection + Face_dqdynomaindirection;
                        dqdzface(i, j, k, m) = nz * Face_dqdmaindirection + Face_dqdznomaindirection;
                    }
                }
            }
        }
    }
}

void NSSolverStructFD::GetGradientAtFaceCorrectionAtPhysicalBoundary(Grid *gridin, int iSurface)
{
    StructGrid *grid = StructGridCast(gridin);

    RDouble4D &xfv  = *(grid->GetFaceVectorX());
    RDouble4D &yfv  = *(grid->GetFaceVectorY());
    RDouble4D &zfv  = *(grid->GetFaceVectorZ());
    RDouble3D &vol  = *(grid->GetCellVolume());

    RDouble4D &dqdx = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &dqdy = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &dqdz = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    RDouble4D &UVWT = UVWTproxy->GetField_STR();
    
    RDouble4D &dqdxface = dqdxFaceProxy->GetField_STR();
    RDouble4D &dqdyface = dqdyFaceProxy->GetField_STR();
    RDouble4D &dqdzface = dqdzFaceProxy->GetField_STR();
    
    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int BCType = bcregion->GetBCType();
        int nsurf_bc = bcregion->GetFaceDirection() + 1;
        if (IsInterface(BCType) || nsurf_bc != iSurface)
        {
            continue;
        }
        else if (IsWall(BCType))
        {
            int ist, ied, jst, jed, kst, ked;
            bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            int *s_lr3d = bcregion->GetFaceDirectionIndex();
            int id, jd, kd;
            GetBCFaceIDX(s_lr3d, id, jd, kd);
            
            int il1, jl1, kl1;
            GetNsurfIndex(iSurface, il1, jl1, kl1);

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        int ibc1 = i + id;
                        int jbc1 = j + jd;
                        int kbc1 = k + kd;
                        
                        RDouble nx  = xfv (ibc1, jbc1, kbc1, iSurface);
                        RDouble ny  = yfv (ibc1, jbc1, kbc1, iSurface);
                        RDouble nz  = zfv (ibc1, jbc1, kbc1, iSurface);  
                        
                        for (int m = 0; m <= 3; ++ m)
                        {
                            RDouble f_kc_et_ct = UVWT(ibc1, jbc1, kbc1, m) - UVWT(ibc1-il1,jbc1-jl1,kbc1-kl1,m);
                            dqdxface(ibc1, jbc1, kbc1, m) = nx * f_kc_et_ct / vol(i, j, k);
                            dqdyface(ibc1, jbc1, kbc1, m) = ny * f_kc_et_ct / vol(i, j, k);
                            dqdzface(ibc1, jbc1, kbc1, m) = nz * f_kc_et_ct / vol(i, j, k);
                        }
                    }
                }
            }
        }
        else
        {
            int ist, ied, jst, jed, kst, ked;
            bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);
            
            int *s_lr3d = bcregion->GetFaceDirectionIndex();
            int id, jd, kd;
            GetBCFaceIDX(s_lr3d, id, jd, kd); 

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {       
                        int ibc1 = i + id;
                        int jbc1 = j + jd;
                        int kbc1 = k + kd;

                        for (int m = 0; m <= 3; ++ m)
                        {
                            dqdxface(ibc1, jbc1, kbc1, m) = dqdx(i, j, k, m);
                            dqdyface(ibc1, jbc1, kbc1, m) = dqdy(i, j, k, m);
                            dqdzface(ibc1, jbc1, kbc1, m) = dqdz(i, j, k, m);
                        }
                    }
                }
            }
        } 
    }
}

void NSSolverStructFD::GetFaceViscousFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();
    
    RDouble4D * UVWTface = 0;
    if (iSurface == 1)
    {
        UVWTface = &(UVWTfaceIproxy->GetField_STR());
    }
    else if (iSurface == 2)
    {
        UVWTface = &(UVWTfaceJproxy->GetField_STR());
    }
    else
    {
        UVWTface = &(UVWTfaceKproxy->GetField_STR());
    }

    RDouble4D &gradqx = dqdxFaceProxy->GetField_STR();
    RDouble4D &gradqy = dqdyFaceProxy->GetField_STR();
    RDouble4D &gradqz = dqdzFaceProxy->GetField_STR();    

    RDouble4D &flux   = fluxProxy->GetField_STR();

    RDouble4D &xfv  = *(grid->GetFaceVectorX());
    RDouble4D &yfv  = *(grid->GetFaceVectorY());
    RDouble4D &zfv  = *(grid->GetFaceVectorZ());
    
    RDouble3D &vist = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist")); 

    Param_NSSolverStruct *parameters = GetControlParameters();
    RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();
    RDouble oPrandtlLaminar    = parameters->GetoPrandtlLaminar();
    RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble specificHeatAtConstantPressure = 1.0 / ((1.4 - 1.0) * refMachNumber * refMachNumber);

    int il1 = 0, jl1 = 0, kl1 = 0;
    GetNsurfIndex(iSurface, il1, jl1, kl1);

    int ist = 0, jst = 0, kst = 0;
    int ied = 0, jed = 0, ked = 0;

    if (solvername.substr(0, 4) == "WCNS")
    {
        ist = 1 - il1;
        jst = 1 - jl1;
        kst = 1 - kl1;
        ied = ni - 1 + 2*il1;
        jed = nj - 1 + 2*jl1;
        ked = nk - 1 + 2*kl1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }
    }
    else if (solvername.substr(0, 4) == "HDCS")
    {
        ist = 1;
        jst = 1;
        kst = 1;
        ied = ni - 1 + il1;
        jed = nj - 1 + jl1;
        ked = nk - 1 + kl1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }
    }

    RDouble fvis[5] = {0};
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble nx = xfv(i, j, k, iSurface);
                RDouble ny = yfv(i, j, k, iSurface);
                RDouble nz = zfv(i, j, k, iSurface);

                RDouble um  = (*UVWTface)(i, j, k, 0);
                RDouble vm  = (*UVWTface)(i, j, k, 1);
                RDouble wm  = (*UVWTface)(i, j, k, 2);
                RDouble tm  = (*UVWTface)(i, j, k, 3);

                RDouble mul = tm * sqrt(tm) * (1.0 + nonDimensionalSutherlandTemperature) / (tm + nonDimensionalSutherlandTemperature);
                RDouble mut = half * (vist(i-il1, j-jl1, k-kl1) + vist(i, j, k));
                RDouble vis = mul + mut;

                RDouble dudx  = gradqx(i, j, k, 0);
                RDouble dudy  = gradqy(i, j, k, 0);
                RDouble dudz  = gradqz(i, j, k, 0);

                RDouble dvdx  = gradqx(i, j, k, 1);
                RDouble dvdy  = gradqy(i, j, k, 1);
                RDouble dvdz  = gradqz(i, j, k, 1);

                RDouble dwdx  = gradqx(i, j, k, 2);
                RDouble dwdy  = gradqy(i, j, k, 2);
                RDouble dwdz  = gradqz(i, j, k, 2);

                RDouble dtdx  = gradqx(i, j, k, 3);
                RDouble dtdy  = gradqy(i, j, k, 3);
                RDouble dtdz  = gradqz(i, j, k, 3);                
                
                RDouble kcp = (mul * oPrandtlLaminar + mut * oPrandtlTurbulence) * specificHeatAtConstantPressure;

                RDouble qx = kcp * dtdx;
                RDouble qy = kcp * dtdy;
                RDouble qz = kcp * dtdz;

                RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

                RDouble txx = vis * (two * dudx - divv2p3);
                RDouble tyy = vis * (two * dvdy - divv2p3);
                RDouble tzz = vis * (two * dwdz - divv2p3);
                RDouble txy = vis * (dudy + dvdx);
                RDouble txz = vis * (dudz + dwdx);
                RDouble tyz = vis * (dvdz + dwdy);

                fvis[0] = 0.0;
                fvis[1] = nx * txx + ny * txy + nz * txz;
                fvis[2] = nx * txy + ny * tyy + nz * tyz;
                fvis[3] = nx * txz + ny * tyz + nz * tzz;
                fvis[4] = um * fvis[1] + vm * fvis[2] + wm * fvis[3] + (nx * qx + ny * qy + nz * qz);

                for (int m = 0; m < 5; ++ m)
                {
                    flux(i, j, k,m) = - fvis[m] / refReNumber;
                }
            }
        }
    }
}

void NSSolverStructFD::GetFaceViscousFluxLES(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();
    
    RDouble4D * UVWTface = 0;
    if (iSurface == 1)
    {
        UVWTface = &(UVWTfaceIproxy->GetField_STR());
    }
    else if (iSurface == 2)
    {
        UVWTface = &(UVWTfaceJproxy->GetField_STR());
    }
    else
    {
        UVWTface = &(UVWTfaceKproxy->GetField_STR());
    }

    RDouble4D &gradqx = dqdxFaceProxy->GetField_STR();
    RDouble4D &gradqy = dqdyFaceProxy->GetField_STR();
    RDouble4D &gradqz = dqdzFaceProxy->GetField_STR();

    RDouble4D &flux = fluxProxy->GetField_STR();

    RDouble4D &xfv = *(grid->GetFaceVectorX());
    RDouble4D &yfv = *(grid->GetFaceVectorY());
    RDouble4D &zfv = *(grid->GetFaceVectorZ());
    
    RDouble3D &vist = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist")); 
    RDouble3D &subgridScaleEnergy = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));
    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));

    Param_NSSolverStruct *parameters = GetControlParameters();
    RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();
    RDouble oPrandtlLaminar    = parameters->GetoPrandtlLaminar();

    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble specificHeatAtConstantPressure = 1.0 / ((1.4 - 1.0) * refMachNumber * refMachNumber);

    int il1 = 0, jl1 = 0, kl1 = 0;
    GetNsurfIndex(iSurface, il1, jl1, kl1);

    int ist = 0, jst = 0, kst = 0;
    int ied = 0, jed = 0, ked = 0;

    if (solvername.substr(0, 4) == "WCNS")
    {
        ist = 1 - il1;
        jst = 1 - jl1;
        kst = 1 - kl1;
        ied = ni - 1 + 2*il1;
        jed = nj - 1 + 2*jl1;
        ked = nk - 1 + 2*kl1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }
    }
    else if (solvername.substr(0, 4) == "HDCS")
    {
        ist = 1;
        jst = 1;
        kst = 1;
        ied = ni - 1 + il1;
        jed = nj - 1 + jl1;
        ked = nk - 1 + kl1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }
    }

    string sgsmodel = " ";
    GlobalDataBase::GetData("sgsmodel", &sgsmodel, PHSTRING, 1);

    RDouble fvis[5] = {0};

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble nx = xfv(i, j, k, iSurface);
                RDouble ny = yfv(i, j, k, iSurface);
                RDouble nz = zfv(i, j, k, iSurface);

                RDouble um  = (*UVWTface)(i, j, k, 0);
                RDouble vm  = (*UVWTface)(i, j, k, 1);
                RDouble wm  = (*UVWTface)(i, j, k, 2);
                RDouble tm  = (*UVWTface)(i, j, k, 3);

                RDouble mul = tm * sqrt(tm) * (1.0 + nonDimensionalSutherlandTemperature) / (tm + nonDimensionalSutherlandTemperature);
                RDouble mut = half * (vist(i-il1, j-jl1, k-kl1) + vist(i, j, k));
                RDouble vis = mul + mut;

                RDouble dudx  = gradqx(i, j, k, 0);
                RDouble dudy  = gradqy(i, j, k, 0);
                RDouble dudz  = gradqz(i, j, k, 0);

                RDouble dvdx  = gradqx(i, j, k, 1);
                RDouble dvdy  = gradqy(i, j, k, 1);
                RDouble dvdz  = gradqz(i, j, k, 1);

                RDouble dwdx  = gradqx(i, j, k, 2);
                RDouble dwdy  = gradqy(i, j, k, 2);
                RDouble dwdz  = gradqz(i, j, k, 2);

                RDouble dtdx  = gradqx(i, j, k, 3);
                RDouble dtdy  = gradqy(i, j, k, 3);
                RDouble dtdz  = gradqz(i, j, k, 3);                

                RDouble oPrandtlTurbulence = 2.0 / (turbulentPrandtlNumber(i-il1, j-jl1, k-kl1) + turbulentPrandtlNumber(i, j, k));
                RDouble kcp = (mul * oPrandtlLaminar + mut * oPrandtlTurbulence) * specificHeatAtConstantPressure;

                RDouble qx = kcp * dtdx;
                RDouble qy = kcp * dtdy;
                RDouble qz = kcp * dtdz;

                RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

                RDouble kSGS = half * (subgridScaleEnergy(i-il1, j-jl1, k-kl1) + subgridScaleEnergy(i, j, k));

                RDouble txx = vis * (two * dudx - divv2p3) - 1.0/3.0 * kSGS;
                RDouble tyy = vis * (two * dvdy - divv2p3) - 1.0/3.0 * kSGS;
                RDouble tzz = vis * (two * dwdz - divv2p3) - 1.0/3.0 * kSGS;
                RDouble txy = vis * (dudy + dvdx);
                RDouble txz = vis * (dudz + dwdx);
                RDouble tyz = vis * (dvdz + dwdy);

                fvis[0] = 0.0;
                fvis[1] = nx * txx + ny * txy + nz * txz;
                fvis[2] = nx * txy + ny * tyy + nz * tyz;
                fvis[3] = nx * txz + ny * tyz + nz * tzz;
                fvis[4] = um * fvis[1] + vm * fvis[2] + wm * fvis[3] + (nx * qx + ny * qy + nz * qz);

                for (int m = 0; m < 5; ++ m)
                {
                    flux(i, j, k,m) = - fvis[m] / refReNumber;
                }
            }
        }
    }
}

void NSSolverStructFD::SpectrumRadius(Grid *grid)
{
    SpectrumRadiusInviscid(grid);

    SpectrumRadiusViscous(grid);    
}

void NSSolverStructFD::SpectrumRadiusInviscid(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist = 0;
    int ied = ni;
    int jst = 0;
    int jed = nj;
    int kst = 0;
    int ked = nk;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    RDouble4D &xfv = *(grid->GetFaceVectorX());
    RDouble4D &yfv = *(grid->GetFaceVectorY());
    RDouble4D &zfv = *(grid->GetFaceVectorZ());

    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));
    RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("rtem"));
    const RDouble gama = 1.4;

    invSpectralRadius = 0.0;

    using namespace IDX;
    int nDim = GetDim();
    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble rm = q(i, j, k, IR);
                    RDouble um = q(i, j, k, IU);
                    RDouble vm = q(i, j, k, IV);
                    RDouble wm = q(i, j, k, IW);
                    RDouble pm = q(i, j, k, IP);
                    RDouble cm = sqrt(gama * pm / rm);

                    RDouble nx = half * (xfv(i, j, k, iSurface) + xfv(il, jl, kl, iSurface));
                    RDouble ny = half * (yfv(i, j, k, iSurface) + yfv(il, jl, kl, iSurface));
                    RDouble nz = half * (zfv(i, j, k, iSurface) + zfv(il, jl, kl, iSurface));
                    RDouble vn = nx * um + ny * vm + nz * wm;
                    RDouble ns = sqrt(nx * nx + ny * ny + nz * nz);                    

                    RDouble eig = ABS(vn) + cm * ns;
                                        
                    RDouble k1 = 0.25;
                    RDouble k2 = 5.0;
                    RDouble kp = rtem(i, j, k);
                    RDouble eig_lim = eig * (k1 + k2 * kp);
                    
                    if (eig < eig_lim)
                    {
                        invSpectralRadius(i, j, k, iSurface - 1) = 0.5 * (eig * eig + eig_lim * eig_lim) / (eig_lim + 1.E-20);
                    }
                    else
                    {
                        invSpectralRadius(i, j, k, iSurface - 1) = eig;
                    }
                }
            }
        }
    }
}

void NSSolverStructFD::SpectrumRadiusViscous(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
    visSpectralRadius = 0.0;

    Param_NSSolverStruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    if (viscousType == INVISCID)
    {
        return;
    }

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist = 0;
    int ied = ni;
    int jst = 0;
    int jed = nj;
    int kst = 0;
    int ked = nk;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    RDouble4D &xfv  = *(grid->GetFaceVectorX());
    RDouble4D &yfv  = *(grid->GetFaceVectorY());
    RDouble4D &zfv  = *(grid->GetFaceVectorZ());
    RDouble3D &vol  = *(grid->GetCellVolume());

    RDouble4D &q    = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));
    RDouble3D &viscousLaminar    = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));

    RDouble refReNumber       = parameters->GetRefReNumber();
    RDouble prandtlLaminar    = parameters->GetPrandtlLaminar();
    RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();
    
    const RDouble gama = 1.4;
    const RDouble foth = 4.0 / 3.0;

    int nDim = GetDim();    

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble nx = half * (xfv(i, j, k, iSurface) + xfv(il, jl, kl, iSurface));
                    RDouble ny = half * (yfv(i, j, k, iSurface) + yfv(il, jl, kl, iSurface));
                    RDouble nz = half * (zfv(i, j, k, iSurface) + zfv(il, jl, kl, iSurface));

                    RDouble ns2 = nx * nx + ny * ny + nz * nz;

                    visSpectralRadius(i, j, k, iSurface - 1) = ns2;
                }

                RDouble rm   = q(i, j, k, 0);
                RDouble vvol = vol(i, j, k);

                RDouble viscLaminar = viscousLaminar(i, j, k);
                RDouble viscTurbulence = viscousTurbulence(i, j, k);

                RDouble coef;
                RDouble coef1 = foth * (viscLaminar + viscTurbulence);
                RDouble coef2 = gama * (viscLaminar / prandtlLaminar +  viscTurbulence / prandtlTurbulence);
                coef  = MAX(coef1, coef2);

                coef *= 2.0 / (refReNumber * rm * vvol);

                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    visSpectralRadius(i, j, k, iSurface - 1) *= coef;
                }
            }
        }
    }    
}

void NSSolverStructFD::TimeStep(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);   
    RDouble3D &dt   = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));
    RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble cfl = ComputeCFL();
    
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK(); 

    int ist = 1;
    int ied = ni - 1;
    int jst = 1;
    int jed = nj - 1;
    int kst = 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    int nDim = GetDim();

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {    
                RDouble rad_euler = zero;
                RDouble rad_ns    = zero;
                for (int iSurface = 0; iSurface < nDim; ++ iSurface)
                {
                    rad_euler += invSpectralRadius(i, j, k, iSurface);
                    rad_ns    += visSpectralRadius(i, j, k, iSurface);
                }
                RDouble rad = rad_euler + rad_ns;

                dt(i, j, k) = cfl / rad;
            }
        }
    }
}

void NSSolverStructFD::ReduceMaxTimeStep(Grid *gridIn, RDouble globalMinDt)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist = 1;
    int ied = ni - 1;
    int jst = 1;
    int jed = nj - 1;
    int kst = 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    RDouble3D &vol = *grid->GetCellVolume();
    RDouble3D &dt  = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));

    RDouble ktmax = GlobalDataBase::GetDoubleParaFromDB("ktmax");

    if (ktmax > 0.0)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    dt(i, j, k) = MIN(dt(i, j, k), static_cast <RDouble> (ktmax * globalMinDt / vol(i, j, k)));
                }
            }
        }
    }
}

void NSSolverStructFD::Diagonal(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    RDouble3D &dt   = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));
    RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble4D &diagonal = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
    RDouble3D &vol = *(grid->GetCellVolume());

    RDouble dualTimeSpectrumC1 = zero;
    RDouble dualTimeSpectrumC2 = one;

    Param_NSSolverStruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble dualTimeCoefficient[7];
        const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
        dualTimeSpectrumC1 = - dualTimeCoefficient[3];
        dualTimeSpectrumC2 =   dualTimeCoefficient[6];
    }
    
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK(); 

    int ist = 1;
    int ied = ni - 1;
    int jst = 1;
    int jed = nj - 1;
    int kst = 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    int nDim = GetDim();

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {    
                RDouble rad_euler = zero;
                RDouble rad_ns    = zero;
                for (int iSurface = 0; iSurface < nDim; ++ iSurface)
                {
                    rad_euler += invSpectralRadius(i, j, k, iSurface);
                    rad_ns    += visSpectralRadius(i, j, k, iSurface);
                }
                RDouble rad = rad_euler + rad_ns;

                RDouble odt =  dualTimeSpectrumC1 * vol(i, j, k) + dualTimeSpectrumC2 / dt(i, j, k);
                for (int m = 0; m < 5; ++ m)
                {
                    diagonal(i, j, k, m) = odt + rad;
                }
            }
        }
    }
}

void NSSolverStructFD::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    if (tscheme == JACOBIAN_ITERATION)
    {
        SolveLUSGSbyJacobiIteration(gridIn, dqProxy, LUplusDQ, sweepNormal);
        return;
    }
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &dq = dqProxy->GetField_STR();
    RDouble4D &dRHS = LUplusDQ->GetField_STR();
    RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble4D &diagonal = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
    RDouble4D &res = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));    
    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));

    RDouble4D &xfv = *(grid->GetFaceVectorX());
    RDouble4D &yfv = *(grid->GetFaceVectorY());
    RDouble4D &zfv = *(grid->GetFaceVectorZ());

    int *iBlank = grid->GetCellTypeContainer();

    Param_NSSolverStruct *parameters = GetControlParameters();
    bool isViscous = parameters->IsViscous();
    const int nEquation  = 5;

    RDouble df[5];
    RDouble df_total[5];
    RDouble prim_former[5];
    RDouble dq_former[5];    

    int ist = 1;
    int ied = ni - 1;
    int jst = 1;
    int jed = nj - 1;
    int kst = 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }    

    int nDim = GetDim();
    int cellIndex = 0;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                EncodeIJK(cellIndex, i, j, k, ni-1, nj-1);
                
                if (iBlank[cellIndex] != ACTIVE)
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        dq(i, j, k, m) = 0.0;
                    }
                    continue;
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    df_total[m] = 0.0;     
                }

                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i - il1;
                    int jl = j - jl1;
                    int kl = k - kl1;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prim_former[m]  =  q(il, jl, kl, m);
                        dq_former[m]= dq(il, jl, kl, m);
                    }
                   
                    RDouble nx = half * (xfv(i, j, k, iSurface) + xfv(il, jl, kl, iSurface));
                    RDouble ny = half * (yfv(i, j, k, iSurface) + yfv(il, jl, kl, iSurface));
                    RDouble nz = half * (zfv(i, j, k, iSurface) + zfv(il, jl, kl, iSurface));

                    RDouble rad = invSpectralRadius(il, jl, kl, iSurface - 1);

                    GetMXDQ(prim_former, nx, ny, nz, dq_former, df, rad, 1);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        df_total[m] += half * df[m];
                    }

                    if (isViscous)
                    {
                        RDouble rad_vis = half * visSpectralRadius(il, jl, kl, iSurface - 1);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            df_total[m] += rad_vis * dq_former[m];
                        }
                    }
                }                

                for (int m = 0; m < nEquation; ++ m)
                {
                    RDouble dqOld = dq(i, j, k, m);

                    RDouble coef = 0.5;

                    dq(i, j, k, m) = (res(i, j, k, m) + coef * dRHS(i, j, k, m) + df_total[m]) / diagonal(i, j, k, m);

                    dRHS(i, j, k, m) = df_total[m];

                    sweepNormal += SQR(dq(i, j, k, m) - dqOld);
                }
            }
        }
    }
}

void NSSolverStructFD::SolveLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    if (tscheme == JACOBIAN_ITERATION)
    {
        return;
    }
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &dq = dqProxy->GetField_STR();
    RDouble4D &dRHS = LUplusDQ->GetField_STR();
    RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble4D &diagonal = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
    RDouble4D &res = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));    
    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));

    RDouble4D &xfv  = *(grid->GetFaceVectorX());
    RDouble4D &yfv  = *(grid->GetFaceVectorY());
    RDouble4D &zfv  = *(grid->GetFaceVectorZ());

    int *iBlank = grid->GetCellTypeContainer();

    Param_NSSolverStruct *parameters = GetControlParameters();
    bool isViscous = parameters->IsViscous();
    const int nEquation  = 5;

    RDouble df[5];
    RDouble df_total[5];
    RDouble prim_after[5];
    RDouble dq_after[5];    

    int ist = 1;
    int ied = ni - 1;
    int jst = 1;
    int jed = nj - 1;
    int kst = 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }    

    int nDim = GetDim();
    int cellIndex = 0;
    for (int k = ked; k >= kst; -- k)
    {
        for (int j = jed; j >= jst; -- j)
        {
            for (int i = ied; i >= ist; -- i)
            {
                EncodeIJK(cellIndex, i, j, k, ni-1, nj-1);

                if (iBlank[cellIndex] != ACTIVE)
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        dq(i, j, k, m) = 0.0;
                    }
                    continue;
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    df_total[m] = 0.0;
                }

                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prim_after[m] =  q(il, jl, kl, m);
                        dq_after[m]   = dq(il, jl, kl, m);
                    }

                    RDouble nx  = half * (xfv (il, jl, kl, iSurface)  + xfv (il + il1, jl + jl1, kl + kl1, iSurface));
                    RDouble ny  = half * (yfv (il, jl, kl, iSurface)  + yfv (il + il1, jl + jl1, kl + kl1, iSurface));
                    RDouble nz  = half * (zfv (il, jl, kl, iSurface)  + zfv (il + il1, jl + jl1, kl + kl1, iSurface));

                    RDouble rad = invSpectralRadius(il, jl, kl, iSurface - 1);
                    
                    GetMXDQ(prim_after, nx, ny, nz, dq_after, df, rad, -1);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        df_total[m] += - half * df[m];
                    }

                    if (isViscous)
                    {
                        RDouble rad_vis = half * visSpectralRadius(il, jl, kl, iSurface - 1);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            df_total[m] += rad_vis * dq_after[m];
                        }
                    }
                }
                for (int m = 0; m < nEquation; ++ m)
                {
                    RDouble dqOld = dq(i, j, k, m);
                    dq(i, j, k, m) = (res(i, j, k, m) + dRHS(i, j, k, m) + df_total[m]) / diagonal(i, j, k, m);
                    dRHS(i, j, k, m) = df_total[m];
                    sweepNormal += SQR(dq(i, j, k, m) - dqOld);
                }
            }
        }
    }
}

void NSSolverStructFD::SolveLUSGSbyJacobiIteration(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *dqtmpProxy, RDouble &sweepNormal)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &dq    = dqProxy->GetField_STR();
    RDouble4D &dqtmp = dqtmpProxy->GetField_STR();
    RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble4D &diagonal= *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
    RDouble4D &res = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));    
    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));

    RDouble4D &xfv = *(grid->GetFaceVectorX());
    RDouble4D &yfv = *(grid->GetFaceVectorY());
    RDouble4D &zfv = *(grid->GetFaceVectorZ());

    Param_NSSolverStruct *parameters = GetControlParameters();
    bool isViscous = parameters->IsViscous();
    const int nEquation  = 5;

    RDouble df[5];
    RDouble df_total[5];
    RDouble prim_former[5];
    RDouble dq_former[5];
    RDouble prim_after[5];
    RDouble dq_after[5]; 

    int ist = 1;
    int ied = ni - 1;
    int jst = 1;
    int jed = nj - 1;
    int kst = 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }    

    int nDim = GetDim();

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    df_total[m] = 0.0;     
                }

                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i - il1;
                    int jl = j - jl1;
                    int kl = k - kl1;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prim_former[m] = q(il, jl, kl, m);
                        dq_former[m] = dq(il, jl, kl, m);
                    }
                   
                    RDouble nx = half * (xfv(i, j, k, iSurface) + xfv(il, jl, kl, iSurface));
                    RDouble ny = half * (yfv(i, j, k, iSurface) + yfv(il, jl, kl, iSurface));
                    RDouble nz = half * (zfv(i, j, k, iSurface) + zfv(il, jl, kl, iSurface));

                    RDouble rad = invSpectralRadius(il, jl, kl, iSurface - 1);

                    GetMXDQ(prim_former, nx, ny, nz, dq_former, df, rad, 1);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        df_total[m] += half * df[m];
                    }

                    if (isViscous)
                    {
                        RDouble rad_vis = half * visSpectralRadius(il, jl, kl, iSurface - 1);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            df_total[m] += rad_vis * dq_former[m];
                        }
                    }
                }

                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prim_after[m] = q(il, jl, kl, m);
                        dq_after[m] = dq(il, jl, kl, m);
                    }

                    RDouble nx  = half * (xfv (il, jl, kl, iSurface) + xfv (il + il1, jl + jl1, kl + kl1, iSurface));
                    RDouble ny  = half * (yfv (il, jl, kl, iSurface) + yfv (il + il1, jl + jl1, kl + kl1, iSurface));
                    RDouble nz  = half * (zfv (il, jl, kl, iSurface) + zfv (il + il1, jl + jl1, kl + kl1, iSurface));

                    RDouble  rad = invSpectralRadius(il, jl, kl, iSurface - 1);
                    
                    GetMXDQ(prim_after, nx, ny, nz, dq_after, df, rad, -1);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        df_total[m] += - half * df[m];
                    }

                    if (isViscous)
                    {
                        RDouble rad_vis = half * visSpectralRadius(il, jl, kl, iSurface - 1);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            df_total[m] += rad_vis * dq_after[m];
                        }
                    }
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    dqtmp(i, j, k, m) = (res(i, j, k, m) + df_total[m]) / diagonal(i, j, k, m);
                }
            }
        }
    }
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    sweepNormal += SQR(dqtmp(i, j, k, m) - dq(i, j, k, m));

                    dq(i, j, k, m) = dqtmp(i, j, k, m);
                }
            }
        }
    }
}

void NSSolverStructFD::LUSGSInitializationStructHighOrder(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ)
{
    StructGrid *grid = StructGridCast(gridIn);   
    RDouble4D &diagonal = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
    RDouble4D &res = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));    
    RDouble4D &dq = dqProxy->GetField_STR();
    RDouble4D &dRHS = LUplusDQ->GetField_STR();
    
    dRHS = 0.0;
    dq = 0.0;
    //! This statement is necessary.
    //! Otherwise debug failed to run for unknown reason.

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist = 1;
    int ied = ni - 1;
    int jst = 1;
    int jed = nj - 1;
    int kst = 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m < 5; ++ m)
                {
                    dq(i, j, k, m) = res(i, j, k, m) / diagonal(i, j, k, m);
                }
            }
        }
    }
}

void NSSolverStructFD::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, 0, -1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

    RDouble4D &q  = qProxy->GetField_STR();
    RDouble4D &dq = dqProxy->GetField_STR();
    RDouble4D &qnew = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_FD"));

    Param_NSSolverStruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

    using namespace IDX;

    RDouble roo = primitiveVarFarfield[IR];
    RDouble poo = primitiveVarFarfield[IP];
    RDouble densitylimit  = 1.0e-6 * roo;
    RDouble pressurelimit = 1.0e-6 * poo;

    int     nNegativeCell = 0;

    RDouble prim[5];
    RDouble qtry[5];
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m < 5; ++ m)
                {
                    prim[m] = q(i, j, k, m);
                }
                
                Primitive2Conservative(prim, qtry);

                if (dq(i, j, k, 0) > 0.0)
                {
                    qtry[0] += dq(i, j, k, 0);
                }
                else
                {
                    qtry[0] = qtry[0] * exp(dq(i, j, k,0) / (qtry[0] + 1.E-6));
                }

                for (int m = 1; m < 5; ++ m)
                {
                    qtry[m] += dq(i, j, k, m);
                }

                Conservative2Primitive(qtry, prim);
                
                prim[IP] = ABS (prim[IP]);

                if (prim[IR] < densitylimit || prim[IP] < pressurelimit)
                {
                    nNegativeCell += 1;
                    SolutionFix(qProxy, prim, i, j, k);
                }

                for (int m = 0; m < 5; ++ m)
                {
                    qnew(i, j, k, m) = prim[m];
                }
            }
        }
    }

    if (nNegativeCell > 0)
    {
        int nTotalCell = (ied - ist + 1) * (jed - jst + 1) * (ked - kst + 1);
        RDouble negativeRatio = (RDouble)nNegativeCell / (RDouble)nTotalCell;
        cout << "Warning: negative pressure or density appears in " << nNegativeCell << " cells ... " << endl;
        cout << "         ratio = " << negativeRatio << ", CFL auto-reduction may be needed if ratio > 1e-3" << endl;
    }

    ObtainBoundaryValue(grid);
    ObtainGamaAndTemperature(grid);

    RDouble4D &tFD = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t_FD"));
    RDouble4D &q_ordinary = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &t_ordinary = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));

    GetRange(ni, nj, nk, -2, 1, I, J, K);
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                t_ordinary(i, j, k, 0) = tFD(i, j, k, 0);
                for (int m = 0; m < 5; ++ m)
                {
                    q_ordinary(i, j, k, m) = qnew(i, j, k, m);
                }
            }
        }
    }
}

void NSSolverStructFD::SolutionFix(FieldProxy *qProxy, RDouble prim[5], int i, int j, int k)
{
    RDouble4D &q  = qProxy->GetField_STR();

    using namespace IDX;
    Param_NSSolverStruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble roo = primitiveVarFarfield[IR];
    RDouble poo = primitiveVarFarfield[IP];
    RDouble densityLimitLoc = 1.0e-6 * roo;
    RDouble pressureLimitLoc = 1.0e-6 * poo;

    int np = 0;
    int nDim = GetDim();

    for (int iEquation = 0; iEquation < 5; ++ iEquation)
    {
        prim[iEquation] = zero;
    }

    for (int ii = -1; ii <= 1; ++ ii)
    {
        for (int jj = -1; jj <= 1; ++ jj)
        {
            for (int kk = -1; kk <= 1; ++ kk)
            {
                int kk2d = kk;
                if (nDim == 2)
                {
                    kk2d = 0;
                }
                if ((ABS(ii) + ABS(jj) + ABS(kk2d)) == 1)
                {
                    RDouble neighRho = q(i + ii, j + jj, k + kk2d, IR);
                    RDouble neighPre = q(i + ii, j + jj, k + kk2d, IP);
                    if (neighRho < densityLimitLoc || neighPre < pressureLimitLoc)
                    {
                        continue;
                    }
                    np = np + 1;
                    for (int iEquation = 0; iEquation < 5; ++ iEquation)
                    {
                        RDouble f = q(i + ii, j + jj, k + kk2d, iEquation);
                        if (iEquation == IR || iEquation == IP)
                        {
                            f = ABS(f);
                        }
                        prim[iEquation] += f;
                    }
                }
            }
        }
    }

    if (np > 0)
    {
        for (int iEquation = 0; iEquation < 5; ++ iEquation)
        {
            prim[iEquation] /= np;
        }
    }
    else
    {
        prim[IR] = roo;
        prim[IU] = primitiveVarFarfield[IU];
        prim[IV] = primitiveVarFarfield[IV];
        prim[IW] = primitiveVarFarfield[IW];
        prim[IP] = poo;
    }
}

void NSSolverStructFD::GetMXDQ(RDouble *prim, RDouble KX, RDouble KY, RDouble KZ, RDouble *dq, RDouble *f, RDouble radius, int ipn)
{
    const RDouble gama = 1.4;
    const RDouble gama1 = gama - 1;

    RDouble rm = prim[0];
    RDouble um = prim[1];
    RDouble vm = prim[2];
    RDouble wm = prim[3];
    RDouble pm = prim[4];
    
    RDouble Vn = KX * um + KY * vm + KZ * wm;
    RDouble halfV2 = half * (um * um + vm * vm + wm * wm);
    RDouble hm = (gama / gama1)* (pm / rm) + halfV2;
    RDouble lmd1 = Vn + ipn * radius;
    
    RDouble dc = Vn * dq[0] - KX * dq[1] - KY * dq[2] - KZ * dq[3];    
    RDouble dh = gama1 * (halfV2 * dq[0] - (um * dq[1] + vm * dq[2] + wm * dq[3]) + dq[4]);

    f[0] = lmd1 * dq[0] -      dc;
    f[1] = lmd1 * dq[1] - um * dc + KX * dh;
    f[2] = lmd1 * dq[2] - vm * dc + KY * dh;
    f[3] = lmd1 * dq[3] - wm * dc + KZ * dh;
    f[4] = lmd1 * dq[4] - hm * dc + Vn * dh;
}

void NSSolverStructFD::UpdateQlQrOnlyforMixGrid(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    total_negative = 0;

    FieldProxy * qlProxy = new FieldProxy();
    FieldProxy * qrProxy = new FieldProxy();

    int nDim = GetDim();
    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        if (iSurface == 1)
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1")));
        }
        else if (iSurface == 2)
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2")));
        }
        else
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3")));
        }

        GetFaceVariable_HDCS_firstlayerExplicit(grid, qlProxy, qrProxy, iSurface);
        GetFaceVariableCorrectionAtPhysicalBoundary(grid, qlProxy, qrProxy, iSurface);
        GetFaceVariable_HDCS_otherlayerImplicit(grid, qlProxy, qrProxy, iSurface);
    }

    /*if (total_negative > 0)
    {
        std::cout << "warning:    " << total_negative << std::endl;
    }*/

    delete qlProxy;
    delete qrProxy;
}

void NSSolverStructFD::UpdateInviscidfluxOnlyforMixGrid(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    FieldProxy * qlProxy = new FieldProxy();
    FieldProxy * qrProxy = new FieldProxy();
    FieldProxy * fluxProxy = new FieldProxy();

    int nDim = GetDim();
    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        if (iSurface == 1)
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1")));
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxI")));
        }
        else if (iSurface == 2)
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2")));
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxJ")));
        }
        else
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3")));
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxK")));
        }

        GetFaceInviscidFlux(grid, qlProxy, qrProxy, fluxProxy, iSurface);
    }

    delete qlProxy;
    delete qrProxy;
    delete fluxProxy;
}

void NSSolverStructFD::UpdateViscousfluxOnlyforMixGrid(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    Param_NSSolverStruct *parameters = GetControlParameters();

    int viscousType = parameters->GetViscousType();
    if (viscousType <= INVISCID)
    {
        return;
    }

    FieldProxy * fluxProxy = new FieldProxy();

    int nDim = GetDim();
    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        if (iSurface == 1)
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxI")));
        }
        else if (iSurface == 2)
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxJ")));
        }
        else
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxK")));
        }

        GetGradientAtFace(gridIn, iSurface);
        GetGradientAtFaceCorrectionAtPhysicalBoundary(gridIn, iSurface);
        GetFaceViscousFlux(gridIn, fluxProxy, iSurface);
    }
    delete fluxProxy;
}

void NSSolverStructFD::LoadGlobalfacefluxtoResOnlyforMixGrid(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    
    FieldProxy * fluxProxy = new FieldProxy();

    int nDim = GetDim();
    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        if (iSurface == 1)
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxI")));
        }
        else if (iSurface == 2)
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxJ")));
        }
        else
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxK")));
        }

        GetCellInviscidFlux(gridIn, iSurface);
        FluxDifference(grid, fluxProxy, iSurface);
    }

    Param_NSSolverStruct *parameters = GetControlParameters();

    int viscousType = parameters->GetViscousType();
    if (viscousType <= INVISCID)
    {
        delete fluxProxy;
        return;
    }

    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        if (iSurface == 1)
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxI")));
        }
        else if (iSurface == 2)
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxJ")));
        }
        else
        {
            fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxK")));
        }

        GetCellViscousFlux(gridIn, iSurface);
        FluxDifference(gridIn, fluxProxy, iSurface);
    }

    delete fluxProxy;
}


}
