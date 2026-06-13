#include "Geo_StructGrid.h"
#include "MatrixLUSGS.h"
#include "Flux_RoeEntropyCorrection.h"
#include "Gas.h"
#include "Glb_Dimension.h"

namespace PHSPACE
{
     using namespace GAS_SPACE;
//! This is for Matrix LU-SGS according NSMB handbook in 1995 year.
//! This is Matrix LU-SGS method.
void SolveMatrixLUSGSForwardSweep(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    StructGrid *grid = StructGridCast(gridIn);

    SolveDiagonalMatrix(grid);

    RDouble4D &dq = dqProxy->GetField_STR();
    RDouble4D &deltaFlux = LUplusDQ->GetField_STR();

    RDouble4D &gxfn  = *(grid->GetFaceNormalX());
    RDouble4D &gyfn  = *(grid->GetFaceNormalY());
    RDouble4D &gzfn  = *(grid->GetFaceNormalZ());
    RDouble4D &gvgn  = *(grid->GetFaceNormalVelocity());
    RDouble4D &garea = *(grid->GetFaceArea());
    /*RDouble3D &gvol  = *(grid->GetCellVolume());*/

    int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

    RDouble *dqNeighbor = new RDouble[nEquation];
    RDouble *prim   = new RDouble[nEquation];
    RDouble *qNeighbor = new RDouble[nEquation];
    RDouble *rhs0   = new RDouble[nEquation];
    RDouble *df     = new RDouble[nEquation];   
    RDouble *dqCell = new RDouble[nEquation];
    RDouble *dqOld  = new RDouble[nEquation];
    RDouble *temperature = new RDouble[nTemperatureModel];
    RDouble *qPrecondition = new RDouble[nEquation];
    RDouble *Ux = new RDouble [nEquation];
    RDouble *deltaQ = new RDouble[nEquation];
    RDouble *rhsMatrix = new RDouble[nEquation];
    RDouble *priml = new RDouble[nEquation];
    RDouble *primr = new RDouble[nEquation];

    Range M(0, nEquation - 1);
    RDouble2D *diagonalMatrixD = new RDouble2D(M, M, fortranArray);
    RDouble2D &diagonalMatrix2D = *reinterpret_cast<RDouble2D *> (diagonalMatrixD);
    RDouble2D *tempMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D &tempMatrix2D = *reinterpret_cast<RDouble2D *> (tempMatrix);
    RDouble2D *jacobianMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D &jacobianMatrix2D = *reinterpret_cast<RDouble2D *> (jacobianMatrix);

    RDouble4D &res = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
    RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble4D &q  = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &gamma = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));
    RDouble5D &diagonalMatrix = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("diagonalMatrix"));
    RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (gridIn->GetDataPtr("rtem"));

    int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
    grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    RDouble3D viscousLaminar, viscousTurbulence;
    if (viscousType > INVISCID)
    {
        viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
    }

    int nDim = GetDim();

    for (int k = kCellStart; k <= kCellEnd; ++ k)
    {
        for (int j = jCellStart; j <= jCellEnd; ++ j)
        {
            for (int i = iCellStart; i <= iCellEnd; ++ i)
            {    
                for (int m = 0; m < nEquation; ++ m)
                {
                    dq(i, j, k, m) = 0.0;
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    rhs0[m] = 0.0;

                    //! Temporary variable, only because it is convenient for in/out in function.
                    qPrecondition[m] = q(i, j, k, m);

                    //! Back up the old dq, to compute the convergence.
                    //! it is initialized by zero or res in the first sweep,
                    //! then it is updated by backward sweep, that is dq*.
                    dqOld[m] = dq(i, j, k, m);

                    //! Here, the deltaFlux is computed in Upper backward sweep!
                    //! res: b      deltaFlux: Ux, which is computed in backward sweep.
                    //! the dq is not the real dq, now.
                    //! Although the 'dq' changed can not the right one, it dosen't matter, since 
                    //! the following only using the lower neighbor cells, whose 'dq' has been updated.
                    //! It is convenient for precondition transform.
                    dq(i, j, k, m) = res(i, j, k, m);

                    Ux[m] = deltaFlux(i, j, k, m);

                    //! Then reset it to zero to store the Lower forward sweep!
                    deltaFlux(i, j, k, m) = 0.0;

                    //! Temporary variable, only because it is convenient for in/out in precondition function.
                    dqCell[m] = dq(i, j, k, m);
                }

                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1, il, jl, kl;
                    grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

                    grid->GetLeftCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prim[m] = q(il, jl, kl, m);
                        qNeighbor[m] = q(il, jl, kl, m);
                        dqNeighbor[m] = dq(il, jl, kl, m);

                        priml[m] = q(il, jl, kl, m);
                        primr[m] = q(i, j, k, m);
                    }
                   
                    RDouble &xfn1  = gxfn (i ,j ,k ,iSurface);
                    RDouble &yfn1  = gyfn (i ,j ,k ,iSurface);
                    RDouble &zfn1  = gzfn (i ,j ,k ,iSurface);
                    RDouble &area1 = garea(i ,j ,k ,iSurface);
                    RDouble &vgn1  = gvgn (i ,j ,k ,iSurface);

                    RDouble &xfn2  = gxfn (il ,jl ,kl ,iSurface);
                    RDouble &yfn2  = gyfn (il ,jl ,kl ,iSurface);
                    RDouble &zfn2  = gzfn (il ,jl ,kl ,iSurface);
                    RDouble &area2 = garea(il ,jl ,kl ,iSurface);
                    RDouble &vgn2  = gvgn (il ,jl ,kl ,iSurface);

                    RDouble xfn  = half * (xfn1  + xfn2);
                    RDouble yfn  = half * (yfn1  + yfn2);
                    RDouble zfn  = half * (zfn1  + zfn2);
                    RDouble area = half * (area1 + area2);
                    RDouble vgn  = half * (vgn1  + vgn2);

                    RDouble gama = gamma(il, jl, kl);

                    ComputeInviscidMatrix(prim, xfn, yfn, zfn, area, vgn, gama, tempMatrix);

                    RDouble rtemL = rtem(il, jl, kl);
                    RDouble rtemR = rtem(i,  j,  k);

                    SolveRoeMatrix(priml, primr, xfn1, yfn1, zfn1, gama, area1, rtemL, rtemR, jacobianMatrix);
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        for (int n = 0; n < nEquation; ++ n)
                        {
                            tempMatrix2D(m, n) += jacobianMatrix2D(m, n);
                        }
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        RDouble fSum = 0.0;
                        for (int n = 0; n < nEquation; ++ n)
                        {
                            fSum += half * tempMatrix2D(m, n) * dqNeighbor[n];
                        }
                        df[m] = fSum;
                    }

                    //! Add Flux together.
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        rhs0[m] -= df[m];
                    }

                    if (viscousType > INVISCID)
                    {
                        RDouble rad_vis = half * visSpectralRadius(il, jl, kl, iSurface - 1);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            rhs0[m] -= rad_vis * dqNeighbor[m];
                        }
                    }
                }
                
                for (int m = 0; m < nEquation; ++ m)
                {
                    rhsMatrix[m] = 0.0;
                    deltaQ[m]    = 0.0;
                    rhsMatrix[m] = dq(i, j, k, m) - Ux[m] - rhs0[m];
                }
                
                for (int m = 0; m < nEquation; ++ m)
                {
                    for (int n = 0; n < nEquation; ++ n)
                    {
                        diagonalMatrix2D(m, n) = diagonalMatrix(i, j, k, m, n);
                    }
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    RDouble qSum = 0.0;
                    for (int n = 0; n < nEquation; ++ n)
                    {
                        qSum += diagonalMatrix2D(m, n) * rhsMatrix[n];
                    }
                    deltaQ[m] = qSum;
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    //! dq(i, j, k, m) = deltQ[m]
                    dq(i, j, k, m) = deltaQ[m];
                    //! Store the lower forward sweep delta-flux, which will be used in the backward sweep.
                    deltaFlux(i, j, k, m) += rhs0[m];
                    sweepNormal += SQR(dq(i, j, k, m) - dqOld[m]);
                }
            }
        }
    }

    delete [] dqNeighbor;
    delete [] prim;
    delete [] qNeighbor;
    delete [] rhs0;
    delete [] priml;
    delete [] primr;

    delete [] df;
    delete [] deltaQ;
    delete [] rhsMatrix;
    delete diagonalMatrixD;
    delete tempMatrix;
    delete jacobianMatrix;

    delete [] qPrecondition;
    delete [] dqCell;
    delete [] temperature;

    delete [] dqOld;
    delete [] Ux;
}

void SolveMatrixLUSGSBackwardSweep(Grid *gridIn, FieldProxy *dq_proxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &dq = dq_proxy->GetField_STR();
    RDouble4D &deltaFlux = LUplusDQ->GetField_STR();

    RDouble4D &gxfn  = *(grid->GetFaceNormalX());
    RDouble4D &gyfn  = *(grid->GetFaceNormalY());
    RDouble4D &gzfn  = *(grid->GetFaceNormalZ());
    RDouble4D &gvgn  = *(grid->GetFaceNormalVelocity());
    RDouble4D &garea = *(grid->GetFaceArea());
    /*RDouble3D &gvol  = *(grid->GetCellVolume());*/

    int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

    RDouble *dqNeighbor = new RDouble[nEquation];
    RDouble *prim = new RDouble[nEquation];
    RDouble *qNeighbor = new RDouble[nEquation];
    RDouble *dqOld = new RDouble[nEquation];
    RDouble *rhs0  = new RDouble[nEquation];
    RDouble *df    = new RDouble[nEquation];
    RDouble *temperature = new RDouble[nTemperatureModel];
    RDouble *qPrecondition = new RDouble[nEquation];
    RDouble *dqCell = new RDouble[nEquation];
    RDouble *Lx = new RDouble[nEquation];
    RDouble *deltaQ = new RDouble[nEquation];
    RDouble *rhsMatrix = new RDouble[nEquation];
    RDouble *priml = new RDouble[nEquation];
    RDouble *primr = new RDouble[nEquation];
    
    //! Allocate the memories.
    Range M(0, nEquation - 1);
    RDouble2D *diagonalMatrixD = new RDouble2D(M, M, fortranArray);
    RDouble2D &diagonalMatrix2D = *reinterpret_cast<RDouble2D *> (diagonalMatrixD);
    RDouble2D *tempMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D &tempMatrix2D = *reinterpret_cast<RDouble2D *> (tempMatrix);
    RDouble2D *jacobianMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D &jacobianMatrix2D = *reinterpret_cast<RDouble2D *> (jacobianMatrix);

    RDouble4D &res = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
    RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble4D &q  = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &gamma = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));
    RDouble5D &diagonalMatrix = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("diagonalMatrix"));
    RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (gridIn->GetDataPtr("rtem"));

    Range II, JJ, KK;
    GetRange(ni, nj, nk, -2, 1, II, JJ, KK);
    Range MM(0, nEquation - 1);

    int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
    grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    RDouble3D viscousLaminar, viscousTurbulence;
    if (viscousType > INVISCID)
    {
        viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
    }

    int nDim = GetDim();

    for (int k = kCellEnd; k >= kCellStart; -- k)
    {
        for (int j = jCellEnd; j >= jCellStart; -- j)
        {
            for (int i = iCellEnd; i >= iCellStart; -- i)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    dq(i, j, k, m) = 0.0;
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    rhs0[m] = 0.0;
                    dqNeighbor[m] = dq(i, j, k, m);

                    //! Temporary variable, only because it is convenient for in/out in function.
                    qPrecondition[m] = q(i, j, k, m);  

                    //! Back up the old dq, to compute the convergence.
                    //! the 'dq' is dq*, which has been updated in forward.
                    dqOld[m] = dq(i, j, k, m);

                    //! the dq is not the real dq, now.
                    //! it is convenient for precondition transform.
                    dq(i, j, k, m) = res(i, j, k, m);

                    Lx[m] = deltaFlux(i, j, k, m);

                    //! Then reset it to zero to store the upper forward sweep!
                    deltaFlux(i, j, k, m) = 0.0;

                    //! Temporary variable, only because it is convenient for in/out in precondition function.
                    dqCell[m] = dq(i, j, k, m);
                }

                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1, il, jl, kl;
                    grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

                    grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prim[m] = q(il, jl, kl, m);
                        qNeighbor[m] = q(il, jl, kl, m);
                        dqNeighbor[m] = dq(il, jl, kl, m);

                        priml[m] = q(i, j, k, m);
                        primr[m] = q(il, jl, kl, m);
                    }

                    //clz begin
                    RDouble &xfn1  = gxfn (il, jl, kl, iSurface);
                    RDouble &yfn1  = gyfn (il, jl, kl, iSurface);
                    RDouble &zfn1  = gzfn (il, jl, kl, iSurface);
                    RDouble &area1 = garea(il, jl, kl, iSurface);
                    RDouble &vgn1  = gvgn (il, jl, kl, iSurface);

                    RDouble &xfn2  = gxfn (il + il1, jl + jl1, kl + kl1, iSurface);
                    RDouble &yfn2  = gyfn (il + il1, jl + jl1, kl + kl1, iSurface);
                    RDouble &zfn2  = gzfn (il + il1, jl + jl1, kl + kl1, iSurface);
                    RDouble &area2 = garea(il + il1, jl + jl1, kl + kl1, iSurface);
                    RDouble &vgn2  = gvgn (il + il1, jl + jl1, kl + kl1, iSurface);

                    RDouble xfn  = half * (xfn1  + xfn2);
                    RDouble yfn  = half * (yfn1  + yfn2);
                    RDouble zfn  = half * (zfn1  + zfn2);
                    RDouble area = half * (area1 + area2);
                    RDouble vgn  = half * (vgn1  + vgn2);

                    RDouble gama = gamma(il, jl, kl);

                    ComputeInviscidMatrix(prim, xfn, yfn, zfn, area, vgn, gama, tempMatrix);

                    RDouble rtemL = rtem(i, j, k);
                    RDouble rtemR = rtem(il, jl, kl);

                    SolveRoeMatrix(priml, primr, xfn1, yfn1, zfn1, gama, area1, rtemL, rtemR, jacobianMatrix);
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        for (int n = 0; n < nEquation; ++ n)
                        {
                            tempMatrix2D(m, n) -= jacobianMatrix2D(m, n);
                        }
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        RDouble fSum = 0.0;
                        for (int n = 0; n < nEquation; ++ n)
                        {
                            fSum += half * tempMatrix2D(m, n) * dqNeighbor[n];
                        }
                        df[m] = fSum;
                    }

                    //! Add Flux together.
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        rhs0[m] += df[m];
                    }

                    if (viscousType > INVISCID)
                    {
                        RDouble rad_vis = half * visSpectralRadius(il, jl, kl, iSurface - 1);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            rhs0[m] += - rad_vis * dqNeighbor[m];
                        }
                    }
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    rhsMatrix[m] = 0.0;
                    deltaQ[m]    = 0.0;
                    rhsMatrix[m] = dq(i, j, k, m) - Lx[m] - rhs0[m];
                }
                
                for (int m = 0; m < nEquation; ++ m)
                {
                    for (int n = 0; n < nEquation; ++ n)
                    {
                        diagonalMatrix2D(m, n) = diagonalMatrix(i, j, k, m, n);
                    }
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    RDouble qSum = 0.0;
                    for (int n = 0; n < nEquation; ++ n)
                    {
                        qSum += diagonalMatrix2D(m, n) * rhsMatrix[n];
                    }
                    deltaQ[m] = qSum;
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    //! dq(i, j, k, m) = deltQ[m]
                    dq(i, j, k, m) = deltaQ[m];

                    //! Store the upper backward sweep delta-flux, which will be used in the forward sweep.
                    deltaFlux(i, j, k, m) += rhs0[m];
                    sweepNormal += SQR(dq(i, j, k, m) - dqOld[m]);
                }
            }
        }
    }

    delete [] dqNeighbor;
    delete [] dqOld;
    delete [] prim;
    delete [] qNeighbor;
    delete [] rhs0;
    delete [] priml;
    delete [] primr;

    delete [] df;
    delete [] deltaQ;
    delete [] rhsMatrix;
    delete diagonalMatrixD;
    delete tempMatrix;
    delete jacobianMatrix;

    delete [] temperature;
    delete [] Lx;
    delete [] qPrecondition;
    delete [] dqCell;
}

void SolveDiagonalMatrix(StructGrid *gridIn)
{
    int iSimplifyViscousTerm = GlobalDataBase::GetIntParaFromDB("iSimplifyViscousTerm");
    //! Compute the number of equations.
    int nEquation = GlobalDataBase::GetIntParaFromDB("nl");

    RDouble4D &invSpectrumRadius = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("invSpectralRadius"));
    RDouble4D &visSpectrumRadius = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("visSpectralRadius"));
    RDouble4D &q = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("q"));
    RDouble4D &t = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("t"));
    RDouble3D &dt = *reinterpret_cast<RDouble3D *> (gridIn->GetDataPtr("dt"));
    RDouble3D &gamma = *reinterpret_cast<RDouble3D *> (gridIn->GetDataPtr("gama"));

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    RDouble3D &viscousLaminar = *reinterpret_cast<RDouble3D *> (gridIn->GetDataPtr("visl"));
    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (gridIn->GetDataPtr("vist"));

    RDouble4D &gxfn  = *(gridIn->GetFaceNormalX());
    RDouble4D &gyfn  = *(gridIn->GetFaceNormalY());
    RDouble4D &gzfn  = *(gridIn->GetFaceNormalZ());
    RDouble3D &gvol  = *(gridIn->GetCellVolume());
    RDouble4D &gvgn  = *(gridIn->GetFaceNormalVelocity());
    RDouble4D &garea = *(gridIn->GetFaceArea());

    RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (gridIn->GetDataPtr("rtem"));

    Range I, J, K, M(0, nEquation - 1);
    RDouble5D &diagonalMatrix = *reinterpret_cast<RDouble5D *> (gridIn->GetDataPtr("diagonalMatrix"));
    RDouble2D *tempMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D &tempMatrix2D = *reinterpret_cast<RDouble2D *> (tempMatrix);
    RDouble2D *diagonalMatrixD = new RDouble2D(M, M, fortranArray);
    RDouble2D &diagonalMatrix2D = *reinterpret_cast<RDouble2D *> (diagonalMatrixD);
    RDouble2D *inverseDiaMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D &inverseDiaMatrix2D = *reinterpret_cast<RDouble2D *> (inverseDiaMatrix);

    RDouble *primitiveVars = new RDouble[nEquation];
    RDouble *temperature = new RDouble[nEquation];
    RDouble *priml = new RDouble[nEquation];
    RDouble *primr = new RDouble[nEquation];
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");

    int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
    gridIn->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

    int nDim = GetDim();
    int iL, jL, kL, iR, jR, kR;
    RDouble gama, visLaminar, visTurbulence;
    RDouble xfn1, yfn1, zfn1, area1, vgn1;
    RDouble xfn2, yfn2, zfn2, area2, vgn2;
    RDouble rad_euler, rad_ns;

    //! To Compute diagonal matrix.
    for (int k = kCellStart; k <= kCellEnd; ++ k)
    {
        for (int j = jCellStart; j <= jCellEnd; ++ j)
        {
            for (int i = iCellStart; i <= iCellEnd; ++ i)
            {
                //! Obtain the primitive variables.
                for (int m = 0; m < nEquation; ++ m)
                {
                    primitiveVars[m] = q(i, j, k, m);
                }
                gama = gamma(i, j, k);
                if (viscousType > INVISCID)
                {
                    visLaminar = viscousLaminar(i, j, k);
                    visTurbulence = viscousTurbulence(i, j, k);
                }
                else
                {
                    visLaminar = 0.0;
                    visTurbulence = 0.0;
                }
                for (int m = 0; m < nTemperatureModel; ++ m)
                {
                    temperature[m] = t(i, j, k, m);
                }

                // initialize.
                for (int m = 0; m < nEquation; ++ m)
                {
                    for (int n = 0; n < nEquation; ++ n)
                    {
                        diagonalMatrix(i, j, k, m, n) = 0.0;
                    }
                }

                //! Compute the scalar coefficient of the diagonal matrix.
                for (int m = 0; m < nEquation; ++ m)
                {
                    diagonalMatrix(i, j, k, m, m) = 1.0 / dt(i, j, k) + 1.0 * gvol(i, j, k);
                }

                rad_euler = 0.0;
                rad_ns = 0.0;
                //! Construct the diagonal matrix called D.
                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    gridIn->GetNsurfIndex(il1, jl1, kl1, iSurface);

                    //! Compute spectrum radius of inviscid flux.
                    if (iSimplifyViscousTerm == 1)// Equivalent method to computing the spectrum radius.
                    {
                        rad_euler += invSpectrumRadius(i, j, k, iSurface - 1);
                        rad_ns += visSpectrumRadius(i, j, k, iSurface - 1);
                    }

                    //! Compute the surface -1/2
                    gridIn->GetLeftCellOfFace(i, j, k, il1, jl1, kl1, iL, jL, kL);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        priml[m]  =  q(iL, jL, kL, m);
                        primr[m]  =  q(i, j, k, m);
                    }

                    xfn1  = gxfn (i ,j ,k ,iSurface);
                    yfn1  = gyfn (i ,j ,k ,iSurface);
                    zfn1  = gzfn (i ,j ,k ,iSurface);
                    area1 = garea(i ,j ,k ,iSurface);
                    vgn1  = gvgn (i ,j ,k ,iSurface);

                    RDouble rtemL = rtem(iL, jL, kL);
                    RDouble rtemR = rtem(i,  j,  k);

                    SolveRoeMatrix(priml, primr, xfn1, yfn1, zfn1, gama, area1, rtemL, rtemR, tempMatrix);
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        for (int n = 0; n < nEquation; ++ n)
                        {
                            diagonalMatrix(i, j, k, m, n) += 0.5 * tempMatrix2D(m, n);
                        }
                    }

                    //! Compute the surface +1/2
                    gridIn->GetRightCellOfFace(i, j, k, il1, jl1, kl1, iR, jR, kR);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        priml[m]  =  q(i, j, k, m);
                        primr[m]  =  q(iR, jR, kR, m);
                    }

                    xfn2  = gxfn (iR, jR, kR, iSurface);
                    yfn2  = gyfn (iR, jR, kR, iSurface);
                    zfn2  = gzfn (iR, jR, kR, iSurface);
                    area2 = garea(iR, jR, kR, iSurface);
                    vgn2  = gvgn (iR, jR, kR, iSurface);

                    rtemL = rtem(i, j, k);
                    rtemR = rtem(iR, jR, kR);

                    SolveRoeMatrix(priml, primr, xfn2, yfn2, zfn2, gama, area2, rtemL, rtemR, tempMatrix);
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        for (int n = 0; n < nEquation; ++ n)
                        {
                            diagonalMatrix(i, j, k, m, n) += 0.5 * tempMatrix2D(m, n);
                        }
                    }
                }

                //! Compute the summation of spectrum radius.
                for (int m = 0; m < nEquation; ++ m)
                {
                    diagonalMatrix(i, j, k, m, m) += rad_ns;
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    for (int n = 0; n < nEquation; ++ n)
                    {
                        diagonalMatrix2D(m, n) = diagonalMatrix(i, j, k, m, n);
                    }
                }

                //! ComputeInverseDiagonalMatrix(diagonalMatrixD, nEquation, inverseDiaMatrix);
                ObtainInverseDiagonalMatrix(diagonalMatrixD, nEquation, inverseDiaMatrix);

                for (int m = 0; m < nEquation; ++ m)
                {
                    for (int n = 0; n < nEquation; ++ n)
                    {
                        diagonalMatrix(i, j, k, m, n) = inverseDiaMatrix2D(m, n);
                    }
                }
            }
        }
    }

    delete [] primitiveVars;
    delete [] temperature;
    delete tempMatrix;
    delete diagonalMatrixD;
    delete inverseDiaMatrix;
    delete [] priml;
    delete [] primr;
}

void ComputeInviscidMatrix(RDouble *primitiveVars, RDouble xfn, RDouble yfn, RDouble zfn, RDouble area, RDouble vgn, RDouble gama, RDouble2D *tempMatrix)
{
    int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    //! Compute the number of equations.
    int neqn = nLaminar + nChemical + nTemperatureModel - 1;

    //! To store the Jacobian matrix of inviscid flux.
    RDouble2D &inviscidMatrix = *reinterpret_cast<RDouble2D *> (tempMatrix);

    using namespace IDX;
    //! Obtain the primitive variables.
    RDouble uVelocity, vVelocity, wVelocity, density, pressure, averageVelocity, rVelocity;
    density   = primitiveVars[IR];
    pressure  = primitiveVars[IP];
    uVelocity = primitiveVars[IU];
    vVelocity = primitiveVars[IV];
    wVelocity = primitiveVars[IW];
    averageVelocity = uVelocity * xfn + vVelocity * yfn + wVelocity * zfn;
    //! The relative velocity.
    rVelocity = averageVelocity - vgn;

    //! To compute the total enthalpy.
    RDouble squareVelocity, totalEnthalpy=0.0;
    squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;
    if (nChemical == 0)
    {
        totalEnthalpy = (gama / (gama - 1)) * pressure / density;
    }
    
    totalEnthalpy += 0.5 * squareVelocity;

    //! Construct the Jacobian matrix.
    RDouble alpha = 0.5 * (gama - 1.0) * squareVelocity;
    inviscidMatrix(0, 0) = -vgn;
    inviscidMatrix(1, 0) = alpha * xfn - uVelocity * averageVelocity;
    inviscidMatrix(2, 0) = alpha * yfn - vVelocity * averageVelocity;
    inviscidMatrix(3, 0) = alpha * zfn - wVelocity * averageVelocity;
    inviscidMatrix(4, 0) = (alpha - totalEnthalpy) * averageVelocity;

    //! The second column.
    inviscidMatrix(0, 1) = xfn;
    inviscidMatrix(1, 1) = (2.0 - gama) * uVelocity * xfn + rVelocity;
    inviscidMatrix(2, 1) = (1.0 - gama) * uVelocity * yfn + vVelocity * xfn;
    inviscidMatrix(3, 1) = (1.0 - gama) * uVelocity * zfn + wVelocity * xfn;
    inviscidMatrix(4, 1) = (1.0 - gama) * uVelocity * averageVelocity + totalEnthalpy * xfn;

    //! The third column.
    inviscidMatrix(0, 2) = yfn;
    inviscidMatrix(1, 2) = (1.0 - gama) * vVelocity * xfn + uVelocity * yfn;
    inviscidMatrix(2, 2) = (2.0 - gama) * vVelocity * yfn + rVelocity;
    inviscidMatrix(3, 2) = (1.0 - gama) * vVelocity * zfn + wVelocity * yfn;
    inviscidMatrix(4, 2) = (1.0 - gama) * vVelocity * averageVelocity + totalEnthalpy * yfn;

    //! The fourth column.
    inviscidMatrix(0, 3) = zfn;
    inviscidMatrix(1, 3) = (1.0 - gama) * wVelocity * xfn + uVelocity * zfn;
    inviscidMatrix(2, 3) = (1.0 - gama) * wVelocity * yfn + vVelocity * zfn;
    inviscidMatrix(3, 3) = (2.0 - gama) * wVelocity * zfn + rVelocity;
    inviscidMatrix(4, 3) = (1.0 - gama) * wVelocity * averageVelocity + totalEnthalpy * zfn;

    //! The fifth column.
    inviscidMatrix(0, 4) = 0.0;
    inviscidMatrix(1, 4) = (gama - 1.0) * xfn;
    inviscidMatrix(2, 4) = (gama - 1.0) * yfn;
    inviscidMatrix(3, 4) = (gama - 1.0) * zfn;
    inviscidMatrix(4, 4) = (gama - 1.0) * averageVelocity + rVelocity;

    //! To multiply the average area located in the center of the cell.
    for (int m = 0; m < neqn; ++ m)
    {
        for (int n = 0; n < neqn; ++ n)
        {
            inviscidMatrix(m, n) *= area;
        }
    }
}

void SolveRoeMatrix(RDouble *priml, RDouble *primr, RDouble xfn, RDouble yfn, RDouble zfn, RDouble gama, RDouble area, RDouble rtemL, RDouble rtemR, RDouble2D *jacobianMatrix)
{
    using namespace IDX;
    RDouble2D &matrix = *reinterpret_cast<RDouble2D *> (jacobianMatrix);

    RDouble rl = priml[IDX::IR];
    RDouble ul = priml[IDX::IU];
    RDouble vl = priml[IDX::IV];
    RDouble wl = priml[IDX::IW];
    RDouble pl = priml[IDX::IP];

    RDouble rr = primr[IDX::IR];
    RDouble ur = primr[IDX::IU];
    RDouble vr = primr[IDX::IV];
    RDouble wr = primr[IDX::IW];
    RDouble pr = primr[IDX::IP];

    RDouble nx = xfn;
    RDouble ny = yfn;
    RDouble nz = zfn;
    RDouble t0 = gama - 1.0;

    if (rl < 1.0e-5) rl = 1.0e-5;
    if (rr < 1.0e-5) rr = 1.0e-5;

    //! Note that not to miss the gama.
    RDouble hl = gama * pl / rl / t0 + half * (ul * ul + vl * vl + wl * wl);
    RDouble hr = gama * pr / rr / t0 + half * (ur * ur + vr * vr + wr * wr);

    RDouble ratio = sqrt(rr / rl);
    RDouble coef = 1.0 / (1.0 + ratio);

    RDouble u = (ul + ur * ratio) * coef;
    RDouble v = (vl + vr * ratio) * coef;
    RDouble w = (wl + wr * ratio) * coef;
    RDouble h = (hl + hr * ratio) * coef;
    RDouble halfV2 = half * (u * u + v * v + w * w);
    RDouble Vn = nx * u + ny * v + nz * w;
    RDouble absVel = sqrt(two * halfV2);

    RDouble c2 = t0 * (h - halfV2);
    RDouble c = sqrt(ABS(c2));

    RDouble eigv1 = ABS(Vn);
    RDouble eigv2 = ABS(Vn + c);
    RDouble eigv3 = ABS(Vn - c);

    //! Entropy fix.
    int RoeEntropyFixMethod = 2;
    RDouble RoeEntropyFixCoef1 = 0.0;
    RDouble RoeEntropyFixCoef2 = 0.0;

    Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, rtemL, rtemR, absVel, Vn, c, RoeEntropyFixMethod, RoeEntropyFixCoef1, RoeEntropyFixCoef2);

    RDouble y0 = eigv1;
    RDouble y1 = half * (2.0 * eigv1 - eigv2 - eigv3) / c2;
    RDouble y2 = half * (eigv2 - eigv3) / c;

    //! First row.
    matrix(0, 0) = y0 - halfV2 * t0 * y1 - Vn * y2;
    matrix(0, 1) = u * t0 * y1 + nx * y2;
    matrix(0, 2) = v * t0 * y1 + ny * y2;
    matrix(0, 3) = w * t0 * y1 + nz * y2;
    matrix(0, 4) = -t0 * y1;

    //! Second row.
    matrix(1, 0) = (nx * Vn * c2 - u * halfV2 * t0) * y1 + (nx * halfV2 * t0 - u * Vn) * y2;
    matrix(1, 1) = (u * u * t0 - c2 * nx * nx) * y1 + (nx * u - nx * u * t0) * y2 + y0;
    matrix(1, 2) = (v * u * t0 - c2 * ny * nx) * y1 + (ny * u - nx * v * t0) * y2;
    matrix(1, 3) = (w * u * t0 - c2 * nz * nx) * y1 + (nz * u - nx * w * t0) * y2;
    matrix(1, 4) = -u * t0 * y1 + nx * t0 * y2;

    //! Third row.
    matrix(2, 0) = (ny * Vn * c2 - v * halfV2 * t0) * y1 + (ny * halfV2 * t0 - v * Vn) * y2;
    matrix(2, 1) = (u * v * t0 - c2 * nx * ny) * y1 + (nx * v - ny * u * t0) * y2;
    matrix(2, 2) = (v * v * t0 - c2 * ny * ny) * y1 + (ny * v - ny * v * t0) * y2 + y0;
    matrix(2, 3) = (w * v * t0 - c2 * nz * ny) * y1 + (nz * v - ny * w * t0) * y2;
    matrix(2, 4) = -v * t0 * y1 + ny * t0 * y2;

    //! Fourth row.
    matrix(3, 0) = (nz * Vn * c2 - w * halfV2 * t0) * y1 + (nz * halfV2 * t0 - w * Vn) * y2;
    matrix(3, 1) = (u * w * t0 - c2 * nx * nz) * y1 + (nx * w - nz * u * t0) * y2;
    matrix(3, 2) = (v * w * t0 - c2 * ny * nz) * y1 + (ny * w - nz * v * t0) * y2;
    matrix(3, 3) = (w * w * t0 - c2 * nz * nz) * y1 + (nz * w - nz * w * t0) * y2 + y0;
    matrix(3, 4) = -w * t0 * y1 + nz * t0 * y2;

    //! Fi(t,r)w.
    matrix(4, 0) = (-halfV2 * halfV2 * t0 + Vn * Vn * c2 - halfV2 * c2) * y1 + (halfV2 * Vn * (t0 - 1.0) - c2 * Vn / t0) * y2;
    matrix(4, 1) = (u * c2 + u * halfV2 * t0 - nx * Vn * c2) * y1 + (nx * halfV2 + nx * c2 / t0 - u * Vn * t0) * y2;
    matrix(4, 2) = (v * c2 + v * halfV2 * t0 - ny * Vn * c2) * y1 + (ny * halfV2 + ny * c2 / t0 - v * Vn * t0) * y2;
    matrix(4, 3) = (w * c2 + w * halfV2 * t0 - nz * Vn * c2) * y1 + (nz * halfV2 + nz * c2 / t0 - w * Vn * t0) * y2;
    matrix(4, 4) = -halfV2 * t0 * y1 + Vn * t0 * y2 + half * (eigv2 + eigv3);

    for (int m = 0; m < 5; ++ m)
    {
        for (int n = 0; n < 5; ++ n)
        {
            matrix(m, n) *= area;
        }
    }
}

void ComputeViscousMatrix(RDouble *primitiveVars, RDouble *temperature, RDouble xfn, RDouble yfn, RDouble zfn, RDouble area, RDouble vol, 
                          RDouble vgn, RDouble gama, RDouble laminarViscosity, RDouble turbulentViscosity, RDouble2D *tempMatrix)
{
    int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    //! Compute the number of equations.
    int neqn = nLaminar + nChemical + nTemperatureModel - 1;

    //! To store the Jacobian matrix of inviscid flux.
    RDouble2D &viscousMatrix = *reinterpret_cast<RDouble2D *> (tempMatrix);

    using namespace IDX;
    //! Obtain the primary variables.
    RDouble uVelocity, vVelocity, wVelocity, density, pressure, averageVelocity;
    density   = primitiveVars[IR];
    pressure  = primitiveVars[IP];
    uVelocity = primitiveVars[IU];
    vVelocity = primitiveVars[IV];
    wVelocity = primitiveVars[IW];
    averageVelocity = uVelocity * xfn + vVelocity * yfn + wVelocity * zfn;

    //! trTemperature is the translation-rotation temperature, vTemperature and eTemperature denote the temperatures \n
    //! of the vibration mode and electron mode, respectively.
    RDouble trTemperature, vTemperature, eTemperature;

    //! To obtain the value of temperature.
    trTemperature = temperature[ITT];
    if (nTemperatureModel > 1)
    {
        vTemperature = temperature[ITV];
        eTemperature = vTemperature;
        if (nTemperatureModel > 2)
        {
            eTemperature = temperature[ITE];
        }
    }
    else
    {
        vTemperature = trTemperature;
        eTemperature = trTemperature;
    }

    //! ktr is the heat conductivity of translation and rotation mode, kv and ke are heat conductivities \n
    //! of vibration mode and electron mode, respectively. kcp denotes the total heat conductivity of mixed gas.
    RDouble ktr, kv, ke, kcp;
    RDouble oprl = GlobalDataBase::GetDoubleParaFromDB("oprl");
    RDouble oprt = GlobalDataBase::GetDoubleParaFromDB("oprt");

    //! Compute the heat conductivity.
    if (nChemical == 0) //The heat conductivity of perfect gas
    {
        RDouble R = gas->GetUniversalGasConstant();
        ktr = (laminarViscosity * oprl + turbulentViscosity * oprt) * gama / (gama - 1.0) * R;
        kv = 0.0;
        ke = 0.0;
    }
    else
    {
        gas->ComputeMixtureGasHeatConductivity(primitiveVars, trTemperature, vTemperature, eTemperature, ktr, kv, ke);
    }
    //! To compute the total heat conductivity and viscosity.
    kcp = ktr + kv + ke;
    RDouble viscosity = laminarViscosity + turbulentViscosity;

    //! To compute the square velocity.
    RDouble squareVelocity;
    squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;

    //! Construct the Jacobian matrix.
    RDouble alpha = viscosity / density;
    RDouble phi = kcp * (gama - 1.0) * pressure / trTemperature;
    viscousMatrix(0, 0) = 0.0;
    viscousMatrix(1, 0) = -alpha * (uVelocity + 1.0/3.0 * averageVelocity * xfn);
    viscousMatrix(2, 0) = -alpha * (vVelocity + 1.0/3.0 * averageVelocity * yfn);
    viscousMatrix(3, 0) = -alpha * (wVelocity + 1.0/3.0 * averageVelocity * zfn);
    viscousMatrix(4, 0) = -kcp * trTemperature / density + 0.5 * phi * squareVelocity
                          -2.0 * alpha * (squareVelocity + 1.0/3.0 * averageVelocity * averageVelocity);

    //! The second column.
    viscousMatrix(0, 1) = 0.0;
    viscousMatrix(1, 1) = alpha * (1.0 + 1.0/3.0 * xfn * xfn);
    viscousMatrix(2, 1) = alpha * (0.0 + 1.0/3.0 * xfn * yfn);
    viscousMatrix(3, 1) = alpha * (0.0 + 1.0/3.0 * xfn * zfn);
    viscousMatrix(4, 1) = 2.0 * alpha * (uVelocity + 1.0/3.0 * averageVelocity * xfn) - phi * uVelocity;

    // The third column.
    viscousMatrix(0, 2) = 0.0;
    viscousMatrix(1, 2) = alpha * (0.0 + 1.0/3.0 * xfn * yfn);
    viscousMatrix(2, 2) = alpha * (1.0 + 1.0/3.0 * yfn * yfn);
    viscousMatrix(3, 2) = alpha * (0.0 + 1.0/3.0 * yfn * zfn);
    viscousMatrix(4, 2) = 2.0 * alpha * (vVelocity + 1.0/3.0 * averageVelocity * yfn) - phi * vVelocity;

    //! The fourth column.
    viscousMatrix(0, 3) = 0.0;
    viscousMatrix(1, 3) = alpha * (0.0 + 1.0/3.0 * xfn * zfn);
    viscousMatrix(2, 3) = alpha * (0.0 + 1.0/3.0 * yfn * zfn);
    viscousMatrix(3, 3) = alpha * (1.0 + 1.0/3.0 * zfn * zfn);
    viscousMatrix(4, 3) = 2.0 * alpha * (wVelocity + 1.0/3.0 * averageVelocity * zfn) - phi * wVelocity;

    //! The fifth column.
    viscousMatrix(0, 4) = 0.0;
    viscousMatrix(1, 4) = 0.0;
    viscousMatrix(2, 4) = 0.0;
    viscousMatrix(3, 4) = 0.0;
    viscousMatrix(4, 4) = phi;

    //! To multiply the average area located in the center of the cell.
    for (int m = 0; m < neqn; ++ m)
    {
        for (int n = 0; n < neqn; ++ n)
        {
            viscousMatrix(m, n) *= area * area / (vol * refReNumber);
        }
    }
}

void ObtainInverseDiagonalMatrix(RDouble2D *diagonalMatrixD, int nEquation, RDouble2D *inverseDiaMatrix)
{
    Range II(0, nEquation - 1);
    RDouble2D *lMatrix = new RDouble2D(II, II, fortranArray);
    RDouble2D &lMatrix2D = *reinterpret_cast<RDouble2D *> (lMatrix);

    RDouble2D *uMatrix = new RDouble2D(II, II, fortranArray);
    RDouble2D &uMatrix2D = *reinterpret_cast<RDouble2D *> (uMatrix);

    RDouble2D *lInverseMatrix = new RDouble2D(II, II, fortranArray);
    RDouble2D &lInverseMatrix2D = *reinterpret_cast<RDouble2D *> (lInverseMatrix);

    RDouble2D *uInverseMatrix = new RDouble2D(II, II, fortranArray);
    RDouble2D &uInverseMatrix2D = *reinterpret_cast<RDouble2D *> (uInverseMatrix);

    RDouble2D &dMatrix = *reinterpret_cast<RDouble2D *>(diagonalMatrixD);
    RDouble2D &inverseMatrix = *reinterpret_cast<RDouble2D *>(inverseDiaMatrix);

    int i, j, k;
    RDouble s, t;

    for (i = 0; i < nEquation; i ++)
    {
        for (j = 0; j < nEquation; j ++)
        {
            lMatrix2D(i, j) = 0.0;
            uMatrix2D(i, j) = 0.0;
            lInverseMatrix2D(i, j) = 0.0;
            uInverseMatrix2D(i, j) = 0.0;
        }
    }

    //! Calculate the diagonal of the l matrix.
    for (i = 0; i < nEquation; i ++)
    {
        lMatrix2D(i, i) = 1;
    }

    for (i = 0; i < nEquation; i ++)
    {
        for (j = i; j < nEquation; j ++)
        {
            s = 0.0;
            for (k = 0; k < i; k ++)
            {
                s += lMatrix2D(i, k) * uMatrix2D(k, j);
            }
            uMatrix2D(i, j) = dMatrix(i, j) - s;    //! Calculate u value by row.
        }

        for (j = i + 1; j < nEquation; j ++)
        {
            s = 0.0;
            for (k = 0; k < i; k ++)
            {
                s += lMatrix2D(j, k) * uMatrix2D(k, i);
            }
            lMatrix2D(j, i) = (dMatrix(j ,i) - s) / uMatrix2D(i, i);    //! Calculate l value by column.
        }
    }

    //! Calculate the inverse matrix of l in row order from top to bottom.
    for (i = 0; i < nEquation; i ++)
    {
        lInverseMatrix2D(i, i) = 1;
    }

    for (i = 1; i < nEquation; i ++)
    {
        for (j = 0; j < i; j ++)
        {
            s = 0.0;
            for (k = 0; k < i; k ++)
            {
                s += lMatrix2D(i, k) * lInverseMatrix2D(k, j);
            }
            lInverseMatrix2D(k, j) = -s;
        }
    }

    //! Calculate the inverse matrix of u in column order from bottom to top.
    for (i = 0; i < nEquation; i ++)
    {
        uInverseMatrix2D(i, i) = 1 / uMatrix2D(i, i);
    }

    for (i = 1; i < nEquation; i ++)
    {
        for (j = i - 1; j >= 0; j --)
        {
            s = 0.0;
            for (k = j + 1; k <= i; k ++)
            {
                s += uMatrix2D(j, k) * uInverseMatrix2D(k, i);
            }
            uInverseMatrix2D(j, i) = -s / uMatrix2D(j, j);
        }
    }

    for (i = 0; i < nEquation; i ++)
    {
        for (j = 0; j < nEquation; j ++)
        {
            t = 0.0;
            for (k = 0; k < nEquation; k ++)
            {
                t += uInverseMatrix2D(i, k) * lInverseMatrix2D(k, j);
            }
            inverseMatrix(i, j) = t;
        }
    }
    delete lMatrix;
    delete uMatrix;
    delete lInverseMatrix;
    delete uInverseMatrix;
}

}