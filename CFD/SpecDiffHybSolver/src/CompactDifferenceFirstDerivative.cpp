#include "DataStruct_Array.h"
#include "GlobalDataBase.h"
#include "TK_Exit.h"
#include "CompactDifferenceFirstDerivative.h"
#include "ExplicitDifferenceBoundary.h"
#include "PHMpi.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include "PHHeader.h"

using namespace std;

namespace PHSPACE
{
CompactDifferenceFirstDerivative::CompactDifferenceFirstDerivative()
{
}

void CompactDifferenceFirstDerivative::InitCompactDifferenceFirstDerivativeData(SpecDiffHybGrid *GridData, ExplicitDifferenceBoundary *DiffBoundaryData)
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex  ;

    int ischeme = GlobalDataBase::GetIntParaFromDB("ischeme");
    if ( ischeme == 1 ) 
    {
        SetCompactDiffData_CarpenterC6(ist, ied);
    }
    else if(ischeme == 3)
    {
        SetCompactDiffData_optC4_InnerOnly(ist, ied, DiffBoundaryData);
    }
    else
    {
        TK_Exit::UnexpectedVarValue( "ischeme", ischeme );
    }
}

void CompactDifferenceFirstDerivative::SetCompactDiffData_CarpenterC6(const int ist, const int ied)
{
    minofLHSMatrixBound = -1; //<--- 左侧模板，左边界最小值  
    maxofLHSMatrixBound =  1; //<--- 左侧模板，右边界最大值  
    minofRHSMatrixBound = -7; //<--- 右侧模板，左边界最小值  
    maxofRHSMatrixBound =  7; //<--- 右侧模板，右边界最大值  

    nBoundary_Neumann = 6;
    AllocCompactDiffData(ist, ied);

    Range I(ist, ied);
    Range IL(minofLHSMatrixBound, maxofLHSMatrixBound);
    Range IR(minofRHSMatrixBound, maxofRHSMatrixBound);
    (*lHSCoefMatrix)(I, IL) = 0.0;
    (*rHSCoefMatrix)(I, IR) = 0.0;

    double alpha_C6 = 1.0 / 3.0;
    double alpha = alpha_C6;               
    double a_C6 = (2.0 * alpha + 4.0) / 3.0;
    double b_C6 = (4.0 * alpha - 1.0) / 3.0;

    //inner points
    (*lowerBoundofLHSCoefMatrix)(I) = -1;
    (*upperBoundofLHSCoefMatrix)(I) = 1;

    (*lHSCoefMatrix)(I,-1) = alpha_C6;
    (*lHSCoefMatrix)(I, 1) = alpha_C6;
    (*lHSCoefMatrix)(I, 0) = 1.0;

    (*lowerBoundofRHSCoefMatrix)(I) = -2;
    (*upperBoundofRHSCoefMatrix)(I) = 2;

    (*rHSCoefMatrix)(I, 2) = b_C6 / 4.0;
    (*rHSCoefMatrix)(I, 1) = a_C6 / 2.0;
    (*rHSCoefMatrix)(I,-2) = -1.0 * (*rHSCoefMatrix)(I, 2);
    (*rHSCoefMatrix)(I,-1) = -1.0 * (*rHSCoefMatrix)(I, 1);
    (*rHSCoefMatrix)(I, 0) = 0.0;

    //left boundary: i = ist
    int i = ist;
    double alp0 = 1809.257;
    double bet0 = -65.1944;

    (*lHSCoefMatrix)(i, IL) = 0.0;
    (*rHSCoefMatrix)(i, IR) = 0.0;

    (*lowerBoundofLHSCoefMatrix)(i) = 0;
    (*upperBoundofLHSCoefMatrix)(i) = 0;
    (*lowerBoundofRHSCoefMatrix)(i) = 0;
    (*upperBoundofRHSCoefMatrix)(i) = 7;

    (*lHSCoefMatrix)(i, 0) = 1.0;

    (*rHSCoefMatrix)(i, 0) = -(alp0 - 28.0 * bet0 +13068.0) / 5040.0;
    (*rHSCoefMatrix)(i, 1) =  (alp0 - 27.0 * bet0 + 5040.0) / 720.0 ;
    (*rHSCoefMatrix)(i, 2) = -(alp0 - 26.0 * bet0 + 2520.0) / 240.0 ;
    (*rHSCoefMatrix)(i, 3) =  (alp0 - 25.0 * bet0 + 1680.0) / 144.0 ;
    (*rHSCoefMatrix)(i, 4) = -(alp0 - 24.0 * bet0 + 1260.0) / 144.0 ;
    (*rHSCoefMatrix)(i, 5) =  (alp0 - 23.0 * bet0 + 1008.0) / 240.0 ;
    (*rHSCoefMatrix)(i, 6) = -(alp0 - 22.0 * bet0 +  840.0) / 720.0 ;
    (*rHSCoefMatrix)(i, 7) =  (alp0 - 21.0 * bet0 +  720.0) / 5040.0;

    //left boundary: i = ist + 1
    i = ist + 1;
    double alp1 = -262.16;
    double bet1 = -26.6742;

    (*lHSCoefMatrix)(i, IL) = 0.0;
    (*rHSCoefMatrix)(i, IR) = 0.0;

    (*lowerBoundofLHSCoefMatrix)(i) = 0;
    (*upperBoundofLHSCoefMatrix)(i) = 0;
    (*lowerBoundofRHSCoefMatrix)(i) = -1;
    (*upperBoundofRHSCoefMatrix)(i) = 6;

    (*lHSCoefMatrix)(i, 0) = 1.0;

    (*rHSCoefMatrix)(i,-1) = -(alp1 - 21.0 * bet1 + 720.0) / 5040.0;
    (*rHSCoefMatrix)(i, 0) =  (alp1 - 20.0 * bet1 -1044.0) / 720.0 ;
    (*rHSCoefMatrix)(i, 1) = -(alp1 - 19.0 * bet1 - 720.0) / 240.0 ;
    (*rHSCoefMatrix)(i, 2) =  (alp1 - 18.0 * bet1 - 360.0) / 144.0 ;
    (*rHSCoefMatrix)(i, 3) = -(alp1 - 17.0 * bet1 - 240.0) / 144.0 ;
    (*rHSCoefMatrix)(i, 4) =  (alp1 - 16.0 * bet1 - 180.0) / 240.0 ;
    (*rHSCoefMatrix)(i, 5) = -(alp1 - 15.0 * bet1 - 144.0) / 720.0 ;
    (*rHSCoefMatrix)(i, 6) =  (alp1 - 14.0 * bet1 - 120.0) / 5040.0;

    //right boundary: i = ied - 1
    i = ied - 1;

    (*lHSCoefMatrix)(i, IL) = 0.0;
    (*rHSCoefMatrix)(i, IR) = 0.0;

    (*lowerBoundofLHSCoefMatrix)(i) = 0;
    (*upperBoundofLHSCoefMatrix)(i) = 0;
    (*lowerBoundofRHSCoefMatrix)(i) = -6;
    (*upperBoundofRHSCoefMatrix)(i) = 1;

    (*lHSCoefMatrix)(i, 0) = 1.0;

    (*rHSCoefMatrix)(i, 1) = -1.0 * (*rHSCoefMatrix)(ist+1,-1);
    (*rHSCoefMatrix)(i, 0) = -1.0 * (*rHSCoefMatrix)(ist+1, 0);
    (*rHSCoefMatrix)(i,-1) = -1.0 * (*rHSCoefMatrix)(ist+1, 1);
    (*rHSCoefMatrix)(i,-2) = -1.0 * (*rHSCoefMatrix)(ist+1, 2);
    (*rHSCoefMatrix)(i,-3) = -1.0 * (*rHSCoefMatrix)(ist+1, 3);
    (*rHSCoefMatrix)(i,-4) = -1.0 * (*rHSCoefMatrix)(ist+1, 4);
    (*rHSCoefMatrix)(i,-5) = -1.0 * (*rHSCoefMatrix)(ist+1, 5);
    (*rHSCoefMatrix)(i,-6) = -1.0 * (*rHSCoefMatrix)(ist+1, 6);

    //right boundary: i = ied
    i = ied;

    (*lHSCoefMatrix)(i, IL) = 0.0;
    (*rHSCoefMatrix)(i, IR) = 0.0;

    (*lowerBoundofLHSCoefMatrix)(i) = 0;
    (*upperBoundofLHSCoefMatrix)(i) = 0;
    (*lowerBoundofRHSCoefMatrix)(i) = -7;
    (*upperBoundofRHSCoefMatrix)(i) = 0;

    (*lHSCoefMatrix)(i, 0) = 1.0;

    (*rHSCoefMatrix)(i, 0) = -1.0 * (*rHSCoefMatrix)(ist, 0);
    (*rHSCoefMatrix)(i,-1) = -1.0 * (*rHSCoefMatrix)(ist, 1);
    (*rHSCoefMatrix)(i,-2) = -1.0 * (*rHSCoefMatrix)(ist, 2);
    (*rHSCoefMatrix)(i,-3) = -1.0 * (*rHSCoefMatrix)(ist, 3);
    (*rHSCoefMatrix)(i,-4) = -1.0 * (*rHSCoefMatrix)(ist, 4);
    (*rHSCoefMatrix)(i,-5) = -1.0 * (*rHSCoefMatrix)(ist, 5);
    (*rHSCoefMatrix)(i,-6) = -1.0 * (*rHSCoefMatrix)(ist, 6);
    (*rHSCoefMatrix)(i,-7) = -1.0 * (*rHSCoefMatrix)(ist, 7);

    Range N(0, nBoundary_Neumann + 1);
    Range NN(-(nBoundary_Neumann + 1), 0);
    (*boundary_Neumann_0)(N) = (*rHSCoefMatrix)(ist, N);
    (*boundary_Neumann_N)(NN) = (*rHSCoefMatrix)(ied, NN);
}

void CompactDifferenceFirstDerivative::SetCompactDiffData_optC4(const int ist, const int ied)
{
}

void CompactDifferenceFirstDerivative::SetCompactDiffData_optC4_InnerOnly(const int ist, const int ied, ExplicitDifferenceBoundary *DiffBoundaryData)
{
    minofLHSMatrixBound = -1; //<--- 左侧模板，左边界最小值  
    maxofLHSMatrixBound =  1; //<--- 左侧模板，右边界最大值  
    minofRHSMatrixBound = -4; //<--- 右侧模板，左边界最小值  
    maxofRHSMatrixBound =  4; //<--- 右侧模板，右边界最大值  

    nBoundary_Neumann = 3;
    AllocCompactDiffData(ist, ied);

    Range I (ist, ied);
    Range IL(minofLHSMatrixBound, maxofLHSMatrixBound);
    Range IR(minofRHSMatrixBound, maxofRHSMatrixBound);

    (*lHSCoefMatrix)(I, IL) = 0.0;
    (*rHSCoefMatrix)(I, IR) = 0.0;

    double alpha_optC4 = 0.3717019586;              
    double a_optC4 = (2.0 * alpha_optC4 + 4.0) / 3.0;
    double b_optC4 = (4.0 * alpha_optC4 - 1.0) / 3.0;

    // !inner points
    (*lowerBoundofLHSCoefMatrix)(I) = -1;
    (*upperBoundofLHSCoefMatrix)(I) = 1;
    (*lHSCoefMatrix)(I, -1) = alpha_optC4;
    (*lHSCoefMatrix)(I, 1) = alpha_optC4;
    (*lHSCoefMatrix)(I, 0) = 1.0;

    (*lowerBoundofRHSCoefMatrix)(I) = -2;
    (*upperBoundofRHSCoefMatrix)(I) = 2;
    (*rHSCoefMatrix)(I, 2) = b_optC4 / 4.0;
    (*rHSCoefMatrix)(I, 1) = a_optC4 / 2.0;
    (*rHSCoefMatrix)(I, -2) = -1.0 * (*rHSCoefMatrix)(I, 2);
    (*rHSCoefMatrix)(I, -1) = -1.0 * (*rHSCoefMatrix)(I, 1);
    (*rHSCoefMatrix)(I, 0) = 0.0;

    // !left boundary, i = ist + 1
    int i = ist + 1;
    (*lHSCoefMatrix)(i, IL) = 0.0;
    (*rHSCoefMatrix)(i, IR) = 0.0;

    (*lowerBoundofLHSCoefMatrix)(i) = 0;
    (*upperBoundofLHSCoefMatrix)(i) = 1;
    (*lowerBoundofRHSCoefMatrix)(i) = -1;
    (*upperBoundofRHSCoefMatrix)(i) = 2;

    double alpha_bc2 = 0.4;

    double bm1_C4BC2 = 1.0 / 6.0 * alpha_bc2 - 1.0 / 3.0; 
    double b00_C4BC2 =      -1.0 * alpha_bc2 - 1.0 / 2.0;
    double bp1_C4BC2 = 1.0 / 2.0 * alpha_bc2 + 1.0;
    double bp2_C4BC2 = 1.0 / 3.0 * alpha_bc2 - 1.0 / 6.0; 

    (*lHSCoefMatrix)(i, 0) = 1.0;
    (*lHSCoefMatrix)(i, 1) = alpha_bc2;

    (*rHSCoefMatrix)(i, -1) = bm1_C4BC2;
    (*rHSCoefMatrix)(i, 0) = b00_C4BC2;
    (*rHSCoefMatrix)(i, 1) = bp1_C4BC2;
    (*rHSCoefMatrix)(i, 2) = bp2_C4BC2;

    // !right boundary, i = ied - 1

    i = ied - 1;

    (*lHSCoefMatrix)(i, IL) = 0.0;
    (*rHSCoefMatrix)(i, IR) = 0.0;

    (*lowerBoundofLHSCoefMatrix)(i) = -1;
    (*upperBoundofLHSCoefMatrix)(i) = 0;
    (*lowerBoundofRHSCoefMatrix)(i) = -2;
    (*upperBoundofRHSCoefMatrix)(i) = 1;

    alpha_bc2 = 0.4;

    (*lHSCoefMatrix)(i, 0) = 1.0;
    (*lHSCoefMatrix)(i, -1) = alpha_bc2;

    (*rHSCoefMatrix)(i, 1) = -1.0 * bm1_C4BC2;
    (*rHSCoefMatrix)(i, 0) = -1.0 * b00_C4BC2;
    (*rHSCoefMatrix)(i, -1) = -1.0 * bp1_C4BC2;
    (*rHSCoefMatrix)(i, -2) = -1.0 * bp2_C4BC2;

    // !left boundary, i = ist
    i = ist;

    (*lHSCoefMatrix)(i, IL) = 0.0;
    (*rHSCoefMatrix)(i, IR) = 0.0;

    (*lowerBoundofLHSCoefMatrix)(i) = 0;
    (*upperBoundofLHSCoefMatrix)(i) = 0;
    (*lowerBoundofRHSCoefMatrix)(i) = 0;
    (*upperBoundofRHSCoefMatrix)(i) = 4;

    (*lHSCoefMatrix)(i, 0) = 1.0;

    RDouble1D *boundary_Neumann_0 = DiffBoundaryData -> boundary_Neumann_0;
    (*rHSCoefMatrix)(i, 0) = (*boundary_Neumann_0)(0);
    (*rHSCoefMatrix)(i, 1) = (*boundary_Neumann_0)(1);
    (*rHSCoefMatrix)(i, 2) = (*boundary_Neumann_0)(2);
    (*rHSCoefMatrix)(i, 3) = (*boundary_Neumann_0)(3);
    (*rHSCoefMatrix)(i, 4) = (*boundary_Neumann_0)(4);

    // !right boundary, i = ied
    i = ied;

    (*lHSCoefMatrix)(i, IL) = 0.0;
    (*rHSCoefMatrix)(i, IR) = 0.0;

    (*lowerBoundofLHSCoefMatrix)(i) = 0;
    (*upperBoundofLHSCoefMatrix)(i) = 0;
    (*lowerBoundofRHSCoefMatrix)(i) = -4;
    (*upperBoundofRHSCoefMatrix)(i) = 0;

    (*lHSCoefMatrix)(i, 0) = 1.0;

    (*rHSCoefMatrix)(i, 0) = -1.0 * (*boundary_Neumann_0)(0);
    (*rHSCoefMatrix)(i, -1) = -1.0 * (*boundary_Neumann_0)(1);
    (*rHSCoefMatrix)(i, -2) = -1.0 * (*boundary_Neumann_0)(2);
    (*rHSCoefMatrix)(i, -3) = -1.0 * (*boundary_Neumann_0)(3);
    (*rHSCoefMatrix)(i, -4) = -1.0 * (*boundary_Neumann_0)(4);
}

void CompactDifferenceFirstDerivative::GetGradientVector( PHComplex *varAddress, PHComplex *dUdX, PHComplex *dUdY, PHComplex *dUdZ, Complex1D *varTemp, Complex1D *dvarTemp, 
                                                         Complex1D *workSpace1, RDouble1D *workSpace2, SpecDiffHybGrid *GridData )
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex  ;
    int jst = GridData -> jStartFourierIndex;
    int jed = GridData -> jEndFourierIndex  ;
    int kst = GridData -> kStartFourierIndex;
    int ked = GridData -> kEndFourierIndex;

    Complex1D *kx = GridData -> complexKX;
    Complex1D *ky = GridData -> complexKY;
    Range I(ist, ied);
    Range J(jst, jed);
    Range K(kst, ked);
    Complex3D var( varAddress, I, J, K, neverDeleteData, fortranArray );
    Complex3D dVardX( dUdX, I, J, K, neverDeleteData, fortranArray );
    Complex3D dVardY( dUdY, I, J, K, neverDeleteData, fortranArray );
    Complex3D dVardZ( dUdZ, I, J, K, neverDeleteData, fortranArray );

    for( int ix = kst; ix <= ked; ++ix)
    {
        for( int iy = jst; iy <= jed; ++iy)
        {
            for( int iz = ist; iz <= ied; ++iz)
            {
                dVardX(iz, iy, ix) = (*kx)(ix) * var(iz, iy, ix);
                dVardY(iz, iy, ix) = (*ky)(iy) * var(iz, iy, ix);
            }
        }
    }

    // !Get gradient in direction z
    GetDvarDz( GridData, &var, &dVardZ, varTemp, dvarTemp, workSpace1, workSpace2 );
}

void CompactDifferenceFirstDerivative::GetDivergence( PHComplex *vecX, PHComplex *vecY, PHComplex *vecZ, Complex3D *divVar, Complex1D *varTemp, Complex1D *dVarTemp, Complex3D *dVardZ, 
                                                     Complex1D *workSpace1, RDouble1D *workSpace2, SpecDiffHybGrid *GridData )
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex  ;
    int jst = GridData -> jStartFourierIndex;
    int jed = GridData -> jEndFourierIndex  ;
    int kst = GridData -> kStartFourierIndex;
    int ked = GridData -> kEndFourierIndex;
    Range I(ist, ied);
    Range J(jst, jed);
    Range K(kst, ked);
    Complex3D vecVarX( vecX, I, J, K, neverDeleteData, fortranArray );
    Complex3D vecVarY( vecY, I, J, K, neverDeleteData, fortranArray );
    Complex3D vecVarZ( vecZ, I, J, K, neverDeleteData, fortranArray );

    Complex1D *kx = GridData -> complexKX;
    Complex1D *ky = GridData -> complexKY;

    (*divVar)(I, J, K) = PHComplex(0.0, 0.0);
    for( int ix = kst; ix <= ked; ++ix)
    {
        for( int iy = jst; iy <= jed; ++iy)
        {
            for( int iz = ist; iz <= ied; ++iz)
            {
                (*divVar)(iz, iy, ix) = (*divVar)(iz, iy, ix) + (*kx)(ix) * vecVarX(iz, iy, ix) + (*ky)(iy) * vecVarY(iz, iy, ix);
            }
        }
    }

    GetDvarDz( GridData, &vecVarZ, dVardZ, varTemp, dVarTemp, workSpace1, workSpace2 );
    (*divVar)(I, J, K) = (*divVar)(I, J, K) + (*dVardZ)(I, J, K);
}

void CompactDifferenceFirstDerivative::GetDvarDz( SpecDiffHybGrid *GridData, Complex3D *var, Complex3D *dVardZ, Complex1D *varTemp, Complex1D *dvarTemp, Complex1D *workSpace1,
                                                 RDouble1D *workSpace2 )
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex  ;
    int jst = GridData -> jStartFourierIndex;
    int jed = GridData -> jEndFourierIndex  ;
    int kst = GridData -> kStartFourierIndex;
    int ked = GridData -> kEndFourierIndex;
    RDouble1D *detaDz = GridData -> realDetaDz;
    Range I(ist, ied);

    for( int ix = kst; ix <= ked; ++ix )
    {
        for( int iy = jst; iy <= jed; ++iy )
        {
            (*varTemp)(I) = (*var)(I, iy, ix);
            GetFirstDerivative(1, ist, ied, varTemp, dvarTemp, lowerBoundofRHSCoefMatrix, upperBoundofRHSCoefMatrix, lHSCoefMatrix, rHSCoefMatrix, workSpace1, workSpace2);
            (*dVardZ)(I, iy, ix) = (*dvarTemp)(I) * (*detaDz)(I);
        }
    }
}

void CompactDifferenceFirstDerivative::GetDvarDz( const int ist, const int ied, RDouble *deta, RDouble *v, RDouble *dV, RDouble *cR, RDouble *rC )
{		
    RDouble1D detaDz( deta, Range(ist, ied), neverDeleteData, fortranArray);
    RDouble1D var( v, Range(ist, ied), neverDeleteData, fortranArray );
    RDouble1D dVardZ( dV, Range(ist, ied), neverDeleteData, fortranArray );
    RDouble1D workSpace1( cR, Range(ist, ied), neverDeleteData, fortranArray );
    RDouble1D workSpace2( rC, Range(ist, ied), neverDeleteData, fortranArray );
    GetFirstDerivative( 1, ist, ied, &var, &dVardZ, lowerBoundofRHSCoefMatrix, upperBoundofRHSCoefMatrix, lHSCoefMatrix, rHSCoefMatrix, &workSpace1, &workSpace2);
    for ( int i = ist; i <= ied; ++i )
    {
        dVardZ(i) = dVardZ(i) * detaDz(i);
    }
}

void CompactDifferenceFirstDerivative::GetDvarDz( const int ist, const int ied, RDouble *deta, PHComplex *v, PHComplex *dV, PHComplex *cR, RDouble *rC )
{		
    RDouble1D detaDz( deta, Range(ist, ied), neverDeleteData, fortranArray);
    Complex1D var( v, Range(ist, ied), neverDeleteData, fortranArray );
    Complex1D dVardZ( dV, Range(ist, ied), neverDeleteData, fortranArray );
    Complex1D workSpace1( cR, Range(ist, ied), neverDeleteData, fortranArray );
    RDouble1D workSpace2( rC, Range(ist, ied), neverDeleteData, fortranArray );

    GetFirstDerivative( 1, ist, ied, &var, &dVardZ, lowerBoundofRHSCoefMatrix, upperBoundofRHSCoefMatrix, lHSCoefMatrix, rHSCoefMatrix, &workSpace1, &workSpace2);

    Range I(ist, ied);
    dVardZ(I) = dVardZ(I) * detaDz(I);
}

void CompactDifferenceFirstDerivative::GetDetaDz( SpecDiffHybGrid *GridData )
{		
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex  ;
    RDouble1D *Z      = GridData -> realZ;
    RDouble1D *detaDz = GridData -> realDetaDz;

    Range I (ist, ied);
    RDouble1D *realWorkSpace1 = new RDouble1D(I, fortranArray);
    RDouble1D *realWorkSpace2 = new RDouble1D(I, fortranArray);

    GetFirstDerivative(1, ist, ied, Z, detaDz, lowerBoundofRHSCoefMatrix, upperBoundofRHSCoefMatrix, lHSCoefMatrix, rHSCoefMatrix, realWorkSpace1, realWorkSpace2);
    for (int i = ist; i <= ied; ++i)
    {
        (*detaDz)(i) = 1.0 / (*detaDz)(i);
    }

    delete realWorkSpace1; realWorkSpace1 = NULL;
    delete realWorkSpace2; realWorkSpace2 = NULL;
}

void CompactDifferenceFirstDerivative::GetFirstDerivative( const int neqn, const int ist, const int ied, RDouble1D *var, RDouble1D *dvar, Int1D *lowerBound_RHS, 
                                                          Int1D *upperBound_RHS, RDouble2D *coeff_LHS, RDouble2D *coeff_RHS, RDouble1D *WorkSpace1,
                                                          RDouble1D *WorkSpace2 )
{
    int n = ied - ist + 1;
    Range I(ist, ied);
    (*dvar)(I) = 0.0;

    for (int i = ist; i <= ied; ++i)
    {
        for ( int idp = (*lowerBound_RHS)(i); idp <= (*upperBound_RHS)(i); ++idp)
        {
            int iip = i + idp;
            (*dvar)(i) = (*dvar)(i) + (*coeff_RHS)(i, idp) * (*var)(iip);
        }
    }

    for (int ieq = 1; ieq <= neqn; ++ieq)
    {
        Tridiag_scalar(n, &(*coeff_LHS)(ist, -1), &(*coeff_LHS)(ist, 0), &(*coeff_LHS)(ist, 1), &(*dvar)(ist), &(*WorkSpace1)(ist), &(*WorkSpace2)(ist));
    }
}

void CompactDifferenceFirstDerivative::GetFirstDerivative( const int neqn, const int ist, const int ied, Complex1D *var, Complex1D *dvar, Int1D *lowerBound_RHS, 
                                                          Int1D *upperBound_RHS, RDouble2D *coeff_LHS, RDouble2D *coeff_RHS, Complex1D *WorkSpace1,
                                                          RDouble1D *WorkSpace2 )
{
    int n = ied - ist + 1;
    Range I(ist, ied);
    (*dvar)(I) = PHComplex(0.0, 0.0);

    for (int i = ist; i <= ied; ++i)
    {
        for ( int idp = (*lowerBound_RHS)(i); idp <= (*upperBound_RHS)(i); ++idp)
        {
            int iip = i + idp;
            (*dvar)(i) = (*dvar)(i) + (*coeff_RHS)(i, idp) * (*var)(iip);
        }
    }

    for (int ieq = 1; ieq <= neqn; ++ieq)
    {
        Tridiag_scalar(n, &(*coeff_LHS)(ist, -1), &(*coeff_LHS)(ist, 0), &(*coeff_LHS)(ist, 1), &(*dvar)(ist), &(*WorkSpace1)(ist), &(*WorkSpace2)(ist) );
    }
}

void CompactDifferenceFirstDerivative::Tridiag_scalar(const int n, RDouble *a_in, RDouble *b_in, RDouble *c_in, RDouble *x_in, RDouble *ry_in, RDouble *rc_in)
{
    RDouble1D a(a_in, Range(1, n), neverDeleteData, fortranArray);
    RDouble1D b(b_in, Range(1, n), neverDeleteData, fortranArray);
    RDouble1D c(c_in, Range(1, n), neverDeleteData, fortranArray);
    RDouble1D x(x_in, Range(1, n), neverDeleteData, fortranArray);
    RDouble1D ry(ry_in, Range(1, n), neverDeleteData, fortranArray);
    RDouble1D rc(rc_in, Range(1, n), neverDeleteData, fortranArray);

    double r = 1.0/b(1);
    rc(1) = r * c(1);
    ry(1) = r * x(1);
    for(int i = 2; i <= n ; ++i)
    {
        r = 1.0/( b(i) - rc(i-1) * a(i) );
        ry(i) = r * ( x(i) - a(i) * ry(i-1) );
        rc(i) = r * c(i);
    }

    x(n) = ry(n);
    for(int i = n-1 ; i >= 1 ; --i)
    {
        x(i) = ry(i) - rc(i) * x(i+1);
    }
}

void CompactDifferenceFirstDerivative::Tridiag_scalar(const int n, RDouble *a_in, RDouble *b_in, RDouble *c_in, PHComplex *cx_in, PHComplex *cy_in, RDouble *rc_in)
{
    RDouble1D a(a_in, Range(1, n), neverDeleteData, fortranArray);
    RDouble1D b(b_in, Range(1, n), neverDeleteData, fortranArray);
    RDouble1D c(c_in, Range(1, n), neverDeleteData, fortranArray);
    Complex1D cx(cx_in, Range(1, n), neverDeleteData, fortranArray);
    Complex1D cy(cy_in, Range(1, n), neverDeleteData, fortranArray);
    RDouble1D rc(rc_in, Range(1, n), neverDeleteData, fortranArray);

    double r = 1.0/b(1);
    rc(1) = r * c(1);
    cy(1) = r * cx(1);
    for(int i = 2; i <= n ; ++i)
    {
        r = 1.0/( b(i) - rc(i-1) * a(i) );
        cy(i) = r * ( cx(i) - a(i) * cy(i-1) );
        rc(i) = r * c(i);
    }

    cx(n) = cy(n);
    for(int i = n-1 ; i >= 1 ; --i)
    {
        cx(i) = cy(i) - rc(i) * cx(i+1);
    }
}

void CompactDifferenceFirstDerivative::AllocCompactDiffData(const int ist, const int ied)
{
    Range N (0, nBoundary_Neumann + 1 );
    Range NN(-(nBoundary_Neumann + 1), 0);
    boundary_Neumann_0 = new RDouble1D(N, fortranArray);
    boundary_Neumann_N = new RDouble1D(NN, fortranArray);

    Range I (ist, ied);
    Range IL(minofLHSMatrixBound, maxofLHSMatrixBound);
    Range IR(minofRHSMatrixBound, maxofRHSMatrixBound);
    lowerBoundofLHSCoefMatrix = new Int1D(I, fortranArray);
    upperBoundofLHSCoefMatrix = new Int1D(I, fortranArray);
    lowerBoundofRHSCoefMatrix = new Int1D(I, fortranArray);
    upperBoundofRHSCoefMatrix = new Int1D(I, fortranArray);
    lHSCoefMatrix = new RDouble2D(I, IL, fortranArray);
    rHSCoefMatrix = new RDouble2D(I, IR, fortranArray);
}

CompactDifferenceFirstDerivative::~CompactDifferenceFirstDerivative()
{
    delete boundary_Neumann_0; boundary_Neumann_0 = NULL;
    delete boundary_Neumann_N; boundary_Neumann_N = NULL;
    delete lowerBoundofLHSCoefMatrix; lowerBoundofLHSCoefMatrix = NULL;
    delete upperBoundofLHSCoefMatrix; upperBoundofLHSCoefMatrix = NULL;
    delete lowerBoundofRHSCoefMatrix; lowerBoundofRHSCoefMatrix = NULL;
    delete upperBoundofRHSCoefMatrix; upperBoundofRHSCoefMatrix = NULL;
    delete lHSCoefMatrix; lHSCoefMatrix = NULL;
    delete rHSCoefMatrix; rHSCoefMatrix = NULL;
}
}

