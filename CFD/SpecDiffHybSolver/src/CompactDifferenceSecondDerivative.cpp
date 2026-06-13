#include "DataStruct_Array.h"
#include "GlobalDataBase.h"
#include <math.h>
#include "CompactDifferenceSecondDerivative.h"
#include "PoissonSolver.h"
#include "TK_Exit.h"
#include "PHHeader.h"
#include "PHMpi.h"
#include <iostream>
#include <sstream>
#include <cstdlib>

using namespace std;

namespace PHSPACE
{
CompactDifferenceSecondDerivative::CompactDifferenceSecondDerivative()
{
}

void CompactDifferenceSecondDerivative::InitCompactDifferenceSecondDerivativeData(SpecDiffHybGrid *GridData)
{
    AllocCompactDifferenceSecondDerivativeData(GridData);		

    SetCompactDifferenceSecondDerivativeData(GridData);   
}

void CompactDifferenceSecondDerivative::AllocCompactDifferenceSecondDerivativeData(SpecDiffHybGrid *GridData)
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;

    int iplst = ist + 1;
    int ipled = ied - 1;

    Range IP (iplst, ipled);
    lowerBoundofLHSCoefMatrix = new Int1D(IP, fortranArray);
    upperBoundofLHSCoefMatrix = new Int1D(IP, fortranArray);
    lowerBoundofRHSCoefMatrix = new Int1D(IP, fortranArray);
    upperBoundofRHSCoefMatrix = new Int1D(IP, fortranArray);

    Range J (-2, 2);
    lHSCoefMatrix = new RDouble2D(IP, J, fortranArray);
    rHSCoefMatrix = new RDouble2D(IP, J, fortranArray);
    coefd2L = new RDouble2D(IP, J, fortranArray);
    coefd2R = new RDouble2D(IP, J, fortranArray);
    coefLUDecp = new RDouble2D(IP, J, fortranArray);

    coefOfLHSBoundary0 = new RDouble1D(J, fortranArray);
    coefOfLHSBoundaryN = new RDouble1D(J, fortranArray);
    coefOfRHSBoundary0 = new RDouble1D(J, fortranArray);
    coefOfRHSBoundaryN = new RDouble1D(J, fortranArray);
}

void CompactDifferenceSecondDerivative::SetCompactDifferenceSecondDerivativeData(SpecDiffHybGrid *GridData)
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;

    int nn = ied - ist;
    RDouble1D *Z = GridData -> realZ;

    Range I(ist, ied);
    Range J(-2, 2);
    RDouble2D *coef_d2L = new RDouble2D(I, J, fortranArray);
    RDouble2D *coef_d2R = new RDouble2D(I, J, fortranArray);

    GetCoefMatrix( nn, &(*coef_d2L)(ist, -2), &(*coef_d2L)(ist, -1), &(*coef_d2L)(ist, 0), &(*coef_d2L)(ist, 1), &(*coef_d2L)(ist, 2),
        &(*coef_d2R)(ist, -2), &(*coef_d2R)(ist, -1), &(*coef_d2R)(ist, 0), &(*coef_d2R)(ist, 1), &(*coef_d2R)(ist, 2), &(*Z)(ist) );

    Range IP(ist + 1, ied - 1);
    (*coefLUDecp)(IP, J) = (*coef_d2L)(IP, J);
    (*coefd2L)(IP, J) = (*coef_d2L)(IP, J);
    (*coefd2R)(IP, J) = (*coef_d2R)(IP, J);

    PoissonSolver::LUDecompose(nn - 1, &(*coefLUDecp)(ist+1, -2));

    (*coefOfLHSBoundary0)(J) = (*coef_d2L)(ist, J);
    (*coefOfLHSBoundaryN)(J) = (*coef_d2L)(ied, J);
    (*coefOfRHSBoundary0)(J) = (*coef_d2R)(ist, J);
    (*coefOfRHSBoundaryN)(J) = (*coef_d2R)(ied, J);

    (*lHSCoefMatrix)(IP, J) = (*coefd2L)(IP, J);
    (*rHSCoefMatrix)(IP, J) = (*coefd2R)(IP, J);

    // !set lowerBoundofLHSCoefMatrix, upperBoundofLHSCoefMatrix, lowerBoundofRHSCoefMatrix, upperBoundofRHSCoefMatrix
    // !inner points
    (*lowerBoundofLHSCoefMatrix)(IP) = -2;
    (*upperBoundofLHSCoefMatrix)(IP) =  2;
    (*lowerBoundofRHSCoefMatrix)(IP) = -2;
    (*upperBoundofRHSCoefMatrix)(IP) =  2;

    // !left boundary: i = ist + 1
    int i = ist + 1;
    (*lowerBoundofLHSCoefMatrix)(i) =  0;
    (*lowerBoundofRHSCoefMatrix)(i) = -1;

    // !left boundary: i = ist + 2
    i = ist + 2;
    (*lowerBoundofLHSCoefMatrix)(i) = -1;

    // !right boundary: i = ied - 2
    i = ied - 2;
    (*upperBoundofLHSCoefMatrix)(i) =  1;

    // !right boundary: i = ied - 1
    i = ied - 1;
    (*upperBoundofLHSCoefMatrix)(i) =  0;
    (*upperBoundofRHSCoefMatrix)(i) =  1;

    delete coef_d2L; coef_d2L = NULL;
    delete coef_d2R; coef_d2R = NULL;
}

void CompactDifferenceSecondDerivative::GetCoefMatrix( const int mz, RDouble *ben , RDouble *aln , RDouble *aaa , RDouble *alp , RDouble *bep ,
                                                      RDouble *an2 , RDouble *an1 , RDouble *a00 , RDouble *ap1 , RDouble *ap2 , RDouble *x )
{
    RDouble1D *rb, *rb0, *rb1, *rb2;
    RDouble2D *rA, *rA0, *rA1, *rA2;

    Range I9 (1, 9);
    Range I8 (1, 8);
    Range I6 (1, 6);
    Range I5 (1, 5);

    rb  = new RDouble1D(I9, fortranArray);
    rA  = new RDouble2D(I9, I9, fortranArray);

    rb0 = new RDouble1D(I5, fortranArray);
    rA0 = new RDouble2D(I5, I5, fortranArray);

    rb1 = new RDouble1D(I6, fortranArray);
    rA1 = new RDouble2D(I6, I6, fortranArray);

    rb2 = new RDouble1D(I8, fortranArray);
    rA2 = new RDouble2D(I8, I8, fortranArray);

    int i = 0;
    for(i = 0; i <= mz; ++i)
    {
        aaa[i] = 1.0;
        ben[i] = 0.0;
        aln[i] = 0.0;
        alp[i] = 0.0;
        bep[i] = 0.0;
        an2[i] = 0.0;
        an1[i] = 0.0;
        a00[i] = 0.0;
        ap1[i] = 0.0;
        ap2[i] = 0.0;
    }

    // left boundary: i = 0
    i = 0;

    (*rA0)(1, 1) = 0.0;
    (*rA0)(1, 2) = 0.0;
    (*rA0)(1, 3) = 1.0;
    (*rA0)(1, 4) = 1.0;
    (*rA0)(1, 5) = 1.0;
    (*rb0)(1) = 0.0;

    (*rA0)(2, 1) = 0.0;
    (*rA0)(2, 2) = 0.0;
    (*rA0)(2, 3) = 0.0;
    (*rA0)(2, 4) = x[i+1] - x[i];
    (*rA0)(2, 5) = x[i+2] - x[i];
    (*rb0)(2) = 0.0;

    (*rA0)(3, 1) = 2.0;
    (*rA0)(3, 2) = 2.0;
    (*rA0)(3, 3) = 0.0;
    (*rA0)(3, 4) = -pow(x[i+1] - x[i], 2.0);
    (*rA0)(3, 5) = -pow(x[i+2] - x[i], 2.0);
    (*rb0)(3) = -2.0;

    for ( int n = 4 ; n <= 5 ; ++n )
    {
        int k = n - 2;
        (*rA0)(n, 1) = k * (k + 1) * pow(x[i+1]-x[i], k-1);
        (*rA0)(n, 2) = k * (k + 1) * pow(x[i+2]-x[i], k-1);
        (*rA0)(n, 3) = 0.0 ;
        (*rA0)(n, 4) = - pow(x[i+1]-x[i], k+1);
        (*rA0)(n, 5) = - pow(x[i+2]-x[i], k+1);
        (*rb0)(n) = 0.0;
    }

    Gaussj(&(*rA0)(1, 1), 5, 5, &(*rb0)(1), 1, 1);

    ben[i] = 0.0; 
    aln[i] = 0.0; 
    alp[i] = (*rb0)(1); 
    bep[i] = (*rb0)(2);	
    an2[i] = 0.0; 
    an1[i] = 0.0; 
    a00[i] = (*rb0)(3); 
    ap1[i] = (*rb0)(4); 
    ap2[i] = (*rb0)(5);

    // left boundary: i = 1
    i = 1;

    (*rA1)(1, 1) = 0.0;
    (*rA1)(1, 2) = 0.0;
    (*rA1)(1, 3) = 1.0;
    (*rA1)(1, 4) = 1.0;
    (*rA1)(1, 5) = 1.0;
    (*rA1)(1, 6) = 1.0;
    (*rb1)(1) = 0.0;

    (*rA1)(2, 1) = 0.0;
    (*rA1)(2, 2) = 0.0;
    (*rA1)(2, 3) = x[i-1] - x[i];
    (*rA1)(2, 4) = 0.0;
    (*rA1)(2, 5) = x[i+1] - x[i];
    (*rA1)(2, 6) = x[i+2] - x[i];
    (*rb1)(2 ) = 0.0;

    (*rA1)(3, 1) = 2.0;
    (*rA1)(3, 2) = 2.0;
    (*rA1)(3, 3) = - pow(x[i-1] - x[i], 2.0);
    (*rA1)(3, 4) = 0.0;
    (*rA1)(3, 5) = - pow(x[i+1] - x[i], 2.0);
    (*rA1)(3, 6) = - pow(x[i+2] - x[i], 2.0);
    (*rb1)(3) = - 2.0;

    for ( int n = 4 ; n <= 6 ; ++n )
    {
        int k = n - 2;
        (*rA1)(n, 1) = k * (k+1) * pow(x[i+1] - x[i], k-1);
        (*rA1)(n, 2) = k * (k+1) * pow(x[i+2] - x[i], k-1);
        (*rA1)(n, 3) = - pow(x[i-1] - x[i], k+1);
        (*rA1)(n, 4) = 0.0;
        (*rA1)(n, 5) = - pow(x[i+1] - x[i], k+1);
        (*rA1)(n, 6) = - pow(x[i+2] - x[i], k+1);
        (*rb1)(n) = 0.0  ;
    }

    Gaussj(&(*rA1)(1, 1), 6, 6, &(*rb1)(1), 1, 1);

    ben[i] = 0.0; 
    aln[i] = 0.0; 
    alp[i] = (*rb1)(1); 
    bep[i] = (*rb1)(2);	
    an2[i] = 0.0; 
    an1[i] = (*rb1)(3); 
    a00[i] = (*rb1)(4); 
    ap1[i] = (*rb1)(5); 
    ap2[i] = (*rb1)(6);

    // left boundary: i = 2
    i = 2;

    (*rA2)(1, 1) = 0.0;
    (*rA2)(1, 2) = 0.0;
    (*rA2)(1, 3) = 0.0;
    (*rA2)(1, 4) = 1.0;
    (*rA2)(1, 5) = 1.0;
    (*rA2)(1, 6) = 1.0;
    (*rA2)(1, 7) = 1.0;
    (*rA2)(1, 8) = 1.0;
    (*rb2)(1) = 0.0;

    (*rA2)(2, 1) = 0.0;
    (*rA2)(2, 2) = 0.0;
    (*rA2)(2, 3) = 0.0;
    (*rA2)(2, 4) = x[i-2] - x[i];
    (*rA2)(2, 5) = x[i-1] - x[i];
    (*rA2)(2, 6) = 0.0;
    (*rA2)(2, 7) = x[i+1] - x[i];
    (*rA2)(2, 8) = x[i+2] - x[i];
    (*rb2)(2) = 0.0;

    (*rA2)(3, 1) = 2.0;
    (*rA2)(3, 2) = 2.0;
    (*rA2)(3, 3) = 2.0;
    (*rA2)(3, 4) = - pow(x[i-2] - x[i], 2.0);
    (*rA2)(3, 5) = - pow(x[i-1] - x[i], 2.0);
    (*rA2)(3, 6) = 0.0;
    (*rA2)(3, 7) = - pow(x[i+1] - x[i], 2.0);
    (*rA2)(3, 8) = - pow(x[i+2] - x[i], 2.0);
    (*rb2)(3) = - 2.0;

    for (int n = 4; n <= 8; ++n)
    {
        int k = n - 2;
        (*rA2)(n, 1) = k * (k+1) * pow(x[i-1] - x[i], k-1);
        (*rA2)(n, 2) = k * (k+1) * pow(x[i+1] - x[i], k-1);
        (*rA2)(n, 3) = k * (k+1) * pow(x[i+2] - x[i], k-1);
        (*rA2)(n, 4) = - pow(x[i-2] - x[i], k+1);
        (*rA2)(n, 5) = - pow(x[i-1] - x[i], k+1);
        (*rA2)(n, 6) = 0.0;
        (*rA2)(n, 7) = - pow(x[i+1] - x[i], k+1);
        (*rA2)(n, 8) = - pow(x[i+2] - x[i], k+1);
        (*rb2)(n) = 0.0;
    }

    Gaussj(&(*rA2)(1, 1), 8, 8, &(*rb2)(1), 1, 1);

    ben[i] = 0.0; 
    aln[i] = (*rb2)(1); 
    alp[i] = (*rb2)(2); 
    bep[i] = (*rb2)(3);	
    an2[i] = (*rb2)(4); 
    an1[i] = (*rb2)(5); 
    a00[i] = (*rb2)(6); 
    ap1[i] = (*rb2)(7); 
    ap2[i] = (*rb2)(8);

    // right boundary: i = mz - 2
    i = mz - 2;

    (*rA2)(1, 1) = 0.0;
    (*rA2)(1, 2) = 0.0;
    (*rA2)(1, 3) = 0.0;
    (*rA2)(1, 4) = 1.0;
    (*rA2)(1, 5) = 1.0;
    (*rA2)(1, 6) = 1.0;
    (*rA2)(1, 7) = 1.0;
    (*rA2)(1, 8) = 1.0;
    (*rb2)(1) = 0.0;

    (*rA2)(2, 1) = 0.0;
    (*rA2)(2, 2) = 0.0;
    (*rA2)(2, 3) = 0.0;
    (*rA2)(2, 4) = x[i-2] - x[i];
    (*rA2)(2, 5) = x[i-1] - x[i];
    (*rA2)(2, 6) = 0.0;
    (*rA2)(2, 7) = x[i+1] - x[i];
    (*rA2)(2, 8) = x[i+2] - x[i];
    (*rb2)(2) = 0.0;

    (*rA2)(3, 1) = 2.0;
    (*rA2)(3, 2) = 2.0;
    (*rA2)(3, 3) = 2.0;
    (*rA2)(3, 4) = - pow(x[i-2] - x[i], 2.0);
    (*rA2)(3, 5) = - pow(x[i-1] - x[i], 2.0);
    (*rA2)(3, 6) = 0.0;
    (*rA2)(3, 7) = - pow(x[i+1] - x[i], 2.0);
    (*rA2)(3, 8) = - pow(x[i+2] - x[i], 2.0);
    (*rb2)(3) = - 2.0;

    for ( int n = 4 ; n <= 8 ; ++n )
    {
        int k = n - 2;
        (*rA2)(n, 1) = k * (k+1) * pow(x[i-2] - x[i], k-1);
        (*rA2)(n, 2) = k * (k+1) * pow(x[i-1] - x[i], k-1);
        (*rA2)(n, 3) = k * (k+1) * pow(x[i+1] - x[i], k-1);
        (*rA2)(n, 4) = - pow(x[i-2] - x[i], k+1);
        (*rA2)(n, 5) = - pow(x[i-1] - x[i], k+1);
        (*rA2)(n, 6) = 0.0;
        (*rA2)(n, 7) = - pow(x[i+1] - x[i], k+1);
        (*rA2)(n, 8) = - pow(x[i+2] - x[i], k+1);
        (*rb2)(n) = 0.0;
    }

    Gaussj(&(*rA2)(1, 1), 8, 8, &(*rb2)(1), 1, 1);

    ben[i] = (*rb2)(1); 
    aln[i] = (*rb2)(2); 
    alp[i] = (*rb2)(3); 
    bep[i] = 0.0;	
    an2[i] = (*rb2)(4); 
    an1[i] = (*rb2)(5); 
    a00[i] = (*rb2)(6); 
    ap1[i] = (*rb2)(7); 
    ap2[i] = (*rb2)(8);

    // right boundary: i = mz - 1
    i = mz - 1;

    (*rA1)(1, 1) = 0.0;
    (*rA1)(1, 2) = 0.0;
    (*rA1)(1, 3) = 1.0;
    (*rA1)(1, 4) = 1.0;
    (*rA1)(1, 5) = 1.0;
    (*rA1)(1, 6) = 1.0;
    (*rb1)(1) = 0.0;

    (*rA1)(2, 1) = 0.0;
    (*rA1)(2, 2) = 0.0;
    (*rA1)(2, 3) = x[i-2] - x[i];
    (*rA1)(2, 4) = x[i-1] - x[i];
    (*rA1)(2, 5) = 0.0;
    (*rA1)(2, 6) = x[i+1] - x[i];
    (*rb1)(2) = 0.0;

    (*rA1)(3, 1) = 2.0;
    (*rA1)(3, 2) = 2.0;
    (*rA1)(3, 3) = - pow(x[i-2] - x[i], 2.0);
    (*rA1)(3, 4) = - pow(x[i-1] - x[i], 2.0);
    (*rA1)(3, 5) = 0.0;
    (*rA1)(3, 6) = - pow(x[i+1] - x[i], 2.0);
    (*rb1)(3) = - 2.0;

    for (int n = 4; n <= 6; ++n)
    {
        int k = n - 2;
        (*rA1)(n, 1) = k * (k+1) * pow(x[i-2] - x[i], k-1);
        (*rA1)(n, 2) = k * (k+1) * pow(x[i-1] - x[i], k-1);
        (*rA1)(n, 3) =           - pow(x[i-2] - x[i], k+1);
        (*rA1)(n, 4) =           - pow(x[i-1] - x[i], k+1);
        (*rA1)(n, 5) = 0.0;
        (*rA1)(n, 6) = - pow(x[i+1] - x[i], k+1);
        (*rb1)(n) = 0.0;
    }

    Gaussj(&(*rA1)(1, 1), 6, 6, &(*rb1)(1), 1, 1);

    ben[i] = (*rb1)(1); 
    aln[i] = (*rb1)(2); 
    alp[i] = 0.0; 
    bep[i] = 0.0;	
    an2[i] = (*rb1)(3); 
    an1[i] = (*rb1)(4); 
    a00[i] = (*rb1)(5); 
    ap1[i] = (*rb1)(6); 
    ap2[i] = 0.0;

    //right boundary: i = mz;
    i = mz;

    (*rA0)(1, 1) = 0.0;
    (*rA0)(1, 2) = 0.0;
    (*rA0)(1, 3) = 1.0;
    (*rA0)(1, 4) = 1.0;
    (*rA0)(1, 5) = 1.0;
    (*rb0)(1) = 0.0;

    (*rA0)(2, 1) = 0.0;
    (*rA0)(2, 2) = 0.0;
    (*rA0)(2, 3) = x[i-2] - x[i];
    (*rA0)(2, 4) = x[i-1] - x[i];
    (*rA0)(2, 5) = 0.0;
    (*rb0)(2) = 0.0;

    (*rA0)(3, 1) = 2.0;
    (*rA0)(3, 2) = 2.0;
    (*rA0)(3, 3) = -pow(x[i-2] - x[i], 2.0);
    (*rA0)(3, 4) = -pow(x[i-1] - x[i], 2.0);
    (*rA0)(3, 5) = 0.0;
    (*rb0)(3) = - 2.0;

    for ( int n = 4 ; n <= 5 ; ++n )
    {
        int k = n - 2;
        (*rA0)(n, 1) = k * (k+1) * pow( x[i-2] - x[i] , k-1 );
        (*rA0)(n, 2) = k * (k+1) * pow( x[i-1] - x[i] , k-1 );
        (*rA0)(n, 3) = - pow( x[i-2] - x[i] , k+1 );
        (*rA0)(n, 4) = - pow( x[i-1] - x[i] , k+1 );
        (*rA0)(n, 5) = 0.0  ;
        (*rb0)(n) = 0.0  ;
    }

    Gaussj(&(*rA0)(1, 1), 5, 5, &(*rb0)(1), 1, 1);

    ben[i] = (*rb0)(1); 
    aln[i] = (*rb0)(2); 
    alp[i] = 0.0; 
    bep[i] = 0.0;	
    an2[i] = (*rb0)(3); 
    an1[i] = (*rb0)(4); 
    a00[i] = (*rb0)(5); 
    ap1[i] = 0.0; 
    ap2[i] = 0.0;

    // inner points
    for (i = 3 ; i <= mz-3 ; ++i)
    {
        (*rA)(1, 1) = 0.0;
        (*rA)(1, 2) = 0.0;
        (*rA)(1, 3) = 0.0;
        (*rA)(1, 4) = 0.0;
        (*rA)(1, 5) = 1.0;
        (*rA)(1, 6) = 1.0;
        (*rA)(1, 7) = 1.0;
        (*rA)(1, 8) = 1.0;
        (*rA)(1, 9) = 1.0;
        (*rb)(1) = 0.0;

        (*rA)(2, 1) = 0.0;
        (*rA)(2, 2) = 0.0;
        (*rA)(2, 3) = 0.0;
        (*rA)(2, 4) = 0.0;
        (*rA)(2, 5) = x[i-2] - x[i];
        (*rA)(2, 6) = x[i-1] - x[i];
        (*rA)(2, 7) = 0.0;
        (*rA)(2, 8) = x[i+1] - x[i];
        (*rA)(2, 9) = x[i+2] - x[i];
        (*rb)(2) = 0.0;

        (*rA)(3, 1) = 2.0;
        (*rA)(3, 2) = 2.0;
        (*rA)(3, 3) = 2.0;
        (*rA)(3, 4) = 2.0;
        (*rA)(3, 5) = - pow(x[i-2] - x[i], 2.0);
        (*rA)(3, 6) = - pow(x[i-1] - x[i], 2.0);
        (*rA)(3, 7) = 0.0;
        (*rA)(3, 8) = - pow(x[i+1] - x[i], 2.0);
        (*rA)(3, 9) = - pow(x[i+2] - x[i], 2.0);
        (*rb)(3) = - 2.0;

        for ( int n = 4 ; n <= 9 ; ++n )
        {
            int k = n - 2;
            (*rA)(n, 1) = k * (k+1) * pow(x[i-2] - x[i], k-1);
            (*rA)(n, 2) = k * (k+1) * pow(x[i-1] - x[i], k-1);
            (*rA)(n, 3) = k * (k+1) * pow(x[i+1] - x[i], k-1);
            (*rA)(n, 4) = k * (k+1) * pow(x[i+2] - x[i], k-1);
            (*rA)(n, 5) = - pow(x[i-2] - x[i], k+1);
            (*rA)(n, 6) = - pow(x[i-1] - x[i], k+1);
            (*rA)(n, 7) = 0.0;
            (*rA)(n, 8) = - pow(x[i+1] - x[i], k+1);
            (*rA)(n, 9) = - pow(x[i+2] - x[i], k+1);
            (*rb)(n) = 0.0;
        }

        Gaussj(&(*rA)(1, 1), 9, 9, &(*rb)(1), 1, 1);

        ben[i] = (*rb)(1); 
        aln[i] = (*rb)(2); 
        alp[i] = (*rb)(3); 
        bep[i] = (*rb)(4); 
        an2[i] = (*rb)(5); 
        an1[i] = (*rb)(6); 
        a00[i] = (*rb)(7); 
        ap1[i] = (*rb)(8); 
        ap2[i] = (*rb)(9);
    }

    delete rb; rb = NULL;
    delete rb0; rb0 = NULL;
    delete rb1; rb1 = NULL;
    delete rb2; rb2 = NULL;
    delete rA; rA = NULL;
    delete rA0; rA0 = NULL;
    delete rA1; rA1 = NULL;
    delete rA2; rA2 = NULL;
}

void CompactDifferenceSecondDerivative::GetLaplacian(SpecDiffHybGrid *GridData, Complex3D *cVar, Complex3D *cLaplacianVar, Complex1D *c1tmp, Complex1D *c2tmp, Complex3D *cd2Var)
{
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;
    Range IF(istf, iedf);
    Range JF(jstf, jedf);
    Range KF(kstf, kedf);

    RDouble1D *realKX = GridData -> realKX;
    RDouble1D *realKY = GridData -> realKY;

    for(int ix = kstf; ix <= kedf; ++ix)
    {
        double k2 = (*realKX)(ix) * (*realKX)(ix);
        (*cLaplacianVar)(IF, JF, ix) = -1.0 * k2 * (*cVar)(IF, JF, ix);
    }

    for(int iy = jstf; iy <= jedf; ++iy)
    {
        double k2 = (*realKY)(iy) * (*realKY)(iy);
        (*cLaplacianVar)(IF, iy, KF) = (*cLaplacianVar)(IF, iy, KF) - k2 * (*cVar)(IF, iy, KF);
    }

    GetD2varDz2(GridData, cVar, cd2Var, c1tmp, c2tmp);

    (*cLaplacianVar)(IF, JF, KF) = (*cLaplacianVar)(IF, JF, KF) + (*cd2Var)(IF, JF, KF);
}

void CompactDifferenceSecondDerivative::GetD2varDz2(SpecDiffHybGrid *GridData, Complex3D *cVar, Complex3D *cdVardx_2rd, Complex1D *cV, Complex1D *cd2Vardz)
{
    // !--- purpose: to get d^2(cVar)/dx^2, by second order derivative evaluated on non-uniform grid directly.
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;
    Range IF(istf, iedf);

    for(int ix = kstf; ix <= kedf; ++ix)
    {
        for(int iy = jstf; iy <= jedf; ++iy)
        {
            (*cV)(IF) = (*cVar)(IF, iy, ix);
            GetSecondDerivative(GridData, cV, cd2Vardz);
            (*cdVardx_2rd)(IF, iy, ix) = (*cd2Vardz)(IF);
        }
    }
}

void CompactDifferenceSecondDerivative::GetSecondDerivative(SpecDiffHybGrid *GridData, Complex1D *cVar, Complex1D *cd2Vardz)
{
    // !--- purpose: to get d(rVar)/dz by compact difference scheme
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    Range IF(istf, iedf);

    (*cd2Vardz)(IF) = PHComplex(0.0, 0.0);
    for(int iz = istf+1; iz <= iedf-1; ++iz)
    {
        for(int idp = (*lowerBoundofRHSCoefMatrix)(iz); idp <= (*upperBoundofRHSCoefMatrix)(iz); ++idp)
        {
            int iip = iz + idp;
            (*cd2Vardz)(iz) = (*cd2Vardz)(iz) + (*rHSCoefMatrix)(iz, idp) * (*cVar)(iip);
        }
    }

    int n = iedf - istf - 1;

    PoissonSolver::LUSolve(n, &(*coefLUDecp)(istf+1, -2), &(*cd2Vardz)(istf+1));

    (*cd2Vardz)(istf) = (*cVar)(istf  ) * (*coefOfRHSBoundary0)(0) 
        + (*cVar)(istf+1) * (*coefOfRHSBoundary0)(1)
        + (*cVar)(istf+2) * (*coefOfRHSBoundary0)(2)
        - (*cd2Vardz)(istf+1) * (*coefOfLHSBoundary0)(1) 
        - (*cd2Vardz)(istf+2) * (*coefOfLHSBoundary0)(2);

    (*cd2Vardz)(iedf) = (*cVar)(iedf  ) * (*coefOfRHSBoundaryN)(0) 
        + (*cVar)(iedf-1) * (*coefOfRHSBoundaryN)(-1)
        + (*cVar)(iedf-2) * (*coefOfRHSBoundaryN)(-2)
        - (*cd2Vardz)(iedf-1) * (*coefOfLHSBoundaryN)(-1) 
        - (*cd2Vardz)(iedf-2) * (*coefOfLHSBoundaryN)(-2);
}

void CompactDifferenceSecondDerivative::Gaussj(RDouble *rA, const int n, const int np, RDouble *rb, const int m, const int mp)
{
    RDouble2D a(rA, Range(1, np), Range(1, np), neverDeleteData, fortranArray);
    RDouble2D b(rb, Range(1, np), Range(1, mp), neverDeleteData, fortranArray);
    Int1D ipiv(Range(1, n), fortranArray);
    Int1D indxr(Range(1, n), fortranArray);
    Int1D indxc(Range(1, n), fortranArray);
    ipiv(Range(1, n)) = 0;
    indxr(Range(1, n)) = 0;
    indxc(Range(1, n)) = 0;
    double rdum = 0.0;

    for(int i = 1; i <= n; ++i)
    {
        double rBig = 0.0;
        int icol = 0;
        int irow = 0;
        for(int j = 1; j <= n; ++j)
        {
            if(ipiv(j) != 1)
            {
                for(int k = 1; k <= n; ++k)
                {
                    if(ipiv(k) == 0)
                    {
                        if(abs(a(j, k)) >= rBig)
                        {
                            rBig = abs(a(j, k));
                            irow = j;
                            icol = k;
                        }
                    }
                    else if(ipiv(k) > 1)
                    {
                        TK_Exit::ExceptionExit("Error: singular matrix in gaussj1\n");
                    } // !ipiv(k) == 0
                }// !for k
            }// !ipiv(j) != 1
        }// !for j

        ipiv(icol) = ipiv(icol) + 1;

        if(irow != icol)
        {
            for(int l = 1; l <= n; ++l)
            {
                rdum = a(irow, l);
                a(irow, l) = a(icol, l);
                a(icol, l) = rdum;
            }

            for(int l = 1; l <= m; ++l)
            {
                rdum = b(irow, l);
                b(irow, l) = b(icol, l);
                b(icol, l) = rdum;
            }
        }// !irow != icol

        indxr(i) = irow;
        indxc(i) = icol;
        if(a(icol,icol) == 0.0)
        {
            TK_Exit::ExceptionExit("Error: singular matrix in gaussj2\n");
        }

        double rpivinv = 1.0 / a(icol, icol);
        a(icol, icol) = 1.0;
        for(int l = 1; l <= n; ++l)
        {
            a(icol, l) = a(icol, l)*rpivinv;
        }
        for(int l = 1; l <= m; ++l)
        {
            b(icol, l) = b(icol, l)*rpivinv;
        }

        for(int ll = 1; ll <= n; ++ll)
        {
            if( ll != icol)
            {
                rdum = a(ll, icol);
                a(ll, icol) = 0.0;
                for(int l = 1; l <= n; ++l)
                {
                    a(ll, l) = a(ll, l) - a(icol, l) * rdum;
                }
                for(int l = 1; l <= m; ++l)
                {
                    b(ll, l) = b(ll, l) - b(icol, l) * rdum;
                }
            }
        }// !for ll
    }// !for i

    for(int l = n; l >= 1; --l)
    {
        if(indxr(l) != indxc(l))
        {
            for(int k = 1; k <= n; ++k)
            {
                rdum = a(k, indxr(l));
                a(k, indxr(l)) = a(k, indxc(l));
                a(k, indxc(l)) = rdum;
            }
        }
    }// !for l
}

CompactDifferenceSecondDerivative::~CompactDifferenceSecondDerivative()
{
    delete lowerBoundofLHSCoefMatrix; lowerBoundofLHSCoefMatrix = NULL;
    delete upperBoundofLHSCoefMatrix; upperBoundofLHSCoefMatrix = NULL;
    delete lowerBoundofRHSCoefMatrix; lowerBoundofRHSCoefMatrix = NULL;
    delete upperBoundofRHSCoefMatrix; upperBoundofRHSCoefMatrix = NULL;

    delete lHSCoefMatrix; lHSCoefMatrix = NULL;
    delete rHSCoefMatrix; rHSCoefMatrix = NULL;
    delete coefd2L; coefd2L = NULL;
    delete coefd2R; coefd2R = NULL;
    delete coefLUDecp; coefLUDecp = NULL;

    delete coefOfLHSBoundary0; coefOfLHSBoundary0 = NULL;
    delete coefOfLHSBoundaryN; coefOfLHSBoundaryN = NULL;
    delete coefOfRHSBoundary0; coefOfRHSBoundary0 = NULL;
    delete coefOfRHSBoundaryN; coefOfRHSBoundaryN = NULL;
}
}

