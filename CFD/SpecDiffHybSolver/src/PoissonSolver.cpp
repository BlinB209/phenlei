#include "PoissonSolver.h"
#include "Precision.h"
#include "GlobalDataBase.h"
#include "CompactDifferenceSecondDerivative.h"
#include "TK_Exit.h"
#include "Complex.h"
#include "PHHeader.h"

using namespace std;

namespace PHSPACE
{
PoissonSolver::PoissonSolver()
{
}

PoissonSolver::~PoissonSolver()
{
    delete coefMatrixPressureNeumann; coefMatrixPressureNeumann = NULL;
    delete coefMatrixVelocityDirichlet; coefMatrixVelocityDirichlet = NULL;
    delete boundaryMatrixPressureNeumann; boundaryMatrixPressureNeumann = NULL;
}

void PoissonSolver::InitPoissonSolverData( SpecDiffHybGrid *GridData, ExplicitDifferenceBoundary *DiffBoundaryData, CompactDifferenceSecondDerivative*   
                                          CompactDifferenceSecondDerivativeData, RDouble1D *realnut, RDouble1D *realnut_00)
{
    AllocPoissonSolverData(GridData, DiffBoundaryData);

    GetCoefMatrixPressureNeumann( GridData, DiffBoundaryData, CompactDifferenceSecondDerivativeData );

    GetCoefMatrixVelocityDirichlet(GridData, DiffBoundaryData, CompactDifferenceSecondDerivativeData, realnut, realnut_00);
}

void PoissonSolver::AllocPoissonSolverData(SpecDiffHybGrid *GridData, ExplicitDifferenceBoundary *DiffBoundaryData)
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex  ;
    int jst = GridData -> jStartFourierIndex;
    int jed = GridData -> jEndFourierIndex  ;
    int kst = GridData -> kStartFourierIndex;
    int ked = GridData -> kEndFourierIndex  ;
    int nbc = DiffBoundaryData -> nBoundary_Neumann;

    Range I(ist, ied);
    Range J(jst, jed);
    Range K(kst, ked);
    Range N(1, nbc);
    Range NN(-2, 2);

    coefMatrixPressureNeumann     = new RDouble4D(I, NN, J, K, fortranArray);
    coefMatrixVelocityDirichlet   = new RDouble4D(I, NN, J, K, fortranArray);
    boundaryMatrixPressureNeumann = new RDouble4D(N, N, J, K, fortranArray);
}

void PoissonSolver::PoissonSolverNeumannBC( const int n, RDouble *coef_LHS, PHComplex *var, int *lowerBound, int *upperBound, const int nbc, const PHComplex boundN, 
                                           const PHComplex bound0, PHComplex *rHS, RDouble *coef, RDouble *boundary )
{
    RDouble2D rA( coef_LHS, Range(1, n-1), Range(-2, 2), neverDeleteData, fortranArray );
    Complex1D cf( var, Range(0, n), neverDeleteData, fortranArray );
    Int1D il( lowerBound, Range(1, n-1), neverDeleteData, fortranArray );
    Int1D ir( upperBound, Range(1, n-1), neverDeleteData, fortranArray );
    Complex1D cR( rHS, Range(0, n), neverDeleteData, fortranArray );
    RDouble2D rC( coef, Range(0, n), Range(-2, 2), neverDeleteData, fortranArray );
    RDouble2D bca( boundary, Range(1, nbc), Range(1, nbc), neverDeleteData, fortranArray );

    int np1 = n + 1;
    cR(Range(0, n)) = PHComplex(0.0, 0.0);
    for ( int i = 1; i <= n-1; ++i )
    {
        for ( int idp = il(i); idp <= ir(i); ++idp )
        {
            int iip = i + idp;
            cR(i) = cR(i) + rA(i, idp) * cf(iip);
        }
    }
    cR(0) = bound0;
    cR(n) = boundN;

    LUSolveNeumannBC( np1, &rC(0, -2), nbc, &bca(1, 1), &cR(0) );

    cf( Range(0, n) ) = cR( Range(0, n) );
}

void PoissonSolver::PoissonSolverDirichletBC( const int n, RDouble *coef_LHS, RDouble *coef_RHS, PHComplex *var, const RDouble rK2, int *lowerBound, int *upperBound, 
                                             const PHComplex boundN, const PHComplex bound0, RDouble *coef, PHComplex *rHS )
{
    RDouble2D rA( coef_LHS, Range(1, n-1), Range(-2, 2), neverDeleteData, fortranArray );
    RDouble2D rB( coef_RHS, Range(1, n-1), Range(-2, 2), neverDeleteData, fortranArray );
    Int1D il( lowerBound, Range(1, n-1), neverDeleteData, fortranArray );
    Int1D ir( upperBound, Range(1, n-1), neverDeleteData, fortranArray );
    Complex1D cf( var, Range(0, n), neverDeleteData, fortranArray );
    RDouble2D rC( coef, Range(0, n), Range(-2, 2), neverDeleteData, fortranArray );
    Complex1D cR( rHS, Range(0, n), neverDeleteData, fortranArray );

    int np1 = n + 1;
    cR( Range(0, n) ) = PHComplex(0.0, 0.0);
    for ( int i = 1; i <= n-1; ++i )
    {
        for ( int idp = il(i); idp <= ir(i); ++idp )
        {
            int iip = i + idp;
            cR(i) = cR(i) + rA(i, idp) * cf(iip);
        }
    }

    rC( Range(0, n) ) = 0.0;
    for ( int i = 1; i <= n-1; ++i )
    {
        for ( int idp = il(i); idp <= ir(i); ++idp )
        {
            int iip = i + idp;
            rC(i, idp) = rB(i, idp) - rA(i, idp) * rK2;
        }
    }

    rC(0, -2) = 0.0;
    rC(0, -1) = 0.0;
    rC(0, 0) = 1.0;
    rC(0, 1) = 0.0;
    rC(0, 2) = 0.0;
    rC(n,-2) = 0.0;
    rC(n,-1) = 0.0;
    rC(n, 0) = 1.0;
    rC(n, 1) = 0.0;
    rC(n, 2) = 0.0;
    cR(0) = bound0;
    cR(n) = boundN;

    LUDecompose(np1, &rC(0, -2));
    LUSolve(np1, &rC(0, -2), &cR(0));

    cf(Range(0, n)) = cR(Range(0, n));
}

void PoissonSolver::PoissonSolverVelocity( const int n, RDouble *coef_LHS, PHComplex *var, int *lowerBound, int *upperBound, const PHComplex boundN, const PHComplex bound0, 
                                          RDouble *realPHI, RDouble *coef, PHComplex *rHS )
{
    RDouble2D rA( coef_LHS, Range(1, n-1), Range(-2, 2), neverDeleteData, fortranArray );
    Complex1D cf( var, Range(0, n), neverDeleteData, fortranArray );
    Int1D il( lowerBound, Range(1, n-1), neverDeleteData, fortranArray );
    Int1D ir( upperBound, Range(1, n-1), neverDeleteData, fortranArray );
    RDouble1D rPHI( realPHI, Range(0, n), neverDeleteData, fortranArray );
    RDouble2D rC( coef, Range(0, n), Range(-2, 2), neverDeleteData, fortranArray );
    Complex1D cR( rHS, Range(0, n), neverDeleteData, fortranArray );

    int np1 = n + 1;
    cR( Range(0, n) ) = PHComplex(0.0, 0.0);
    for ( int i = 1; i <= n-1; ++i )
    {
        for ( int idp = il(i); idp <= ir(i); ++idp )
        {
            int iip = i + idp;
            cR(i) = cR(i) + rA(i, idp) * cf(iip) * rPHI(iip);
        }
    }
    cR(0) = bound0;
    cR(n) = boundN;

    LUSolve( np1, &rC(0, -2), &cR(0) );

    cf( Range(0, n) ) = cR( Range(0, n) );
}

void PoissonSolver::PoissonSolverVelocity(PHComplex &uLocal, const RDouble timeStep, const RDouble k2, const RDouble phi)
{
    RDouble psi = phi / timeStep + k2;

    uLocal = phi * uLocal / psi;
}

void PoissonSolver::GetCoefMatrixPressureNeumann( SpecDiffHybGrid *GridData, ExplicitDifferenceBoundary *DiffBoundaryData, CompactDifferenceSecondDerivative*   
                                                 CompactDifferenceSecondDerivativeData)
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex  ;
    int jst = GridData -> jStartFourierIndex;
    int jed = GridData -> jEndFourierIndex  ;
    int kst = GridData -> kStartFourierIndex;
    int ked = GridData -> kEndFourierIndex  ;
    int nsol = ied - ist;

    (*coefMatrixPressureNeumann)( Range(ist, ied), Range(-2, 2), Range(jst, jed), Range(kst, ked) ) = 0.0;

    int nbc_neumann = DiffBoundaryData -> nBoundary_Neumann;
    (*boundaryMatrixPressureNeumann)( Range(1, nbc_neumann), Range(1, nbc_neumann), Range(jst, jed), Range(kst, ked) ) = 0.0;

    RDouble1D *realKX = GridData -> realKX;
    RDouble1D *realKY = GridData -> realKY;

    RDouble2D *coef_LHS = CompactDifferenceSecondDerivativeData -> lHSCoefMatrix;
    RDouble2D *coef_RHS = CompactDifferenceSecondDerivativeData -> rHSCoefMatrix;

    Int1D *lowerBound = CompactDifferenceSecondDerivativeData -> lowerBoundofRHSCoefMatrix;
    Int1D *upperBound = CompactDifferenceSecondDerivativeData -> upperBoundofRHSCoefMatrix;

    RDouble1D *boundaryNeumann0 = DiffBoundaryData -> boundary_Neumann_0;
    RDouble1D *boundaryNeumannN = DiffBoundaryData -> boundary_Neumann_N;

    for ( int ix = kst; ix <= ked; ++ix )
    {
        for (int iy = jst; iy <= jed; ++iy)
        {
            RDouble kx = (*realKX)(ix);
            RDouble ky = (*realKY)(iy);
            RDouble k2  = kx * kx + ky * ky;
            if( k2 < TINY ) continue;

            DecomCoefMatrixPressureNeumannBC( nsol, &(*coef_LHS)(ist+1, -2), &(*coef_RHS)(ist+1, -2), k2, &(*lowerBound)(ist+1), &(*upperBound)(ist+1), nbc_neumann, 
                &(*boundaryNeumann0)(0), &(*boundaryNeumannN)(-nbc_neumann-1), &(*coefMatrixPressureNeumann)(ist, -2, iy, ix), 
                &(*boundaryMatrixPressureNeumann)(1, 1, iy, ix) );
        }
    }
}

void PoissonSolver::GetCoefMatrixVelocityDirichlet(SpecDiffHybGrid *GridData, ExplicitDifferenceBoundary *DiffBoundaryData, CompactDifferenceSecondDerivative*   
                                                   CompactDifferenceSecondDerivativeData, RDouble1D *realnut, RDouble1D *realnut_00)
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex  ;
    int jst = GridData -> jStartFourierIndex;
    int jed = GridData -> jEndFourierIndex  ;
    int kst = GridData -> kStartFourierIndex;
    int ked = GridData -> kEndFourierIndex  ;
    int nsol = ied - ist;
    Range I(ist, ied);
    Range J(jst, jed);
    Range K(kst, ked);

    RDouble reynolds = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble timeStep = GlobalDataBase::GetDoubleParaFromDB("timeStep");

    (*coefMatrixVelocityDirichlet)(I, Range(-2, 2), J, K) = 0.0;

    RDouble1D *realPHI = new RDouble1D(I, fortranArray);
    (*realPHI)(I) = 2.0 * reynolds/( 1.0 + (*realnut)(I) );

    RDouble1D *realKX = GridData -> realKX;
    RDouble1D *realKY = GridData -> realKY;

    RDouble2D *coef_LHS = CompactDifferenceSecondDerivativeData -> lHSCoefMatrix;
    RDouble2D *coef_RHS = CompactDifferenceSecondDerivativeData -> rHSCoefMatrix;

    Int1D *lowerBound = CompactDifferenceSecondDerivativeData -> lowerBoundofRHSCoefMatrix;
    Int1D *upperBound = CompactDifferenceSecondDerivativeData -> upperBoundofRHSCoefMatrix;

    for ( int ix = kst; ix <= ked; ++ix )
    {
        for (int iy = jst; iy <= jed; ++iy)
        {
            RDouble kx = (*realKX)(ix);
            RDouble ky = (*realKY)(iy);
            RDouble k2  = kx * kx + ky * ky;

            DecomCoefMatrixVelocityDirichletBC( nsol, &(*coef_LHS)(ist+1, -2), &(*coef_RHS)(ist+1, -2), k2, &(*lowerBound)(ist+1), &(*upperBound)(ist+1), timeStep, &(*realPHI)(ist),
                &(*coefMatrixVelocityDirichlet)(ist, -2, iy, ix) );
        }
    }

    int imode00fix = GlobalDataBase::GetIntParaFromDB("imode00fix");
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if(imode00fix > 0 && myid ==0)
    {
        (*realPHI)(I) = 2.0 * reynolds/( 1.0 + (*realnut_00)(Range(ist, ied)) );

        int ix = kst;
        int iy = jst;

        RDouble kx = (*realKX)(ix);
        RDouble ky = (*realKY)(iy);
        RDouble k2  = kx * kx + ky * ky;

        DecomCoefMatrixVelocityDirichletBC( nsol, &(*coef_LHS)(ist+1, -2), &(*coef_RHS)(ist+1, -2), k2, &(*lowerBound)(ist+1), &(*upperBound)(ist+1), timeStep, &(*realPHI)(ist),
            &(*coefMatrixVelocityDirichlet)(ist, -2, iy, ix) );
    }

    delete realPHI; realPHI = NULL;
}

void PoissonSolver::DecomCoefMatrixPressureNeumannBC( const int n, RDouble *coef_LHS, RDouble *coef_RHS, const double rK2, int *lowerBound, int *upperBound, const int nbc, 
                                                     RDouble *boundaryNeumann0, RDouble *boundaryNeumannN, RDouble *coef, RDouble *boundary )
{
    RDouble2D rA( coef_LHS, Range(1, n-1), Range(-2, 2), neverDeleteData, fortranArray );
    RDouble2D rB( coef_RHS, Range(1, n-1), Range(-2, 2), neverDeleteData, fortranArray );
    Int1D il( lowerBound, Range(1, n-1), neverDeleteData, fortranArray );
    Int1D ir( upperBound, Range(1, n-1), neverDeleteData, fortranArray );
    RDouble1D bca_0( boundaryNeumann0, Range(0, nbc+1), neverDeleteData, fortranArray );
    RDouble1D bca_n( boundaryNeumannN, Range(-(nbc+1), 0), neverDeleteData, fortranArray );
    RDouble2D rC( coef, Range(0, n), Range(-2, 2), neverDeleteData, fortranArray );
    RDouble2D bca( boundary, Range(1, nbc), Range(1, nbc), neverDeleteData, fortranArray );

    int np1 = n + 1;
    rC( Range(0, n), Range(-2, 2) ) = 0.0;
    for ( int i = 1; i <= n-1; ++i )
    {
        for( int idp = il(i); idp <= ir(i); ++idp )
        {
            rC(i,idp) = rB(i,idp) - rA(i,idp) * rK2;
        }
    }

    rC(0, 0) = bca_0(0);
    rC(0, 1) = bca_0(1);
    rC(0, 2) = bca_0(2);

    rC(n, 0) = bca_n( 0);
    rC(n,-1) = bca_n(-1);
    rC(n,-2) = bca_n(-2);

    bca( Range(1, nbc), Range(1, nbc) ) = 0.0;
    for ( int i = 2; i <= nbc; ++i )
    {
        bca(1,i) = bca_0(i+1);
    }

    for ( int i = 1; i <= nbc-1; ++i )
    {
        bca(nbc,i) = bca_n(i-nbc-2);
    }

    LUDecomposeNeumannBC( np1, &rC(0, -2), nbc, &bca(1, 1) );
}

void PoissonSolver::DecomCoefMatrixVelocityDirichletBC( const int n, RDouble *coef_LHS, RDouble *coef_RHS, const double rK2, int *lowerBound, int *upperBound, const double rdt, 
                                                       RDouble *realPHI, RDouble *coef )
{
    RDouble2D rA( coef_LHS, Range(1, n-1), Range(-2, 2), neverDeleteData, fortranArray );
    RDouble2D rB( coef_RHS, Range(1, n-1), Range(-2, 2), neverDeleteData, fortranArray );
    Int1D il( lowerBound, Range(1, n-1), neverDeleteData, fortranArray );
    Int1D ir( upperBound, Range(1, n-1), neverDeleteData, fortranArray );
    RDouble1D rPHI( realPHI, Range(0, n), neverDeleteData, fortranArray );
    RDouble2D rC( coef, Range(0, n), Range(-2, 2), neverDeleteData, fortranArray );

    rC( Range(0, n), Range(-2, 2) ) = 0.0;
    for( int i = 1; i <= n-1; ++i )
    {
        for( int idp = il(i); idp <= ir(i); ++idp)
        {
            int iip = i + idp;
            rC(i,idp) = -rB(i,idp) + rA(i,idp) * ( rPHI(iip)/rdt + rK2 );
        }
    }

    //--- Dirichelt B.C. : u(+1)=0,u(-1)=0
    rC(0,-2) = 0.0;
    rC(0,-1) = 0.0;
    rC(0, 0) = 1.0;
    rC(0, 1) = 0.0;
    rC(0, 2) = 0.0;
    rC(n,-2) = 0.0;
    rC(n,-1) = 0.0;
    rC(n, 0) = 1.0;
    rC(n, 1) = 0.0;
    rC(n, 2) = 0.0;

    int np1 = n + 1;

    LUDecompose( np1, &rC(0, -2) );
}

void PoissonSolver::LUDecompose( const int n, RDouble *coef )
{
    RDouble2D a( coef, Range(1, n), Range(1, 5), neverDeleteData, fortranArray );
    int i = 1;

    a(2, 2) = a(2, 2) / a(1, 3);
    a(2, 3) = a(2, 3) - a(2, 2) * a(1, 4);
    a(2, 4) = a(2, 4) - a(2, 2) * a(1, 5);

    for ( i = 3; i <= n-1; ++i)
    {
        a(i, 1) = a(i, 1) / a(i-2, 3);
        a(i, 2) = ( a(i, 2) - a(i, 1) * a(i-2, 4) ) / a(i-1, 3);
        a(i, 3) = a(i, 3) - a(i, 1) * a(i-2, 5) - a(i, 2) * a(i-1, 4);
        a(i, 4) = a(i, 4) - a(i, 2) * a(i-1, 5);
    }

    i = n;
    a(i, 1) = a(i, 1) / a(i-2, 3);
    a(i, 2) = ( a(i, 2) - a(i, 1) * a(i-2, 4) ) / a(i-1, 3);
    a(i, 3) = a(i, 3) - a(i, 1) * a(i-2, 5) - a(i, 2) * a(i-1, 4);
}

void PoissonSolver::LUDecomposeNeumannBC( const int n, RDouble *coef, int nbc, RDouble *boundary )
{
    RDouble2D a( coef, Range(1, n), Range(1, 5), neverDeleteData, fortranArray );
    RDouble2D bca( boundary, Range(1, nbc), Range(1, nbc), neverDeleteData, fortranArray );
    int i, j;

    if ( n <= (2*nbc+2) )
    {
        TK_Exit::ExceptionExit("n<=2*nbc+2 in LUDecomposeNeumannBC");
    }

    a(2, 2) = a(2, 2) / a(1, 3);
    a(2, 3) = a(2, 3) - a(2, 2) * a(1, 4);
    a(2, 4) = a(2, 4) - a(2, 2) * a(1, 5);

    //bca 对角线不用
    //  bca(1,2:6) 包含 Neumann 边界系数，分解后不变.
    //  bca(6,1:5) 包含 Neumann 边界系数，分解后不变
    //  其余初始值为零

    a(2, 5) = a(2, 5) - a(2, 2) * bca(1, 2);
    for ( i = 3; i <= nbc; ++i )
    {
        bca(2, i) = bca(2, i) - a(2, 2) * bca(1, i);
    }

    // i=3 ! 对 nbc=3也适用

    for ( i = 3; i <= nbc-1; ++i )
    {
        a(i, 1) = a(i, 1) / a(i-2, 3);
        a(i, 2) = ( a(i, 2) - a(i, 1) * a(i-2, 4) ) / a(i-1, 3);
        a(i, 3) = a(i, 3) - a(i, 1) * a(i-2, 5) - a(i, 2) * a(i-1, 4);
        //different parts:
        a(i, 4) = a(i, 4) - a(i, 2) * a  (i-1, 5) - a(i, 1) * bca(i-2, i-1);
        a(i, 5) = a(i, 5) - a(i, 2) * bca(i-1, i) - a(i, 1) * bca(i-2, i  );
        for ( j = i+1; j <= nbc; ++j )
        {
            bca(i, j) = bca(i, j) - a(i, 2) * bca(i-1, j) - a(i, 1) * bca(i-2, j);
        }
    }

    i = nbc;
    a(i, 1) = a(i, 1) / a(i-2, 3);
    a(i, 2) = ( a(i, 2) - a(i, 1) * a(i-2, 4) ) / a(i-1, 3);
    a(i, 3) = a(i, 3) - a(i, 1) * a(i-2, 5) - a(i, 2) * a(i-1, 4);
    // different parts:
    a(i, 4) = a(i, 4) - a(i, 2) * a  (i-1, 5) - a(i, 1) * bca(i-2, i-1);
    a(i, 5) = a(i, 5) - a(i, 2) * bca(i-1, i) - a(i, 1) * bca(i-2, i  );

    i = nbc + 1;
    a(i, 1) = a(i, 1) / a(i-2, 3);
    a(i, 2) = ( a(i, 2) - a(i, 1) * a(i-2, 4) ) / a(i-1, 3);
    a(i, 3) = a(i, 3) - a(i, 1) * a(i-2, 5) - a(i, 2) * a(i-1, 4);
    // different parts:
    a(i, 4) = a(i, 4) - a(i, 2) * a(i-1, 5) - a(i, 1) * bca(i-2, i-1);

    //inner points:
    for ( i = nbc+2; i <= n-1; ++i )
    {
        a(i, 1) = a(i, 1) / a(i-2, 3);
        a(i, 2) = ( a(i, 2) - a(i, 1) * a(i-2, 4) ) / a(i-1, 3);
        a(i, 3) = a(i, 3) - a(i, 1) * a(i-2, 5) - a(i, 2) * a(i-1, 4);
        a(i, 4) = a(i, 4) - a(i, 2) * a(i-1, 5);
    }

    i = n;
    bca(nbc, 1) = bca(nbc, 1) / a(n-nbc-1, 3);
    bca(nbc, 2) = ( bca(nbc, 2) - bca(nbc, 1) * a(n-nbc-1, 4) ) / a(n-nbc,3);
    for ( j = 3; j <= nbc-1; ++j )
    {
        bca(nbc, j) = ( bca(nbc, j) - bca(nbc, j-2) * a(n-nbc-4+j, 5) - bca(nbc, j-1) * a(n-nbc-3+j, 4) ) / a(n-nbc-2+j, 3);
    }

    a(n, 1) = ( a(n, 1) - bca(nbc, nbc-2) * a(n-4, 5) - bca(nbc, nbc-1) * a(n-3, 4) ) / a(n-2, 3);
    a(n, 2) = ( a(n, 2) - bca(nbc, nbc-1) * a(n-3, 5) - a  (n  ,     1) * a(n-2, 4) ) / a(n-1, 3);
    a(n, 3) = a(n, 3) - a(n, 1) * a(n-2, 5) - a(n, 2) * a(n-1, 4);
}

void PoissonSolver::LUSolve( const int n, RDouble *coef, RDouble *rHS )
{
    RDouble2D a( coef, Range(1, n), Range(1, 5), neverDeleteData, fortranArray );
    RDouble1D b( rHS, Range(1, n), neverDeleteData, fortranArray );

    b(2) = b(2) - a(2, 2) * b(1);
    for ( int i = 3; i <= n; ++i )
    {
        b(i) = b(i) - a(i, 1) * b(i-2) - a(i, 2) * b(i-1);
    }

    // !back substitution
    b(n) = b(n) / a(n, 3);
    b(n-1) = ( b(n-1) - a(n-1, 4) * b(n) ) / a(n-1, 3);
    for ( int i = n-2; i >= 1; --i )
    {
        b(i) = ( b(i) - a(i, 4) * b(i+1) - a(i, 5) * b(i+2) ) / a(i, 3);
    }
}

void PoissonSolver::LUSolve( const int n, RDouble *coef, PHComplex *rHS )
{
    RDouble2D a( coef, Range(1, n), Range(1, 5), neverDeleteData, fortranArray );
    Complex1D b( rHS, Range(1, n), neverDeleteData, fortranArray );

    // !forward substitution
    b(2) = b(2) - a(2, 2) * b(1);
    for ( int i = 3; i <= n; ++i )
    {
        b(i) = b(i) - a(i, 1) * b(i-2) - a(i, 2) * b(i-1);
    }

    // !back substitution
    b(n) = b(n) / a(n, 3);
    b(n-1) = ( b(n-1) - a(n-1, 4) * b(n) ) / a(n-1, 3);
    for ( int i = n-2; i >= 1; --i )
    {
        b(i) = ( b(i) - a(i, 4) * b(i+1) - a(i, 5) * b(i+2) ) / a(i, 3);
    }
}

void PoissonSolver::LUSolveNeumannBC( const int n, RDouble *coef, const int nbc, RDouble *boundary, PHComplex *rHS )
{
    RDouble2D a( coef, Range(1, n), Range(1, 5), neverDeleteData, fortranArray );
    RDouble2D bca( boundary, Range(1, nbc), Range(1, nbc), neverDeleteData, fortranArray);
    Complex1D cb( rHS, Range(1, n), neverDeleteData, fortranArray);

    if ( n <= (2*nbc+2) )
    {
        TK_Exit::ExceptionExit("n<=2*nbc+2 in LUSolveNeumannBC");
    }

    // !forward substitution
    cb(2) = cb(2) - a(2, 2) * cb(1);
    for ( int i = 3; i <= n-1; ++i )
    {
        cb(i) = cb(i) - a(i, 1) * cb(i-2) - a(i, 2) * cb(i-1);
    }
    cb(n) = cb(n) - a(n, 1) * cb(n-2) - a(n, 2) * cb(n-1);
    for ( int j = 1; j <= (nbc-1); ++j )
    {
        cb(n) = cb(n) - bca(nbc, j) * cb(n-nbc-2+j);
    }

    // !back substitution
    cb(n) = cb(n) / a(n, 3);
    cb(n-1) = ( cb(n-1) - a(n-1, 4) * cb(n) ) / a(n-1, 3);
    for ( int i = (n-2); i >= nbc; --i )
    {
        cb(i) = ( cb(i) - a(i, 4) * cb(i+1) - a(i, 5) * cb(i+2) ) / a(i ,3);
    }
    for ( int i = (nbc-1); i >= 1; --i )
    {
        cb(i) = ( cb(i) - a(i, 4) * cb(i+1) - a(i, 5) * cb(i+2) );
        for ( int j = (i+1); j <= nbc; ++j )
        {
            cb(i) = cb(i) - bca(i, j) * cb(j+2);
        }
        cb(i) = cb(i) / a(i, 3);
    }
}

}

