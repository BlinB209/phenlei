#include "CorrectWallDivergence.h"
#include "GlobalDataBase.h"
#include "CompactDifferenceFirstDerivative.h"
#include "CompactDifferenceSecondDerivative.h"
#include "PoissonSolver.h"
#include "PHMpi.h"
#include "TK_Exit.h"
#include "PHHeader.h"

using namespace std;

namespace PHSPACE
{
CorrectWallDivergence::CorrectWallDivergence()
{
}

CorrectWallDivergence::~CorrectWallDivergence()
{
    delete complexDeltaU1; complexDeltaU1 = NULL;
    delete complexDeltaV1; complexDeltaV1 = NULL;
    delete complexDeltaW1; complexDeltaW1 = NULL;
    delete complexDeltaP1; complexDeltaP1 = NULL;
    delete complexDeltaU2; complexDeltaU2 = NULL;
    delete complexDeltaV2; complexDeltaV2 = NULL;
    delete complexDeltaW2; complexDeltaW2 = NULL;
    delete complexDeltaP2; complexDeltaP2 = NULL;
    delete complexDivDeltaU1LowerBound; complexDivDeltaU1LowerBound = NULL;
    delete complexDivDeltaU1UpperBound; complexDivDeltaU1UpperBound = NULL;
    delete complexDivDeltaU2LowerBound; complexDivDeltaU2LowerBound = NULL;
    delete complexDivDeltaU2UpperBound; complexDivDeltaU2UpperBound = NULL;
}

void CorrectWallDivergence::InitCorrectWallDivergenceData( SpecDiffHybGrid *GridData, CompactDifferenceFirstDerivative *CompactDifferenceFirstDerivativeData,
                                                          CompactDifferenceSecondDerivative *CompactDifferenceSecondDerivativeData, PoissonSolver *PoissonSolverData, 
                                                          RDouble1D *realnut )
{
    AllocCorrectWallDivergenceData(GridData);

    SetCorrectWallDivergenceData(GridData, CompactDifferenceFirstDerivativeData, CompactDifferenceSecondDerivativeData, PoissonSolverData, realnut);
}

void CorrectWallDivergence::AllocCorrectWallDivergenceData( SpecDiffHybGrid *GridData )
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex  ;
    int jst = GridData -> jStartFourierIndex;
    int jed = GridData -> jEndFourierIndex  ;
    int kst = GridData -> kStartFourierIndex;
    int ked = GridData -> kEndFourierIndex  ;

    Range I(ist, ied);
    Range J(jst, jed);
    Range K(kst, ked);

    complexDeltaU1 = new Complex3D(I, J, K, fortranArray);
    complexDeltaV1 = new Complex3D(I, J, K, fortranArray);
    complexDeltaW1 = new Complex3D(I, J, K, fortranArray);
    complexDeltaP1 = new Complex3D(I, J, K, fortranArray);

    complexDeltaU2 = new Complex3D(I, J, K, fortranArray);
    complexDeltaV2 = new Complex3D(I, J, K, fortranArray);
    complexDeltaW2 = new Complex3D(I, J, K, fortranArray);
    complexDeltaP2 = new Complex3D(I, J, K, fortranArray);

    complexDivDeltaU1LowerBound = new Complex2D(J, K, fortranArray);
    complexDivDeltaU1UpperBound = new Complex2D(J, K, fortranArray);
    complexDivDeltaU2LowerBound = new Complex2D(J, K, fortranArray);
    complexDivDeltaU2UpperBound = new Complex2D(J, K, fortranArray);
}

void CorrectWallDivergence::SetCorrectWallDivergenceData( SpecDiffHybGrid *GridData, CompactDifferenceFirstDerivative *CompactDifferenceFirstDerivativeData,
                                                         CompactDifferenceSecondDerivative *CompactDifferenceSecondDerivativeData, PoissonSolver *PoissonSolverData, 
                                                         RDouble1D *realnut )
{
    const PHComplex cOne = PHComplex(1.0, 0.0);
    const PHComplex cZero = PHComplex(0.0, 0.0);

    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex  ;
    int jst = GridData -> jStartFourierIndex;
    int jed = GridData -> jEndFourierIndex  ;
    int kst = GridData -> kStartFourierIndex;
    int ked = GridData -> kEndFourierIndex  ;

    RDouble reynolds = GlobalDataBase::GetDoubleParaFromDB("refReNumber");

    int nSolve = ied - ist;

    Range I(ist , ied);
    RDouble1D *realPHI = new RDouble1D( I , fortranArray );
    (*realPHI)(I) = 2.0 * reynolds/( 1.0 + (*realnut)(I) );

    RDouble1D *rkx = GridData -> realKX;
    RDouble1D *rky = GridData -> realKY;
    Complex1D *ckx = GridData -> complexKX;
    Complex1D *cky = GridData -> complexKY;
    RDouble1D *detaDz = GridData -> realDetaDz;

    Complex1D *complexP = new Complex1D(I, fortranArray);
    Complex1D *complexU = new Complex1D(I, fortranArray);
    Complex1D *complexV = new Complex1D(I, fortranArray);
    Complex1D *complexW = new Complex1D(I, fortranArray);
    Complex1D *complexdPdZ = new Complex1D(I, fortranArray);
    Complex1D *cR = new Complex1D(I, fortranArray);

    RDouble2D *coef_LHS = CompactDifferenceSecondDerivativeData -> lHSCoefMatrix;
    RDouble2D *coef_RHS = CompactDifferenceSecondDerivativeData -> rHSCoefMatrix;
    Int1D *lowerBound = CompactDifferenceSecondDerivativeData -> lowerBoundofRHSCoefMatrix;
    Int1D *upperBound = CompactDifferenceSecondDerivativeData -> upperBoundofRHSCoefMatrix;

    RDouble4D *coefVelocity = PoissonSolverData -> coefMatrixVelocityDirichlet;
    Range J(-2 , 2);
    RDouble2D *rC = new RDouble2D(I, J, fortranArray);

    for ( int ix = kst; ix <= ked; ++ix )
    {
        for ( int iy = jst; iy <= jed; ++iy )
        {
            RDouble kxx = (*rkx)(ix);
            RDouble kyy = (*rky)(iy);
            RDouble k2  = kxx * kxx + kyy * kyy;
            if( k2 < TINY ) continue;

            //! To correct wall divergence using influence-matrix technique
            //! To solve p1 and u1
            (*complexP)(I) = cZero;
            PoissonSolver::PoissonSolverDirichletBC( nSolve, &(*coef_LHS)(ist+1, -2), &(*coef_RHS)(ist+1, -2), &(*complexP)(ist), k2, &(*lowerBound)(ist+1), &(*upperBound)(ist+1), 
                cOne, cZero, &(*rC)(ist, -2), &(*cR)(ist) );
            CompactDifferenceFirstDerivativeData -> GetDvarDz( ist, ied, &(*detaDz)(ist), &(*complexP)(ist), &(*complexdPdZ)(ist), &(*cR)(ist), &(*rC)(ist, -2) );

            for( int iz = ist; iz <= ied; ++iz )
            {
                (*complexU)(iz) = -(*ckx)(ix) * (*complexP)(iz);
                (*complexV)(iz) = -(*cky)(iy) * (*complexP)(iz);
                (*complexW)(iz) = -(*complexdPdZ)(iz);
            }

            PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(ist+1, -2), &(*complexU)(ist),  &(*lowerBound)(ist+1), &(*upperBound)(ist+1), cZero, cZero, &(*realPHI)(ist),
                &(*coefVelocity)(ist, -2, iy, ix), &(*cR)(ist) );
            PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(ist+1, -2), &(*complexV)(ist),  &(*lowerBound)(ist+1), &(*upperBound)(ist+1), cZero, cZero, &(*realPHI)(ist),
                &(*coefVelocity)(ist, -2, iy, ix), &(*cR)(ist) );
            PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(ist+1, -2), &(*complexW)(ist),  &(*lowerBound)(ist+1), &(*upperBound)(ist+1), cZero, cZero, &(*realPHI)(ist),
                &(*coefVelocity)(ist, -2, iy, ix), &(*cR)(ist) );

            (*complexW)(ist) = cZero;
            (*complexW)(ied) = cZero;
            (*complexU)(ist) = cZero;
            (*complexU)(ied) = cZero;
            (*complexV)(ist) = cZero;
            (*complexV)(ied) = cZero;

            (*complexDeltaP1)(I, iy, ix) = (*complexP)(I);
            (*complexDeltaU1)(I, iy, ix) = (*complexU)(I);
            (*complexDeltaV1)(I, iy, ix) = (*complexV)(I);
            (*complexDeltaW1)(I, iy, ix) = (*complexW)(I);

            CompactDifferenceFirstDerivativeData -> GetDvarDz( ist, ied, &(*detaDz)(ist), &(*complexW)(ist), &(*complexdPdZ)(ist), &(*cR)(ist), &(*rC)(ist, -2) );
            (*complexDivDeltaU1LowerBound)(iy, ix) = (*complexdPdZ)(ist);
            (*complexDivDeltaU1UpperBound)(iy, ix) = (*complexdPdZ)(ied);

            // !To solve p2 and u2
            (*complexP)(I) = cZero;
            PoissonSolver::PoissonSolverDirichletBC( nSolve, &(*coef_LHS)(ist+1, -2), &(*coef_RHS)(ist+1, -2), &(*complexP)(ist), k2, &(*lowerBound)(ist+1), &(*upperBound)(ist+1), 
                cZero, cOne, &(*rC)(ist, -2), &(*cR)(ist) );
            CompactDifferenceFirstDerivativeData -> GetDvarDz( ist, ied, &(*detaDz)(ist), &(*complexP)(ist), &(*complexdPdZ)(ist), &(*cR)(ist), &(*rC)(ist, -2) );
            for( int iz = ist; iz <= ied; ++iz )
            {
                (*complexU)(iz) = -(*ckx)(ix) * (*complexP)(iz);
                (*complexV)(iz) = -(*cky)(iy) * (*complexP)(iz);
                (*complexW)(iz) = -(*complexdPdZ)(iz);
            }

            PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(ist+1, -2), &(*complexU)(ist), &(*lowerBound)(ist+1), &(*upperBound)(ist+1), cZero, cZero, &(*realPHI)(ist),
                &(*coefVelocity)(ist, -2, iy, ix), &(*cR)(ist) );
            PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(ist+1, -2), &(*complexV)(ist), &(*lowerBound)(ist+1), &(*upperBound)(ist+1), cZero, cZero, &(*realPHI)(ist),
                &(*coefVelocity)(ist, -2, iy, ix), &(*cR)(ist) );
            PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(ist+1, -2), &(*complexW)(ist), &(*lowerBound)(ist+1), &(*upperBound)(ist+1), cZero, cZero, &(*realPHI)(ist),
                &(*coefVelocity)(ist, -2, iy, ix), &(*cR)(ist) );

            (*complexW)(ist) = cZero;
            (*complexW)(ied) = cZero;
            (*complexU)(ist) = cZero;
            (*complexU)(ied) = cZero;
            (*complexV)(ist) = cZero;
            (*complexV)(ied) = cZero;

            (*complexDeltaP2)(I, iy, ix) = (*complexP)(I);
            (*complexDeltaU2)(I, iy, ix) = (*complexU)(I);
            (*complexDeltaV2)(I, iy, ix) = (*complexV)(I);
            (*complexDeltaW2)(I, iy, ix) = (*complexW)(I);

            CompactDifferenceFirstDerivativeData -> GetDvarDz( ist, ied, &(*detaDz)(ist), &(*complexW)(ist), &(*complexdPdZ)(ist), &(*cR)(ist), &(*rC)(ist, -2) );
            (*complexDivDeltaU2LowerBound)(iy, ix) = (*complexdPdZ)(ist);
            (*complexDivDeltaU2UpperBound)(iy, ix) = (*complexdPdZ)(ied);
        }
    }

    delete complexP; complexP = NULL;
    delete complexU; complexU = NULL;
    delete complexV; complexV = NULL;
    delete complexW; complexW = NULL;
    delete complexdPdZ; complexdPdZ = NULL;
    delete cR; cR = NULL;
    delete rC; rC = NULL;
    delete realPHI; realPHI = NULL;
}

void CorrectWallDivergence::rotProjection()
{
}

}

