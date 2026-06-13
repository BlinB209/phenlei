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
//! @file      CompacDiff2ndNonUnifInnerData.h
//! @brief     Compact difference scheme for approximation of second derivatives on non-uniform meshes for SpecDiffHybSolver.
//! @author    ZhangZipei.

#pragma once
#include "TypeDefine.h"
#include "SpecDiffHybGrid.h" 
using namespace std;

namespace PHSPACE
{

	class CompactDifferenceSecondDerivative 
	{
	public:
		CompactDifferenceSecondDerivative();
		~CompactDifferenceSecondDerivative();

	public:
		RDouble2D *lHSCoefMatrix, *rHSCoefMatrix, *coefd2L, *coefd2R, *coefLUDecp;
        RDouble1D *coefOfLHSBoundary0, *coefOfLHSBoundaryN, *coefOfRHSBoundary0, *coefOfRHSBoundaryN;
		Int1D *lowerBoundofLHSCoefMatrix, *upperBoundofLHSCoefMatrix, *lowerBoundofRHSCoefMatrix, *upperBoundofRHSCoefMatrix;

	private:   

//		double **fiveDiagonalMatrixLUDecp;
//		double *leftBoundofLHSCoefMatrix, *rightBoundofLHSCoefMatrix, *leftBoundoRHSCoefMatrix, *rightBoundofRHSCoefMatrix;

	public:
		void InitCompactDifferenceSecondDerivativeData( SpecDiffHybGrid *GridData );
        void GetLaplacian(SpecDiffHybGrid *GridData, Complex3D *cVar, Complex3D *cLaplacianVar, Complex1D *c1tmp, Complex1D *c2tmp, Complex3D *cd2Var);
        void GetD2varDz2(SpecDiffHybGrid *GridData, Complex3D *cVar, Complex3D *cdVardx_2rd, Complex1D *cV, Complex1D *cd2Vardz);
        void GetSecondDerivative(SpecDiffHybGrid *GridData, Complex1D *cVar, Complex1D *cd2Vardz);

	private:
        void AllocCompactDifferenceSecondDerivativeData(SpecDiffHybGrid *);
        void SetCompactDifferenceSecondDerivativeData(SpecDiffHybGrid *);
		void GetCoefMatrix( const int mz, RDouble * ben , RDouble * aln , RDouble * aaa , RDouble * alp , RDouble * bep ,
		                                                                 RDouble * an2 , RDouble * an1 , RDouble * a00 , RDouble * ap1 , RDouble * ap2 , RDouble * x );
        void Gaussj(RDouble *, const int, const int, RDouble *, const int, const int);

	};


}



