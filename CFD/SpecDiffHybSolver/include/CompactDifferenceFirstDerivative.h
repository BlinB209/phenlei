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
//! @file      CompactDiffData.h
//! @brief     Compact difference data for SpecDiffHybSolver.
//! @author    ZhangZipei.

#pragma once
#include "TypeDefine.h"
#include "SpecDiffHybGrid.h"
#include "ExplicitDifferenceBoundary.h"
#include "Complex.h"

using namespace std;

namespace PHSPACE
{

	class CompactDifferenceFirstDerivative
	{
	public:
		CompactDifferenceFirstDerivative();
		~CompactDifferenceFirstDerivative();

	private:   
		int nBoundary_Neumann;
		int minofLHSMatrixBound, maxofLHSMatrixBound, minofRHSMatrixBound, maxofRHSMatrixBound;
		Int1D *lowerBoundofLHSCoefMatrix, *upperBoundofLHSCoefMatrix, *lowerBoundofRHSCoefMatrix, *upperBoundofRHSCoefMatrix;
		RDouble2D *lHSCoefMatrix, *rHSCoefMatrix;
		RDouble1D *boundary_Neumann_0, *boundary_Neumann_N;

	public:
        void InitCompactDifferenceFirstDerivativeData(SpecDiffHybGrid *GridData, ExplicitDifferenceBoundary *DiffBoundaryData);
        void GetDetaDz( SpecDiffHybGrid * );
        void GetDvarDz( SpecDiffHybGrid *, Complex3D *, Complex3D *, Complex1D *, Complex1D *, Complex1D *, RDouble1D *);
        void GetDvarDz( const int, const int, RDouble *, RDouble *, RDouble *, RDouble *, RDouble * );
        void GetDvarDz( const int, const int, RDouble *, PHComplex *, PHComplex *, PHComplex *, RDouble * );
        void GetGradientVector( PHComplex *, PHComplex *, PHComplex *, PHComplex *, Complex1D *, Complex1D *, Complex1D *, RDouble1D *, SpecDiffHybGrid * );
        void GetDivergence( PHComplex *, PHComplex *, PHComplex *, Complex3D *, Complex1D *, Complex1D *, Complex3D *, Complex1D *, RDouble1D *, SpecDiffHybGrid * );

	private:	
		void SetCompactDiffData_CarpenterC6(const int, const int);
		void SetCompactDiffData_optC4(const int, const int);
		void SetCompactDiffData_optC4_InnerOnly(const int, const int, ExplicitDifferenceBoundary *);
		void AllocCompactDiffData(const int, const int);
		void GetFirstDerivative(const int, const int, const int, RDouble1D *, RDouble1D *, Int1D *, Int1D *, RDouble2D *, RDouble2D *, RDouble1D *,
			                           RDouble1D *);
		void GetFirstDerivative(const int, const int, const int, Complex1D *, Complex1D *, Int1D *, Int1D *, RDouble2D *, RDouble2D *, Complex1D *,
			                           RDouble1D *);
		static void Tridiag_scalar(const int, RDouble *, RDouble *, RDouble *, RDouble*, RDouble *, RDouble *);
        static void Tridiag_scalar(const int, RDouble *, RDouble *, RDouble *, PHComplex*, PHComplex *, RDouble *);

	};


}



