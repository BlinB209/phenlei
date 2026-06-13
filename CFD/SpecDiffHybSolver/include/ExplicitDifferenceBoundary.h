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
//! @file      ExplicitDifferenceBoundary.h
//! @brief     Explicit difference bondary condition data for SpecDiffHybSolver.
//! @author    ZhangZipei.

#pragma once
#include "TypeDefine.h"
using namespace std;

namespace PHSPACE
{

	class ExplicitDifferenceBoundary 
	{
	public:
		ExplicitDifferenceBoundary();
		~ExplicitDifferenceBoundary();

	public:
		int nBoundary_Neumann;
		RDouble1D *boundary_Neumann_0, *boundary_Neumann_N;

	private:   		
//		double rb0, rb1, rb2, rb3, rb4;
		
	public:
		void GetSchemeCoefDiffBoundaryData();
	};

}



