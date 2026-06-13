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
//! @file      Complex.h
//! @brief     Complex array.
//! @author    Zhang Zipei.

#pragma once
#include "Precision.h"
#include "DataStruct_Array.h"
#include "DataStruct_AdtTree.h"
#include <string>
#include <vector>
#include <map>
#include <numeric> 
#include <complex> 

using namespace std;

namespace PHSPACE
{
	typedef complex<double> PHComplex;

	typedef PHArray<complex<double>, 1> Complex1D;
	typedef PHArray<complex<double>, 2> Complex2D;
	typedef PHArray<complex<double>, 3> Complex3D;
	typedef PHArray<complex<double>, 4> Complex4D;
	typedef PHArray<complex<double>, 5> Complex5D;
	typedef PHArray<complex<double>, 6> Complex6D;
}

