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
//! @file      SpecDiffHybFieldView.h
//! @brief     Field view of Spectral-Difference hybrid method.
//! @author    Zhang Zipei.

#pragma once

#include "TypeDefine.h"
#include "Complex.h"
using namespace std;

namespace PHSPACE
{
	class SpecDiffHybFieldView 
	{
	public:
		SpecDiffHybFieldView();
		~SpecDiffHybFieldView();

	private:
		Complex3D *complexU, *complexV, *complexW, *complexP;
        RDouble3D *realU, *realV, *realW, *realP;
        Complex2D *complexIn;
        RDouble2D *realOut;
        RDouble1D *realX, *realY, *realZ;
        Int2D *fourierStartIndex, *fourierEndIndex;
        int nTProcessor;
        int nx, ny, nz;

	public:
		void Run();

	private:
		void Init();
        void AllocFieldData();
        void InputGridDistribute();
        void InputVelocity();
        void InputPressure();
        void FieldView();
        void OutputGrid();
        void OutputField();
	};

}