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
//! @file      SpecGrid.h
//! @brief     Grid data for SpecSolver.
//! @author    ZhangZipei.

#pragma once
#include "TypeDefine.h"
#include "Complex.h"
using namespace std;

namespace PHSPACE
{
class SpecGrid
{
public:
    SpecGrid();
    ~SpecGrid();
protected:   
    int nyquistX, nyquistY, nyquistZ;		
    Int1D *physicalStartIndex, *physicalEndIndex, *physicalIndexSize;
    //		Int1D *fourierStartIndexGlobal, *fourierEndIndexGlobal;
    Int1D *dimensionOfP3DFFT;		
    string gridFileName, zgridFileName;
    int procIDWithKx0Ky0;

public:
    int nX, nY, nZ, nX2, nY2, nZ2, nXYAll, nXYZAll, nXYZFourier, nXYZPhysical;
    int iStartFourierIndex , jStartFourierIndex , kStartFourierIndex , iEndFourierIndex , jEndFourierIndex , kEndFourierIndex ;
    int iStartPhysicalIndex, jStartPhysicalIndex, kStartPhysicalIndex, iEndPhysicalIndex, jEndPhysicalIndex, kEndPhysicalIndex;
    //int x0Fall, y0Fall, z0Fall, xNFall, yNFall, zNFall;
    int iStartFourierIndexGlobal, jStartFourierIndexGlobal, kStartFourierIndexGlobal, iEndFourierIndexGlobal, jEndFourierIndexGlobal, kEndFourierIndexGlobal;
    RDouble realKNyquistX, realKNyquistY, realKNyquistZ;

    Int1D *fourierStartIndex, *fourierEndIndex, *fourierIndexSize;
    Int2D *startAndEndFourierIndexOfAll, *startAndEndPhysicalIndexOfAll;
    //		int **startAndEndAll, **startAndEndPhysical;
    //				Int1D * physicalStart, ;
    RDouble1D *realKX, *realKY, *realKZ;
    Complex1D *complexKX, *complexKY, *complexKZ;
    RDouble1D *realZ, *realDetaDz;

protected:
    virtual void InitGridFileName();

    virtual void InitP3DFFT();

    virtual void AllocGridData();

    virtual void InitWavenumber();

public:
    virtual void InitGridData();

    void InitGridSize();

    int GetProcIDWithKx0Ky0();
};

}
