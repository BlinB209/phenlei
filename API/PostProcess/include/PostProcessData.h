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
//! @file      PostProcessData.h
//! @brief     Flow post-processing data store.
//! @author    Baka.

#pragma once
#include "Precision.h"

using namespace std;

namespace PHSPACE
{

class PostProcessData
{
private:
    int numberofSolidWallPart;

    //! Number of total part, is numberofSolidWallPart plus ONE.
    //! 'ONE' represent the total solid wall force, stored on the first position.
    int numberofTotalPart;

    //! The corresponding Boundary Condition ID of each part.
    int *partBCID;

    RDouble *cfx, *cfy, *cfz;
    RDouble *cpx, *cpy, *cpz;
    RDouble *cmx, *cmy, *cmz;
    RDouble *clTot, *cdTot, *cd_pr, *cd_sf, *cd_cl2pa, *xcp, *side;
    RDouble *hingeMoment;

public:
    PostProcessData();
    ~PostProcessData();

public:
    void SetNWallPart(int numberofSolidWallPartIn);
    void SetNPart(int numberofTotalPartIn);
    void SetPartBCID(int *partBCIDIn);
    void SetCfx(RDouble *cfxIn);
    void SetCfy(RDouble *cfyIn);
    void SetCfz(RDouble *cfzIn);
    void SetCpx(RDouble *cpxIn);
    void SetCpy(RDouble *cpyIn);
    void SetCpz(RDouble *cpzIn);
    void SetCmx(RDouble *cmxIn);
    void SetCmy(RDouble *cmyIn);
    void SetCmz(RDouble *cmzIn);
    void SetClTot(RDouble *clTotIn);
    void SetCdTot(RDouble *cdTotIn);
    void SetCdPr(RDouble *cd_prIn);
    void SetCdSf(RDouble *cd_sfIn);
    void SetCdCl2Pa(RDouble *cd_cl2paIn);
    void SetXcp(RDouble *xcpIn);
    void SetSide(RDouble *sideIn);
    void SetHingeMoment(RDouble *hingeMomentIn);

    int GetNWallPart();
    int GetNPart();
    int * GetPartBCID();
    RDouble * GetCfx();
    RDouble * GetCfy();
    RDouble * GetCfz();
    RDouble * GetCpx();
    RDouble * GetCpy();
    RDouble * GetCpz();
    RDouble * GetCmx();
    RDouble * GetCmy();
    RDouble * GetCmz();
    RDouble * GetClTot();
    RDouble * GetCdTot();
    RDouble * GetCdPr();
    RDouble * GetCdSf();
    RDouble * GetCdCl2Pa();
    RDouble * GetXcp();
    RDouble * GetSide();
    RDouble * GetHingeMoment();
};

}