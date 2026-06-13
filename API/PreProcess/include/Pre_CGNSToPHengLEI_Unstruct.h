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
//! @file      Pre_CGNSToPHengLEI_Unstruct.h
//! @brief     Grid conversion from CGNS to PHengLEI.
//! @author    Bell, He Xin.

#pragma once
#include "Pre_GridConversion.h"
#include "Pre_GridBase.h"
using namespace std;

namespace PHSPACE
{
//! @brief Pre_CGNSToPHengLEI_Unstruct class defines the method of 
//! Grid conversion from CGNS format to PHengLEI.
class Pre_CGNSToPHengLEI_Unstruct : public Pre_GridConversion
{
public:
    //! @param[in] Grid file name of field view grid.
    LIB_EXPORT Pre_CGNSToPHengLEI_Unstruct(const string &gridFileName);

    LIB_EXPORT ~Pre_CGNSToPHengLEI_Unstruct();

public:
    //! Dump additional information.
    LIB_EXPORT void WriteAdditionalInformation(const string &targetGridFileName, Grid **grids_in);

private:
    //! Read the information of original grid.
    void ReadGrid();

    //! Perform the grid conversion.
    void Conversion();

private:
    CGNSFactory *factoryCGNS;
};

namespace FIELDVIEW_SPACE
{

}

}