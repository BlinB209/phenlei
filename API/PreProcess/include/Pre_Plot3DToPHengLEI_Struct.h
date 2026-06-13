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
//! @file      Pre_Plot3DToPHengLEI_Struct.h
//! @brief     Grid conversion from Plot3D to PHengLEI.
//! @author    Xu Gang.

#pragma once
#include "Pre_GridConversion.h"
#include "Geo_StructGrid.h"

using namespace std;

namespace PHSPACE
{

class Pre_Plot3DToPHengLEI_Struct : public Pre_GridConversion
{
public:
    LIB_EXPORT Pre_Plot3DToPHengLEI_Struct(const string &gridFileName);
    LIB_EXPORT ~Pre_Plot3DToPHengLEI_Struct();

private:
    void ReadGrid();
    void Conversion();

private:
    int  CheckPlot3DFileIfASCII();
    int  CheckGridgenCoorASCIIFileType();
    void ReadGridgenCoorASCIIPlot3D();
    void ReadGridgenCoorASCIIGridgen();
    void ReadGridgenCoor();
    void ReadGridgenBC();
    void CheckMeshMultigrid();

};

}