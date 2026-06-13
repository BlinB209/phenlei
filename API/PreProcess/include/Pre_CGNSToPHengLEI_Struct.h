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
//! @file      Pre_CGNSToPHengLEI_Struct.h
//! @brief     Grid conversion from CGNS to PHengLEI.
//! @author    Xu Gang.

#pragma once
#include "LIB_Macro.h"
#include "Pre_CGNSConversion_Struct.h"
#include "cgnslib.h"

using namespace std;

namespace PHSPACE
{
class CGNS_Str_Data;

class Pre_CGNSToPHengLEI_Struct : public Pre_CGNSConversion_Struct
{
public:
    LIB_EXPORT Pre_CGNSToPHengLEI_Struct(const string &gridFileName);
    LIB_EXPORT ~Pre_CGNSToPHengLEI_Struct();

private:
    void Conversion();

private:
    void SetCoordinate();
    void SetVolumeInfo();
    void SetBCInfo();
    void CheckMeshMultigrid();

};

}