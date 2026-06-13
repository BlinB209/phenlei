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
//! @file      FieldProxy.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "TypeDefine.h"

namespace PHSPACE
{
class FieldProxy
{
private:
    RDouble ***field_uho;
    RDouble  **field_uns;
    RDouble4D *field_str;
    bool del_field;
public:
    FieldProxy();
    ~FieldProxy();
public:
    RDouble *** GetField_UHO();
    RDouble  ** GetField_UNS();
    RDouble4D & GetField_STR();

    void SetField_UHO(RDouble ***field_uho, bool del_field = false);
    void SetField_UNS(RDouble  **field_uns, bool del_field = false);
    void SetField_STR(RDouble4D *field_str, bool del_field = false);
};

}
