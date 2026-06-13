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
//! @file      GeomProxy.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "TypeDefine.h"

namespace PHSPACE
{

class GeomProxy
{
public:
    // type definitions
    typedef RDouble value_type;
    typedef RDouble *iterator;
    typedef const RDouble *const_iterator;
    typedef RDouble &reference;
private:
    iterator xfn, yfn, zfn, vgn, area;
public:
    GeomProxy();
    ~GeomProxy();
public:
    void Create(int nsize);

    //! Return the X of face normal.
    iterator GetFaceNormalX();

    //! Return the Y of face normal.
    iterator GetFaceNormalY();

    //! Return the Z of face normal.
    iterator GetFaceNormalZ();

    //! Return the face area.
    iterator GetFaceArea();

    //! Return the face velocity.
    iterator GetFaceVelocity(); 
};

}
