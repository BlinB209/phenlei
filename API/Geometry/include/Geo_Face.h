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
//! @file      Geo_Face.h
//! @brief     It defines the face topology, with left and right cell.
//! @author    Bell.
#pragma once
using namespace std;

namespace PHSPACE
{
class Geo_Face
{
public:
    //! Left cell id.
    int le;

    //! Right cell id.
    int re;

    //! Face id.
    int faceID;

    //! Maximum face id in global, it would be a estimated appropriate value.
    int maxFaceID;

public:
    Geo_Face();
    Geo_Face(const int &le, const int &re, const int &id, int maxFaceID = 102400);
    Geo_Face(const Geo_Face &rhs);
    Geo_Face & operator = (const Geo_Face &rhs);
    ~Geo_Face();

public:
    //! Get face id.
    int GetID() const { return faceID; }
    void SetID(const int &id) { this->faceID = id; }

    int hash_key() const;

    bool operator < (const Geo_Face &rhs) const;
    friend bool operator == (const Geo_Face &face1, const Geo_Face &face2)
    {
        if (&face1 == &face2) return true;

        return (face1.le == face2.le && face1.re == face2.re);
        //return (face1.le == face2.le && face1.re == face2.re && face1.faceID == face2.faceID);
    }
};

}

#include "Geo_Face.hxx"
