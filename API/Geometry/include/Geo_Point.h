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
//! @file      Geo_Point.h
//! @brief     Class of grid point.
//! @author    Bell, He Xin.

#pragma once
#include "Math_BasisFunction.h"
using namespace std;

namespace PHSPACE
{
template < typename T >
class Geo_Point
{
public:
    typedef Geo_Point<T> point_type;
private:
    T x, y, z;
    int id;
    RDouble tolerance;
public:
    Geo_Point();
    Geo_Point(const T &x, const T &y, const T &z, int id = 0);
    Geo_Point(const point_type &rhs);
    Geo_Point & operator = (const point_type &rhs);
    ~Geo_Point();
public:
    T X() const { return x; }
    T Y() const { return y; }
    T Z() const { return z; }
    int ID() const { return id; }

    void SetTolerance(RDouble dataIn) { this->tolerance = dataIn; }

    void SetX(const T &x) { this->x  = x;  }
    void SetY(const T &y) { this->y  = y;  }
    void SetZ(const T &z) { this->z  = z;  }
    void SetID(int id)    { this->id = id; }

    int hash_key() const;
public:
    bool operator < (const point_type &rhs) const;
    friend bool operator == (const Geo_Point &pt1, const Geo_Point &pt2)
    {
        RDouble epsilon = 1.0e-13;

        if (&pt1 == &pt2) return true;

        RDouble dx = ABS(pt1.X() - pt2.X());
        RDouble dy = ABS(pt1.Y() - pt2.Y());
        RDouble dz = ABS(pt1.Z() - pt2.Z());

        return (dx <= epsilon && dy <= epsilon && dz <= epsilon);
    }
};

typedef Geo_Point< RDouble > Point3D;

}

#include "Geo_Point.hxx"
