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
//! @file      Geo_FaceMetrics_Struct.h
//! @brief     It defines the face metrics of the structured grid,
//!            such as the face area and face normal.
//! @author    Bell, He Xin, Zhang Jian.
#pragma once
#include "TypeDefine.h"

namespace PHSPACE
{
//! @brief Geo_FaceMetrics_Struct class defines the face metrics of structured grid.
//! 1: face normal.
//! 2: face area.
//! 3: face vector.
class Geo_FaceMetrics_Struct
{
private:
    //! UNIT face normal.
    RDouble4D *xfn, *yfn, *zfn;

    //! Face area.
    RDouble4D *area;

    //! Face vector.
    RDouble4D *xFaceVector, *yFaceVector, *zFaceVector;

    RDouble5D *xyzFaceVector_FD;
    RDouble5D *xyzFaceNormal_FD;

public:
    Geo_FaceMetrics_Struct();
    ~Geo_FaceMetrics_Struct();

public:
    //! Return the X component of UNIT face normal (xfn).
    RDouble4D * GetFaceNormalX() const;

    //! Return the Y component of UNIT face normal (yfn).
    RDouble4D * GetFaceNormalY() const;

    //! Return the Z component of UNIT face normal (zfn).
    RDouble4D * GetFaceNormalZ() const;

    //! Return the face area (area).
    RDouble4D * GetFaceArea() const;

    //! Return the X component of face vector (xFaceVector).
    RDouble4D * GetFaceVectorX() const;

    //! Return the Y component of face vector (yFaceVector).
    RDouble4D * GetFaceVectorY() const;

    //! Return the Z component of face vector (zFaceVector).
    RDouble4D * GetFaceVectorZ() const;

    RDouble5D * GetFaceVector_FD() const;
    RDouble5D * GetFaceNormal_FD() const;

    //! Init all the metrics pointers (allocate memory).
    //! @param[in] I.
    //! @param[in] J.
    //! @param[in] K.
    //! @param[in] D.
    //! @note construct the four dimensional array with Range I, J, K, D.
    void AllocateMetrics(Range &I, Range &J, Range &K, Range &D);

private:
    //! Deallocate memory of all the metrics array pointers.
    void DeallocateAll();
};

#include "Geo_FaceMetrics_Struct.hxx"
}