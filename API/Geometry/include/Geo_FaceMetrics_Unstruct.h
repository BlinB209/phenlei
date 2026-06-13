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
//! @file      Geo_FaceMetrics_Unstruct.h
//! @brief     It defines the face metrics of the unstructured grid,
//!            such as the face center and face normal.
//! @author    Bell (Modified).

#pragma once
#include "LIB_Macro.h"
#include "Precision.h"
namespace PHSPACE
{
//! @brief Geo_FaceMetrics_Unstruct class defines the face metrics of unstructured grid.\n
//! 1: face center.\n
//! 2: face normal.\n
//! 3: face area.
class Geo_FaceMetrics_Unstruct
{
private:
    //! Face center data.
    RDouble *xfc, *yfc, *zfc;

    //! Face UNIT normal data.
    RDouble *xfn, *yfn, *zfn;

    //! Face area.
    RDouble *area;
public:
    LIB_EXPORT Geo_FaceMetrics_Unstruct();
    LIB_EXPORT ~Geo_FaceMetrics_Unstruct();

public:
    //! Assign the given pointer to the X Face center data.
    void SetFaceCenterX(RDouble *xfc);

    //! Assign the given pointer to the Y Face center data.
    void SetFaceCenterY(RDouble *yfc);

    //! Assign the given pointer to the Z Face center data.
    void SetFaceCenterZ(RDouble *zfc);

    //! Assign the given pointer to the X Face UNIT normal data.
    void SetFaceNormalX(RDouble *xfn);

    //! Assign the given pointer to the Y Face UNIT normal data.
    void SetFaceNormalY(RDouble *yfn);

    //! Assign the given pointer to the Z Face UNIT normal data.
    void SetFaceNormalZ(RDouble *zfn);

    //! Assign the given pointer to the Face area.
    void SetFaceArea(RDouble *area);

    //! Return the X pointer of face center.
    RDouble * GetFaceCenterX() const;

    //! Return the Y pointer of face center.
    RDouble * GetFaceCenterY() const;

    //! Return the Z pointer of face center.
    RDouble * GetFaceCenterZ() const;

    //! Return the X pointer of face UNIT normal.
    RDouble * GetFaceNormalX() const;

    //! Return the Y pointer of face UNIT normal.
    RDouble * GetFaceNormalY() const;

    //! Return the Z pointer of face UNIT normal.
    RDouble * GetFaceNormalZ() const;

    //! Return the face area.
    RDouble * GetFaceArea() const;
};

#include "Geo_FaceMetrics_Unstruct.hxx"
}