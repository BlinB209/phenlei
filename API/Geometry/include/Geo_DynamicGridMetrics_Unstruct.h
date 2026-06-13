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
//! @file      Geo_DynamicGridMetrics_Unstruct.h
//! @brief     It defines the dynamic grid metrics of the unstructured grid,
//!            such as the face normal velocity.
//! @author    Bell (Modified).

#pragma once
#include "LIB_Macro.h"
#include "Precision.h"
namespace PHSPACE
{
//! @brief Geo_DynamicGridMetrics_Unstruct class defines dynamic mesh metrics of unstructured grid.\n
//! 1: face normal velocity.\n
//! 2: face velocity.
class Geo_DynamicGridMetrics_Unstruct
{
private:
    //! Cell volume of the last time step.
    RDouble *voln;

    //! Face normal velocity.
    RDouble *vgn;

    //! face velocity.
    RDouble *xfv, *yfv, *zfv;

    //! cell velocity.
    RDouble *xcv, *ycv, *zcv;

public:
    LIB_EXPORT Geo_DynamicGridMetrics_Unstruct();
    LIB_EXPORT ~Geo_DynamicGridMetrics_Unstruct();

public:
    //! Assign the given pointer to the X face velocity.
    void SetFaceVelocityX(RDouble *xfv);

    //! Assign the given pointer to the Y face velocity.
    void SetFaceVelocityY(RDouble *yfv);

    //! Assign the given pointer to the Z face velocity.
    void SetFaceVelocityZ(RDouble *zfv);
    
    //! Assign the given pointer to the X cell velocity.
    void SetCellVelocityX(RDouble *xcv);

    //! Assign the given pointer to the Y cell velocity.
    void SetCellVelocityY(RDouble *ycv);

    //! Assign the given pointer to the Z cell velocity.
    void SetCellVelocityZ(RDouble *zcv);

    //! Assign the given pointer to the face normal velocity.
    void SetFaceNormalVelocity(RDouble *vgn);

    //! Assign the given pointer to the Cell volume of the last time step.
    void SetCellVolumeOld(RDouble *voln);

    //! Return the X face velocity pointer.
    RDouble * GetFaceVelocityX() const;

    //! Return the Y face velocity pointer.
    RDouble * GetFaceVelocityY() const;

    //! Return the Z face velocity pointer.
    RDouble * GetFaceVelocityZ() const;

    //! Return the X cell velocity pointer.
    RDouble * GetCellVelocityX() const;

    //! Return the Y cell velocity pointer.
    RDouble * GetCellVelocityY() const;

    //! Return the Z cell velocity pointer.
    RDouble * GetCellVelocityZ() const;

    //! Return the face normal velocity pointer.
    RDouble * GetFaceNormalVelocity() const;

    //! Return the cell volume of the last time step.
    RDouble * GetCellVolumeOld() const;
};

#include "Geo_DynamicGridMetrics_Unstruct.hxx"
}