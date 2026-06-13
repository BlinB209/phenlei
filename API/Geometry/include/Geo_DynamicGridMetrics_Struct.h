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
//! @file      Geo_DynamicGridMetrics_Struct.h
//! @brief     It defines the dynamic grid metrics of the structured grid,
//!            such as the face normal velocity.
//! @author    Bell, Chang Xinghua, Zhang Jian, He Xin.

#pragma once
#include "TypeDefine.h"
namespace PHSPACE
{
//! @brief Geo_DynamicGridMetrics_Struct class defines the dynamic grid metrics of the structured grid.
//! 1: Cell volume of the last time step.
//! 2: Face normal velocity.
//! 3: Face velocity.
class Geo_DynamicGridMetrics_Struct
{
private:
    //! Cell volume of the last time step.
    RDouble3D *voln;

    //! Face normal velocity.
    RDouble4D *vgn;

    //! Face velocity.
    RDouble4D *xfv, *yfv, *zfv;

    //! cell velocity.
    RDouble3D *xcv, *ycv, *zcv;

public:
    Geo_DynamicGridMetrics_Struct();
    ~Geo_DynamicGridMetrics_Struct();

public:
    //! Return the cell volume of the last time step (voln).
    RDouble3D * GetCellVolumeOld() const;

    //! Return the Face normal velocity (vgn).
    RDouble4D * GetFaceNormalVelocity() const;

    //! Return the x component of face velocity (xfv).
    RDouble4D * GetFaceVelocityX() const;

    //! Return the y component of face velocity (yfv).
    RDouble4D * GetFaceVelocityY() const;

    //! Return the z component of face velocity (zfv).
    RDouble4D * GetFaceVelocityZ() const;

    //! Return the X cell velocity pointer.
    RDouble3D * GetCellVelocityX() const;

    //! Return the Y cell velocity pointer.
    RDouble3D * GetCellVelocityY() const;

    //! Return the Z cell velocity pointer.
    RDouble3D * GetCellVelocityZ() const;

    //! Assign the cell volume of the last time step (voln).
    void SetCellVolumeOld(RDouble3D *voln);

    //! Assign the Face normal velocity (vgn).
    void SetFaceNormalVelocity(RDouble4D *vgn);

    //! Assign the x component of face velocity (xfv).
    void SetFaceVelocityX(RDouble4D *xfv);

    //! Assign the y component of face velocity (yfv).
    void SetFaceVelocityY(RDouble4D *yfv);

    //! Assign the z component of face velocity (zfv).
    void SetFaceVelocityZ(RDouble4D *zfv);

    //! Assign the given pointer to the X cell velocity.
    void SetCellVelocityX(RDouble3D *xcv);

    //! Assign the given pointer to the Y cell velocity.
    void SetCellVelocityY(RDouble3D *ycv);

    //! Assign the given pointer to the Z cell velocity.
    void SetCellVelocityZ(RDouble3D *zcv);

private:
    void DeallocateAll();

};

#include "Geo_DynamicGridMetrics_Struct.hxx"
}