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
//! @file      Geo_CellMetrics_Struct.h
//! @brief     It defines the cell metrics of the structured grid,
//!            such as the cell center and cell volume.
//! @author    Bell, He Xin.
#pragma once
#include "TypeDefine.h"

namespace PHSPACE
{
//! @brief Geo_CellMetrics_Struct class defines the cell metrics of structured grid.
//! 1: cell center.
//! 2: cell volume.
class Geo_CellMetrics_Struct
{
private:
    //! Cell center.
    RDouble3D *xcc, *ycc, *zcc;

    //! Cell volume.
    RDouble3D *vol;

    //! Cell length.
    RDouble3D *xlen, *ylen, *zlen;

    RDouble3D *jacobian;

    //! Used for DES simulation ONLY: largest local grid length, distance of neighbor cell centers.
    RDouble3D *largestLocalGridLength;

    //! Used for DES simulation ONLY: smallest local grid length, distance of neighbor cell centers.
    RDouble3D *smallestLocalGridLength;

    //! Used for DES simulation ONLY: sub-grid length.
    RDouble3D *subgridLength;

public:
    Geo_CellMetrics_Struct();
    ~Geo_CellMetrics_Struct();

public:
    //! Return the X component of cell center (xcc).
    RDouble3D * GetCellCenterX() const;

    //! Return the Y component of cell center (ycc).
    RDouble3D * GetCellCenterY() const;

    //! Return the Z component of cell center (zcc).
    RDouble3D * GetCellCenterZ() const;

    //! Return the cell volume (vol).
    RDouble3D * GetCellVolume() const;

    //! Return the X direction cell length.
    RDouble3D * GetCellLengthX() const;

    //! Return the Y direction cell length.
    RDouble3D * GetCellLengthY() const;

    //! Return the Z direction cell length.
    RDouble3D * GetCellLengthZ() const;

    //! 
    RDouble3D * GetCellJacobian() const;

    //! Return the largest local grid spacing.
    RDouble3D * GetLargestLocalGridLength();

    //! Return the smallest local grid spacing.
    RDouble3D * GetSmallestLocalGridLength();

    //! Return the sub-grid length, used for DES simulation ONLY.
    RDouble3D * GetSubgridLength();

    //! Assign the largest local grid spacing.
    void SetLargestLocalGridLength(RDouble3D *length);

    //! Assign the smallest local grid spacing.
    void SetSmallestLocalGridLength(RDouble3D *length);

    //! Assign the sub-grid length, used for DES simulation ONLY.
    void SetSubgridLength(RDouble3D *length);

    //! Init all the metrics pointers (allocate memory).
    //! @param[in] I
    //! @param[in] J
    //! @param[in] K
    //! @note construct the three dimensional array with Range I, J, K.
    void AllocateMetrics(Range &I, Range &J, Range &K);

private:
    //! Deallocate memory of all the metrics array pointers.
    void DeallocateAll();
};

#include "Geo_CellMetrics_Struct.hxx"
}