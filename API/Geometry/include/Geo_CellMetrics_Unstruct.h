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
//! @file      Geo_CellMetrics_Unstruct.h
//! @brief     It defines the cell metrics of the unstructured grid,
//!            such as the cell center and cell volume.
//! @author    Bell, He Xin.

#pragma once
#include "LIB_Macro.h"
#include "Precision.h"
namespace PHSPACE
{
//! @brief Geo_CellMetrics_Unstruct class defines the cell metrics of unstructured grid.\n
//! 1: cell center.\n
//! 2: cell volume.
class Geo_CellMetrics_Unstruct
{
private:
    //! Cell center data, including ghost cells.
    RDouble *xcc, *ycc, *zcc;

    //! Cell volume data, excluding ghosts.
    RDouble *vol; 

    //! Cell skewness data, excluding ghosts.
    RDouble *cellSkewness;

    //! Used for DES simulation ONLY: largest local grid length, distance of neighbor cell centers.
    RDouble *largestLocalGridLength;

    //! Used for DES simulation ONLY: sub-grid length.
    RDouble *subgridLength;
public:
    LIB_EXPORT Geo_CellMetrics_Unstruct();
    LIB_EXPORT ~Geo_CellMetrics_Unstruct();

public:
    //! Assign the given pointer to the X Cell center data.
    void SetCellCenterX(RDouble *xcc);

    //! Assign the given pointer to the Y Cell center data.
    void SetCellCenterY(RDouble *ycc);

    //! Assign the given pointer to the Z Cell center data.
    void SetCellCenterZ(RDouble *zcc);

    //! Assign the given pointer to the Cell volume data.
    void SetCellVolume(RDouble *vol);

    //! Assign the given pointer to the Cell skewness data.
    void SetCellSkewness(RDouble *skewness);

    //! Assign the largest local grid spacing.
    void SetLargestLocalGridLength(RDouble *length);

    //! Assign the sub-grid length, used for DES simulation ONLY.
    void SetSubgridLength(RDouble *length);

    //! Return the X pointer of cell center.
    RDouble * GetCellCenterX() const;

    //! Return the X pointer of cell center.
    RDouble * GetCellCenterY() const;

    //! Return the X pointer of cell center.
    RDouble * GetCellCenterZ() const;

    //! Return the cell volume pointer.
    RDouble * GetCellVolume()  const;

    //! Return the cell skewness pointer.
    RDouble * GetCellSkewness()  const;

    //! Return the largest local grid spacing.
    RDouble * GetLargestLocalGridLength();

    //! Return the sub-grid length, used for DES simulation ONLY.
    RDouble * GetSubgridLength();
};

#include "Geo_CellMetrics_Unstruct.hxx"

}