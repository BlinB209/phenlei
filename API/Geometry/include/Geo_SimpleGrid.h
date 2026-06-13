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
//! @file      Geo_SimpleGrid.h
//! @brief     It is the base class of geometry grid.
//!            The inheriting order is: SimpleGrid -> Grid -> StructuredGrid/UnstructuredGrid.
//! @author    He Xin, Bell.

#pragma once
#include "LIB_Macro.h"
#include "TypeDefine.h"

namespace PHSPACE
{
//! @brief It defines the class 'SimpleGrid', which is the base class of geometry grid.
//!        The inheriting order is: SimpleGrid -> Grid -> StructuredGrid/UnstructuredGrid.
class SimpleGrid
{
protected:
    //! Number of total nodes.
    int nTotalNode;

    //! Coordinates.
    RDouble *x, *y, *z;

    //! Min-max box of the nodes coordinates.
    RDouble *pmin, *pmax;

    //! The minimum edge length.
    RDouble minEdgeLength;

    //! The maximum edge length.
    RDouble maxEdgeLength;

    //! If the min-max box has been calculated.
    bool has_cal_minmax;

    //! The material type of grid, such as water and air.
    int materialType;

public:
    SimpleGrid();
    SimpleGrid(const SimpleGrid &rhs);
    virtual ~SimpleGrid();

public:
    SimpleGrid & operator = (const SimpleGrid &rhs);

public:
    //! Assign the given pointer to the X coordinate pointer of this point.
    void SetX(RDouble *x);

    //! Assign the given pointer to the Y coordinate pointer of this point.
    void SetY(RDouble *y);

    //! Assign the given pointer to the Z coordinate pointer of this point.
    void SetZ(RDouble *z);

    //! Return its X coordinate pointer.
    RDouble * GetX() const;

    //! Return its Y coordinate pointer.
    RDouble * GetY() const;

    //! Return its Z coordinate pointer.
    RDouble * GetZ() const;

    //! Assign the given value to the number of the nodes.
    void SetNTotalNode(int nTotalNode);

    //! Return the number of the nodes.
    int GetNTotalNode() const;

    //! Return the minimum box of the grid zone.
    RDouble * GetMinBox();
    
    //! Return the maximum box of the grid zone.
    RDouble * GetMaxBox();

    //! Compute the minimum and maximum box of the grid zone.
    LIB_EXPORT void ComputeMinMaxBox();

    //! Compute the minimum and maximum box of the grid zone.
    LIB_EXPORT bool IfBoxOverset(RDouble *pmin, RDouble *pmax, int nDim, RDouble eps = 1.0e-6);

    //! Compute the minimum distance of the grid edges.
    LIB_EXPORT virtual RDouble CalMinEdgeLength();

    //! Compute the maximum distance of the grid edges.
    LIB_EXPORT virtual RDouble CalMaxEdgeLength();

private:
    //! Free memory.
    void FreeMemory();

    //! Copy object.
    void Copy(const SimpleGrid &rhs);
};

#include "Geo_SimpleGrid.hxx"

}