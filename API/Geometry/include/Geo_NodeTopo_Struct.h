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
//! @file      Geo_NodeTopo_Struct.h
//! @brief     It defines the node topology of the structured grid,
//!            such as the node dimension(ni, nj, nk) and three dimensional 
//!            array of node coordinates(structx[i][j][k], structy[i][j][k], structz[i][j][k]).
//! @author    Bell, Zhang Jian, He Xin.

#pragma once
#include "TypeDefine.h"
namespace PHSPACE
{
//! @brief Geo_NodeTopo_Struct class defines the node dimension and coordinates.
class Geo_NodeTopo_Struct
{
private:
    //! The dimension of three direction.
    int ni, nj, nk;

    //! The node coordinates array of three direction.
    RDouble3D *structx, *structy, *structz;

public:
    Geo_NodeTopo_Struct();
    ~Geo_NodeTopo_Struct();

public:
    //! Set the node dimension of I direction (ni).
    void SetNI(int ni);

    //! Set the node dimension of J direction (nj).
    void SetNJ(int nj);

    //! Set the node dimension of K direction (nk).
    void SetNK(int nk);

    //! Get the node dimension of I direction (ni).
    int GetNI() const;

    //! Get the node dimension of J direction (nj).
    int GetNJ() const;

    //! Get the node dimension of K direction (nk).
    int GetNK() const;

    //! Get the all node's x coordinates array (structx).
    RDouble3D * GetStructX() const;

    //! Get the all node's y coordinates array (structy).
    RDouble3D * GetStructY() const;

    //! Get the all node's z coordinates array (structz).
    RDouble3D * GetStructZ() const;

    //! Set the three dimensional array of coordinates from one dimensional coordinates.
    //! @param[in] x    x coordinates of all nodes.
    //! @param[in] y    y coordinates of all nodes.
    //! @param[in] z    z coordinates of all nodes.
    //! @note           use x to create a three dimensional array structx, with range I(1, ni), J(1, nj), K(1, nk).
    void SetArrayLayout(RDouble *x, RDouble *y, RDouble *z);

    //! Set the three dimensional array of coordinates from one dimensional coordinates for HighOrder structured solver.
    //! @param[in] x    x coordinates of all nodes.
    //! @param[in] y    y coordinates of all nodes.
    //! @param[in] z    z coordinates of all nodes.
    void SetArrayLayoutStructHighOrder(RDouble *x, RDouble *y, RDouble *z);

private:
    void DeallocateAll();
};

#include "Geo_NodeTopo_Struct.hxx"
}