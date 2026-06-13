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
//! @file      GeometryUnit.h 
//! @brief     It defines the hybrid grid structures, for arbitrary mixed cells.
//! @author    Bell, Baka.

#pragma once

namespace PHSPACE
{
class DYNode
{
public:
    //! Coor before deform.
    double x, y, z;

    //! Coor after deform.
    double xnew, ynew, znew;

    int moveType;

    //! Number of node around per node.
    int nodeNumberAround;

    //! Node to node.
    int *node2node;

    //! Spring 'k' of each node.
    double *node2nodeK;

public:
    DYNode();
    ~DYNode();
};


class DYFace
{
public:
    //! Face boundary type.
    int bcType;

    //! Left and right cell index.
    int leftCell, rightCell;

    //! Number of node in per face.
    int nodeNumber;

    //! Node index.
    int *face2node;

    double area;
    double xfc, yfc, zfc;
    double xfn, yfn, zfn;
public:
    DYFace();
    ~DYFace();
};


class DYCell
{
public:
    //! Number of face in per cell.
    int faceNumber;

    //! Face index.
    int *cell2face;

    //! Number of node in per cell.
    int nodeNumber;

    //! Node index.
    int *cell2node;

    //! Number of cell around per cell.
    int cellNumberAround;

    //! Cell index.
    int *cell2cell;

    double volume;
    double xCenter, yCenter, zCenter;
    int celltype, material;
public:
    DYCell();
    ~DYCell();
};

}