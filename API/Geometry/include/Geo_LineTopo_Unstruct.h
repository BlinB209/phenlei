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
//! @file      Geo_LineTopo_Unstruct.h
//! @brief     It defines the face topology of the unstructured grid,
//!              such as the face-cell and face-node connected relationship.
//! @author    WanYunBo.

#pragma once
#include "LIB_Macro.h"
#include "Pointer.h"
#include "PHVector3.h"
namespace PHSPACE
{
//! @brief Geo_LineTopo_Unstruct class defines the line for line implicit LUSGS.\n
class Geo_LineTopo_Unstruct
{
private:
    //! the number of lines.
    int nLine;

    //! the number of cells in all lines.
    int *nCellsInLines;

    //! the cell No. of all lines.
    int **CellOfLine;

    //! the face No. of all lines.
    int **FaceOfLine;

    //! the line No. of all cells.
    int *LineOfCell;

    //! the node No. of Line.
    int starIndex, endIndex;

    //! the node No. of all lines.
    int *Line2Node;

    //! the node No. of all lines.
    int **Face2Line;

public:
    LIB_EXPORT Geo_LineTopo_Unstruct();
    LIB_EXPORT ~Geo_LineTopo_Unstruct();

    Vector3D *startPoint;

    Vector3D *endPoint;

public:
    //! set the number of lines.
    void SetLI_nLine(int nLine);

    //! set the number of cells in all lines.
    void SetLI_nCellsInLine(int *nCellsInLines);

    //! set the cell No. of all lines.
    void SetLI_CellOfLine(int **CellOfLine);

    //! set the face No. of all lines.
    void SetLI_FaceOfLine(int **FaceOfLine);

    //! set the line No. of all cells.
    void SetLI_LineOfCell(int *LineOfCell);

    //! set the node No. of all lines.
    void SetLine2Node(int *Line2NodeIn);

    //! set the line No. of all faces.
    void SetFace2Line(int **Face2LineIn);

    //! get the number of lines.
    int GetLI_nLine() const;

    //! get the number of cells in all lines.
    int * GetLI_nCellsInLine() const;

    //! get the cell No. of all lines.
    int ** GetLI_CellOfLine() const;

    //! get the face No. of all lines.
    int ** GetLI_FaceOfLine() const;

    //! get the line No. of all cells.
    int * GetLI_LineOfCell() const;

    void SetPointIndex(const int starIndexIn, const int endIndexIn){starIndex = starIndexIn; endIndex = endIndexIn;}

    int  GetStarIndex() const { return starIndex; }

    int  GetEndIndex() const { return endIndex; }

    int  hash_key() const { if(starIndex>endIndex) return starIndex+23*endIndex; else return endIndex+23*starIndex; }

    friend bool operator==(const Geo_LineTopo_Unstruct &Line1, const Geo_LineTopo_Unstruct &Line2)
    {
        return (Line1.starIndex==Line2.starIndex&&Line1.endIndex==Line2.endIndex) || (Line1.starIndex==Line2.endIndex&&Line1.endIndex==Line2.starIndex);
}

    //! get the node No. of all lines.
    int * GetLine2Node() const { return Line2Node; }

};

#include "Geo_LineTopo_Unstruct.hxx"
}
