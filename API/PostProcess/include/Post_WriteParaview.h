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
//! @file      Post_WriteParaview.h
//! @brief     Write flow field into paraview file.
//! @author    Xu Gang, Yu Yaojie(External contributor).

#pragma once
#include "Post_WriteVisualFile.h"

namespace PHSPACE
{

//! @brief Post_WriteParaview class realize paraview visual file output function.\n
class Post_WriteParaview : public Post_WriteVisualFile
{
public:
    Post_WriteParaview();
    ~Post_WriteParaview();

private:
    void Initialize();

    //! Collect data for paraview pvd file.
    void StorePVDData();
    void StorePVDData(int zoneIndex, DataContainer *cdata);

    //! Write variable file.
    void WriteFieldData();
    void WriteUnsBoundaryFieldData(fstream &flowFile, int iBoundary);
    void WriteStrBoundaryFieldData(fstream &flowFile, int iBoundary);
    void WriteUnsBlockFieldData(fstream &flowFile, int iZone);
    void WriteStrBlockFieldData(fstream &flowFile, int iZone);

    void WritePVDFile();

    //! Dump data into paraview file.
    void WriteFile();

    void ClearFieldData();

private:
    vector<DataContainer *> pvdDataList;
    string pvdFileName;
    int iter;

};

typedef enum
{
    VTK_VERTEX = 1,
    VTK_POLY_VERTEX,
    VTK_LINE,
    VTK_POLY_LINE,
    VTK_TRIANGLE,
    VTK_TRIANGLE_STRIP,
    VTK_POLYGON,
    VTK_PIXEL,
    VTK_QUAD,
    VTK_TETRA,
    VTK_VOXEL,
    VTK_HEXAHEDRON,
    VTK_WEDGE,
    VTK_PYRAMID
} VTKElementType;

}




