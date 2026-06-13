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
//! @file      GPUFaceColor.h
//! @brief     Compute face color.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include <iostream>
#include "Geo_Grid.h"
#include "Geo_UnstructGrid.h"
#include "Region.h"
#include "Zone.h"
#include "cudaErrorHandle.h"

using namespace PHSPACE;
namespace GPUFaceColor
{
    extern int *BoundFaceConflict;
    extern int *BoundFaceConflictPosi;
    extern int *BoundFaceConflictNum;
    extern int *InteriorFaceConflict;
    extern int *InteriorFaceConflictPosi;
    extern int *InteriorFaceConflictNum;
    //! The order of face label in different colors
    extern int *BoundFaceGroup;
    //! The offset of the first face label of one face group in BoundFaceGroup
    extern int *BoundFaceGroupPosi;
    //! the number of every face group (the number of faces in one color)
    extern int *BoundFaceGroupNum;
    //! the total number of color (the total number of face groups)
    extern int BoundFaceColorNum;
    //! the device variable of BoundFaceGroup
    extern int *d_BoundFaceGroup;
    extern int *InteriorFaceGroup;
    extern int *InteriorFaceGroupPosi;
    extern int *InteriorFaceGroupNum;
    extern int  InteriorFaceColorNum;
    extern int *d_InteriorFaceGroup;

    void GPUFaceColorMain();
    void CreateFaceConflict(Grid *grid);
    void ColorFaces(Grid *grid);
    void FaceGroupHostToDevice(Grid *grid);
    void FreeFaceColorVars();

} //! namespace GPUFaceColor
