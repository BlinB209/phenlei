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
//! @file      Geo_FaceTopo_Unstruct.h
//! @brief     It defines the face topology of the unstructured grid,
//!            such as the face-cell and face-node connected relationship.
//! @author    Bell (Modified).

#pragma once
#include "LIB_Macro.h"
namespace PHSPACE
{
//! @brief Geo_FaceTopo_Unstruct class defines the face-node connected relationship of unstructured grid.\n
//! 1: right and left cell of each face.\n
//! 2: nodes list of each face.
class Geo_FaceTopo_Unstruct
{
private:
    //! Left cell index of each face.\n
    //! Left Cell: the cell that the INVERSE face normal point to.
    int *left_cell_of_face;

    //! Right cell index of each face.\n
    //! Right Cell: the cell that the face normal point to.
    int *right_cell_of_face;

    //! Number of nodes per face.  
    int *node_number_of_each_face;

    //! Face to node: node index of each face.
    int *face2node;  

    //! The face2node is a 1D array which record every face's node index.
    //! face2nodeSubscript records each face's starting subcript of the face2node array.
    long long int *face2nodeSubscript;

    //! Face to node: node index of each face.
    int **face2nodeArray;

public:
    LIB_EXPORT Geo_FaceTopo_Unstruct();
    LIB_EXPORT ~Geo_FaceTopo_Unstruct();

public:
    //! Set the left cell list of each face.
    void SetLeftCellOfFace(int *left_cell_of_face_in);

    //! Set the right cell list of each face.
    void SetRightCellOfFace(int *right_cell_of_face_in); 

    //! Set the number of nodes per face.
    void SetNodeNumberOfEachFace(int *node_number_of_each_face);

    //! Set the face to node, node index of each face.
    void SetFace2Node(int *face2node);

    //! Set each face's starting subscript in the face2node array.
    void SetFace2NodeSubscript(long long int *face2nodeSubscript);

    //!
    void SetFace2NodeArray(int **face2NodeArrayIn);

    //! Get the left cell list of each face.
    int * GetLeftCellOfFace() const;

    //! Get the right cell list of each face.
    int * GetRightCellOfFace() const;

    //! Get the number of nodes per face.
    int * GetNodeNumberOfEachFace() const;

    //! Get the face to node, node index of each face.
    int * GetFace2Node() const;

    //! Get each face's starting subscript in the face2node array.
    long long int * GetFace2NodeSubscript() const;

    //! Get the face to node, node index of each face.
    int ** GetFace2NodeArray() const;
};

#include "Geo_FaceTopo_Unstruct.hxx"
}