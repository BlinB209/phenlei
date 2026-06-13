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
//! @file      Geo_Element.h
//! @brief     Explain this file briefly.
//! @author    xxx

#pragma once
#include "Math_BasisFunction.h"
#include "cgnslib.h"
using namespace std;

namespace PHSPACE
{
class BasicElement
{
public:
    typedef vector< int > int_array;
public:
    BasicElement();
    ~BasicElement();
private:
    int elem_type;
    vector< int_array > mp_struct;     //! Midpoint data structure.
    vector< int_array > child_element_index;
    vector< int_array > child_element_index_case2;
    vector< int > child_element_type;
    vector< int > face_element_type;
    vector< int_array > face_list;

    vector< int >       composite_face_element_type;
    vector< int_array > composite_face_list;
public:
    void Init(int elem_type);
    //!
    int_array & GetElementPhysicsFace(int iFace) { return composite_face_list[iFace]; }
    //!
    int GetElementPhysicsFaceType(int iFace) const { return composite_face_element_type[iFace]; }
    //! Get the face number of element.
    uint_t GetElementFaceNumber() const { return face_list.size(); }
    //!
    uint_t GetElementPhysicsFaceNumber() const { return composite_face_list.size(); }
    //! Get the element type.
    int GetElementType() const { return elem_type; }
    //! Get the node number of element.
    int GetElementNodeNumber(int elem_type);
    //! Get the refine element type.
    int GetRefineElemType(int elem_type, const bool &anisotropicAdapt);
    //! Get the simple element type.
    int GetSimpleElemType(int elem_type);
    //! Whether the type is basic element or not.
    bool ISBasicElement(int elem_type);
    //! Whether the type is face element or not.
    bool IsFaceElementType(int elementType);
    //!
    bool IsFaceElementTypeAtLeast(int elementType);
    //! Whether the type is volume element or not.
    bool IsBasicVolumeElementType(int elementType);
    //! Get the faces of element.
    int_array & GetElementFace(int iFace) { return face_list[iFace]; }
    //! Get the type of the element faces.
    int GetFaceElementType(int iFace) const { return face_element_type[iFace]; }

    int_array & GetMPlist(int i) { return mp_struct[i]; }
    //! Get the index of child element.
    int_array & GetChildIndex(int i) { return child_element_index[i]; }
    //! Get the number of child element.
    uint_t GetChildElementNumber() const { return child_element_type.size(); }
    //! Get the type of child element.
    int GetChildElementType(int i) { return child_element_type[i]; }
    //! Get the relative node index of child element.
    vector< int > & GetChildElementRelativeNodeIndex(int iChildElement) { return child_element_index[iChildElement]; }
    vector< int > & GetChildElementRelativeNodeIndexCase2(int iChildElement) { return child_element_index_case2[iChildElement]; }
};

class ElementProxy
{
public:
    ElementProxy();
    ~ElementProxy();
private: 
    BasicElement **basic_element;
    int basic_element_number;
public:
    BasicElement * GetBasicElement(int elem_type) { return basic_element[elem_type]; }
};

class IndexArray
{
public:
    IndexArray();
    ~IndexArray();
private:
    vector<int> index;
public:
    uint_t size() { return index.size(); }
    void resize(int n) { index.resize(n); }

    int& operator[](int i)
    {
        return index[i];
    }

    const int& operator[](int i) const
    {
        return index[i];
    }
};

//! The following sections is copy from HyperFLOW 3.0.
//! Copy start.
class ElementInitializingClass
{
public:
    ElementInitializingClass();
    ~ElementInitializingClass();
public:
    void InitConstants();
};

const int INVALID_INDEX = -1;
const int PENTA_12      = NofValidElementTypes;
const int QUAD_6        = NofValidElementTypes+1;
//! Get the refine element type.
int GetRefineElemType(int elementType, const bool &anisotropicAdapt);
//! Get the simple element type.
int GetSimpleElementType(int elementType);
//! Whether the type is basic element or not.
bool ISBasicElement(int elementType);
//! Whether the type is face element or not.
bool IsFaceElementType(int elementType);
//!
bool IsFaceElementTypeAtLeast(int elementType);
//! Whether the type is volume element or not.
bool IsBasicVolumeElementType(int elementType);
//! Get the basic element type.
BasicElement * GetBasicElement(int elementType);
//! Get the number of child element.
uint_t GetChildElementNumbers(int elementType);
//! Get the number of child element.
uint_t GetChildElementNumber(int elementType);
//! Get the node number of element.
int GetElementNodeNumbers(int elementType);
//!
vector<int> & GetRelatedPointListForMiddlePointComputation(int elem_type, int i);
//!//! If the right cell is boundary face.
bool IsBoundaryFace(const int &rightCellIndex);
//! If the right cell is not exist.
bool IsElementNotExist(const int &rightCellIndex);
//! Copy end.
//! Get the number of basic element.
int GetNumberOfBasicElement();

typedef enum
{

    TYPE_EDGE,            /* Edge */
    TYPE_FACE_TRIA,       /* Triangle */
    TYPE_FACE_QUAD,       /* Quadrangle */
    TYPE_FACE_POLY,       /* Simple Polygon */
    TYPE_CELL_TETRA,      /* Tetrahedron */
    TYPE_CELL_PYRAM,      /* Pyramid */
    TYPE_CELL_PRISM,      /* Prism (pentahedron) */
    TYPE_CELL_HEXA,       /* Hexahedron (brick) */
    TYPE_CELL_POLY,       /* Simple Polyhedron (convex or quasi-convex) */
    TYPE_N_ELEMENT_TYPES, /* Number of element types */
    TYPE_NOT_CONSIDERED
} ElementType;
}
