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
//! @file      Pre_FieldViewToPHengLEI_Unstruct.h
//! @brief     Grid conversion from fieldview format to PHengLEI.
//! @author    Bell, He Xin.

#pragma once
#include "cgnslib.h"
#include "Pre_GridConversion.h"
#include "Pre_GridBase.h"
using namespace std;

namespace PHSPACE
{
//! @brief Pre_FieldViewToPHengLEI_Unstruct class defines the method of 
//! Grid conversion from fieldview format to PHengLEI.
class Pre_FieldViewToPHengLEI_Unstruct : public Pre_GridConversion
{
    static const int FIELDVIEW     = 1;
    static const int FIELDVIEW_2_4 = 2;
    static const int FIELDVIEW_3_0 = 3;
public:
    //! @param[in] Grid file name of field view grid.
    LIB_EXPORT Pre_FieldViewToPHengLEI_Unstruct(const string &gridFileName);

    LIB_EXPORT ~Pre_FieldViewToPHengLEI_Unstruct();

public:
    //! Dump additional information.
    LIB_EXPORT void WriteAdditionalInformation(const string &out_grid_file, Grid **grids_in);

private:
    //! Read the original grid.
    void ReadGrid();

    //! Get the version of Field view grid.
    void GetFieldViewVersion(fstream &file, int &nversion);

    //! Get the block number of Field view grid.
    void GetFieldViewNBlocks(fstream &file, int nversion);

    //! Get the BC type number of Field view grid.
    void GetFieldViewNBCType(fstream &file, int nversion, int &nBCType);

    //!
    void GetFieldViewBCNameTable(fstream &file, int nversion, int nBCType, string *bctype_name_list);

    //!
    void GetFieldViewBCTypeTable(int nBCType, int nversion, string *bctype_name_list, int *bctype_list);

    //!
    void Gridgen2CGNSBCTypeTable(int nBCType, int *bctype_list);

    //! Get the total node number of Field view grid.
    void GetFieldViewNTNode(fstream &file, int nversion, int &nTotalNode);

    //! Get the boundary face number of Field view grid.
    void GetFieldViewNBFace(fstream &file, int nversion, int &nBoundFace);

    //! Read the coordinates of Field view grid.
    void ReadFieldViewCoor(fstream &file, int nTotalNode, RDouble *x, RDouble *y, RDouble *z);

    //!
    void GetFieldViewBCCell(fstream &file, int nBoundFace, int *bctype_list,
                            int &nBC_Sections, int *&elem_type, int *&istart, int *&iend, cgsize_t **&conn_list,
                            int &nBCRegions, cgsize_t **&bc_conn_list, int *&nBCElem, int *&bc_type);

    //!
    void GetFieldViewMainCell(fstream &file, int &nTotalCell, int &nSections, int *&elem_type,
                              int *&istart, int *&iend, cgsize_t **&conn_list);

    //!
    void MergeFieldViewCell(CGNSFactory *factory_cgns, int iZone, int nBoundFace, int nTotalNode, int nTotalCell,
                            int nElem_Sections, int *elem_type, int *elem_istart, int *elem_iend, cgsize_t **elem_conn_list,
                            int nBC_Sections, int *bc_elem_type, int *bc_istart, int *bc_iend, cgsize_t **bc_elem_conn_list,
                            int nBCRegions, cgsize_t **bc_conn_list, int *nBCElem, int *bc_type);

    void ElemReorderFieldView2CGNS(int etype, cgsize_t *node_swap, cgsize_t *node_results);

    void Conversion();

private:
    CGNSFactory *factory_cgns;
};

namespace FIELDVIEW_SPACE
{
    //! Cell type definition in field view.
    const int TETR = 1;    //! Tetrahedron.
    const int PYRA = 4;    //! Pyramid.
    const int PRIS = 3;    //! Prism.
    const int HEXA = 2;    //! Hexahedron.

    //! The following two types are added for convenience and have not been confirmed.
    const int TRIA = 5;    //! Surface triangle.
    const int QUAD = 6;    //! Surface quadrilateral.
}

}