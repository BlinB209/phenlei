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
//! @file      LinkStruct.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "LIB_Macro.h"
#include "DataStruct_Sort.h"
#include "Geo_Point.h"

namespace PHSPACE
{

//! LinkStruct only works when to generate and convert grids, and it will not be used the rest of time.
class LinkStruct
{
public:
    typedef DataStruct_AdtTree<int,RDouble> AdtTree;
private:
    //! coor_tree is created to find the point column of interfaces.
    AdtTree coor_tree;
    //! facelist is mainly used as a search list of global face.
    set < DataStruct_Sort< VInt > > facelist;
    //! In general, each interface consists of two faces of different blocks.
    //! This requires the block index number of each face and the faces index number in current block.
    //! That is what zoneid and faceid do, and consistent with the order of global faces.
    VVInt zoneid;
    VVInt faceid;
    VVInt facemap;
    RDouble tol;
    vector< RDouble > x, y, z;
    set< Geo_Point< RDouble > > point_set;
public:
    AdtTree & GetCoordinateTree() { return coor_tree; }
    set < DataStruct_Sort< VInt > > & GetFaceList() { return facelist; }
    VVInt & GetZoneID() { return zoneid; }
    VVInt & GetFaceID() { return faceid; }
    VVInt & GetFaceMap() { return facemap; }
    RDouble GetTolerance() const { return tol; }
    vector< RDouble > & GetX() { return x; }
    vector< RDouble > & GetY() { return y; }
    vector< RDouble > & GetZ() { return z; }
    set< Geo_Point< RDouble >  > & GetPointSet() { return point_set; }
public:
    LinkStruct(RDouble *pmin, RDouble *pmax, RDouble tol, int nblocks) : coor_tree(3, pmin, pmax)
    {
        this->tol = tol;
        facemap.resize(nblocks);
    };

    ~LinkStruct()
    {
    };
};

}
