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
//! @file      DataStruct_AdtNode.h
//! @brief     DataStruct_AdtNode defines a "Node" data structure used by ADT Tree.
//!            The "node" has its coordinates(x, y, z), level in the tree, and pointer to its "son" node.
//!            Data is stored in node's item.
//! @author    Bell, He Xin, reference from Dr. Wang.

#pragma once
#include <vector>
#include <queue>
#include <string.h>
#include "Precision.h"
using namespace std;

namespace PHSPACE 
{

template < typename T, typename U >
class DataStruct_AdtNode 
{
public:
    typedef DataStruct_AdtNode<T,U> AdtNode;
    typedef vector<AdtNode *> AdtNodeList;
    typedef typename AdtNodeList::iterator AdtNodeListIter;

public:
    //! The coordinate of the node.
    U *point;

    //! The level in the tree.
    int level; 

    //! The left tree.
    AdtNode *left;

    //! The right tree.
    AdtNode *right;

    //! Any data stored.
    T item;

    //! The dimension of point's coordinates.
    int ndim;

public:
    //! Constructor, default dimension is three.
    DataStruct_AdtNode(int ndim = 3);

    //! Constructor, set with coordinates and data.
    DataStruct_AdtNode(int ndim, U *coor, T data);

    //! Destructor.
    ~DataStruct_AdtNode();

    //! Add an Adt node under the current node .
    void AddNode(AdtNode *node, U *nwmin, U *nwmax, const int &ndim);

    //! Judge if its the current node inside region (pmin,pmax).
    bool IsInRegion(U *pmin, U *pmax, const int &ndim);

    //! Find all the nodes inside region (pmin,pmax) among all descendant nodes, (including node self).
    void FindNodesInRegion(U *pmin, U *pmax, U *nwmin, U *nwmax, const int &ndim, AdtNodeList &ld, const uint_t &sizeLimit = 0);

    //! Count all the descendant nodes, including self.
    int GetNodeNum();

    //! Get the date stored in item.
    T GetData() const { return item; }
};

#include "DataStruct_AdtNode.hxx"

}