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
//! @file      DataStruct_AdtTree.h
//! @brief     DataStruct_AdtTree defines a ADT Tree data structure.
//!            The tree starts from root node, and has a coordinates range limit.
//!            The node's definition refers to DataStruct_AdtNode.
//! @author    Bell, He Xin, reference from Dr. Wang.

#pragma once
#include "DataStruct_AdtNode.h"

namespace PHSPACE
{

template < typename T, typename U >
class DataStruct_AdtTree
{
public:
    typedef typename DataStruct_AdtNode<T,U>::AdtNode         AdtNode;
    typedef typename DataStruct_AdtNode<T,U>::AdtNodeList     AdtNodeList;
    typedef typename DataStruct_AdtNode<T,U>::AdtNodeListIter AdtNodeListIter;
    typedef          DataStruct_AdtTree<T,U>                  AdtTree;

private:
    int ndim;
    U *pmin, *pmax;
    AdtNode *root;

public:
    //! Construct an Adt tree.
    //! @param[in] ndim the dimension of element's coordinate in node, default is 3.
    //! The element's numerical range is default with 0 to 1.
    DataStruct_AdtTree(int ndim = 3);

    //! Construct an Adt tree with the region's elements' coordinate range.
    //! pmin and pmax is a pointer to a ndim array which stores each dimension's elements' coordinate range.
    //! @param[in] ndim the dimension of elements coordinate in node
    //! @param[in] pmin the minimum value of all elements' coordinate
    //! @param[in] pmax the maximum value of all elements' coordinate
    DataStruct_AdtTree(int ndim, U *pmin, U *pmax);

    //! Destructor.
    ~DataStruct_AdtTree();

    //! Add an Adt node to the AdtTree.
    void AddNode(AdtNode *node);

    //! Find All nodes inside the region (pmin,pmax) from the tree.
    void FindNodesInRegion(U *pmin, U *pmax, AdtNodeList &ld, const size_t &sizeLimit = 0);

    //! Get the total number of all nodes.
    int GetNodeNum();

    //! Get the min coordinates of the tree.
    U * GetMin() const;

    //! Get the max coordinates of the tree.
    U * GetMax() const;

    //! Get the root node of this tree.
    AdtNode * GetRoot() const;

    //! Traverse all the node from root.
    void LevelTraverse() const;

};
#include "DataStruct_AdtTree.hxx"

}