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
//! @file      DataStruct_KDNode.h
//! @brief     DataStruct_KDNode defines a "Node" data structure used by KDTree
//!            The "node" has its coordinates(x, y, z) range, level in the tree, and pointer to its "son" node.
//!            Data is stored in node's item.
//! @author    He Kun.
#include "Precision.h"

#pragma once
using namespace std;

namespace PHSPACE 
{

class DataStruct_KDNode 
{
public:
    typedef DataStruct_KDNode  KDNode;
public:

    //! The coordinate of the node.
    RDouble *pos;

    //! The split Axis.
    int dir;

    //! The left tree.
    KDNode *left;

    //! The right tree.
    KDNode *right;

    int data;

public:
    //! Constructor, default dimension is three.
    DataStruct_KDNode(int dir)
    {
        this->dir = dir;
        pos = new RDouble[dir];
        left  = nullptr;
        right = nullptr;
        data = 0;
    };

    //! Constructor, set with coordinates and data.
    DataStruct_KDNode(RDouble *pos, int dir)
    {
        this->dir = dir;
        this->pos = new RDouble[this->dir];
        for (int i = 0; i < dir; ++ i)
        {
            this->pos[i] = pos[i];
        }

        left  = nullptr;
        right = nullptr;
        data = 0;
    };

    //! Constructor, set with another node.
    DataStruct_KDNode(const KDNode &node)
    {
        this->dir = node.dir;
        this->pos = new RDouble[this->dir];
        for (int i = 0; i < this->dir; ++ i)
        {
            this->pos[i] = node.pos[i];
        }

        this->left  = node.left;
        this->right = node.right;
        this->data  = node.data;
    };

    void operator = (const KDNode &node)
    {
        this->dir = node.dir;
        this->pos = new RDouble[this->dir];
        for (int i = 0; i < this->dir; ++ i)
        {
            this->pos[i] = node.pos[i];
        }

        this->left  = node.left;
        this->right = node.right;
        this->data  = node.data;
    };

    //! Destructor.
    ~DataStruct_KDNode()
    {
        if (pos) 
        {
            delete [] pos;
            pos = nullptr;
        }
        if (left)
        {
            delete left;
            left = nullptr;
        }
        if (right)
        {
            delete right;
            right = nullptr;
        }
    };

    //! Add an Adt node under the current node.

    //! Count all the descendant nodes, including self.
    int GetNodeNum()
    {
        int count = 0;
        count += 1;
        if (this->left)
        {
            count += left->GetNodeNum();
        }
        if (this->right)
        {
            count += right->GetNodeNum();
        }
        return count;
    };

    //! Get the date stored in item.
    RDouble * GetData() const
    {
        return pos;
    };
};

}