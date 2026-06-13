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
//! @file      DataStruct_KDTree.h
//! @brief     DataStruct_KDTree defines a k-dimension Tree data structure.
//!            The tree starts from root node, and has a coordinates range limit.
//!            The node's definition refers to DataStruct_KDNode.
//! @author    He Kun.

#pragma once
#include "DataStruct_KDNode.h"
#include "Precision.h"
#include "string.h"

namespace PHSPACE
{

class DataStruct_KDHyperRect
{
public:
    int dim;
    RDouble *min, *max;
public:
    DataStruct_KDHyperRect()
    {
        dim = 0;
        min = NULL;
        max = NULL;
    };

    DataStruct_KDHyperRect(int dim)
    {
        this->dim = dim;
        min = new RDouble[dim];
        max = new RDouble[dim];
        for (int i = 0; i < dim; ++ i)
        { 
            min[i] = 0.0;
            max[i] = 1.0;
        }
    };

    ~DataStruct_KDHyperRect()
    {
        if (this->min)
        {
            delete [] this->min;
        }
        if (this->max)
        {
            delete [] this->max;
        }
    };

    void Free()
    {
        if (this->min)
        {
            delete [] this->min;
        }
        if (this->max)
        {
            delete [] this->max;
        }
    };

    void Extend(const RDouble *pos);

    RDouble distSQ(const RDouble *pos);
};

class DataStruct_ResNode
{
public:
    typedef DataStruct_KDNode KDNode;
    DataStruct_ResNode()
    {
        item = NULL;
        distSQ = 0.0;
        next = NULL;
    }
public:
    KDNode *item;
    RDouble distSQ;
    DataStruct_ResNode *next;
};

typedef DataStruct_KDNode      KDNode;
typedef DataStruct_KDHyperRect KDHyperRect;
typedef DataStruct_ResNode     KDResNode;

int ListInsert(KDResNode *list, KDNode *item, RDouble distSQ);

class DataStruct_KDTree
{
public:
    int dim;
    KDNode *root;
    KDHyperRect *rect;
    void (*destr)(int);

public:
    //! Construct an K-Dimension tree.
    //! @param[in] ndim the dimension of element's coordinate in node, default is 3.
    //! The element's numerical range is default with 0 to 1.
    DataStruct_KDTree()
    {
        dim  = 0;
        root = NULL;
        rect = NULL;
        destr = NULL;
    };

    DataStruct_KDTree(int dim)
    {
        this->dim = dim;
        root = NULL;
        rect = new KDHyperRect(dim);
        destr = NULL;
    };

    ~DataStruct_KDTree()
    {
        if (this)
        {
            this->Clear();
        }
    };

    void Free();

    void Clear();

    void ClearRec(KDNode *node, void(*destr)(int));

    void DestructorKDData(void(*destr)(int));

    void ExtendHyperRect(int dim, const RDouble *pos);

    KDResNode * KDNearest(DataStruct_KDTree *kdtree, const RDouble *pos);
};

typedef DataStruct_KDTree KDTree;

class DataStruct_KDRes
{
public:
    KDTree *kdTree;
    KDResNode *rlist, *riter;
    int size;
    DataStruct_KDRes()
    {
        size   = 0;
        kdTree = NULL;
        rlist  = NULL;
        riter  = NULL;
    };
    ~DataStruct_KDRes()
    {
        if (kdTree) delete kdTree;
        if (rlist)  delete rlist;
        //if (riter)  delete riter;
    }


public:
    int Size()
    {
        return size;
    };

    void Rewind()
    {
        riter = rlist->next;
    };

    int End()
    {
        return riter == NULL;
    };

    int Next()
    {
        riter = riter->next;
        return riter != NULL;
    };

    int Item(RDouble *pos)
    {
        if (riter)
        {
            if (pos)
            {
                memcpy(pos, riter->item->pos, kdTree->dim * sizeof(RDouble));
            }
            return riter->item->data;
        }
        return 0;
    };

    RDouble * Position()
    {
        if (riter)
        {
            if (riter->item->pos)
            {
                return riter->item->pos;
            }
        }
        return NULL;
    };

    int itemData()
    {
        return Item(0);
    };
};
typedef DataStruct_KDRes KDRes;

KDTree * CreatKDTree(int k);

void ClearRec(KDNode *node, void(*destr)(void *));

KDRes * NearestRange(KDTree *kdTree, const RDouble *pos, RDouble range);

KDRes * KDNearest(KDTree *kdTree, const RDouble *pos);

void KDNearestI(KDNode *node, const RDouble *pos, const RDouble **result, RDouble *resultDistSQ, KDHyperRect *rect);

KDHyperRect * DuplicateHyperRect(const KDHyperRect *rect);

int KDInsert(KDTree *kdTree, const RDouble *pos, int data);

int InsertRList(KDResNode *list, KDNode *item, RDouble distsq);

void FreeKDRes(KDRes *rset);

int FindNearest(KDNode *node, const RDouble *pos, RDouble range, KDResNode *list, int ordered, int dim);

int InsertRec(KDNode **nptr, const RDouble *pos, int data, int dir, int dim);

KDHyperRect * CreatHyperRect(int dim, const RDouble *min, const RDouble *max);

void clearResults(KDRes *rset);

}