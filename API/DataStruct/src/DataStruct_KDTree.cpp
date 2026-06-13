#pragma once
#include "DataStruct_KDTree.h"
#include <cstring>
#include <cmath>
#include <new>
using namespace std;

#define SQ(x) ((x)*(x))

namespace PHSPACE
{
    KDHyperRect * CreatHyperRect(int dim, const RDouble *min, const RDouble *max)
    {
        size_t size = dim * sizeof(RDouble);
        KDHyperRect *rect = new KDHyperRect(dim);
        if (!rect)
        {
            return NULL;
        }

        /*rect->dim = dim;
        
        rect->min = new(std::nothrow) RDouble[dim];
        if (rect->min == NULL)
        {
            delete rect;
            return NULL;
        }

        rect->max = new(std::nothrow) RDouble[dim];
        if (rect->max == NULL)
        {
            delete rect->min;
            delete rect;
            return NULL;
        }*/

        memcpy(rect->min, min, size);
        memcpy(rect->max, max, size);

        return rect;
    }

    void KDHyperRect::Extend(const RDouble *pos)
    {
        for (int i = 0; i < dim; ++ i)
        {
            if (pos[i] < min[i])
            {
                min[i] = pos[i];
            }

            if (pos[i] > max[i])
            {
                max[i] = pos[i];
            }
        }
    }

    RDouble KDHyperRect::distSQ(const RDouble *pos)
    {
        RDouble result = 0.0;
        for (int i = 0; i < dim; ++ i)
        {
            if (pos[i] < min[i])
            {
                result += SQ(min[i] - pos[i]);
            }

            if (pos[i] > max[i])
            {
                result += SQ(max[i] - pos[i]);
            }
        }
        return result;
    }

    KDTree * CreatKDTree(int k)
    {
        KDTree *kdTree = new KDTree(k);
        if (!kdTree)
        {
            return NULL;
        }

        return kdTree;
    };

    void KDTree::Free()
    {
        if (this)
        {
            this->Clear();
        }
    }

    void KDTree::Clear()
    {
        ClearRec(root, destr);
        root = NULL;

        if (rect)
        {
            delete rect;
            rect = NULL;
        }
    }

    void KDTree::ClearRec(KDNode *node, void(*destr)(int))
    {
        if (!node)
        {
            return;
        }
        ClearRec(node->left, destr);
        ClearRec(node->right, destr);

        if (destr)
        {
            destr(node->data);
        }
        if (node->right)
        {
            node->right = nullptr;
        }
        if (node->left)
        {
            node->left = nullptr;
        }
        //delete [] node->pos;
        delete node;
    }

    void KDTree::DestructorKDData(void(*destr)(int))
    {
        this->destr = destr;
    }

    int InsertRec(KDNode **nptr, const RDouble *pos, int data, int dir, int dim)
    {
        int newDir;
        KDNode *node;

        if (!*nptr)
        {
            try {
                node = new KDNode(dir);
            }
            catch (bad_alloc) {
                return -1;
            }

            node->pos = new(std::nothrow) RDouble[dim];
            if ( node->pos == NULL )
            {
                delete node;
                return -1;
            }

            memcpy(node->pos, pos, dim * sizeof(node->pos));
            node->data = data;
            node->dir = dir;
            node->left = NULL;
            node->right = NULL;
            *nptr = node;
            return 0;
        }

        node = *nptr;
        newDir = (node->dir + 1) % dim;
        if (pos[node->dir] < node->pos[node->dir])
        {
            return InsertRec(&(*nptr)->left, pos, data, newDir, dim);
        }
        return InsertRec(&(*nptr)->right, pos, data, newDir, dim);
    }

    int KDInsert(KDTree *kdTree, const RDouble *pos, int data)
    {
        if (InsertRec(&kdTree->root, pos, data, 0, kdTree->dim))
        {
            return -1;
        }

        if (kdTree->rect == NULL)
        {
            kdTree->rect = CreatHyperRect(kdTree->dim, pos, pos);
        }
        else
        {
            kdTree->rect->Extend(pos);
        }

        return 0;
    }

    int ListInsert(KDResNode *list, KDNode *item, RDouble distSQ)
    {
        KDResNode *rnode = new KDResNode();
        if (!rnode)
        {
            return -1;
        }
        rnode->item = item;
        rnode->distSQ = distSQ;

        if (distSQ >= 0.0)
        {
            while (list->next && list->next->distSQ < distSQ)
            {
                list = list->next;
            }
        }
        rnode->next = list->next;
        list->next = rnode;
        return 0;
    }

    KDHyperRect * DuplicateHyperRect(const KDHyperRect *rect)
    {
        return CreatHyperRect(rect->dim, rect->min, rect->max);
    };

    int FindNearest(KDNode *node, const RDouble *pos, RDouble range, KDResNode *list, int ordered, int dim)
    {
        if (!node)
        {
            return 0;
        }

        RDouble distSQ = 0.0;
        int addedRes = 0;
        for (int i = 0; i < dim; i ++)
        {
            distSQ += SQ(node->pos[i] - pos[i]);
        }
        if (distSQ <= SQ(range))
        {
            if (ListInsert(list, node, ordered ? distSQ: -1.0) == -1)
            {
                return -1;
            }
            addedRes = 1;
        }

        RDouble dx = pos[node->dir] - node->pos[node->dir];
        int ret = FindNearest(dx <= 0.0 ? node->left: node->right, pos, range, list, ordered, dim);
        if (ret >= 0 && fabs(dx) < range)
        {
            addedRes += ret;
            ret = FindNearest(dx <= 0.0 ? node->right: node->left, pos, range, list, ordered, dim);
        }
        if (ret == -1)
        {
            return -1;
        }
        addedRes += ret;
        return addedRes;
    }

    void FreeKDRes(KDRes *rset)
    {
        clearResults(rset);
        //delete rset->rlist;
        delete rset;
    }

    void clearResults(KDRes *rset)
    {
        KDResNode *node = rset->rlist->next;
        while (node)
        {
            KDResNode *tmp = node;
            node = node->next;
            delete tmp;
        }
        rset->rlist->next = NULL;
    }

    void KDNearestI(KDNode *node, const RDouble *pos, KDNode **result, RDouble *resultDistSQ, KDHyperRect *rect)
    {
        int dir = node->dir;
        RDouble distSQ;
        KDNode *nearestSubTree, *fatherSubTree;
        RDouble *nearerHyperRectCoord, *fatherHyperRectCoord;

        RDouble dummy = pos[dir] - node->pos[dir];
        if (dummy <= 0.0)
        {
            nearestSubTree = node->left;
            fatherSubTree  = node->right;
            nearerHyperRectCoord = rect->max + dir;
            fatherHyperRectCoord = rect->min + dir;
        }
        else
        {
            nearestSubTree = node->right;
            fatherSubTree  = node->left;
            nearerHyperRectCoord = rect->min + dir;
            fatherHyperRectCoord = rect->max + dir;
        }

        if (nearestSubTree)
        {
            dummy = *nearerHyperRectCoord;
            *nearerHyperRectCoord = node->pos[dir];

            KDNearestI(nearestSubTree, pos, result, resultDistSQ, rect);

            *nearerHyperRectCoord = dummy;
        }

        distSQ = 0.0;
        for (int i = 0; i < rect->dim; ++ i)
        {
            distSQ += SQ(node->pos[i] - pos[i]);
        }
        if (distSQ < *resultDistSQ)
        {
            *result = node;
            *resultDistSQ = distSQ;
        }

        if (fatherSubTree)
        {
            dummy = *fatherHyperRectCoord;
            *fatherHyperRectCoord = node->pos[dir];
            if (rect->distSQ(pos) < *resultDistSQ)
            {
                KDNearestI(fatherSubTree, pos, result, resultDistSQ, rect);
            }
            *fatherHyperRectCoord = dummy;
        }
    }

    KDRes * KDNearest(KDTree *kdTree, const RDouble *pos)
    {
        KDHyperRect *rect;
        KDNode *result;
        KDRes *rset;
        RDouble distSQ;

        if (!kdTree)
        {
            return NULL;
        }

        if (!kdTree->rect)
        {
            return NULL;
        }

        // Allocate result set
        try {
            rset = new KDRes;
        }
        catch (bad_alloc) {
            return NULL;
        }

        try {
            rset->rlist = new KDResNode;
        }
        catch (bad_alloc) {
            delete rset;
            return NULL;
        }

        rset->rlist->next = NULL;
        //rset->kdTree = kdTree;

        rect = DuplicateHyperRect(kdTree->rect);
        if (rect == NULL)
        {
            FreeKDRes(rset);
            return NULL;
        }

        result = kdTree->root;
        distSQ = 0.0;
        for (int i = 0; i < kdTree->dim; ++ i)
        {
            distSQ += SQ(result->pos[i] - pos[i]);
        }

        KDNearestI(kdTree->root, pos, &result, &distSQ, rect);

        //rect->Free();
        delete rect;

        if (result)
        {
            if (InsertRList(rset->rlist, result, -1.0) == -1)
            {
                FreeKDRes(rset);
                return NULL;
            }
            rset->size = 1;
            rset->Rewind();
            return rset;
        }
        else
        {
            FreeKDRes(rset);
            return NULL;
        }
    }
    int InsertRList(KDResNode *list, KDNode *item, RDouble distSQ)
    {
        KDResNode *rNode = new KDResNode;
        if (!rNode)
        {
            return -1;
        }
        rNode->item = item;
        rNode->distSQ = distSQ;

        if (distSQ >= 0.0)
        {
            while (list->next && list->next->distSQ < distSQ)
            {
                list = list->next;
            }
        }

        rNode->next = list->next;
        list->next = rNode;
        return 0;

    }

    KDRes * NearestRange(KDTree *kdTree, const RDouble *pos, RDouble range)
    {
        KDRes *rset = new KDRes;
        if (!rset)
        {
            return NULL;
        }

        try {
            rset->rlist = new KDResNode;
        }
        catch (bad_alloc) {
            delete rset;
            return NULL;
        }

        rset->rlist->next = NULL;
        rset->kdTree = kdTree;
        int ret = 0;
        if ((ret = FindNearest(kdTree->root, pos, range, rset->rlist, 0, kdTree->dim)) == -1)
        {
            FreeKDRes(rset);
            return NULL;
        }
        rset->size = ret;
        rset->Rewind();
        return rset;
    }

}