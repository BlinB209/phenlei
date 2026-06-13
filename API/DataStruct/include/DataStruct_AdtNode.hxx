template < typename T, typename U >
DataStruct_AdtNode<T, U>::DataStruct_AdtNode(int ndim)
{
    this->ndim = ndim;
    point = new U[ndim];
    level = 0;
    left  = 0;
    right = 0;
}

template < typename T, typename U >
DataStruct_AdtNode<T, U>::DataStruct_AdtNode(int ndim, U *coor, T data)
{
    this->ndim = ndim;
    point = new U[ndim];
    memcpy(point, coor, ndim*sizeof(U));

    level = 0;
    left  = 0;
    right = 0;
    item  = data;
}

template < typename T, typename U >
DataStruct_AdtNode<T, U>::~DataStruct_AdtNode()
{
    delete [] point;
    delete left;
    delete right;
}

template < typename T, typename U >
int DataStruct_AdtNode<T, U>::GetNodeNum()
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
}

//! Add an Adt node under the current node.
template < typename T, typename U >
void DataStruct_AdtNode<T, U>::AddNode(AdtNode *node, U *nwmin, U *nwmax, const int &ndim)
{
    int axis = level%ndim;
    U mid = 0.5 * (nwmin[axis] + nwmax[axis]);

    if (node->point[axis] <= mid)
    {
        if (left)
        {
            nwmax[axis] = mid;
            left->AddNode(node, nwmin, nwmax, ndim);
        }
        else
        {
            left = node;
            node->level = level + 1;
        }
    }
    else
    {
        if (right)
        {
            nwmin[axis] = mid;
            right->AddNode(node, nwmin, nwmax, ndim);
        }
        else
        {
            right = node;
            node->level = level + 1;
        }
    }
}

// is the current node inside region (pmin,pmax)?
template < typename T, typename U >
bool DataStruct_AdtNode<T, U>::IsInRegion(U *pmin, U *pmax, const int &ndim)
{
    for (int i = 0; i < ndim; ++ i)
    {
        if (point[i] < pmin[i] || point[i] > pmax[i])
        {
            return false;
        }
    }

    return true;
}

// ld carries all the nodes inside region (pmin,pmax)
template < typename T, typename U >
void DataStruct_AdtNode<T, U>::FindNodesInRegion(U *pmin, U *pmax, U *nwmin, U *nwmax, const int &ndim, AdtNodeList &ld, const uint_t & sizeLimit)
{
    int axis;
    U mid, temp;

    if (IsInRegion(pmin, pmax, ndim))
    {
        ld.push_back(this);
    }

    uint_t sizeInList = ld.size();

    if (sizeLimit > 0 && sizeInList >= sizeLimit)
    {
        return;
    }

    axis = level%ndim;
    mid = 0.5 * (nwmin[axis] + nwmax[axis]);

    if (left)
    {
        if (pmin[axis] <= mid && pmax[axis] >= nwmin[axis])
        {
            temp        = nwmax[axis];
            nwmax[axis] = mid;
            left->FindNodesInRegion(pmin, pmax, nwmin, nwmax, ndim, ld, sizeLimit);
            nwmax[axis] = temp;
        }
    }

    if (right)
    {
        if (pmax[axis] >= mid && pmin[axis] <= nwmax[axis])
        {
            temp        = nwmin[axis];
            nwmin[axis] = mid;
            right->FindNodesInRegion(pmin, pmax, nwmin, nwmax, ndim, ld, sizeLimit);
            nwmin[axis] = temp;
        }
    }

    return;
}