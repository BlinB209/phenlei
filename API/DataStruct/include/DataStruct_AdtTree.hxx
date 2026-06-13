template < typename T, typename U >
DataStruct_AdtTree<T, U>::DataStruct_AdtTree(int ndim)
{
    this->ndim = ndim;
    pmin = new U[ndim];
    pmax = new U[ndim];
    for (int i = 0; i < ndim; ++ i)
    { 
        pmin[i] = 0.0;
        pmax[i] = 1.0;
    }
    root = 0;
}

template < typename T, typename U >
DataStruct_AdtTree<T, U>::DataStruct_AdtTree(int ndim, U *pmin, U *pmax)
{
    this->ndim = ndim;
    this->pmin = new U[ndim];
    this->pmax = new U[ndim];
    for (int i = 0; i < ndim; ++ i)
    { 
        this->pmin[i] = pmin[i];
        this->pmax[i] = pmax[i];
    }
    root = 0;
}

template < typename T, typename U >
DataStruct_AdtTree<T, U>::~DataStruct_AdtTree()
{  
    delete [] pmin;
    delete [] pmax;
    delete root;
}

// Add an Adt node to the AdtTree 
template < typename T, typename U >
void DataStruct_AdtTree<T, U>::AddNode(AdtNode *node)
{
    //这里似乎可以加速,不要老new来new去的
    U *nwmin = new U[ndim];
    U *nwmax = new U[ndim];
    memcpy(nwmin, this->pmin, ndim*sizeof(U));
    memcpy(nwmax, this->pmax, ndim*sizeof(U));

    if (root == 0)
    {
        root = node;
    }
    else
    {
        root->AddNode(node, nwmin, nwmax, ndim);
    }
        
    delete [] nwmin;
    delete [] nwmax;
}

// Find All nodes inside the region (pmin,pmax) from the tree
template < typename T, typename U >
void DataStruct_AdtTree<T, U>::FindNodesInRegion(U *pmin, U *pmax, AdtNodeList &ld, const size_t &sizeLimit)
{
    U *nwmin = new U[ndim];
    U *nwmax = new U[ndim];
    memcpy(nwmin, this->pmin, ndim*sizeof(U));
    memcpy(nwmax, this->pmax, ndim*sizeof(U));

    if (root)
    {
        root->FindNodesInRegion(pmin, pmax, nwmin, nwmax, ndim, ld, sizeLimit);
    }
    delete [] nwmin;
    delete [] nwmax;
}

// Get the min coordinates of the tree
template < typename T, typename U >
U * DataStruct_AdtTree<T, U>::GetMin() const
{
    return pmin;
}

// Get the max coordinates of the tree
template < typename T, typename U >
U * DataStruct_AdtTree<T, U>::GetMax() const
{
    return pmax;
}

template < typename T, typename U >
void DataStruct_AdtTree<T, U>::LevelTraverse() const
{
    queue< AdtNode * > nodeQueue;
    AdtNode *start = root;
    if (start == 0)
    {
        return;
    }

    nodeQueue.push(start);
    while (!nodeQueue.empty())
    {
        start = nodeQueue.front();
        nodeQueue.pop();

        //cout<<start->item<<" ";
        if (start->left)
        {
            nodeQueue.push(start->left);
        }
        if (start->right)
        {
            nodeQueue.push(start->right);
        }
    }
}

template < typename T, typename U >
int  DataStruct_AdtTree<T, U>::GetNodeNum()
{ 
    if (root)
    {
        return root->GetNodeNum();
    }
    else
    {
        return 0;
    }
}

template < typename T, typename U >
typename DataStruct_AdtTree<T, U>::AdtNode * DataStruct_AdtTree<T,U>::GetRoot() const
{
    return root;
}

