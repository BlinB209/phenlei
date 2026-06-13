inline Grid ** Pre_GridConversion::GetGrids()
{
    return this->grids;
}

inline int Pre_GridConversion::GetNumberofBlocks()
{
    return this->nBlocks;
}

inline void Pre_GridConversion::SetNumberofBlocks(int dataIn)
{
    this->nBlocks = dataIn;
}