inline int * Geo_NodeTopo_Unstruct::GetCellNumberOfEachNode() const
{
    return cell_number_of_each_node;
}

inline int * Geo_NodeTopo_Unstruct::GetNode2Cell() const
{
    return node2cell;
}

inline void Geo_NodeTopo_Unstruct::SetCellNumberOfEachNode(int *nCPN)
{
    if (this->cell_number_of_each_node && this->cell_number_of_each_node != nCPN)
    {
        delete [] this->cell_number_of_each_node;
    }
    this->cell_number_of_each_node = nCPN;
}

inline void Geo_NodeTopo_Unstruct::SetNode2Cell(int *n2c)
{
    if (this->node2cell && this->node2cell != n2c)
    {
        delete [] this->node2cell;
    }
    this->node2cell = n2c; 
}

inline void Geo_NodeTopo_Unstruct::SetNode2CellArray(int **node2Cell)
{
    if (this->node2CellArray && this->node2CellArray != node2Cell)
    {
        delete [] this->node2CellArray;
    }
    this->node2CellArray = node2Cell;
}

inline int ** Geo_NodeTopo_Unstruct::GetNode2CellArray() const
{
    return this->node2CellArray;
}