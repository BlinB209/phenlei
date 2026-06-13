inline void Geo_CellTopo_Unstruct::SetFaceNumberOfEachCell(int *face_number_of_each_cell)
{
    if (this->face_number_of_each_cell && this->face_number_of_each_cell != face_number_of_each_cell)
    {
        delete [] this->face_number_of_each_cell;
    }
    this->face_number_of_each_cell = face_number_of_each_cell;
}

inline void Geo_CellTopo_Unstruct::SetCell2Face(int **cell2faceIn)
{
    if (this->cell2face && this->cell2face != cell2faceIn)
    {
        delete [] cell2face[0];
        delete [] cell2face;
    }
    cell2face = cell2faceIn;
}

inline void Geo_CellTopo_Unstruct::SetNodeNumberOfEachCell(int *node_number_of_each_cell)
{
    if (this->node_number_of_each_cell && this->node_number_of_each_cell != node_number_of_each_cell)
    {
        delete [] this->node_number_of_each_cell; 
    }
    this->node_number_of_each_cell = node_number_of_each_cell;
}

inline void Geo_CellTopo_Unstruct::SetCell2Node(int *cell2node)
{
    if (this->cell2node && this->cell2node != cell2node)
    {
        delete [] this->cell2node;
    }
    this->cell2node = cell2node;
}

inline void Geo_CellTopo_Unstruct::SetCell2NodeArray(int **cell2nodeArrayIn)
{
    if (this->cell2nodeArray && this->cell2nodeArray != cell2nodeArray)
    {
        delete [] cell2nodeArray[0];
        delete [] cell2nodeArray;
    }
    this->cell2nodeArray = cell2nodeArrayIn;
}

inline void Geo_CellTopo_Unstruct::SetNumberOfNeighborCell(int *number_of_neighbor_cell)
{
    if (this->number_of_neighbor_cell && this->number_of_neighbor_cell != number_of_neighbor_cell)
    {
        delete [] this->number_of_neighbor_cell;
    }
    this->number_of_neighbor_cell = number_of_neighbor_cell;
}

inline void Geo_CellTopo_Unstruct::SetCell2Cell(vector<int> *cell2cell)
{
    if (this->cell2cell && this->cell2cell != cell2cell)
    {
        delete this->cell2cell;
    }
    this->cell2cell = cell2cell;
}

// GMRESVis
inline void Geo_CellTopo_Unstruct::GMRESSetNeighborCells(vector<int> * neighborCells)
{
	if(this->neighborCells && this->neighborCells != neighborCells)
	{
		delete this->neighborCells;
	}
	this->neighborCells = neighborCells;
}

// GMRESVis
inline void Geo_CellTopo_Unstruct::GMRESSetNeighborFaces(vector<int> * neighborFaces)
{
	if(this->neighborFaces && this->neighborFaces != neighborFaces)
	{
		delete this->neighborFaces;
	}
	this->neighborFaces = neighborFaces;
}

// GMRESVis
inline void Geo_CellTopo_Unstruct::GMRESSetNeighborLR(vector<int> * neighborLR)
{
	if(this->neighborLR && this->neighborLR != neighborLR)
	{
		delete this->neighborLR;
	}
	this->neighborLR = neighborLR;
}

inline int * Geo_CellTopo_Unstruct::GetFaceNumberOfEachCell() 
{
    return face_number_of_each_cell;
}

inline int ** Geo_CellTopo_Unstruct::GetCell2Face()
{
    return cell2face;
}

inline int * Geo_CellTopo_Unstruct::GetNodeNumberOfEachCell()
{
    return node_number_of_each_cell;
}

inline int * Geo_CellTopo_Unstruct::GetCell2Node()
{
    return cell2node;
}

inline int ** Geo_CellTopo_Unstruct::GetCell2NodeArray()
{
    return cell2nodeArray;
}

inline int * Geo_CellTopo_Unstruct::GetNumberOfNeighborCell()
{
    return number_of_neighbor_cell;
}

inline vector<int> * Geo_CellTopo_Unstruct::GetCell2Cell()
{
    return cell2cell;
}

// GMRESVis
inline vector<int> * Geo_CellTopo_Unstruct::GMRESGetNeighborCells()
{
	return neighborCells;
}

// GMRESVis
inline vector<int> * Geo_CellTopo_Unstruct::GMRESGetNeighborFaces()
{
	return neighborFaces;
}

// GMRESVis
inline vector<int> * Geo_CellTopo_Unstruct::GMRESGetNeighborLR()
{
	return neighborLR;
}
