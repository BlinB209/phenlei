inline void Geo_FaceTopo_Unstruct::SetLeftCellOfFace(int *left_cell_of_face_in) 
{ 
//  if ( this->left_cell_of_face && this->left_cell_of_face != left_cell_of_face_in ) 
//  {
//      delete [] this->left_cell_of_face;
//      this->left_cell_of_face = 0;
//  }
    this->left_cell_of_face = left_cell_of_face_in;
}

inline void Geo_FaceTopo_Unstruct::SetRightCellOfFace(int *right_cell_of_face_in)
{
//  if (this->right_cell_of_face && this->right_cell_of_face != right_cell_of_face_in)
//  {
//      delete [] this->right_cell_of_face;
//      this->right_cell_of_face = 0;
//  }
    this->right_cell_of_face = right_cell_of_face_in;
}

inline void Geo_FaceTopo_Unstruct::SetNodeNumberOfEachFace(int *node_number_of_each_face)
{
//  if (this->node_number_of_each_face && this->node_number_of_each_face != node_number_of_each_face)
//  {
//      delete [] this->node_number_of_each_face;
//      this->node_number_of_each_face = 0;
//  }
    this->node_number_of_each_face = node_number_of_each_face;
}

inline void Geo_FaceTopo_Unstruct::SetFace2Node(int *face2node)
{
//  if (this->face2node && this->face2node != face2node)
//  {
//      delete [] this->face2node;
//      this->face2node = 0;
//  }
    this->face2node = face2node;
}

inline void Geo_FaceTopo_Unstruct::SetFace2NodeSubscript(long long int *face2nodeSubscript)
{
    this->face2nodeSubscript = face2nodeSubscript;
}

inline int * Geo_FaceTopo_Unstruct::GetLeftCellOfFace() const
{
    return left_cell_of_face;
}

inline int * Geo_FaceTopo_Unstruct::GetRightCellOfFace() const
{
    return right_cell_of_face;
}

inline int * Geo_FaceTopo_Unstruct::GetNodeNumberOfEachFace() const
{
    return node_number_of_each_face;
}

inline int * Geo_FaceTopo_Unstruct::GetFace2Node() const
{
    return face2node;
}

inline int ** Geo_FaceTopo_Unstruct::GetFace2NodeArray() const
{
    return face2nodeArray;
}

inline long long int * Geo_FaceTopo_Unstruct::GetFace2NodeSubscript() const
{
    return face2nodeSubscript;
}

inline void Geo_FaceTopo_Unstruct::SetFace2NodeArray(int **face2NodeArrayIn)
{
    this->face2nodeArray = face2NodeArrayIn;
}