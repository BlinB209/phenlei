#include "Geo_NodeTopo_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_NodeTopo_Unstruct::Geo_NodeTopo_Unstruct()
{
    cell_number_of_each_node = 0;
    node2cell                = 0;
    node2CellArray           = 0;
    nTotalFacesOfEachNode    = 0;
    node2Face                = 0;
    node2FaceArray           = 0;
}

LIB_EXPORT Geo_NodeTopo_Unstruct::~Geo_NodeTopo_Unstruct()
{
    if (cell_number_of_each_node != 0)
    {
        delete [] cell_number_of_each_node;
        cell_number_of_each_node = 0;
    }

    if (node2cell != 0)
    {
        delete [] node2cell;
        node2cell = 0;
    }

    if (node2CellArray != 0)
    {
        delete [] node2CellArray;
        node2CellArray = 0;
    }

    if (!nTotalFacesOfEachNode)
    {
        delete [] nTotalFacesOfEachNode;    nTotalFacesOfEachNode = nullptr;
    }

    if (!node2Face)
    {
        delete [] node2Face;    node2Face = nullptr;
    }

    if (!node2FaceArray)
    {
        delete [] node2FaceArray;    node2FaceArray = nullptr;
    }
}

void Geo_NodeTopo_Unstruct::SetNTotalFacesOfEachNode(int *nFPN)
{
    if (this->nTotalFacesOfEachNode && this->nTotalFacesOfEachNode != nFPN)
    {
        delete [] this->nTotalFacesOfEachNode;    nTotalFacesOfEachNode = nullptr;
    }
    this->nTotalFacesOfEachNode = nFPN;
}

void Geo_NodeTopo_Unstruct::SetNode2Face(int *n2F)
{
    if (this->node2Face && this->node2Face != n2F)
    {
        delete [] this->node2Face;    this->node2Face = nullptr;
    }
    this->node2Face = n2F;
}

int * Geo_NodeTopo_Unstruct::GetNTotalFacesOfEachNode() const
{
    return nTotalFacesOfEachNode;
}

int * Geo_NodeTopo_Unstruct::GetNode2Face() const
{
    return node2Face;
}

void Geo_NodeTopo_Unstruct::SetNode2FaceArray(int **node2FaceArrayIn)
{
    if (this->node2FaceArray && this->node2FaceArray != node2FaceArrayIn)
    {
        delete [] this->node2FaceArray;    node2FaceArray = nullptr;
    }
    this->node2FaceArray = node2FaceArrayIn;
}

int ** Geo_NodeTopo_Unstruct::GetNode2faceArray() const
{
    return this->node2FaceArray;
}
}