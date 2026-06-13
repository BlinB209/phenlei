#include "Geo_FaceTopo_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_FaceTopo_Unstruct::Geo_FaceTopo_Unstruct()
{
    left_cell_of_face  = NULL;
    right_cell_of_face = NULL;
    node_number_of_each_face = NULL;
    face2node = NULL;

    face2nodeSubscript = NULL;
    face2nodeArray = NULL;
}

LIB_EXPORT Geo_FaceTopo_Unstruct::~Geo_FaceTopo_Unstruct()
{
    delete [] left_cell_of_face;    left_cell_of_face  = NULL;
    delete [] right_cell_of_face;    right_cell_of_face = NULL;
    delete [] node_number_of_each_face;    node_number_of_each_face = NULL;
    delete [] face2node;    face2node = NULL;

    if (face2nodeSubscript != 0)
    {
        delete [] face2nodeSubscript;    face2nodeSubscript = NULL;
    }

    if (face2nodeArray != 0)
    {
        delete [] face2nodeArray;    face2nodeArray = NULL;
    }
}

}