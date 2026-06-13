#include "GeometryUnit.h"
#include <cstddef>

namespace PHSPACE
{
DYNode::DYNode()
{
    node2node  = NULL;
    node2nodeK = NULL;
}

DYNode::~DYNode()
{
    delete [] node2node;
    delete [] node2nodeK;
}

DYFace::DYFace()
{
    face2node = NULL;
}
DYFace::~DYFace()
{
    delete [] face2node;
}

DYCell::DYCell()
{
    cell2face = NULL;
    cell2node = NULL;
    cell2cell = NULL;
}
DYCell::~DYCell()
{
    delete [] cell2face;
    delete [] cell2node;
    delete [] cell2cell;
}

}