#include "Geometry.h"
#include "Geo_Grid.h"

namespace PHSPACE
{

PHGeometry::PHGeometry(vector < Grid * > *grids)
{
    this->grids = grids;
}

PHGeometry::~PHGeometry()
{
    for (std::size_t i = 0; i < grids->size(); ++ i)
    {
        delete (*grids)[i];
    }
    delete grids;
}

Grid * PHGeometry::GetGrid(int level)
{
    return (*grids)[level];
}

void PHGeometry::AddGrid(Grid *grid)
{
    this->grids->push_back(grid);
}

vector < Grid * > * PHGeometry::GetVectorGrid()
{
    return grids;
}

}
