#include "GridType.h"
#include "Geo_StructGrid.h"
#include "Geo_UnstructGrid.h"

namespace PHSPACE
{

Grid * CreateGridGeneral(int gridtype, GridID *index, int level, int nDim)
{
    if (gridtype == UNSTRUCTGRID)
    {
        Grid *grid = new UnstructGrid();
        grid->InitGrid(index,level,nDim,UNSTRUCTGRID);
        return grid;
    }
    else if (gridtype == STRUCTGRID)
    {
        Grid *grid = new StructGrid();
        grid->InitGrid(index,level,nDim,STRUCTGRID);
        return grid;
    }
    cout << "no this type grid\n";
    return 0;
}


int GetSystemGridType(set<int> &grid_type)
{
    if (grid_type.size() > 1)
    {
        return MIXGRID;
    }
    else
    {
        return *grid_type.begin();
    }
}

}
