#include "GlobalDataBase.h"
#include "Geo_CellMetrics_Struct.h"

namespace PHSPACE
{
Geo_CellMetrics_Struct::Geo_CellMetrics_Struct()
{
    xcc = NULL;
    ycc = NULL;
    zcc = NULL;
    vol = NULL;

    xlen = NULL;
    ylen = NULL;
    zlen = NULL;

    jacobian = NULL;

    largestLocalGridLength = NULL;
    smallestLocalGridLength = NULL;
    subgridLength = NULL;
}

Geo_CellMetrics_Struct::~Geo_CellMetrics_Struct()
{
    DeallocateAll();

    delete largestLocalGridLength; largestLocalGridLength = NULL;
    delete smallestLocalGridLength; smallestLocalGridLength = NULL;
    delete subgridLength; subgridLength = NULL;
}

void Geo_CellMetrics_Struct::AllocateMetrics(Range &I, Range &J, Range &K)
{
    DeallocateAll();
    xcc = new RDouble3D(I, J, K, fortranArray);
    ycc = new RDouble3D(I, J, K, fortranArray);
    zcc = new RDouble3D(I, J, K, fortranArray);
    vol = new RDouble3D(I, J, K, fortranArray);

    xlen = new RDouble3D(I, J, K, fortranArray);
    ylen = new RDouble3D(I, J, K, fortranArray);
    zlen = new RDouble3D(I, J, K, fortranArray);

    int SolverStructOrder = PHSPACE::GlobalDataBase::GetIntParaFromDB("SolverStructOrder");
    if (SolverStructOrder > 2)
    {
        jacobian = new RDouble3D(I, J, K, fortranArray);
    }
}
}