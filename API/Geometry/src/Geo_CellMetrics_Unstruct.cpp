#include "Geo_CellMetrics_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_CellMetrics_Unstruct::Geo_CellMetrics_Unstruct()
{
    xcc = NULL;
    ycc = NULL;
    zcc = NULL;
    vol = NULL;

    cellSkewness = NULL;

    largestLocalGridLength = NULL;
    subgridLength = NULL;
}

LIB_EXPORT Geo_CellMetrics_Unstruct::~Geo_CellMetrics_Unstruct()
{
    delete [] xcc; xcc = NULL;
    delete [] ycc; ycc = NULL;
    delete [] zcc; zcc = NULL;
    delete [] vol; vol = NULL;

    delete [] cellSkewness; cellSkewness = NULL;

    delete [] largestLocalGridLength; largestLocalGridLength = NULL;
    delete [] subgridLength; subgridLength = NULL;
}

}