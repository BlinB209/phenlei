#include "Geo_DynamicGridMetrics_Struct.h"

namespace PHSPACE
{

Geo_DynamicGridMetrics_Struct::Geo_DynamicGridMetrics_Struct()
{
    voln = NULL;
    xfv = NULL;
    yfv = NULL;
    zfv = NULL;
    xcv = NULL;
    ycv = NULL;
    zcv = NULL;
    vgn = NULL;
}

Geo_DynamicGridMetrics_Struct::~Geo_DynamicGridMetrics_Struct()
{
    DeallocateAll();
}

}