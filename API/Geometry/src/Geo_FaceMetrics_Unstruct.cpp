#include "Geo_FaceMetrics_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_FaceMetrics_Unstruct::Geo_FaceMetrics_Unstruct()
{
    xfc = NULL;
    yfc = NULL;
    zfc = NULL;

    xfn = NULL;
    yfn = NULL;
    zfn = NULL;

    area = NULL;
}

LIB_EXPORT Geo_FaceMetrics_Unstruct::~Geo_FaceMetrics_Unstruct()
{
    delete [] xfc;    xfc = NULL;
    delete [] yfc;    yfc = NULL;
    delete [] zfc;    zfc = NULL;

    delete [] xfn;    xfn = NULL;
    delete [] yfn;    yfn = NULL;
    delete [] zfn;    zfn = NULL;

    delete [] area;    area = NULL;
}

}