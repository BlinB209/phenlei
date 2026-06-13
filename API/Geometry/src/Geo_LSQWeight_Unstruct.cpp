#include "Geo_LSQWeight_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_LSQWeight_Unstruct::Geo_LSQWeight_Unstruct()
{
    iwt = NULL;

    ixx = NULL;
    iyy = NULL;
    izz = NULL;
    ixy = NULL;
    ixz = NULL;
    iyz = NULL;

    fMark = NULL;
}

LIB_EXPORT Geo_LSQWeight_Unstruct::~Geo_LSQWeight_Unstruct()
{
    delete [] iwt;

    delete [] ixx;
    delete [] iyy;
    delete [] izz;
    delete [] ixy;
    delete [] ixz;
    delete [] iyz;

    delete [] fMark;
}

}