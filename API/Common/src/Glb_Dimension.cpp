#include "Glb_Dimension.h"

namespace PHSPACE
{
int dimension = THREE_D;

LIB_EXPORT int GetDim()
{
    return dimension;
}

LIB_EXPORT void SetDim(int dim)
{
    dimension = dim;
}

}