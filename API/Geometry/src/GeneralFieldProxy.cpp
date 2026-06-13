#include "GeneralFieldProxy.h"
#include "Pointer.h"

namespace PHSPACE
{

GeneralFieldProxy::GeneralFieldProxy()
{
    ndim      = 0;
    nsize     = 0;
    field     = 0;
    del_field = false;
}

GeneralFieldProxy::~GeneralFieldProxy()
{
    if (del_field)
    {
        DelPointer2(field);
    }
}

void GeneralFieldProxy::Create(int ndim, int nsize)
{
    this->ndim  = ndim;
    this->nsize = nsize;
    del_field   = true;
    field = NewPointer2<RDouble>(ndim, nsize);
}

RDouble ** GeneralFieldProxy::GetField()
{
    return field;
}
}

