#include "FieldProxy.h"
#include "Pointer.h"

namespace PHSPACE
{
FieldProxy::FieldProxy()
{
    field_uns  = 0;
    field_str  = 0;
    del_field  = false;
}

FieldProxy::~FieldProxy()
{
    if (del_field)
    {
        delete field_str;
        DelPointer2(field_uns);
    }
}

RDouble ***FieldProxy::GetField_UHO()
{
    return field_uho;
}

RDouble ** FieldProxy::GetField_UNS()
{
    return field_uns;
}

RDouble4D & FieldProxy::GetField_STR()
{
    return *field_str;
}

void FieldProxy::SetField_UHO(RDouble ***field_uho, bool del_field)
{
    this->field_uho = field_uho;
    this->del_field = del_field;
}

void FieldProxy::SetField_UNS(RDouble **field_uns, bool del_field)
{
    this->field_uns = field_uns;
    this->del_field = del_field;
}

void FieldProxy::SetField_STR(RDouble4D *field_str, bool del_field)
{
    this->field_str = field_str;
    this->del_field = del_field;
}
}

