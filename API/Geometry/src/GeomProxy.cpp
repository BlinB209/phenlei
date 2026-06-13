#include "GeomProxy.h"

using namespace std;
namespace PHSPACE
{
    GeomProxy::GeomProxy()
    {
        xfn  = NULL;
        yfn  = NULL;
        zfn  = NULL;
        vgn  = NULL;
        area = NULL;
    }

    GeomProxy::~GeomProxy()
    {
        delete [] xfn;
        delete [] yfn;
        delete [] zfn;
        delete [] vgn;
        delete [] area;
    }

    void GeomProxy::Create(int nsize)
    {
        xfn  = new value_type[nsize];
        yfn  = new value_type[nsize];
        zfn  = new value_type[nsize];
        vgn  = new value_type[nsize];
        area = new value_type[nsize];
    }

    GeomProxy::iterator GeomProxy::GetFaceNormalX()
    {
        return xfn;
    }

    GeomProxy::iterator GeomProxy::GetFaceNormalY()
    {
        return yfn;
    }

    GeomProxy::iterator GeomProxy::GetFaceNormalZ()
    {
        return zfn;
    }

    GeomProxy::iterator GeomProxy::GetFaceArea()
    {
        return area;
    }

    GeomProxy::iterator GeomProxy::GetFaceVelocity()
    {
        return vgn;
    }
}

