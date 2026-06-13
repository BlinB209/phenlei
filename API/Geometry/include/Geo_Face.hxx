#include "Math_BasisFunction.h"
using namespace std;

namespace PHSPACE
{

inline Geo_Face::Geo_Face()
{

}

inline Geo_Face::Geo_Face(const int &le, const int &re, const int &id, int maxFaceID)
{
    this->le     = le;
    this->re     = re;
    this->faceID = id;
    this->maxFaceID = maxFaceID;
}

inline Geo_Face::Geo_Face(const Geo_Face &rhs)
{
    this->le     = rhs.le;
    this->re     = rhs.re;
    this->faceID = rhs.faceID;
    this->maxFaceID = rhs.maxFaceID;
}

inline Geo_Face::~Geo_Face(void)
{

}

inline Geo_Face & Geo_Face::operator = (const Geo_Face &rhs)
{
    if (this == &rhs) return *this;

    this->le     = rhs.le;
    this->re     = rhs.re;
    this->faceID = rhs.faceID;
    this->maxFaceID = rhs.maxFaceID;

    return *this;
}

inline bool Geo_Face::operator < (const Geo_Face &rhs) const
{
    int dx = le - rhs.le;
    int dy = re - rhs.re;

    if (ABS(dx) > 0) return le < rhs.le;
    if (ABS(dy) > 0) return re < rhs.re;

    return false;
}

inline int Geo_Face::hash_key() const
{
    unsigned long temp1, temp2;
    int hashNumber;

    temp1 = (le + re) * 23;
    temp2 = le * re;
    hashNumber = (temp1 + temp2) % maxFaceID;
    return hashNumber;
}

}
