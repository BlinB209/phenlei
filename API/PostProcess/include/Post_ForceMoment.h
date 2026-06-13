//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Post_ForceMoment.h
//! @brief     Offer some subroutine to calculate moment of force
//! @author    Zhang Jian

#include "Precision.h"

namespace PHSPACE
{

class Post_ForceMoment
{
public:
    Post_ForceMoment(){};
    ~Post_ForceMoment(){};

public:
    //! Calculate the moment from a point normal to a line.
    //! @param[in] point    coordinates of the point.
    //! @param[in] linePoint1/2    two end point of the line.
    //! @param[in] force    force at the point.
    RDouble ComputeMoment(RDouble *point, RDouble *linePoint1, RDouble *linePoint2, RDouble *force);

private:
    //! Calculate the distance from a point normal to a line.
    //! @param[in] point    coordinates of the point.
    //! @param[in] linePoint1/2    two end point of the line.
    RDouble DistPoint2Line(RDouble *point, RDouble *linePoint1, RDouble *linePoint2);

    //! Give a point and a line, calculate the plane's normal.
    //! @param[out] normalX/Y/Z.
    void PlaneNormal(RDouble *point, RDouble *linePoint1, RDouble *linePoint2, RDouble &normalX, RDouble &normalY, RDouble &normalZ);
};

}