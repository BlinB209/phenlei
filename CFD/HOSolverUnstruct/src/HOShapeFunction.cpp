#include "HODefine.h"
#include "HOShapeFunction.h"

using PHSPACE::RDouble;

namespace HOUnstruct
{
//! calculate shape function
//! @param[in]   refCoord   reference coordinates
//! @param[out]  ret        shape function
//!
//! @return        void

void ShapeFunctionOfTri(const RDouble * refCoord, RDouble * ret)
{
    RDouble xi = refCoord[0];
    RDouble et = refCoord[1];

    ret[0] = 1.0 - xi - et;
    ret[1] = xi;
    ret[2] = et;
}

void ShapeFunctionOfTriGrad(RDouble * ret)
{
    ret[0] = -1.0;
    ret[1] = -1.0;

    ret[2] =  1.0;
    ret[3] =  0.0;

    ret[4] =  0.0;
    ret[5] =  1.0;
}

void ShapeFunctionOfQuad(const RDouble * refCoord, RDouble * ret)
{
    RDouble xi = refCoord[0];
    RDouble et = refCoord[1];

    ret[0] = (1.0 - et) * (1.0 - xi);
    ret[1] = xi * (1.0 - et);
    ret[2] = xi * et;
    ret[3] = (1.0 - xi) * et;
}

void ShapeFunctionOfQuadGrad(const RDouble * refCoord, RDouble * ret)
{
    RDouble xi = refCoord[0];
    RDouble et = refCoord[1];

    ret[0] = et - 1.0;
    ret[1] = xi - 1.0;

    ret[2] =  1.0 - et;
    ret[3] = -xi;

    ret[4] =  et;
    ret[5] =  xi;

    ret[6] = -et;
    ret[7] =  1.0 - xi;
}

void ShapeFunctionOfTetr(const RDouble * refCoord, RDouble * ret)
{
    RDouble xi = refCoord[0];
    RDouble et = refCoord[1];
    RDouble zt = refCoord[2];

    ret[0] = 1.0 - xi - et - zt;
    ret[1] =       xi;
    ret[2] =            et;
    ret[3] =                 zt;
}

void ShapeFunctionOfTetrGrad(RDouble * ret)
{
    ret[ 0] = -1.0;
    ret[ 1] = -1.0;
    ret[ 2] = -1.0;

    ret[ 3] =  1.0;
    ret[ 4] =  0.0;
    ret[ 5] =  0.0;

    ret[ 6] =  0.0;
    ret[ 7] =  1.0;
    ret[ 8] =  0.0;

    ret[ 9] =  0.0;
    ret[10] =  0.0;
    ret[11] =  1.0;
}

void ShapeFunctionOfHex(const RDouble * refCoord, RDouble * ret)
{
    RDouble xi = refCoord[0];
    RDouble et = refCoord[1];
    RDouble zt = refCoord[2];

    ret[0] = (1 - xi) * (1 - et) * (1 - zt);
    ret[1] = xi  * (1 - et) * (1 - zt);
    ret[2] = xi  * et  * (1 - zt);
    ret[3] = (1 - xi) * et  * (1 - zt);
    ret[4] = (1 - xi) * (1 - et) * zt;
    ret[5] = xi  * (1 - et) * zt;
    ret[6] = xi  * et  * zt;
    ret[7] = (1 - xi) * et  * zt;
}

void ShapeFunctionOfHexGrad(const RDouble * refCoord, RDouble * ret)
{
    RDouble xi = refCoord[0];
    RDouble et = refCoord[1];
    RDouble zt = refCoord[2];

    ret[0] = -(1 - et) * (1 - zt);
    ret[1] = -(1 - xi) * (1 - zt);
    ret[2] = -(1 - xi) * (1 - et);

    ret[3] =  (1 - et) * (1 - zt);
    ret[4] = -xi  * (1 - zt);
    ret[5] = -xi  * (1 - et);

    ret[6] =  et  * (1 - zt);
    ret[7] =  xi  * (1 - zt);
    ret[8] = -xi  * et;

    ret[9 ] = -et  * (1 - zt);
    ret[10] =  (1 - xi) * (1 - zt);
    ret[11] = -(1 - xi) * et;

    ret[12] = -(1 - et) * zt;
    ret[13] = -(1 - xi) * zt;
    ret[14] =  (1 - xi) * (1 - et);

    ret[15] =  (1 - et) * zt;
    ret[16] = -xi  * zt;
    ret[17] =  xi  * (1 - et);

    ret[18] =  et  * zt;
    ret[19] =  xi  * zt;
    ret[20] =  xi  * et;

    ret[21] = -et  * zt;
    ret[22] =  (1 - xi) * zt;
    ret[23] =  (1 - xi) * et;
}

void ShapeFunctionOfPrism(const RDouble * refCoord, RDouble * ret)
{
    RDouble xi = refCoord[0];
    RDouble et = refCoord[1];
    RDouble zt = refCoord[2];

    ret[0] = (1.0 - xi - et) * (1.0 - zt);
    ret[1] = xi * (1.0 - zt);
    ret[2] = et * (1.0 - zt);
    ret[3] = (1.0 - xi - et) * zt;
    ret[4] = xi * zt;
    ret[5] = et * zt;
}

void ShapeFunctionOfPrismGrad(const RDouble * refCoord, RDouble * ret)
{
    RDouble xi = refCoord[0];
    RDouble et = refCoord[1];
    RDouble zt = refCoord[2];

    ret[0] =  zt-1.0;
    ret[1] =  zt-1.0;
    ret[2] =  xi+et-1.0;

    ret[3] =  1.0-zt;
    ret[4] =  0.0;
    ret[5] = -xi;

    ret[6] =  0.0;
    ret[7] =  1.0-zt;
    ret[8] = -et;

    ret[9 ] = -zt;
    ret[10] = -zt;
    ret[11] =  1.0-xi-et;

    ret[12] =  zt;
    ret[13] =  0.0;
    ret[14] =  xi;

    ret[15] =  0.0;
    ret[16] =  zt;
    ret[17] =  et;
}

void ShapeFunctionOfPyra(const RDouble * refCoord, RDouble * ret)
{
    RDouble xi = refCoord[0];
    RDouble et = refCoord[1];
    RDouble zt = refCoord[2];

    ret[0] = (1.0 - xi) * (1.0 - et) * (1.0 - zt);
    ret[1] = xi * (1.0 - et) * (1.0 - zt);
    ret[2] = xi * et * (1.0 - zt);
    ret[3] = (1.0 - xi) * et * (1.0 - zt);
    ret[4] = zt;
}

void ShapeFunctionOfPyraGrad(const RDouble * refCoord, RDouble * ret)
{
    RDouble xi = refCoord[0];
    RDouble et = refCoord[1];
    RDouble zt = refCoord[2];

    ret[0] = -(1.0 - et) * (1.0 - zt);
    ret[1] = -(1.0 - xi) * (1.0 - zt);
    ret[2] = -(1.0 - xi) * (1.0 - et);

    ret[3] =  (1.0 - et) * (1.0 - zt);
    ret[4] = -xi  * (1.0 - zt);
    ret[5] = -xi  * (1.0 - et);

    ret[6] =  et  * (1.0 - zt);
    ret[7] =  xi  * (1.0 - zt);
    ret[8] = -xi  * et;

    ret[9 ] = -et  * (1.0 - zt);
    ret[10] =  (1.0 - xi) * (1.0 - zt);
    ret[11] = -(1.0 - xi) * et;

    ret[12] =  0.0;
    ret[13] =  0.0;
    ret[14] =  1.0;
}

}

#include <iostream>
#include <cmath>
//using namespace std;

//using namespace HOUnstruct;

namespace HOUnstruct
{
//! test input 0 or 1
//! @param[in]   valToTest
//! @param[in]   eps = HO_SMALL
//! @param[out]  
//!
//! @return        int , result of test
//! @retval        £¨return£©0  valToTest near 0
//! @retval        £¨return£©1  valToTest near 1
//! @retval        £¨return£©2  unknown

int Test0Or1(RDouble & valToTest, RDouble eps = HO_SMALL)
{
    if (abs(valToTest) < eps)
    {
        if (abs(valToTest-1) < eps) return 2;
        else return 0;
    }
    else
    {
        if (abs(valToTest-1) < eps) return 1;
        else return 2;
    }
}

}