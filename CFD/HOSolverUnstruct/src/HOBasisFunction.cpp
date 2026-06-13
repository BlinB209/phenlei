#include "HODefine.h"
#include "HOShapeFunction.h"

using PHSPACE::RDouble;

namespace HOUnstruct
{
//! calculate monomial basis function
//! @param[in]   refCoord   coordinates
//! @param[in]  i        the i-th basis
//!
//! @return      basis function value
RDouble MonomialBasis(const PHSPACE::RDouble * refCoord, unsigned int i)
{
    const RDouble xi = refCoord[0];
    const RDouble et = refCoord[1];
    const RDouble zt = refCoord[2];

    //(i < ((order)+1)*((order)+2)*((order)+3)/6);

    //! monomials. since they are hierarchic we only need one case block.
    switch (i)
    {
        //! constant
        case 0:
            return 1.0;
        
        //! linears
        case 1:
            return xi;
        
        case 2:
            return et;
        
        case 3:
            return zt;
        
        //! quadratics
        case 4:
            return xi*xi;
        
        case 5:
            return xi*et;
        
        case 6:
            return et*et;
        
        case 7:
            return xi*zt;
        
        case 8:
            return zt*et;
        
        case 9:
            return zt*zt;
        
        default:
            return 0.0;
    }
}

//! calculate monomial basis function 's grad
//! @param[in]   refCoord   coordinates
//! @param[in]  i        the i-th basis
//! @param[in]  j        the j-th grad 0,1,2 means dxi, dzt, dzt
//!
//! @return      basis function grad value
RDouble MonomialBasisGrad(const PHSPACE::RDouble * refCoord, const unsigned int i, const unsigned int j)
{
    const RDouble xi = refCoord[0];
    const RDouble et = refCoord[1];
    const RDouble zt = refCoord[2];

    //! monomials. since they are hierarchic we only need one case block.
    switch (j)
    {
        //! d()/dxi
        case 0:
            {
                switch (i)
                {
                    //! constant
                    case 0:
                        return 0.;

                    //! linear
                    case 1:
                        return 1.;

                    case 2:
                        return 0.;

                    case 3:
                        return 0.;

                    //! quadratic
                    case 4:
                        return 2.*xi;

                    case 5:
                        return et;

                    case 6:
                        return 0.;

                    case 7:
                        return zt;

                    case 8:
                        return 0.;

                    case 9:
                        return 0.;

                    default:
                        return 0.0;
                }
            }
        
        //! d()/deta
        case 1:
            {
                switch (i)
                {
                    //! constant
                    case 0:
                        return 0.;

                    //! linear
                    case 1:
                        return 0.;

                    case 2:
                        return 1.;

                    case 3:
                        return 0.;

                    //! quadratic
                    case 4:
                        return 0.;

                    case 5:
                        return xi;

                    case 6:
                        return 2.*et;

                    case 7:
                        return 0.;

                    case 8:
                        return zt;

                    case 9:
                        return 0.;

                    default:
                        return 0.0;
                }
            }

        //! d()/dzeta
        case 2:
            {
                switch (i)
                {
                    //! constant
                    case 0:
                        return 0.;

                    //! linear
                    case 1:
                        return 0.;

                    case 2:
                        return 0.;

                    case 3:
                        return 1.;

                    //! quadratic
                    case 4:
                        return 0.;

                    case 5:
                        return 0.;

                    case 6:
                        return 0.;

                    case 7:
                        return xi;

                    case 8:
                        return et;
 
                    case 9:
                        return 2.*zt;

                    default:
                        return 0.0;
                }
            }

        default:
            return 0.0;    //! ("Invalid shape function derivative j = " << j);
    }
}

//! calculate monomial basis function
//! @param[in]   refCoord   coordinates
//! @param[in]  i        the i-th basis
//!
//! @return      basis function value
RDouble MonomialBasis2d(const PHSPACE::RDouble * refCoord, const unsigned int i)
{
    const RDouble xi = refCoord[0];
    const RDouble et = refCoord[1];

    //(i < ((order)+1)*((order)+2)/2);

    // monomials. since they are hierarchic we only need one case block.
    switch (i)
    {
    //! constant
    case 0:
        return 1.0;

    //! linears
    case 1:
        return xi;

    case 2:
        return et;

    //! quadratics
    case 3:
        return xi*xi;

    case 4:
        return xi*et;

    case 5:
        return et*et;

    //! cubics
    case 6:
        return xi*xi*xi;

    case 7:
        return xi*xi*et;

    case 8:
        return xi*et*et;

    case 9:
        return et*et*et;

    default:
        return 0.0;
    }
}
}

#include <iostream>
#include <cmath>
//using namespace std;

//using namespace HOUnstruct;

namespace HOUnstruct
{
//! test shape function, print some info.
//! @param[in]   
//! @param[out]  
//!
//! @return        void
void TestBasisFunction()
{
    const RDouble refCoord[3] = { 0.1, 0.2, 0.3 };

    std::cout << "monomial basis function:" << std::endl;
    for (int i = 0; i < 10; ++i)
    {
        std::cout <<  MonomialBasis(refCoord, i) << " ";
    }
    std::cout << std::endl;

    for (int j = 0; j < 3; ++j)
    {
        std::cout << "monomial basis function grad:" << j << std::endl;

        for (int i = 0; i < 10; ++i)
        {
            std::cout <<  MonomialBasisGrad(refCoord, i, j) << " ";
        }
        std::cout << std::endl;
    }

    /*/ out
    monomial basis function:
    1 0.1 0.2 0.3 0.01 0.02 0.04 0.03 0.06 0.09
    monomial basis function grad:0
    0 1 0 0 0.2 0.2 0 0.3 0 0
    monomial basis function grad:1
    0 0 1 0 0 0.1 0.4 0 0.3 0
    monomial basis function grad:2
    0 0 0 1 0 0 0 0.1 0.2 0.6
    */

    //char c;
    //cin.get(c);
}

}