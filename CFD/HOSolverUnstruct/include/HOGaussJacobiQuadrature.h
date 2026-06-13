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
//! @file      HOGaussJacobiQuadrature.h
//! @brief     compute Gauss Jacobi Quadrature point for High Order Solver.
//! @author    Li Ming, Gong Xiaoquan, Zhang Jian, Wan yunbo, Xu gang.
#pragma once

#include <vector>

using namespace PHSPACE;
using namespace std;

namespace HOUnstruct
{

class CGaussJacobiQuadrature;
/*!
* \class CGaussJacobiQuadrature
* \brief Class used to determine the quadrature points of the Gauss Jacobi integration rules.
*/
class CGaussJacobiQuadrature
{
public:
    CGaussJacobiQuadrature() {}
    ~CGaussJacobiQuadrature() {}

    /*!
    * \brief Function, which serves as the API to compute the integration points
    and weights.
    * \param[in]     alpha     Parameter in the weighting function (b-x)^alpha*(x-a)^beta
    in the Gauss Jacobi rule.
    * \param[in]     beta      Parameter in the weighting function (b-x)^alpha*(x-a)^beta
    in the Gauss Jacobi rule.
    * \param[in]     a         Lower bound of the integration interval, usually -1.0.
    * \param[in]     b         Upper bound of the integration interval, usually  1.0.
    * \param[in,out] GJPoints  Location of the Gauss-Jacobi integration points.
    * \param[in,out] GJWeights Weights of the Gauss-Jacobi integration points.
    */
    void GetQuadraturePoints(const RDouble alpha,const RDouble beta,const RDouble a,const RDouble b,
                            vector<RDouble> &GJPoints, vector<RDouble> &GJWeights);
private:
    /*!
    * \brief Function in the original implementation of John Burkardt to compute
    the integration points of the Gauss-Jacobi quadrature rule.
    */
    void cdgqf(int nt, int kind, RDouble alpha, RDouble beta, RDouble t[],
        RDouble wts[]);

    /*!
    * \brief Function in the original implementation of John Burkardt to compute
    the integration points of the Gauss-Jacobi quadrature rule.
    */
    void cgqf(int nt, int kind, RDouble alpha, RDouble beta, RDouble a,
        RDouble b, RDouble t[], RDouble wts[]);

    /*!
    * \brief Function in the original implementation of John Burkardt to compute
    the integration points of the Gauss-Jacobi quadrature rule.
    */
    RDouble class_matrix(int kind, int m, RDouble alpha, RDouble beta,
        RDouble aj[], RDouble bj[]);

    /*!
    * \brief Function in the original implementation of John Burkardt to compute
    the integration points of the Gauss-Jacobi quadrature rule.
    */
    void imtqlx(int n, RDouble d[], RDouble e[], RDouble z[]);

    /*!
    * \brief Function in the original implementation of John Burkardt to compute
    the integration points of the Gauss-Jacobi quadrature rule.
    */
    void parchk(int kind, int m, RDouble alpha, RDouble beta);

    /*!
    * \brief Function in the original implementation of John Burkardt to compute
    the integration points of the Gauss-Jacobi quadrature rule.
    */
    RDouble r8_epsilon();

    /*!
    * \brief Function in the original implementation of John Burkardt to compute
    the integration points of the Gauss-Jacobi quadrature rule.
    */
    RDouble r8_sign(RDouble x);

    /*!
    * \brief Function in the original implementation of John Burkardt to compute
    the integration points of the Gauss-Jacobi quadrature rule.
    */
    void scqf(int nt, RDouble t[], int mlt[], RDouble wts[], int ndx[],
        RDouble swts[], RDouble st[], int kind, RDouble alpha,
        RDouble beta, RDouble a, RDouble b);

    /*!
    * \brief Function in the original implementation of John Burkardt to compute
    the integration points of the Gauss-Jacobi quadrature rule.
    */
    void sgqf(int nt, RDouble aj[], RDouble bj[], RDouble zemu, RDouble t[],
        RDouble wts[]);
};

void TestGaussLegendrePoints1D();

}