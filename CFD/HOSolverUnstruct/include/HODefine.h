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
//! @file      HODefine.h
//! @brief     Constant variables definition for High Order Solver.
//! @author    xxx

#pragma once
#include <vector>
#include "Precision.h"
#include "Constants.h"
#include "cgnslib.h"    //! ElementType_t
#include "../../../PHengLEI/3rdparty/eigen/Eigen/Eigen"

using namespace std;
using namespace PHSPACE;

namespace HOUnstruct
{
#define HO_UNIT_TEST

//#define HODebug     true
#define HODebug     false

#define HO_USE_REAL_DOUBLE

//! ----------------------------------------------------------------------------*/
//!
//!   Basic constants, copy from PETSC
//!

#if defined(HO_USE_REAL___FLOAT128)
#define HO_PI                 M_PIq
#elif defined(M_PI)
#define HO_PI                 M_PI
#else
#define HO_PI                 3.14159265358979323846264338327950288419716939937510582
#endif

#if !defined(HO_USE_64BIT_INDICES)
#define HO_MAX_INT            2147483647
#define HO_MIN_INT            (-HO_MAX_INT - 1)
#else
#define HO_MAX_INT            9223372036854775807L
#define HO_MIN_INT            (-HO_MAX_INT - 1)
#endif


#if defined(HO_USE_REAL_SINGLE)
#define HO_MAX_REAL                3.40282346638528860e+38F
#define HO_MIN_REAL                -HO_MAX_REAL
#define HO_MACHINE_EPSILON         1.19209290e-07F
#define HO_SQRT_MACHINE_EPSILON    3.45266983e-04F
#define HO_SMALL                   1.e-5
#define HO_TOLERANCE               2.5e-3
#elif defined(HO_USE_REAL_DOUBLE)
#define HO_MAX_REAL                1.7976931348623157e+308
#define HO_MIN_REAL                -HO_MAX_REAL
#define HO_MACHINE_EPSILON         2.2204460492503131e-16
#define HO_SQRT_MACHINE_EPSILON    1.490116119384766e-08
#define HO_SMALL                   1.e-10
#define HO_TOLERANCE               1.e-6
#elif defined(HO_USE_REAL___FLOAT128)
#define HO_MAX_REAL                FLT128_MAX
#define HO_MIN_REAL                -FLT128_MAX
#define HO_MACHINE_EPSILON         FLT128_EPSILON
#define HO_SQRT_MACHINE_EPSILON    1.38777878078e-17
#define HO_SMALL                   1.e-20
#define HO_TOLERANCE               1.e-8
#endif

#define HO_INFINITY                HO_MAX_REAL/4.0
#define HO_NINFINITY              -HO_INFINITY

#define HO_GAMA                    1.4
#define HO_GAMAM1                  0.4

//! ----------------------------------------------------------------------------*/
//!
//!   Basic constants, copy from MOOSE
//!
/**
 * \enum ::Order defines an \p enum for polynomial orders.
 * Fixing each label to a specific int, since \p InfFE and p refinement
 * may cast between them
 */
enum Order
{
    CONSTANT     =  0,
    FIRST        =  1,
    SECOND       =  2,
    THIRD        =  3,
    FOURTH       =  4,
    FIFTH        =  5,
    SIXTH        =  6,
    SEVENTH      =  7,
    EIGHTH       =  8,
    NINTH        =  9,
    TENTH        = 10,

    ELEVENTH     = 11,
    TWELFTH      = 12,
    THIRTEENTH   = 13,
    FOURTEENTH   = 14,
    FIFTEENTH    = 15,
    SIXTEENTH    = 16,
    SEVENTEENTH  = 17,
    EIGHTTEENTH  = 18,
    NINTEENTH    = 19, // deprecated
    NINETEENTH   = 19,
    TWENTIETH    = 20,

    TWENTYFIRST   = 21,
    TWENTYSECOND  = 22,
    TWENTYTHIRD   = 23,
    TWENTYFOURTH  = 24,
    TWENTYFIFTH   = 25,
    TWENTYSIXTH   = 26,
    TWENTYSEVENTH = 27,
    TWENTYEIGHTH  = 28,
    TWENTYNINTH   = 29,
    THIRTIETH     = 30,

    THIRTYFIRST   = 31,
    THIRTYSECOND  = 32,
    THIRTYTHIRD   = 33,
    THIRTYFOURTH  = 34,
    THIRTYFIFTH   = 35,
    THIRTYSIXTH   = 36,
    THIRTYSEVENTH = 37,
    THIRTYEIGHTH  = 38,
    THIRTYNINTH   = 39,
    FORTIETH      = 40,

    FORTYFIRST   = 41,
    FORTYSECOND  = 42,
    FORTYTHIRD   = 43,

    INVALID_ORDER
};

enum 
{
    HO_STEADY   =  0,
    HO_UNSTEADY =  1
};

// ------------------------------------------------------------
// enum FEFamily definition

/**
 * \enum ::FEFamily defines an \p enum for finite element families.
 */
// vanilla C0
enum FEFamily
{
    LAGRANGE     = 0,
    HIERARCHIC   = 1,

    // discontinuous, in local coordinates
    MONOMIAL      = 2,
    L2_HIERARCHIC = 6,
    L2_LAGRANGE   = 7,

    // higher-order
    BERNSTEIN    = 3,
    SZABAB       = 4,

    // discontinuous, in global coordinates
    XYZ          = 5,

    // infinite element stuff
    INFINITE_MAP = 11,     //   for 1/r-map
    JACOBI_20_00 = 12,     //   i_max = 19
    JACOBI_30_00 = 13,     //   i_max = 19
    LEGENDRE     = 14,     //   i_max = 19

    // C1 elements
    CLOUGH       = 21,
    HERMITE      = 22,
    SUBDIVISION  = 23,

    // A scalar variable that couples to
    // all other DOFs in the system
    SCALAR       = 31,

    // Vector-valued elements
    LAGRANGE_VEC = 41,
    NEDELEC_ONE  = 42,

    INVALID_FE   = 99
};

/**
 * \enum ::FEContinuity defines an \p enum for finite element
 * types to _assert a certain level (or type? Hcurl?) of continuity.
 */
enum FEContinuity {DISCONTINUOUS,
                   C_ZERO,
                   C_ONE,
                   H_CURL};

/**
 * \enum ::FEFieldType defines an \p enum for finite element
 * field types - i.e. is it a scalar element, vector, tensor, etc.
 */
enum FEFieldType {TYPE_SCALAR = 0,
                  TYPE_VECTOR};


/*!
 * \brief types of geometric entities based on VTK nomenclature
 */
enum GEO_TYPE
{
    TRIANGLE = 5,       /*!< \brief VTK nomenclature for defining a triangle element. */
    QUADRILATERAL = 8,      /*!< \brief VTK nomenclature for defining a quadrilateral element. */
    TETRAHEDRON = 11,       /*!< \brief VTK nomenclature for defining a tetrahedron element. */
    HEXAHEDRON = 14,        /*!< \brief VTK nomenclature for defining a hexahedron element. */
    PRISM = 17,             /*!< \brief VTK nomenclature for defining a prism element. */
    PYRAMID = 20        /*!< \brief VTK nomenclature for defining a pyramid element. */
};

//enum GEO_TYPE {
//  VERTEX = 1,           /*!< \brief VTK nomenclature for defining a vertex element. */
//  LINE = 3,            /*!< \brief VTK nomenclature for defining a line element. */
//  TRIANGLE = 5,         /*!< \brief VTK nomenclature for defining a triangle element. */
//  QUADRILATERAL = 9,        /*!< \brief VTK nomenclature for defining a quadrilateral element. */
//  TETRAHEDRON = 10,         /*!< \brief VTK nomenclature for defining a tetrahedron element. */
//  HEXAHEDRON = 12,          /*!< \brief VTK nomenclature for defining a hexahedron element. */
//  PRISM = 13,             /*!< \brief VTK nomenclature for defining a prism element. */
//  PYRAMID = 14          /*!< \brief VTK nomenclature for defining a pyramid element. */
//};

template < typename T >
T ComputePressure(T * conserv)
{
    return HO_GAMAM1 * (conserv[4] - 0.5 * (pow(conserv[1],2) + pow(conserv[2],2) + pow(conserv[3],2))/conserv[0]);
}

class HOMemoryHandler
{
public:
    HOMemoryHandler()
    {

    }

private:
    int d_level,*d_orderP;  //level:p-multigrid,order:orderP,
    RDouble *** u_h, *** w_h; //u_h[iCell][iEqu][iDof], w_h[0] -> high order rho 's dof for DG/FV

    RDouble * massMatrix;
    RDouble * invMassMatrix; // dof X dof //
    RDouble * reconstruMatrix; //DG/FV reconstruction massmatrix

    int * nGaussPtsOfCell;  //number of gauss points for each Cell,nGaussPtsOfCell[nTCell]

    RDouble (* gaussPtsOfCell)[3]; //volume gauss point coordinate, gaussPtsOfCell[Cell][0-2]
    RDouble *** basisGaussPtsOfCell; //volume gauss point coordinate, basisGaussPtsOfCell[iCell][iGauss][iDof]
    RDouble *** basisGradGaussPtsOfCell; //volume gauss point coordinate, basisGradGaussPtsOfCell[iCell][iGauss][3*iDof]

    int * nGaussPtsOfFace;   //number of gauss points for each Face,nGaussPtsOfFace[nTFace]

    RDouble (* gaussPtsOfFace)[3]; //face gauss points coordinate xyz, gaussPtsOfFace[iFace][0-2]
    RDouble *** basisGaussPtsOfFace; //face gauss point coordinate basis, basisGaussPtsOfFace[2*iFace][iGauss][iDof]
    RDouble *** basisGradGaussPtsOfFace; //face gauss point coordinate, basisGradGaussPtsOfFace[2*iFace][iGauss][3*iDof]

};

}