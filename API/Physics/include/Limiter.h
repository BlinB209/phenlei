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
//! @file      Limiter.h
//! @brief     Spatial inviscid flux limiter.
//! @author    Bell.

#pragma once
#include "Geo_UnstructGrid.h"
using namespace std;
namespace PHSPACE
{
template < typename T >
class Data2D;
class Gradient;
//! @brief Limiter class defines the method of limiter function during inviscid flux calculation.\n
//! calculate limiter with inter faces communication.
class Limiter
{
public:
    typedef RDouble * iterator;

private:
    //! The unstructured grid that the variables load on.
    UnstructGrid * grid;

    //! The variables that need to calculate limiter.
    RDouble ** q;

    //! Number of variables, which is the above 'q', that need to calculate limiter.
    int nLimiterVariables;

    //! Gradient of the variables, which is the above 'q', that need to calculate limiter.
    //Gradient *gradient;

    //! Two-dimensional limiter pointer.\n
    //! 1th dim point to the variable limiter.\n
    //! 2nd dim point to the limiter value on each cell.
    Data2D<RDouble> * limit;

    //! Limiter variable name, used to identity interface limiter when parallel communication.
    string limiterVarName;

    //! Limiter type.
    //! - # ILMT_NOLIM (1): Without limiter.
    //! - # ILMT_FIRST (2): First order.
    //! - # ILMT_BARTH (3): Barth limiter, for unstructured solver.
    //! - # ILMT_VENCAT (4): Vencat limiter, for unstructured solver.
    //! - # ILMT_STRUCT (10):
    //! - # ILMT_MINMOD (11): MINMOD limiter, for structured solver.
    //! - # ILMT_VAN_ALBADA (12): Van Albada limiter, for structured solver.
    //! - # ILMT_VAN_ALBADA_CLZ (13):
    //! - # ILMT_VANLEER (14): VanLeer limiter, for structured solver.
    //! - # ILMT_SMOOTH (15):
    //! - # ILMT_MINMODTEST (16):
    //! - # ILMT_MIN_VAN (17):
    int limiterTypeID;

    //! Limiter model:\n
    //! -# 0: limit only for pressure and density, then get the min value\n
    //! -# 1: limit for every variables, then get the min value
    int limitModel;

    //! Scalar or vector limiter:\n
    //! -# 0: Each variable use the same limiter coefficient.\n
    //! -# 1: Each variable use the respective limiter coefficients.
    int limitVector;

    //! Vencat limiter type.\n
    //! -# 0: The modification of Dr. Wang, AIAA-96-2091. Independent of grid scale.
    //!       venkatCoeff = [0.01, 0.2].
    //! -# 2: This limiter type was proposed by CXH. Which has nothing to do with grid scale, mach number.
    //!       But it need to be V&V in more numerical experiments.
    //! -# 3: This limiter type was modified by NeverMore.
    //!       Which the limiter coefficient of eps^2 is influenced by cell volume,
    //!       other than distance between cell center and face center, that was used in VencatLimiter4.
    //! -# 4: Original format vencat limiter, depend on grid scale.
    int vencatLimiterType;

public:
    //! @param[in] gridin     The unstructured grid that the variables load on.\n
    //! @param[in] q          The variables that need to calculate limiter.\n
    //! @param[in] gradient   Gradient of the variables, which is the above 'q', 
    //!                       that need to calculate limiter.\n
    //! @param[in] nVariable  Number of variables, which is the above 'q', 
    //!                       that need to calculate limiter.\n
    //! @param[in] limiterVar Name Limiter variable name, 
    //!                       used to identity interface limiter when parallel communication.
    LIB_EXPORT Limiter(UnstructGrid *gridin, RDouble **q, int nVariable, const string &limiterVarName);
    LIB_EXPORT ~Limiter();

public:
    //! Set limiter type.\n
    //! - # ILMT_NOLIM (1): Without limiter.\n
    //! - # ILMT_FIRST (2): First order.\n
    //! - # ILMT_BARTH (3): Barth limiter, for unstructured solver.\n
    //! - # ILMT_VENCAT (4): Vencat limiter, for unstructured solver.\n
    //! - # ILMT_STRUCT (10):\n
    //! - # ILMT_MINMOD (11): MINMOD limiter, for structured solver.\n
    //! - # ILMT_VAN_ALBADA (12): Van Albada limiter, for structured solver.\n
    //! - # ILMT_VAN_ALBADA_CLZ (13):\n
    //! - # ILMT_VANLEER (14): VanLeer limiter, for structured solver.\n
    //! - # ILMT_SMOOTH (15):\n
    //! - # ILMT_MINMODTEST (16):\n
    //! - # ILMT_MIN_VAN (17):\n
    LIB_EXPORT void SetLimiterTypeID(int uns_limiter); 

    //! Set limiter model:\n
    //! -# 0: limit only for pressure and density, then get the min value\n
    //! -# 1: limit for every variables, then get the min value
    LIB_EXPORT void SetLimiterModel(int limit_mode);   

    //! Set scalar or vector limiter:\n
    //! -# 0: Each variable use the same limiter coefficient.\n
    //! -# 1: Each variable use the respective limiter coefficients.
    LIB_EXPORT void SetLimiterVector(int limit_vec);

    //! Set Vencat limiter type.\n
    //! -# 0: The modification of Dr. Wang, AIAA-96-2091. Independent of grid scale.
    //!       venkatCoeff = [0.01, 0.2].
    //! -# 2: This limiter type was proposed by CXH. Which has nothing to do with grid scale, mach number.
    //!       But it need to be V&V in more numerical experiments.\n
    //! -# 3: This limiter type was modified by NeverMore.
    //!       Which the limiter coefficient of eps^2 is influenced by cell volume,
    //!       other than distance between cell center and face center, that was used in VencatLimiter4.\n
    //! -# 4: Original format vencat limiter, depend on grid scale.
    LIB_EXPORT void SetVencatLimiterType(int VenkatLimiterType);

    //! Return limiter of ith-variable.
    LIB_EXPORT iterator GetLimiter(int m = 0);

    //! Run the limiter calculation.
    LIB_EXPORT void Calculation(Grid *grid);

    //! Initialzie the Two-dimensional limiter pointer(limit)'s data
    LIB_EXPORT void InitLimitData();

public:
    void CreateLimiter(int nsize, int neqn = 1);

    void BarthLimiter  (RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void VencatLimiter (RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int iVariable = 0);
    void VencatLimiter0(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void VencatLimiter1(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void VencatLimiter2(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int iVariable = 0);
    void VencatLimiter3(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int iVariable = 0);
    void VencatLimiter4(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);    // Bell 20121028 add
    void VencatLimiter5(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void VencatLimiter6(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void VencatLimiter7(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int iVariable = 0);

    void MinMaxDiff(RDouble *dmin, RDouble *dmax, RDouble *q);
};

#include "Limiter.hxx"

template < typename T >
class Data2D
{
private:
    T *data;
    int nsize, neqn, ntotal_size;
public:
    //! type definitions
    typedef T  value_type;
    typedef T* iterator;
public:
    Data2D()
    {
        data = 0;
    }
    ~Data2D()
    {
        delete [] data;
    }
public:
    void CreateData(int nsize, int neqn = 1)
    {
        this->nsize = nsize;
        this->neqn  = neqn;
        ntotal_size = nsize * neqn;
        data = new T [ ntotal_size ];
        std::fill_n(data, ntotal_size, 0);
    }
    iterator GetData(int m)
    {
        return &data[ m * nsize ];
    }

    template < typename T2 >
    void SetData(const T2 & value)
    {
        T t_value = value;
        std::fill_n(data, ntotal_size, t_value);
    }
};

void Init_Limiters(Grid *grid_in, RDouble *limit);

int GetLimiterID(const string &limiter_name);

}
