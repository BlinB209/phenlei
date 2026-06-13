/// \file unap.hpp
/// \author Hanfeng GU
/// \version 1.0
/// \date 2019-01-30
/// \modified 2022-04-15
/// \brief 基础头文件，从这里引入utilities的头文件

#ifndef UNAP_HPP
#define UNAP_HPP

#if defined(_WIN32)
#define NOMINMAX
#endif

#include <cmath>
#include <limits>

#include "mpiWrapper.hpp"
#include "utilities.h"

namespace UNAP  ///< 命名空间 UNAP，本模块所有类均在UNAP命名空间下
{
/// \brief sleep for the specified number of seconds
unsigned int sleep(const unsigned int);

/// \brief small scalar for the use in solvers
// #define SMALL (1.0e-37)
// #define VSMALL (1.0e-128)
#define SMALL (std::numeric_limits<scalar>::epsilon())
#define VSMALL (std::numeric_limits<scalar>::min())

/// \brief large scalar for the use in solvers
#define GREAT (1.0 / std::numeric_limits<scalar>::epsilon())
#define VGREAT (std::numeric_limits<scalar>::max())

/// \brief check if pointer exists
#define CHECK_POINTER(ptr)                                             \
  if (!ptr) {                                                          \
    POUT << "ERROR in " << __FILE__ << " " << __LINE__ << ": " << #ptr \
         << " is NULL!" << ENDL;                                       \
    ERROR_EXIT;                                                        \
  }

#define ALLOCATE_POINTER(ptr, oldObj, T)                  \
  {                                                       \
    DELETE_POINTER(ptr)                                   \
    ptr = new T(oldObj.size(), oldObj.getCommunicator()); \
    T &newObj = *ptr;                                     \
    newObj = oldObj;                                      \
  }

/// \brief 返回绝对值
/// 名称比较普通，可能会与其他库同名
#define ABS(v) std::abs(v)

#ifdef SW_SLAVE

/// \brief 判断矩阵的行数rowSize小于神威加速临界值accUsingSize
/// accUsingSize目前定义在swArrays vectorOps_struct.h中
#define IFNOT_SWACC(size) if (size < accUsingSize)

/// \brief 判断矩阵的行数rowSize小于神威加速临界值2500，一般用于操作数较多的计算
#define IFNOT_SWACC_SMALL(size) if (size < 2500)

#else
#define IFNOT_SWACC(size)

#define IFNOT_SWACC_SMALL(size)
#endif

/*! partitioning size for 2D block-cyclic process grid */
#define PROCESSGRID_DEFAULT_BLOCKSIZE 64

/*! 判断是否为复数 */
template <typename Cmpt>
inline bool isComplex() {
  return false;
}
template <>
inline bool isComplex<std::complex<scalar> >() {
  return true;
}
/*! 实数时返回实数，复数时返回复数 */
template <typename Cmpt>
Cmpt getScalar(scalar vr, scalar vi) {
  return Cmpt(vr);
}
template <>
inline std::complex<scalar> getScalar(scalar vr, scalar vi) {
  return std::complex<scalar>(vr, vi);
}
/*! 返回某数的共轭 */
inline scalar myConj(scalar a) { return a; }
inline std::complex<scalar> myConj(std::complex<scalar> a) {
  return std::conj(a);
}

/*! 返回实部 */
template <typename Cmpt>
inline scalar real(Cmpt a) {
  return a;
}
template <>
inline scalar real(std::complex<scalar> a) {
  return a.real();
}

}  // namespace UNAP

#endif
