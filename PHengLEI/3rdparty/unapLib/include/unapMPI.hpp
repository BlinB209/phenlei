/// \file unapMPI.hpp
/// \author Hanfeng GU
/// \version 1.0
/// \date 2019-01-30
/// \brief function reduceSum

#ifndef UNAPMPI_HPP
#define UNAPMPI_HPP

#include <mpi.h>

#include "unap.hpp"

namespace UNAP {

/// \brief 数组和规约，提供数组每个位置上的和规约
/// \param T 数组类型
/// \param v 数组首地址
/// \param commcator 规约操作通信域指针
/// \param n 数组大小，默认大小为1
template <typename T>
void reduceSum(T *v, Communicator *commcator, label n = 1) {
  if (commcator->getMySize() > 1) {
    T *vLocal = new T[n];
    memcpy(vLocal, v, n * sizeof(T));

    CommData myType;
    if (typeid(T) == typeid(label32)) {
      myType = COMM_INT;
    } else if (typeid(T) == typeid(label64)) {
      myType = COMM_LONG;
    } else if (typeid(T) == typeid(scalar32)) {
      myType = COMM_FLOAT;
    } else if (typeid(T) == typeid(scalar64)) {
      myType = COMM_DOUBLE;
    }

    commcator->allReduce("sum", &vLocal[0], v, n, myType, COMM_SUM);
    commcator->finishTask("sum");
    DELETE_POINTERS(vLocal);
  }
}

}  // namespace UNAP

#endif  //- UNAPMPI_HPP
