/*! \file mpiWrapper.hpp
 *  \author Zhao Chengpeng (chengpeng_zhao@foxmail.com)
 *  \date 2022-02-21
 *  \modified 2022-04-18
 *  \brief utilities不能满足的mpi相关功能在此实现
 */
#ifndef TOOLS_WRAPPER_MPIWRAPPER_HPP
#define TOOLS_WRAPPER_MPIWRAPPER_HPP
#include <mpi.h>

#include <complex>
#include <limits>
#include <numeric>
#include <vector>

#if defined(_WIN32)
#define NOMINMAX
#endif
#include "utilities.h"

namespace UNAP {

/*!
 * Return the corresponding MPI_Datatype for simple C++ data types.
 * \tparam T C++ type for which to return the corresponding MPI_Datatype
 */
template <typename T>
MPI_Datatype mpiType() {
  return T::mpiType();
}
/*! return MPI datatype for C++ char */
template <>
inline MPI_Datatype mpiType<char>() {
  return MPI_CHAR;
}
/*! return MPI datatype for C++ bool */
// template <>
// inline MPI_Datatype mpiType<bool>() {
//  return MPI_CXX_BOOL;
//}
/*! return MPI datatype for C++ int */
template <>
inline MPI_Datatype mpiType<int>() {
  return MPI_INT;
}
/*! return MPI datatype for C++ long */
template <>
inline MPI_Datatype mpiType<long>() {
  return MPI_LONG;
}
/*! return MPI datatype for C++ unsigned long */
template <>
inline MPI_Datatype mpiType<unsigned long>() {
  return MPI_UNSIGNED_LONG;
}
/*! return MPI datatype for C++ long long int */
template <>
inline MPI_Datatype mpiType<long long int>() {
  return MPI_LONG_LONG_INT;
}
/*! return MPI datatype for C++ float */
template <>
inline MPI_Datatype mpiType<float>() {
  return MPI_FLOAT;
}
/*! return MPI datatype for C++ double */
template <>
inline MPI_Datatype mpiType<double>() {
  return MPI_DOUBLE;
}
/*! return MPI datatype for C++ std::complex<float> */
template <>
inline MPI_Datatype mpiType<std::complex<float> >() {
#ifndef MPI_CXX_FLOAT_COMPLEX
#define MPI_CXX_FLOAT_COMPLEX MPI_COMPLEX
#endif
  return MPI_CXX_FLOAT_COMPLEX;
}
/*! return MPI datatype for C++ std::complex<double> */
template <>
inline MPI_Datatype mpiType<std::complex<double> >() {
#ifndef MPI_CXX_DOUBLE_COMPLEX
#define MPI_CXX_DOUBLE_COMPLEX MPI_DOUBLE_COMPLEX
#endif
  return MPI_CXX_DOUBLE_COMPLEX;
}
/*! return MPI datatype for C++ std::pair<int,int> */
template <>
inline MPI_Datatype mpiType<std::pair<int, int> >() {
  return MPI_2INT;
}
/*! return MPI datatype for C++ std::pair<long int,long int> */
template <>
inline MPI_Datatype mpiType<std::pair<long int, long int> >() {
  static MPI_Datatype llMpiType = MPI_DATATYPE_NULL;
  if (llMpiType == MPI_DATATYPE_NULL) {
    /*! creates an MPI datatype by replicating an existing one a certain number
     * of times. */
    MPI_Type_contiguous(2, mpiType<long int>(), &llMpiType);
    MPI_Type_commit(&llMpiType);
  }
  return llMpiType;
}

/*!
 * Perform an MPI_Alltoallv. Each rank sends sbuf[i] to process i
 * The results are received in a single contiguous vector rbuf
 * pbuf has pointers into rbuf, with pbuf[i] pointing to the
 * data received from rank i.
 * \tparam T type of data to send
 * \param sbuf send buffers (should be size P)
 * \param rbuf receive buffer, can be empty, will be allocated
 * \param pbuf pointers (to positions in rbuf) to where data
 * received from different ranks start
 * \param Ttype MPI_Datatype corresponding to the template parameter T
 */
template <typename T>
void allToallv(MPI_Comm mpiComm, std::vector<std::vector<T> >& sbuf,
               std::vector<T>& rbuf, std::vector<T*>& pbuf,
               const MPI_Datatype Ttype);

/*!
 * Perform an MPI_Alltoallv. Each rank sends sbuf[i] to process
 * i. The results are received in a single contiguous vector which
 * is returned.
 * \tparam T type of data to send, this should have a
 * corresponding mpiType<T>() implementation or should define
 * T::mpiType()
 * \param sbuf send buffers (should be size P)
 * \rbuf return, receive buffer
 */
template <typename T>
void allToallv(std::vector<T>& rbuf, MPI_Comm comm,
               std::vector<std::vector<T> >& sbuf) {
  std::vector<T*> pbuf;
  allToallv(comm, sbuf, rbuf, pbuf, mpiType<T>());
}

/*!
 * Receive a vector of T's from process src, with tag. The message
 * size does not need to be known in advance.
 * \tparam T template parameter of vector to receive, should have
 * a corresponding mpiType<T>() implementation
 * \param src process to receive from
 * \param tag tag to match the message
 * \param rbuf return, std::vector<T> with the data to be received.
 */
template <typename T>
void recv(std::vector<T>& rbuf, label src, label tag, MPI_Comm comm_) {
  MPI_Status stat;
  MPI_Probe(src, tag, comm_, &stat);
  label msgsize;
  MPI_Get_count(&stat, mpiType<T>(), &msgsize);
  rbuf.resize(msgsize);
  MPI_Recv(rbuf.data(), msgsize, mpiType<T>(), src, tag, comm_,
           MPI_STATUS_IGNORE);
}
/*!
 * Receive a vector of T's from any processes, with tag. The message
 * size does not need to be known in advance.
 * \tparam T template parameter of vector to receive, should have
 * a corresponding mpiType<T>() implementation
 * \param tag tag to match the message
 * \param rPair a pair that contains the rank of the process that sent the
 * message and a rvalue reference of std::vector<T> with the data to be
 * received.
 */
template <typename T>
void recvAnySrc(std::pair<label, std::vector<T> >& rPair, label tag,
                MPI_Comm comm_) {
  MPI_Status stat;
  /*! Tell MPI to receive msg without restricting the rank of the sender. */
  MPI_Probe(MPI_ANY_SOURCE, tag, comm_, &stat);
  label msgsize;
  MPI_Get_count(&stat, mpiType<T>(), &msgsize);
  rPair.second.resize(msgsize);
  MPI_Recv(rPair.second.data(), msgsize, mpiType<T>(), stat.MPI_SOURCE, tag,
           comm_, MPI_STATUS_IGNORE);
  rPair.first = stat.MPI_SOURCE;
}
/*!
 * Receive a value of T from any process src, with tag
 * \tparam T template parameter of vector to receive, should have
 * a corresponding mpiType<T>() implementation
 * \param tag tag to match the message
 */
template <typename T>
T recvOne(label src, label tag, MPI_Comm comm_) {
  T t;
  MPI_Recv(&t, 1, mpiType<T>(), src, tag, comm_, MPI_STATUS_IGNORE);
  return t;
}

/*!
 * Return a subcommunicator with P ranks, starting from rank P0
 * [P0:stride:P0+stride*P) will be included in the new communicator
 * This operation is collective on all the processes in this communicator
 * \param P0 first rank in the new communicator
 * \param P number of ranks in the new communicator
 * \param stride stride between ranks in this communicator
 * determining which ranks will go into new communicator
 * \return new communicator containing P ranks from this
 * communicator [P0:stride:P0+stride*P)
 */
void mpiSubComm(MPI_Comm& newComm, const MPI_Comm& comm, label P0, label P,
                label stride = 1);
/*!
 * \class Triplet
 * \brief helper class to communicate matrix data
 */
template <typename Cmpt>
class Triplet {
 public:
  label r_, c_;
  Cmpt v_;
  Triplet() {}
  Triplet(label row, label col, Cmpt value) : r_(row), c_(col), v_(value) {}
  static MPI_Datatype tripletMpiType_;
  static MPI_Datatype mpiType();
  static void freeMpiType();
};

/*!
 * \class IdxIJ
 * \brief helper class to get symm pattern matrix after pivoting
 */
class IdxIJ {
 public:
  label i_, j_;
  IdxIJ() {}
  IdxIJ(label ii, label jj) : i_(ii), j_(jj) {}
  static MPI_Datatype idxijMpiType_;
  static MPI_Datatype mpiType();
  static void freeMpiType();
};
/*!
 * \class IdxVal
 * \brief helper class to permute solved result when using pivoting
 */
template <typename Cmpt>
class IdxVal {
 public:
  label i_;
  Cmpt v_;
  IdxVal() {}
  IdxVal(label idx, Cmpt value) : i_(idx), v_(value) {}
  static MPI_Datatype idxvalMpiType_;
  static MPI_Datatype mpiType();
  static void freeMpiType();
};

} /* End namespace UNAP */

#endif /* end of include guard: TOOLS_WRAPPER_MPIWRAPPER_HPP */
