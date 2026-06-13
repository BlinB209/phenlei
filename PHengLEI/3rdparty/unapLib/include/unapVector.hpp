/// \file unapVector.hpp
/// \author Hanfeng GU
/// \version 1.0
/// \date 2019-01-30
/// \brief base data structure in UNAP

#ifndef UNAPVECTOR_HPP
#define UNAPVECTOR_HPP

#include <string.h>

#include "unap.hpp"
#include "unapMPI.hpp"

namespace UNAP {
/// \brief store array in continuous
template <typename T>
class Vector {
 private:
  /// \brief size of vector
  label length_;
  /// \brief data pointer of vector
  T *values_;
  /// \brief communication domain
  Communicator *commcator_;
  /// \brief don't delete when destructing if true
  bool dontDel_;

 public:
  /*! * * * * * * * * * * * * assign communicator * * * * * * * * * * */
  // vector default value is 0

  /// \brief constructor
  /// \param comm communication domain
  Vector(Communicator *comm);
  /// \brief constructor
  /// \param length size of vector
  /// \param comm communication domain
  Vector(const label length, Communicator *comm);
  /// \brief copy from an existing array with a giving length
  /// \param val data pointer
  /// \param length size of vector
  /// \param comm communication domain
  /// \param reUse if the data will be reused
  /// \param dontDel won't delete pointer if true
  Vector(const T *val, const label &length, Communicator *comm,
         const bool reUse = false, const int dontDel = 0);
  /// \brief build a vector with a given length and same value
  /// \param length size of vector
  /// \param value data in vector is assigned to the same value
  /// \param comm communication domain
  Vector(const label length, const T value, Communicator *comm);

  /*! * * * * * * * * * * * * use global communicator * * * * * * * * * * */
  /// \brief constructor
  Vector();
  /// \brief constructor
  /// \param length size of vector
  Vector(const label length);
  /// \brief copy from an existing array with a giving length
  /// \param val data pointer
  /// \param length size of vector
  /// \param reUse if the data will be reused
  /// \param dontDel won't delete pointer if true
  Vector(const T *val, const label &length, const bool reUse = false,
         const int dontDel = 0);
  /// \brief build a vector with a given length and same value
  /// \param length size of vector
  /// \param value data in vector is assigned to the same value
  Vector(const label length, const T value);

  /*! * * * * * * * * * * * * others * * * * * * * * * * */
  /// \brief copy from an existing vector
  /// \param v an existing vector
  Vector(const Vector<T> &v);
  /// \brief copy from an existing stl vector
  /// \param v an existing stl vector
  Vector(const Array<T> &v);

  /// \brief destructor
  virtual ~Vector();

  /// \brief set communicator
  /// \param comm communication domain
  void setCommunicator(Communicator *comm) { commcator_ = comm; }

  /// \brief get communicator
  Communicator *getCommunicator() const { return commcator_; }

  /// \brief return value of giving index
  /// \param i giving index
  inline T &operator[](const label i) { return (this->values_[i]); }

  /// \brief return constant value of giving index
  /// \param i giving index
  inline const T operator[](const label i) const { return (this->values_[i]); }

  /// \brief copy from an existing vector
  /// \param v an existing vector
  Vector &operator=(const Vector<T> &v);

  /// \brief set all values to the scalar a
  Vector &operator=(const T &a);

  /// \brief both vectors must be of equal length
  Vector &operator+=(const Vector<T> &v);
  Vector &operator-=(const Vector<T> &v);

  Vector &operator+=(const T &a);
  Vector &operator-=(const T &a);
  Vector &operator*=(const T &a);
  Vector &operator/=(const T &a);

  const Vector operator+(const Vector<T> &v) const;
  const Vector operator-(const Vector<T> &v) const;
  const Vector operator*(const T &a) const;

  /// \brief return size of vector
  label size() const { return (this->length_); }

  /// \brief set size of vector
  void SET_size(const label newSize);

  /// \brief scale the  vector
  void scale(T a);

  /// \brief return data pointer of vector
  T *data() const { return (this->values_); }

  /// \brief sum of all values in vector
  /// sum = v1 + v2 + v3 + ...
  T Sum() const;

  /// \brief absolute sum of all values in vector
  /// sum = abs(v1) + abs(v2) + abs(v3) + ...
  T SumMag() const;

  /// \brief square sum of all values in vector
  /// sum = v1*v1 + v2*v2 + v3*v3 + ...
  T SumSqr() const;

  /// \brief squart sum of all values in vector
  /// sum = (v1*v1 + v2*v2 + v3*v3 + ...)^0.5
  scalar SumSqrt() const;

  /// \brief set all values to be zero in vector
  void SET_zero();
};

/*! * * * * * * * * * * * * * * *  * * * * * * * * * * * * * */
/// \brief y=1/y
template <typename T>
void reciprocal(Vector<T> &y);
/// \brief y=1/x
template <typename T>
void reciprocal(Vector<T> &y, const Vector<T> &x);
/// \brief res=v1*v2
template <typename T>
T dot(const Vector<T> &v1, const Vector<T> &v2);
/// \brief y=a*x+b*y
template <typename T>
void axpby(Vector<T> &y, T a, T b, const Vector<T> &x);
/// \brief y=y+a*x
template <typename T>
void axpy(Vector<T> &y, T a, const Vector<T> &x);
/// \brief y=a*y+x
template <typename T>
void aypx(Vector<T> &y, T a, const Vector<T> &x);
/// \brief w=ax+by
template <typename T>
void waxpy(Vector<T> &w, T a, const Vector<T> &x, const Vector<T> &y);
/// \brief y=ax+by+z
template <typename T>
void axpbypz(Vector<T> &y, T a, T b, const Vector<T> &x, const Vector<T> &z);
/// \brief wi=xi*yi
template <typename T>
void pointwiseMul(Vector<T> &w, const Vector<T> &x, const Vector<T> &y);
/// \brief wi=xi/yi
template <typename T>
void pointwiseDiv(Vector<T> &w, const Vector<T> &x, const Vector<T> &y);

typedef Vector<label> labelVector;
typedef Vector<scalar> scalarVector;
}  // namespace UNAP
#endif  // end UNAPVECTOR_HPP
