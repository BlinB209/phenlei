#ifndef UNAPMATRIX_HPP
#define UNAPMATRIX_HPP

#include "interfaces.hpp"
#include "unapVector.hpp"

#ifdef SW_SLAVE
// #include "iterator.h"
// #include "swMacro.h"
#endif

namespace UNAP {
/// \brief matrix base class
class Matrix {
 protected:
  Communicator *commcator_;

 public:
  /// \brief constructor
  /// \param other_comm communication domain
  Matrix(Communicator *other_comm) : commcator_(other_comm) {}

  /// \brief destructor
  virtual ~Matrix() {}

  /// \brief sparse matrix-vector multiplication: y=Ax
  /// \param Apsi the result vector y
  /// \param psi vector x
  virtual void spMV(scalarVector &Apsi, const scalarVector &psi) const = 0;

  /// \brief number of equations( = number of rows)
  virtual label size() const = 0;

  /// \brief access to diagonal data
  virtual scalarVector &diag() const = 0;

  /// \brief return if matrix is symmetric or asymmetric
  virtual bool symm() const = 0;

  /// \brief access to interfaces
  virtual Interfaces &matrixInterfaces() const = 0;

  /// \brief initialize interfaces
  virtual void initInterfaces(const scalarVector &psi) const = 0;

  /// \brief update interfaces
  virtual void updateInterfaces(scalarVector &Apsi) const = 0;

  /// \brief set communicator
  void setCommunicator(Communicator *other_comm) { commcator_ = other_comm; }

  /// \brief communicator
  Communicator *getCommunicator() const { return commcator_; }

  /// \brief return local nonzero number
  virtual label lnnz() const = 0;
};

}  // namespace UNAP
#endif  //- UNAPMATRIX_HPP
