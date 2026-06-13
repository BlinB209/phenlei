#ifndef LDUDILUPRECOND_HPP
#define LDUDILUPRECOND_HPP

#include "unapMatrix.hpp"
#include "unapPreconditioner.hpp"

namespace UNAP {
class LduMatrix;

/// \brief Simplified Diagonal-based Incomplete LU preconditioner
/// support symmetric and asymmetric, same as DIC if matrix is symmetric
class LduDILUPrecond : public Preconditioner {
 private:
  /// \brief the reciprocal diagonal
  scalarVector rD_;

  /// \brief LduMatrix
  const LduMatrix *APtr_;

 public:
  /// \brief constructor
  /// \param A matrix input
  LduDILUPrecond(const LduMatrix &A);

  /// \brief destructor
  virtual ~LduDILUPrecond() {}

  /// \brief calculate the reciprocal of the preconditioned diagonal
  /// \param rD the resulting reciprocal diagonal
  /// \param A matrix input
  static void calcReciprocalD(scalarVector &rD, const LduMatrix &A);

  /// \brief return wA the preconditioned form of residual rA
  /// \param w resulting vector, preconditioned by rA
  /// \param r Preconditioner vector
  virtual void precondition(scalarVector &w, const scalarVector &r) const;

  /// \brief return the reciprocal diagonal
  virtual const scalarVector &rD() const { return rD_; }
};

}  // namespace UNAP

#endif  // LDUDILUPRECOND_HPP
