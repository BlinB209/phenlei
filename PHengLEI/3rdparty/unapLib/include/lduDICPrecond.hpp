#ifndef LDUDICPRECOND_HPP
#define LDUDICPRECOND_HPP

#include "unapMatrix.hpp"
#include "unapPreconditioner.hpp"

namespace UNAP {
class LduMatrix;

/// \brief Simplified Diagonal-based Incomplete Cholesky preconditioner
/// used if matrix is symmetric
class LduDICPrecond : public Preconditioner {
 private:
  /// \brief the reciprocal diagonal
  scalarVector rD_;

  /// \brief LduMatrix
  const LduMatrix *APtr_;

 public:
  /// \brief constructor
  /// \param A matrix input
  LduDICPrecond(const LduMatrix &A);

  /// \brief destructor
  virtual ~LduDICPrecond() {}

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

#endif  // LDUDICPRECOND_HPP
