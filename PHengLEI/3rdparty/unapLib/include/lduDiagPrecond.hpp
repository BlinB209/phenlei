#ifndef LDUDIAGPRECOND_HPP
#define LDUDIAGPRECOND_HPP

#include "unapMatrix.hpp"
#include "unapPreconditioner.hpp"

namespace UNAP {
class LduMatrix;

/// \brief diagonal Preconditioner in LDU type
/// also seen as Jacobi Preconditioner
class LduDiagPrecond : public Preconditioner {
 private:
  /// \brief the reciprocal diagonal
  scalarVector rD_;

 public:
  /// \brief constructor
  /// \param A matrix input
  LduDiagPrecond(const LduMatrix &A);

  /// \brief destructor
  virtual ~LduDiagPrecond() {}

  /// \brief return wA the preconditioned form of residual rA
  /// \param w resulting vector, preconditioned by rA
  /// \param r Preconditioner vector
  virtual void precondition(scalarVector &w, const scalarVector &r) const;

  /// \brief return the reciprocal diagonal
  virtual const scalarVector &rD() const { return rD_; }
};

}  // namespace UNAP

#endif  // LDUDIAGPRECOND_HPP
