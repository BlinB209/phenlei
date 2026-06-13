#ifndef SOLVERS_PRECONDITIONERS_LDUMGPRECOND_HPP
#define SOLVERS_PRECONDITIONERS_LDUMGPRECOND_HPP
#include "lduAgglomeration.hpp"
#include "unapMatrix.hpp"
#include "unapMultigrid.hpp"
#include "unapPreconditioner.hpp"

namespace UNAP {
class LduMatrix;

/// \brief Simplified Diagonal-based Incomplete Cholesky preconditioner
/// used if matrix is symmetric
class LduMGPrecond : public Preconditioner {
 private:
  /// \brief the reciprocal diagonal
  scalarVector rD_;

  /// Number of V-cycles to perform
  label nVcycles_;

  MGSolver *mgSolver_;
  LduAgglomeration *aggl_;
  /// \brief coarse grid correction fields, Ae=r, e
  mutable PtrList<scalarVector> coarseCorrFields_;
  /// \brief coarse grid sources, Ae=r, r
  mutable PtrList<scalarVector> coarseSources_;

 public:
  /// \brief constructor
  /// \param A matrix input
  LduMGPrecond(const LduMatrix &A, label nCycles = 1);

  /// \brief destructor
  virtual ~LduMGPrecond() {
    DELETE_POINTER(mgSolver_);
    DELETE_POINTER(aggl_);
  }

  /// \brief return wA the preconditioned form of residual rA
  /// \param w resulting vector, preconditioned by rA
  /// \param r Preconditioner vector
  virtual void precondition(scalarVector &w, const scalarVector &r) const;

  /// \brief return the reciprocal diagonal
  virtual const scalarVector &rD() const { return rD_; }
};

}  // namespace UNAP

#endif /* end of include guard: SOLVERS_PRECONDITIONERS_LDUMGPRECOND_HPP */
