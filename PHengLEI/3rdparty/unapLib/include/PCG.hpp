#ifndef PCG_HPP
#define PCG_HPP

#include "solverPerformance.hpp"
#include "unapMatrix.hpp"
#include "unapMatrixSolver.hpp"
#include "unapPreconditioner.hpp"

namespace UNAP {
/// \brief preconditioned Conjugate Gradient Method
/// see Barrett, R. et al., Templates for the solution of linear systems:
/// building blocks for iterative methods.
class PCG : public Solver {
 private:
  /// \brief use if the CG solver needs to create a null Preconditioner
  bool deletePrecondPtr_;

  /// \brief the Preconditioner
  Preconditioner *precondPtr_;

 public:
  /// \brief constructor
  /// \param other_comm communication domain
  PCG(Communicator *other_comm);

  /// \brief constructor
  /// \param precond Preconditioner used
  PCG(Preconditioner &precond);

  /// \brief destructor
  virtual ~PCG() {
    if (deletePrecondPtr_) {
      delete precondPtr_;
      precondPtr_ = NULL;
      deletePrecondPtr_ = false;
    }
  }

  /// \brief solve the matrix with this solver: Ax = b
  /// \param x vector needed to be solved or updated
  /// \param A Matrix
  /// \param b rhs
  virtual SolverPerformance solve(scalarVector &x, const Matrix &A,
                                  const scalarVector &b);
};
}  // namespace UNAP

#endif  // PCG_HPP
