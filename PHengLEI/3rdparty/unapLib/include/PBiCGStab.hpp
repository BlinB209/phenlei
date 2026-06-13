#ifndef PBICGSTAB_HPP
#define PBICGSTAB_HPP

#include "solverPerformance.hpp"
#include "unapMatrix.hpp"
#include "unapMatrixSolver.hpp"
#include "unapPreconditioner.hpp"

namespace UNAP {
/// \brief Preconditioned BiConjugate Gradient Stabilized method
/// see Barrett, R. et al., Templates for the solution of linear systems:
/// building blocks for iterative methods.
class PBiCGStab : public Solver {
 private:
  /// \brief use if the CG solver needs to create a null Preconditioner
  bool deletePrecondPtr_;

  /// \brief the Preconditioner
  Preconditioner *precondPtr_;

 public:
  /// \brief constructor
  /// \param other_comm communication domain
  PBiCGStab(Communicator *other_comm);

  /// \brief constructor
  /// \param precond Preconditioner used
  /// 采用的是右预处理
  PBiCGStab(Preconditioner &precond);

  /// \brief destructor
  virtual ~PBiCGStab() {
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

  /// \brief set the Preconditioner
  /// \param precond Preconditioner used
  void SET_Preconditioner(Preconditioner &precond) {
    DELETE_POINTER(precondPtr_);
    precondPtr_ = &precond;
  }
};
}  // namespace UNAP

#endif  // PBICGSTAB_HPP
