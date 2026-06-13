#ifndef GMRES_HPP
#define GMRES_HPP

#include "solverPerformance.hpp"
#include "unapMatrix.hpp"
#include "unapMatrixSolver.hpp"
#include "unapPreconditioner.hpp"

namespace UNAP {
/// \brief generalized minimal residual method solver
/// see Barrett, R. et al., Templates for the solution of linear systems:
/// building blocks for iterative methods.
class GMRES : public Solver {
 private:
  /// \brief use if the CG solver needs to create a null Preconditioner
  bool deletePrecondPtr_;

  /// \brief the Preconditioner
  Preconditioner *precondPtr_;

  /// \brief Krylov space dimension, GMRES(m) m=nDirs_, 为固定的重启动参数
  /// nDirs_过小会导致计算x与r的次数增大，从而明显降低求解效率，nDirs_过大则程序需要的空间越多
  label nDirs_;

  /// \brief Givens rotation
  void givensRotation(const scalar &H, const scalar &beta, scalar &c,
                      scalar &s) const;
  /// \brief apply rotation
  void applyRotation(scalar &dx, scalar &dy, const scalar &c,
                     const scalar &s) const;

 public:
  /// \brief constructor
  /// \param other_comm communication domain
  GMRES(Communicator *other_comm);

  /// \brief constructor
  /// \param precond Preconditioner used
  GMRES(Preconditioner &precond);

  /// \brief constructor
  /// \param m 子空间维数
  GMRES(Preconditioner &precond, label m);

  /// \brief destructor
  virtual ~GMRES() {
    if (deletePrecondPtr_) {
      delete precondPtr_;
      precondPtr_ = NULL;
      deletePrecondPtr_ = false;
    }
  }

  /// \brief solve the matrix with this solver: Ax = b, 右预处理
  /// \param x vector needed to be solved or updated
  /// \param A Matrix
  /// \param b rhs
  virtual SolverPerformance solve(scalarVector &x, const Matrix &A,
                                  const scalarVector &b);
};

}  // namespace UNAP

#endif  // GMRES_HPP
