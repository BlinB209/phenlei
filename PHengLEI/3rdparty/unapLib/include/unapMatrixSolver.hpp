#ifndef UNAP_MATRIX_SOLVER_HPP
#define UNAP_MATRIX_SOLVER_HPP

#include "solverPerformance.hpp"
#include "unapMatrix.hpp"
#include "unapVector.hpp"

namespace UNAP {
/// \brief abstract base-class for matrix solvers
class Solver {
 protected:
  /// \brief maximum number of iterations in the iteration
  label maxIter_;

  /// \brief minimum number of iterations in the iteration
  label minIter_;

  /// \brief convergence tolerance relative to the initial
  scalar relTol_;

  /// \brief final convergence tolerance
  scalar tolerance_;

  /// \brief print converging process
  bool ifPrint_;

  /// \brief communicator
  Communicator *commcator_;

 public:
  /// \brief constructors
  /// \param other_comm communication domain
  Solver(Communicator *other_comm);

  /// \brief destructor
  virtual ~Solver() {}

  /// \brief solve the equation Ax=b
  /// \param x unkown vector
  /// \param A Matrix
  /// \param b rhs
  /// \return return the solver results, including residual, iteration numbers
  virtual SolverPerformance solve(scalarVector &x, const Matrix &A,
                                  const scalarVector &b) = 0;

  /// \brief return the matrix norm used to normalize the residual for the
  /// stopping criterion
  /// \param soure input source
  scalar normFactor(const scalarVector &source) const;

  /// \brief return the maximum iteration numbers
  inline label maxIter() const { return maxIter_; }

  /// \brief return the minmum iteration numbers
  inline label minIter() const { return minIter_; }

  /// \brief return the relative tolerance
  inline scalar relTol() const { return relTol_; }

  /// \brief return the absolute tolerance
  inline scalar tolerance() const { return tolerance_; }

  /// \brief return if it needs to print verbosely
  inline bool ifPrint() const { return ifPrint_; }

  /// \brief set the maximum iteration numbers
  /// \param maxIter the maximum iteration numbers
  void SET_maxIter(label maxIter);

  /// \brief set the minmum iteration numbers
  /// \param minIter the minmum iteration numbers
  void SET_minIter(label minIter);

  /// \brief set the relative tolerance
  /// \param relTol the relative tolerance
  void SET_relTol(scalar relTol);

  /// \brief set the absolute tolerance
  /// \param relTol the absolute tolerance
  void SET_tolerance(scalar tolerance);

  /// \brief set the value if it needs to print verbosely
  /// \param ifPrint the value if it needs to print verbosely
  void SET_ifPrint(bool ifPrint);

  /// \brief set communicator

  void setCommunicator(Communicator *other_comm) { commcator_ = other_comm; }

  /// \brief get communicator
  /// \param other_comm communication domain
  Communicator *getCommunicator() const { return commcator_; }
};

}  // namespace UNAP
#endif
