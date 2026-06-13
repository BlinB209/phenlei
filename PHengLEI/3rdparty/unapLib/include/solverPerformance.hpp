/*
 * @Author: guhf
 * @Email: hanfenggu@gmail.com
 * @Date: 2021-07-13 11:09:03
 * @Brief:
 */

#ifndef SOLVER_PERFORMANCE_HPP
#define SOLVER_PERFORMANCE_HPP

#include "unap.hpp"

namespace UNAP {
/// \brief Class returned by the solver, containing performance statistics
class SolverPerformance {
 private:
  /// \brief the initial residual
  scalar initialResidual_;

  /// \brief the final residual
  scalar finalResidual_;

  /// \brief the tempororary residual
  scalar previousResidual_;

  /// \brief number of iterations
  label nIterations_;

  /// \brief if results converged
  bool converged_;

  /// \brief if it is too small or not
  bool singular_;

 public:
  /// \brief Construct null
  SolverPerformance()
      : initialResidual_(0),
        finalResidual_(0),
        previousResidual_(0),
        nIterations_(0),
        converged_(false),
        singular_(false) {}

  /// \brief Construct from components
  /// \param iRes the initial residual
  /// \param fRes the final residual
  /// \param nIter number of iterations
  /// \param converged if the results meet the relTol or tolerance
  /// \param singular if it is too small or not
  SolverPerformance(const scalar iRes, const scalar fRes, const label nIter,
                    const bool converged, const bool singular)
      : initialResidual_(iRes),
        finalResidual_(fRes),
        previousResidual_(iRes),
        nIterations_(nIter),
        converged_(converged),
        singular_(singular) {}

  // Member functions

  /// \brief Return initial residual
  scalar initialResidual() const { return initialResidual_; }

  /// \brief return initial residual
  scalar &initialResidual() { return initialResidual_; }

  /// \brief set initial residual
  void initialResidual(scalar init) {
    initialResidual_ = init;
    finalResidual_ = init;
    previousResidual_ = init;
  }

  /// \brief return final residual
  scalar finalResidual() const { return finalResidual_; }

  /// \brief return final residual
  scalar &finalResidual() { return finalResidual_; }

  /// \brief set final residual
  void finalResidual(scalar init) {
    previousResidual_ = finalResidual_;
    finalResidual_ = init;
  }

  /// \brief return previous residual
  scalar previousResidual() const { return previousResidual_; }

  /// \brief return previous residual
  scalar &previousResidual() { return previousResidual_; }

  /// \brief return number of iterations
  label nIterations() const { return nIterations_; }

  /// \brief return number of iterations
  label &nIterations() { return nIterations_; }

  /// \brief has the solver converged?
  bool converged() const { return converged_; }

  /// \brief is the matrix singular?
  bool singular() const { return singular_; }

  /// \brief convergence check
  /// \param tolerance the absolute tolerance
  /// \param relTolerance the relative tolerance
  /// \param minIter minimum number of iterations
  /// \return if solver needs to finish
  bool checkConvergence(const scalar tolerance, const scalar relTolerance,
                        const label minIter);

  /// \brief convergence check
  /// \param tolerance the absolute tolerance
  /// \param minIter minimum number of iterations
  /// \return if solver needs to finish
  bool checkConvergence(const scalar tolerance, const label minIter);

  /// \brief convergence check
  /// \param tolerance the absolute tolerance
  /// \param relTolerance the relative tolerance
  /// \return if solver needs to finish
  bool checkConvergence(const scalar tolerance, const scalar relTolerance);

  /// \brief singularity check
  /// \param residual compared with SMALL(1e-37)
  bool checkSingularity(const scalar residual);

  /// \brief print residual for every step
  /// \param nIter iteration number
  /// \param precision print precision for scalar
  /// \param comm Communicator pointer
  void printResidual(label precision = 8);
};

}  // namespace UNAP

#endif
