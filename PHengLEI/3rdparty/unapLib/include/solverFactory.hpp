#ifndef SOLVERS_SOLVERFACTORY_HPP
#define SOLVERS_SOLVERFACTORY_HPP
/*! \file solverFactory.hpp
 *  \author Zhao Chengpeng (chengpeng_zhao@foxmail.com)
 *  \date 2022-11-21
 *  \brief factory for solvers and preconditioners
 */

#include "unapMatrixSolver.hpp"
#include "unapPreconditioner.hpp"

namespace UNAP {
class LduMatrix;

typedef Preconditioner *(*precondCreatorPtr)(const LduMatrix &A);
class PrecondFactory {
 public:
  static HASHMAP<string, precondCreatorPtr> precondMap_;
  static Preconditioner *buildPrecond(const LduMatrix &A, string name) {
    return precondMap_.find(name)->second(A);
  }
};

typedef Solver *(*solverCreatorPtr)(Preconditioner &precond);
class SolverFactory {
 public:
  static HASHMAP<string, solverCreatorPtr> solverMap_;
  static Solver *buildSolver(Preconditioner &precond, string name) {
    return solverMap_.find(name)->second(precond);
  }
};

}  // namespace UNAP
#endif /* end of include guard: SOLVERS_SOLVERFACTORY_HPP */
