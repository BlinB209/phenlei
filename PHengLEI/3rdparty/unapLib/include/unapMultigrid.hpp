#ifndef UNAP_MULTIGRID_HPP
#define UNAP_MULTIGRID_HPP

#include "solverPerformance.hpp"
#include "unapAgglomeration.hpp"
#include "unapMatrix.hpp"
#include "unapMatrixSolver.hpp"
#include "unapSmoother.hpp"

namespace UNAP {
class LduMGPrecond;

/// \brief Agglomeration based algebraic multigrid solver(concurrently)
class MGSolver : public Solver {
  friend LduMGPrecond;

 private:
  const Matrix &finestMatrix_;
  /// \brief number of pre-smoothing sweeps at first coarser level
  label nPreSweeps_;
  /// \brief number of post-smoothing sweeps at first coarser level
  label nPostSweeps_;
  /// \brief max number of pre-smoothing sweeps(cut-off)
  label maxPreSweeps_;
  /// \brief max number of post-smoothing sweeps(cut-off)
  label maxPostSweeps_;
  /// \brief number of smoothing sweeps on finest mesh
  label nFinestSweeps_;
  /// \brief level multiplier for additional pre-sweeps at subsequent coarser
  /// levels;
  label preSweepsLevelMultiplier_;
  /// \brief level multiplier for additional post-sweeps at subsequent coarser
  /// levels;
  label postSweepsLevelMultiplier_;
  /// \brief choose if the corrections should be scaled.
  /// by default corrections for symmetric matrices are scaled
  /// but not for asymmetric matrices.
  bool scaleCorrection_;
  /// \brief the Agglomeration
  const Agglomeration *Agglomeration_;
  /// \brief Smoother type
  const char *smootherType_;
  /// \brief Smoother type set
  static std::set<string> smTypeSet_;
  /// \brief Smoother in each level
  PtrList<Smoother> *Smoothers_;
  /// \brief if Smoothers_ is built by MGSolver
  bool hasSmoothers_;
  /// \brief coarsest preconditioner name
  const char *coarsestPrecond_;
  /// \brief coarsest solver name
  const char *coarsestSolver_;
  /// \brief grid levels
  label levels_;
  /// \brief if no coarse grid
  bool noCoarse_;
  /// \bried run Kcycle rather than Vcycle when leveli>=levelKcycle_
  label levelKcycle_;

  /// \brief Acceleration of convergence
  /// Calculate and return the scaling factor from Acf, coarseSource and
  /// coarseField. At the same time do a Jacobi iteration on the coarseField
  /// using the Acf provided after the coarseField values are used for the
  /// scaling factor.
  /// \param field  vector will be updated
  /// \param source rhs, used to updated field using Jacobi iteration
  /// \param A Matrix
  scalar scalingFactor(scalarVector &field, const scalarVector &source,
                       const Matrix &A) const;

  /// \brief perform a single AMG V-cycle with pre, post and finest smoothing.
  /// \param psi unknown field in the finest level, will be updated
  /// \param source rhs in the finest level
  /// \param finestCorrection correction field in the finest level
  /// \param finestResidual residual in the finest level
  /// \param coarseCorrFields correction fields in coarse levels
  /// \param coarseSources rhs fields in coarse levels
  void Vcycle(label level, scalarVector &psi, const scalarVector &source,
              PtrList<scalarVector> &coarseCorrFields,
              PtrList<scalarVector> &coarseSources, label k);
  void Kcycle(label level, scalarVector &psi, const scalarVector &source,
              PtrList<scalarVector> &coarseCorrFields,
              PtrList<scalarVector> &coarseSources, label k);

  /// \brief solve the coarsest level with an iterative solver
  /// \param coarsestCorrField correction field in the coarsest level
  /// \param coarsestSource rhs field in the coarsest level
  void solveCoarsestLevel(scalarVector &coarsestCorrField,
                          const scalarVector &coarsestSource) const;

 public:
  /// \brief constructors
  /// \param A Matrix
  /// \param agglomerator matrix Topology in coarse levels
  MGSolver(const Matrix &A, const Agglomeration &agglomerator);

  /// \brief constructors
  /// \param A Matrix
  /// \param agglomerator matrix Topology in coarse levels
  /// \param Smoothers Smoothers in all levels
  MGSolver(const Matrix &A, const Agglomeration &agglomerator,
           PtrList<Smoother> &Smoothers);

  /// \brief destructor
  virtual ~MGSolver();

  /// \brief solve the matrix with this solver: Ax = b
  /// \param x vector needed to be solved or updated
  /// \param A Matrix
  /// \param b rhs filed
  virtual SolverPerformance solve(scalarVector &x, const Matrix &A,
                                  const scalarVector &b);

  /// \brief choose a smoother type
  /// Jacobi relaxedJacobi Smoother
  /// GS gaussSeidel Smoother
  /// Cheby chebyshev Smoother
  void SET_smootherType(string type);
  void SET_smootherType(const char *type);
  /// \brief initialize the Smoother in all levels
  /// should be done after SET_smootherType
  void initSmoothers(PtrList<Smoother> *smPtr = NULL);
  /// \brief initialize MGSolver
  void setUp();
  /// \brief set the number of pre sweeps in each level
  /// \param nPreSweeps number of pre sweeps in each level
  void SET_nPreSweeps(const label nPreSweeps);
  /// \brief set the number of post sweeps in each level
  /// \param nPostSweeps number of post sweeps in each level
  void SET_nPostSweeps(const label nPostSweeps);
  /// \brief set the number of sweeps in the finest level
  /// \param nFinestSweeps number of sweeps in the finest level
  void SET_nFinestSweeps(const label nFinestSweeps);
  /// \brief set the max number of pre-smoothing sweeps(cut-off)
  void SET_maxPreSweeps(const label maxPreSweeps);
  /// \brief set the max number of post-smoothing sweeps(cut-off)
  void SET_maxPostSweeps(const label maxPostSweeps);
  /// \brief set the level multiplier for additional pre-sweeps at subsequent
  /// coarser levels;
  void SET_preSweepsLevelMultiplier(const label multiplier);
  /// \brief set the level multiplier for additional post-sweeps at subsequent
  /// coarser levels;
  void SET_postSweepsLevelMultiplier(const label multiplier);
  /// \brief set if using scaling factor to accelerate the convergence
  /// only work if matrix is symmetric
  /// \param scaleCorrection true(use scaling factor) or false
  void SET_scaleCorrection(bool scaleCorrection);

  /// \brief not implemented yet
  void SET_coarsestSolver(const string &solverName) const;
  void SET_coarsestSolver(const char *solverName) const;

  void SET_KcycleLevel(const label li) { levelKcycle_ = li; }
  label KcycleLevel() { return levelKcycle_; }

  /// \brief restrict (integrate by summation) face field
  /// \tparam T vector type
  /// \param cf coarse level vector
  /// \param ff fine level vector
  /// \param fineLevelIndex fine level index
  template <typename T>
  void restrictField(Vector<T> &cf, const Vector<T> &ff,
                     const label fineLevelIndex) const;
  /// \brief prolong (interpolate by injection) cell field
  /// \tparam T vector type
  /// \param cf coarse level vector
  /// \param ff fine level vector
  /// \param coarseLevelIndex coarse level index
  template <typename T>
  void prolongField(Vector<T> &ff, const Vector<T> &cf,
                    const label coarseLevelIndex) const;
  /// \brief return matrix of given level
  /// \param levlei level index
  const Matrix &matrix(const label leveli) const;
  /// \brief return a pointer to agglomeration
  const Agglomeration *agglomeration() const { return Agglomeration_; }

  void printStat(bool ifPrintProc = false);
};

}  // namespace UNAP
#endif  // UNAP_MULTIGRID_HPP
