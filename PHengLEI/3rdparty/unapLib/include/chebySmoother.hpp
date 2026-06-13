#ifndef CHEBYSMOOTHER_CPP
#define CHEBYSMOOTHER_CPP

#include "unapLduMatrix.hpp"
#include "unapSmoother.hpp"

namespace UNAP {
/// \brief chebyshev Smoother
/// see Barrett, R. et al., Templates for the solution of linear systems:
/// building blocks for iterative methods.
class ChebySmoother : public Smoother {
 private:
  /// \brief Number of PCGs to get alphas and betas, used to calculate
  /// eigenvalue
  label nDiagPCGs_;

  /// \brief Upper bound of the bounding ellipse of the eigenvalues of the
  /// matrix A
  mutable scalar maxEigPCG_;

  /// \brief Detect if eigenvalue has been calculated,
  /// if not, Chebyshev Smoother will call PCGs
  mutable bool eigFirstTimeComputed_;

  /// \brief Estimate the minEig by maxEig divided by eigRatio
  /// minEig = maxEig / eigRatioCheby_
  scalar eigRatioCheby_;

  /// \brief factor to enlarge the maxEig, because
  /// the maxEig obtained is usually underestimated
  /// maxEig *= boostFactorCheby_
  scalar boostFactorCheby_;

  /// \brief not used here
  scalar eigRatioCoarest_;

  //-detect if using preSmooth_
  //-if using, first Amul operation can be ignored
  label preSmoothUsing_;

 public:
  static const string typeName_;
  /// \brief constructor
  /// \param other_comm communication domain
  ChebySmoother(Communicator *other_comm)
      : Smoother(other_comm),
        nDiagPCGs_(10),
        maxEigPCG_(0.0),
        eigFirstTimeComputed_(true),
        eigRatioCheby_(1.5),     // 30
        boostFactorCheby_(1.1),  // 1.05
        eigRatioCoarest_(1.0),
        preSmoothUsing_(0) {}

  /// \brief smooth the solution for a given number of sweeps
  /// \param x vector needed to be updated
  /// \param A Matrix input
  /// \param b rhs vector
  /// \param nSweeps number of sweeps
  virtual void smooth(scalarVector &x, const Matrix &A, const scalarVector &b,
                      const label nSweeps) const;

  /// \brief smooth the solution for a given number of sweeps
  /// \param x vector needed to be updated
  /// \param A matrix input
  /// \param b rhs vector
  /// \param nSweeps number of sweeps
  void smooth(scalarVector &x, const LduMatrix &A, const scalarVector &b,
              const label nSweeps) const;

  /// \brief eigFirstTimeComputed_ is used to determine whether the eigenvalue
  /// needs to be calculated
  virtual void init() const { eigFirstTimeComputed_ = true; }
};
}  // namespace UNAP

#endif  // CHEBYSMOOTHER_CPP
