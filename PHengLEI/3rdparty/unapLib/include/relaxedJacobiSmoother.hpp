#ifndef MULTIGRID_SMOOTHERS_RELAXEDJACOBISMOOTHER_HPP
#define MULTIGRID_SMOOTHERS_RELAXEDJACOBISMOOTHER_HPP

/*! \file relaxedJacobiSmoother.hpp
 *  \author Zhao Chengpeng (chengpeng_zhao@foxmail.com)
 *  \date 2022-10-24
 *  \brief relaxed jacobi smoother
 */

#include "unapLduMatrix.hpp"
#include "unapSmoother.hpp"

namespace UNAP {
/// \brief relaxed Jacobi Smoother
class RelaxedJacobiSmoother : public Smoother {
 private:
  scalar omega_;

 public:
  static const string typeName_;
  /// \brief constructors
  /// \param other_comm communication domain
  RelaxedJacobiSmoother(Communicator *other_comm, scalar omega = 2. / 3)
      : Smoother(other_comm), omega_(omega) {}

  /// \brief destructor
  virtual ~RelaxedJacobiSmoother() {}

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

  /// \brief init interface, do nothing in Gauss-Seidel
  virtual void init() const {}
};

}  // namespace UNAP

#endif /* end of include guard: MULTIGRID_SMOOTHERS_RELAXEDJACOBISMOOTHER_HPP \
        */
