#ifndef LDUGAUSSSEIDELSMOOTHER_HPP
#define LDUGAUSSSEIDELSMOOTHER_HPP

#include "unapLduMatrix.hpp"
#include "unapSmoother.hpp"

namespace UNAP {
/// \brief Gauss-Seidel Smoother
class LduGaussSeidelSmoother : public Smoother {
 private:
 public:
  static const string typeName_;
  /// \brief constructors
  /// \param other_comm communication domain
  LduGaussSeidelSmoother(Communicator *other_comm) : Smoother(other_comm) {}

  /// \brief destructor
  virtual ~LduGaussSeidelSmoother() {}

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

#endif  // LDUGAUSSSEIDELSMOOTHER_HPP
