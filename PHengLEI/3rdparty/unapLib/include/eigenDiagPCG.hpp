#ifndef EIGENDIAGPCG_HPP
#define EIGENDIAGPCG_HPP

#include <vector>

#include "unapMatrix.hpp"
#include "unapPreconditioner.hpp"

namespace UNAP {
/// \brief compute maximum eigenvalue of the matrix by PCG loops
class EigenDiagPCG {
 private:
  /// \brief the maximum eigenvalue
  mutable scalar maxEigenValue_;

  /// \brief Sturm Sequence Method，find polynomial vector p(i) to find eigen
  /// value \param TruMatrix tridiagonal matrix \param lamb guess eigen value
  /// (used in bisection method) \param p a pvector polynomials (pi_1 to pi_n)
  /// \param s the number of sign change
  /// \param nPCGs nDiagPCGs_
  void h14Sturm(scalar **TriMatrix, const scalar lamb, Array<scalar> &p,
                label &s, const label nPCGs) const;

  /// \brief form the matrix from alphas and betas, used to calculate eigenvalue
  /// algorithm \param alphas alphas from nPCGs times of PCG \param betas betas
  /// from nPCGs times of PCG \param TriMatrix tridiagonal matrix \param nPCGs
  /// nDiagPCGs_
  void computeValueForMatrix(const scalarVector &alphas,
                             const scalarVector &betas, scalar **TriMatrix,
                             const label nPCGs) const;

  /// \brief estimate the range of eigenvalues
  /// \param TriMatrix tridiagonal matrix
  /// \param xBegin left bound of eigenvalues
  /// \param xEnd right bound of eigenvalues
  /// \param nPCGs nDiagPCGs_
  void determineEigRange(scalar **TriMatrix, scalar &xBegin, scalar &xEnd,
                         const label nPCGs) const;

  /// \brief compute eigenvalue
  /// \param alphas alphas from nPCGs times of PCG
  /// \param betas betas from nPCGs times of PCG
  /// \param nPCGs nDiagPCGs_
  /// \param k kth largest eigenvalue, default k=1 for the maximum eigenvalue
  void computeMaxEig(const scalarVector &alphas, const scalarVector &betas,
                     const label nPCGs, const label k) const;

  /// \brief diagPCG loops
  /// \param A input Matrix
  /// \param x input and output vector, x will be updated
  /// \param b input vector rhs
  /// \param precond Preconditioner
  /// \param nDiagPCGs loop count
  /// \param alphas alphas from nPCGs times of PCG
  /// \param betas betas from nPCGs times of PCG
  void diagPCGLoops(const Matrix &A, scalarVector &x, const scalarVector &b,
                    const Preconditioner &precond, const label nDiagPCGs,
                    scalarVector &alphas, scalarVector &betas) const;

  /// \brief allocate memory space of the 2D matrix
  /// \tparam T type of the 2D matrix
  /// \param nPCGs loop count, number of row in 2D matrix
  /// \return return the pointer of the 2D matrix
  template <typename T>
  T **allocateSym2D(label nPCGs) const;

  /// \brief free memory space of the 2D matrix
  /// \tparam T type of the 2D matrix
  /// \param arr pointer of the 2D matrix
  /// \param nPCGs number of rows in 2D matrix
  template <typename T>
  void deleteSym2D(T **arr, label nPCGs) const;

 public:
  /// \brief constructor
  /// \param A input Matrix
  /// \param x input and output vector, x will be updated
  /// \param b input vector rhs
  /// \param precond Preconditioner
  /// \param nDiagPCGs loop count
  EigenDiagPCG(const Matrix &A, scalarVector &x, const scalarVector &b,
               const Preconditioner &precond, const label nDiagPCGs);

  /// \brief return the maximum eigenvalue
  /// \return return value of the maximum eigenvalue
  scalar maxEigenValue() const { return maxEigenValue_; }
};

// allocate memory
template <typename T>
T **UNAP::EigenDiagPCG::allocateSym2D(label n) const {
  T **arr2D;
  arr2D = new T *[n];

  for (label i = 0; i < n; i++) {
    arr2D[i] = new T[n];
  }

  return (T **)arr2D;
}

// free memory
template <typename T>
void UNAP::EigenDiagPCG::deleteSym2D(T **arr, label n) const {
  if (arr != NULL) {
    for (label i = 0; i < n; i++) {
      delete[] arr[i];
    }

    delete[] arr;
    arr = NULL;
  }
}

}  // namespace UNAP
#endif  //-EIGENDIAGCG_HPP
