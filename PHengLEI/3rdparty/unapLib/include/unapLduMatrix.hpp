#ifndef UNAP_LDUMATRIX_HPP
#define UNAP_LDUMATRIX_HPP

#include "unapMatrix.hpp"

#ifdef SW_SLAVE
#include "iterator.hpp"
#include "multiLevelBlockIterator.hpp"
#include "rowSubsectionIterator.hpp"
#endif

namespace UNAP {
template <typename T>
class CSRMatrix;

/// \brief LDU type matrix
class LduMatrix : public Matrix {
 private:
  /// \brief number of cells
  label rowSize_;

  /// \brief addressing
  labelVector *lowerAddrPtr_, *upperAddrPtr_;

  /// \brief coefficients (not including interfaces)
  scalarVector *lowerPtr_, *diagPtr_, *upperPtr_;

  /// \brief interfaces
  Interfaces *interfacesPtr_;

  /// \brief losort addressing
  mutable labelVector *losortPtr_;

  /// \brief owner start addressing
  mutable labelVector *ownerStartPtr_;

  /// \brief losort start addressing
  mutable labelVector *losortStartPtr_;

  label *startRow_;

  /// \brief calculate losort
  void calcLosort() const;

  /// \brief calculate owner start
  void calcOwnerStart() const;

  /// \brief calculate losort start
  void calcLosortStart() const;

#ifdef SW_SLAVE
  /// \brief multi-level block iterator
  UNAT::MultiLevelBlockIterator *mlbIter_;

  /// \brief start from 1, can be negative, which stands for lower part
  label *unatEdgeMap_;

  /// \brief start from 0
  label *unatCellMap_;

  /// \brief row subsection iterator
  mutable UNAT::RowSubsectionIterator *rssIter_;

  /// \brief mlbIter_ and unatIter_ may be replaced with unatIter_
  /// according to the new interfaces from unat
  mutable UNAT::Iterator *unatIter_;

  mutable UNAT::Topology *topo_;

  enum swSpeedUpType { USE_MLB, USE_RSS, NONE } swSpeedUp_;

#endif

 public:
  /// \brief constructor
  /// \param other_comm communication domain
  LduMatrix(Communicator *other_comm);

  /// \brief constructor
  /// \param rowSize local row size
  /// \param other_comm communication domain
  LduMatrix(label rowSize, Communicator *other_comm);

  /// \brief constructor
  /// \param rowSize number of rows in the matrix
  /// \param lowerAddr the row index of the upper part of the matrix
  /// \param upperAddr the col index of the upper part of the matrix
  /// \param lower lower part of the matrix
  /// \param diag diagonal part of the matrix
  /// \param upper upper part of the matrix
  /// \param other_comm communication domain
  LduMatrix(const label &rowSize, const labelVector &lowerAddr,
            const labelVector &upperAddr, const scalarVector &lower,
            const scalarVector &diag, const scalarVector &upper,
            Communicator *other_comm);

  /// \brief constructor
  /// \param rowSize number of rows in the matrix
  /// \param lowerAddr the row index of the upper part of the matrix
  /// \param upperAddr the col index of the upper part of the matrix
  /// \param lower lower part of the matrix
  /// \param diag diagonal part of the matrix
  /// \param upper upper part of the matrix
  /// \param reUse reUse the vectors, so the matrix will allocate new space
  /// \param other_comm communication domain
  LduMatrix(const label &rowSize, const labelVector &lowerAddr,
            const labelVector &upperAddr, const scalarVector &lower,
            const scalarVector &diag, const scalarVector &upper,
            const bool reUse, Communicator *other_comm);

  /// \brief constructor only Topology
  /// \param rowSize number of rows in the matrix
  /// \param lowerAddr the row index of the upper part of the matrix
  /// \param upperAddr the col index of the upper part of the matrix
  /// \param other_comm communication domain
  LduMatrix(const label &rowSize, const labelVector &lowerAddr,
            const labelVector &upperAddr, Communicator *other_comm);

  /// \brief constructor only Topology
  /// \param rowSize number of rows in the matrix
  /// \param lowerAddr the row index of the upper part of the matrix
  /// \param upperAddr the col index of the upper part of the matrix
  /// \param reUse reUse the vectors, so the matrix will allocate new space
  /// \param other_comm communication domain
  LduMatrix(const label &rowSize, const labelVector &lowerAddr,
            const labelVector &upperAddr, const bool reUse,
            Communicator *other_comm);

  LduMatrix(const CSRMatrix<scalar> &A);

  /// \brief destructor
  ~LduMatrix();

  /// \brief access to lower addressing
  virtual labelVector &lowerAddr() const;

  /// \brief access to upper addressing
  virtual labelVector &upperAddr() const;

  /// \brief access to lower coefficients
  virtual scalarVector &lower() const;

  /// \brief access to diagonal coefficients
  virtual scalarVector &diag() const;

  /// \brief access to upper coefficients
  virtual scalarVector &upper() const;

  label *startRow() const;

  label startRow(label i) const;

  void SET_startRow(label *startRow) {
    DELETE_POINTERS(startRow_)
    startRow_ = startRow;
  }

  /// \brief set lower addressing
  /// \param newLowerAddr the new lower addressing
  void SET_lowerAddr(labelVector &newLowerAddr) {
    ALLOCATE_POINTER(lowerAddrPtr_, newLowerAddr, labelVector)
  }

  /// \brief set upper addressing
  /// \param newUpperAddr the new upper addressing
  void SET_upperAddr(labelVector &newUpperAddr) {
    ALLOCATE_POINTER(upperAddrPtr_, newUpperAddr, labelVector)
  }

  /// \brief set lower coefficients
  /// \param newLower the new lower coefficients
  void SET_lower(scalarVector &newLower) {
    if (this->symm()) {
      lowerPtr_ = NULL;
    }

    ALLOCATE_POINTER(lowerPtr_, newLower, scalarVector)
  }

  /// \brief set lower size
  /// \param newSize the new lower size
  void SET_lower(label newSize) {
    DELETE_POINTER(lowerPtr_)

    lowerPtr_ = new scalarVector(newSize, this->commcator_);
  }

  /// \brief set upper coefficients
  /// \param newUpper the new upper coefficients
  void SET_upper(scalarVector &newUpper) {
    ALLOCATE_POINTER(upperPtr_, newUpper, scalarVector)
  }

  /// \brief set upper size
  /// \param newSize the new upper size
  void SET_upper(label newSize) {
    DELETE_POINTER(upperPtr_)

    upperPtr_ = new scalarVector(newSize, this->commcator_);
  }

  /// \brief set diagonal coefficients
  /// \param newDiag the new diagonal coefficients
  void SET_diag(scalarVector &newDiag) {
    ALLOCATE_POINTER(diagPtr_, newDiag, scalarVector)
  }

  /// \brief set diagonal size
  /// \param newSize the new diagonal size
  void SET_diag(label newSize) {
    DELETE_POINTER(diagPtr_)

    diagPtr_ = new scalarVector(newSize, this->commcator_);
  }

  /// \brief return the diagonal coefficients
  scalarVector &diag() { return *diagPtr_; }

  /// \brief free the lower coefficients
  void freeLower() {
    if (!this->symm()) {
      DELETE_POINTER(lowerPtr_);
      lowerPtr_ = upperPtr_;
    }
  }

  /// \brief make upperPtr_ = lowerPtr_
  void setSymm() {
    if (upperPtr_) {
      DELETE_POINTER(lowerPtr_);
      lowerPtr_ = upperPtr_;
    } else if (lowerPtr_) {
      DELETE_POINTER(upperPtr_);
      upperPtr_ = lowerPtr_;
    }
  }
  void setSymm(bool sym) {
    if (sym) setSymm();
  }

  /// \brief return if the matrix is symmetric or asymmetric
  virtual bool symm() const {
    if (lowerPtr_ == upperPtr_) {
      return true;
    } else {
      return false;
    }
  }

  /// \brief retrun row size the matrix
  virtual label size() const { return rowSize_; }

  /// \brief  set the row size the matrix
  void size(label size);

  /// \brief sparse matrix-vector multiplication: y=Ax
  /// \param Apsi the result vector y
  /// \param psi vector x
  virtual void spMV(scalarVector &Apsi, const scalarVector &psi) const;

  /// \brief return losort addressing
  const labelVector &losortAddr() const;

  /// \brief return owner start addressing
  const labelVector &ownerStartAddr() const;

  /// \brief return losort start addressing
  const labelVector &losortStartAddr() const;

  /// \brief return the reference of the matrix interface
  virtual Interfaces &matrixInterfaces() const { return *interfacesPtr_; }

  /// \brief set the reference of the matrix interface
  /// \param a the new matrix interface
  virtual void matrixInterfaces(Interfaces &a) { interfacesPtr_ = &a; }

  /// \brief set the size of patches
  /// \param size number of the patches
  void matrixInterfaces(const label size) {
    PtrList<Patch> *patchesPtr = new PtrList<Patch>(size);
    interfacesPtr_ = new Interfaces(*patchesPtr, this->commcator_);
  }

  /// Not use any more
  /// \brief create the interface Topology
  /// \param nNeiProcs number of neighbor processors, equal to number of patches
  /// \param destRank array pointer of the neighbor processor ID
  /// \param offDiagRows array pointer of the local addressing of the
  /// inter-processor faces
  //  \param offDiagStarts start index for every patch in
  /// offDiagRow array
  // void createInterfacesTopology(const label nNeiProcs, const label *destRank,
  //                              const label *offDiagRows,
  //                              const label *offDiagStarts);

  /// \brief create the interface Topology
  /// \param nNeiProcs number of neighbor processors, equal to number of patches
  /// \param destRank array pointer of the neighbor processor ID
  /// \param offDiagRows array pointer of the local row addressing of the
  /// inter-processor faces
  /// \param offDiagCols array pointer of the local col addressing of the
  /// inter-processor faces
  //  \param offDiagStarts start index for every patch in
  /// offDiagRow array
  void createInterfacesTopology(const label nNeiProcs, const label *destRank,
                                const label *offDiagRows,
                                const label *offDiagCols,
                                const label *offDiagStarts);

  /// \brief fill the coefficients in the interface
  /// \param offDiagStarts start index for every patch in offDiagCoeffs array
  /// \param offDiagCoeffs the coefficients array
  void fillInterfacesCofficients(const label *offDiagStarts,
                                 const scalar *offDiagCoeffs);

  /// \brief initialize interfaces, start to sending
  /// \param psi vector needed to be updated
  virtual void initInterfaces(const scalarVector &psi) const;

  /// \brief update interfaces, finish receiving
  /// \param Apsi the result array, which will be calculated in this function
  virtual void updateInterfaces(scalarVector &Apsi) const;

  /// \brief fill coefficients, allocate new memory space
  /// \param diag diagonal part
  /// \param upper upper part
  /// \param lower lower part
  void setMatrixCoeffients(const scalarVector &diag, const scalarVector &upper,
                           const scalarVector &lower);

  /// \brief fill coefficients, reuse vectors input
  /// \param diag diagonal part
  /// \param upper upper part
  /// \param lower lower part
  /// \param reUse reuse the vectors input
  void setMatrixCoeffients(const scalarVector &diag, const scalarVector &upper,
                           const scalarVector &lower, const bool reuse);

  /// \brief fill Topology, reuse the vectors input
  /// \param rowSize matrix local rowSize
  /// \param upperAddr column index for upper part
  /// \param lowerAddr row index for upper part
  /// \param reUse reuse the vectors input
  void setMatrixTopology(const label rowSize, const labelVector &upperAddr,
                         const labelVector &lowerAddr,
                         const bool reUse = false);
  /*! calc local nonzero number
   *  \return scalar, lnnz which includes interfaces entities
   */
  virtual label lnnz() const;

#ifdef SW_SLAVE
/// \brief spMV神威加速临界值，指的是矩阵上三角系数的数量
/// 目前用于判断MLB加速下粗层是否需要重排
#define SpMVAccSize 10000
  // #define SpMVAccSize 500

  /// \brief construct the MLB iterator
  void constructMLBIterator();

  /// \brief return the MLB iterator
  UNAT::MultiLevelBlockIterator *mlbIter() const { return mlbIter_; }

  /// \brief reorder the coefficients in lower, diagonal, upper
  void reorderLDUValues();

  /// \brief restore the vector
  /// \param vv vector input
  void restoreVector(scalarVector &vv);

  /// \brief reorder the vector
  /// \param vv vector input
  void reorderVector(scalarVector &vv);

  /// \brief return the edge map in MLB
  const label *unatEdgeMap() const { return unatEdgeMap_; }

  /// \brief return the vertex map in MLB
  const label *unatCellMap() const { return unatCellMap_; }

  /// \brief set the giving MLB iterator
  void setMlbIter(UNAT::MultiLevelBlockIterator *mlbIter);

  /// \brief free edge and vertex map
  /// TODO:temperary
  void unatMapFree() {
    DELETE_POINTERS(unatEdgeMap_);
    DELETE_POINTERS(unatCellMap_);
  }

  /// \brief construct the RSS iterator
  void constructRSSIterator();

#endif
};

}  // namespace UNAP

#endif  // UNAP_LDUMATRIX_HPP
