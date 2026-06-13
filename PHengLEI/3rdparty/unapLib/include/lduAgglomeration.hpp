#ifndef LDUAGGLOMERATION_HPP
#define LDUAGGLOMERATION_HPP
/*! \file lduAgglomeration.hpp
 *  \author Zhao Chengpeng (chengpeng_zhao@foxmail.com)
 *  \modified 2023-02-20
 *  \brief simple aggregation(which OpenFOAM uses) has been replaced by multiple
 * pairwise aggregation
 */

#include "unapAgglomeration.hpp"
#include "unapLduMatrix.hpp"

namespace UNAP {
/// \brief Pairwise aggregation
/// NOTAY Y. AN AGGREGATION-BASED ALGEBRAIC MULTIGRID METHOD[J]. 2010: 24.
class LduAgglomeration : public Agglomeration {
 protected:
  /// \brief max number of levels
  label maxLevels_;

  /// \brief number of cells in coarsest level
  label rowSizeInCoarsestLevel_;

  /// \brief the number of cells in each level
  labelVector rowSize_;

  /// \brief number of levels created
  label nCreatedLevels_;

  /// \brief number of levels to merge, 1 = don't merge, 2 = merge pairs etc.
  label mergeLevels_;

  /// \brief finest matrix
  const LduMatrix &finestMatrix_;

  /// \brief hierarchy of coarse matrix levels
  PtrList<LduMatrix> coarseMatrixLevels_;

  /// \brief cell restriction addressing array.
  /// maps from the finer to the coarser level.
  PtrList<labelVector> restrictAddressing_;

  /// \brief face restriction addressing array.
  /// maps from the finer to the coarser level.
  /// positive indices map the finer faces which form part of the boundary
  /// of the coarser cells to the corresponding coarser cell face.
  /// negative indices map the finer faces which are internal to the
  /// coarser cells to minus the corresponding coarser cell index minus 1.
  PtrList<labelVector> faceRestrictAddressing_;

  /// \brief face restriction addressing in interface patch
  /// maps from the finer to the coarser level.
  PtrList<labelVector> patchFaceRestrictAddressing_;

  /// \brief do rcm reordering for coarse level i if ifRCMforLevels[i]==true
  bool *ifRCMforLevels_;

  /// \brief assemble coarse mesh addressing
  /// \param fineLevelIndex fine level specified
  void agglomerateLduAddressing(const label fineLevelIndex);

  /// \brief shrink the number of levels to that specified
  /// \param nCreatedLevels number of resulting levels specified
  void compactLevels(const label nCreatedLevels);

  /// \brief check the need for further Agglomeration
  /// \param nCoarseCells the number of rows in the coarsest level
  bool continueAgglomerating(const label nCoarseCells,
                             const label nCoarseCellsOld) const;

  /// \brief combine levels
  /// \param vurLevel current level
  void combineLevels(const label curLevel);

  /// \brief calculate and return Agglomeration of given level
  /// \param nCoarseCells number of rows in coarse level
  /// \param fineA fine level matrix
  /// \param
  virtual labelVector &agglomerate(label &nCoarseCells, const LduMatrix &fineA,
                                   label level);

  /// \brief agglomerate coarse matrix
  /// \param fineLevelIndex fine level index
  virtual void agglomerateMatrix(const label fineLevelIndex);

  /// \brief return the matrix of giving level
  /// \param leveli index of giving level
  const LduMatrix &matrixLevel(const label leveli) const;

 public:
  /// \brief constructors
  /// \param A input matrix to produce coarse level matrix
  /// \param rowSizeInCoarsestLevel
  LduAgglomeration(const LduMatrix &A, label rowSizeInCoarsestLevel = 100);

  /// \brief destructor
  virtual ~LduAgglomeration();

  /// \brief agglomerate all levels starting from default face weights
  virtual void agglomerate();

  /// \brief return the number of total levels
  virtual label size() const { return nCreatedLevels_; }

  /// \brief access to coarse matrix level of giving level
  /// \param leveli index of giving level
  virtual const LduMatrix &coarseMatrixLevels(const label leveli) const {
    return coarseMatrixLevels_[leveli];
  }

  /// \brief access to coarse matrix level of giving level
  /// \param leveli index of giving level
  LduMatrix &coarseMatrix(const label leveli) const {
    return coarseMatrixLevels_[leveli];
  }

  /// \brief access to restrictAddressing_
  /// return cell restrict addressing of given level
  /// \param leveli index of giving level
  virtual const labelVector &restrictAddressing(const label leveli) const {
    return restrictAddressing_[leveli];
  }

  /// \brief access to faceRestrictAddressing_
  /// return face restrict addressing of given level
  /// \param leveli index of giving level
  virtual const labelVector &faceRestrictAddressing(const label leveli) const {
    return faceRestrictAddressing_[leveli];
  }

  /// \brief set the number of rows in the coarsest level
  /// \param i number of rows in the coarsest level
  virtual void SET_rowSizeInCoarsestLevel(const label i);

  /// \brief set the maximum number of levles
  /// \param i maximum number of levles
  virtual void SET_maxLevels(const label i);

  /// \brief reorder the Topology according to the results in MLB
  void AgglomerationReorderTopo();

  /// \brief construct the RSS iterator in coarse levels [0,maxL)
  /// \param maxL max allowed coarse levels to create RSS iterator, default is
  /// nCreatedLevels_
  void agglConstructRSSIterators(const label maxL = 1000);

  /// \brief construct the RSS iterator in coarse levels i
  void agglConstructRSSIterator(const label i);

  /// \brief do rcm reordering for coarse level [0,maxL)
  void SET_rcmLevels(const label maxL = 0);
  /// \brief do rcm reordering for coarse level i
  void SET_rcmLevel(const label i);
};

}  // namespace UNAP

#endif  // LDUAGGLOMERATION_HPP
