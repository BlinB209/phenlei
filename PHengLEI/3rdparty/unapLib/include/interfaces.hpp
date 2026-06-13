#ifndef INTERFACES_HPP
#define INTERFACES_HPP

#include "PtrList.hpp"
#include "patch.hpp"

namespace UNAP {
/// \brief integration of patches, including all neighbors for current processor
class Interfaces {
 private:
  /// \brief patch list
  PtrList<Patch> &patches_;

  /// \brief send buffer for initializing and updating matrix interfaces
  mutable scalar *sendBuffer_;

  /// \brief receive buffer for initializing and updating matrix interfaces
  mutable scalar *recvBuffer_;

  /// \brief locPosition
  mutable label *sendLocPosition_;

  /// \brief locPosition
  mutable label *recvLocPosition_;

  /// \brief destRank
  mutable label *destRank_;

  /// \brief communicator
  mutable Communicator *commcator_;

  /// \brief total interfaces size
  label size_;
  /// \brief  compressed row size
  label rowSize_;
  /// \brief compressed row index
  label *compRow_;
  /// \brief  row index for every row (rowSize_ in total)
  label *rowIndex_;
  /// \brief relative col index
  label *coloffset_;
  /// \brief matrix coeff
  scalar *coef_;

 public:
  /// \brief constructor
  /// \param size number of patches
  /// \param other_comm communication domain
  Interfaces(const label size, Communicator *other_comm);

  /// \brief constructor
  /// \param patches list of patches
  /// \param other_comm communication domain
  Interfaces(PtrList<Patch> &patches, Communicator *other_comm);

  virtual ~Interfaces();

  /// \brief Initialize the update of interfaced interfaces for matrix
  /// operations \param psi scalar vector needed to be updated
  virtual void initMatrixInterfaces(const scalarVector &psi) const;

  /// \brief update interfaced interfaces for matrix operations
  /// \param psi scalar vector needed to be updated
  virtual void updateMatrixInterfaces(scalarVector &result) const;

  /// \brief return the patch from giving index
  /// \param i patch index
  virtual Patch &patchList(const label i) const { return patches_[i]; }

  /// \brief return the list of patches
  PtrList<Patch> &patchList() { return patches_; }

  /// \brief return the number of patches in this interface
  virtual label size() const { return patches_.size(); }

  /// \brief renumber the face-cell addressing for giving cell map
  /// \param cellMap the renumbered cell index
  void reorderIntFaceCells(const label *cellMap);

  /// \brief set communicator
  /// \param other_comm communication domain
  void setCommunicator(Communicator *other_comm) { commcator_ = other_comm; }

  /// \brief return communicator
  Communicator *getCommunicator() const { return commcator_; }

  /// \brief setup interfaces, must be called before first spmv
  /// \param startRow array that indicates matrix local row size
  /// length is nProcs+1
  void setUp(label *startRow);
};

}  // namespace UNAP
#endif  // INTERFACES_HPP
