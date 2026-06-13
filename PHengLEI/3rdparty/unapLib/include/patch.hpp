#ifndef PATCH_HPP
#define PATCH_HPP
/*! \file patch.hpp
 *  \modified 2022-11-14
 */

#include "unapVector.hpp"

namespace UNAP {
/// \brief store the inter-faces between two neighboring processors
class Patch {
 private:
  /// \brief number of faces in this patch
  label size_;

  /// \brief current processor number
  label myProcNo_;

  /// \brief neighbor processor number
  label neighbProcNo_;

  /// \brief face-cell addressing, local row index
  mutable labelVector *faceCellsPtr_;

  /// \brief face-cell2 addressing, global col index
  mutable labelVector *faceCells2Ptr_;

  /// \brief face-cell2 addressing, relative col index
  mutable labelVector *faceCells2OffsetPtr_;

  /// \brief coefficients
  mutable scalarVector *patchCoeffsPtr_;

  /// \brief face restriction addressing if needed
  mutable labelVector *faceRestrictAddressingPtr_;

  /// \brief row index size without duplicate
  label sendIndSize_;

  /// \brief row index size without duplicate
  label recvIndSize_;

  /// \brief row index array without duplicate, build from faceCellsPtr_
  /// length is sendIndSize_
  label *rowindices_;

 public:
  /// \brief constructor
  /// \param size size of the patch
  /// \param myProcNo current processor ID
  /// \param neighbProcNo the neighbor processor ID of the patch
  Patch(label size, label myProcNo, label neighbProcNo);

  virtual ~Patch();

  /// \brief return neighbor processor number
  inline label neighbProcNo() const { return neighbProcNo_; }

  /// \brief set the neighbor processor ID
  /// \param i the neighbor processor ID
  inline void neighbProcNo(const label i) { neighbProcNo_ = i; }

  /// \brief return current processor number
  inline label myProcNo() const { return myProcNo_; }

  /// \brief set the current processor ID
  /// \param i the current processor ID
  inline void myProcNo(const label i) { myProcNo_ = i; }

  /// \brief return faceCells
  inline labelVector &faceCells() const {
    CHECK_POINTER(faceCellsPtr_)
    return *faceCellsPtr_;
  }
  /// \brief return faceCells2
  inline labelVector &faceCells2() const {
    CHECK_POINTER(faceCells2Ptr_)
    return *faceCells2Ptr_;
  }
  /// \brief return faceCells2Offset
  inline labelVector &faceCells2Offset() const {
    CHECK_POINTER(faceCells2OffsetPtr_)
    return *faceCells2OffsetPtr_;
  }

  /// \brief set faceCells
  /// \param a the new faceCells vector
  inline void faceCells(labelVector &a) const {
    DELETE_POINTER(faceCellsPtr_)
    faceCellsPtr_ = &a;
  }
  /// \brief set faceCells2
  /// \param a the new faceCells2 vector
  inline void faceCells2(labelVector &a) const {
    DELETE_POINTER(faceCells2Ptr_)
    faceCells2Ptr_ = &a;
  }
  /// \brief set faceCells2Offset
  /// \param a the new faceCells2Offset vector
  inline void faceCells2Offset(labelVector &a) const {
    DELETE_POINTER(faceCells2OffsetPtr_)
    faceCells2OffsetPtr_ = &a;
  }

  /// \brief return patch coefficients
  inline scalar patchCoeffs(const label faceI) const {
    return (*patchCoeffsPtr_)[faceI];
  }

  /// \brief return the data in patch
  inline scalarVector &patchCoeffs() const {
    CHECK_POINTER(patchCoeffsPtr_)
    return *patchCoeffsPtr_;
  }

  /// \brief set the data in patch
  inline void patchCoeffs(scalarVector &a) {
    DELETE_POINTER(patchCoeffsPtr_)
    patchCoeffsPtr_ = &a;
  }

  /// \brief return the size of patch data
  inline label size() const { return size_; }

  /// \brief set the size of patch data
  inline void size(const label i) { size_ = i; }

  /// \brief set the face restriction addressing
  inline void faceRestrictAddressing(labelVector &a) const {
    DELETE_POINTER(faceRestrictAddressingPtr_)
    faceRestrictAddressingPtr_ = &a;
  }

  /// \brief return the face restriction addressing
  inline labelVector &faceRestrictAddressing() const {
    CHECK_POINTER(faceRestrictAddressingPtr_)
    return *faceRestrictAddressingPtr_;
  }

  /// \brief renumber the face-cell addressing using cellMap
  /// \param cellMap new cell map
  void reorderPatchFaceCells(const label *cellMap);

  /// \brief return the array size to send by this patch
  label sendIndSize() const { return sendIndSize_; }
  /// \brief return the array size need to recv by this patch
  label recvIndSize() const { return recvIndSize_; }

  /// \brief return pointer to row index array without duplicate
  label *rowindices() const { return rowindices_; }
  label rowindices(label i) const { return rowindices_[i]; }

  /// \brief setup the patch, must be called before use
  void setUp();
};

}  // namespace UNAP
#endif  // PATCH_HPP
