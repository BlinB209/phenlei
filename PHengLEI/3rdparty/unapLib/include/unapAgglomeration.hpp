#ifndef UNAP_AGGLOMERATION_HPP
#define UNAP_AGGLOMERATION_HPP

#include "unapMatrix.hpp"

// #define AMG_RCM

namespace UNAP {
/// \brief abstract base-class for matrix Agglomeration
class Agglomeration {
 protected:
  /// \brief communication domain
  Communicator *commcator_;
  bool noCoarse_;

 public:
  /// \brief constructor
  /// \param other_comm communication domain
  Agglomeration(Communicator *other_comm)
      : commcator_(other_comm), noCoarse_(false) {}

  bool ifNoCoarse() const { return noCoarse_; }

  /// \brief destructor
  virtual ~Agglomeration() {}

  /// \brief number of coarse matrix levels
  virtual label size() const = 0;

  /// \brief return coarse matrix of given level
  /// \param levlei coarse level index
  virtual const Matrix &coarseMatrixLevels(const label leveli) const = 0;

  /// \brief return cell restrict addressing of given level
  /// \param levlei coarse level index
  virtual const labelVector &restrictAddressing(const label leveli) const = 0;

  /// \brief return face restrict addressing of given level
  /// \param levlei coarse level index
  virtual const labelVector &faceRestrictAddressing(
      const label leveli) const = 0;

  /// \brief agglomerate coarse matrix
  virtual void agglomerateMatrix(const label fineLevelIndex) = 0;

  /// \brief calculate the coarse level matrix using given weights
  /// \param weights the given weights
  virtual void agglomerate() = 0;

  /// \brief set the row size in the coarsest level
  /// \param i row size in the coarsest level
  virtual void SET_rowSizeInCoarsestLevel(const label i) = 0;

  /// \brief set the maximum number of coarse levels
  /// \param i the maximum number of coarse levels
  virtual void SET_maxLevels(const label i) = 0;

  // /// \brief temporary
  // virtual void AgglomerationReorderTopo() = 0;

  /// \brief restrict (integrate by summation) cell field
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

  /// \brief restrict (integrate by summation) face field
  /// \tparam T vector type
  /// \param cf coarse level vector
  /// \param ff fine level vector
  /// \param fineLevelIndex fine level index
  template <typename T>
  void restrictFaceField(Vector<T> &cf, const Vector<T> &ff,
                         const label fineLevelIndex) const;
};

template <typename T>
void Agglomeration::restrictField(Vector<T> &cf, const Vector<T> &ff,
                                  const label fineLevelIndex) const {
  const labelVector &fineToCoarse = restrictAddressing(fineLevelIndex);

#ifdef DEBUG
  if (ff.size() != fineToCoarse.size()) {
    commcator_->log()
        << "Error in restrictField: field does not correspond to level "
        << fineLevelIndex << " sizes: field = " << ff.size()
        << " level = " << fineToCoarse.size() << ENDL;
    ERROR_EXIT;
  }
#endif

  cf.SET_zero();
  label len = ff.size();

  forAll(i, len) { cf[fineToCoarse[i]] += ff[i]; }
}

template <typename T>
void Agglomeration::prolongField(Vector<T> &ff, const Vector<T> &cf,
                                 const label coarseLevelIndex) const {
  const labelVector &fineToCoarse = restrictAddressing(coarseLevelIndex);

#ifdef DEBUG
  if (ff.size() != fineToCoarse.size()) {
    commcator_->log()
        << "Error in prolongField: field does not correspond to level "
        << coarseLevelIndex << " sizes: field = " << ff.size()
        << " level = " << fineToCoarse.size() << ENDL;
    ERROR_EXIT;
  }
#endif

  label len = ff.size();

  forAll(i, len) { ff[i] = cf[fineToCoarse[i]]; }
}

template <typename T>
void Agglomeration::restrictFaceField(Vector<T> &cf, const Vector<T> &ff,
                                      const label fineLevelIndex) const {
  const labelVector &fineToCoarse = faceRestrictAddressing(fineLevelIndex);

#ifdef DEBUG
  if (ff.size() != fineToCoarse.size()) {
    commcator_->log()
        << "Error in restrictFaceField: field does not correspond to level "
        << fineLevelIndex << " sizes: field = " << ff.size()
        << " level = " << fineToCoarse.size() << ENDL;
    ERROR_EXIT;
  }
#endif

  cf.SET_zero();
  label len = ff.size();

  forAll(i, len) {
    label cFace = fineToCoarse[i];

    if (cFace >= 0) {
      cf[cFace] += ff[i];
    }
  }
}

template <typename T>
int sign(T val) {
  return (T(0) < val) - (val < T(0));
}

}  // namespace UNAP

#endif
