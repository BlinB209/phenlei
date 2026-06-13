#ifndef SWAgglomeration_HPP
#define SWAgglomeration_HPP

#include "swRestInterStruct.h"
#include "unapAgglomeration.hpp"
#include "unapMatrix.hpp"

namespace UNAP {
/// \brief strore the struct of accessing address for slave cores
/// used in restriction and interpolation in multigrid
class SwRestInterMap {
 private:
  const Agglomeration &aggl_;

 public:
  friend class Agglomeration;

  SwRestInterMap(const Agglomeration &aggl) : aggl_(aggl) {}

  static label bandSize_;
  static label minRowSizeUsingSW_;

  /// \brief restrict
  static bool *restFirstUse_;
  static restStruct *restStructLevels_;

  /// \brief face restrict
  static bool *faceRestFirstUse_;
  static restStruct *faceRestStructLevels_;

  /// \brief interpolate
  static bool *interFirstUse_;
  static interStruct *interStructLevels_;

  /// \brief agglomerate matrix upper
  static bool *aggMatrixUpperFirstUse_;
  static aggMatrixUpperStruct *aggMatrixUpperStructLevels_;

  void initRestInterSize();

  void initRestStruct(Vector<scalar> &cf, const Vector<scalar> &ff,
                      const label fineLevelIndex);

  void initInterStruct(Vector<scalar> &ff, const Vector<scalar> &cf,
                       const label levelIndex);

  void agglomerateMatrixUpper(scalarVector &coarseUpper,
                              scalarVector &coarseDiag,
                              const scalarVector &fineUpper,
                              const label fineLevelIndex);

  void initAggMatrixUpperStruct(scalarVector &coarseUpper,
                                scalarVector &coarseDiag,
                                const scalarVector &fineUpper,
                                const label &fineLevelIndex);
};

}  // namespace UNAP

#endif
