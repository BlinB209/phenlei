#ifndef MULTIGRID_AGGLOMERATION_DIRECTIONALAGGLOMERATION_HPP
#define MULTIGRID_AGGLOMERATION_DIRECTIONALAGGLOMERATION_HPP
/*! \file directionalAgglomeration.hpp
 *  \author Zhao Chengpeng (chengpeng_zhao@foxmail.com)
 *  \date 2022-11-02
 *  \brief an agglomeration method, usually better than LduAgglomeration and
 * MetisAgglomeration.
 *  ref: <<Directional agglomeration multigrid techniques for
 * high-Reynolds-number viscous flows>>
 */

#include "lduAgglomeration.hpp"

namespace UNAP {
/// \brief abstract base-class for matrix Agglomeration in LDU type
class DirectionalAgglomeration : public LduAgglomeration {
 private:
  /// \brief calculate and return Agglomeration of given level
  /// \param nCoarseCells number of rows in coarse level
  /// \param fineA fine level matrix
  /// \param weights vecttor of weights to agglomerate matrix
  virtual labelVector &agglomerate(label &nCoarseCells, const LduMatrix &fineA,
                                   label level) override;

 public:
  /// \brief constructors
  /// \param A input matrix to produce coarse level matrix
  /// \param rowSizeInCoarsestLevel
  DirectionalAgglomeration(const LduMatrix &A,
                           label rowSizeInCoarsestLevel = 50);

  /// \brief agglomerate all levels starting from default face weights
  virtual void agglomerate() override;
};

}  // namespace UNAP

#endif /* end of include guard: \
          MULTIGRID_AGGLOMERATION_DIRECTIONALAGGLOMERATION_HPP */
