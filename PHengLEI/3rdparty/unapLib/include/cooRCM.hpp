/// \brief C++ program for Implementation of Reverse Cuthill Mckee Algorithm
#ifndef cooRCM_HPP_
#define cooRCM_HPP_

#include "utilityContainer.hpp"

using UTILITY::Array;

namespace COORCM {
class ReorderingSSM {
 private:
  /// \brief number of non-zeros
  const label _nnz;

  /// \brief number of rows
  const label _nRows;

  /// \brief array of rows
  const label *_row;

  /// \brief array of columns
  const label *_col;

  /// \brief number of columns in each row
  static Array<label> _globalDegree;

  /// \brief vertex order after reordering
  Array<label> *_vetexOrderPtr;

  /// \brief edge order after reordering
  Array<label> *_edgeOrderPtr;

  // label* _rowsOffset;

  void degreeGenerator();

  static bool compareDegree(label i, label j);

  void CuthillMckee();

  // Implementation of reverse Cuthill-Mckee algorithm
  void ReverseCuthillMckee();

 public:
  /// \brief constructor
  /// \param nnz number of non zeros
  /// \param nRows number of rows
  /// \param row row index array
  /// \param col column index array
  ReorderingSSM(const label nnz, const label nRows, const label *row,
                const label *col);

  ~ReorderingSSM();

  /// \brief retrun array of vertex order
  Array<label> *getVetexOrder();

  /// \brief return array of edge order
  Array<label> *getEdgeOrder();
};

}  // namespace COORCM

// Driver Code
#endif
