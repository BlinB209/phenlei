#ifndef MATRIXCONVERSION_HPP
#define MATRIXCONVERSION_HPP

#include "unapLduMatrix.hpp"

namespace UNAP {

/// \brief sort data form giving order
/// \param data
/// \param order giving order
/// \param cellFaces arrray of number of faces for each cell
/// \param other_comm communication domain
void sortData(scalarVector &data, const labelVector &order,
              const labelVector &cellFaces, Communicator *other_comm);

/// \brief reorder COO matrix, to make it ordered in columns in every row
/// to use this, the input COO matrix must be ordered in rows
/// row, col must be started from 0
/// \param dataPtr pointer of coefficients array
/// \param rowsPtr row index
/// \param columnPtr column index
/// \param rowSize size of rows
/// \param size size of non-zeros
/// \param other_comm communication domain
void reorderCOO(scalar *dataPtr, label *rowsPtr, label *columnPtr,
                const label rowSize, const label size,
                Communicator *other_comm);

/// \brief reorder values by newOrder
/// \param val pointer of values array
/// \param newOrder new ordering of values
/// \param size size of non-zeros
/// \param other_comm communication domain
void reorderValue(scalar *val, const label *newOrder, const label size,
                  Communicator *other_comm);

/// \brief reorder upper part as row-wise
/// \param rowsPtr row index
/// \param columnPtr column index
/// \param rowSize size of rows
/// \param size size of non-zeros
/// \param newOrder position index in the new upper array
/// \param other_comm communication domain
void reorderUFace(label *rowsPtr, label *columnPtr, const label rowSize,
                  const label size, label *newOrder, Communicator *other_comm);

/// \brief reorder lower part as row-wise
/// \param rowsPtr row index
/// \param columnPtr column index
/// \param rowSize size of rows
/// \param size size of non-zeros
/// \param newOrder position index in the new lower array
/// \param other_comm communication domain
void reorderLFace(label *rowsPtr, label *columnPtr, const label rowSize,
                  const label size, label *newOrder, Communicator *other_comm);

/// \brief convert COO matrix to LDU matrix
/// only for diagonal(inner) part
/// \param dataPtr values array
/// \param rowsPtr row index array
/// \param columnPtr column index array
/// \param rowSize number of rows
/// \param size number of non-zeros
/// \param symm if matrix is symmetric or not
/// \param other_comm communication domain
/// \return return a LduMatrix
LduMatrix &coo2Ldu(const scalar *dataPtr, const label *rowsPtr,
                   const label *columnPtr, const label rowSize,
                   const label size,
                   const bool symm,  // symm refers to the data
                   Communicator *other_comm);

/// \brief convert CSR matrix to LDU matrix
/// only for diagonal(inner) part
/// \param dataPtr values array
/// \param compRowsPtr array of compressed row
/// \param columnPtr column index array
/// \param rowSize number of rows
/// \param size number of non-zeros
/// \param symm if matrix is symmetric or not
/// \param other_comm communication domain
/// \return return a LduMatrix
LduMatrix &csr2Ldu(const scalar *dataPtr, const label *compRowsPtr,
                   const label *columnPtr, const label rowSize,
                   const label size,
                   const bool symm,  // symm refers to the data
                   Communicator *other_comm);

}  // namespace UNAP

#endif  // MATRIXCONVERSION_HPP
