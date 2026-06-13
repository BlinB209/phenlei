#ifndef READFROMHYPRE_HPP
#define READFROMHYPRE_HPP

#include <map>

#include "unapLduMatrix.hpp"

namespace UNAP {

/// \brief construct LduMatrix(inner) from Hypre output files
/// \param LduA LduMatrix
/// \param fileName file name
void constructLDUMatrixFromHypre(LduMatrix &LduA, const char *fileName);

/// \brief sort inter-faces as neighbor processors
/// \param val value in inter-faces
/// \param row row index of inter-faces
/// \param col column index of inter-faces
/// \param faceToProcNO neighbor processor ID of inter-faces
/// \param faceSize size of inter-faces
/// \param mapProcNO map: <key, value> - <processor ID, local processor ID
/// count(from 0)> \param neiProcNO vector of neighbor processor ID \param
/// procSize number of neighbor processors \param globalRowStart global index of
/// first row in current processor \param globalRowEnd global index of last row
/// in current processor
void sortInterFaces(scalarVector &val, labelVector &row, labelVector &col,
                    const labelVector &faceToProcNO, const label faceSize,
                    const labelVector &faceStart, Table<int, int> &mapProcNO,
                    const labelVector &neiProcNo, const label procSize,
                    const labelVector &globalRowStart,
                    const labelVector &globalRowEnd);

/// \brief construct LduMatrix(interfaces) from Hypre output files
/// \param LduA LduMatrix
/// \param nNeiProcs number of neighbor processors
/// \param destRank vector of neighbor processor ID, size is nNeiProcs
/// \param locPosition vector of accumulated number of inter-faces in one
/// processor \param faceCells local row index of inter-faces \param data values
/// of inter-faces
void constructLDUInterfacesFromHypre(LduMatrix &LduA, const label nNeiProcs,
                                     const labelVector &destRank,
                                     const labelVector &locPosition,
                                     const labelVector &faceCells,
                                     const scalarVector &data);

/// \brief construct vector from Hypre output files
/// \param b vector
/// \param fileName file name
void constructVectorFromHypre(scalarVector &b, const char *fileName);
}  // namespace UNAP
#endif
