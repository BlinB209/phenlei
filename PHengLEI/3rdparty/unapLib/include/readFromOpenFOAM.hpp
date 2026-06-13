#ifndef READFROMOPENFOAM_HPP
#define READFROMOPENFOAM_HPP

#include "unapLduMatrix.hpp"

namespace UNAP {
/// \brief construct LduMatrix(inner) from OpenFOAM output files
/// \param LduA LduMatrix
/// \param fileName file name
void constructLDUMatrixFromOpenFOAM(LduMatrix &LduA, const char *fileName);

/// \brief construct vector from OpenFOAM output files
/// \param b vector
/// \param fileName file name
void constructVectorFromOpenFOAM(scalarVector &b, const char *fileName);

/// \brief construct LduMatrix(interfaces) from OpenFOAM output files
/// \param LduA LduMatrix
/// \param fileName file name
void constructLDUInterfacesFromOpenFOAM(LduMatrix &LduA, const char *fileName);

}  // namespace UNAP
#endif
