#ifndef UTILITY_WRITETOFILE_HPP
#define UTILITY_WRITETOFILE_HPP

#include "utilityCommunicator.hpp"

namespace UTILITY {

/// \brief add a processor id after file name, e.g. fieldA_1.txt
/// \param comm communicator
/// \param fileName fileName will be edited with a processor ID
void addFileNameWithProcID(Communicator *comm, char *fileName);

/// \brief write label arrays into files, number of arrays is not limited
/// \param comm communicator
/// \param fileName output file name
/// \param n number of arrays
/// \param len size of an array
/// \param ... pointer of array, whose size equals to number of arrays
void writeToFileLabel(Communicator *comm, const char *fileName, int n, int len,
                      ...);

/// \brief write scalar arrays into files, number of arrays is not limited
/// \param comm communicator
/// \param fileName output file name
/// \param n number of arrays
/// \param len size of an array
/// \param ... pointer of array, whose size equals to number of arrays
void writeToFileScalar(Communicator *comm, const char *fileName, int n, int len,
                       ...);

}  // namespace UTILITY
#endif
