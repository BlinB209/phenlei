#ifndef PRINTUNAP_HPP
#define PRINTUNAP_HPP

#include <fstream>
#include <iomanip>
#include <sstream>

#include "unapLduMatrix.hpp"

#define FILEOPEN(fout, fileName, MYID)     \
  std::ostringstream os;                   \
  os << fileName << "_" << MYID << ".txt"; \
  fout.open(os.str().c_str())

#define FILECLOSE(fout, fileName) fout.close()

namespace UNAP {
/// \brief print LduMatrix in UNAP
/// \param A LduMatrix
/// \param name output file name
void printLDUMatrix(const LduMatrix &A, const char *name);

void printMTXMatrix(const LduMatrix &A, const char *name);

/// \brief print vector in UNAP
/// \tparam T vector type
/// \param b vector
/// \param name output file name
template <typename T>
void printVector(const T &b, const char *name);

template <typename T>
void printVector(const T &b, const char *name,
                 std::vector<label> distGlobalOrder);

/// \brief print processor interfaces of LduMatrix in UNAP
/// \param A LduMatrix
/// \param name output file name
void printInterfaces(const LduMatrix &A, const char *name);
}  // namespace UNAP

template <typename T>
void UNAP::printVector(const T &v, const char *fileName) {
  const label size = v.size();
  std::ofstream fout;

  FILEOPEN(fout, fileName, v.getCommunicator()->getMyId());

  fout << "rowSize: " << size << std::endl;
  forAll(i, size) {
    fout << std::setiosflags(std::ios::scientific) << std::setprecision(15)
         << v[i] << std::endl;
  }

  FILECLOSE(fout, fileName);
}
template <typename T>
void UNAP::printVector(const T &v, const char *fileName,
                       std::vector<label> distGlobalOrder) {
  const label size = v.size();
  if ((label)distGlobalOrder.size() != size) exit(1);
  std::ofstream fout;

  FILEOPEN(fout, fileName, v.getCommunicator()->getMyId());

  fout << "rowSize: " << size << std::endl;
  forAll(i, size) {
    fout << distGlobalOrder[i] << " " << std::setiosflags(std::ios::scientific)
         << std::setprecision(15) << v[i] << std::endl;
  }

  FILECLOSE(fout, fileName);
}

#endif  //- PRINTUNAP_HPP
