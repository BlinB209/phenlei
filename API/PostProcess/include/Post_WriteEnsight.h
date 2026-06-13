//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Post_WriteEnsight.h
//! @brief     Write flow field into ensight file.
//! @author    Xu Gang, Yan Yong(External contributor).

#pragma once
#include "Post_WriteVisualFile.h"

namespace PHSPACE
{

//! File acess method
typedef enum
{
    STDIO_SERIAL,
    MPI_PARALLEL_IO,
} file_access_t;


//! @brief Post_WriteEnsight class realize ensight visual file output function.\n
class Post_WriteEnsight : public Post_WriteVisualFile
{
public:
    Post_WriteEnsight();
    ~Post_WriteEnsight();

private:
    void Initialize();

    //! Dump data into visual file.
    void WriteFile();

    //! Write geometry file.
    void WriteGeometryFile();
    void WriteGeometryFileSerial();
    void WriteGeometryFileParallel();

    //! Write variable file.
    void WriteFieldData();
    void WriteFieldDataSerial();
    void WriteFieldDataParallel();

    //! Write and update Ensight case file
    void UpdateCaseFile();

    //! Clear data.
    void ClearFieldData();

private:
    //! Open file.
    void GetFileHandle(const char *fileName);

    //! Close file.
    void ClearFileHandle();

    //! Write string data to file.
    void WriteString(const char *str);

    //! Write int data to file
    void WriteInt(int num);

    //! Write string data to file use mpi.
    void MpiWriteString(const char *str, MPI_Offset offset);

    //! Write int data to file use mpi.
    void MpiWriteInt(int num, MPI_Offset offset);

private:
    //! Ensight file type.
    //!     0 -- binary.
    //!     1 -- ASCII.
    int fileType;

    //! Visual data location.
    //!     0 -- On node, default.
    //!     1 -- On cell center.
    bool dataOnElement;

    //! If preserve visual file each iter.
    //!     0 -- Not, only preserve file last iter.
    //!     1 -- Yes, for instance, unsteady simulation.
    bool writeFieldEachIter;

    //! Current dump iter.
    int    currentIter;
    string currentIterStr;

    //! Geometry file name.
    string meshFileName;

    //! File write access type.
    //!     0 -- serial.
    //!     1 -- parallel.
    file_access_t fileMethod;

    //! File handle, when serial write.
    FILE *serialFilehandle;

    //! MPI file handle.
    MPI_File mpiFilehandle;

    //! MPI file offset.
    MPI_Offset fileOffset;

};



}
