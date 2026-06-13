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
//! @file      IO_VirtualFile.h
//! @brief     Define the class VirtualFile.
//! @author    Refactored by He Xianyao, Bell.

#pragma once
#include <fstream>
#include "DataContainer.h"

using namespace std;

namespace PHSPACE
{

//! @brief The virtual file class.
class VirtualFile
{
private:
    int key;
    fstream *file;
    DataContainer *cdata;

public:
    VirtualFile(fstream *file) { this->file = file; key = 0; }
    VirtualFile(DataContainer *cdata) { this->cdata = cdata; key = 1; }

public:
    //! Read data from file or datacontainer.
    //! @param[in] data    data address.
    //! @param[in] size    data size.
    void read(void *data, streamsize size);

    //! Write data to file or datacontainer.
    //! @param[in] data    data address.
    //! @param[in] size    data size.
    void write(void *data, streamsize size);

    //! Begin to read.
    void BeginReadWork();

    //! Mark section end only for key = 0.
    void EndReadWork();

    //! Mark and write section_begin only for key = 0.
    void BeginWriteWork();
    
    //! Do something to end writing.
    void EndWriteWork();

    //! Read data length.
    void ReadDataLength();
    
    //! Write section_begin.
    void ReservePlaceholder();
    
    //! Return 0 or file pointer location.
    streamsize GetCurrentPosition();

private:
    //! section_begin assignment.
    void MarkSectionBegin();

    //! section_end assignment.
    void MarkSectionEnd();
    
    //! Return section_begin.
    streamsize GetSectionBegin();
    
    //! Return section_end.
    streamsize GetSectionEnd();
    
    //! Move to some position.
    void MoveToPosition(streamsize section_pos);
    
    //! Move to the begin of section.
    void MoveToSectionBegin();
    
    //! Move to the end of section.
    void MoveToSectionEnd();
    
    //! Calculate the data_len and write it.
    void WriteDataLength();

private:
    streamsize section_begin, section_end, data_len;
};

#include "IO_VirtualFile.hxx"
}