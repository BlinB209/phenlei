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
//! @file      DataContainer.h
//! @brief     It defines kernel data structure of PHengLEI, named 'DataContainer'.
//!            All types of data are standardized into DataChar. "Three-In-One" data is designed:
//!            1 - MPI communication; 2 - File I/O; 3 - Communication for hybrid parallel.
//! @author    Bell, Xu Qingxin, He Xin.

#pragma once
#include <fstream>
#include <string>
#include "TypeDefine.h"
#include "AMRDef.h"
using namespace std;

namespace PHSPACE
{
class DataChar;
typedef vector< char >::size_type CharVecSizeType;
typedef vector< DataChar * >::size_type ContainerSizeType;
typedef string::size_type StringSizeType;
//! the purpose of designing DataContainer is to make up the deficiencies of the original DataChar.
//! the interface should look like to be the same as the DataChar.
//! the design concept is as this at least.
class DataContainer
{
public:
    typedef int iterator;
public:
    DataContainer(DataContainer *rhs);    // Bell add, but i am not sure it's perfect right.
    DataContainer();
    ~DataContainer();
private:
    vector< DataChar * > *data;
    iterator start;
    iterator finish;
    CharVecSizeType pos, elemPos, maxUnitSize;

public:

    //! get data from vector *data.
    DataChar * GetCurr();

    //! return vector *data's value at the site of i.
    DataChar * GetIter(ContainerSizeType i);

    //! return start.
    iterator Begin() { return start; }

    //! return finish.
    iterator End() { return finish; };

    //! invoking delete function to delete *iter.
    void Destroy(DataChar *iter);

    //!invoking destroy function to delete *iter data,start with ist and end with ied.
    //!@param[in ]: ist  begin position.
    //!@param[in ]: ied  end position
    void Erase(ContainerSizeType ist, ContainerSizeType ied);

    //! return the number of vector *data 
    ContainerSizeType ElementSize();

    //!read DataContainer vector<char> * data,and receive it by void *data.
    //!@param[in ]: *data , target object.
    //!@param[in ]: size
    void Read(void *data, CharVecSizeType size);

    //!write *data into DataContainer vector<char> *data
    //!@param[in ]: *data ,source object.
    //!@param[in ]: size
    void Write(void *data, CharVecSizeType size);

    //!initialize dataContainer with char 
    void Write(DataContainer *dataContainer);

    //!read data from file.
    void ReadFile(fstream &file);

    //!write data into file.
    void WriteFile(fstream &file);

    //!read DataContainer vector<char> *data, and put it in cs.
    void ReadString(string &cs);

    //!write cs into DataContainer vector<char> *data.
    void WriteString(string &cs);

    //!put oss data into vector<char> *data.
    void Write(ostringstream *oss);

    //!return size of sum of all data object.
    CharVecSizeType Size();

    //!vector<char> *data object resize itself.
    void Resize(CharVecSizeType nlen);

    //!char transform to string
    void ToString(string &str);

    //!link another data
    void Append(void *data, CharVecSizeType size);

    //!link sting object
    void AppendString(string &cs);

    //!before invoking resize function, calculate it size.
    void SecureRelativeSpace(CharVecSizeType size);

    //! invoking resize function.
    void SecureAbsoluteSpace(CharVecSizeType needSize);

    //! invoking moveto function.
    void MoveToBegin();

    //! invoking moveto function.
    void MoveToEnd();

    //! calculate: pos += size;
    void ForwardPosition(ContainerSizeType size);

    void BackwardPosition(ContainerSizeType size);

private:

    //!use new DataChar object fill up ,until the size equal new_size.
    void ElementResize(ContainerSizeType newSize);

    //!return max_unit_size - remainder;
    CharVecSizeType Remain();

public:

    //!invoking DataChar::read function
    //! @param[in ]:   *data :  target object.\n
    //! @param[in ]:
    void StaticRead(void *data, CharVecSizeType size);

    //!invoking DataChar::write function
    //! @param[in ]: *data :source object.\n
    //! @param[in ]:
    void StaticWrite(void *data, CharVecSizeType size);


    //!invoking DataChar::read function with two parameter
    //! @param[in ]:   *data :  target object.\n
    //! @param[in ]:
    //! @param[in ]:
    void StaticRead(CharVecSizeType begin, void *data, CharVecSizeType size);

    //!invoking DataChar::write function with two parameter
    //! @param[in ]:  *data  :source object.\n
    //! @param[in ]:
    //! @param[in ]:
    void StaticWrite(CharVecSizeType begin, void *data, CharVecSizeType size);
public:

    //!return size of sum of all data object
    CharVecSizeType MySize(streamsize size);

    //!
    streamsize GetRemainder(streamsize position) 
        const { return position % maxUnitSize; }

    //!
    int GetElementIndex(streamsize position) 
        const { return (int)(position / maxUnitSize); }

    //!
    streamsize GetElementPos(streamsize elem) 
        const { return elem * maxUnitSize; }
public:

    //!from elem_beg to elem_num.use  moveto(0) function initial vector<char> * data .
    //! @param[in ]: bufBeg   : \n
    //! @param[in ]: elemBeg  : the number of invoking moveto_begin() function.\n
    //! @param[in ]: elemNum  : the number of invoking moveto_begin() function.
    void InitBuffer(CharVecSizeType bufBeg, int elemBeg, int elemNum);

    //!
    //! @param[in ]: ptrElemLen  : use to resize vector<char> *data.\n
    //! @param[in ]: ptrElemBeg  : use to store which object to be resize.\n
    //! @param[in ]: ptrElemNum  : use to store which object to be resize.\n
    //! @param[in ]: elemNum
    //! @param[in ]: ptrDim
    void InitBuffer(int elemNum, int ptrDim, ContainerSizeType *ptrElemBeg, 
         ContainerSizeType *ptrElemLen, ContainerSizeType *ptrElemNum);
};

class DataChar
{
public:
    DataChar();
    ~DataChar();
public:

    //!
    CharVecSizeType Size();

    //! invoking function memcpy
    //! @param[in ]:   *data : target object.\n
    //! @param[in ]:        
    void Read(void *data, CharVecSizeType size);

    //! invoking function read with two parameter
    //! @param[in ]:   *data :  target object.\n
    //! @param[in ]:   pos:     moveto(pos);\n
    //! @param[in ]:
    void Read(void *data, CharVecSizeType size, int pos);

    //! invoking function memcpy
    //! @param[in ]:   *data :source object.\n
    //! @param[in ]:            
    void Write(void *data, CharVecSizeType size);

    //! invoking function write with two parameter
    //! @param[in ]:   *data:    source object.\n
    //! @param[in ]:   pos:      moveto(pos);\n
    //! @param[in ]:
    void Write(void *data, CharVecSizeType size, int pos);

    //!file.read(data, nlen);
    void ReadFile(fstream &file);

    //!file.write(data, nlen);
    void WriteFile(fstream &file);

    //!str.append(begin(),size());
    void ToString(string &str);

    //!return &((*data)[0]);
    char * Begin();

    //!
    void MoveToBegin() { MoveTo(0); }

    //!
    void MoveToEnd() { pos = Size(); }

    //!invoking data->resize(size);
    void Resize(CharVecSizeType size);

    //!return &((*data)[pos]);
    char * DataPointer();

    //! calculate pos += size;
    void ForwardPosition(CharVecSizeType size);

    void MoveTo(CharVecSizeType pos);

private:
    CharVecSizeType pos;
    vector<char> *data;

public:

    //!
    char * DataPointer(int beg) { return &((*data)[beg]); }
};

//!return reinterpret_cast<char *>(data) + size;
char * MovePointer(void *data, streamsize size);
}
