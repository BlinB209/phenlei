#include "DataContainer.h"
using namespace std;

namespace PHSPACE
{

DataChar::DataChar()
{
    data = new vector<char>;
    pos  = 0;
}

DataChar::~DataChar()
{
    delete data;
}

void DataChar::MoveTo(CharVecSizeType pos)
{
    if ((0 < pos) && (pos < Size()))
    {
        this->pos = pos;
    }
    else if (0 == pos)
    {
        this->pos = pos;
    }
    else
    {
        ostringstream oss;
        oss << "Out of Range: pos  " << pos << endl;
        oss << "Out of Range: Size " << Size() << endl;

        //PHMPI::FreeBlockData();
        cout << oss.str() << endl;
        exit(0);
    }
}

void DataChar::ForwardPosition(CharVecSizeType size)
{
    pos += size;
}

CharVecSizeType DataChar::Size()
{
    return data->size();
}

char * DataChar::Begin()
{ 
    //! should consider that the size = 0;
    if (this->Size() == 0) return 0;
    return &((*data)[0]);
}

char * DataChar::DataPointer()
{ 
    //! should consider that the size = 0;
    if (this->Size() == 0) return 0;
    return &((*data)[pos]);
}

void DataChar::ToString(string &str)
{
    if (Size())
    {
        str.append(Begin(), Size());
    }
}

void DataChar::Write(void *data, CharVecSizeType size)
{
    if (size <= 0) return;

    memcpy(DataPointer(), data, size);
    ForwardPosition(size);
}

void DataChar::Read(void *data, CharVecSizeType size)
{
    if (size <= 0) return;
    memcpy(data, DataPointer(), size);
    ForwardPosition(size);
}

void DataChar::Write(void *data, CharVecSizeType size, int pos)
{
    MoveTo(pos);
    Write(data, size);
}

void DataChar::Read(void *data, CharVecSizeType size, int pos)
{
    MoveTo(pos);
    Read(data, size);
}

void DataChar::Resize(CharVecSizeType size)
{
    this->data->resize(size);
}

void DataChar::ReadFile(fstream &file)
{
    CharVecSizeType nlen = this->Size();
    if (nlen <= 0) return;

    char *data = this->Begin();

    file.read(data, nlen);
}

void DataChar::WriteFile(fstream &file)
{
    CharVecSizeType nlen = this->Size();
    if (nlen <= 0) return;

    char *data = this->Begin();

    file.write(data, nlen);
}

DataContainer::DataContainer(DataContainer *rhs)
{
    maxUnitSize = 102400000;
    start  = 0;
    finish = 1;
    pos = 0;
    elemPos = 0;
    data = new vector< DataChar * >;
    ContainerSizeType nelement = rhs->ElementSize();
    data->resize(nelement);
    for (ContainerSizeType i = 0; i < nelement; ++ i)
    {
        DataChar *rhsData = rhs->GetIter(i);
        CharVecSizeType nlen = rhsData->Size();

        DataChar *sectionData = new DataChar;
        sectionData->Resize(nlen);

        sectionData->MoveToBegin();
        rhsData->MoveToBegin();
        memcpy(sectionData->DataPointer(), rhsData->DataPointer(), nlen);
        (*data)[i] = sectionData;
    }
}

DataContainer::DataContainer()
{
    maxUnitSize = 102400000;
    start  = 0;
    finish = 1;
    data = new vector< DataChar * >;
    data->push_back(new DataChar());
    pos = 0;
    elemPos = 0;
}

DataContainer::~DataContainer()
{
    for (ContainerSizeType i = 0; i < data->size(); ++ i)
    {
        delete (*data)[i];
    }

    delete data;
}

ContainerSizeType DataContainer::ElementSize()
{
    return data->size();
}

DataChar * DataContainer::GetCurr()
{
    //! concern that the current element is self consistent with the pos.
    elemPos = pos / maxUnitSize;
    return (*data)[elemPos];
}

DataChar * DataContainer::GetIter(ContainerSizeType i)
{
    return (*data)[i];
}

void DataContainer::ForwardPosition(ContainerSizeType size)
{
    pos += size;
}

void DataContainer::BackwardPosition(ContainerSizeType size)
{
    elemPos = pos / maxUnitSize;

    CharVecSizeType iElem_pos = pos % maxUnitSize;
    if (iElem_pos == 0)
    {
        elemPos -= 1;
        iElem_pos = maxUnitSize;
    }

    if (iElem_pos > size)
    {
        GetIter(elemPos)->MoveTo(iElem_pos - size);
    }
    else
    {
        CharVecSizeType new_size = size - iElem_pos;
        GetIter(elemPos)->MoveToBegin();

        elemPos -= 1;
        if (new_size != 0)
        {
            GetIter(elemPos)->MoveTo(maxUnitSize - new_size);
        }
    }

    pos -= size;
}

void DataContainer::Write(DataContainer *dataContainer)
{
    dataContainer->MoveToBegin();
    CharVecSizeType length = dataContainer->Size();

    char *data;

    data = new char[length];

    dataContainer->Read(data, length);

    this->Write(data, length);

    delete [] data;
}

void DataContainer::Read(void *data, CharVecSizeType size)
{
    //! note. reading from DataContainer require the DataContainer is not empty.
    if (size <= 0) return;

    CharVecSizeType remainSize = Remain();

    if (remainSize >= size)
    {
        //! if the sub data length is not small than 'size', it can read directly.
        GetCurr()->Read(data, size);
        ForwardPosition(size);
    }
    else
    {
        //! read the data by segmentation.
        //! first read the remained data int this DataChar.
        GetCurr()->Read(data, remainSize);
        ForwardPosition(remainSize);

        void *newData = MovePointer(data, remainSize);
        CharVecSizeType newSize = size - remainSize;

        Read(newData, newSize);
    }
}

void DataContainer::Write(void *data, CharVecSizeType size)
{
    if (size <= 0) return;

    SecureRelativeSpace(size);

    CharVecSizeType remainSize = 0;
    CharVecSizeType newSize = size;
    void *newData = data;

    CharVecSizeType sizeSum = this->Size();

    while (sizeSum > this->pos)
    {
        remainSize = Remain();
        if (remainSize >= newSize)
        {
            GetCurr()->Write(newData, newSize);
            ForwardPosition(newSize);
        }
        else
        {
            GetCurr()->Write(newData, remainSize);
            ForwardPosition(remainSize);
            newData = MovePointer(newData, remainSize);
            newSize = newSize - remainSize;
        }
    }
}

void DataContainer::ReadString(string &cs)
{
    int nlen = 0;
    Read(&nlen, sizeof(int));

    char *data = new char[nlen+1];

    Read(data, nlen+1);

    cs = data;

    delete [] data;
}

void DataContainer::WriteString(string &cs)
{
    int nlen = static_cast<int>(cs.length());

    Write(&nlen, sizeof(int));

    char *data = new char[nlen+1];

    cs.copy(data, nlen);
    data[nlen] = '\0';

    Write(data, nlen+1);

    delete [] data;
}

void DataContainer::AppendString(string &cs)
{
    MoveToEnd();
    WriteString(cs);
}

void DataContainer::Write(ostringstream *oss)
{
    string str = oss->str();
    CharVecSizeType strSize = static_cast<CharVecSizeType>(str.size());
    Write(const_cast<char *>(str.c_str()), strSize * sizeof(char));
}

CharVecSizeType DataContainer::Size()
{
    CharVecSizeType sum = 0;
    for (unsigned int i = 0; i < ElementSize(); ++ i)
    {
        sum += GetIter(i)->Size();
    }
    return sum;
}

void DataContainer::Resize(CharVecSizeType nlen)
{
    if (nlen <= 0) return;

    //23 divided by 3 is 7, remainder 2.

    ContainerSizeType basicElem = nlen / maxUnitSize;
    ContainerSizeType remainder  = nlen % maxUnitSize;

    ContainerSizeType additionalElem = 0;
    if (remainder) additionalElem = 1;

    ContainerSizeType newSize = basicElem + additionalElem;
    ElementResize(newSize);

    for (ContainerSizeType i = 0; i < ElementSize(); ++ i)
    {
        CharVecSizeType needSize = maxUnitSize;
        if (i == basicElem)
        {
            needSize = remainder;
        }
        GetIter(i)->Resize(needSize);
    }
}

void DataContainer::ElementResize(ContainerSizeType newSize)
{
    if (newSize <= ElementSize())
    {
        Erase(newSize, ElementSize());
        data->resize(newSize);
    }
    else
    {
        ContainerSizeType ist = ElementSize();
        ContainerSizeType ied = newSize;

        for (ContainerSizeType i = ist; i != ied; ++ i)
        {
            data->push_back(new DataChar());
        }
    }
}

void DataContainer::SecureRelativeSpace(CharVecSizeType size)
{
    CharVecSizeType needSize = pos + size;

    SecureAbsoluteSpace(needSize);
}

void DataContainer::SecureAbsoluteSpace(CharVecSizeType needSize)
{
    //! if the space is enough, then it does not need to resize.
    //! it may cause some questions of these chars, if we do not resize, there may be some extra chars not been cleared in the memory.
    //if (needSize <= size()) return;

    Resize(needSize);
}

void DataContainer::MoveToBegin()
{
    pos = 0;
    for (ContainerSizeType i = 0; i < ElementSize(); ++ i)
    {
        GetIter(i)->MoveToBegin();
    }
}

void DataContainer::MoveToEnd()
{
    pos = Size();
    for (ContainerSizeType i = 0; i < ElementSize(); ++ i)
    {
        GetIter(i)->MoveToEnd();
    }
}

CharVecSizeType DataContainer::Remain()
{
    CharVecSizeType remainder = pos % maxUnitSize;
    return maxUnitSize - remainder;
}

void DataContainer::ReadFile(fstream &file)
{
    //! write the data to DataContainer from the file.
    streamsize nlen = 0;
    file.read(reinterpret_cast<char *>(&nlen), sizeof(streamsize));
    if (nlen <= 0) return;

    SecureAbsoluteSpace(static_cast< CharVecSizeType >(nlen));

    for (ContainerSizeType i = 0; i < ElementSize(); ++ i)
    {
        GetIter(i)->ReadFile(file);
    }
}

void DataContainer::WriteFile(fstream &file)
{
    streamsize nlen = Size();

    //! it need to be written to the file whenever the nlen is less than 0 or not.
    file.write(reinterpret_cast<char *>(&nlen), sizeof(streamsize));
    if (nlen <= 0) return;

    for (ContainerSizeType i = 0; i < ElementSize(); ++ i)
    {
        GetIter(i)->WriteFile(file);
    }
}

void DataContainer::ToString(string &str)
{
    for (ContainerSizeType i = 0; i < ElementSize(); ++ i)
    {
        GetIter(i)->ToString(str);
    }
}

void DataContainer::Append(void *data, CharVecSizeType size)
{
    MoveToEnd();
    CharVecSizeType needSize = this->Size() + size;
    SecureAbsoluteSpace(needSize);

    CharVecSizeType remainSize = 0;
    CharVecSizeType newSize = size;
    void *newData = data;

    CharVecSizeType sizeSum = this->Size();

    while (sizeSum > this->pos)
    {
        remainSize = Remain();
        if (remainSize >= newSize)
        {
            GetCurr()->Write(newData, newSize);
            ForwardPosition(newSize);
        }
        else
        {
            GetCurr()->Write(newData, remainSize);
            ForwardPosition(remainSize);
            newData = MovePointer(newData, remainSize);
            newSize = newSize - remainSize;
        }
    }
}

void DataContainer::Destroy(DataChar *iter)
{
    delete iter;
}

void DataContainer::Erase(ContainerSizeType ist, ContainerSizeType ied)
{
    for (ContainerSizeType i = ist; i != ied; ++ i)
    {
        Destroy(GetIter(i));
    }
}

//start to add by Xu Qingxin 2013.03.20., modified 2013.05.09.
CharVecSizeType DataContainer::MySize(streamsize size)
{
    CharVecSizeType sum = 0;
    ContainerSizeType elemBeg = static_cast< ContainerSizeType >(pos / maxUnitSize);
    ContainerSizeType elemEnd = static_cast< ContainerSizeType >((pos + size - 1) / maxUnitSize);

    for (ContainerSizeType i = elemBeg; i <= elemEnd; ++ i)
    {
        sum += GetIter(i)->Size();
    }

    return sum;
}

void DataContainer::InitBuffer(CharVecSizeType bufBeg, int elemBeg, int elemNum)
{
    pos = elemBeg * maxUnitSize + bufBeg;

    for (int i = elemBeg; i < elemBeg + elemNum; ++ i)
    {
        GetIter(i)->MoveToBegin();
    }
}

void DataContainer::InitBuffer(int elemNum, int ptrDim, 
     ContainerSizeType *ptrElemBeg, ContainerSizeType *ptrElemLen, 
     ContainerSizeType *ptrElemNum)
{
    ElementResize(elemNum);

    for (int i = 0; i < ptrDim; i ++)
    {
        CharVecSizeType remainder  = ptrElemLen[i] % maxUnitSize;
        ContainerSizeType elemEnd = ptrElemBeg[i] + ptrElemNum[i];

        for (ContainerSizeType inum = ptrElemBeg[i]; inum < elemEnd; inum ++)
        {
            CharVecSizeType needSize;
            if (inum == elemEnd - 1)
            {
                needSize = remainder;
            }
            else
            {
                needSize = maxUnitSize;
            }

            GetIter(inum)->Resize(needSize);
        }
    }

    MoveToBegin();
}

void DataContainer::StaticWrite(CharVecSizeType begin, void *data, CharVecSizeType size)
{
    //using namespace PHMPI;

    pos += begin;

    //check:
    if (size > this->MySize(size))
    {
        ostringstream oss;

        oss << "Error: Proc ? " //<< GetCurrentProcessorID() 
            << ", The write buffer size not enough!\n"
            << "   ---- begin =" << pos << ", actual size = " << size 
            << ", buffer size = " << this->MySize(size) << " ----\n"
            << " (In StaticWrite, from ActionKey.cpp)\n";

        string errmsg = oss.str();

        //PHMPI::FreeBlockData();
        cout << errmsg << endl;
        exit(0);
    }

    begin = pos % maxUnitSize;
    CharVecSizeType remainSize = Remain();

    if (remainSize >= size)
    {
        GetCurr()->Write(data, size, static_cast<int>(begin));
        ForwardPosition(size);
    }
    else 
    {
        //! write the data by segmentation.
        //! write some data for the rest space.
        GetCurr()->Write(data, remainSize, static_cast<int>(begin));

        ForwardPosition(remainSize);

        void *newData = MovePointer(data, remainSize);
        CharVecSizeType newSize = size - remainSize;

        //! repeat this function for the rest data that has node been written.
        StaticWrite(newData, newSize);
    }
}

void DataContainer::StaticWrite(void *data, CharVecSizeType size)
{
    //using namespace PHMPI;

    //check:
    if (size > this->MySize(size))
    {
        ostringstream oss;

        oss << "Error: Proc ? " //<< GetCurrentProcessorID() 
            << ", The write buffer size not enough!\n"
            << "   ---- pos =" << pos << ", actual size = " << size 
            << ", buffer size = " << this->MySize(size) << " ----\n"
            << " (In StaticWrite, from ActionKey.cpp)\n";

        string errmsg = oss.str();

        //PHMPI::FreeBlockData();
        cout << errmsg << endl;
        exit(0);
    }

    CharVecSizeType remainSize = Remain();

    if (remainSize >= size)
    {
        GetCurr()->Write(data, size);
        ForwardPosition(size);
    }
    else 
    {
        //! write the data by segmentation.
        //! write some data for the rest space.
        GetCurr()->Write(data, remainSize);

        ForwardPosition(remainSize);

        void *newData = MovePointer(data, remainSize);
        CharVecSizeType newSize = size - remainSize;

        //! repeat this function for the rest data that has node been written.
        StaticWrite(newData, newSize);
    }
}

void DataContainer::StaticRead(CharVecSizeType begin, void *data, CharVecSizeType size)
{
    //using namespace PHMPI;

    pos += begin;

    //check:
    if (size > this->MySize(size))
    {
        ostringstream oss;

        oss << "Error: Proc ? " //<< GetCurrentProcessorID() 
            << ", The read buffer size not enough!\n"
            << "   ---- begin =" << pos << ", actual size = " << size 
            << ", buffer size = " << this->MySize(size) << " ----\n"
            << " (In StaticRead, from ActionKey.cpp)\n";

        string errmsg = oss.str();

        //PHMPI::FreeBlockData();
        cout << errmsg << endl;
        exit(0);
    }

    CharVecSizeType remainSize = Remain();
    begin = pos % maxUnitSize;

    if (remainSize >= size)
    {
        //! if the sub data length is not small than 'size', it can read directly.
        GetCurr()->Read(data, size, static_cast<int>(begin));
        ForwardPosition(size);
    }
    else
    {
        //! read the data by segmentation.
        //! first read the remained data int this superchar.
        GetCurr()->Read(data, remainSize, static_cast<int>(begin));
        ForwardPosition(remainSize);

        void *newData = MovePointer(data, remainSize);
        CharVecSizeType newSize = size - remainSize;

        StaticRead(newData, newSize);
    }
}

void DataContainer::StaticRead(void *data, CharVecSizeType size)
{
    //using namespace PHMPI;

    streamsize needSize = pos + size;

    //check:
    if (needSize > static_cast<streamsize>(this->Size()))
    {
        ostringstream oss;

        oss << "Error: Proc ? " //<< GetCurrentProcessorID() 
            << ", The read buffer size not enough!\n"
            << "   ---- pos =" << pos << ", actual size = " << needSize 
            << ", buffer size = " << this->Size() << " ----\n"
            << " (In StaticRead, from ActionKey.cpp)\n";

        string errmsg = oss.str();

        //PHMPI::FreeBlockData();
        cout << errmsg << endl;
        exit(0);
    }

    CharVecSizeType remainSize = Remain();

    if (remainSize >= size)
    {
        //! if the sub data legth is not small than 'size', it can read directly.
        GetCurr()->Read(data, size);
        ForwardPosition(size);
    }
    else
    {
        //! read the data by segmentation.
        //! first read the remained data int this superchar.
        GetCurr()->Read(data, remainSize);
        ForwardPosition(remainSize);

        void *newData = MovePointer(data, remainSize);
        CharVecSizeType newSize = size - remainSize;

        StaticRead(newData, newSize);
    }
}

char * MovePointer(void *data, streamsize size)
{
    return reinterpret_cast<char *>(data) + size;
}
}
