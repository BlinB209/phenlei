#include "IO_VirtualFile.h"

namespace PHSPACE
{

void VirtualFile::read(void *data, streamsize size)
{
    if (key == 0)
    {
        file->read(reinterpret_cast<char *>(data), size);
    }
    else
    {
        cdata->Read(data, static_cast<CharVecSizeType>(size));
    }
}

void VirtualFile::write(void *data, streamsize size)
{
    if (key == 0)
    {
        file->write(reinterpret_cast<char *>(data), size);
    }
    else
    {
        cdata->Write(data, static_cast<CharVecSizeType>(size));
    }
}

void VirtualFile::BeginReadWork()
{
    if (key == 0)
    {
        MarkSectionBegin();
        ReadDataLength();
    }
    else
    {
        cdata->MoveToBegin();
    }
}

void VirtualFile::EndReadWork()
{
    if (key == 0)
    {
        MarkSectionEnd();
    }
}

void VirtualFile::BeginWriteWork()
{
    if (key == 0)
    {
        MarkSectionBegin();
        ReservePlaceholder();
    }
}

void VirtualFile::EndWriteWork()
{
    if (key == 0)
    {
        MarkSectionEnd();
        MoveToSectionBegin();
        WriteDataLength();
        MoveToSectionEnd();
    }
}

void VirtualFile::ReadDataLength()
{
    read(&data_len, sizeof(streamsize));
}

void VirtualFile::ReservePlaceholder()
{
    write(&section_begin, sizeof(streamsize));
}

streamsize VirtualFile::GetCurrentPosition()
{
    if (key == 0)
    {
        return file->tellp();
    }
    else
    {
        return 0;
    }
}

void VirtualFile::MarkSectionBegin()
{
    this->section_begin = GetCurrentPosition();
}

void VirtualFile::MarkSectionEnd()
{
    this->section_end = GetCurrentPosition();
}

void VirtualFile::MoveToPosition(streamsize section_pos)
{
    file->seekp(section_pos);
}

void VirtualFile::MoveToSectionBegin()
{
    MoveToPosition(this->GetSectionBegin());
}

void VirtualFile::MoveToSectionEnd()
{
    MoveToPosition(this->GetSectionEnd());
}

void VirtualFile::WriteDataLength()
{
    data_len = GetSectionEnd() - GetSectionBegin() - sizeof(streamsize);

    write(&data_len, sizeof(streamsize));
}

}