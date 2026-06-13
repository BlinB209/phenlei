//! Read a value from file.
//! @param[in] file     file stream flow.
//! @param[in] value    variable name.
template < typename T >
LIB_EXPORT void PHRead(fstream *file, T &value)
{
    file->read(reinterpret_cast<char *>(&value), sizeof(T));
}

//! Read a value from file.
//! @param[in] file     file stream flow.
//! @param[in] value    variable name.
template < typename T >
LIB_EXPORT void PHRead(fstream &file, T &value)
{
    PHRead(&file, value);
}

//! Read a field from file.
//! @param[in] file     file stream flow.
//! @param[in] field    variable field.
//! @param[in] numberOfElements    number of field elements.
template < typename T >
LIB_EXPORT void PHRead(fstream *file, T *field, int numberOfElements)
{
    file->read(reinterpret_cast<char *>(field), numberOfElements * sizeof(T));
}

//! Read a field from file.
//! @param[in] file     file stream flow.
//! @param[in] field    variable field.
//! @param[in] numberOfElements    number of field elements.
template < typename T >
LIB_EXPORT void PHRead(fstream &file, T *field, int numberOfElements)
{
    PHRead(&file, field, numberOfElements);
}

//! Read a field from file.
//! @param[in] file     file stream flow.
//! @param[in] field    variable field.
//! @param[in] numberOfElements    number of field elements.
template < typename T >
LIB_EXPORT void PHRead(fstream *file, vector< T > &field, int numberOfElements)
{
    PHRead(*file, field, numberOfElements);
}

//! Read a field from file.
//! @param[in] file     file stream flow.
//! @param[in] field    variable field.
//! @param[in] numberOfElements    number of field elements.
template < typename T >
LIB_EXPORT void PHRead(fstream &file, vector< T > &field, int numberOfElements)
{
    PHRead(&file, &field[0], numberOfElements);
}

//! Read a value from virtual file.
//! @param[in] virtualFile    virtual file.
//! @param[in] value          variable name.
template < typename T >
LIB_EXPORT void PHRead(VirtualFile *virtualFile, T &value)
{
    virtualFile->read(&value, sizeof(T));
}

//! Read a field from virtual file.
//! @param[in] virtualFile    virtual file.
//! @param[in] field          variable field.
//! @param[in] numberOfElements    number of field elements.
template < typename T >
LIB_EXPORT void PHRead(VirtualFile *virtualFile, T *field, int numberOfElements)
{
    virtualFile->read(field, numberOfElements * sizeof(T));
}

template < typename T >
void PHSkipRead(fstream *file, int numberOfElements)
{
    int bufferSize = sizeof(T) * numberOfElements;
    file->seekg(bufferSize, ios::cur);
}

template < typename T >
void PHSkipRead(fstream &file, int numberOfElements)
{
    PHSkipRead< T >(&file, numberOfElements);
}

//! Write a value to file.
//! @param[in] file     file stream flow.
//! @param[in] value    variable name.
template < typename T >
LIB_EXPORT void PHWrite(fstream *file, T &value)
{
    file->write(reinterpret_cast<char *>(&value), sizeof(T));
}

//! Write a value to file.
//! @param[in] file     file stream flow.
//! @param[in] value    variable name.
template < typename T >
LIB_EXPORT void PHWrite(fstream &file, T &value)
{
    PHWrite(&file, value);
}

//! Write a field to file.
//! @param[in] file     file stream flow.
//! @param[in] field    variable field.
//! @param[in] numberOfElements    number of field elements.
template < typename T >
LIB_EXPORT void PHWrite(fstream *file, T *field, int numberOfElements)
{
    file->write(reinterpret_cast<char *>(field), numberOfElements * sizeof(T));
}

//! Write a field to file.
//! @param[in] file     file stream flow.
//! @param[in] field    variable field.
//! @param[in] numberOfElements    number of field elements.
template < typename T >
LIB_EXPORT void PHWrite(fstream &file, T *field, int numberOfElements)
{
    PHWrite(&file, field, numberOfElements);
}

//! Write a field to file.
//! @param[in] file     file stream flow.
//! @param[in] field    variable field.
//! @param[in] numberOfElements    number of field elements.
template < typename T >
LIB_EXPORT void PHWrite(fstream &file, vector< T > &field, int numberOfElements)
{
    PHWrite(&file, &field[0], numberOfElements);
}

//! Write a value to virtual file.
//! @param[in] virtualFile    virtual file.
//! @param[in] value          variable name.
template < typename T >
LIB_EXPORT void PHWrite(VirtualFile *virtualFile, T &value)
{
    virtualFile->write(&value, sizeof(T));
}

//! Write a field to virtual file.
//! @param[in] virtualFile    virtual file.
//! @param[in] field          variable field.
//! @param[in] numberOfElements    number of field elements.
template < typename T >
LIB_EXPORT void PHWrite(VirtualFile *virtualFile, T *field, int numberOfElements)
{
    virtualFile->write(field, numberOfElements * sizeof(T));
}