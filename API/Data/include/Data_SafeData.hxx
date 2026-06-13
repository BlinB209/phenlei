namespace PHSPACE
{

inline string Data_SafeData::GetName() const
{
    return name;
}

inline int Data_SafeData::GetType() const
{
    return type;
}

inline int Data_SafeData::GetSize() const
{
    return size;
}

inline void * Data_SafeData::GetData() const
{
    return data;
}

}