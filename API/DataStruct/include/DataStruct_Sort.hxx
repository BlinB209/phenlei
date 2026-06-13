template<typename T>
DataStruct_Sort<T>::DataStruct_Sort()
{

}

template<typename T>
DataStruct_Sort<T>::DataStruct_Sort(T &value, int index)
{
    this->value = value; 
    this->index = index;
}

template<typename T>
bool DataStruct_Sort<T>::operator < (const DataStruct_Sort &rhs) const
{
    return value < rhs.value;
}

template<typename T>
bool DataStruct_Sort<T>::operator > (const DataStruct_Sort &rhs) const
{
    return value > rhs.value;
}