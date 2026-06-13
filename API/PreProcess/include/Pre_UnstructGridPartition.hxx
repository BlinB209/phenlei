inline Grid ** Pre_UnstructGridPartition::GetPartitionGrids()
{
    return this->grids;
}

inline void Pre_UnstructGridPartition::SetParallelPartitionMethod(int parallelPartitionMethod)
{
    this->parallelPartitionMethod = parallelPartitionMethod;
}

inline bool Pre_UnstructGridPartition::IsBelongtoThisZone(int faceMark)
{
    if (faceMark == NOT_IN_THIS_ZONE)
    {
        return false;
    }
    else if (faceMark == IN_THIS_ZONE)
    {
        return true;
    }
    else
    {
        TK_Exit::UnexpectedVarValue("faceMark = ", faceMark);
        return true;
    }
}