LIB_EXPORT inline std::size_t IO_GridReader::GetNumberOfGridFile() const
{
	return gridGroupNameList.size();
}

LIB_EXPORT inline int IO_GridReader::GetNumberOfGrids()
{
	return gridGroup->GetNumberofGrid();
}

LIB_EXPORT inline int IO_GridReader::GetNumberofGlobalZones()
{
	return gridGroup->GetNZones();
}

LIB_EXPORT inline Grid * IO_GridReader::GetGrid(int iGrid)
{
	return gridGroup->GetGrid(iGrid);
}

LIB_EXPORT inline GridGroup * IO_GridReader::GetGridGroup()
{
	return gridGroup;
}

LIB_EXPORT inline vector< GridGroup * > & IO_GridReader::GetGridGroups()
{
	return gridGroups;
}