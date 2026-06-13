inline int Post_Visual::GetNumberofVisualVariables()
{
    return nVisualVariables;
}

inline bool Post_Visual::IsNeedVisualization(int variableType)
{
    return visualVariables.find(variableType) != visualVariables.end();
}

inline string Post_Visual::GetVariableName(int variableNumber) const
{
    map<int, string>::const_iterator it = visualVariablesMap.find(variableNumber);
    return it->second;
}

inline const set<int> &Post_Visual::GetVisualVariables() const
{
    return visualVariables;
}

inline void * Post_Visual::GetVisualNodeVarPtr(const string &name)
{
    return visualVariablesPtr->GetDataPtr(name);
}

inline void Post_Visual::DeleteVisualNodeVarPtr(const string &name)
{
    visualVariablesPtr->DeleteDataPtr(name);
}

inline Grid * Post_Visual::GetGrid() const
{
    return grid;
}

inline string Post_VisualWall::GetVariableName(int variableNumber) const
{
    map<int, string>::const_iterator it = visualVariablesMap.find(variableNumber);
    return it->second;
}