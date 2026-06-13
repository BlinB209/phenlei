inline int Post_Probes::GetNumberofProbeVariables()
{
    return nProbeVariables;
}

inline bool Post_Probes::IsNeedMonitoring(int variableType)
{
    return probeVariables.find(variableType) != probeVariables.end();
}

inline string Post_Probes::GetProbesVariableName(int variableNumber) const
{
    map<int, string>::const_iterator it = probesVariablesMap.find(variableNumber);
    return it->second;
}

inline const set<int> & Post_Probes::GetProbeVariables() const
{
    return probeVariables;
}

inline void * Post_Probes::GetProbeVarPtr(const string &name)
{
    return probeVariablesPtr->GetDataPtr(name);
}

inline void Post_Probes::DeleteProbeVarPtr(const string &name)
{
    probeVariablesPtr->DeleteDataPtr(name);
}

inline Grid * Post_Probes::GetGrid() const
{
    return grid;
}