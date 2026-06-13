inline SpecDiffHybGrid *SpecDiffHybSolver::GetGrid()
{
    return this->GridData;
}

inline void SpecDiffHybSolver::SetGrid(SpecDiffHybGrid *GridData)
{ 
    this->GridData = GridData;
}

inline RescaleField *SpecDiffHybSolver::GetRescaleFieldData()
{
    return this->RescaleFieldData;
}

inline void SpecDiffHybSolver::SetRescaleFieldData(RescaleField *RescaleFieldData)
{ 
    this->RescaleFieldData = RescaleFieldData;
}

inline Statistics *SpecDiffHybSolver::GetStatisticData()
{
    return this->StatisticData;
}

inline void SpecDiffHybSolver::SetStatisticData(Statistics *StatisticData)
{
    this->StatisticData = StatisticData;
}

inline LES *SpecDiffHybSolver::GetLESData()
{
    return this->LESData;
}

inline void SpecDiffHybSolver::SetLESData(LES *LESData)
{
    this->LESData = LESData;
}

