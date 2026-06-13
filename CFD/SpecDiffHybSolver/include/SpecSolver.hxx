inline SpecGrid *SpecSolver::GetGrid()
{
    return this->GridData;
}

inline void SpecSolver::SetGrid(SpecGrid *GridData)
{ 
    this->GridData = GridData;
}

inline Spectrum *SpecSolver::GetSpectrumData()
{
    return this->SpectrumData;
}

inline void SpecSolver::SetSpectrumData(Spectrum *SpectrumData)
{
    this->SpectrumData = SpectrumData;
}

inline RescaleFieldHIT *SpecSolver::GetRescaleFieldData()
{
    return this->RescaleFieldData;
}

inline void SpecSolver::SetRescaleFieldData(RescaleFieldHIT *RescaleFieldData)
{ 
    this->RescaleFieldData = RescaleFieldData;
}

inline StatisticsHIT *SpecSolver::GetStatisticData()
{
    return this->StatisticData;
}

inline void SpecSolver::SetStatisticData(StatisticsHIT *StatisticData)
{
    this->StatisticData = StatisticData;
}

inline LESHIT *SpecSolver::GetLESData()
{
    return this->LESData;
}

inline void SpecSolver::SetLESData(LESHIT *LESData)
{
    this->LESData = LESData;
}
