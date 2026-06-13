inline int CFDSolver::GetNPostSolve() 
{
    return 1; 
}

inline int CFDSolver::GetZoneLocalID(void)
{
    int g = 0;
    Grid * grid = this->GetGeometry()->GetGrid(g);
    return grid->GetZoneLocalID();
}

inline uint_t CFDSolver::GetNumberOfMultiGrid()
{
    return this->GetGeometry()->GetNumberOfMultiGrid(); 
}

inline Grid *CFDSolver::GetGrid(int level)
{
    return this->GetGeometry()->GetGrid(level);
}

inline string CFDSolver::CastName(const string &name)
{
    return name; 
}

inline int CFDSolver::GetNumberOfEquations() 
{
    return 1; 
}

inline void CFDSolver::SetPostVisualization(Post_Visual *postVisualization)
{
    this->postVisualization = postVisualization;
}

inline Post_Visual* CFDSolver::GetPostVisualization()
{
    return this->postVisualization;
}

inline void CFDSolver::SetPostVisualizationWall(Post_VisualWall *postVisualization_in)
{
    this->postVisualWall = postVisualization_in;
}

inline Post_VisualWall* CFDSolver::GetPostVisualizationWall()
{
    return this->postVisualWall;
}

