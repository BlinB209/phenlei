inline void PHSolver::SetKey(int key) 
{ 
    this->key = key;
}

inline int PHSolver::GetKey() const 
{ 
    return key;
}

inline void PHSolver::SetName(string name)
{
    this->solverName = name;
}

inline string PHSolver::GetName() const
{
    return this->solverName;
}

inline void PHSolver::SetIndex(int index) 
{ 
    this->index = index; 
}

inline int PHSolver::GetIndex() const 
{ 
    return index; 
}

inline PHGeometry * PHSolver::GetGeometry() 
{ 
    return geometry; 
}

inline void PHSolver::SetGeometry(PHGeometry *geometry)
{ 
    this->geometry = geometry;
}