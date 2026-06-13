inline RDouble Param_TurbSolver::GetTurbCFLScale() const
{
    return turbCFLScale;
}

inline int Param_TurbSolver::GetDESType() const
{
    return DESType;
}

inline int Param_TurbSolver::GetTransitionType() const
{
    return transitionType;
}

inline int Param_TurbSolver::GetSAProductType() const
{
    return SAProductType;
}

inline int Param_TurbSolver::GetSSTProductType() const
{
    return SSTProductType;
}

inline int Param_TurbSolver::GetNTurbulenceEquation() const
{
    return nTurbulenceEquation;
}

inline RDouble Param_TurbSolver::GetEddyViscosityLimit() const
{
    return eddyViscosityLimit;
}

inline int Param_TurbSolver::GetWallFunctionType() const
{
    return wallFunctionType;
}

inline RDouble Param_TurbSolver::GetSST_a1() const
{
    return SST_a1;
}

inline RDouble Param_TurbSolver::GetSST_alphaw1() const
{
    return SST_alphaw1;
}

inline RDouble Param_TurbSolver::GetSST_alphaw2() const
{
    return SST_alphaw2;
}

inline RDouble Param_TurbSolver::GetSST_beta1() const
{
    return SST_beta1;
}

inline RDouble Param_TurbSolver::GetSST_beta2() const
{
    return SST_beta2;
}

inline RDouble Param_TurbSolver::GetSST_betaStar() const
{
    return SST_betaStar;
}

inline RDouble Param_TurbSolver::GetSA_cv1_cube() const
{
    return SA_cv1_cube;
}

inline RDouble Param_TurbSolver::GetSSTProductLimit() const
{
    return SSTProductLimit;
}

inline RDouble Param_TurbSolver::GetKW_sigma() const
{
    return KW_sigma;
}

inline RDouble Param_TurbSolver::GetKW_sigmaW() const
{
    return KW_sigmaW;
}

inline RDouble Param_TurbSolver::GetKW_sigmaW1() const
{
    return KW_sigmaW1;
}

inline RDouble Param_TurbSolver::GetKW_sigmaW2() const
{
    return KW_sigmaW2;
}

inline RDouble Param_TurbSolver::GetKW_sigmaK() const
{
    return KW_sigmaK;
}

inline RDouble Param_TurbSolver::GetKW_sigmaK1() const
{
    return KW_sigmaK1;
}

inline RDouble Param_TurbSolver::GetKW_sigmaK2() const
{
    return KW_sigmaK2;
}

inline RDouble Param_TurbSolver::GetFreeStreamViscosity() const
{
    return freeStreamViscosity;
}

inline RDouble * Param_TurbSolver::GetFreeStreamTurbVar() const
{
    return freeStreamTurbVar;
}

inline int Param_TurbSolver::GetKindOfTurbSource() const
{
    return kindOfTurbSource;
}

inline int Param_TurbSolver::GetWennSchemeFlag() const
{
    return isWennScheme;
}

inline int Param_TurbSolver::UsingConservativeForm() const
{
    return rhoFlag;
}
inline RDouble Param_TurbSolver::GetInviscidSpectrumRadiusCoef() const
{
    return inviscidSpectrumRadiusCoef;
}
