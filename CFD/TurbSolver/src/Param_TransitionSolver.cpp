#include "Param_TransitionSolver.h"
#include "Transition.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#pragma warning(disable:6386)

namespace PHSPACE
{

LIB_EXPORT Param_TransitionSolver::~Param_TransitionSolver()
{
    if (freeStreamTransitionVar != NULL)
    {
        delete [] freeStreamTransitionVar;
        freeStreamTransitionVar = NULL;
    }
}

LIB_EXPORT void Param_TransitionSolver::Init()
{
    Param_CFDSolver::Init();

    transitionCFLScale = GlobalDataBase::GetDoubleParaFromDB("turbCFLScale");

    transitionType = GlobalDataBase::GetIntParaFromDB("transitionType");

    nTransitionEquation = GlobalDataBase::GetIntParaFromDB("n_transition");

    freeStreamViscosity = GlobalDataBase::GetDoubleParaFromDB("freeStreamViscosity");

    using namespace IDX;
    RDouble *prim_inf = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble reference_density_farfield = prim_inf[IR];

    string viscousName;
    GlobalDataBase::GetData("viscousName", &viscousName, PHSTRING, 1);

    RDouble refReNumber = 1.0e5;
    GlobalDataBase::GetData("refReNumber", &refReNumber, PHDOUBLE, 1);

    RDouble turb_cmu=1.0;
    if (transitionType == IREGAMA)
    {
        ce1 = 1.0;
        ca1 = 2.0;
        ce2 = 50.0;
        ca2 = 0.06;
        dct = 2.0;
        df = 1.0;
        cct = 0.03;
        s1 = 2.0;
        ccf = 0.6;
    }

    freeStreamTransitionVar = new RDouble[nTransitionEquation];

    RDouble TUoo = 0.18;
    RDouble Fcta = 1.0;
    freeStreamTransitionVar[IGAMA] = 1.0;

    RDouble mutoo = GlobalDataBase::GetDoubleParaFromDB("freeStreamViscosity");
    RDouble kwoo = GlobalDataBase::GetDoubleParaFromDB("kwoo") / refReNumber;
    RDouble keoo = mutoo * kwoo / (turb_cmu * reference_density_farfield);

    turbulenceIntensity = GlobalDataBase::GetDoubleParaFromDB("turbIntensity");

    if (turbulenceIntensity >= SMALL)
    {
        keoo = 1.5 * (0.01 * turbulenceIntensity) * (0.01 * turbulenceIntensity);
    }

    TUoo = TurbulenceIntensity(1.0, keoo);
    freeStreamTransitionVar[IRECT] = EmpiricalCorrelationOfRectat(TUoo, Fcta);

}

}