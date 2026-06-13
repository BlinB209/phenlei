//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Solver.h
//! @brief     Define the base solver class of all other solvers.
//! @author    Bell, He Xin, Xu Qingxin.

#pragma once
#include <sstream>
#include <vector>
#include "Geo_Grid.h"
#pragma warning(disable:4100)

using namespace std;

namespace PHSPACE
{
typedef enum
{
    NS_SOLVER,
    TURB_SOLVER,
    TRANSITION_SOLVER,
    INCOMPRESSIBLE_SOLVER
}solverEnum;

class PHNsTag
{
public:
    static int GetID() { return NS_SOLVER; }
};

class PHTurbulentTag
{
public:
    static int GetID() { return TURB_SOLVER; }
};

class PHTransitionTag
{
public:
    static int GetID() { return TRANSITION_SOLVER; }
};

class ActionKey;
class ActionTag;
class PHSolver;
class Post_Visual;
class Post_Probes;
class PHGeometry;
class Data_ParamFieldSuite;
class InterfaceDataProxy;
class InterpointDataProxy;

//! @brief The base class of all other solvers. Although different types of solver are create,
//!        they are all handled in the form of 'Solver' during simulation, this let them look same.
//!        The inherit order is: Solver->NSSolver->NSSolverStructFV/NSSolverStructFD/NSSolverUnstruct.
//!                                    ->TurbSolver->TurbSolverStr/TurbSolverUnstr.
//!                                    ->INCompSolver->INCompSolverStruct/INCompSolverUnstruct.
//!                                    ->UGKSSimulation, and other solvers.
//!        Note that, each zone can load several different solvers.
class PHSolver
{
private:
    //! The index of the solver.
    int index;

    //! The name of the solver.
    string solverName;

    //! The geometry, which is grid actually, bind to the solver.
    PHGeometry *geometry;

    //!
    int key; 

public:
    PHSolver(){};
    virtual ~PHSolver(){};

public:
    //! Assign the solver index, the default index is '0' if without given.
    void SetIndex(int index = 0);

    //! Assign the solver name, this function is for incompressible solver. 
    void SetName(string name);

    //! Assign the geometry(grid) to the solver.
    void SetGeometry(PHGeometry *geometry);

    //! Return the solver name,
    inline string GetName() const;

    //! Return the solver index, 
    int  GetIndex() const;

    //! Return the geometry(grid).
    PHGeometry *GetGeometry();
    
    //!
    void SetKey(int key);

    //!
    int  GetKey() const;
    
    //! Initialize the data, 
    //! here it is a pure virtual function, it should be re-defined in the real solver.
    virtual void InitMemory() = 0;

    virtual void ReleaseMemory() = 0;

    //! Initialize the flow, here it is a pure virtual function, 
    //! here it is a pure virtual function, it should be re-defined in the real solver.
    virtual void InitFlow() = 0;

    virtual void InitProtectedFlow() {};

    //! Interface data communications.    
    virtual void CommunicationInterfaceData() {};

    //! Interpoint data communications.
    virtual void CommunicationInterpointData() {};

    //! Overset data communications.
    virtual void CommunicationOversetSlipData(){};

    //! To get the number of equations to be solve in the solver.
    //! here it is a pure virtual function, it should be re-defined in the real solver.
    virtual int GetNumberOfEquations() = 0;

    //! For multi-grid technique, initialize the flow on the coarse grid.
    virtual void InitCoarseGridsFlow(){};

    //! The main entrance of solving the equation.
    //! here it is a pure virtual function, it should be re-defined in the real solver.
    virtual void Solve() = 0;

    //! Do some left work after solving the equation.
    //! here it is a pure virtual function, it should be re-defined in the real solver.
    virtual void Post() = 0;

    //! The following function is old version of task, should be removed latter.
    virtual void ActionReflect(ActionTag *acttag){};
    virtual void Action(ActionKey *actkey){};
    virtual void TranslateAction(ActionKey *actkey){};
    virtual streamsize TranslateActionLength(ActionKey *actkey){ return 0; };
    virtual void PrepareInterfaceData(Grid *grid, Data_ParamFieldSuite *datastore, InterfaceDataProxy *interfaceDataProxy){};
    virtual void ComputePostVisualVariables(Post_Visual * postVisualization) = 0;
    virtual void ComputePostProbesVariables(Post_Probes * postProbesVar){};
    virtual void AirForceCoefBodies(ActionKey *actkey){};
    virtual void RegisterOversetField(){};
    virtual void DeAllocateOversetInterfaceVar(Grid *grid){};
};
#include "Solver.hxx"

//! @brief GlobalSolvers stores the global solver list on each zone.
class GlobalSolvers
{
private:
    //! The global solvers list, in which the solver list is stored on each global zone.
    //! The first dimension is the zone index, and the second is the solver index.
    static vector< vector< PHSolver *> * > *globalSolversList; 
               
public:    
    //! Return the 'solverID'-th solver on 'zoneID'-th zone from global solvers list.
    //! @param[in] zoneID    The global zone index.
    //! @param[in] solverID  The solver index on the zone of 'zoneID'.
    static PHSolver *GetSolver(int zoneID, int solverID = 0);

    //! Return the total number of solvers on the zone of 'iZone'.
    //! @param[in] iZone    The zone index that you want to know the number of solvers.
    static int GetNumberOfSolvers(int iZone = 0);

    //! Add the solvers list onto the zone of 'iZone'.
    //! This function is used after all solvers have been created on the zone of 'iZone'.
    //! @param[in] iZone          The zone index that you want to add the solvers.
    //! @param[in] solversList    The solvers list that you want to add to the zone of 'iZone'.
    static void AddSolverToGlobal(int iZone, vector < PHSolver * > *solversList);
};

}
