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
//! @file      Task.h
//! @brief     Task defines several design patterns, such as communication pattern,
//!            read and write pattern, and so on. It is not recommended for the learners 
//!            in the first stage.
//! @author    He Xin.
#pragma once
#include "PHMpi.h"
using namespace std;

namespace PHSPACE
{

typedef void (*ServerFUN)(ActionKey *actkey);
typedef void (*ActionKeyFUN)(ActionKey *actkey);
void NOAction (ActionKey *actkey);

//! @brief The base class of all the design patterns
class LIB_EXPORT Task
{
private:
    ActionKeyFUN main_action, pre_action, post_action;
    ActionKey *actkey;

public:
    static vector<Task *> *ootasklist;
//  static const int SOLVE_FIELD     ;
//  static const int CREATE_GRID     ;
//  static const int CAL_WALL_DIST   ;
//  static const int PARTITION_GRID  ;
//  static const int OVERSETGRID_VIEW;
//  static const int OVERSET_CONFIG;

public:
    Task();
    virtual~Task();

public:
    //! Return the main_action
    //! @return main_action
    ActionKeyFUN GetMainAction();

    //! Return the pre_action
    //! @return pre_action
    ActionKeyFUN GetPreAction ();

    //! Return the post_action
    //! @return post_action
    ActionKeyFUN GetPostAction();

    //! Set main_action with action
    //! @param[in] action    An input action
    void SetMainAction(ActionKeyFUN action);

    //! Set pre_action with action
    //! @param[in] action    An input action
    void SetPreAction (ActionKeyFUN action);

    //! Set post_action with action
    //! @param[in] action    An input action
    void SetPostAction(ActionKeyFUN action);

    //! Return the actkey
    //! @return actkey
    ActionKey * GetActionKey();

    //! Set this->actkey with actkey
    //! @param[in] action    An input actkey
    void SetActionKey(ActionKey *actkey);

    //! Run all the task(PreTask, MainTask, PostTask)
    void RunTask();

private:
    virtual void PreTask(ActionKey *actkey) = 0;
    virtual void PostTask(ActionKey *actkey) = 0;
    virtual void MainTask(ActionKey *actkey) = 0;
};

void RunSimuTask(Task *task, bool run_now);

void GlobalSolve();
void AddTasktoSimu(Task *task);
void ClearTaskList();

void PH_Interface(ActionKey *actkey, int iZone, int jZone);
void PH_InterfaceNonBlocking(ActionKey *actkey, DataContainer * sendBuffer, DataContainer * receivedBuffer, int iZone, int jZone, vector< PH_Request > & requestContainer);
void GeneralTranslateAction(ActionKey *actkey, int iZone);
streamsize GeneralTranslateActionLength(ActionKey *actkey, int iZone);
CharVecSizeType PH_Interface_GetLength(ActionKey *actkey, int iZone, int jZone);

void PH_ReadPattern (ActionKey *actkey, ActionKeyFUN main_action, int send_proc, int recv_proc, int iZone);
void PH_WritePattern(ActionKey *actkey, ActionKeyFUN main_action, int send_proc, int recv_proc, int iZone);

void GeneralAction(ActionKey *actkey, int iZone);

#include "Task.hxx"

}