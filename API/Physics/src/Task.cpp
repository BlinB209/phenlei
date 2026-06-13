#include "Task.h"
#include "Solver.h"

#include "Math_BasisFunction.h"
#include "Geo_Grid.h"
#include "Constants.h"
using namespace std;

namespace PHSPACE
{
vector<Task *> * Task::ootasklist = new vector<Task *>;

// const int Task::SOLVE_FIELD      = 0;
// const int Task::CREATE_GRID      = 1;
// const int Task::CAL_WALL_DIST    = 2;
// const int Task::PARTITION_GRID   = 3;
// const int Task::OVERSETGRID_VIEW = 5;
// const int Task::OVERSET_CONFIG   = 6;

Task::Task()
{
    this->main_action = NOAction;
    this->pre_action  = NOAction;
    this->post_action = NOAction;
    this->actkey      = 0;
}

Task::~Task()
{
    delete this->actkey;
}

void Task::RunTask()
{
    ActionKey *actkey = this->GetActionKey();

    this->PreTask(actkey);

    this->MainTask(actkey);

    this->PostTask(actkey);
}

void NOAction (ActionKey *actkey)
{

}

void AddTasktoSimu(Task *task)
{
    Task::ootasklist->push_back(task);
}

void RunSimuTask(Task *task, bool run_now)
{
    if (run_now)
    {
        task->RunTask();
    }
    else
    {
        AddTasktoSimu(task);
    }
}

vector<RDouble> store_value;

RDouble MaxStoreValue()
{
    RDouble value = store_value[0];

    for (size_t i = 0; i < store_value.size(); ++ i)
    {
        value = MAX(store_value[i],value);
    }
    return value;
}

void GlobalSolve()
{
    uint_t taskNumber = Task::ootasklist->size();

    for (int iTask = 0; iTask < taskNumber; ++ iTask)
    {
        Task *task = (*Task::ootasklist)[iTask];

        task->RunTask();

        delete task;
    }

    ClearTaskList();
}

void ClearTaskList()
{
    Task::ootasklist->resize(0);
}

void GeneralAction(ActionKey *actkey, int iZone)
{
    if (actkey->kind == SOLVER_BASED)
    {
        PHSolver *solver = GlobalSolvers::GetSolver(iZone, actkey->solverID);
        solver->Action(actkey);
    }
    else
    {
        Grid *grid = GetGrid(iZone, actkey->level);
        grid->Action(actkey);
    }
}

void PH_Interface(ActionKey *actkey, int iZone, int jZone)
{
    using namespace PHMPI;

    //! iZone is the current zone, jZone is the target zone.
    int send_proc = GetZoneProcessorIDSepMode(iZone);
    int recv_proc = GetZoneProcessorIDSepMode(jZone);

    int tag = GetSendRecvTag(actkey, iZone);
    int myid = GetCurrentProcessorID();

    if (myid == send_proc)
    {
        //! If this process is a send process.
        actkey->ipos = jZone;
        GeneralAction(actkey, iZone);
    }

    
    PH_Trade(actkey, send_proc, recv_proc, tag);   

    if (myid == recv_proc)
    {
        //! If this process is a receive process,there are two situations:
        //! One is that this process is a receive process and also a send process.
        //! The other is that this process is a receive process not a send process.
        //! In both situation,the target zone data needs to be translate.
        actkey->ipos = iZone;
        GeneralTranslateAction(actkey, jZone);
    }
}

void PH_InterfaceNonBlocking(ActionKey *actkey, DataContainer *sendBuffer, DataContainer *receivedBuffer, int iZone, int jZone, vector <PH_Request> &requestContainer)
{
    using namespace PHMPI;

    int currentProcessor  = GetCurrentProcessorID();
    int sendProcessor     = GetZoneProcessorIDSepMode(iZone);
    int receiveProcessor  = GetZoneProcessorIDSepMode(jZone);

    int tag = GetSendRecvTag(actkey, iZone);

    //! Compress the send information into the actkey.
    if (currentProcessor == sendProcessor)
    {
        DataContainer *originalData = actkey->GetData();

        DataContainer *compressData = sendBuffer;
        if (sendProcessor == receiveProcessor)
        {
            compressData = receivedBuffer;
        }
        //! Fill the send buffer by actkey.
        actkey->SetData(compressData);
        actkey->ipos = jZone;
        GeneralAction(actkey, iZone);

        actkey->SetData(originalData);
    }

    if (sendProcessor == receiveProcessor)
    {
        return;
    }
    //! Communication.
    if (currentProcessor == sendProcessor)
    {
        streamsize nlen = sendBuffer->Size();
        //! Send the data to neighbors.
        send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
    }
    else if (currentProcessor == receiveProcessor)
    {
        //! Pre-calculating the data length of the received data.
        CharVecSizeType nlen = PH_Interface_GetLength(actkey, iZone, jZone);
        //! Receive data from neighbors.
        receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
    }
}

void GeneralTranslateAction(ActionKey *actkey, int iZone)
{
    if (actkey->kind == SOLVER_BASED)
    {
        GlobalSolvers::GetSolver(iZone, actkey->solverID)->TranslateAction(actkey);
    }
    else
    {
        GetGrid(iZone, actkey->level)->TranslateAction(actkey);
    }
}

streamsize GeneralTranslateActionLength(ActionKey *actkey, int iZone)
{
    if (actkey->kind == SOLVER_BASED)
    {
        return GlobalSolvers::GetSolver(iZone, actkey->solverID)->TranslateActionLength(actkey);
    }
    else
    {
        return GetGrid(iZone, actkey->level)->TranslateActionLength(actkey);
    }
}


CharVecSizeType PH_Interface_GetLength(ActionKey *actkey, int iZone, int jZone)
{
    using namespace PHMPI;
    CharVecSizeType nlen = 0;

    int send_proc = GetZoneProcessorIDSepMode(iZone);
    int recv_proc = GetZoneProcessorIDSepMode(jZone);

    int myid = GetCurrentProcessorID();

    if (myid == send_proc)
    {
        DataContainer *cdata = actkey->GetData();
        nlen = cdata->Size();
    }
    if (myid == recv_proc)
    {
        actkey->ipos = iZone;
        nlen = static_cast< CharVecSizeType >(GeneralTranslateActionLength(actkey, jZone));
    }
    return nlen;
}


void PH_ReadPattern(ActionKey *actkey, ActionKeyFUN main_action, int send_proc, int recv_proc, int iZone)
{
    using namespace PHMPI;

    int myid = GetCurrentProcessorID();
    int tag = GetSendRecvTag(actkey, iZone);

    if (myid == send_proc)
    {
        //! If this process is a send process.
        main_action(actkey);
    }

    PH_Trade(actkey, send_proc, recv_proc, tag);

    if (myid == recv_proc)
    {
        //! If this process is a receive process,there are two situations:
        //! One is that this process is a receive process and also a send process.
        //! The other is that this process is a receive process not a send process.
        //! In both situation,the target zone data needs to be translate.
        GeneralAction(actkey, iZone);
    }
}

void PH_WritePattern(ActionKey *actkey, ActionKeyFUN main_action, int send_proc, int recv_proc, int iZone)
{
    using namespace PHMPI;

    int myid = GetCurrentProcessorID();
    int tag = GetSendRecvTag(actkey, iZone);

    if (myid == send_proc)
    {
        //! If the process zone i belongs to is the current process.
        GeneralAction(actkey, iZone);
    }

    PH_Trade(actkey, send_proc, recv_proc, tag);

    if (myid == recv_proc)
    {
        //! If this process is server process.
        main_action(actkey);
    }
}

}
