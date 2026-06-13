#include "Task_WriteBinaryTecplot.h"
#include "PHIO.h"

#include "Post_Visual.h"
using namespace std;
#pragma warning(disable:4100)
namespace PHSPACE
{


void Task_WriteBinaryTecplot::PostTask(ActionKey *actkey)
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if (actkey->VTKvisual == true && myid == server)
    {
        ActionKeyFUN post_action = this->GetPostAction();
        post_action(actkey);
    }

    if (datalist.size() > 0)
    {
        ActionKeyFUN post_action = this->GetPostAction();
        post_action(actkey);

        #ifdef USE_TecplotLib
        DumpToTecio dumptotecio(datalist, actkey->tecfilename, this->flowType);
        dumptotecio.Run();
        #endif
    }
}

LIB_EXPORT Task_WriteBinaryTecplot::Task_WriteBinaryTecplot(int flowType)
{
    this->flowType = flowType;
}

LIB_EXPORT Task_WriteBinaryTecplot::~Task_WriteBinaryTecplot(){};

void Task_WriteBinaryTecplot::PreTask(ActionKey *actkey){};

void Task_WriteBinaryTecplot::MainTask(ActionKey *actkey)
{
    using namespace PHMPI;
    ServerWrite(actkey);
}

void Task_WriteBinaryTecplot::ServerWrite(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    
    int recv_proc;
    int myid;

    if (actkey->VTKvisual == true)
    {
        if (! PHMPI::IsParallelRun())
        {
            ParallelOpenFile(actkey);
        }
        else
        {
            ParallelOpenFile2(actkey);
        }
    }

    for (int iZone =  0; iZone < nZones; ++ iZone)
    {
        DataContainer *cdata = actkey->GetData();
        cdata->MoveToBegin();
        
        DataContainer *tecdata = actkey->GetTecData();
        tecdata->MoveToBegin();

        int send_proc = GetZoneProcessorID(iZone);
        recv_proc = GetZoneFileID(iZone);

        myid = GetCurrentProcessorID();
        int tag = GetSendRecvTag(actkey, iZone);

        if (myid == send_proc)
        {
            GeneralAction(actkey, iZone);
        }
        
        PH_Trade(actkey, send_proc, recv_proc, tag);

        if (myid == recv_proc)
        {
            if (actkey->VTKvisual == true)
            {
                WriteASCIIFile(actkey);
            }

            datalist.push_back(tecdata);
            DataContainer *cdata1 = new DataContainer();
            actkey->SetTecData(cdata1);
        }
    }

    if (actkey->VTKvisual == true)
    {
        ParallelCloseFile(actkey);
    }
}

}
