#include "Task_Simple.h"

namespace PHSPACE
{
//Task_Simple::Task_Simple(){}
//Task_Simple::~Task_Simple(){}

void Task_Simple::PreTask(ActionKey *actkey){}

void Task_Simple::PostTask(ActionKey *actkey){}

void Task_Simple::MainTask(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int myid = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = GetZoneProcessorIDSepMode(iZone);

        if (myid == proc)
        {
            GeneralAction(actkey, iZone);
        }
    }
}

}