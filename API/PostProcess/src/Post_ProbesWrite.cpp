#include "Post_ProbesWrite.h"
#include "IO_FileName.h"
#include "TK_Parse.h"
#include "PHIO.h"
#include "HyList.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"
#include "Geo_UnstructBC.h"
#include "Geo_StructBC.h"
#include "Post_Probes.h"
#include "Pointer.h"
#include "Solver.h"
#include "Geo_UnstructGrid.h"
#include "TK_Log.h"

namespace PHSPACE
{

LIB_EXPORT Post_ProbesWrite::Post_ProbesWrite(int flowType)
{
    this->flowType = flowType;
}

LIB_EXPORT Post_ProbesWrite::~Post_ProbesWrite()
{

}

LIB_EXPORT void Post_ProbesWrite::Run()
{
    string samplefile = "sample.dat";
    GlobalDataBase::GetData("samplefile", &samplefile, PHSTRING, 1);
    ActionKey *actkey = new ActionKey();
    actkey->level    = 0;
    actkey->filename = samplefile;
    actkey->openmode = ios_base::out|ios_base::app;

    ServerProbesWrite(actkey);
    PostProbesWrite(actkey);
    delete actkey;    actkey = nullptr;
    this->Clear();
}

void Post_ProbesWrite::ServerProbesWrite(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = PHMPI::GetNumberofGlobalZones();

    int recv_proc;
    int myid;

    for (int iZone =  0; iZone < nZones; ++ iZone)
    {
        DataContainer *cdata = actkey->GetData();
        cdata->MoveToBegin();

        int send_proc = GetZoneProcessorID(iZone);
        recv_proc = GetServerProcessorID();

        myid = GetCurrentProcessorID();
        int tag = GetSendRecvTag(actkey, iZone);

        if (myid == send_proc)
        {
            Grid *grid = GetGrid(iZone, actkey->level);
            int nProbeVariables = GlobalDataBase::GetIntParaFromDB("nProbeVariables");
            int probeVariables[100] = {0};

            GlobalDataBase::GetData("probeVariables", probeVariables, PHINT, nProbeVariables);
            Post_Probes *postProbesVar = new Post_Probes(nProbeVariables, probeVariables, grid);
            int Nsolver = GlobalSolvers::GetNumberOfSolvers(iZone);
            for (int iSolver = 0; iSolver < Nsolver; ++ iSolver)
            {
                GlobalSolvers::GetSolver(iZone, iSolver)->ComputePostProbesVariables(postProbesVar);
            }

            StoreProbesVarariables(actkey, postProbesVar);
            FreePointer(postProbesVar);
        }

        PH_Trade(actkey, send_proc, recv_proc, tag);

        if (myid == recv_proc)
        {
            datalist.push_back(cdata);
            DataContainer *cdata1 = new DataContainer();
            actkey->SetData(cdata1);
        }
    }
}

void Post_ProbesWrite::StoreProbesVarariables(ActionKey *actkey, Post_Probes *postProbesVar)
{
    Grid *grid_in = postProbesVar->GetGrid();

    int nProbeVariables = postProbesVar->GetNumberofProbeVariables();
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nProbeVariables = nProbeVariables + 2 * numberOfSpecies;
    }
    RDouble **q_probes = new RDouble *[nProbeVariables];
    postProbesVar->GetAllProbesVarPtr(q_probes);

    if (actkey->filename != "")
    {
        int probesLocalNumber = grid_in->GetZoneProbesNumber();
        PHWrite(actkey->GetData(), &probesLocalNumber, 1);
        vector<int> probesGlobalID = grid_in->GetZoneProbesGlobalID();
        vector<int> probesLineID = grid_in->GetZoneProbesLineID();
        vector<int> probesSurfaceID = grid_in->GetZoneProbesSurfaceID();
        vector<vector<RDouble> > probesCoordinates = grid_in->GetZoneProbesCoordinates();
        int nDim;

        for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
        {
            nDim = static_cast<int>(probesCoordinates[iProbe].size());

            PHWrite(actkey->GetData(), &nDim, 1);
            PHWrite(actkey->GetData(), probesGlobalID[iProbe]);
            PHWrite(actkey->GetData(), probesLineID[iProbe]);
            PHWrite(actkey->GetData(), probesSurfaceID[iProbe]);
            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHWrite(actkey->GetData(), probesCoordinates[iProbe][iDim]);
            }
            for (int iVariable = 0; iVariable < nProbeVariables; ++ iVariable)
            {
                PHWrite(actkey->GetData(), q_probes[iVariable][iProbe]);
            }
        }
    }
    delete [] q_probes;    q_probes = NULL;
}

void Post_ProbesWrite::PostProbesWrite(ActionKey *actkey)
{
    if (datalist.size() > 0)
    {
        DumpToFileASCII dumpToFileASCII(datalist, actkey->filename);
        dumpToFileASCII.Run();
    }

    for (unsigned int iZone = 0; iZone < datalist.size(); ++ iZone)
    {
        delete datalist[iZone];
    }
}

void Post_ProbesWrite::Clear()
{
    this->datalist.clear();
}

}