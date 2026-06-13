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
//! @file      ActionKey.h
//! @brief     Explain this file briefly.
//! @author    He Xin, Bell.

#pragma once
#include "DataContainer.h"
#include "cstdint"
#pragma warning( disable : 26495 )
using namespace std;

namespace PHSPACE
{
typedef enum
{
    READ_RESTART = 1,
    READ_INTERPOLATE,
    BUILD_WALL_STRUCT,
    CALC_WALL_DIST,
    FILL_WALL_STRUCT,
    READ_WALL_DIST,
    WRITE_WALL_DIST,
    DUMP_RESTART,
    UPDATE_INTERFACE_DATA,
    UPLOAD_INTERFACE_DATA,
    DOWNLOAD_INTERFACE_DATA,
    UPDATE_INTERPOINT_DATA,
    UPLOAD_INTERPOINT_DATA,
    DOWNLOAD_INTERPOINT_DATA,
    UPDATE_OVERSET_DATA,
    UPLOAD_OVERSET_DATA,
    DOWNLOAD_OVERSET_DATA,
    DUMP_COMPONENT_AIR_FORCE_COEF,
    DUMP_AIR_FORCE_COEF,
    DUMP_CP_DISTRI,
    VISUALIZATION,
    DUMP_RESIDUAL,
    ALLOCATEINTERFACE,
    COMPUTEMETRICS,
    COMMCELLCENTERDATA,
    COMMCELLIBLANK,
    ALLOCATE_WALL_DIST,
    TEST_RECONSTRUCTION,
    TEST_INTERFACE,
    SIMPLE_ACTION,
    SHOW_VISTMAX,
    INIT_FIRST,
    INIT_SECOND,
    INIT_FINAL,
    DUMP_AVERAGE_FLOW,
    READ_AVERAGE_FLOW,
    VISUALIZATION_AVERAGE_FLOW,
    DUMP_SURFACE_INFO,
    DUMP_PROBES_FLOW,
    VISUALIZATION_AVERAGE_ReynoldsStress,
    DUMP_HEATFLUX,
    COMPUTE_AERODYNAMIC_FORCE
}ActionEnum;
const int NofValidAction = 41;
extern const char * ActionName[NofValidAction];

class ActionKey
{
public:
    int action;
    int subact;
    int solver;
    int solverID;
    int level;
    int ipos;
    int kind;    //! 0 solver, 1 grid.
    int64_t filepos;    //! new version of HDF5 library requires int64_t type of hid_t, instead of int type.
    int gridfilepos;
    int format;
    fstream *file;
    bool del_file;
    bool VTKvisual;
    ios_base::openmode openmode;
    string filename;
    string gridfilename;
    string tecfilename;
    DataContainer *data;
    DataContainer *tecdata;
    string taskname;
public:
    ActionKey();
    ~ActionKey();
    DataContainer * GetData() const { return data; }
    DataContainer * GetTecData() const { return tecdata; }
    fstream * GetFile() const { return file; }
    string GetFileName() const { return filename; }
    bool GetVTKVisual() const { return VTKvisual; }
    int GetAction() { return action; }
    void SetData(DataContainer *dataIn) { this->data = dataIn; }
    void SetTecData(DataContainer *dataIn) { this->tecdata = dataIn; }
    bool IsNeedOpenFile() const { return this->filename != ""; }
};

class ActionTag
{
public:
    ActionTag(){};
    ActionTag(int action, int level, bool run_now = true)
    {
        this->action  = action;
        this->level   = level;
        this->run_now = run_now;
    };
    ~ActionTag(){};
public:
    int  action;
    int  level;
    bool run_now;
};

void SetCurrentActionKey(ActionKey *ActionKeyIn);
ActionKey * GetCurrentActionKey();

}

