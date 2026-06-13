#include "ActionKey.h"
using namespace std;

namespace PHSPACE
{
const char * ActionName[NofValidAction] =
{
    "READ_RESTART",
    "BUILD_WALL_STRUCT",
    "CALC_WALL_DIST",
    "FILL_WALL_STRUCT",
    "READ_WALL_DIST",
    "WRITE_WALL_DIST",
    "DUMP_RESTART",
    "UPDATE_INTERFACE_DATA",
    "UPLOAD_INTERFACE_DATA",
    "DOWNLOAD_INTERFACE_DATA",
    "UPDATE_INTERPOINT_DATA",
    "UPLOAD_INTERPOINT_DATA",
    "DOWNLOAD_INTERPOINT_DATA",
    "UPDATE_OVERSET_DATA",
    "UPLOAD_OVERSET_DATA",
    "DOWNLOAD_OVERSET_DATA",
    "DUMP_COMPONENT_AIR_FORCE_COEF",
    "DUMP_AIR_FORCE_COEF",
    "DUMP_CP_DISTRI",
    "VISUALIZATION",
    "DUMP_RESIDUAL",
    "ALLOCATEINTERFACE",
    "COMPUTEMETRICS",
    "COMMCELLCENTERDATA",
    "COMMCELLIBLANK",
    "ALLOCATE_WALL_DIST",
    "TEST_RECONSTRUCTION",
    "TEST_INTERFACE",
    "SIMPLE_ACTION",
    "SHOW_VISTMAX",
    "INIT_FIRST",
    "INIT_SECOND",
    "INIT_FINAL",
    "DUMP_AVERAGE_FLOW",
    "READ_AVERAGE_FLOW",
    "VISUALIZATION_AVERAGE_FLOW",
    "DUMP_SURFACE_INFO",
    "DUMP_PROBES_FLOW",
    "DUMP_HEATFLUX",
    "COMPUTE_AERODYNAMIC_FORCE"
};

ActionKey::ActionKey()
{
    data    = new DataContainer();
    tecdata = new DataContainer();

    action      = 0;
    subact      = 0;
    solver      = 0;
    solverID    = 0;
    level       = 0;
    ipos        = 0;
    kind        = 0;
    filepos     = 0;
    gridfilepos = 0;
    format      = 0;
    file        = 0;
    del_file    = false;
    VTKvisual   = false;
    openmode    = ios_base::out;
    filename  = "";
    gridfilename  = "";
    tecfilename = "";
    taskname  = "";
}

ActionKey::~ActionKey()
{
    delete data;    data = NULL;
    delete tecdata;    tecdata = NULL;
}

ActionKey *currentActionKey = 0;

void SetCurrentActionKey(ActionKey *currentActionKeyIn)
{
    currentActionKey = currentActionKeyIn;
}

ActionKey * GetCurrentActionKey()
{
    return currentActionKey;
}

}
