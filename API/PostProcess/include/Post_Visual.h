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
//! @file      Post_Visual.h
//! @brief     Post_Visual defines the post visualization variables.
//! @author    Bell.

#pragma once
#include "GridType.h"
#include "Geo_UnstructGrid.h"
#include "Geo_StructGrid.h"

using namespace std;

namespace PHSPACE
{

class Grid;
//! @brief Post_Visual class is used to manage post visualization variables data.
class Post_Visual
{
public:
    //! @param[in] nVisualVariables_in    Number of visual variables, parsed from parameter file usually.
    //! @param[in] visualVariables_in     Type of visual variables, parsed from parameter file usually.
    LIB_EXPORT Post_Visual(int nVisualVariables_in, int *visualVariables_in, Grid *grid_in = NULL, int flowType_in = 0);

    LIB_EXPORT ~Post_Visual();

public:
    //! Return the ASCII type of tecplot title header.
    LIB_EXPORT void GetASCIITecplotHeader(vector<string> &title_tecplot);

    //! Return the BINARY type of tecplot title header.
    LIB_EXPORT void GetBinaryTecplotHeader(vector<string> &title_tecplot);

    //! Return number of visual variables.
    int GetNumberofVisualVariables();

    //! Return if the variableType should be visualized.
    bool IsNeedVisualization(int variableType);

    //! Compute variables on node through variables on cell center and or update to database.
    //! Unstructured grid data type.
    //! @param[in] name              Source data's name.
    //! @param[in] cellCenterData    Pointer the Source cell center data.
    LIB_EXPORT void UpdateVisualNodeVarPtr(const string &name, RDouble *cellCenterData);
    LIB_EXPORT void UpdateVisualNodeVarPtr(const string &name, int *cellCenterData);

    //! Compute variables on node through variables on cell center and or update to database.
    //! Structured grid data type.
    //! @param[in] name              Source data's name.
    //! @param[in] cellCenterData    Pointer the Source cell center data.
    //! @param[in] index             The array index from RDouble4D array,
    //!                              i.e we only need the cellCenterData[...][...][...][index].
    LIB_EXPORT void UpdateVisualNodeVarPtr(const string &name, RDouble4D &cellCenterData, int index);

    //! Compute variables on node through variables on cell center and or update to database.
    //! Structured grid data type.
    //! @param[in] name              Source data's name.
    //! @param[in] cellCenterData    Pointer the Source cell center data.
    LIB_EXPORT void UpdateVisualNodeVarPtr(const string &name, RDouble3D &cellCenterData);

    //! Compute overset information on node through variables on cell center and or update to database.
    //! Structured grid data type.
    //! @param[in] name              Source data's name.
    //! @param[in] cellCenterData    Pointer the Source cell center data.
    LIB_EXPORT void UpdateVisualNodeVarPtr(const string &name, Int3D &oversetcellCenterdata);

    LIB_EXPORT void UpdateVisualNodeNotInterpolation(const string &name, int *cellCenterdata);

    //! Obtain data's pointer(address) from database by name.
    //! If the dataSet doesn't has a same name data, returns null.
    //! @param[in] name    Source data's name.
    void * GetVisualNodeVarPtr(const string &name);

    //! Get all visualization varialbes' pointer (RDouble * for unstruct grid).
    LIB_EXPORT void GetAllVisualNodeVarPtr(RDouble **qn);

    //! Get all visualization varialbes' pointer (RDouble4D * for struct grid).
    LIB_EXPORT void GetAllVisualNodeVarPtr(RDouble4D **qn);

    //! Get variables' name through variable's tag index.
    //! @param[in] variableNumber    variable's index.
    string GetVariableName(int variableNumber) const;

    //! Get set of variables.
    const set<int> & GetVisualVariables() const;

    Grid * GetGrid() const;
    int GetFlowType() { return this->flowType; }
    void UpdateVisualNodeVarPtr(const string &name, void *data);

    void UpdateDataPtr(const string &name, void *data);

private:
    //! Geometry information.
    Grid *grid;

    //! Grid type.
    int gridType;

    //! Number of visual variables.
    int nVisualVariables;

    //! Type of visual variables, each number represents ONE variable.
    set<int> visualVariables;

    //! Variable names, displayed in Tecplot.
    map<int, string> visualVariablesMap;

    //! Store visual variables.
    Data_Field *visualVariablesPtr;

    int flowType;

private:
    //! Free all allocated memory by visual node variables' pointer.
    void FreeAllVisualNodeVarPtr();

    //! Free visual variables memory.
    //! @param[in] name    variables' name.
    void DeleteVisualNodeVarPtr(const string &name);
};

const int VISUAL_DENSITY               = 0;
const int VISUAL_U                     = 1;
const int VISUAL_V                     = 2;
const int VISUAL_W                     = 3;
const int VISUAL_PRESSURE              = 4;
const int VISUAL_TEMPERATURE           = 5;
const int VISUAL_MACH                  = 6;
const int VISUAL_VISCOSITY_LAMINAR     = 7;
const int VISUAL_VISCOSITY_TURBULENT   = 8;
const int VISUAL_VORTICITY_X           = 9;
const int VISUAL_VORTICITY_Y           = 10;
const int VISUAL_VORTICITY_Z           = 11;
const int VISUAL_VORTICITY_MAGNITUDE   = 12;
const int VISUAL_STRAIN_RATE           = 13;
const int VISUAL_Q_CRITERIA            = 14;
const int VISUAL_CP                    = 15;
const int VISUAL_TIME_STEP             = 16;
const int VISUAL_VOLUME                = 17;
const int VISUAL_MODELED_TKE           = 18;
const int VISUAL_MODELED_DISSIPATION   = 19;
const int VISUAL_SST_F1                = 20;
const int VISUAL_SST_F2                = 21;
const int VISUAL_INCOMP_U              = 22;
const int VISUAL_INCOMP_V              = 23;
const int VISUAL_INCOMP_W              = 24;
const int VISUAL_INCOMP_P              = 25;
const int VISUAL_INCOMP_CP             = 26;
const int VISUAL_INCOMP_T              = 27;
const int VISUAL_INCOMP_ENTHALPY       = 70;
const int VISUAL_INCOMP_DEN            = 28;
const int VISUAL_INCOMP_VIS            = 29;
const int VISUAL_KN                    = 30;  //! Useless.
const int VISUAL_INCOMP_TE             = 31;
const int VISUAL_INCOMP_ED             = 32;
const int VISUAL_INCOMP_WD             = 63;
const int VISUAL_INCOMP_SPECIES        = 71;
const int VISUAL_TEMPERATURE_VIBRATION = 33;
const int VISUAL_TEMPERATURE_ELECTRON  = 34;
const int VISUAL_ENERGY_VIBRATION      = 35;
const int VISUAL_ENERGY_ELECTRON       = 36;
const int VISUAL_ELECTRON_NUMBER       = 37;
const int VISUAL_DIMENSIONAL_DENSITY   = 38;
const int VISUAL_DIMENSIONAL_PRESSURE  = 39;
const int VISUAL_DIMENSIONAL_TEMPERATURE = 40;
const int VISUAL_GRADIENT_UX           = 41;
const int VISUAL_GRADIENT_UY           = 42;
const int VISUAL_GRADIENT_VX           = 43;
const int VISUAL_GRADIENT_VY           = 44;
const int VISUAL_STREAMLINE_U          = 45;
const int VISUAL_STREAMLINE_V          = 46;
const int VISUAL_STREAMLINE_W          = 47;
const int VISUAL_TRANSITION_GAMAEFF    = 48;
const int VISUAL_TRANSITION_RESCF      = 49;
const int VISUAL_WALL_DIST             = 50;
const int VISUAL_INTERMITTENCY         = 51;
const int VISUAL_MOMENTUMTHICKREYNOLDS = 52;
const int VISUAL_DIMENSIONAL_U         = 53;
const int VISUAL_DIMENSIONAL_V         = 54;
const int VISUAL_DIMENSIONAL_W         = 55;
const int VISUAL_GAMA                  = 56;
const int VISUAL_CFL1                  = 57;
const int VISUAL_CFL2                  = 58;
const int VISUAL_STREAMLINE_MACH       = 59;
const int VISUAL_KNUDSEN_NUMBER        = 60;  //! Useless.
const int VISUAL_DAMKOHLER_NUMBER      = 61;
const int VISUAL_VIBNONEQ_NUMBER       = 62;
const int VISUAL_VELOCITY_MAGNITUDE    = 63;

const int VISUAL_IBLANK                = 81;
//const int VISUAL_O2                    = 15;
//const int VISUAL_NO                    = 16;
//const int VISUAL_N                     = 17;
//const int VISUAL_N2                    = 18;

const int VISUAL_TAU_XX = 90;
const int VISUAL_TAU_YY = 91;
const int VISUAL_TAU_ZZ = 92;
const int VISUAL_TAU_XY = 93;
const int VISUAL_TAU_XZ = 94;
const int VISUAL_TAU_YZ = 95;

//Visual Wall Variables.
const int VISUAL_SPECIES_START = 30;
const int VISUAL_WALL_CP       = 0;
const int VISUAL_WALL_CF       = 1;
const int VISUAL_WALL_YPLUS    = 2;
const int VISUAL_WALL_QNONDIM  = 3;
const int VISUAL_WALL_QDIM     = 4;
const int VISUAL_WALL_PW       = 5;
const int VISUAL_WALL_TW       = 6;
const int VISUAL_WALL_RHOW     = 7;
const int VISUAL_WALL_QT       = 8;
const int VISUAL_WALL_QS       = 9;
const int VISUAL_WALL_QV       = 10;
const int VISUAL_WALL_QE       = 11;
const int VISUAL_WALL_NS       = 12;
const int VISUAL_WALL_VX       = 13;
const int VISUAL_WALL_VY       = 14;
const int VISUAL_WALL_VZ       = 15;
const int VISUAL_SLIP_TS       = 16;
const int VISUAL_SLIP_TV       = 17;
const int VISUAL_SLIP_TE       = 18;
const int VISUAL_WALL_VS       = 19;
const int VISUAL_WALL_ST       = 20;
const int VISUAL_WALL_CH       = 21;
const int VISUAL_SLIP_DTS      = 22;
const int VISUAL_WALL_RE       = 23;
const int VISUAL_WALL_KN       = 24; //The Knudsen number.
//! The code number must be less than VISUAL_SPECIES_START.

//! A controller to select BLK data or surface data as output data generating tecplot binary files
//! in the condition of THREE_D mode and FieldVisual type.
const int OUTPUT_SURFACE       = 0;
const int OUTPUT_BLK           = 1;

class Post_VisualWall
{
public:
    LIB_EXPORT Post_VisualWall();

    LIB_EXPORT Post_VisualWall(int nVisualNumber, int *visualVariablesType);

    LIB_EXPORT ~Post_VisualWall();

public:
    int GetVisualVariablesNumber() const { return nVisualVariables; }

    //! Get variables' name through variable's tag index.
    //! @param[in] variableNumber    variable's index.
    string GetVariableName(int variableNumber) const;

    int * GetVisualVariablesType() const { return visualVariables; }

    int GetVisualVariablesType(int nIndex) const { return visualVariables[nIndex]; }

private:
    //! Number of visual variables.
    int nVisualVariables;

    //! Type of visual variables, each number represents ONE variable.
    int *visualVariables;

    //! Variable names, displayed in Tecplot.
    map<int, string> visualVariablesMap;
};

const int AverageFlow = 1;
const int AverageReynoldsStress = 2;

#include "Post_Visual.hxx"

class DumpToVTK
{
private:
    vector<DataContainer *> vtkdatalist;
    string vtkfilename;

    //! Visual flow type: 0 -- main flow;  1 -- average flow.
    int visualFlowType;

public:
    DumpToVTK(vector<DataContainer *> datalist_in, string filename_in, bool CharacteristicBoundary_in, int visualFlowType = 0);
    ~DumpToVTK();

public:
    void Run();

private:
    void InitVariables();

    //! Set number of variables in tecplot.
    void SetNvarplot(int nvarplot_in) { this->nvarplot = nvarplot_in; }

    //! Return number of variables in tecplot.
    int GetNvarplot() { return nvarplot; }

    //! Get the file name in VTK format.
    string GetFilename() { return vtkfilename; }

    void WriteVTKFilm();

private:
    int nvarplot;
    bool CharacteristicBoundary;
};

#ifdef USE_TecplotLib
class DumpToTecio
{
private:
    //ActionKey *actkey;
    vector<DataContainer *> datalist;
    string filename;

    //! Visual flow type: 0 -- main flow;  1 -- average flow.
    int visualFlowType;

public:
    DumpToTecio(vector<DataContainer *> datalist_in, string filename_in, int visualFlowType = 0);
    ~DumpToTecio();

public:
    void Run();

private:
    void WriteGridFilm_old();
    void WriteGridFilm(int isOutputBLK = OUTPUT_SURFACE);
    void WriteSliceFilm();
    void InitVariables();

    //! Set the variables.
    void SetVariables(string Variables_in) { this->Variables = Variables_in; }

    //! Get the variables.
    string GetVariables() { return Variables; }

    //! Set number of variables in tecplot.
    void SetNvarplot(int nvarplot_in) { this->nvarplot = nvarplot_in; }

    //! Return number of variables in tecplot.
    int GetNvarplot() { return nvarplot; }

    //! Get the file name.
    string GetFilename() { return filename; }
private:
    string Variables;
    int nvarplot;
};
#endif

class DumpToTecplotASCII
{
private:
    vector<DataContainer *> dataList;
    string fileName;

    //! Visual flow type: 0 -- main flow;  1 -- average flow.
    int visualFlowType;

public:
    DumpToTecplotASCII(vector<DataContainer *> datalist_in, string filename_in, int visualFlowType = 0);
    ~DumpToTecplotASCII();

public:
    void Run();

    //! Set the variables.
    void SetVariables(string Variables_in) { this->Variables = Variables_in; }

    //! Get the variables.
    string GetVariables() { return Variables; }

    //! Set number of variables in tecplot.
    void SetNvarplot(int nVariable_in) { this->nVariable = nVariable_in; }

    //! Return number of variables in tecplot.
    int GetNvarplot() { return nVariable; }

    //! Get the file name.
    string GetFilename() { return fileName; }

private:
    void InitVariables();
    void WriteFlowFilm();
    void DumpOrderedFlow(fstream &flowFile, DataContainer *data, string zoneTitle);
    void DumpFlowByCellTopo(fstream &flowFile, DataContainer *data, string zoneTitle);
    void DumpFlowByFaceTopo(fstream &flowFile, DataContainer *data, string zoneTitle);

private:
    int nVariable;
    string Variables;
};

class PHCutPlane
{
public:
    char **title;
    int nTitle;
    int neqn;
    int nTotalNode, nTotalFace;
    RDouble **data;
    RDouble *x, *y, *z;
    int *faceIndex[4];
public:
    PHCutPlane();
    ~PHCutPlane();
    void SortData(RDouble *value);
public:
};

template < typename T >
void ReorderList(T *x, vector<int> &index_list)
{
    T *tmp = new T [index_list.size()];

    for (std::size_t i = 0; i < index_list.size(); ++ i)
    {
        tmp[i] = x[i];
    }

    for (std::size_t i = 0; i < index_list.size(); ++ i)
    {
        x[i] = tmp[index_list[i]];
    }

    delete [] tmp;
}

void CutPlane(UnstructGrid *grid, RDouble **node_var, int nvar, RDouble cut, PHCutPlane *cutpl, int ndir);
void CutPlaneCellCenter(UnstructGrid *grid, RDouble **q_var, int nvar, RDouble cut, PHCutPlane *cutpl, int ndir);

void write_string(fstream &file, string &cs);
void read_string(fstream &file, string &cs);
}