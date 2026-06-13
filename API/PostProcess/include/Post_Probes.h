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
//! @file      Post_Probes.h
//! @brief     Post_Probes defines the post variables of probes that need to be monitored.
//! @author    Meng Liyuan.

#pragma once
#include "GridType.h"
#include "Geo_Grid.h"
#include "Geo_UnstructGrid.h"
#include "Geo_StructGrid.h"

using namespace std;

namespace PHSPACE
{

class Grid;
class SampleLocationInfo;

//! @brief Post_Probes class is used to manage post variables data of probes.
class Post_Probes
{
public:
    //! @param[in] nVariables_in    Number of probe variables, parsed from parameter file usually.
    //! @param[in] Variables_in     Type of probe variables, parsed from parameter file usually.
    LIB_EXPORT Post_Probes(int nProbeVariables_in, int *probeVariables_in, Grid *grid_in = NULL);

    LIB_EXPORT ~Post_Probes();

public:
    //! Return the ASCII type file title header.
    LIB_EXPORT void GetASCIIFileHeader(vector<string> &title_tecplot);

    //! Return the probes ASCII type file title header.
    LIB_EXPORT void GetProbesASCIIFileHeader(vector<string> &title_tecplot);

    //! Return the lines or surfaces ASCII type file title header.
    LIB_EXPORT void GetLinesOrSurfacesASCIIFileHeader(vector<string> &title_tecplot);

    //! Return number of probe variables.
    int GetNumberofProbeVariables();

    //! Return if the variableType should be Monitored.
    bool IsNeedMonitoring(int variableType);

    //! Compute variables data of probe through variables on probes cell center.
    //! Unstructured grid data type.
    //! @param[in] name              Variable data's name.
    //! @param[in] cellCenterData    Pointer the cell center data.
    LIB_EXPORT void UpdateProbesVarPtr(const string &name, RDouble *cellCenterData);

    //! Compute variables data of probe by using the value of probes cell as probe real value.
    LIB_EXPORT void UpdateProbesVarPtrM1(const string &name, RDouble *cellCenterData);

    //! Interpolation from probes and neighboring cell to probe.
    LIB_EXPORT void UpdateProbesVarPtrM2(const string &name, RDouble *cellCenterData);

    //! Interpolation from probes cell nodes to probe.
    LIB_EXPORT void UpdateProbesVarPtrM3(const string &name, RDouble *cellCenterData);

    //! Compute the necessary variables for interpolation.
    LIB_EXPORT void NecessaryVarForInterpolation(UnstructGrid *grid, int iprobe, int index, vector<RDouble> &numerator, vector<int> &interpolationIndexs, bool flag);

    //! Compute variables data of probe through variables on probes cell center.
    //! Structured grid data type.
    //! @param[in] name              Variable data's name.
    //! @param[in] cellCenterData    Pointer the nearest cell center data.
    //! @param[in] index             The array index from RDouble4D array,
    //!                              i.e we only need the cellCenterData[...][...][...][index].
    LIB_EXPORT void UpdateProbesVarPtr(const string &name, RDouble4D &cellCenterData, int index);

    //! Compute variables data of probe by using the value of probes cell as probe real value.
    LIB_EXPORT void UpdateProbesVarPtrM1(const string &name, RDouble4D &cellCenterData, int index);

    //! Interpolation from probes and neighboring cell to probe.
    LIB_EXPORT void UpdateProbesVarPtrM2(const string &name, RDouble4D &cellCenterData, int index);

    //! Interpolation from probes cell nodes to probe.
    LIB_EXPORT void UpdateProbesVarPtrM3(const string &name, RDouble4D &cellCenterData, int index);

    //! Compute the necessary variables for interpolation.
    LIB_EXPORT void NecessaryVarForInterpolation(StructGrid *grid, int iprobe, int i, int j, int k, vector<RDouble> &numerator, vector<vector<int> > &interpolationIndexs, bool flag);

    //! Obtain probe variable data's pointer(address) from database by name.
    //! If the dataSet doesn't has a same name data, returns null.
    //! @param[in] name    Variable data's name.
    void * GetProbeVarPtr(const string &name);

    //! Get all probe variables' pointer.
    LIB_EXPORT void GetAllProbesVarPtr(RDouble **q_probes);

    //! Get probe variables' name through variable's index.
    //! @param[in] variableNumber    variable's index.
    string GetProbesVariableName(int variableNumber) const;

    //! Get set of variables.
    const set<int> & GetProbeVariables() const;

    map<int, string> probesVariablesMap;

    Grid * GetGrid() const;

private:
    //! Geometry information.
    Grid *grid;

    //! Grid type;
    int gridType;

    //! Number of probe variables.
    int nProbeVariables;

    //! Type of probe variables, each number represents ONE variable.
    set<int> probeVariables;

    //! Store probe variables.
    Data_Field *probeVariablesPtr;

private:
    //! Free all allocated memory by probe variables' pointer.
    void FreeAllProbeVarPtr();

    //! Free probe variables memory.
    //! @param[in] name    variables' name.
    void DeleteProbeVarPtr(const string &name);
};

const int PROBE_DENSITY     = 0;
const int PROBE_U           = 1;
const int PROBE_V           = 2;
const int PROBE_W           = 3;
const int PROBE_PRESSURE    = 4;
const int PROBE_TEMPERATURE = 5;
const int PROBE_MACH        = 6;
const int PROBE_DIMENSIONAL_DENSITY     = 7;
const int PROBE_DIMENSIONAL_U           = 8;
const int PROBE_DIMENSIONAL_V           = 9;
const int PROBE_DIMENSIONAL_W           = 10;
const int PROBE_DIMENSIONAL_PRESSURE    = 11;
const int PROBE_DIMENSIONAL_TEMPERATURE = 12;
const int PROBE_VELOCITY_MAGNITUDE      = 13;
const int PROBE_VISCOSITY_LAMINAR     = 14;
const int PROBE_VISCOSITY_TURBULENT   = 15;
const int PROBE_MODELED_TKE           = 16;
const int PROBE_MODELED_DISSIPATION   = 17;

#include "Post_Probes.hxx"

class DumpToFileASCII
{
private:
    vector <DataContainer *> dataList;
    string fileName;

public:
    DumpToFileASCII(vector<DataContainer *> datalist_in, string filename_in);
    ~DumpToFileASCII();

public:
    void Run();

    //! Set the variables.
    void SetVariables(string Variables_in) { this->variables = Variables_in; }

    //! Get the variables.
    string GetVariables() { return variables; }

    //! Set number of variables in file.
    void SetnVariablesToDump(int nVariables_in) { this->nVariables = nVariables_in; }

    //! Return number of variables in file.
    int GetnVariablesToDump() { return nVariables; }

    //! Get the file name.
    string GetFilename() { return fileName; }

private:
    //! Init the probes variables need to dump.
    void InitProbesVariables();

    void WriteProbeFlowFile();

    //! Write the variables data into file.
    void WriteFlowFile();

    //! Write the probes variables data into file.
    void WriteProbesFlowFile();
    //! Write the probes variables data into file.
    void WriteProbesFlowFile(int outerStep);
    //! Write the probes of lines variables data into file.
    void WriteLinesFlowFile();
    //! Write the probes of lines variables data into file.
    void WriteLinesFlowFile(int outerStep);
    //! Write the probes of surfaces variables data into file.
    void WriteSurfacesFlowFile();
    //! Write the probes of surfaces variables data into file.
    void WriteSurfacesFlowFile(int outerStep);

private:
    int nVariables;
    string variables;
};

}