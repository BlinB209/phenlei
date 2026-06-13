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
//! @file      Post_WriteVisualFile.h
//! @brief     Write flow field into post visual file.
//! @author    Xu Gang.

#pragma once
#include "Geo_Grid.h"
#include "Post_Visual.h"

namespace PHSPACE
{

//! Dump results into visual file.
//! @param[in]  flowTypeIn         Flow field data type.
//!             flowTypeIn = 0,    Default.
//!             flowTypeIn = 1,    AverageFlow
//!             flowTypeIn = 2,    AverageReynoldsStress
void WriteVisualFile(int flowTypeIn = 0);


//! @brief Post_WriteVisualFile class realize visual file output function.\n
class Post_WriteVisualFile
{
public:
    Post_WriteVisualFile();
    ~Post_WriteVisualFile();

public:
    //! @param[in]  flowTypeIn         Flow field data type.
    //!             flowTypeIn = 0,    Default.
    //!             flowTypeIn = 1,    AverageFlow
    //!             flowTypeIn = 2,    AverageReynoldsStress
    void Run(int flowTypeIn = 0);

    //! Boundary face grid will be reconstruct next output, used when mesh changed, such as overset grid.
    void BoundaryGridNeedConstruct();

public:
    //! Store visual variables on node after interpolation.
    vector <Post_Visual * > flowFieldDataOnNode;

    //! Boundary face grid.
    vector <Grid * > boundaryGrid;

    //! The grid index correspondence between boundary face grid and original grid.
    vector <int> gridMap;

    //! The node index correspondence between boundary face grid and original grid.
    vector <int **> nodeMap;

    //! Boundary face grid name.
    vector <string> boundaryName;

    //! Boundary face grid type.
    vector <int> boundaryType;

private:
    //! Construct boundary face grid.
    void ConstructBoundaryGrid();
    void ConstructStrBoundaryGrid(Grid *gridIn);
    void ConstructUnsBoundaryGrid(Grid *gridIn);

    //! Compute node data.
    void ComputeNodeData();

    //! Collect visual data for software interface.
    void StoreVTKData();
    void StoreVTKData(int zoneIndex, DataContainer *cdata);
    void StoreStrBoundaryVTKData(int zoneIndex, DataContainer *cdata);
    void StoreStrFieldVTKData(int zoneIndex, DataContainer *cdata);
    void StoreUnsBoundaryVTKData(int zoneIndex, DataContainer *cdata);
    void StoreUnsFieldVTKData(int zoneIndex, DataContainer *cdata);
    void StoreUnsFieldVTKData2D(int zoneIndex, DataContainer *cData);
    void StoreUnsFieldVTKData3D(int zoneIndex, DataContainer *cdata);

private:
    virtual void Initialize() = 0;

    //! Collect visual data for tecplot visual file.
    virtual void StoreVisualizationData() {};

    //! Dump data into visual file.
    virtual void WriteFile() = 0;

    //! Clear data.
    virtual void ClearFieldData() = 0;

public:
    int    flowType;
    bool   boundaryGridNeedConstruct;
    bool   gridChanged;
    bool   VTKvisual;
    string VTKFileName;
    string visualFileName;
    vector<DataContainer *> VTKDataList;
};

bool WantVisualField(Grid *grid);


}