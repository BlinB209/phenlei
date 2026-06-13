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
//! @file      Post_WriteTecplot.h
//! @brief     Write flow field into tecplot file.
//! @author    Xu Gang.

#pragma once
#include "Post_WriteVisualFile.h"

namespace PHSPACE
{

//! @brief Post_WriteTecplot class realize tecpolt visual file output function.\n
class Post_WriteTecplot : public Post_WriteVisualFile
{
public:
    Post_WriteTecplot();
    ~Post_WriteTecplot();

private:
    void Initialize();

    //! Collect visual data for tecplot visual file.
    void StoreVisualizationData();
    void StoreVisualizationData(int zoneIndex, DataContainer *cdata, DataContainer *cBlkData);
    void StoreStrBoundaryVisualData(int zoneIndex, DataContainer *cdata);
    void StoreStrFieldVisualData(int zoneIndex, DataContainer *cdata);
    void StoreUnsBoundaryVisualData(int zoneIndex, DataContainer *cdata);
    void StoreUnsFieldVisualData(int zoneIndex, DataContainer *cdata);

    //! Dump data into visual file.
    void WriteFile();

    //! Clear data.
    void ClearFieldData();

private:
    vector<DataContainer *> visualDataList;

};

class Post_WriteTecplotByOriginalGrid 
{
private:
    int  nOriginalGrids;
    Grid **OriginalGrid;

public:
    Post_WriteTecplotByOriginalGrid();
    ~Post_WriteTecplotByOriginalGrid();

public:
    void Run();

private:
    void ServerWrite(ActionKey *actKey);
    void ComputerNodeDataOnPartitionGrid(int iZone, ActionKey *actKey);
    void CollectNodeData(ActionKey *actKey);
    void CollectStructedData(ActionKey *actKey, Post_Visual *postVisualization);
    void CollectUnstructedData(ActionKey *actKey, Post_Visual *postVisualization);

    void StoreVisualizationVar(ActionKey *actKey);
    void BoundaryVTKVisualization(Grid *gridIn, DataContainer *cData, RDouble4D **qn, int nVarPlot);
    void FieldVisualization(Grid *gridIn, DataContainer *cData, RDouble4D **qn, int nVarPlot);
    void BoundaryVisualization(Grid *gridIn, DataContainer *cData, RDouble4D **qn, int nVarPlot);
    void SaveDataForTecio(Grid *gridIn, DataContainer *cData, RDouble4D **qn, int nl);

    void BoundaryVTKVisualization(Grid *gridIn, DataContainer *cData, RDouble **qn, int nVarPlot);
    void FieldVisualization(Grid *gridIn, DataContainer *cData, RDouble **qn, int nVarPlot);
    void FieldVisualizationForVTK(Grid *gridIn, DataContainer *cData, RDouble **qn, int nVarPlot);
    void BoundaryVisualization(Grid *gridIn, DataContainer *cData, RDouble **qn, int nVarPlot);
    void SaveDataForTecio(Grid *gridIn, DataContainer *cData, RDouble **qn, int nVarPlot);

    void PostWrite(ActionKey *actKey);

    void UpdataNodeData(ActionKey *actKey, Grid *grid, RDouble4D **qn, int nVisualVariables);
    void UpdataNodeData(ActionKey *actKey, Grid *grid, RDouble **qn, int nVisualVariables);
    int GetOriginalGridProc(int iZone);

public:
    void Initialize();
    void SetOriginalGrid(Grid **gridsIn) { this->OriginalGrid = gridsIn; };
    void SetOriginalGridProc(int nOriginalGridsIn, int *originalGridProcOut);

private:
    int *originalGridProc;
    int *originalGridIndexOfEachGrid;
    vector <Post_Visual* > originalGridData;
    vector<DataContainer*> dataList;
    vector<DataContainer*> vtkDataList;

};

Post_WriteTecplotByOriginalGrid *GetWriteTecplotByOriginalGrid();


}
