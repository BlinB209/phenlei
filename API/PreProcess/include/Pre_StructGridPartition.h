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
//! @file      Pre_StructGridPartition.h
//! @brief     Structured grid partition: manager of block splitting.
//! @author    Baka, Guo Yongheng.

#include "Region.h"
#include "Pre_Block_Struct.h"

#pragma once

namespace PHSPACE
{
//! @brief Task_WriteBinaryTecplot class finish the structured grid partition.
//!
//! Pre_StrGridPartition
class Pre_StrGridPartition
{
public:
    LIB_EXPORT Pre_StrGridPartition(int numberOfProcessors);
    LIB_EXPORT ~Pre_StrGridPartition();

    //! Partition structured grid.
    LIB_EXPORT void Run();

private:
    //! Number of original blocks.
    int numberOfOriginalBlocks;

    //! Number of general blocks.
    int numberOfGeneralBlocks;

    //! Number of processors.
    int numberOfProcessors;

    //! Geometric dimension.
    int geometricDimension;

    int numberOfBasicBlocks;

private:
    void InitSplittingParameter();
    void BinarySplitBlocks();
    void GreedySplitBlocks();
    void SetBlockIndex();
    void OutputProcessorInformation();
    void OutputPartitionedGrid();

    void ReadBoundaryData();
    void ReadLnkFile(const string &fileName_in);
    void ReadLnkInfor(const string &fileName_in);
    void InitParititionGrids();
    void ProcessBCInfo();
    void WriteGridCoordinate(StructGrid *grid, Pre_Block_Struct *simpleBlock, vector<int> &nst, vector<int> &ned);

    void CheckMultigrid(vector<int> nst, vector<int> ned, int numberOfUnitCells);
    void CreateInterfaceBoundaryFace(Pre_Block_Struct *leftBlock, Pre_Block_Struct *rightBlock);
    void CreateInterfaceBoundaryFace(Pre_Block_Struct *leftBlock, Pre_Block_Struct *rightBlock, vector<int> &st, vector<int> &ed);
    bool GetInterfaceCoordinate(vector<int> &l_st, vector<int> &l_ed, vector<int> &r_st, vector<int> &r_ed, vector<int> &st, vector<int> &ed);
    bool IsSameDomain(vector<int> &l_st, vector<int> &l_ed, vector<int> &r_st, vector<int> &r_ed);
    void MarkPatchFace(vector<int> &st, vector<int> &ed, int geometricDimension);
    void OutputBoundaryInformation(Pre_Block_Struct *Pre_Block_Struct_in, StructGrid *grid, int iZone);
    void OutputBoundaryInformation(Pre_BcFace_Struct *leafBoundaryFace);
    void ExtractBoundaryInformation();
    void ExtractOriginalGridInformation();
    void OutputPartitionedGridFile();
    void TraceMark(Pre_Block_Struct *leftBlock, Pre_Block_Struct *rightBlock, Pre_Block_Struct *parentBlock);
    void DumpMarkInfo();

    void SetWallDist();
    void WriteGridWallDist(StructGrid *grid, StructGrid *Ordinary_grid, vector<int> &nst);
    void WriteWalldist(const string &out_grid_file, Grid **grids_in);
    void WriteWalldist(ActionKey *actkey, StructGrid *grid);

    void ProbeLargestBlock(Pre_Block_Struct *&simpleBlock, RDouble &numberOfCells);
    void SplitSingleBlock(Pre_Block_Struct *simpleBlock, RDouble ntarget, int ip, bool allotted);
    void AssignSingleBlock(Pre_Block_Struct *simpleBlock, int ip, RDouble numberOfCells);

    int GetDimension() { return this->geometricDimension; }
    int GetNumberOfProcessors() { return this->numberOfProcessors; }

    int GetMinimumElementIndex(vector<RDouble> &x);
    Grid ** GetOrdinaryGrid() { return this->OrdinaryGrid; }
    Grid ** GetParititionGrids() { return this->parititionGrids; }

private:
    RDouble delta, idealTime, epsilon;
    vector<RDouble> reserves;
    DataContainer *boundaryContainer;
    vector<Pre_Block_Struct *> referencedBlockContainer;
    set<Pre_Block_Struct *> unallottedBlockGroup, generalBlockGroup;
    vector<vector<Pre_Block_Struct *> > specialBlockTable;

    Grid **OrdinaryGrid;
    Grid **parititionGrids;

    bool LnkFileExist;
    bool partitionWallDist;

    int traceMark;
    int blockIndexOfMark;
    int *cellIndexOfMark;
};

int GetIntegerDigitWidth(int integerData);
}
