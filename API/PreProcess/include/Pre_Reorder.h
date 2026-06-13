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
//! @file      Pre_Reorder.h
//! @brief     Unstructured grid reorder.
//! @author    Xi Zhang, Yong Zhang.

//! [1] Zhang Yong, Zhang Xi, Wan Yunbo, He Xianyao, Zhao Zhong, Lu Yutong.
//!     Optimizations of mesh renumbering for unstructured finite volume Computational Fluid Dynamics[C]. 
//!     Journal of Computer Research and Development,2021.10,p.123-131.
#include "Pre_GridBase.h"
using namespace std;
using namespace PHSPACE;
namespace Reorder
{
    //! Reorder face label for reducing GPU global memory accessing overheads
    //! In CGNSBase::ComputeGridConn after CallReorderCellLabel
    void CallReorderFaceLabel(Base_Grid_Conn *gconn);
    
    //! Reorder cell label for reducing bandwidth
    //! In CGNSBase::ComputeGridConn after PostProcess
    void CallReorderCellLabel(Base_Grid_Conn *gconn);
    
    //! Evaluate the bandwidth in CGNSBase::ComputeGridConn
    void CallBandwidthEvaluate(Base_Grid_Conn *gconn);
    
    //! Evaluate the bandwidth in GPUGeomInfoAllocCopy of BasicDeviceVariables
    void CallBandwidthFTS(const int nTotalCell, vector<int> *cell2Cell, const int *nNeighborCellOfEachCell);
    
    //! Get the current level of vertex rStart
    //! Input: rStart, nFaceOfEachCell, cell2cell, offsetCell2Cell, reorderCell
    //! Output: reorderCell, numSortCur, cellIDSortCur
    void GetCurrentLevel(int rStart, const int *nFaceOfEachCell, const int *cell2Cell, const int *offsetCell2Cell, int *reorderCell, int &numSortCur, int *cellIDSortCur);
    
    //! Get the next level of current level
    //! Input:
    //! Output:
    void GetNextLevel(const int &numSortCur, const int *cellIDSortCur, const int *nFaceOfEachCell, const int *cell2Cell, const int *offsetCell2Cell, int *reorderCell, int &numSortNext, int *cellIDSortNext);
    
    //! Copy cellIDSortNext into cellIDSortCur
    void SetCellIDSortCurByCellIDSortNext(const int &numSortNext, const int *cellIDSortNext, int &numSortCur, int *cellIDSortCur);
    
    //! Reset cellIDSortCur or cellIDSortNext by -1
    void ResetCellSort(int &numSort, int *cellIDSort);
    
    //! Output of cellIDSortCur or cellIDSortNext
    void TestCellSort(const int numSort, const int *cellIDSort);
    
    //! Find the smallest degree in cellIDSortCur
    int FindMinCellIDSortCur(const int numSortCur, const int *cellIDSortCur, const int *nFaceOfEachCell);
    
    //! Get current level of rStart for reOrder, the difference from  GetCurrentLevel    is that cancel the set of reorderCell by 1.
    void GetCurrentLevelForReOrder(int rStart, const int *nFaceOfEachCell, const int *cell2Cell, const int *offsetCell2Cell, int *reorderCell, int &numSortCur, int *cellIDSortCur);
    
    //! Arrange cellIDSortCur by Ascending order
    void AscendingCellIDSortCurByDegree(const int numSortCur, const int *nFaceOfEachCell, int *cellIDSortCur);
    
    //! Test degree of cellIDSortCur from 0 to numSortCur
    void TestCellIDSortCurDegree(const int numSortCur, const int *cellIDSortCur, const int *nFaceOfEachCell);
    
    //! Check cell2cell by 3 rules:
    //! rule1 one cell's neighbor cells should not contain itself.
    //! rule2 0<= one cell's neighbor cells < nTotalCell
    //! rule3 one cell's neighbor cells should be different.
    void CheckCell2Cell(const int nTotalCell, const int *nFaceOfEachCell, const int *offsetCell2Cell, const int * cell2Cell);
    
    //! Reorder cells in cellIDSortCur
    void ReOrderLevelCur(const int numSortCur, const int *cellIDSortCur, int &numVertex, int *reorderCell);

    void GetNextLevelForReOrder(const int &numSortCur, const int *cellIDSortCur, const int *nFaceOfEachCell, const int *cell2Cell, const int *offsetCell2Cell, int *reorderCell, int &numSortNext, int *cellIDSortNext);
    
    //! Check rules for cell index
    //! rule1 0 <= newCell < nTotalCell
    //! rule2 cell index should be different
    void CheckNewCell(const int nTotalCell, const int *newCell);
    
    //! Compute nFaceOfEachCell and nTotalCellFaces
    void SetFaceNumberOfEachCell(const int nTotalFace, const int nTotalCell, const int *leftCellOfFace, const int *rightCellOfFace, int *nFaceOfEachCell, int &nTotalCellFaces);
    
    //! Set offsetCell2Face by nFaceOfEachCell
    void SetOffSetCell2Face(const int nTotalCell, const int *nFaceOfEachCell, int *offsetCell2Face);
    
    //! Set cell2face by leftCellOfFace, rightCellOfFace, offsetCell2Face
    //! It should be noted that boundary faces are not considered here.
    void SetCell2Face(const int nTotalFace, const int nTotalCell, const int *leftCellOfFace, const int *rightCellOfFace, int *offsetCell2Face, int *cell2Face);
    
    //! Set cell2Cell
    void SetCell2Cell(const int nTotalCell, const int *nFaceOfEachCell, const int *offsetCell2Face, const int *cell2Face, const int *leftCellOfFace, const int *rightCellOfFace, int *cell2Cell);
    
    //! Compute bandwidth of each cell
    void ComputeBandwidth(const int nTotalCell, const int *nFaceOfEachCell, const int *offsetCell2Face, const int *cell2Cell);
    
    //! Update face_number_of_each_face by swap method
    //! Intput: nTotalCell, face_number_of_each_face, reflectCell
    //! Output:nSwapFaceOfEachCell
    //! The function is used to provide a swap method to update a cell related varible
    //! After cell reorder. In the future, cell2node can be update by the method
    void SwapFaceNumberOfEachFace(const int nTotalCell, const int *nFaceOfEachCell, const int *reflectCell, const int *nNewFaceOfEachCell);
    
    //! Check swap method for face_number_of_each_face
    //! After swap, OrgCell=reflectCell, SwapCell=newCell, nSwapFaceOfEachCell=nNewFaceOfEachCell
    void CheckSwapFaceNumber(const int nTotalCell, const int *orgCell, const int *reflectCell, const int *nSwapFaceOfEachCell, const int *nNewFaceOfEachCell);
    
    //! Reorder cell label from nBoundFace to nTotalFace just with the order of faces in cell2face
    void CellReorderByCell2FaceOrder(const int nTotalCell, const int nTotalFace, const int *offsetCell2Face, const int *nFaceOfEachCell, const int *cell2Face, const int nBoundFace, int labelFace, int *newFaces);
    
    //! A caller for MetricLeftRightCellOfFace
    void CallMetricLeftRightCellOfFace(Base_Grid_Conn *gconn);
    
    //! Metric for computing leftCellOfFace and rightCellOfFace for GPU 
    void MetricLeftRightCellOfFace(const int nTotalFace, const int nBoundFace, const int *leftCellOfFace, const int *rightCellOfFace);
}
