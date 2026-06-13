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
//! @file      GPUComputeNodeValue.h
//! @brief     Compute nodevalue.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUComputeNodeValueInterpolation(const int nTotalNode, const int nBoundFace, const int nTotalFace,
                                              const int nTotalCell, const int nEquation, const int nIPoint);
    void CallGPUMultiProcessReset(const int nIPoint, const int nTotalNode, const int nEquation);
    void CallGPUComputeNodeValue(Grid *gridIn, int nEquation);
    void CallGPUModifyNodeValue(Grid *gridIn, int nEquation);

    __global__ void GPUModifyqNodetNode(const int nTotalNode, const int nEquation, RFloat *qNode, RFloat *tNode,
                                        RDouble *nodeWeight);

    __global__ void GPUnodeBCForSoildSurface(const int nBoundFace, const int *node_number_of_each_face,
                                             const int *face2node, const long long int *nodepos, int *d_nodeBC,
                                             const int *boundaryType, const int SOLID_SURFACE);
    __global__ void GPUnodeBCForFarField(const int nBoundFace, const int *node_number_of_each_face,
                                         const int *face2node, const long long int *nodepos, int *d_nodeBC,
                                         const int *boundaryType, const int SOLID_SURFACE, const int FARFIELD);
    __global__ void GPUnodeBCForOtherBC(const int nBoundFace, const int *node_number_of_each_face, const int *face2node,
                                        const long long int *nodepos, int *d_nodeBC, const int *boundaryType,
                                        const int SOLID_SURFACE, const int FARFIELD, const int SYMMETRY,
                                        const int INTERFACE);
    __global__ void GPUqOnBCFaceIni(const int nEquation, const int nBoundFace, const int nTotalCell,
                                    const int *left_cell_of_face, const int *right_cell_of_face, RFloat *d_qOnBCFace,
                                    const RFloat *q, const RFloat *t, const int *boundaryType, const int SYMMETRY,
                                    const int INTERFACE);
    __global__ void GPUqtWeightNodeForSolidSurface(const int nTotalNode, const int nEquation, const int nBoundFace,
                                                   const int *node_number_of_each_face, RFloat *d_qOnBCFace,
                                                   const int *face2node, const long long int *nodepos, RFloat *qNode,
                                                   RFloat *tNode, RDouble *nodeWeight, const int *boundaryType,
                                                   const int SOLID_SURFACE, const RDouble *d_x, const RDouble *d_y,
                                                   const RDouble *d_z, const RDouble *d_xfc, const RDouble *d_yfc,
                                                   const RDouble *d_zfc);
    __global__ void GPUqtWeightNodeForFarField(const int nTotalNode, const int nEquation, const int nBoundFace,
                                               const int *node_number_of_each_face, RFloat *d_qOnBCFace,
                                               const int *face2node, const long long int *nodepos, RFloat *qNode,
                                               RFloat *tNode, RDouble *nodeWeight, const int *d_nodeBC,
                                               const int *boundaryType, const int SOLID_SURFACE, const int FARFIELD,
                                               const RDouble *d_x, const RDouble *d_y, const RDouble *d_z,
                                               const RDouble *d_xfc, const RDouble *d_yfc, const RDouble *d_zfc);

    __global__ void GPUqtWeightNodeForOtherBC(const int nTotalNode, const int nEquation, const int nBoundFace,
                                              const int *node_number_of_each_face, RFloat *d_qOnBCFace,
                                              const int *face2node, const long long int *nodepos, RFloat *qNode,
                                              RFloat *tNode, RDouble *nodeWeight, const int *d_nodeBC,
                                              const int *boundaryType, const int SOLID_SURFACE, const int FARFIELD,
                                              const int SYMMETRY, const int INTERFACE, const RDouble *d_x,
                                              const RDouble *d_y, const RDouble *d_z, const RDouble *d_xfc,
                                              const RDouble *d_yfc, const RDouble *d_zfc);
    __global__ void GPUqtWeightNodeForInterioFace(const int nTotalNode, const int nEquation, const int nBoundFace,
                                                  const int nTotalCell, RFloat *qNode, RFloat *tNode,
                                                  RDouble *nodeWeight, const int *d_nodeBC, RFloat *q, RFloat *t,
                                                  const int *nodeNumberOfEachCell, const int *cell2Node,
                                                  const int *cell2NodePosition, const RDouble *d_x, const RDouble *d_y,
                                                  const RDouble *d_z, const RDouble *d_xcc, const RDouble *d_ycc,
                                                  const RDouble *d_zcc);

    __global__ void GPUNCountReset(const int nIPoint, int *nCount, const int *interPoint2GlobalPoint,
                                   const int *cellNumberOfInterPoint);
    __global__ void GPUqNodetNodeReset(const int nIPoint, const int nEquation, const int nTotalNode, RFloat *qNode,
                                       RFloat *tNode, const int *interPoint2GlobalPoint,
                                       const int *cellNumberOfInterPoint, const int *labelOfInterPoint);
    __global__ void GPUqIPtIPSetZero(const int nIPoint, const int nEquation, const int nTotalNode, RFloat *qInterPoint,
                                     RFloat *tInterPoint, const int *interPoint2GlobalPoint,
                                     const int *cellNumberOfInterPoint, const int *labelOfInterPoint);
    __global__ void GPUNcountQNodeTNodeInit(const int nTotalNode, const int nEquation, RFloat *qNode, RFloat *tNode,
                                            int *nCount);
    __global__ void GPUNodeWeightQNodeTNodeInit(const int nTotalNode, const int nEquation, RFloat *qNode, RFloat *tNode,
                                                RFloat *nodeWeight);

    void CallGPUComputeNodeValueBoundaryTreatMent(const int nEquation, const int nTotalNode, const int nBoundFace,
                                                  const int nTotalCell);

    __global__ void GPUComputeNodeValueBoundaryInit(const int nEquation, const int nBoundFace, const int nTotalNode,
                                                    const int *face2node, const int *node_number_of_each_face,
                                                    const int *boundaryType, const long long int *nodepos,
                                                    RFloat *qNode, RFloat *tNode, int *nCount, const int bctype_in);

    __global__ void GPUComputeNodeValueBoundaryCal(const int nTotalNode, const int nEquation, const int nBoundFace,
                                                   const int nTotalCell, const int *left_cell_of_face,
                                                   const int *face2node, const int *node_number_of_each_face,
                                                   const int *boundaryType, const long long int *nodepos,
                                                   const RFloat *q, const RFloat *t, RFloat *qNode, RFloat *tNode,
                                                   int *nCount, const int bctype_in);

    __global__ void GPUBoundaryFaceNCountQNodeTNodeCal(const int nTotalNode, const int nEquation, const int nBoundFace,
                                                       const int nTotalCell, const int *left_cell_of_face,
                                                       const int *right_cell_of_face,
                                                       const int *node_number_of_each_face, const int *face2node,
                                                       const long long int *nodepos, const RFloat *q, const RFloat *t,
                                                       RFloat *qNode, RFloat *tNode, int *nCount,
                                                       const int *boundaryType);

    __global__ void GPUCellLoopQNodeCalFinalSep(const int nTotalNode, const int equationID, const int nTotal,
                                                const int nTotalCell, const int nTotalFace,
                                                const int *nodeNumberOfEachCell, const int *cell2Node,
                                                const int *cell2NodePosition, const int *cell2NodeCount,
                                                const RFloat *qNS, RFloat *qNode);
    __global__ void GPUCellLoopTNodeCalFinalSep(const int nTotalNode, const int nTotal, const int nTotalCell,
                                                const int nTotalFace, const int *nodeNumberOfEachCell,
                                                const int *cell2Node, const int *cell2NodePosition,
                                                const int *cell2NodeCount, const RFloat *tCell, RFloat *tNode);
    __global__ void GPUCellLoopNCountCalFinalSep(const int nTotalNode, const int nTotal, const int nTotalCell,
                                                 const int nTotalFace, const int *nodeNumberOfEachCell,
                                                 const int *cell2Node, const int *cell2NodePosition,
                                                 const int *cell2NodeCount, int *nCount);

    __global__ void GPUCellLoopNCountQNodeTNodeCalFinal(const int nTotalNode, const int nEquation, const int nTotal,
                                                        const int nTotalCell, const int nTotalFace,
                                                        const int *nodeNumberOfEachCell, const int *cell2Node,
                                                        const int *cell2NodePosition, const int *cell2NodeCount,
                                                        const RFloat *qNS, const RFloat *tCell, RFloat *qNode,
                                                        RFloat *tNode, int *nCount);

    __global__ void GPUInteriorFaceNCountQNodeTNodeCal(const int nTotalNode, const int nEquation, const int nBoundFace,
                                                       const int nTotalCell, const int nTotalFace,
                                                       const int *left_cell_of_face, const int *right_cell_of_face,
                                                       const int *node_number_of_each_face, const int *face2node,
                                                       const long long int *nodepos, const RFloat *q, const RFloat *t,
                                                       RFloat *qNode, RFloat *tNode, int *nCount);

    __global__ void GPUComputeNodeValueAvr(const int nEquation, const int nTotalNode, RFloat *qNode, RFloat *tNode,
                                           int *nCount);

    void CallGPUModifyNodeValue(const int nIPoint, const int nEquation, const int nTotalNode);

    __global__ void GPUModifyNodeValue(const int nIPoint, const int nEquation, const int nTotalNode, RFloat *qNode,
                                       RFloat *tNode, RFloat *qInterPoint, RFloat *tInterPoint,
                                       const int *interPoint2GlobalPoint);

    void CallGPUModifyQTurbNodeValue(const int nEquation, const int nTotalNode);

    __global__ void GPUModifyQTurbNodeValue(const int nEquation, const int nTotalNode,
                                            RFloat *qTurbNode, const int *ncount);
} //！ namespace GPUKernels
