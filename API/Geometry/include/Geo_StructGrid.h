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
//! @file      Geo_StructGrid.h
//! @brief     It is type of Structured geometry grid.
//!            The inheriting order is: SimpleGrid->Grid->StructuredGrid/UnstructuredGrid
//! @author    He Xin, Bell, Min Yaobing.

#pragma once
#include "Geo_Grid.h"
#include "Geo_UnstructGrid.h"
#include "Geo_OversetGridTopo_Struct.h"
#include "Geo_NodeTopo_Struct.h"
#include "Geo_FaceMetrics_Struct.h"
#include "Geo_CellMetrics_Struct.h"
#include "Geo_MultiGridInfo_Struct.h"
#include "Geo_DynamicGridMetrics_Struct.h"
#include "IO_VirtualFile.h"
#include "Geo_StructBC.h"

namespace PHSPACE
{
class HoleGrid;
class BackgroundTree;
class OversetCellCollector;
class UnstructGrid;

class VirtualFile;
class StructBCSet;
class StructBC;

class StructGrid : public Grid
{
private:
    //! Geo_NodeTopo_Struct class defines the node dimension and coordinates.
    //! 1: The dimension of three direction (ni, nj ,nk).
    //! 2: The node coordinates array of three direction (structx, structy, structz).
    Geo_NodeTopo_Struct *nodeTopology;

    //! Geo_FaceMetrics_Struct class defines the face metrics of structured grid.
    //! 1: unit face normal (xfn, yfn, zfn).
    //! 2: face area (area).
    //! 3: face vector (xSurfaceVector, ySurfaceVector, zSurfaceVector).
    Geo_FaceMetrics_Struct *faceMetrics;

    //! Geo_CellMetrics_Struct class defines the cell metrics of structured grid.
    //! 1: cell center (xcc).
    //! 2: cell volume (vol).
    //! 3: cell length (xlen).
    Geo_CellMetrics_Struct *cellMetrics;

    //! Multi-grid step information.
    //! 1:istp, jstp, kstp.
    Geo_MultiGridInfo_Struct *multigridInfo;

    //! Geo_DynamicGridMetrics_Struct class defines the dynamic grid metrics of the structured grid.
    //! 1: Cell volume of the last time step (voln).
    //! 2: Face normal velocity (vgn).
    //! 3: Face velocity (xfv, yfv, zfv).
    Geo_DynamicGridMetrics_Struct *dynamicGridMetrics;

    //! Distance to wall
    RDouble3D *walldist;

    //! Normal direction of the wall face nearest to a cell.
    RDouble3D *nearestwallfacenormalx;
    RDouble3D *nearestwallfacenormaly;
    RDouble3D *nearestwallfacenormalz;

    //! overset iblank
    Int3D *iBlank3D;

    //! If the cell is the firse cell near wall.
    Int3D *wallCellLabel;

    Int3D *cellBoundaryType;

    //! Boundary information.
    //! 1: Number of boundary face (nBoundFace).
    //! 2: Number of interface face (nIFace).
    StructBCSet *structBCSet;

    //! Topology information of the overset grid.
    Geo_OversetGridTopo_Struct *oversetGridTopology;

    //! 
    PHDouble2D localFrame;

    //! Probes nearest node dimension of i,j,k directions.
    vector <int> zoneProbesCellNI;
    vector <int> zoneProbesCellNJ;
    vector <int> zoneProbesCellNK;

    UnstructGrid *gridUnstr;

    //! Ordinary grid information.
    int ordinaryGridIndex;
    int ordinaryDimStartIndex[3], ordinaryDimEndIndex[3];

    RDouble partialCFL;

public:
    LIB_EXPORT StructGrid();
    LIB_EXPORT ~StructGrid();

public:
    //! Initialize the Grid, create nodeTopology, faceMetrics, et al.
    LIB_EXPORT void InitGrid(GridID *index, int level, int dim, int type);

    //! Return the wall distance (walldist) pointer.
    RDouble3D * GetWallDist() const;

    //! Return the nearest wall face normal x pointer.
    RDouble3D * GetNearestWallFaceNormalX() const;

    //! Return the nearest wall face normal y pointer.
    RDouble3D * GetNearestWallFaceNormalY() const;

    //! Return the nearest wall face normal z pointer.
    RDouble3D * GetNearestWallFaceNormalZ() const;

    //! Return the wall cell label for the first layer of wall boundary.
    Int3D * GetWallCellLabel();

    void InitWallDist();

    void InitIBlank3D();

    Int3D * GetIBlank3D();

    void CopyIBlank3D(int *iBlank);
    void SetIBlank(int *iBlank);

    UnstructGrid * GetUnstrGrid() { return this->gridUnstr; }
    void SetUnstrGrid(UnstructGrid *gridIn) { this->gridUnstr = gridIn; }
    //! 
    void DecodeWallDist(DataContainer *cdata);

    // ! Get the i,j,k when map current grid to the level layer coarse grid.
    // ! @param[in]  level    level layer coarse grid.
    // ! @param[out] i,j,k    number of the grid point in  three direction.
    LIB_EXPORT void RemapMultigridIJK(int level, int &i, int &j, int &k);

    //! Return the cell center (xcc) pointer.
    RDouble3D * GetCellCenterX() const;

    //! Return the cell center (ycc) pointer.
    RDouble3D * GetCellCenterY() const;

    //! Return the cell center (zcc) pointer.
    RDouble3D * GetCellCenterZ() const;

    //! Return the UNIT face normal (xfn) pointer.
    RDouble4D * GetFaceNormalX() const;

    //! Return the UNIT face normal (yfn) pointer.
    RDouble4D * GetFaceNormalY() const;

    //! Return the UNIT face normal (zfn) pointer.
    RDouble4D * GetFaceNormalZ() const;

    //! Return the face area (area) pointer.
    RDouble4D * GetFaceArea() const;

    //! Return the face vector (xFaceVector) pointer.
    RDouble4D * GetFaceVectorX() const;

    //! Return the face vector (yFaceVector) pointer.
    RDouble4D * GetFaceVectorY() const;

    //! Return the face vector (zFaceVector) pointer.
    RDouble4D * GetFaceVectorZ() const;

    //! Return the cell volume (vol) pointer.
    RDouble3D * GetCellVolume() const;

    //! 
    RDouble3D * GetCellLengthX() const;

    //! 
    RDouble3D * GetCellLengthY() const;

    //!
    RDouble3D * GetCellLengthZ() const;

    //! Return the cell volume pointer of the last time step (voln).
    RDouble3D * GetCellVolumeOld() const;

    //! Return the X velocity pointer of cell center.
    RDouble3D *GetCellVelocityX() const;

    //! Return the Y velocity pointer of cell center.
    RDouble3D *GetCellVelocityY() const;

    //! Return the Z velocity pointer of cell center.
    RDouble3D *GetCellVelocityZ() const;

    //! Return the face normal velocity (vgn) pointer.
    RDouble4D * GetFaceNormalVelocity() const;

    RDouble5D * GetFaceVector_FD() const;
    RDouble5D * GetFaceNormal_FD() const;
    RDouble3D * GetCellJacobian() const;

    //! Return the node's x coordinates array (structx) pointer.
    RDouble3D * GetStructX() const;

    //! Return the node's y coordinates array (structx) pointer.
    RDouble3D * GetStructY() const;

    //! Return the node's z coordinates array (structx) pointer.
    RDouble3D * GetStructZ() const;

    //! Set the node dimension of I direction (ni).
    void SetNI(int ni);

    //! Set the node dimension of J direction (nj).
    void SetNJ(int nj);

    //! Set the node dimension of K direction (nk).
    void SetNK(int nk);

    //! Return the node dimension of I direction (ni).
    int GetNI() const;

    //! Return the node dimension of J direction (nj).
    int GetNJ() const;

    //! Return the node dimension of K direction (nk).
    int GetNK() const;

    //! Set the probes cell dimension of I direction (ni).
    void SetZoneProbesCellNI(vector <int> zoneProbesCellNI);

    //! Set the probes cell dimension of J direction (nj).
    void SetZoneProbesCellNJ(vector <int> zoneProbesCellNJ);

    //! Set the probes cell dimension of K direction (nk).
    void SetZoneProbesCellNK(vector <int> zoneProbesCellNK);

    //! Return the probes cell dimension of I direction (ni).
    vector <int> GetZoneProbesCellNI() const;

    //! Return the probes cell dimension of J direction (nj).
    vector <int> GetZoneProbesCellNJ() const;

    //! Return the probes cell dimension of K direction (nk).
    vector <int> GetZoneProbesCellNK() const;

    //! Return the cell iteration index, for example: if layer = 0, then iCellStart = 1, iCellEnd = ni-1.
    void GetCellIterationIndex(int &iCellStart, int &iCellEnd, int &jCellStart, int &jCellEnd, int &kCellStart, int &kCellEnd, int layer = 0);

    //! Return the node iteration index, for example: iNodeStart = 1, iNodeEnd = ni.
    void GetNodeIterationIndex(int &iNodeStart, int &iNodeEnd, int &jNodeStart, int &jNodeEnd, int &kNodeStart, int &kNodeEnd);

    //! Return the face iteration index, for example: if nSurface = 1, then iFaceStart = 0, iFaceEnd = ni-1.
    void GetFaceIterationIndex(int &iFaceStart, int &iFaceEnd, int &jFaceStart, int &jFaceEnd, int &kFaceStart, int &kFaceEnd, int iSurface);

    //! Get the face index.
    void GetNsurfIndex(int &il1, int &jl1, int &kl1, int iDimension);

    //! Get the left cell list of each face.
    void GetLeftCellOfFace(int &i, int &j, int &k, int &il1, int &jl1, int &kl1, int &ile, int &jle, int &kle);

    //! Get the right cell list of each face.
    void GetRightCellOfFace(int &i, int &j, int &k, int &il1, int &jl1, int &kl1, int &ire, int &jre, int &kre);

    //! Assign the node dimension of I, J, K direction to si, sj ,sk.
    void GetND(int &si, int &sj, int &sk) const;

    //! Return multigrid step I (istp);
    int GetMultigridStepI() const;

    //! Return multigrid step J (jstp);
    int GetMultigridStepJ() const;

    //! Return multigrid step K (kstp);
    int GetMultigridStepK() const;

    //! Sets multigrid step I (istp);
    void SetMultigridStepI(int istp);

    //! Sets multigrid step J (jstp);
    void SetMultigridStepJ(int jstp);

    //! Sets multigrid step K (kstp);
    void SetMultigridStepK(int kstp);

    //! Copy coordinates(X, Y, Z) and boundary information from another grid.
    LIB_EXPORT void CopyGrid(StructGrid *grid);

    //! Copy boundary information from another grid.
    LIB_EXPORT void CopyStructBCSet(StructBCSet *structBCSet);

    //! Copy coordinates(X, Y, Z) information from another grid.
    LIB_EXPORT void CopyXYZ(int ni, int nj, int nk, RDouble *xx, RDouble *yy, RDouble *zz);

    //! Return cell volume reference.
    LIB_EXPORT RDouble3D & GetVolume(int istep = 0) const;

    //! Set the three dimensional array of coordinates from one dimensional coordinates.
    LIB_EXPORT void SetArrayLayout();

    //! Allocate memory for face metrics, cell metrics (dynamic grid metrics alternatively).
    LIB_EXPORT void AllocateMetrics(ActionKey *actkey);
    
    //! Compute metrics (face normal, cell center, et al.).
    LIB_EXPORT void ComputeMetrics(ActionKey *actkey = 0);
    
    //! Create composite boundary condition region (allocate memory for structBCSet).
    LIB_EXPORT void CreateCompositeBCRegion(int nBCRegion);

    //! Return composite boundary condition region pointer (structBCSet).
    StructBCSet * GetStructBCSet() const;

    //! Set boundary face information.
    LIB_EXPORT void SetBCFaceInfo();

    //! Get source boundary grid's I, J, K index.
    LIB_EXPORT void GetSourceIndexIJK(int iFace, int ipos, int &i, int &j, int &k);

    //! Get source boundary grid's I, J, K index.
    LIB_EXPORT void GetSourceIndexIJK_Nsurf_LR(int iFace, int ipos, int &i, int &j, int &k, int &nsurf, int &s_lr);

    //! Get target boundary grid's I, J, K index.
    LIB_EXPORT void GetTargetIndexIJK(int iFace, int ipos, int &i, int &j, int &k);

    LIB_EXPORT void GetSourceIndexIJK_fornodetransfer(int iFace, int &id, int &jd, int &kd, int &i_lr, int &j_lr, int &k_lr);
    LIB_EXPORT void GetTargetIndexIJK_fornodetransfer(int iFace, int &id, int &jd, int &kd, int &i_lr, int &j_lr, int &k_lr);
    void ComputeMetrics_BC();
    void ComputeMetrics_BC_2D();
    void ComputeMetrics_BC_3D();
    void ComputeMetricsWithBCforStructHighOrder();
    void ComputeMetricsWithBCforStructHighOrder2D();
    void ComputeMetricsWithBCforStructHighOrder3D();
    void ComputeJacobianforStructHighOrder();
    void ComputeJacobianforStructHighOrder2D();
    void ComputeJacobianforStructHighOrder3D();

    void UploadGridInfoOnlyforStructHighOrder(DataContainer *&dataContainer, const int &neighborZoneIndex, int mission);
    void DownloadGridInfoOnlyforStructHighOrder(DataContainer *&dataContainer, const int &neighborZoneIndex, int mission);
    void UploadGridInfoOnlyforStructHighOrder_mission1(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void UploadGridInfoOnlyforStructHighOrder_mission2(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void UploadGridInfoOnlyforStructHighOrder_mission3(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void UploadGridInfoOnlyforStructHighOrder_mission4(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void UploadGridInfoOnlyforStructHighOrder_mission5(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void UploadGridInfoOnlyforStructHighOrder_mission6(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void UploadGridInfoOnlyforStructHighOrder_mission7(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void UploadGridInfoOnlyforStructHighOrder_mission8(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void UploadGridInfoOnlyforStructHighOrder_mission9(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void UploadGridInfoOnlyforStructHighOrder_mission10(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void UploadGridInfoOnlyforStructHighOrder_mission11(DataContainer *&dataContainer, const int &neighborZoneIndex);

    void DownloadGridInfoOnlyforStructHighOrder_mission1(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void DownloadGridInfoOnlyforStructHighOrder_mission2(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void DownloadGridInfoOnlyforStructHighOrder_mission3(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void DownloadGridInfoOnlyforStructHighOrder_mission4(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void DownloadGridInfoOnlyforStructHighOrder_mission5(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void DownloadGridInfoOnlyforStructHighOrder_mission6(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void DownloadGridInfoOnlyforStructHighOrder_mission7(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void DownloadGridInfoOnlyforStructHighOrder_mission8(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void DownloadGridInfoOnlyforStructHighOrder_mission9(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void DownloadGridInfoOnlyforStructHighOrder_mission10(DataContainer *&dataContainer, const int &neighborZoneIndex);
    void DownloadGridInfoOnlyforStructHighOrder_mission11(DataContainer *&dataContainer, const int &neighborZoneIndex);
         
    Int1D *DiscretePrecisionKXI;
    Int1D *DiscretePrecisionETA;
    Int1D *DiscretePrecisionCTA;

    void InitialDiffOrderforStructHighOrder();
    void CorrectDiffOrderNearUnstructBlockforStructHighOrder(const int &neighborZoneIndex);

    //! Read x,y,z coordinates from file to construct node topology (structx, structy, structz).
    LIB_EXPORT void ReadXYZ(fstream &file);

    //! Write node topology (structx, structy, structz) from range i(j, k)st to i(j, k)ed to file.
    LIB_EXPORT void WriteXYZ(fstream &file, int ist, int ied, int jst, int jed, int kst, int ked);

    //! Compute and set nTotalNode, nTotalFacenTFace, nTCell.
    LIB_EXPORT void SetBasicDimension();

    LIB_EXPORT void SetBasicDimension(int ni, int nj, int nk);

    //! Compute face center ? Yes, Verified by Li Peng on Mar 13,2019.
    LIB_EXPORT void FaceCoor(int i, int j, int k, int nsurf, RDouble &xc, RDouble &yc, RDouble &zc);

    //! Compute cell center ? Yes, Verified by Li Peng on Mar 13,2019.
    LIB_EXPORT void CenterCoor(int i, int j, int k, RDouble &xc, RDouble &yc, RDouble &zc);

    //! Compute average norm of q_0_proxy and q_1_proxy field's difference.
    LIB_EXPORT void CVGNorm(FieldProxy *q_1_proxy, FieldProxy *q_0_proxy, int neqn, RDouble &norm);

    //! 
    LIB_EXPORT void ActionReflect(ActionTag *acttag);

    //! Grid based actions.
    LIB_EXPORT void Action(ActionKey *actkey);

    //! Action of trading cell center.
    LIB_EXPORT void TranslateAction(ActionKey *actkey);

    //! Total length of data need to be translated.
    LIB_EXPORT streamsize TranslateActionLength(ActionKey *actkey);

    //! CoarseGrids.
    LIB_EXPORT void CoarseGrids(int nlmax);

    //! Process boundary condition information.
    LIB_EXPORT void ProcessBCInfo();
    LIB_EXPORT void SetLnkInfo();

    //! Read grid information from file.
    LIB_EXPORT void ReadGrid(fstream &file);

    //! Write grid information to file
    LIB_EXPORT void WriteGrid(fstream &file);

    //! Encode data from datacontainer to virtual file.
    LIB_EXPORT void Encode(DataContainer *cdata, int flag);
    
    //! Decode data from virtual file to datacontainer.
    LIB_EXPORT void Decode(DataContainer *cdata, int flag);

    //! Revise the zone index of neighbor zones.
    LIB_EXPORT void ReviseNeighborZoneIndex(int zoneStart);
    
    //! Update coordinate data to datacontainer.
    LIB_EXPORT void UpdateCoordinate(DataContainer *cdata);

    //! Get grid's min and max distance.
    LIB_EXPORT void GetMinMaxDS(RDouble &dismin, RDouble &dismax);

    //! Get boundary grid's min and max distance.
    LIB_EXPORT void GetBoundaryMinMax(RDouble *pmin, RDouble *pmax);

    //! Establish grid link information (parallel computation).
    LIB_EXPORT void BuildGridLink(LinkStruct *link);

    //! Compute total number of Wall Type cell.
    LIB_EXPORT int GetNumberOfWallCell();

    //! To obtain the information of wall boundary conditions.
    //! param[out]: wallBCNum indicates the number of the wall boundary conditions;
    //! param[out]: maxCellNum denotes the maximum number of the cells in the wall boundary conditions.
    LIB_EXPORT void GetWallBCInfo(int &wallBCNum, int &maxCellNum);

    //! Dump to tecplot.
    LIB_EXPORT void ToTecplot();

    //! Get surface velocity.
    LIB_EXPORT void GridSurfaceVelocity(RDouble *xt, RDouble *yt, RDouble *zt);

    //! Get cell velocity.
    LIB_EXPORT void GridCellVelocity(RDouble *xt, RDouble *yt, RDouble *zt);
    
    LIB_EXPORT void ComputeGridFaceVelocity();
    //! Change BC type.
    LIB_EXPORT void ChangeBCType(int from_bctype, int to_bctype);

    //! Compute variables at nodes.
    LIB_EXPORT void CompNodeVar(RDouble4D &q_n, RDouble4D &q);
    
    //! Allocate memory for walldist.
    LIB_EXPORT void AllocateWalldist();

    //! Get the largest local grid length, distance of neighbor cell centers.
    LIB_EXPORT RDouble3D * GetLargestLocalGridLength();

    //! Get the smallest local grid length, distance of neighbor cell centers.
    LIB_EXPORT RDouble3D * GetSmallestLocalGridLength();

    //! Get the sub-grid length, used for DES simulation ONLY.
    LIB_EXPORT RDouble3D * GetSubgridLength();
    
    //! Compute the minimum distance of the grid edges.
    LIB_EXPORT RDouble CalMinEdgeLength();

    //! Compute the maximum distance of the grid edges.
    LIB_EXPORT RDouble CalMaxEdgeLength();

    //! Write face boundary condition to file stream.
    LIB_EXPORT void WriteFaceBC(fstream &file);
#pragma warning(disable:4100)
    LIB_EXPORT RDouble SkewnessSummary(ActionKey *actkey = 0) { return 1.0; }
#pragma warning(default:4100)
    
    void GetLocalIndex(int gp, int &ip, int &jp, int &kp) const;
    void GetLocalCenter(int gp, int &ip, int &jp, int &kp) const;
    LIB_EXPORT void EncodeOversetGrid(DataContainer *cdata);
    LIB_EXPORT void DecodeOversetGrid(DataContainer *cdata);
    LIB_EXPORT void AllocateOversetGrid();
    LIB_EXPORT void CopyOversetFromUnstructGrid();
    LIB_EXPORT int * GetCellTypeContainer();
    LIB_EXPORT void CreateLinkPointStructure();
    LIB_EXPORT void SetLinkPointLabel(int *iLabel, int iDirection);
    //LIB_EXPORT void DeleteLinkPointStructure();
    LIB_EXPORT void SetZoneStartPointLabel(int zoneStartPointLabel);
    //LIB_EXPORT void DeleteCellContainer();
    int GetZoneStartPointLabel() const;
    int * GetILinkPointLabel() const;
    int * GetJLinkPointLabel() const;
    int * GetKLinkPointLabel() const;
    LIB_EXPORT void GatherCellsByColor(int &pCell, int color);
    LIB_EXPORT void ResetCellType();
    LIB_EXPORT void ProbeExternalOverlapCells(int &counter);
    LIB_EXPORT void ProbeHoleCell(HoleGrid *holeGrid);
    LIB_EXPORT void ProbeBorderCell(int iType, int jType, vector< int > *friendPoint, vector< int > *friendZone);
    LIB_EXPORT void ColorOverlapCell(int gp, int type);
    LIB_EXPORT void SetZoneStartCenterLabel(int zoneStartCenterLabel);
    int GetNumberOfCores() const;
    LIB_EXPORT int * GetHingedPointContainer() const;
    LIB_EXPORT RDouble * GetXCoreContainer() const;
    LIB_EXPORT RDouble * GetYCoreContainer() const;
    LIB_EXPORT RDouble * GetZCoreContainer() const;
    //LIB_EXPORT void InitCoreParameterContainer();
    LIB_EXPORT void CreateCoreParameterContainer();
    //LIB_EXPORT void DeleteCoreParameterContainer();
    LIB_EXPORT void SetHingedPointContainer(int *iLabel);
    LIB_EXPORT void SetCoreCoordinateContainer(RDouble *sLabel, int iDirection);
    LIB_EXPORT PHDouble2D & GetLocalFrame();    //! Frequent calls will reduce the computational efficiency, it will be optimized later.
    LIB_EXPORT void ComputeWeight();

    LIB_EXPORT void ComputeUnitCellCapsule(BackgroundTree *backgroundTree, int &iUnitCell);
    LIB_EXPORT void CollectOversetCells(OversetCellCollector *oversetCellCollector, int myid);
    LIB_EXPORT void ProbeNearestPoint(RDouble xx, RDouble yy, RDouble zz, vector< RDouble > &ss, vector< int > &gp);
    LIB_EXPORT void ComputeLocalCoordinateByNewtonMethod(int iCore, int jCore, int kCore, RDouble xx, RDouble yy, RDouble zz, vector< RDouble > &ss, vector< int > &node);
    
    void DecodeBCFace(DataContainer *cdata);
    void DecodeBcDirFace(DataContainer *cdata);

    void GhostCell3DExceptInterface(RDouble3D &TargetArray3D);
    void GhostCell3DExceptInterface(RDouble4D &TargetArray4D, const int &nVariables);

    void InitVariableWallTemperature();

    RDouble4D * GetFaceVelocityX() const;
    RDouble4D * GetFaceVelocityY() const;
    RDouble4D * GetFaceVelocityZ() const;

    void GetBcNamePair(set< pair<int, string> > &bcNamePairSet);
    void DumpWallFaceCenter(ActionKey *actkey);
    void GetCellIBlank(ActionKey *actkey);
    void TranslateCellIBlank(ActionKey *actkey);
    void ComputeCellBoundaryType();
    Int3D *GetCellBoundaryType();

    void SetOrdinaryGridIndex(int ordinaryGridIndexIn);
    void SetOrdinaryDimStartIndex(int *ordinaryDimStartIndexIn);
    void SetOrdinaryDimEndIndex  (int *ordinaryDimEndIndexIn);
    int  GetOrdinaryGridIndex();    
    int *GetOrdinaryDimStartIndex();
    int *GetOrdinaryDimEndIndex();
    void GetOrdinaryNodeIndex(int iIndexIn, int jIndexIn,int kIndexIn,int &iOrdinaryIndexOut,int &jOrdinaryIndexOut,int &kOrdinaryIndexOut);
    void GetOrdinaryCellIndex(int iIndexIn, int jIndexIn,int kIndexIn,int &iOrdinaryIndexOut,int &jOrdinaryIndexOut,int &kOrdinaryIndexOut);

    RDouble  GetPartialCFL();
    void SetPartialCFL(RDouble PartialCFL){ partialCFL = PartialCFL;}

    void CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, const int &neighborZoneIndex, const int &neqn);
    void DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, const int &neighborZoneIndex, const int &neqn);

    void ComputeBCOriginal();

private:
    void ReadXYZ(fstream &file, RDouble3D &coor, int ist, int ied, int jst, int jed, int kst, int ked);
    void ReadXYZASCII(fstream &file, RDouble3D &coor, int ist, int ied, int jst, int jed, int kst, int ked);
    void WriteXYZ(fstream &file, RDouble3D &coor, int ist, int ied, int jst, int jed, int kst, int ked);
    void AllocateMetricsALE(ActionKey *actkey);
    void ComputeMetrics3D(ActionKey *actkey);
    void ComputeMetrics2D(ActionKey *actkey);
    void ComputeCellCenter();

    //! @brief grid metrics computation for struct high order solver.
    void AllocateMetricsStructHighOrder(ActionKey *actkey);
    void ComputeMetricsStructHighOrder2D(ActionKey *actkey);
    void ComputeMetricsStructHighOrder3D(ActionKey *actkey);
    void ComputeMetricsStructHighOrder2D_HDCS();
    void ComputeMetricsStructHighOrder3D_HDCS();
    void GhostMetricsStructHighOrder2D();
    void GhostMetricsStructHighOrder3D();

    //! @brief grid difference operator for a whole line.
    void GridDelta(const Range &NP, const Range &NI, RDouble2D *LineData1, RDouble2D *LineData2);

    void ReadMovingGrid();
    void UpdateVolold();
    void ReadGrid(VirtualFile *vfile);
    void WriteGrid(VirtualFile *vfile);
    void EncodeGrid(DataContainer *cdata);
    void DecodeGrid(DataContainer *cdata);
    void DumpWalldist(ActionKey *actkey);
    void ReadWalldist(RDouble *wallDistIn);
    void ReadNearestwallfacenormalx(RDouble *nearestwallfaceNormalxIn);
    void ReadNearestwallfacenormaly(RDouble *nearestwallfaceNormalyIn);
    void ReadNearestwallfacenormalz(RDouble *nearestwallfaceNormalzIn);
    void GetCellCenter(ActionKey *actkey);
    void TranslateCellCenter(ActionKey *actkey);
    bool CommunicateNodeLocation();
    void GetNodeIJKfromSourceIndex(StructGrid *fgrid, int iFace, int layer, int IJKofNodes[4][3]);
    void GetNodeIJKfromTargetIndex(StructGrid *fgrid, int iFace, int layer, int IJKofNodes[4][3]);
    void GetRelationBetweenTwoPairsofNodes(int RelationBetweenTwoPairsofNodes[4], RDouble Nodes_ABCD[4][3], RDouble Nodes_abcd[4][3]);
    void GetNodeIJKfromSourceIndexghostofghost(StructGrid *fgrid, int iFace, int layer, int IJKofNodes[12][3]);
    void GetNodeIJKfromTargetIndexghostofghost(StructGrid *fgrid, int iFace, int layer, int IJKofNodes[12][3]);
    int GetRelationshipBetweenTwoPairsofNodes(RDouble NodesofABCD[4][3], RDouble Nodesofabcd[4][3]);
    streamsize TranslateCellCenterLength(ActionKey *actkey);
    void GhostMetrics3D();
    void GhostMetrics2D();
    void SimpleAction(ActionKey *actkey);

    int CompNIFace();
    void FillCoreIndex(vector< int > &core_index, int iCore, int jCore, int kCore);
    void ComputePhysicalCoordinate(vector< vector< RDouble > > &amer, vector< int > &core_index, vector< RDouble > &xx, vector< RDouble > &yy,  vector< RDouble > &zz);
    void ComputePhysicalCoordinate(vector< vector< RDouble > > &localFrame, RDouble xx, RDouble yy, RDouble zz, vector< RDouble > &globalCoordinate);
    void FillCoefficientContainer(vector< vector< RDouble > > &source, vector< RDouble > &target, int i);
    void ComputeTempContainer(vector< RDouble > &ksai, vector< RDouble > &temp);
    void ComputeJacobiMatrix(vector< RDouble > &a, vector< RDouble > &b, vector< RDouble > &c, vector< RDouble > &ksai, vector< vector< RDouble > > &amer);
    void GaussElimination(vector< vector< RDouble > > &a, vector< RDouble > &b, int n, vector< RDouble > &x, int &flag);

    void ComputeLargestLocalGridLength();
    void ComputeSmallestLocalGridLength();
    void ComputeSubgridLength();
};

#include "Geo_StructGrid.hxx"

//! Communicate the interface variables at the boundary faces between different zones.
//! @param[in] grid    computational mesh.
//! @param[in] variable    variable pointer of which need to be communicated.
//! @param[in] variableName    variable name of which need to be communicated.
//! @param[in] variableDimension    variable dimension of which need to be communicated.
//! -# 1: 1-D array.
//! -# 2: 2-D array.
//! -# n: n-D array.
void CommunicateInterfaceValue(StructGrid *grid, RDouble4D *variable, const string &variableName, int variableDimension);
void CommunicateInterfaceValue(StructGrid *grid, RDouble3D *variable, const string &variableName);

void UploadInterfaceValue(StructGrid *grid, RDouble4D *f, const string &name, int neqn);
void UploadInterfaceValue(StructGrid *grid, RDouble3D *f, const string &name);
void DownloadInterfaceValue(StructGrid *grid, RDouble4D *f, const string &name, int neqn);
void DownloadInterfaceValue(StructGrid *grid, RDouble3D *f, const string &name);

void UploadOversetValue(StructGrid *grid, RDouble4D *field, const string &name, int neqn);
void DownloadOversetValue(StructGrid *grid, RDouble4D *field, const string &name, int neqn);
void UploadOversetValue(StructGrid *grid, RDouble3D *field, const string &name);
void DownloadOversetValue(StructGrid *grid, RDouble3D *field, const string &name);

void MinMaxDiffINF(StructGrid *grid, int ig, RDouble *dmin, RDouble *dmax, RDouble4D *q, int m);
void BarthLimiterINF (StructGrid *grid, int ig, RDouble *limit, RDouble4D *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int m);
void VencatLimiterINF(StructGrid *grid, int ig, RDouble *limit, RDouble4D *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int m);

void EncodeIJK(int &index, int i, int j, int k, int ni, int nj);
void DecodeIJK(int index, int &i, int &j, int &k, int ni, int nj);

void GetIJKRegion(const Range &I, const Range &J, const Range &K, int &ist, int &ied, int &jst, int &jed, int &kst, int &ked);
void GetNsurfIndex(int nsurf, int &i, int &j, int &k);

void GetIJKROT(int ijk[4][3], int *irot, int *jrot, int *krot);

extern int mide[3][3];

StructGrid * StructGridCast(Grid *grid_in);

void GetRange(int ni, int nj, int nk, int startShift, int endShift, Range &I, Range &J, Range &K);
void GetRange(Range &I, Range &J, Range &K, int ist, int ied, int jst, int jed, int kst, int ked);

//! Used by Post_Visual class.
void CompNodeVar(Grid *grid_in, RDouble4D &qn, int m, RDouble4D &q, int n, bool isVelocityForPostVisual = false);
void CompNodeVar(Grid *grid_in, RDouble4D &qn, int m, RDouble3D &q);

void ModifyWallTemperature(const string& name, Grid *grid_in, RDouble4D &qn);
void ModifyBoundaryNodeValue(const string &name, Grid *grid_in, RDouble4D &qn, RDouble4D &cellCenterData, int index, int bcType);
void ModifyBoundaryNodeValue(const string &name, Grid *grid_in, RDouble4D &qn, RDouble3D &cellCenterData, int bcType);
}
