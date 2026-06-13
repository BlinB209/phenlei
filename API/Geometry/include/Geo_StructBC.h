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
//! @file      Geo_StructBC.h
//! @brief     It is the class 'StructBC', which is type of Structured geometry grid operation for boundary condition.
//!            The inheriting order is: SimpleBC -> StructBC/UnstructBC.
//! @author    Zhang Yong, Bell, He Xin.

#pragma once
#include "LIB_Macro.h"
#include "TypeDefine.h"
#include "Geo_SimpleBC.h"
#include "PHMpi.h"
#include "Constants.h"
#include "TK_Log.h"
using namespace std;

namespace PHSPACE
{
//! @brief It defines the class 'StructBC', which is the base class of structured geometry grid operation for boundary condition.
//!        The inheriting order is: SimpleBC -> StructBC/UnstructBC.
class StructBC : public SimpleBC
{
public:
    //! Construct the StructBC.
    //! @param[in] zoneID      Zone Index.
    //! @param[in] regionID    Region Index.
    StructBC(int zoneID, uint_t regionID);

    ~StructBC();

private:
    //! Region index (the number of region).
    uint_t nr;

    //! s_nd: Dimensional number in the current  zone. s - source, nd - the number of dimension.
    //! t_nd: Dimensional number in the neighbor zone. t - target, nd - the number of dimension.
    //! -# 0: i direction.
    //! -# 1: j direction.
    //! -# 2: k direction.
    int s_nd, t_nd;

    //! s_lr : This BC region is on the left or right side in the current zone.  s - source, lr - Left/Right.
    //! t_lr : This BC region is on the left or right side in the neighbor zone. t - target, lr - Left/Right.
    //! -# -1: Left side.
    //! -#  1: Right side.
    int s_lr, t_lr;

    //! s_lr3d && t_lr3d stores the composite information of both s_nd/t_nd && s_lr/t_lr.
    //! It is used to uniform the loop format only. The values are equal to the following:
    //!     if (i == s_nd) s_lr3d[i] = s_lr;  else  s_lr3d[i] = 0;
    //!     if (i == t_nd) t_lr3d[i] = t_lr;  else  t_lr3d[i] = 0;
    int s_lr3d[3], t_lr3d[3];

    //! nbs: The current  zone index of this BC region. nb - the number of block, s - source.
    //! nbt: The neighbor zone index of this BC region. nb - the number of block, t - target.
    int nbs, nbt;

    //! s_st: Start local point index in the current zone. s - source, st - start.
    //! s_ed: End   local point index in the current zone. s - source, ed - end.
    int s_st[3], s_ed[3];

    //!
    int lnkInfo[12];

    //! s_dir3d: the i\j\k direction between the current face and matching target face if the face is a interface.
    int s_dir3d[3];

    //! WhetherMerge:the symbol of whether has been merged.
    int WhetherMerge;
    
    //! stOrignal: Start local point index in the current block. s - source, st - start.
    //! edOrignal: End   local point index in the current block. s - source, ed - end.
    int stOrignal[3], edOrignal[3];

public:
    StructBC & operator = (const StructBC &rhs);

public:
    //! Copy struct boundary.
    void CopyBC(const StructBC *bcNew);

    //!Set the current  zone index of this BC region.
    void SetCurrentZoneIndex(int nbs) { this->nbs = nbs; }

    //! Set the symbol of whether has been merged.
    void SetWhetherMerge(int WhetherMerge);

    //! Get the symbol of whether has been merged.
    int GetWhetherMerge();

    //! Set the start local point index in the current block. s - source, st - start.
    //! Set the end   local point index in the current block. s - source, ed - end.
    void SetStartOrignal(int *st);
    void SetEndOrignal(int *end);

    //! Get the start local point index in the current block. s - source, st - start.
    //! Get the end   local point index in the current block. s - source, ed - end.
    int GetStartOrignal(int st);
    int GetEndOrignal(int ed);

    // Calculate stOrignal and edOrignal.
    void ComputeOriginal();
    
    void SetFaceMatchingTargetDirIndex(int *dir3dIn);

    //! Set the current source zone index.
    void SetRegionBlock(int currentSourceZoneID);

    //! Set the target zone index.
    void SetTargetRegionBlock(int targetZoneID);
    
    //! Set the start/end range index in the current zone in i/j/k dimensions.
    void SetIJKRegion(int ist, int ied, int jst, int jed, int kst, int ked);

    //! Set the start/end range index in the neighbor zone in i/j/k dimensions.
    void SetTargetIJKRegion(int ist, int ied, int jst, int jed, int kst, int ked);

    //! Set this BC region on the left(-1) or right(1) in the current zone.
    void SetFaceLeftOrRightIndex(int s_lr);

    //! Set this BC region on the left(-1) or right(1) in the neighbor target zone.
    void SetFaceLeftOrRightIndexOfTargetBlock(int t_lr);

    //! Set the dimension number of this BC region on the current zone.
    //! 0-1-2, represents i, j, k direction, respectively.
    void SetFaceDirection(int s_nd);

    //! Set the dimension number of this BC region on the neighbor target zone.
    //! 0-1-2, represents i, j, k direction, respectively.
    void SetFaceDirectionOfTargetBlock(int t_nd);

    //! Return the current source zone index.
    int GetRegionBlock() const;

    //! Return the neighbor target zone index.
    int GetTargetRegionBlock() const;

    //! Return face direction.
    int * GetFaceDirectionIndex();

    //! Return direction between current face and matching target face.
    int * GetFaceMatchingTargetDirIndex();

    //! Return face direction in the neighbor zone.
    int * GetFaceDirectionIndexOfTargetBlock();

    //! Process to construct start/end local point index in the current/neighbor zones.
    void ProcessBCInfo();
    void SetLnkInfo();

    //! Init Face direction, all to zero.
    void InitFaceDirectionIndex();

    //! Return the start/end i/j/k range in the current zone.
    void GetIJKRegion(int &ist, int &ied, int &jst, int &jed, int &kst, int &ked);

    //! Return the three dimensional start index in the current zone.
    int * GetStartPoint();

    //! Return the m-th dimensional start index in the current zone.
    int GetStartPoint(int m);

    //! Return the three dimensional end index in the current zone.
    int * GetEndPoint(void);

    //! Return the m-th dimensional end index in the current zone.
    int GetEndPoint(int m);

    int * GetInkInfo();

    //! Return the m-th dimensional start index in the neighbor target zone.
    int GetTargetStart(int m);

    //! Return the m-th dimensional end index in the neighbor target zone.
    int GetTargetEnd(int m);
    
    //! Return the three dimensional start index in the neighbor target zone.
    int * GetTargetStart();

    //! Return the three dimensional end index in the neighbor target zone.
    int * GetTargetEnd();
    
    //! Return this BC region on the left(-1) or right(1) in the current zone.
    int GetFaceLeftOrRightIndex() const;

    //! Return this BC region on the left(-1) or right(1) in the neighbor target zone.
    int GetFaceLeftOrRightIndexOfTargetBlock() const;

    //! Return the dimension number of this BC region on the current zone.
    //! 0-1-2, represents i, j, k direction, respectively.
    int GetFaceDirection() const;

    //! Return the dimension number of this BC region on the neighbor target zone.
    //! 0-1-2, represents i, j, k direction, respectively.
    int GetFaceDirectionOfTargetBlock() const;

    //! Compute the coarse grid i/j/k range for multi-grid.
    void ComputeMultiGridIJKRegion(int imin, int imax, int jmin, int jmax, int kmin, int kmax, int istp, int jstp, int kstp, int nsurf);

    //! The followings used to accelerate the massive grid converting, or for overset grid,
    //! Do not care about it.
    void ComputeRelativeParameters();
    int * GetStam();
    int * GetRate();
    void ProbeExternalOverlapCells(int *iBlank, int nci, int ncj, int &counter);

    //! Return the i/j/k range in the current source zone, it is rarely used.
    void GetNormalizeIJKRegion(int &ist, int &ied, int &jst, int &jed, int &kst, int &ked);

    //! Return the cell index inside flowfield.
    void GetInsideCellIndex(int i, int j, int k, int &is, int &js, int &ks, int nss);

    //! Return the cell index of ghostcell.
    void GetGhostCellIndex(int i, int j, int k, int &it, int &jt, int &kt, int ntt);

    //! Return the index of boundary face.
    void GetBoundaryFaceIndex(int i, int j, int k, int &ib, int &jb, int &kb);

private:
    //! t_st: Start local point index in the neighbor zone. t - target, st - start.
    //! t_ed: End   local point index in the neighbor zone. t - target, ed - end.
    int t_st[3], t_ed[3];

    //! stam[3] && rate[3] are used to accelerate the massive grid converting.
    //! Do not care about the two.
    int stam[3], rate[3];
};

//! @brief StructBCSet class defines the composites of StructBC informations.
//! In one zone (block), several StructBCs may be exist, so all StructBCs in one zone
//! are defined as a 'StructBCSet'.
class StructBCSet
{
private:
    //! Number of zone (block) index this BC in.
    int zoneID;

    //! The BCRegions list in the current StructBCSet.
    vector<StructBC *> *bcRegions;


    //! The I/J/K id of interface faces.
    int *iID, *jID, *kID;

    //! The bcRegions id of interface faces.
    int *bcRegionIDofIFace; 

    //! Boundary name of given type;
    string boundaryName[100];

public:
    StructBCSet(int zoneID = 0);
    ~StructBCSet();

public:
    void SetzoneID(int ZONEID) { zoneID = ZONEID; }
    
    string GetBCName(int bcType) { return boundaryName[bcType]; }

    void CopyStructBCSet(StructBCSet *rightHandSide);

    StructBC * GetBCRegion(uint_t iBCRegion) const;

    void CreateBCRegion(uint_t nBCRegion);

    void SetBCRegion(uint_t iBCRegion, StructBC *bcregion);

    int GetnBCRegion() const;

    void GetSourceIndexIJK(int iFace, int ipos, int &i, int &j, int &k);
    void GetSourceIndexIJK_Nsurf_LR(int iFace, int ipos, int &i, int &j, int &k, int &nsruf, int &s_lr);
    void GetTargetIndexIJK(int iFace, int ipos, int &i, int &j, int &k);

    void GetSourceIndexIJK_fornodetransfer(int iFace, int &id, int &jd, int &kd, int &i_lr, int &j_lr, int &k_lr);
    void GetTargetIndexIJK_fornodetransfer(int iFace, int &id, int &jd, int &kd, int &i_lr, int &j_lr, int &k_lr);

    void ProcessBCInfo();
    void SetLnkInfo();

    void SetIFaceInfo();

    //! Get the bcregion id by interfaceIndexContainerForReceive value.
    int * GetIFaceInfo();

private:
    void DecodeIJK(int index, int &i, int &j, int &k);
    void DeleteIndex();
    void InitBoundaryName();
};

const int MERGED = 1;

#include "Geo_StructBC.hxx"

void GetBCFaceIDX(int *s_lr3d, int &id, int &jd, int &kd);

template < typename T >
void GhostCell2D(PHArray<T, 2> &w, int ni, int nj)
{
    //! w(-1:ni+1, -1:nj+1).

    for (int j = 1; j <= nj - 1; ++ j)
    {
        w(0 , j) = w(1   , j);
        w(ni, j) = w(ni-1, j);
    }

    for (int i = 0; i <= ni; ++ i)
    {
        w(i, 0) = w(i, 1  );
        w(i, nj) = w(i, nj-1);
    }

    //! Initialize volume in second row of ghost cells.
    for (int j = 0; j <= nj; ++ j)
    {
        w(-1  , j) = w(0 , j);
        w(ni+1, j) = w(ni, j);
    }

    for (int i = -1; i <= ni + 1; ++ i)
    {
        w(i, -1 ) = w(i, 0);
        w(i, nj+1) = w(i, nj);
    }
}

template < typename T >
void GhostCell2D(PHArray<T, 3> &w, int ni, int nj, int nm)
{
    //! w(-1:ni+1, -1:nj+1, 0:nm-1).

    for (int m = 0; m < nm; ++ m)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            w(0 , j, m) = w(1   , j, m);
            w(ni, j, m) = w(ni-1, j, m);
        }

        for (int i = 0; i <= ni; ++ i)
        {
            w(i, 0 , m) = w(i, 1   , m);
            w(i, nj, m) = w(i, nj-1, m);
        }

        //! Initialize volume in second row of ghost cells.
        for (int j = 0; j <= nj; ++ j)
        {
            w(-1  , j, m) = w(0 , j, m);
            w(ni+1, j, m) = w(ni, j, m);
        }

        for (int i = - 1; i <= ni + 1; ++ i)
        {
            w(i, -1  , m) = w(i, 0 , m);
            w(i, nj+1, m) = w(i, nj, m);
        }
    }
}

template < typename T >
void GhostCell3D(PHArray<T, 3> &w, int ni, int nj, int nk)
{
    //! w(-1:ni+1, -1:nj+1, -1:nk+1).
    int kst, ked;
    kst = 1;
    ked = nk;
    if (nk == 1) ked = 1;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            w(0 , j, k) = w(1   , j, k);
            w(ni, j, k) = w(ni-1, j, k);
        }

        for (int i = 0; i <= ni; ++ i)
        {
            w(i, 0 , k) = w(i, 1   , k);
            w(i, nj, k) = w(i, nj-1, k);
        }
        //! Initialize volume in second row of ghost cells.
        for (int j = 0; j <= nj; ++ j)
        {
            w(-1  , j, k) = w(0 , j, k);
            w(ni+1, j, k) = w(ni, j, k);
        }

        for (int i = -1; i <= ni + 1; ++ i)
        {
            w(i, -1  , k) = w(i, 0 , k);
            w(i, nj+1, k) = w(i, nj, k);
        }
    }

    if (nk == 1) return;

    for (int j = -1; j <= nj + 1; ++ j)
    {
        for (int i = -1; i <= ni + 1; ++ i)
        {
            w(i, j, 0) = w(i, j, 1);
            w(i, j, -1) = w(i, j, 1);

            w(i, j, nk ) = w(i, j, nk-1);
            w(i, j, nk+1) = w(i, j, nk-1);
        }
    }
}

template < typename T >
void GhostCell3D(PHArray<T, 4> &w, int ni, int nj, int nk, int nm)
{
    //! w(-1:ni+1, -1:nj+1, -1:nk+1, 0:nm-1).
    if (nk != 1)
    {
        for (int m = 0; m < nm; ++ m)
        {
            for (int k = 1; k <= nk - 1; ++ k)
            {
                for (int j = 1; j <= nj - 1; ++ j)
                {
                    w(0 , j, k, m) = w(1   , j, k, m);
                    w(ni, j, k, m) = w(ni-1, j, k, m);
                }

                for (int i = 0; i <= ni; ++ i)
                {
                    w(i, 0 , k, m) = w(i, 1   , k, m);
                    w(i, nj, k, m) = w(i, nj-1, k, m);
                }

                //! Initialize volume in second row of ghost cells.
                for (int j = 0; j <= nj; ++ j)
                {
                    w(-1  , j, k, m) = w(0 , j, k, m);
                    w(ni+1, j, k, m) = w(ni, j, k, m);
                }

                for (int i = -1; i <= ni + 1; ++ i)
                {
                    w(i, -1  , k, m) = w(i, 0 , k, m);
                    w(i, nj+1, k, m) = w(i, nj, k, m);
                }
            }

            for (int j = -1; j <= nj + 1; ++ j)
            {
                for (int i = -1; i <= ni + 1; ++ i)
                {
                    w(i, j, 0 , m) = w(i, j, 1, m);
                    w(i, j, -1, m) = w(i, j, 1, m);

                    w(i, j, nk  , m) = w(i, j, nk-1, m);
                    w(i, j, nk+1, m) = w(i, j, nk-1, m);
                }
            }
        }
    }
    else
    {
        int k = 1;
        for (int m = 0; m < nm; ++ m)
        {
            for (int j = 1; j <= nj - 1; ++ j)
            {
                w(0 , j, k, m) = w(1   , j, k, m);
                w(ni, j, k, m) = w(ni-1, j, k, m);
            }

            for (int i = 0; i <= ni; ++ i)
            {
                w(i, 0 , k, m) = w(i, 1   , k, m);
                w(i, nj, k, m) = w(i, nj-1, k, m);
            }

            //! Initialize volume in second row of ghost cells.
            for (int j = 0; j <= nj; ++ j)
            {
                w(-1  , j, k, m) = w(0 , j, k, m);
                w(ni+1, j, k, m) = w(ni, j, k, m);
            }

            for (int i = -1; i <= ni + 1; ++ i)
            {
                w(i, -1  , k, m) = w(i, 0 , k, m);
                w(i, nj+1, k, m) = w(i, nj, k, m);
            }
        }
    }
}

template < typename T >
void FillCornerPoint3D(PHArray<T, 3> &w, int ni, int nj, int nk)
{
    //! This subroutine is developed to fill the ghost points which are 
    //! shared by two planes. this is necessary for the calculation of the 
    //! viscous fluxes.

    //! The figure below shows a 2d plane.

    //!      nsurf -->

    //!         $(a)  |   x   |                x  = internal cell
    //!               |       |
    //!      .........-----------------        $  = ghost cell filled by
    //!               .       .                     bcin or bcvis
    //!         ?     .   $(b).                ?  = ghost cell filled by 
    //!      ..........................             this routine

    //! The algorithm used to fill the points ? is the following. in the
    //! algorithm for the calculation of the viscous fluxes, the average
    //! is taken from the 4 points shown in the figure. this average is
    //! replaced by taking the average of the points $(a) and x, hence
    //! the value of ? is set to accomplish this.
    //! 3d: fill the 8 vertex lines for each direction nsurf.

    if (nk != 1)
    {
        for (int k = 1; k <= nk - 1; ++ k)
        {
            w(0,  0, k) = w(0 ,    1, k) + w(  1,  0, k) - w(  1,    1, k);
            w(0, nj, k) = w(0 , nj-1, k) + w(  1, nj, k) - w(  1, nj-1, k);
            w(ni,  0, k) = w(ni,    1, k) + w(ni-1,  0, k) - w(ni-1,    1, k);
            w(ni, nj, k) = w(ni, nj-1, k) + w(ni-1, nj, k) - w(ni-1, nj-1, k);
        }

        for (int j = 1; j <= nj - 1; ++ j)
        {
            w(0, j,  0) = w(0, j,    1) + w(  1, j,  0) - w(  1, j,    1);
            w(0, j, nk) = w(0, j, nk-1) + w(  1, j, nk) - w(  1, j, nk-1);
            w(ni, j,  0) = w(ni, j,    1) + w(ni-1, j,  0) - w(ni-1, j,    1);
            w(ni, j, nk) = w(ni, j, nk-1) + w(ni-1, j, nk) - w(ni-1, j, nk-1);
        }

        for (int i = 1; i <= ni - 1; ++ i)
        {
            w(i,  0,  0) = w(i,  0,    1) + w(i,    1,  0) - w(i,    1,    1);
            w(i,  0, nk) = w(i,  0, nk-1) + w(i,    1, nk) - w(i,    1, nk-1);
            w(i, nj,  0) = w(i, nj,    1) + w(i, nj-1,  0) - w(i, nj-1,    1);
            w(i, nj, nk) = w(i, nj, nk-1) + w(i, nj-1, nk) - w(i, nj-1, nk-1);
        }

        //! Fill the 8 corner points.
        w(0,  0,  0) = third * (w(  1,  0,  0) + w(0,    1,  0) + w(0,  0,    1));
        w(0, nj,  0) = third * (w(  1, nj,  0) + w(0, nj-1,  0) + w(0, nj,    1));
        w(0,  0, nk) = third * (w(  1,  0, nk) + w(0,    1, nk) + w(0,  0, nk-1));
        w(0, nj, nk) = third * (w(  1, nj, nk) + w(0, nj-1, nk) + w(0, nj, nk-1));
        w(ni,  0,  0) = third * (w(ni-1,  0,  0) + w(ni,    1,  0) + w(ni,  0,    1));
        w(ni, nj,  0) = third * (w(ni-1, nj,  0) + w(ni, nj-1,  0) + w(ni, nj,    1));
        w(ni,  0, nk) = third * (w(ni-1,  0, nk) + w(ni,    1, nk) + w(ni,  0, nk-1));
        w(ni, nj, nk) = third * (w(ni-1, nj, nk) + w(ni, nj-1, nk) + w(ni, nj, nk-1));
    }
    else
    {
        int k = 1;
        w(0,  0, k) = w(0 ,    1, k) + w(  1,  0, k) - w(  1,    1, k);
        w(0, nj, k) = w(0 , nj-1, k) + w(  1, nj, k) - w(  1, nj-1, k);
        w(ni,  0, k) = w(ni,    1, k) + w(ni-1,  0, k) - w(ni-1,    1, k);
        w(ni, nj, k) = w(ni, nj-1, k) + w(ni-1, nj, k) - w(ni-1, nj-1, k);
    }
}

template < typename T >
void FillCornerPoint3D(PHArray<T, 4 > &w, int ni, int nj, int nk, int nl)
{
    //! This subroutine is developed to fill the ghost points which are 
    //! shared by two planes. this is necessary for the calculation of the 
    //! viscous fluxes.

    //! The figure below shows a 2d plane.

    //!      nsurf -->

    //!         $(a)  |   x   |                x  = internal cell
    //!               |       |
    //!      .........-----------------        $  = ghost cell filled by
    //!               .       .                     bcin or bcvis
    //!         ?     .   $(b).                ?  = ghost cell filled by 
    //!      ..........................             this routine

    //! The algorithm used to fill the points ? is the following. in the
    //! algorithm for the calculation of the viscous fluxes, the average
    //! is taken from the 4 points shown in the figure. this average is
    //! replaced by taking the average of the points $(a) and x, hence
    //! the value of ? is set to accomplish this.
    //! 3d: fill the 8 vertex lines for each direction nsurf.

    if (nk != 1)
    {
        for (int m = 0; m < nl; ++ m)
        {
            for (int k = 1; k <= nk - 1; ++ k)
            {
                w(0,  0, k, m) = w(0 ,    1, k, m) + w(  1,  0, k, m) - w(  1,    1, k, m);
                w(0, nj, k, m) = w(0 , nj-1, k, m) + w(  1, nj, k, m) - w(  1, nj-1, k, m);
                w(ni,  0, k, m) = w(ni,    1, k, m) + w(ni-1,  0, k, m) - w(ni-1,    1, k, m);
                w(ni, nj, k, m) = w(ni, nj-1, k, m) + w(ni-1, nj, k, m) - w(ni-1, nj-1, k, m);
            }

            for (int j = 1; j <= nj - 1; ++ j)
            {
                w(0, j,  0, m) = w(0, j,    1, m) + w(  1, j,  0, m) - w(  1, j,    1, m);
                w(0, j, nk, m) = w(0, j, nk-1, m) + w(  1, j, nk, m) - w(  1, j, nk-1, m);
                w(ni, j,  0, m) = w(ni, j,    1, m) + w(ni-1, j,  0, m) - w(ni-1, j,    1, m);
                w(ni, j, nk, m) = w(ni, j, nk-1, m) + w(ni-1, j, nk, m) - w(ni-1, j, nk-1, m);
            }

            for (int i = 1; i <= ni - 1; ++ i)
            {
                w(i,  0,  0, m) = w(i,  0,    1, m) + w(i,    1,  0, m) - w(i,    1,    1, m);
                w(i,  0, nk, m) = w(i,  0, nk-1, m) + w(i,    1, nk, m) - w(i,    1, nk-1, m);
                w(i, nj,  0, m) = w(i, nj,    1, m) + w(i, nj-1,  0, m) - w(i, nj-1,    1, m);
                w(i, nj, nk, m) = w(i, nj, nk-1, m) + w(i, nj-1, nk, m) - w(i, nj-1, nk-1, m);
            }
        }

        //! Fill the 8 corner points.
        for (int m = 0; m < nl; ++ m)
        {
            w(0,  0,  0, m) = third * (w(  1,  0,  0, m) + w(0,    1,  0, m) + w(0,  0,    1, m));
            w(0, nj,  0, m) = third * (w(  1, nj,  0, m) + w(0, nj-1,  0, m) + w(0, nj,    1, m));
            w(0,  0, nk, m) = third * (w(  1,  0, nk, m) + w(0,    1, nk, m) + w(0,  0, nk-1, m));
            w(0, nj, nk, m) = third * (w(  1, nj, nk, m) + w(0, nj-1, nk, m) + w(0, nj, nk-1, m));
            w(ni,  0,  0, m) = third * (w(ni-1,  0,  0, m) + w(ni,    1,  0, m) + w(ni,  0,    1, m));
            w(ni, nj,  0, m) = third * (w(ni-1, nj,  0, m) + w(ni, nj-1,  0, m) + w(ni, nj,    1, m));
            w(ni,  0, nk, m) = third * (w(ni-1,  0, nk, m) + w(ni,    1, nk, m) + w(ni,  0, nk-1, m));
            w(ni, nj, nk, m) = third * (w(ni-1, nj, nk, m) + w(ni, nj-1, nk, m) + w(ni, nj, nk-1, m));
        }
    }
    else
    {
        int k = 1;
        for (int m = 0; m < nl; ++ m)
        {
            w(0,  0, k, m) = w(0 ,    1, k, m) + w(  1,  0, k, m) - w(  1,    1, k, m);
            w(0, nj, k, m) = w(0 , nj-1, k, m) + w(  1, nj, k, m) - w(  1, nj-1, k, m);
            w(ni,  0, k, m) = w(ni,    1, k, m) + w(ni-1,  0, k, m) - w(ni-1,    1, k, m);
            w(ni, nj, k, m) = w(ni, nj-1, k, m) + w(ni-1, nj, k, m) - w(ni-1, nj-1, k, m);
        }
    }
}

}