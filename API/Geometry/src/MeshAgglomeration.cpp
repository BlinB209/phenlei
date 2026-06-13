#include <list>
#include "MeshAgglomeration.h"
#include "GridType.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "Geo_FaceMetrics_Unstruct.h"
#include "Geo_CellMetrics_Unstruct.h"
#include "Geo_DynamicGridMetrics_Unstruct.h"
#include "Geo_LSQWeight_Unstruct.h"
#include "Geo_UnstructBC.h"
#include "Geo_Face.h"
#include "GRHash.h"
#include "Geo_FaceTopo_Unstruct.h"
#include "mgridgen.h"
#include "Glb_Dimension.h"

namespace PHSPACE{
LIB_EXPORT Mesh_Agglomeration::Mesh_Agglomeration(UnstructGrid * fineGridIn, int fineGridLevelIn)
{
    this->fineGrid = fineGridIn;
    this->fineGridLevel = fineGridLevelIn;

    this->coarseGrid = 0;
}

LIB_EXPORT Mesh_Agglomeration::~Mesh_Agglomeration()
{
}

LIB_EXPORT void Mesh_Agglomeration::CoarsenGridOnce()
{
    int nTotalFace = fineGrid->GetNTotalFace();
    int nBoundFace = fineGrid->GetNBoundFace();
    int nTotalCell = fineGrid->GetNTotalCell();
    int *left_cell_of_face = fineGrid->GetLeftCellOfFace();
    int *right_cell_of_face = fineGrid->GetRightCellOfFace();

    RDouble * vol = fineGrid->GetCellVolume();
    RDouble * ns  = fineGrid->GetFaceArea();

    RDouble * xcc = fineGrid->GetCellCenterX();
    RDouble * ycc = fineGrid->GetCellCenterY();
    RDouble * zcc = fineGrid->GetCellCenterZ();

    RDouble * virtualCellCenterX = new RDouble [nTotalFace];
    RDouble * virtualCellCenterY = new RDouble [nTotalFace];
    RDouble * virtualCellCenterZ = new RDouble [nTotalFace];
    PHSPACE::SetField(virtualCellCenterX, 0.0, nTotalFace);
    PHSPACE::SetField(virtualCellCenterY, 0.0, nTotalFace);
    PHSPACE::SetField(virtualCellCenterZ, 0.0, nTotalFace);


    RDouble *tarea = new RDouble[nTotalCell];   // total area of each cell
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        tarea[iCell] = 0.0;
    }

    int le = 0, re = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = left_cell_of_face[iFace];
        tarea[le] += ns[iFace];
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = left_cell_of_face [iFace];
        re = right_cell_of_face[iFace];
        tarea[le] += ns[iFace];
        tarea[re] += ns[iFace];
    }

    int *fmark = new int[nTotalFace];

    SetField(fmark, -1, nTotalFace);
    int cnTFace = 0;

    //! Prepare the inner-faces that have same left&&right cell, on coarse grid.
    //! fmark: mapping face on fine grid to face on coarse grid.
    //! fmark = -1: the two neighbor cells would not been agglomerated.
    using namespace GRIDER_SPACE;
    int nInteriorFace = nTotalFace - nBoundFace;
    FaceHash faceHash;
    vector<Geo_Face *> faceList;
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = left_cell_of_face [iFace];
        re = right_cell_of_face[iFace];

        Geo_Face *face = new Geo_Face(le, re, iFace, nInteriorFace);

        bool isFaceExist;

        // Insert the node into the tetNode, which is a hash table.
        int faceIndexinHash = faceHash.insert(face, isFaceExist) - 1;

        // Push the node into node list of GRAdvanceFront.
        if (isFaceExist)
        {
            // the coarse face have been built.
            delete face;
            face = faceHash[faceIndexinHash];
            fmark[iFace] = face->GetID();
        }
        else
        {
            // build new coarse face.
            face->SetID(cnTFace);
            fmark[iFace] = cnTFace;
            cnTFace ++;

            faceList.push_back(face);
        }
    }
    faceHash.clear();
    vector<Geo_Face *>::iterator faceIter;
    for (faceIter = faceList.begin(); faceIter != faceList.end(); ++ faceIter)
    {
        delete (* faceIter);
    }
    faceList.clear();


    // barea: total area between two cells, may be more than one face between two cells on coarse grid.
    RDouble *barea = new RDouble[ cnTFace ];
    for (int i = 0; i < cnTFace; ++ i)
    {
        barea[i] = 0.0;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        barea[fmark[iFace]] += ns[iFace];
    }

    RDouble *ratio = new RDouble[nTotalFace];
    int fno;
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = left_cell_of_face [iFace];
        re = right_cell_of_face[iFace];
        fno = fmark[iFace];

        ratio[iFace] = pow((tarea[le] + tarea[re] - 2.0 * barea[fno]), 1.5) / (vol[le] + vol[re]);

        virtualCellCenterX[iFace] = (vol[le] * xcc[le] + vol[re] * xcc[re]) / (vol[le] + vol[re]);
        virtualCellCenterY[iFace] = (vol[le] * ycc[le] + vol[re] * ycc[re]) / (vol[le] + vol[re]);
        virtualCellCenterZ[iFace] = (vol[le] * zcc[le] + vol[re] * zcc[re]) / (vol[le] + vol[re]);
    }


    delete [] barea;
    delete [] tarea;

    RDouble limitAngle = -90.0;
    AspectRatioModifyBySkewness(ratio, virtualCellCenterX, virtualCellCenterY, virtualCellCenterZ, limitAngle);

    typedef pair<RDouble,int> TRatioIndex;
    typedef set<TRatioIndex,less<TRatioIndex> > Tset;

    Tset myset;

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        myset.insert(TRatioIndex(ratio[iFace], iFace));
    }

    Tset::iterator curr = myset.begin();

    typedef list<int> Tlist;
    Tlist mylist;

    bool * unMerged = new bool [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        unMerged[iCell] = true;
    }

    for (curr = myset.begin(); curr != myset.end(); ++ curr)
    {
        fno = (*curr).second;
        le = left_cell_of_face[ fno ];
        re = right_cell_of_face[ fno ];
        if (unMerged[le] && unMerged[re] && ratio[fno] >= 0.0)
        {
            mylist.push_back(fno);
            unMerged[le] = false;
            unMerged[re] = false;
        }
    }

    delete [] virtualCellCenterX;
    delete [] virtualCellCenterY;
    delete [] virtualCellCenterZ;

    // now merging
    int *cell2coarsegridcell = new int[ nTotalCell ];
    fineGrid->SetCell2CoarseGridCell(cell2coarsegridcell);
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2coarsegridcell[iCell] = - 1;
    }

    int cnTCell = 0;
    Tlist::iterator lcurr;
    for (lcurr = mylist.begin(); lcurr != mylist.end(); ++ lcurr)
    {
        fno =  * lcurr;
        le = left_cell_of_face[ fno ];
        re = right_cell_of_face[ fno ];
        cell2coarsegridcell[le] = cnTCell;
        cell2coarsegridcell[re] = cnTCell ++;
    }

    // some cells may not have been merged, merge them now
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (cell2coarsegridcell[iCell] == -1)
        {
            cell2coarsegridcell[iCell] = cnTCell ++;
        }
    }

    delete [] ratio;
    delete [] unMerged;
    delete [] fmark;

    coarseGrid = new UnstructGrid();
    coarseGrid->InitGrid(new GridID(*fineGrid->GetGridID()), fineGridLevel + 1, fineGrid->GetDim(), UNSTRUCTGRID);

    coarseGrid->SetNTotalCell(cnTCell);
    int nTotalNode = fineGrid->GetNTotalNode();
    coarseGrid->SetNTotalNode(nTotalNode);
    RDouble *x, *y, *z;
    RDouble * xx = fineGrid->GetX();
    RDouble * yy = fineGrid->GetY();
    RDouble * zz = fineGrid->GetZ();

    x = new RDouble[ nTotalNode ];
    y = new RDouble[ nTotalNode ];
    z = new RDouble[ nTotalNode ];

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        x[iNode] = xx[iNode];
        y[iNode] = yy[iNode];
        z[iNode] = zz[iNode];
    }
    coarseGrid->SetX(x);
    coarseGrid->SetY(y);
    coarseGrid->SetZ(z);
}

LIB_EXPORT void Mesh_Agglomeration::CoarsenGridOncebyMGridgen()
{
    int nTotalFace = fineGrid->GetNTotalFace();
    int nBoundFace = fineGrid->GetNBoundFace();
    int nTotalCell = fineGrid->GetNTotalCell();
    int *left_cell_of_face = fineGrid->GetLeftCellOfFace();

    RDouble * vol = fineGrid->GetCellVolume();
    RDouble * area= fineGrid->GetFaceArea();

    //! bcArea: BC face area of each cell.
    RDouble *bcArea = new RDouble [nTotalCell];
    PHSPACE::SetField(bcArea, 0.0, nTotalCell);

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = left_cell_of_face[iFace];
        bcArea[le] += area[iFace];
    }

    // now merging
    const int NON_COARSED = -1;
    int *cell2coarsegridcell = new int[ nTotalCell ];
    fineGrid->SetCell2CoarseGridCell(cell2coarsegridcell);
    PHSPACE::SetField(cell2coarsegridcell, NON_COARSED, nTotalCell);
    
    //! From the view of graph theory, each face connection is defined by two neighbor cells.
    //! so, number of connection equal to the two times of interior faces.
    int nConnection = 2 * (nTotalFace - nBoundFace);

    //! Adjacent connection of cells, using CSR format which is same to METIS. 
    idx_t *xadj   = new idx_t[nTotalCell + 1];
    idx_t *adjncy = new idx_t[nConnection];
    Get_Xadj_Adjncy(fineGrid, xadj, adjncy, 0, 0);

    RDouble *adjncywgt = new RDouble [nConnection];
    Get_Xadj_Adjncy_Weight(fineGrid, xadj, adjncy, adjncywgt);

    //! Convert the metis type to mgrid type.
    idxtype *xadjNew   = new idxtype[nTotalCell + 1];
    idxtype *adjncyNew = new idxtype[nConnection];
    for (int iCell = 0; iCell < nTotalCell + 1; ++ iCell)
    {
        xadjNew[iCell] = static_cast<idxtype>(xadj[iCell]);
    }
    for (int iFace = 0; iFace < nConnection; ++ iFace)
    {
        adjncyNew[iFace] = static_cast<idxtype>(adjncy[iFace]);
    }
    delete [] xadj;
    delete [] adjncy;

    int nDim = PHSPACE::GetDim();
    int minCoarseCellNumber = 1;
    int maxCoarseCellNumber = 6;
    if (nDim == PHSPACE::THREE_D)
    {
        minCoarseCellNumber = 4;
        maxCoarseCellNumber = 10;
    }

    int options[4] = {4, 6, 0, nDim};
    int nmove, nparts;
    MGridGen(nTotalCell, xadjNew, vol, bcArea, adjncyNew, adjncywgt, minCoarseCellNumber, maxCoarseCellNumber,
        options, &nmove, &nparts, cell2coarsegridcell);

    //! For better grid quality in parallel,
    //! ParMGridgen should be used in future, then, the parMetis should be followed.
    //! For now, the serial version of mgridgen is used even in parallel.

    delete [] xadjNew;
    delete [] adjncyNew;   
    delete [] adjncywgt; 
    delete [] bcArea;

    //! Now, construct the basic geometry information of coarse grid.
    int cnTCell = nparts;

    coarseGrid = new UnstructGrid();
    coarseGrid->InitGrid(new GridID(*fineGrid->GetGridID()), fineGridLevel + 1, fineGrid->GetDim(), UNSTRUCTGRID);

    coarseGrid->SetNTotalCell(cnTCell);
    int nTotalNode = fineGrid->GetNTotalNode();
    coarseGrid->SetNTotalNode(nTotalNode);
    RDouble *x, *y, *z;
    RDouble * xx = fineGrid->GetX();
    RDouble * yy = fineGrid->GetY();
    RDouble * zz = fineGrid->GetZ();

    x = new RDouble[ nTotalNode ];
    y = new RDouble[ nTotalNode ];
    z = new RDouble[ nTotalNode ];

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        x[iNode] = xx[iNode];
        y[iNode] = yy[iNode];
        z[iNode] = zz[iNode];
    }
    coarseGrid->SetX(x);
    coarseGrid->SetY(y);
    coarseGrid->SetZ(z);
}


//global agglomeration method used to coarsen the grid to generate cell2coarsegridcell only
LIB_EXPORT void Mesh_Agglomeration::CoarsenGridOnceWithoutSkewnessLimit()
{
    int nTotalFace = fineGrid->GetNTotalFace();
    int nBoundFace = fineGrid->GetNBoundFace();
    int nTotalCell = fineGrid->GetNTotalCell();
    int *left_cell_of_face = fineGrid->GetLeftCellOfFace();
    int *right_cell_of_face = fineGrid->GetRightCellOfFace();

    RDouble * vol = fineGrid->GetCellVolume();
    RDouble * ns  = fineGrid->GetFaceArea();

    int le = 0, re = 0;

    RDouble *tarea = new RDouble[nTotalCell];   // total area of each cell
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        tarea[iCell] = 0.0;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = left_cell_of_face[iFace];
        tarea[le] += ns[iFace];
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = left_cell_of_face [iFace];
        re = right_cell_of_face[iFace];
        tarea[le] += ns[iFace];
        tarea[re] += ns[iFace];
    }

    int *fmark = new int[nTotalFace];
    SetField(fmark, -1, nTotalFace);
    int cnTFace = 0;

    //! Prepare the inner-faces that have same left&&right cell, on coarse grid.
    //! fmark: mapping face on fine grid to face on coarse grid.
    //! fmark = -1: the two neighbor cells would not been agglomerated.
    using namespace GRIDER_SPACE;
    int nInteriorFace = nTotalFace - nBoundFace;
    FaceHash faceHash;
    vector<Geo_Face *> faceList;
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = left_cell_of_face [iFace];
        re = right_cell_of_face[iFace];

        Geo_Face *face = new Geo_Face(le, re, iFace, nInteriorFace);

        bool isFaceExist;

        // Insert the node into the tetNode, which is a hash table.
        int faceIndexinHash = faceHash.insert(face, isFaceExist) - 1;

        // Push the node into node list of GRAdvanceFront.
        if (isFaceExist)
        {
            // the coarse face have been built.
            delete face;
            face = faceHash[faceIndexinHash];
            fmark[iFace] = face->GetID();
        }
        else
        {
            // build new coarse face.
            face->SetID(cnTFace);
            fmark[iFace] = cnTFace;
            cnTFace ++;
        }
    }
    faceHash.clear();
    vector<Geo_Face *>::iterator faceIter;
    for (faceIter = faceList.begin(); faceIter != faceList.end(); ++ faceIter)
    {
        delete (*faceIter);
    }
    faceList.clear();

    // barea: total area between two cells, may be more than one face between two cells on coarse grid.
    RDouble *barea = new RDouble[ cnTFace ];
    for (int i = 0; i < cnTFace; ++ i)
    {
        barea[i] = 0.0;
    }
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        barea[fmark[iFace]] += ns[iFace];
    }

    RDouble * ratio = new RDouble[nTotalFace];
    int fno;
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = left_cell_of_face [iFace];
        re = right_cell_of_face[iFace];
        fno = fmark[iFace];
        ratio[iFace] = pow((tarea[le] + tarea[re] - 2.0 * barea[fno]), 1.5) / (vol[le] + vol[re]);
    }
    delete [] barea;
    delete [] tarea;

    typedef pair<RDouble,int> TRatioIndex;
    typedef set<TRatioIndex,less<TRatioIndex> > Tset;
    Tset myset;
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        myset.insert(TRatioIndex(ratio[iFace], iFace));
    }

    Tset::iterator curr = myset.begin();
    typedef list<int> Tlist;
    Tlist mylist;
    bool * unMerged = new bool [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        unMerged[iCell] = true;
    }
    for (curr = myset.begin(); curr != myset.end(); ++ curr)
    {
        fno = (*curr).second;
        le = left_cell_of_face[ fno ];
        re = right_cell_of_face[ fno ];
        if (unMerged[le] && unMerged[re])
        {
            mylist.push_back(fno);
            unMerged[le] = false;
            unMerged[re] = false;
        }
    }

    // now merging
    int *cell2coarsegridcell = new int[ nTotalCell ];
    fineGrid->SetCell2CoarseGridCell(cell2coarsegridcell);
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2coarsegridcell[iCell] = - 1;
    }

    int cnTCell = 0;
    Tlist::iterator lcurr;
    for (lcurr = mylist.begin(); lcurr != mylist.end(); ++ lcurr)
    {
        fno =  * lcurr;
        le = left_cell_of_face[ fno ];
        re = right_cell_of_face[ fno ];
        cell2coarsegridcell[le] = cnTCell;
        cell2coarsegridcell[re] = cnTCell ++;
    }

    // some cells may not have been merged, merge them now
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (cell2coarsegridcell[iCell] == -1)
        {
            cell2coarsegridcell[iCell] = cnTCell ++;
        }
    }

    delete [] ratio;
    delete [] unMerged;
    delete [] fmark;

    coarseGrid = new UnstructGrid();
    coarseGrid->InitGrid(new GridID(*fineGrid->GetGridID()), fineGridLevel + 1, fineGrid->GetDim(), UNSTRUCTGRID);

    coarseGrid->SetNTotalCell(cnTCell);
    int nTotalNode = fineGrid->GetNTotalNode();
    coarseGrid->SetNTotalNode(nTotalNode);
    RDouble *x, *y, *z;
    RDouble * xx = fineGrid->GetX();
    RDouble * yy = fineGrid->GetY();
    RDouble * zz = fineGrid->GetZ();

    x = new RDouble[ nTotalNode ];
    y = new RDouble[ nTotalNode ];
    z = new RDouble[ nTotalNode ];

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        x[iNode] = xx[iNode];
        y[iNode] = yy[iNode];
        z[iNode] = zz[iNode];
    }
    coarseGrid->SetX(x);
    coarseGrid->SetY(y);
    coarseGrid->SetZ(z);
}


void Mesh_Agglomeration::AspectRatioModifyBySkewness(RDouble *ratio, RDouble * xcc, RDouble * ycc, RDouble * zcc, const RDouble & limitAngle)
{
    int nTotalFace = fineGrid->GetNTotalFace();
    int nBoundFace = fineGrid->GetNBoundFace();
    int *left_cell_of_face = fineGrid->GetLeftCellOfFace();
    int *right_cell_of_face = fineGrid->GetRightCellOfFace();

    RDouble * xfc = fineGrid->GetFaceCenterX();
    RDouble * yfc = fineGrid->GetFaceCenterY();
    RDouble * zfc = fineGrid->GetFaceCenterZ();

    RDouble * nxs  = fineGrid->GetFaceNormalX();
    RDouble * nys  = fineGrid->GetFaceNormalY();
    RDouble * nzs  = fineGrid->GetFaceNormalZ();

    int * face_number_of_each_cell = fineGrid->GetFaceNumberOfEachCell();
    int ** cell2face = fineGrid->GetCell2Face();

    int LRCell[2];
    RDouble dx, dy, dz, dist, dot;
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        LRCell[0] = left_cell_of_face[iFace];
        LRCell[1] = right_cell_of_face[iFace];


        RDouble minAngle = 100.0;
        for (int iCell = 0; iCell < 2; ++ iCell)
        {
            int cellID = LRCell[iCell];

            int nFace = face_number_of_each_cell[cellID];
            for (int jFace = 0; jFace < nFace; ++ jFace)
            {
                int faceID = cell2face[cellID][jFace];
                if (faceID == iFace) continue;

                dx = xfc[faceID] - xcc[iFace];
                dy = yfc[faceID] - ycc[iFace];
                dz = zfc[faceID] - zcc[iFace];
                dist  = DISTANCE(dx, dy, dz);
                dot = (nxs[faceID] * dx + nys[faceID] * dy + nzs[faceID] * dz) / (dist + TINY);
                if (dot >  1.0) dot =   1.0;
                if (dot < -1.0) dot = - 1.0;

                if (right_cell_of_face[faceID] == cellID)
                {
                    dot *= -1;
                }

                RDouble angle  = asin(dot) * 180.0 / PI;
                minAngle = MIN(angle, minAngle);
            }
        }

        RDouble factor = minAngle / 90.0;

        factor = 1.0 - factor;

        ratio[iFace] *= factor;

        if (minAngle < limitAngle) ratio[iFace] = -1;
    }
}

}