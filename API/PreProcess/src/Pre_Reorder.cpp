#include "Pre_Reorder.h"
namespace Reorder
{
    void CallBandwidthFTS(const int nTotalCell, vector<int> *cell2Cell, const int *nNeighborCellOfEachCell)
    {
        printf("Program is running in CallBandwidthEvaluate\n");
        int *bandwidthOfEachCell = new int[nTotalCell];
        for (int iCell = 0; iCell < nTotalCell; iCell ++)
        {
            int numCells = nNeighborCellOfEachCell[iCell];
            int localMaxBand = 0;
            for (int ngbCell = 0; ngbCell < numCells; ++ ngbCell)
            {
                int bandDiff = abs(iCell - cell2Cell[iCell][ngbCell]);
                if (bandDiff > localMaxBand) localMaxBand = bandDiff;
            }
            bandwidthOfEachCell[iCell] = localMaxBand;
        }

        //statistics of bandwidthOfEachCell
        int maxGlobalBand = 0; //maximum bandwidth
        int maxCellID = 0; //cellID of maximum bandwidth
        int bandDist[33] = {0}; //distribution of bandwidth from 1 to 32 and more
        for (int i = 0; i < 33; ++ i)
        {
            bandDist[i] = 0;
        }

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            int bandCell = bandwidthOfEachCell[iCell];
            if (bandCell <= 0) 
            {
                printf("Error: bandwidth<=0 on cell %d, bandwidth=%d\n", iCell, bandCell);
                exit(1);
            }
            if (bandCell > maxGlobalBand)
            {
                maxGlobalBand = bandCell;
                maxCellID = iCell;
            }
            if (bandCell > 32) bandDist[32]++;
            else bandDist[bandCell-1]++;
        }
        
        printf("maxGlobalBand = %d on cellID %d\n", maxGlobalBand, maxCellID);
        int nTotalDist = 0;
        for (int i = 0; i < 33; ++ i)
        {
            nTotalDist += bandDist[i];
            printf("bandDist %d: %d in nTotalCells\n", i+1, bandDist[i]);
        }

        printf("nTotalDist = %d, nTotalCell = %d\n", nTotalDist, nTotalCell);

    }

    void CallBandwidthEvaluate(Base_Grid_Conn *gconn)
    {
        printf("Program is running in CallBandwidthEvaluate\n");
        cgsize_t nTotalFace = gconn->GetNTotalFace();
        cgsize_t nBoundFace = gconn->GetNBoundFace();
        cgsize_t nTotalCell = gconn->GetNTotalCell();
        printf("nTotalFace = %d, nBoundFace = %d, nTotalCell = %d\n", static_cast<int>(nTotalFace), static_cast<int>(nBoundFace), static_cast<int>(nTotalCell));
        //How to get leftCellOfFace and rightCellOfFace is from FYGridBase.cpp
        vector<cgsize_t> &leftCellOfFace  = *gconn->GetLeftCellOfFace();
        vector<cgsize_t> &rightCellOfFace = *gconn->GetRightCellOfFace();
        //It is difficult to create cell2Face, because the number of faces belonging to a face is not determined. So, I decide to use the familiar schemes
        int *cell2Face;
        int *cell2Cell;

        //offset of iCell's cell faces in cell2Face
        int *offsetCell2Face = new int[nTotalCell];
        int nTotalCellFaces = 0;

        int *nFaceOfEachCell = new int[nTotalCell];
        int *bandwidthOfEachCell = new int[nTotalCell];

        //Initialization of nFaceOfEachCell
        printf("Initialization\n");
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            nFaceOfEachCell[iCell] = 0;
            bandwidthOfEachCell[iCell] = 0;
            offsetCell2Face[iCell] = 0;
        }
        //set nFaceOfEachCell and nTotalCellFaces
        printf("set nFaceOfEachCell and nTotalCellFaces\n");
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            //for leftCellOfFace
            int cellLabelLeft = static_cast<int>(leftCellOfFace[iFace]);
            if ((cellLabelLeft < 0) || (cellLabelLeft >= nTotalCell))
            {
                printf("Error: cellLabelLeft = %d\n", cellLabelLeft);
                exit(1);
            }
            //nFaceOfEachCell[cellLabel]++;
            //nTotalCellFaces++;
            //for rightCellOfFace
            int cellLabelRight = static_cast<int>(rightCellOfFace[iFace]);
            if (cellLabelRight == -1) continue; //boundy face's right neighbor

            if ((cellLabelRight < 0) || (cellLabelRight >= nTotalCell)) 
            {
                printf("Error: cellLabelRight = %d\n", cellLabelRight);
                exit(1);
            }
            nFaceOfEachCell[cellLabelLeft]++;
            nTotalCellFaces++;
            nFaceOfEachCell[cellLabelRight]++;
            nTotalCellFaces++;
        }
        printf("set the size of cell2Cell and cell2Face\n");
        //test begin
            //printf("On cell %d, nFaceOfEachCell[%d]=%d\n", 0, 0, nFaceOfEachCell[0]);
        //test end
        //set the size of cell2Cell and cell2Face
        cell2Face = new int[nTotalCellFaces];
        cell2Cell = new int[nTotalCellFaces];
        //set offsetCell2Face
        offsetCell2Face[0] = 0;
        for (int iCell = 1; iCell < nTotalCell; ++ iCell)
        {
            offsetCell2Face[iCell] = offsetCell2Face[iCell-1] + nFaceOfEachCell[iCell-1];
        }
        //set cell2Face
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            //for leftCellOfFace
            int cellLabelLeft = static_cast<int>(leftCellOfFace[iFace]);
            if ((cellLabelLeft < 0) || (cellLabelLeft >= nTotalCell))
            {
                printf("Error: cellLabelLeft = %d\n", cellLabelLeft);
                exit(1);
            }
            //int cell2faceOffset = offsetCell2Face[cellLabel];
            //cell2Face[cell2faceOffset] = iFace;
            //offsetCell2Face update for the next set of face
            //offsetCell2Face[cellLabel]++;
            //for rightCellOfFace
            int cellLabelRight = static_cast<int>(rightCellOfFace[iFace]);
            if (cellLabelRight == -1) continue; //boundy face's right neighbor

            if ((cellLabelRight < 0) || (cellLabelRight >= nTotalCell))
            {
                printf("Error: cellLabelRight = %d\n", cellLabelRight);
                exit(1);
            }
            int cell2faceOffset = offsetCell2Face[cellLabelLeft];
            cell2Face[cell2faceOffset] = iFace;
            //offsetCell2Face update for the next set of face
            offsetCell2Face[cellLabelLeft]++;
            cell2faceOffset = offsetCell2Face[cellLabelRight];
            cell2Face[cell2faceOffset] = iFace;
            offsetCell2Face[cellLabelRight]++;
        }

        //Reset offsetCell2Face
        printf("Before reset offsetCell2Face[%d]=%d, nTotalCellFaces = %d\n", static_cast<int>(nTotalCell-1), offsetCell2Face[nTotalCell-1], nTotalCellFaces);
        offsetCell2Face[0] = 0;
        for (int iCell = 1; iCell < nTotalCell; ++ iCell)
        {
            //offsetCell2Face[iCell] = nFaceOfEachCell[iCell - 1];
            offsetCell2Face[iCell] = offsetCell2Face[iCell-1] + nFaceOfEachCell[iCell-1];
        }
        printf("After reset offsetCell2Face[%d]=%d, nFaceOfEachCell[%d] = %d, nTotalCellFaces = %d\n", static_cast<int>(nTotalCell-1), offsetCell2Face[nTotalCell-1], static_cast<int>(nTotalCell-1), nFaceOfEachCell[nTotalCell-1], nTotalCellFaces);
        //set cell2Cell
        //test begin
            //printf("rightCellOfFace[%d]=%d\n", 23221, rightCellOfFace[23221]);
        //test end
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            int sizeCellFace = nFaceOfEachCell[iCell];
            int offsetCellFace = offsetCell2Face[iCell];
            for (int iFace = 0; iFace < sizeCellFace; ++ iFace)
            {
                int faceID = cell2Face[offsetCellFace + iFace];
                int cellLabel = static_cast<int>(leftCellOfFace[faceID]);
                if ((cellLabel < 0) || (cellLabel >= nTotalCell)) 
                {
                    printf("Error: cellLabel = %d\n", cellLabel);
                    exit(1);
                }

                //test begin
                    //if (iCell == 0) printf("On the left, sizeCellFace =%d, iFace = %d, faceID = %d, cellLabel = %d, rightCell = %d\n", sizeCellFace, iFace, faceID, cellLabel, rightCellOfFace[faceID]);
                //test end

                if (cellLabel != iCell) cell2Cell[offsetCellFace + iFace] = cellLabel;

                cellLabel = static_cast<int>(rightCellOfFace[faceID]);
                //test begin
                    //if (iCell == 0) printf("On the right, faceID = %d, cellLabel = %d\n", faceID, cellLabel);
                //test end
                if (cellLabel == -1) continue; //boundy face's right neighbor
                if ((cellLabel < 0) || (cellLabel >= nTotalCell))
                {
                    printf("Error: cellLabel = %d\n", cellLabel);
                    exit(1);
                }

                if (cellLabel != iCell) cell2Cell[offsetCellFace + iFace] = cellLabel;
                
            }
        }

        //set bandwidthOfEachCell
        ComputeBandwidth(static_cast<int>(nTotalCell), nFaceOfEachCell, offsetCell2Face, cell2Cell);
        /*
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            int sizeCellFace = nFaceOfEachCell[iCell];
            int offsetCellFace = offsetCell2Face[iCell];

            //The max band between iCell and neighbor cells
            int maxLocalBand = 0;
            for (int iFace = 0; iFace < sizeCellFace; iFace++)
            {
                int cellID = cell2Cell[offsetCellFace+iFace];
                int bandDiff = abs(iCell - cellID);
                if (bandDiff > maxLocalBand) maxLocalBand = bandDiff;
                //test begin
                    //if (iCell == 0) printf("offsetCellFace = %d, iFace = %d, cellID = %d, bandDiff = %d\n", offsetCellFace, iFace, cellID, bandDiff);
                //test end
            }
            bandwidthOfEachCell[iCell] = maxLocalBand;
        }
        //statistics of bandwidthOfEachCell
        int maxGlobalBand = 0; //maximum bandwidth
        int maxCellID = 0; //cellID of maximum bandwidth
        int bandDist[33]; //distribution of bandwidth from 1 to 32 and more
        for (int i = 0; i < 33; i++)
        {
            bandDist[i] = 0;
        }

        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            int bandCell = bandwidthOfEachCell[iCell];
            if (bandCell <= 0) {
                printf("Error: bandwidth<=0 on cell %d, bandwidth=%d\n", iCell, bandCell);
                exit(1);
            }
            if (bandCell > maxGlobalBand) 
            {
                maxGlobalBand = bandCell;
                maxCellID = iCell;
            }
            if (bandCell > 32) bandDist[32]++;
            else bandDist[bandCell-1]++;
        }
        
        printf("maxGlobalBand = %d on cellID %d\n", maxGlobalBand, maxCellID);
        int nTotalDist = 0;
        for (int i = 0; i < 33; i++)
        {
            nTotalDist += bandDist[i];
            printf("bandDist %d: %d in nTotalCells\n", i+1, bandDist[i]);
        }

        printf("nTotalDist = %d, nTotalCell = %d\n", nTotalDist, nTotalCell);
        */
        delete [] bandwidthOfEachCell;    bandwidthOfEachCell = NULL;
    }
    
    void CallReorderFaceLabel(Base_Grid_Conn *gconn)
    {
        printf("Program is running in CallReorderFaceLabel\n");
        cgsize_t nTotalFace = gconn->GetNTotalFace();
        cgsize_t nBoundFace = gconn->GetNBoundFace();
        cgsize_t nTotalCell = gconn->GetNTotalCell();
        printf("nTotalFace = %d, nBoundFace = %d, nTotalCell = %d\n", static_cast<int>(nTotalFace), static_cast<int>(nBoundFace), static_cast<int>(nTotalCell));

        vector<cgsize_t> &leftCellOfFace  = *gconn->GetLeftCellOfFace();
        vector<cgsize_t> &rightCellOfFace = *gconn->GetRightCellOfFace();
        int *oldLeftCellOfFace = new int[nTotalFace];
        int *oldRightCellOfFace = new int[nTotalFace];
        
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            oldLeftCellOfFace[iFace] = static_cast<int>(leftCellOfFace[iFace]);
            oldRightCellOfFace[iFace] = static_cast<int>(rightCellOfFace[iFace]);
        }
        int *cell2Face;
        int *leftRightCellFaces;//corresponding to cell2Face, leftRightCellFaces record the cell is the left or right cell of one face. 

        //int * cell2Cell;
        //offset of iCell's cell faces in cell2Face
        int *offsetCell2Face = new int[nTotalCell];
        int nTotalCellFaces = 0;

        int *nFaceOfEachCell = new int[nTotalCell];

        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            nFaceOfEachCell[iCell] = 0;
            offsetCell2Face[iCell] = 0;
        }
        //set nFaceOfEachCell and nTotalCellFaces
        printf("set nFaceOfEachCell and nTotalCellFaces\n");
        SetFaceNumberOfEachCell(static_cast<int>(nTotalFace), static_cast<int>(nTotalCell), oldLeftCellOfFace, oldRightCellOfFace, nFaceOfEachCell, nTotalCellFaces);
        //set the size of cell2Cell and cell2Face
        printf("Before create of cell2Cell nTotalCellFaces = %d\n", nTotalCellFaces);
        cell2Face = new int[nTotalCellFaces];
        leftRightCellFaces = new int[nTotalCellFaces];
        //cell2Cell = new int[nTotalCellFaces];

        //set offsetCell2Face
        SetOffSetCell2Face(static_cast<int>(nTotalCell), nFaceOfEachCell, offsetCell2Face);
        //set cell2Face
        SetCell2Face(static_cast<int>(nTotalFace), static_cast<int>(nTotalCell), oldLeftCellOfFace, oldRightCellOfFace, offsetCell2Face, cell2Face);
        //Reset offsetCell2Face
        printf("Before reset offsetCell2Face[%d]=%d, nTotalCellFaces = %d\n", static_cast<int>(nTotalCell-1), offsetCell2Face[nTotalCell-1], nTotalCellFaces);
        //offsetCell2Face is modified due to update of cell2Face
        //set offsetCell2Face again
        SetOffSetCell2Face(static_cast<int>(nTotalCell), nFaceOfEachCell, offsetCell2Face);
        printf("After reset offsetCell2Face[%d]=%d, nFaceOfEachCell[%d] = %d, nTotalCellFaces = %d\n", static_cast<int>(nTotalCell-1), offsetCell2Face[nTotalCell-1], static_cast<int>(nTotalCell-1), nFaceOfEachCell[nTotalCell-1], nTotalCellFaces);
        //check cell2Face
        //rule 1: nBoundFace-1<cell2Face<nTotalFace
        //rule 2: in cell2Face, every value owns at least 1 item.
        printf("Checking cell2Face ...\n");
        int *faceIDRecord = new int[nTotalFace-nBoundFace];
        for (int iFace = 0; iFace < (nTotalFace-nBoundFace); ++ iFace)
        {
            faceIDRecord[iFace] = 0;
        }
        for (int cellID = 0; cellID < nTotalCell; ++ cellID)
        {
            int numFaces = nFaceOfEachCell[cellID];
            int offset = offsetCell2Face[cellID];
            for (int iFace=0; iFace < numFaces; ++ iFace)
            {
                int faceID = cell2Face[offset+iFace];
                //rule 1
                if ((faceID < nBoundFace) || (faceID > nTotalFace))
                {
                    printf("Error: faceID %d on cellID %d and iFace %d is outside of nBoundFace %d and nTotalFace-1 %d\n", faceID, cellID, iFace, static_cast<int>(nBoundFace), static_cast<int>(nTotalFace-1));
                    exit(1);
                }
                int remains = faceID - static_cast<int>(nBoundFace);
                faceIDRecord[remains]++;
            }
        }
        for (int iFace = 0; iFace < (nTotalFace-nBoundFace); ++ iFace)
        {
            int numFaces = faceIDRecord[iFace];
            //rule 2
            if (numFaces < 1) 
            {
                printf("Error: the number of faceID %d is %d\n", iFace+ static_cast<int>(nBoundFace), numFaces);
                exit(1);
            }
        }
        
        //set leftRightCellFaces, which stores face
        //SetLeftRightCellFaces();
        int numBoundaryFaces = 0;
        for (int cellID = 0; cellID < nTotalCell; ++ cellID)
        {
            int offset = offsetCell2Face[cellID];
            int numFaces = nFaceOfEachCell[cellID];
            for (int iFace = 0; iFace < numFaces; iFace++)
            {
                int faceID = cell2Face[offset+iFace];
                if (oldLeftCellOfFace[faceID] == -1)
                {
                    printf("Error: on faceID=%d, leftCellOfFace is -1\n", faceID);
                    exit(1);
                }
                if (oldRightCellOfFace[faceID] == -1) numBoundaryFaces ++;
                //if (oldRightCellOfFace[faceID] == -1) numBoundaryFaces++;
                if (oldLeftCellOfFace[faceID] == cellID) leftRightCellFaces[offset+iFace] = 1;
                if (oldRightCellOfFace[faceID] == cellID) leftRightCellFaces[offset+iFace] = -1;
                
            }
        }
        //if (numBoundaryFaces != nBoundFace) {
        //Because cell2Face does not contain boudnary faces
        if (numBoundaryFaces != 0)
        {
            printf("Error: numBoundaryFaces=%d is not equal to 0\n", numBoundaryFaces);
            exit(1);
        }

        //Reorder faceID from nBoundFace to nTotalFace
        int *newFaces = new int[nTotalFace]; //Old faceID and containerID
        //set the initial value of labelFace as nBoundFace-1
        int labelFace = static_cast<int>(nBoundFace)-1;
        //Initialize newFaces by -1, which means not changed.
        for (int faceID = 0; faceID < nTotalFace; ++ faceID)
        {
            newFaces[faceID] = -1;
        }
        
        int faceReorderMethod = PHSPACE::GlobalDataBase::GetIntParaFromDB("faceReorderMethod");
        if (faceReorderMethod == 0)
        {
            //BSFCELLFACE means Broad Search First: faceID is distributed face by face in one cell with the ascending order. One cell starts before faceID distribution in one cell is finished. In fact, it is based on original faceID order. Before ReOrder, cell2Face has been reordered with old faceID in ascending order.
            //Ascending cell2Face for each cell
            printf("BSFCELLFACEORG: Ascending cell2Face order...\n");
            for (int cellID = 0; cellID < nTotalCell; ++ cellID)
            {
                int offset = offsetCell2Face[cellID];
                int numFaces = nFaceOfEachCell[cellID];
                int loopNum = numFaces - 1;
                //bubble sort
                for (int jFace = 0; jFace < numFaces-1; ++ jFace)
                {
                    for (int iFace = 0; iFace < loopNum; ++ iFace)
                    {
                        int faceID = cell2Face[offset+iFace];
                        int faceIDNext = cell2Face[offset+iFace+1];
                        if (faceID > faceIDNext)
                        {
                            cell2Face[offset+iFace] = faceIDNext;
                            cell2Face[offset+iFace+1] = faceID;
                            int temp = leftRightCellFaces[offset+iFace];
                            leftRightCellFaces[offset+iFace] = leftRightCellFaces[offset+iFace+1];
                            leftRightCellFaces[offset+iFace+1] = temp;
                        }
                    }

                    loopNum --;

                }
            }
            printf("BSFCELLFACEORG: checking cell2Face order...\n");
            for (int cellID = 0; cellID < nTotalCell; ++ cellID)
            {
                int offset = offsetCell2Face[cellID];
                int numFaces = nFaceOfEachCell[cellID];
                for (int iFace = 0; iFace < numFaces-1; iFace++)
                {
                    int faceID = cell2Face[offset+iFace];
                    int faceIDNext = cell2Face[offset+iFace+1];
                    if (faceID > faceIDNext)
                    {
                        printf("Error: faceID %d on cellID %d iFace %d is larger than faceIDNext %d\n", faceID, cellID, iFace, faceIDNext);
                        exit(1);
                    }
                }
            }
            printf("BSFCELLFACEORG: Reorder cell labels ...\n");
            CellReorderByCell2FaceOrder(static_cast<int>(nTotalCell), static_cast<int>(nTotalFace), offsetCell2Face, nFaceOfEachCell, cell2Face, static_cast<int>(nBoundFace), labelFace, newFaces);
        }else if (faceReorderMethod == 1)
        {
            //It is quite similar as BSFCELLFACEORG. Firstly, cell2Face is reorded with old faceID with leftCellOfFace in ascending order. Then, set those faceID from nBoundFace to nTotalFace
            //Ascending cell2Face for each cell
            printf("BSFCELLFACELEFT: Ascending cell2Face order by leftCellOfFace ...\n");
            for (int cellID = 0; cellID < nTotalCell; ++ cellID)
            {
                int offset = offsetCell2Face[cellID];
                int numFaces = nFaceOfEachCell[cellID];
                int loopNum = numFaces-1;
                //bubble sort
                for (int jFace = 0; jFace < numFaces-1; ++ jFace)
                {
                    for (int iFace = 0; iFace < loopNum; ++ iFace)
                    {
                        int faceID = cell2Face[offset+iFace];
                        int leftCell = oldLeftCellOfFace[faceID];
                        int faceIDNext = cell2Face[offset+iFace+1];
                        int leftCellNext = oldLeftCellOfFace[faceIDNext];
                        if (leftCell > leftCellNext){
                            cell2Face[offset+iFace] = faceIDNext;
                            cell2Face[offset+iFace+1] = faceID;
                            int temp = leftRightCellFaces[offset+iFace];
                            leftRightCellFaces[offset+iFace] = leftRightCellFaces[offset+iFace+1];
                            leftRightCellFaces[offset+iFace+1] = temp;
                        }
                    }
                    loopNum --;
                }
            }
            printf("BSFCELLFACELEFT: checking cell2Face order...\n");
            for (int cellID = 0; cellID < nTotalCell; ++ cellID)
            {
                int offset = offsetCell2Face[cellID];
                int numFaces = nFaceOfEachCell[cellID];
                for (int iFace = 0; iFace < numFaces-1; ++ iFace)
                {
                    int faceID = cell2Face[offset+iFace];
                    int faceIDNext = cell2Face[offset+iFace+1];
                    int leftCell = oldLeftCellOfFace[faceID];
                    int leftCellNext = oldLeftCellOfFace[faceIDNext];
                    if (leftCell > leftCellNext){
                        printf("Error: faceID %d on cellID %d iFace %d with leftCell %d is larger than leftCellNext %d\n", faceID, cellID, iFace, leftCell, leftCellNext);
                        exit(1);
                    }
                }
            }
            //Reorder cell label by faces in cell2Face
            CellReorderByCell2FaceOrder(static_cast<int>(nTotalCell), static_cast<int>(nTotalFace), offsetCell2Face, nFaceOfEachCell, cell2Face, static_cast<int>(nBoundFace), labelFace, newFaces);
        }else if (faceReorderMethod == 2)
        {
            //It is quite similar as BSFCELLFACEORG. Firstly, cell2Face is reorded with old faceID with rightCellOfFace in ascending order. Then, set those faceID from nBoundFace to nTotalFace
            //Ascending cell2Face for each cell
            printf("BSFCELLFACERIGHT: Ascending cell2Face order by rightCellOfFace ...\n");
            for (int cellID = 0; cellID < nTotalCell; ++ cellID)
            {
                int offset = offsetCell2Face[cellID];
                int numFaces = nFaceOfEachCell[cellID];
                int loopNum = numFaces-1;
                //bubble sort
                for (int jFace = 0; jFace < numFaces-1; ++ jFace)
                {
                    for (int iFace = 0; iFace < loopNum; ++ iFace)
                    {
                        int faceID = cell2Face[offset+iFace];
                        int rightCell = oldRightCellOfFace[faceID];
                        int faceIDNext = cell2Face[offset+iFace+1];
                        int rightCellNext = oldRightCellOfFace[faceIDNext];
                        if (rightCell > rightCellNext)
                        {
                            cell2Face[offset+iFace] = faceIDNext;
                            cell2Face[offset+iFace+1] = faceID;
                            int temp = leftRightCellFaces[offset+iFace];
                            leftRightCellFaces[offset+iFace] = leftRightCellFaces[offset+iFace+1];
                            leftRightCellFaces[offset+iFace+1] = temp;
                        }
                    }
                    loopNum --;
                }
            }
            printf("BSFCELLFACERIGHT: checking cell2Face order...\n");
            for (int cellID = 0; cellID < nTotalCell; ++ cellID)
            {
                int offset = offsetCell2Face[cellID];
                int numFaces = nFaceOfEachCell[cellID];
                for (int iFace = 0; iFace < numFaces-1; ++ iFace)
                {
                    int faceID = cell2Face[offset+iFace];
                    int faceIDNext = cell2Face[offset+iFace+1];
                    int rightCell = oldRightCellOfFace[faceID];
                    int rightCellNext = oldRightCellOfFace[faceIDNext];
                    if (rightCell > rightCellNext)
                    {
                        printf("Error: faceID %d on cellID %d iFace %d with rightCell %d is larger than rightCellNext %d\n", faceID, cellID, iFace, rightCell, rightCellNext);
                        exit(1);
                    }
                }
            }
            //Reorder cell label by faces in cell2Face
            CellReorderByCell2FaceOrder(static_cast<int>(nTotalCell), static_cast<int>(nTotalFace), offsetCell2Face, nFaceOfEachCell, cell2Face, static_cast<int>(nBoundFace), labelFace, newFaces);
        }

        //check reorder results
        //rule 1 Before nBoundFace, newFaces=-1
        //rule 2 nBoundFace <= newFaces < nTotalFace
        //rule 3 After nBoundFace, no elements are the same in newFaces
        printf("Checking newFaces which is Reorder cell labels results...\n");
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            //rule 1
            if (newFaces[iFace]!=-1) {
                printf("Error in rule 1: on iFace %d, newFaces is -1\n", iFace);
                exit(1);
            }
        }
        int *recordFaceID = new int[nTotalFace-nBoundFace];
        for (int iFace = 0; iFace < nTotalFace-nBoundFace; ++ iFace)
        {
            recordFaceID[iFace] = 0;
        }
        for (int iFace = static_cast<int>(nBoundFace); iFace < nTotalFace-nBoundFace; ++ iFace)
        {
            int faceID = newFaces[iFace];
            //rule 2
            if ((faceID < nBoundFace) || (faceID>(nTotalFace-1)))
            {
                printf("Error in rule 2: on iFace %d of newFaces, faceID = %d is outside of nBoundFace %d and nTotalFace - 1 %d\n", iFace, faceID, static_cast<int>(nBoundFace), static_cast<int>(nTotalFace-1));
                exit(1);
            }
            int remains = (static_cast<int>(nTotalFace)-1) - faceID;
            recordFaceID[remains] ++;
        }
        //rule 3
        for (int iFace = 0; iFace < nTotalFace-nBoundFace; ++ iFace)
        {
            if (recordFaceID[iFace] > 1) 
            {
                printf("Error: more than two faces own the same faceID %d\n", static_cast<int>(nTotalFace) - 1 - iFace);
                exit(1);
            }
        }
        
        //set newFaces bwtween 0 and nBoundFace-1
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            newFaces[iFace] = iFace;
        }
    
        //create reflectFaces for recording relationship of containerID and old cellID
        int *reflectFaces = new int[nTotalFace];
        for (int oldFaceID = 0; oldFaceID < nTotalFace; ++ oldFaceID)
        {
            int containerID = newFaces[oldFaceID];
            reflectFaces[containerID] = oldFaceID;
        }
        //check reflectFaceID
        //rule 1: 0~nBoundFace, reflectFaces[iFace]=iFace
        //rule 2: nBoundFace~nTotalFace, nBoundFace<=reflectFaces<nTotalFace
        //rule 3: 0~nTotalFace, reflectFaces does not own the same old faceID
        printf("Checking reflectFaceID  results...\n");
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            //rule 1
            if (reflectFaces[iFace] != iFace) 
            {
                printf("Error: on containerID %d, is not the same faceID\n", iFace);
                exit(1);
            }
        }
        for (int iFace = 0; iFace < nTotalFace-nBoundFace; ++ iFace)
        {
            recordFaceID[iFace] = 0;
        }
        for (int iFace = static_cast<int>(nBoundFace); iFace < nTotalFace; ++ iFace)
        {
            int oldFaceID = reflectFaces[iFace];
            if ((oldFaceID < nBoundFace)||(oldFaceID > nTotalFace))
            {
                //rule 2
                printf("Error: on containerID %d, oldFaceID %d is outside nBoundFace and nTotalFace\n", iFace, oldFaceID);
                exit(1);
            }
            int remains = oldFaceID - static_cast<int>(nBoundFace);
            recordFaceID[remains] ++;
        }
        for (int iFace = 0; iFace < nTotalFace-nBoundFace; ++ iFace)
        {
            if (recordFaceID[iFace] != 1) 
            {
                //rule 3
                printf("Error: on iFace %d, face number %d is not 1\n", iFace+ static_cast<int>(nBoundFace), recordFaceID[iFace]);
                exit(1);
            }
        }
        
        //swap method for updating nNodeOfEachFace, face2Node,
        //leftCellOfFace and rightCellOfFace
        vector< vector<cgsize_t> > & face2Node = *gconn->GetFace2Node();
        vector<int> & nNodeOfEachFace = *gconn->GetNodeNumberOfEachFace();
                
        int *orgFace = new int[nTotalFace]; //containerID and old faceID reflection
        int *swapFace = new int[nTotalFace];//old faceID and containerID reflection
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            orgFace[iFace] = iFace;
            swapFace[iFace] = iFace;
        }

        for (int faceID = 0; faceID < nTotalFace; ++ faceID)
        {
            if (orgFace[faceID] == reflectFaces[faceID]) continue;
            int targetID = reflectFaces[faceID];
            int currentID = orgFace[faceID];
            int swapID = swapFace[targetID];
            //update nNodeOfEachFace on faceID and swapID
            int tempNumNode = nNodeOfEachFace[swapID];
            nNodeOfEachFace[swapID] = nNodeOfEachFace[faceID];
            nNodeOfEachFace[faceID] = tempNumNode;
            //update face2Node on faceID and swapID
            vector<cgsize_t> tempFace = face2Node[swapID];
            face2Node[swapID] = face2Node[faceID];
            face2Node[faceID] = tempFace;
            //update leftCellOfFace on faceID and swapID
            int tempCell = static_cast<int>(leftCellOfFace[swapID]);
            leftCellOfFace[swapID] = leftCellOfFace[faceID];
            leftCellOfFace[faceID] = tempCell;
            //update rightCellOfFace on faceID and swapID
            tempCell = static_cast<int>(rightCellOfFace[swapID]);
            rightCellOfFace[swapID] = rightCellOfFace[faceID];
            rightCellOfFace[faceID] = tempCell;
            //update orgFace on faceID and swapID
            orgFace[faceID] = targetID;
            orgFace[swapID] = currentID;
            //update swapID on targetID and currentID
            swapFace[targetID] = faceID;
            swapFace[currentID] = swapID;
        }
        //check orgFace and swapFace
        printf("Check orgFace and swapFace ...\n");
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            if (orgFace[iFace] != reflectFaces[iFace])
            {
                printf("Error: on iFace %d, orgFace %d is not equal to reflectFaces %d\n", iFace, orgFace[iFace], reflectFaces[iFace]);
                exit(1);
            }

            if (swapFace[iFace] != newFaces[iFace])
            {
                printf("Error:  on iFace %d, swapFace %d is not equal to newFaces %d\n", iFace, swapFace[iFace], newFaces[iFace]);
                exit(1);
            }
        }
        
    }

    void CallReorderCellLabel(Base_Grid_Conn * gconn)
    {
        printf("Program is running in CallReorderCellLabel\n");
        cgsize_t nTotalFace = gconn->GetNTotalFace();
        cgsize_t nBoundFace = gconn->GetNBoundFace();
        cgsize_t nTotalCell = gconn->GetNTotalCell();
        printf("nTotalFace = %d, nBoundFace = %d, nTotalCell = %d\n", static_cast<int>(nTotalFace), static_cast<int>(nBoundFace), static_cast<int>(nTotalCell));

        vector<cgsize_t> &leftCellOfFace  = *gconn->GetLeftCellOfFace();
        vector<cgsize_t> &rightCellOfFace = *gconn->GetRightCellOfFace();
        int *oldLeftCellOfFace = new int[nTotalFace];
        int *oldRightCellOfFace = new int[nTotalFace];
        
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            oldLeftCellOfFace[iFace] = static_cast<int>(leftCellOfFace[iFace]);
            oldRightCellOfFace[iFace] = static_cast<int>(rightCellOfFace[iFace]);
        }

        int *cell2Face;
        int *cell2Cell;
        //offset of iCell's cell faces in cell2Face
        int *offsetCell2Face = new int[nTotalCell];
        int nTotalCellFaces = 0;

        int *nFaceOfEachCell = new int[nTotalCell];

        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            nFaceOfEachCell[iCell] = 0;
            offsetCell2Face[iCell] = 0;
        }

        //set nFaceOfEachCell and nTotalCellFaces
        printf("set nFaceOfEachCell and nTotalCellFaces\n");
        SetFaceNumberOfEachCell(static_cast<int>(nTotalFace), static_cast<int>(nTotalCell), oldLeftCellOfFace, oldRightCellOfFace, nFaceOfEachCell, nTotalCellFaces);
        
        //set the size of cell2Cell and cell2Face
        printf("Before create of cell2Cell nTotalCellFaces = %d\n", nTotalCellFaces);
        cell2Face = new int[nTotalCellFaces];
        cell2Cell = new int[nTotalCellFaces];

        //set offsetCell2Face
        SetOffSetCell2Face(static_cast<int>(nTotalCell), nFaceOfEachCell, offsetCell2Face);
        
        //set cell2Face
        SetCell2Face(static_cast<int>(nTotalFace), static_cast<int>(nTotalCell), oldLeftCellOfFace, oldRightCellOfFace, offsetCell2Face, cell2Face);
        
        //Reset offsetCell2Face
        printf("Before reset offsetCell2Face[%d]=%d, nTotalCellFaces = %d\n", static_cast<int>(nTotalCell)-1, offsetCell2Face[nTotalCell-1], nTotalCellFaces);
        
        //offsetCell2Face is modified due to update of cell2Face
        //set offsetCell2Face again
        SetOffSetCell2Face(static_cast<int>(nTotalCell), nFaceOfEachCell, offsetCell2Face);
        
        printf("After reset offsetCell2Face[%d]=%d, nFaceOfEachCell[%d] = %d, nTotalCellFaces = %d\n", static_cast<int>(nTotalCell)-1, offsetCell2Face[nTotalCell-1], static_cast<int>(nTotalCell)-1, nFaceOfEachCell[nTotalCell-1], nTotalCellFaces);

        //set cell2Cell
        SetCell2Cell(static_cast<int>(nTotalCell), nFaceOfEachCell, offsetCell2Face, cell2Face, oldLeftCellOfFace, oldRightCellOfFace, cell2Cell);
        
        //check cell2Cell, to see some rules:
        //rule 1 one cell's neighbor cells should not contain itself.
        //rule 2 0<= one cell's neighbor cells < nTotalCell
        //rule 3 one cell's neighbor cells should be different.
        //It should be noted that for O type grid, break 3 may be broken.
        CheckCell2Cell(static_cast<int>(nTotalCell), nFaceOfEachCell, offsetCell2Face, cell2Cell);
        //By above process, cell2Cell is established, which stores one cell's neighbour cells' labels. offsetCell2Face can be used by access of cell2Cell. nFaceOfEachCell stores face number owned by one cell, which can also be used for neighbour cell number owned by one cell. Boundary faces are not considered in both cell2Cell, offsetCell2Face and nFaceOfEachCell.
        //For reording cell index, cell2Cell with offsetCell2Face can supply cell connection information. nFaceOfEachCell can supply degree of cell. 
        int maxDegree = 0;
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            int degreeCell = nFaceOfEachCell[iCell];
            if (maxDegree < degreeCell) maxDegree = degreeCell;
        }
        printf("The maximum degree of cell is %d\n", maxDegree);
        //establish containers
        int *reorderCell = new int[nTotalCell]; //reorder cell labels
        int *newCell = new int[nTotalCell]; //reorder cell labels
        //int maxSizeLevel = maxDegree * maxDegree; //maximum number of one level
        //It may be very large, not just maxDegree*maxDegree
        //int maxSizeLevel = 5000; //maximum number of one level
        //The worst condition
        int maxSizeLevel = static_cast<int>(nTotalCell); //maximum number of one level
        int *cellIDSortCur = new int[maxSizeLevel];
        int numSortCur = 0;
        int *cellIDSortNext = new int[maxSizeLevel];
        int numSortNext = 0;
        //Initialization of reorderCell by -1
        for (int iCell = 0; iCell < nTotalCell; ++iCell) reorderCell[iCell] = -1;
        //Initialization of newCell by -1
        for (int iCell = 0; iCell < nTotalCell; ++ iCell) newCell[iCell] = -1;
        //Initialization of cellIDSortCur and cellIDSortNext by -1
        for (int iCell = 0; iCell < maxSizeLevel; ++ iCell) 
        {
            cellIDSortCur[iCell] = -1;
            cellIDSortNext[iCell] = -1;
        }
        
        printf("Searching the root vertex ...\n");
        //start from vetex 0
        int rStart = 0;
        printf("At the beginning, rStart = %d\n", rStart);
        //Here, reorderCell is used for marking a cell exsiting in root level structure
        reorderCell[rStart] = 1;
        int lCur = 0;
        int lOld = -1;
        //Get current level of vertex r
        GetCurrentLevel(rStart, nFaceOfEachCell, cell2Cell, offsetCell2Face, reorderCell, numSortCur, cellIDSortCur);
        
        //Get the next level of vertex r
        GetNextLevel(numSortCur, cellIDSortCur, nFaceOfEachCell, cell2Cell, offsetCell2Face, reorderCell, numSortNext, cellIDSortNext);
        
        while (lCur > lOld)
        {
            lOld = lCur;
            lCur = 0;
            while (numSortNext > 0)
            {
                lCur ++;
                
                ResetCellSort(numSortCur, cellIDSortCur);
                SetCellIDSortCurByCellIDSortNext(numSortNext, cellIDSortNext, numSortCur, cellIDSortCur);
                ResetCellSort(numSortNext, cellIDSortNext);
                GetNextLevel(numSortCur, cellIDSortCur, nFaceOfEachCell, cell2Cell, offsetCell2Face, reorderCell, numSortNext, cellIDSortNext);
            }
            if (lCur > lOld)
            {
                
                //get the minimum degree in cellIDSortCur, start the next iteration
                rStart = FindMinCellIDSortCur(numSortCur, cellIDSortCur, nFaceOfEachCell);
                //printf("rStart = %d in numSortCur %d\n", rStart, numSortCur);
                //Reset R by -1
                for (int iCell = 0; iCell < nTotalCell; ++ iCell) reorderCell[iCell] = -1;
                //Reset numSortCur by -1
                ResetCellSort(numSortCur, cellIDSortCur);
                //Reset numSortNext by -1
                ResetCellSort(numSortNext, cellIDSortNext);

                reorderCell[rStart] = 1;

                GetCurrentLevel(rStart, nFaceOfEachCell, cell2Cell, offsetCell2Face, reorderCell, numSortCur, cellIDSortCur);

                GetNextLevel(numSortCur, cellIDSortCur, nFaceOfEachCell, cell2Cell, offsetCell2Face, reorderCell, numSortNext, cellIDSortNext);
                    
                
            }
            printf("lCur = %d, lOld = %d, numSortNext=%d, rStart = %d\n", lCur, lOld, numSortNext ,rStart);
        }

        printf("Reorder cell index ...\n");
        printf("Start from %d\n", rStart);
        //Reset reorderCell by -1
        for (int iCell = 0; iCell < nTotalCell; ++ iCell) reorderCell[iCell] = -1;
        //Reset numSortCur by -1
        ResetCellSort(numSortCur, cellIDSortCur);

        //Reset numSortNext by -1
        ResetCellSort(numSortNext, cellIDSortNext);
    
        //set rStart as 0, and numVertex as 1
        reorderCell[rStart] = 1;
        int numVertex = 0;
        newCell[rStart] = numVertex;
        int nLevel = 0;
        
        //Get current level of rStart
        //GetCurrentLevelForReOrder(rStart, nFaceOfEachCell, cell2Cell, offsetCell2Face, reorderCell, numSortCur, cellIDSortCur);
        GetCurrentLevel(rStart, nFaceOfEachCell, cell2Cell, offsetCell2Face, reorderCell, numSortCur, cellIDSortCur);

        while(numSortCur > 0)
        {
            printf("nLevel = %d, numSortCur = %d\n", nLevel, numSortCur);
         
            AscendingCellIDSortCurByDegree(numSortCur, nFaceOfEachCell, cellIDSortCur);
            ReOrderLevelCur(numSortCur, cellIDSortCur, numVertex, newCell);

            GetNextLevel(numSortCur, cellIDSortCur, nFaceOfEachCell, cell2Cell, offsetCell2Face, reorderCell, numSortNext, cellIDSortNext);

            ResetCellSort(numSortCur, cellIDSortCur);
            SetCellIDSortCurByCellIDSortNext(numSortNext, cellIDSortNext, numSortCur, cellIDSortCur);
            ResetCellSort(numSortNext, cellIDSortNext);
            nLevel ++;

        }
        printf("Finish cell reorder by RCM scheme...\n");
        printf("nToatalCell = %d, numVertex = %d\n", static_cast<int>(nTotalCell), numVertex);
        //reverse newCell
        //Not sure
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            newCell[iCell] = static_cast<int>(nTotalCell) - 1 - newCell[iCell];
        }
        //! Check rules for cell index
        CheckNewCell(static_cast<int>(nTotalCell), newCell);

        int *reflectCell = new int[nTotalCell];
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            int cellID = newCell[iCell];
            reflectCell[cellID] = iCell;
        }
        //define new nFaceOfEachCell and new offsetCell2Face
        int *nNewFaceOfEachCell = new int[nTotalCell];
        int *newOffsetCell2Face = new int[nTotalCell];
        //set nNewFaceOfEachCell
        for (int cellID = 0; cellID < nTotalCell; ++ cellID)
        {
            int reflectCellID = reflectCell[cellID];
            nNewFaceOfEachCell[cellID] = nFaceOfEachCell[reflectCellID];
        }

        SetOffSetCell2Face(static_cast<int>(nTotalCell), nNewFaceOfEachCell, newOffsetCell2Face);
        //! Define new leftCellOfFace and rightCellOfFace
        //! Define new cell2Cell and cell2Face
        printf("Before create of newCell2Cell, nTotalCellFaces = %d\n", nTotalCellFaces);
        int *newCell2Cell = new int[nTotalCellFaces];
        
        //! Set newCell2Cell
        for (int cellID = 0; cellID < nTotalCell; ++ cellID)
        {
            int reflectCellID = reflectCell[cellID];
            int numCellCellsOld = nFaceOfEachCell[reflectCellID];
            int numCellCellsNew = nNewFaceOfEachCell[cellID];
            if (numCellCellsOld != numCellCellsNew) 
            {
                printf("Error: on cellID %d, numCellCellsNew %d is not equal to numCellCellsOld %d on reflectCellID %d\n", cellID, numCellCellsNew, numCellCellsOld, reflectCellID);
                exit(1);
            }
            int offsetNew = newOffsetCell2Face[cellID];
            int offsetOld = offsetCell2Face[reflectCellID];
            for (int iCell = 0; iCell < numCellCellsNew; ++ iCell)
            {
                int newCellID = newCell[cell2Cell[offsetOld+iCell]];
                newCell2Cell[offsetNew+iCell] = newCellID;
            }
        }
    
        //! Check newCell2Cell
        CheckCell2Cell(static_cast<int>(nTotalCell), nNewFaceOfEachCell, newOffsetCell2Face, newCell2Cell);
        
        //! Calculate bandwidth again    
        ComputeBandwidth(static_cast<int>(nTotalCell), nNewFaceOfEachCell, newOffsetCell2Face, newCell2Cell);
        
        //! According to cell reorder, Update face_number_of_each_face by swap
        //! Input: face_nuber_of_each_cell and reflectCell
        SwapFaceNumberOfEachFace(static_cast<int>(nTotalCell), nFaceOfEachCell, reflectCell, nNewFaceOfEachCell);
        
        //! Update leftCellOfFace and rightCellOfFace
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            //! For leftCellOfFace
            int oldCellID = static_cast<int>(leftCellOfFace[iFace]);
            //! Not bounday faces, treatment is requierd
            if (oldCellID != -1)
            {
                int NewCellID = newCell[oldCellID];
                leftCellOfFace[iFace] = NewCellID;
            }
            //! For rightCellOfFace
            oldCellID = static_cast<int>(rightCellOfFace[iFace]);
            if (oldCellID != -1)
            {
                int NewCellID = newCell[oldCellID];
                rightCellOfFace[iFace] = NewCellID;

            }
        }
        //! Update cell2Node and nNodeOfEachCell by swap method
        vector< vector<cgsize_t> > &cell2Node = *gconn->GetCell2Node();
        vector<int> &nNodeOfEachCell = *gconn->GetNodeNumberOfEachCell();
        int *orgCell = new int[nTotalCell];//! Contianer index and old cellID reflection
        int *swapCell = new int[nTotalCell]; //! Old cellID and container index reflection
        for (int cellID = 0; cellID < nTotalCell; ++ cellID)
        {
            //! Before swap, orgCell is equal to swapCell
            orgCell[cellID] = cellID; 
            swapCell[cellID] = cellID;
        }

        for (int cellID = 0; cellID < nTotalCell; ++ cellID)
        {
            if (orgCell[cellID] == reflectCell[cellID]) continue;
            int targetID = reflectCell[cellID];
            int currentID = orgCell[cellID];
            int swapID = swapCell[targetID];
            //! Update correspondings of cellID and swapID for nNodeOfEachCell
            int tempNodeNumber = nNodeOfEachCell[swapID];

            nNodeOfEachCell[swapID] = nNodeOfEachCell[cellID];
            nNodeOfEachCell[cellID] = tempNodeNumber;
            //! Update correspondings of cellID and swapID for cell2Node
            vector<cgsize_t> tempCell = cell2Node[swapID];
            cell2Node[swapID] = cell2Node[cellID];
            cell2Node[cellID] = tempCell;
            //! Update correspondings of cellID and swapID for orgCell
            orgCell[cellID] = targetID;
            orgCell[swapID] = currentID;
            //! Update correspondings of targerID and currentID for swapID
            swapCell[targetID] = cellID;
            swapCell[currentID] = swapID;
        }
    }

    void GetCurrentLevel(int rStart, const int *nFaceOfEachCell, const int *cell2Cell, const int *offsetCell2Cell, int *reorderCell, int &numSortCur, int *cellIDSortCur)
    {
        numSortCur = 0;
        int numCellCells = nFaceOfEachCell[rStart];
        int offset = offsetCell2Cell[rStart];

        for (int iCell = 0; iCell < numCellCells; ++ iCell)
        {
            int cellID = cell2Cell[offset + iCell];
            //! cellID is not in the root level structure
            if (reorderCell[cellID] == -1)
            {
                cellIDSortCur[numSortCur] = cellID;
                numSortCur ++;
                reorderCell[cellID] = 1;
            }
        }
    }

    void GetNextLevel(const int &numSortCur, const int *cellIDSortCur, const int *nFaceOfEachCell, const int *cell2Cell, const int *offsetCell2Cell, int *reorderCell, int &numSortNext, int *cellIDSortNext)
    {
        numSortNext = 0;
        for (int iCell = 0; iCell < numSortCur; ++ iCell)
        {
            int cellCellID = cellIDSortCur[iCell];
            int numCellCells = nFaceOfEachCell[cellCellID];
            int offset = offsetCell2Cell[cellCellID];
            for (int iCellCell = 0; iCellCell < numCellCells; iCellCell++)
            {
                int cellID = cell2Cell[offset + iCellCell];
                if (reorderCell[cellID] == -1)
                {
                    cellIDSortNext[numSortNext] = cellID;
                    numSortNext ++;
                    reorderCell[cellID] = 1;
                }
            }
        }
    }

    void SetCellIDSortCurByCellIDSortNext(const int &numSortNext, const int *cellIDSortNext, int &numSortCur, int *cellIDSortCur)
    {
        numSortCur = numSortNext;

        for (int iSort=0; iSort < numSortNext; ++ iSort)
        {
            cellIDSortCur[iSort] = cellIDSortNext[iSort];
        }
    }

    void ResetCellSort(int &numSort, int *cellIDSort)
    {
        for (int iSort=0; iSort<numSort; ++ iSort)
        {
            cellIDSort[iSort] = -1;
        }

        numSort = 0;
    }

    void TestCellSort(const int numSort, const int *cellIDSort)
    {
        for (int iSort = 0; iSort<numSort; ++ iSort)
        {
            printf("iSort = %d of numSort %d, cellID = %d\n", iSort, numSort, cellIDSort[iSort]);
        }
    }

    int FindMinCellIDSortCur(const int numSortCur, const int *cellIDSortCur, const int *nFaceOfEachCell)
    {
        int minCellID = cellIDSortCur[0];
        int minDegree = nFaceOfEachCell[minCellID];

        for (int iSort = 1; iSort < numSortCur; ++ iSort)
        {
            int cellID = cellIDSortCur[iSort];
            if (nFaceOfEachCell[iSort] < minDegree)
            {
                minDegree = nFaceOfEachCell[iSort];
                minCellID = cellID;
            }

        }

        return minCellID;
    }

    void GetCurrentLevelForReOrder(int rStart, const int *nFaceOfEachCell, const int *cell2Cell, const int *offsetCell2Cell, int *reorderCell, int &numSortCur, int *cellIDSortCur)
    {
        numSortCur = 0;
        int numCellCells = nFaceOfEachCell[rStart];
        int offset = offsetCell2Cell[rStart];
       
        for (int iCell = 0; iCell < numCellCells; ++ iCell)
        {
            int cellID = cell2Cell[offset + iCell];
            if (reorderCell[cellID] == -1)
            {
                cellIDSortCur[numSortCur] = cellID;
                numSortCur ++;
            }
        }
    }
    
    void AscendingCellIDSortCurByDegree(const int numSortCur, const int *nFaceOfEachCell, int *cellIDSortCur)
    {
        int numVertexInSort = numSortCur;
        //! At least 2 vertices for sort
        while(numVertexInSort > 1)
        {
            for (int iSort=0; iSort < numVertexInSort-1; ++ iSort)
            {
                int cellID = cellIDSortCur[iSort];
                int cellIDNext = cellIDSortCur[iSort+1];
                int degree = nFaceOfEachCell[cellID];
                int degreeNext = nFaceOfEachCell[cellIDNext];
                if (degree > degreeNext)
                {
                    //! Exchange vertex order
                    cellIDSortCur[iSort] = cellIDNext;
                    cellIDSortCur[iSort+1] = cellID;
                }
            }
            numVertexInSort --;
        }
    }

    void TestCellIDSortCurDegree(const int numSortCur, const int *cellIDSortCur, const int *nFaceOfEachCell)
    {
        //! Test Ascending order of vertex 0
        for (int iSort = 0; iSort < numSortCur; ++ iSort)
        {
            printf("iSort = %d, cellID = %d, degree = %d\n", iSort, cellIDSortCur[iSort], nFaceOfEachCell[cellIDSortCur[iSort]]);
        }
    }

    //! Check cell2Cell, to see some rules:
    //! rule1 one cell's neighbor cells should not contain itself.
    //! rule2 0<= one cell's neighbor cells < nTotalCell
    //! rule3 one cell's neighbor cells should be different.
    //! It should be noted that for O type grid, break 3 may be broken.
    void CheckCell2Cell(const int nTotalCell, const int *nFaceOfEachCell, const int *offsetCell2Cell, const int *cell2Cell)
    {
        printf("Check cell2Cell...\n");

        for (int cellID = 0; cellID < nTotalCell; ++ cellID)
        {
            int offset = offsetCell2Cell[cellID];
            int numCellCells = nFaceOfEachCell[cellID];

            for (int iCell = 0; iCell < numCellCells; ++ iCell)
            {
                int cellCellID = cell2Cell[offset + iCell];
                //! Check rule 1
                if (cellID == cellCellID) printf("rule 1 breaks in cellID %d in %d\n", cellID, iCell);
                //! Check rule 2
                if ((cellCellID < 0)||(cellCellID >= nTotalCell)) printf("rule 2 breaks in cellID %d in %d for %d\n", cellID, iCell, cellCellID);
            }

            //! Check rule 3
            int numComp = numCellCells;
            //! At least two cells in cell2Cell
            while (numComp > 1)
            {
                for (int jCell = 0; jCell < numComp-1; ++ jCell)
                {
                    int cellIDLoop = cell2Cell[offset+jCell];
                    int cellIDComp = cell2Cell[offset+numComp-1];
                    if (cellIDLoop == cellIDComp) printf("rule 3 breaks in cellID %d in %d and %d are both %d\n", cellID, jCell, numComp-1, cellIDLoop);
                }
                numComp --;
            }
        }

    }

    void ReOrderLevelCur(const int numSortCur, const int *cellIDSortCur, int &numVertex, int *newCell)
    {
        for (int iSort = 0; iSort < numSortCur; ++ iSort)
        {
            int cellID = cellIDSortCur[iSort];
            //! Avoid the same cell index in cellIDSortCur on the condition of O type grid
            if (newCell[cellID] == -1)
            {
                numVertex ++;
                newCell[cellID] = numVertex;
            }
        }

    }

    void GetNextLevelForReOrder(const int &numSortCur, const int *cellIDSortCur, const int *nFaceOfEachCell, const int *cell2Cell, const int *offsetCell2Cell, int *reorderCell, int & numSortNext, int *cellIDSortNext)
    {
        numSortNext = 0;

        for (int iCell = 0; iCell < numSortCur; ++ iCell)
        {
            int cellCellID = cellIDSortCur[iCell];
            int numCellCells = nFaceOfEachCell[cellCellID];
            int offset = offsetCell2Cell[cellCellID];
            for (int iCellCell = 0; iCellCell < numCellCells; ++ iCellCell)
            {
                int cellID = cell2Cell[offset + iCellCell];

                if (reorderCell[cellID] == -1)
                {
                    cellIDSortNext[numSortNext] = cellID;
                    numSortNext ++;
                }
            }
        }
    }
    //! Check rules for cell index
    //! rule1 0 <= newCell < nTotalCell
    //! rule2 cell index should be different
    void CheckNewCell(const int nTotalCell, const int * newCell)
    {
        printf("Check newCell ...\n");

        int *numCells = new int[nTotalCell];

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            int cellID = newCell[iCell];
            if ((cellID < 0) || (cellID > nTotalCell)) printf("CheckNewCell: cellID %d breaks rule1 on %d of newCell\n", cellID, iCell);
            int remain = cellID % nTotalCell;
            numCells[remain] ++;
        }
    
        for (int iCell = 0; iCell < nTotalCell; ++ iCell) 
        {
            if (numCells[iCell] != 1) printf("CheckNewCell, cellID %d in newCell breaks rule 2, has %d cells\n", iCell, numCells[iCell]);
        }

        delete [] numCells;    numCells = NULL;
    }
    //! Compute nFaceOfEachCell and nTotalCellFaces
    //! It should be noted that boundary faces are not considered here.
    void SetFaceNumberOfEachCell(const int nTotalFace, const int nTotalCell, const int *leftCellOfFace, const int *rightCellOfFace, int *nFaceOfEachCell, int &nTotalCellFaces)
    {
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            //! For leftCellOfFace
            int cellLabelLeft = leftCellOfFace[iFace];

            if ((cellLabelLeft < 0)||(cellLabelLeft >= nTotalCell)) 
            {
                printf("Error: cellLabelLeft = %d\n", cellLabelLeft);
                exit(1);
            }
            int cellLabelRight = rightCellOfFace[iFace];
            //! Boundary faces are not considered. In fact, on boundary face, cellLabelLeft is cell itself. cellLabelRight is -1, which is not necessary stored in cell2Cell.
            if (cellLabelRight == -1) continue; //! Boundy face's right neighbor

            if ((cellLabelRight < 0) || (cellLabelRight >= nTotalCell))
            {
                printf("Error: cellLabelRight = %d\n", cellLabelRight);
                exit(1);
            }
            nFaceOfEachCell[cellLabelLeft] ++;
            nTotalCellFaces ++;
            nFaceOfEachCell[cellLabelRight] ++;
            nTotalCellFaces ++;
        }
    }
    //! Set offsetCell2Face by nFaceOfEachCell
    void SetOffSetCell2Face(const int nTotalCell, const int *nFaceOfEachCell, int *offsetCell2Face)
    {
        offsetCell2Face[0] = 0;
        for (int iCell = 1; iCell < nTotalCell; ++ iCell)
        {
            offsetCell2Face[iCell] = offsetCell2Face[iCell-1] + nFaceOfEachCell[iCell-1];
        }
    }
    //! Set cell2Face by leftCellOfFace, rightCellOfFace, offsetCell2Face
    //! It should be noted that boundary faces are not considered here.
    void SetCell2Face(const int nTotalFace, const int nTotalCell, const int *leftCellOfFace, const int *rightCellOfFace, int *offsetCell2Face, int *cell2Face)
    {
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            //! For leftCellOfFace
            int cellLabelLeft = leftCellOfFace[iFace];
            if ((cellLabelLeft < 0)||(cellLabelLeft >= nTotalCell)) 
            {
                printf("Error: cellLabelLeft = %d\n", cellLabelLeft);
                exit(1);
            }

            int cellLabelRight = rightCellOfFace[iFace];
            //! Boundary face is not considered.
            if (cellLabelRight == -1) continue; //! Boundy face's right neighbor

            if ((cellLabelRight < 0) || (cellLabelRight >= nTotalCell))
            {
                printf("Error: cellLabelRight = %d\n", cellLabelRight);
                exit(1);
            }
            int cell2faceOffset = offsetCell2Face[cellLabelLeft];
            cell2Face[cell2faceOffset] = iFace;
            //! OffsetCell2Face update for the next set of face
            offsetCell2Face[cellLabelLeft] ++;
            cell2faceOffset = offsetCell2Face[cellLabelRight];
            cell2Face[cell2faceOffset] = iFace;
            offsetCell2Face[cellLabelRight] ++;
        }
    }

    //! Set cell2Cell
    void SetCell2Cell(const int nTotalCell, const int *nFaceOfEachCell, const int *offsetCell2Face, const int *cell2Face, const int *leftCellOfFace, const int *rightCellOfFace, int *cell2Cell){
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            int sizeCellFace = nFaceOfEachCell[iCell];
            int offsetCellFace = offsetCell2Face[iCell];

            for (int iFace = 0; iFace < sizeCellFace; ++ iFace)
            {
                int faceID = cell2Face[offsetCellFace + iFace];
                int cellLabel = leftCellOfFace[faceID];
                if ((cellLabel < 0)||(cellLabel >= nTotalCell))
                {
                    printf("Error: cellLabel = %d\n", cellLabel);
                    exit(1);
                }

                if (cellLabel != iCell) cell2Cell[offsetCellFace + iFace] = cellLabel;

                cellLabel = rightCellOfFace[faceID];

                if (cellLabel == -1) continue; //! Boundy face's right neighbor
                if ((cellLabel < 0) || (cellLabel >= nTotalCell))
                {
                    printf("Error: cellLabel = %d\n", cellLabel);
                    exit(1);
                }

                if (cellLabel != iCell) cell2Cell[offsetCellFace + iFace] = cellLabel;
                
            }
        }
    }

    //! Compute bandwidth of each cell
    void ComputeBandwidth(const int nTotalCell, const int *nFaceOfEachCell, const int *offsetCell2Face, const int *cell2Cell)
    {
        int *bandwidthOfEachCell = new int[nTotalCell];
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            int sizeCellFace = nFaceOfEachCell[iCell];
            int offsetCellFace = offsetCell2Face[iCell];

            //!The max band between iCell and neighbor cells
            int maxLocalBand = 0;
            for (int iFace = 0; iFace < sizeCellFace; ++ iFace)
            {
                int cellID = cell2Cell[offsetCellFace+iFace];
                int bandDiff = abs(iCell - cellID);
                if (bandDiff > maxLocalBand) maxLocalBand = bandDiff;
            }
            bandwidthOfEachCell[iCell] = maxLocalBand;
        }
        //! Statistics of bandwidthOfEachCell
        int maxGlobalBand = 0; //! Maximum bandwidth
        int maxCellID = 0; //! CellID of maximum bandwidth
        int level = 20;
        int gapLevel = 6000/level;
        printf("gapLevel = %d\n", gapLevel);
        int bandDist[20] = {0}; //! Distribution of bandwidth from 1 to 32 and more
        for (int i = 0; i < level; ++ i)
        {
            bandDist[i] = 0;
        }

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            int bandCell = bandwidthOfEachCell[iCell];
            if (bandCell <= 0) 
            {
                printf("Error: bandwidth<=0 on cell %d, bandwidth=%d\n", iCell, bandCell);
                exit(1);
            }
            if (bandCell > maxGlobalBand)
            {
                maxGlobalBand = bandCell;
                maxCellID = iCell;
            }
            int ratio = static_cast<int>(floor(bandCell/gapLevel));
            if (ratio < 0) 
            {
                printf("Error: ratio = %d, iCell = %d, bandCell = %d, gapLevel=%d\n", ratio, iCell, bandCell, gapLevel);
                exit(0);
            }
            if (ratio > (level-1)) bandDist[level-1] ++;
            else bandDist[ratio] ++;
        }
        
        printf("maxGlobalBand = %d on cellID %d\n", maxGlobalBand, maxCellID);
        int nTotalDist = 0;
        for (int i = 0; i < level; ++ i)
        {
            nTotalDist += bandDist[i];
            printf("bandDist %d: %d in nTotalCells\n", i, bandDist[i]);
        }
        printf("nTotalDist = %d, nTotalCell = %d\n", nTotalDist, nTotalCell);
        
        delete [] bandwidthOfEachCell;   bandwidthOfEachCell = NULL;
    }

    //! Update face_number_of_each_face by swap method
    //! Intput: nTotalCell, face_number_of_each_face, reflectCell
    //! Output:nSwapFaceOfEachCell
    //! The function is used to provide a swap method to update a cell related varible
    //! After cell reorder. In the future, cell2Node can be update by the method
    void SwapFaceNumberOfEachFace(const int nTotalCell, const int *nFaceOfEachCell, const int *reflectCell, const int *nNewFaceOfEachCell)
    {
        //! Try swap method
        int *nSwapFaceOfEachCell = new int[nTotalCell];
        //!It should be notes that orgCell and swapCell are updated during swap
        //! Finally, orgCell should be equal to reflectCell and swapCell should be equal to newCell. 
        int *orgCell = new int[nTotalCell];//! Contianer index and old cellID reflection
        int *swapCell = new int[nTotalCell]; //! Old cellID and container index reflection
        for (int cellID = 0; cellID < nTotalCell; ++ cellID)
        {
            nSwapFaceOfEachCell[cellID] = nFaceOfEachCell[cellID];
            //! Before swap, orgCell is equal to swapCell
            orgCell[cellID] = cellID; 
            swapCell[cellID] = cellID;
        }
        //! After swap, orgCell should be equal to reflectCell
        //!After swap, swapCell should be equal to newCell
        for (int cellID = 0; cellID < nTotalCell; ++ cellID)
        {
            //! In one container, the corresponding old cellID is the same
            //! In orgCell and reflectCell. It means that no swap is required.
            if (orgCell[cellID] == reflectCell[cellID]) continue;
            //! targetID is the target of orgCell[cellID]
            //! swapID is the container index corresponds to targetID
            //! currentID is the real value of orgCell[cellID]
            //! Obviously, orgCell[cellID] and orgCell[swapID] should be swap.
            //! After swap, swapCell should be updated. swapCell[targetID] 
            //! should be cellID. and swapCell[currentID] should be swapID.
            //! In fact, Swap[targerID] and Swap[currentID] is also swap.
            int targetID = reflectCell[cellID];
            int swapID = swapCell[targetID];
            int currentID = orgCell[cellID];
            if (swapCell[currentID] != cellID)
            {
                //! It means that something wrong in  swap of orgCell and swapCell
                printf("Error: Old cellID %d is not in the containerID %d\n", currentID, cellID);
                printf("Old cellID %d is in the containerID %d\n", currentID, swapCell[currentID]);
                printf("Container %d owns old cellID %d\n", cellID, orgCell[cellID]);
                exit(1);
            }
            //! swap nFaceOfEachCell of cellID and swapID
            int temp = nSwapFaceOfEachCell[swapID];
            nSwapFaceOfEachCell[swapID] = nSwapFaceOfEachCell[cellID];
            nSwapFaceOfEachCell[cellID] = temp;
            //! swap orgCell of cellID and swapID
            orgCell[cellID] = targetID;
            orgCell[swapID] = currentID;
            //! swap swapID of targerID and currentID
            swapCell[targetID] = cellID;
            swapCell[currentID] = swapID;
        }
        //! Check results
        CheckSwapFaceNumber(nTotalCell, orgCell, reflectCell, nSwapFaceOfEachCell, nNewFaceOfEachCell);

        
    }

    //! Check results
    void CheckSwapFaceNumber(const int nTotalCell, const int *orgCell, const int *reflectCell, const int *nSwapFaceOfEachCell, const int *nNewFaceOfEachCell)
    {
        printf("Checking swap method for update nFaceOfEachCell just by itself...\n");
        for (int cellID = 0; cellID < nTotalCell; ++ cellID)
        {
            if (orgCell[cellID] != reflectCell[cellID]) 
            {
                printf("Error: cellID %d is not set yet\n", cellID);
                exit(1);
            }
            if (nSwapFaceOfEachCell[cellID]!= nNewFaceOfEachCell[cellID])
            {
                printf("Error: on cellID %d, swap faces %d is not equal  to new faces %d\n", cellID, nSwapFaceOfEachCell[cellID], nNewFaceOfEachCell[cellID]);
                exit(1);
            }
        }
    }
    //! Reorder cell label from nBoundFace to nTotalFace just with the order of faces in cell2Face
    void CellReorderByCell2FaceOrder(const int nTotalCell, const int nTotalFace, const int *offsetCell2Face, const int *nFaceOfEachCell, const int *cell2Face, const int nBoundFace, int labelFace, int *newFaces)
    {
        for (int cellID = 0; cellID < nTotalCell; ++ cellID)
        {
            int offset = offsetCell2Face[cellID];
            int numFaces = nFaceOfEachCell[cellID];
            for (int iFace = 0; iFace < numFaces; ++ iFace)
            {
                int faceID = cell2Face[offset+iFace];
                //! Do not consider boundary faces
                if (faceID < nBoundFace) 
                {
                    printf("Error: faceID %d on cellID %d, iFace %d is boundary face\n", faceID, cellID, iFace);
                    exit(1);
                }
                if (faceID >= nBoundFace)
                {
                    //! faceID does not have label yet.
                    if (newFaces[faceID] == -1) 
                    {
                        labelFace ++;
                        newFaces[faceID] = labelFace;
                    }
                }
            }
        }

        if (labelFace != nTotalFace-1) 
        {
            printf("Error: labelFace %d is not equal to nTotalFace-1 %d\n", labelFace, nTotalFace-1);
            exit(1);
        }
    }
    //! Metric for computing leftCellOfFace and rightCellOfFace for GPU 
    void MetricLeftRightCellOfFace(const int nTotalFace, const int nBoundFace, const int *leftCellOfFace, const int *rightCellOfFace)
    {
        int *diffCount = new int[129]; //! 0~127, 128 for larger
        //! Left part
        //! Start from nBoundFace, boundary faces are not considered
        int segAddress = nBoundFace;
        for (int i = 0; i < 129; ++ i)
        {
            diffCount[i] = 0;
        }

        while (segAddress < nTotalFace)
        {
                int localMax = leftCellOfFace[segAddress];
                int localMin = leftCellOfFace[segAddress];
                for (int offset = 0; offset < 32; ++ offset)
                {
                    int le = leftCellOfFace[segAddress];
                    if (le < localMin) localMin = le;
                    if (le > localMax) localMax = le;
                    segAddress++;
                }
                int diff = localMax - localMin;
                if (diff > 128) diffCount[128] ++;
                else diffCount[diff] ++;
            }
            //! Compute metric elements
            int totalNum = 0;
            for (int i = 0; i < 129; ++ i)
            {
                totalNum += diffCount[i];
            }
            for (int i = 0; i < 129; ++ i)
            {
                double percent = double(diffCount[i])/totalNum;
                printf("i = %d, percent = %f\n", i, percent);
            }

            //! Right part
            segAddress = nBoundFace;

            for (int i = 0; i < 129; ++ i)
            {
                diffCount[i] = 0;
            }
            while (segAddress < nTotalFace)
            {
                int localMax = rightCellOfFace[segAddress];
                int localMin = rightCellOfFace[segAddress];
                for (int offset = 0; offset < 32; offset++)
                {
                    int re = rightCellOfFace[segAddress];
                    if (re < localMin) localMin = re;
                    if (re > localMax) localMax = re;
                    segAddress++;
                }
                int diff = localMax - localMin;
                if (diff > 128) diffCount[128] ++;
                else diffCount[diff] ++;

            }
            //! Test
            totalNum = 0;
            for (int i = 0; i < 129; ++ i)
            {
                totalNum += diffCount[i];
            }
            for (int i = 0; i < 129; ++ i)
            {
                double percent = double(diffCount[i])/totalNum;
                printf("i = %d, percent = %f\n", i, percent);
            }
            delete [] diffCount;   diffCount = NULL;
        
    }
    //!A caller for MetricLeftRightCellOfFace
    void CallMetricLeftRightCellOfFace(Base_Grid_Conn *gconn)
    {
        printf("Program is running in CallMetricLeftRightCellOfFace\n");
        cgsize_t nTotalFace = gconn->GetNTotalFace();
        cgsize_t nBoundFace = gconn->GetNBoundFace();
        vector<cgsize_t> &leftCellOfFace  = *gconn->GetLeftCellOfFace();
        vector<cgsize_t> &rightCellOfFace = *gconn->GetRightCellOfFace();
        int *oldLeftCellOfFace = new int[nTotalFace];
        int *oldRightCellOfFace = new int[nTotalFace];
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            oldLeftCellOfFace[iFace] = static_cast<int>(leftCellOfFace[iFace]);
            oldRightCellOfFace[iFace] = static_cast<int>(rightCellOfFace[iFace]);
        }

        MetricLeftRightCellOfFace(static_cast<int>(nTotalFace), static_cast<int>(nBoundFace), oldLeftCellOfFace, oldRightCellOfFace);

    }

}
