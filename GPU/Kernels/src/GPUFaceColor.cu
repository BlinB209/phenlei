#include "GPUFaceColor.h"
using namespace PHSPACE;
using namespace PHMPI;
using namespace std;
namespace GPUFaceColor
{

    int *BoundFaceConflict;
    int *BoundFaceConflictPosi;
    int *BoundFaceConflictNum;
    int *InteriorFaceConflict;
    int *InteriorFaceConflictPosi;
    int *InteriorFaceConflictNum;
    int *BoundFaceGroup;
    int *BoundFaceGroupPosi;
    int *BoundFaceGroupNum;
    int  BoundFaceColorNum = 0;
    int *d_BoundFaceGroup;
    int *InteriorFaceGroup;
    int *InteriorFaceGroupPosi;
    int *InteriorFaceGroupNum;
    int  InteriorFaceColorNum = 0;
    int *d_InteriorFaceGroup;
    void GPUFaceColorMain()
    {
#ifdef CUDAUNITTEST
#ifdef UNITTESTOUTPUT
        cout << "program is in GPUFaceColorMain" << endl;
#endif
#endif
        int nLocalZones = GetNumberofLocalZones();
        for (int izone = 0; izone < nLocalZones; izone++)
        {
            int           zoneID = GetLocalZoneIDToGlobalZoneID(izone);
            int           level  = 0;
            UnstructGrid *grid   = UnstructGridCast(GetGrid(zoneID, level));
            CreateFaceConflict(grid);
            ColorFaces(grid);
            FaceGroupHostToDevice(grid);
        }

        FreeFaceColorVars();
    }

    void CreateFaceConflict(Grid *grid_in)
    {
#ifdef CUDAUNITTEST
#ifdef UNITTESTOUTPUT
        cout << "program is running in CreateFaceConflict" << endl;
#endif
#endif
        UnstructGrid *grid = UnstructGridCast(grid_in);

        const int nBoundFace              = grid->GetNBoundFace();
        const int nTotalFace              = grid->GetNTotalFace();
        int       sumBoundFaceConflict    = 0;
        int       sumInteriorFaceConflict = 0;
        int     **cell2face               = grid->GetCell2Face();

        const int *left_cell_of_face        = grid->GetLeftCellOfFace();
        const int *right_cell_of_face       = grid->GetRightCellOfFace();
        const int *face_number_of_each_cell = grid->GetFaceNumberOfEachCell();

        BoundFaceConflictPosi    = new int[nBoundFace];
        BoundFaceConflictNum     = new int[nBoundFace];
        InteriorFaceConflictPosi = new int[nTotalFace - nBoundFace];
        InteriorFaceConflictNum  = new int[nTotalFace - nBoundFace];
        int le, re;
        for (int i = 0; i < nBoundFace; i++)
        {
            le = left_cell_of_face[i];
            //! re = right_cell_of_face[i];
            //! re is not 0. However, face_number_of_each_cell[re] = 0. Do not consider
            //! about ghost cells.
            BoundFaceConflictPosi[i] = sumBoundFaceConflict;
            //! add the number of face in left and right cell, the total number should
            //! except itself.
            sumBoundFaceConflict += face_number_of_each_cell[le] - 1;
        }

        for (int i = nBoundFace; i < nTotalFace; i++)
        {
            le = left_cell_of_face[i];
            re = right_cell_of_face[i];

            InteriorFaceConflictPosi[i - nBoundFace] = sumInteriorFaceConflict;
            //! add the number of face in left and right cell, the total number should
            //! except itself.
            sumInteriorFaceConflict += face_number_of_each_cell[le] - 1;
            sumInteriorFaceConflict += face_number_of_each_cell[re] - 1;
        }
//! just for test
#ifdef CUDAUNITTEST
#ifdef UNITTESTOUTPUT
        cout << "sumBoundFaceConflict = " << sumBoundFaceConflict << endl;
        cout << "sumInteriorFaceConflict = " << sumInteriorFaceConflict << endl;
#endif
#endif
        BoundFaceConflict    = new int[sumBoundFaceConflict];
        InteriorFaceConflict = new int[sumInteriorFaceConflict];

        for (int i = 0; i < nBoundFace; i++)
        {
            le = left_cell_of_face[i];
            //! Initializaion of BoundFaceConflictNum
            BoundFaceConflictNum[i] = 0;
            for (int j = 0; j < face_number_of_each_cell[le]; j++)
            {
                if (i != cell2face[le][j])
                {
                    BoundFaceConflict[BoundFaceConflictPosi[i] + BoundFaceConflictNum[i]] = cell2face[le][j];
                    BoundFaceConflictNum[i]++;
                }
            }
        }

        for (int i = nBoundFace; i < nTotalFace; i++)
        {
            le            = left_cell_of_face[i];
            re            = right_cell_of_face[i];
            int localFace = i - nBoundFace;
            //! Initializaion of InteriorFaceConflictNum
            InteriorFaceConflictNum[localFace] = 0;
            //!
            for (int j = 0; j < face_number_of_each_cell[le]; j++)
            {
                if (i != cell2face[le][j])
                {
                    InteriorFaceConflict[InteriorFaceConflictPosi[localFace] + InteriorFaceConflictNum[localFace]] =
                        cell2face[le][j];
                    InteriorFaceConflictNum[localFace]++;
                }
            }

            for (int j = 0; j < face_number_of_each_cell[re]; j++)
            {
                if (i != cell2face[re][j])
                {
                    InteriorFaceConflict[InteriorFaceConflictPosi[localFace] + InteriorFaceConflictNum[localFace]] =
                        cell2face[re][j];
                    InteriorFaceConflictNum[localFace]++;
                }
            }
        }
    }

    void ColorFaces(Grid *grid_in)
    {
#ifdef CUDAUNITTEST
#ifdef UNITTESTOUTPUT
        cout << "program is running in ColorFaces" << endl;
#endif
#endif
        UnstructGrid *grid       = UnstructGridCast(grid_in);
        const int     nBoundFace = grid->GetNBoundFace();
        const int     nTotalFace = grid->GetNTotalFace();
        int           colorMax   = 0;
        BoundFaceGroup           = new int[nBoundFace];
        InteriorFaceGroup        = new int[nTotalFace - nBoundFace];
        //! int * faceColor = new int[nBoundFace];
        int *faceColor = new int[nTotalFace];

        //! coloring boundary faces
        //! for (int i = 0; i < nBoundFace; i++) {
        //! firstly, color the boundary faces. nTotalFace should be used, because
        //! boundary face conflict relationship contain interior faces. After face
        //! coloring on boundary, from i=0 to i=nBoundFace-1, faceColor owns different
        //! value. On the contrary, from i=nBoundFace to i=nTotalFace-1, faceColor is
        //! always -1.
        for (int i = 0; i < nTotalFace; i++)
        {
            //! Initialization of faceColor with 0, which means no color
            faceColor[i] = -1;
        }

        for (int i = 0; i < nBoundFace; i++)
        {
            int color     = 0;
            int colorSame = 0;
            while (faceColor[i] == -1)
            {
                for (int j = 0; j < BoundFaceConflictNum[i]; j++)
                {
                    int faceConflict = BoundFaceConflict[BoundFaceConflictPosi[i] + j];
                    if (color == faceColor[faceConflict])
                    {
                        colorSame = 1;
                        break;
                    }
                }
                if (colorSame == 0)
                    faceColor[i] = color;
                else
                {
                    color++;
                    colorSame = 0;
                }
            }
            //! record the maximum color
            if (faceColor[i] > colorMax) colorMax = faceColor[i];
        }
        BoundFaceColorNum = colorMax + 1;
#ifdef CUDAUNITTEST
#ifdef UNITTESTOUTPUT
        printf("Boundary faces own %d colors\n", BoundFaceColorNum);
#endif
#endif
        BoundFaceGroupNum         = new int[BoundFaceColorNum];
        BoundFaceGroupPosi        = new int[BoundFaceColorNum];
        int *BoundFaceColorOffset = new int[BoundFaceColorNum];
        BoundFaceGroupPosi[0]     = 0;
        for (int i = 0; i < BoundFaceColorNum; i++)
        {
            //! Initializaiton with zero
            BoundFaceGroupNum[i]    = 0;
            BoundFaceColorOffset[i] = 0;
        }

        for (int i = 0; i < nBoundFace; i++)
        {
            int color = faceColor[i];
            BoundFaceGroupNum[color]++;
        }

        for (int i = 1; i < BoundFaceColorNum; i++)
        {
            BoundFaceGroupPosi[i] = BoundFaceGroupPosi[i - 1] + BoundFaceGroupNum[i - 1];
        }

        for (int i = 0; i < BoundFaceColorNum; i++)
        {
            //! Initializaiton with zero
            BoundFaceColorOffset[i] = 0;
        }

        for (int i = 0; i < nBoundFace; i++)
        {
            int color                 = faceColor[i];
            int colorPosi             = BoundFaceGroupPosi[color] + BoundFaceColorOffset[color];
            BoundFaceGroup[colorPosi] = i;
            BoundFaceColorOffset[color]++;
        }
        
        //! coloring interior faces
        colorMax = 0; //! reset colorMax for interior faces

        //! nTotalFace is used, because in conflict relationship, nTotalFace is used.
        //! However, only interior faces are colored. Thus, by coloring, from i=0 to
        //! i=nBoundFace-1, faceColor is always -1. From i=nBoundFace to
        //! i=nTotalFace-1, faceColor should be ranged from 0 to InteriorFaceColorNum.
        for (int i = 0; i < nTotalFace; i++)
        {
            //! Initialization of faceColor with 0, which means no color
            faceColor[i] = -1;
        }

        for (int i = nBoundFace; i < nTotalFace; i++)
        {
            int color     = 0;
            int colorSame = 0;
            int localFace = i - nBoundFace;
            while (faceColor[i] == -1)
            {
                for (int j = 0; j < InteriorFaceConflictNum[localFace]; j++)
                {
                    int faceConflict = InteriorFaceConflict[InteriorFaceConflictPosi[localFace] + j];
                    if (color == faceColor[faceConflict])
                    {
                        colorSame = 1;
                        break;
                    }
                }
                if (colorSame == 0)
                    faceColor[i] = color;
                else
                {
                    color++;
                    colorSame = 0;
                }
            }
            //! record the maximum color
            if (faceColor[i] > colorMax) colorMax = faceColor[i];
        }
        InteriorFaceColorNum = colorMax + 1;
#ifdef CUDAUNITTEST
#ifdef UNITTESTOUTPUT
        printf("The interior faces own %d colors\n", InteriorFaceColorNum);
#endif
#endif
        InteriorFaceGroupNum         = new int[InteriorFaceColorNum];
        InteriorFaceGroupPosi        = new int[InteriorFaceColorNum];
        int *InteriorFaceColorOffset = new int[InteriorFaceColorNum];
        InteriorFaceGroupPosi[0]     = 0;
        for (int i = 0; i < InteriorFaceColorNum; i++)
        {
            //! Initializaiton with zero
            InteriorFaceGroupNum[i]    = 0;
            InteriorFaceColorOffset[i] = 0;
        }

        for (int i = nBoundFace; i < nTotalFace; i++)
        {
            int color = faceColor[i];
            InteriorFaceGroupNum[color]++;
        }

        for (int i = 1; i < InteriorFaceColorNum; i++)
        {
            InteriorFaceGroupPosi[i] = InteriorFaceGroupPosi[i - 1] + InteriorFaceGroupNum[i - 1];
        }

        for (int i = 0; i < InteriorFaceColorNum; i++)
        {
            //! Initializaiton with zero
            InteriorFaceColorOffset[i] = 0;
        }

        for (int i = nBoundFace; i < nTotalFace; i++)
        {
            int color                    = faceColor[i];
            int colorPosi                = InteriorFaceGroupPosi[color] + InteriorFaceColorOffset[color];
            InteriorFaceGroup[colorPosi] = i;
            InteriorFaceColorOffset[color]++;
        }
    }

    void FaceGroupHostToDevice(Grid *grid_in)
    {
#ifdef CUDAUNITTEST
#ifdef UNITTESTOUTPUT
        cout << "program is running in FaceGroupHostToDevice" << endl;
#endif
#endif
        UnstructGrid *grid = UnstructGridCast(grid_in);
        //! copy BoundFaceGroup into d_BoundFaceGroup
        const int nBoundFace    = grid->GetNBoundFace();
        size_t    sizeBoundFace = nBoundFace * sizeof(int);
        HANDLE_API_ERR(cudaMalloc((void **)&d_BoundFaceGroup, sizeBoundFace));
        HANDLE_API_ERR(cudaMemcpy(d_BoundFaceGroup, BoundFaceGroup, sizeBoundFace, cudaMemcpyHostToDevice));

        const int nTotalFace       = grid->GetNTotalFace();
        size_t    sizeInteriorFace = (nTotalFace - nBoundFace) * sizeof(int);
        HANDLE_API_ERR(cudaMalloc((void **)&d_InteriorFaceGroup, sizeInteriorFace));
        HANDLE_API_ERR(cudaMemcpy(d_InteriorFaceGroup, InteriorFaceGroup, sizeInteriorFace, cudaMemcpyHostToDevice));
    }

    void FreeFaceColorVars()
    {
#ifdef CUDAUNITTEST
#ifdef UNITTESTOUTPUT
        cout << "program is running in FreeFaceColorVars" << endl;
#endif
#endif
        delete[] BoundFaceConflictPosi;
        delete[] BoundFaceConflict;
        delete[] BoundFaceConflictNum;
        delete[] BoundFaceGroup;
        delete[] InteriorFaceConflictPosi;
        delete[] InteriorFaceConflict;
        delete[] InteriorFaceConflictNum;
        delete[] InteriorFaceGroup;
    }
} //! namespace GPUFaceColor
