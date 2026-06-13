#include "Pre_CGNSToPHengLEI_Unstruct.h"
#include "Pre_GridBase.h"
#include "Math_BasisFunction.h"
#include "Geo_UnstructBC.h"
#include "PHIO.h"
#include "Glb_Dimension.h"
#include "TK_Exit.h"
using namespace std;
#pragma warning(disable:6385)
#pragma warning(disable:26812)
#pragma warning(disable:4701)

namespace PHSPACE
{
LIB_EXPORT Pre_CGNSToPHengLEI_Unstruct::Pre_CGNSToPHengLEI_Unstruct(const string &gridFileName) : 
    Pre_GridConversion(gridFileName)
{
    factoryCGNS = 0;
}

LIB_EXPORT Pre_CGNSToPHengLEI_Unstruct::~Pre_CGNSToPHengLEI_Unstruct()
{
    delete factoryCGNS;
}

void Pre_CGNSToPHengLEI_Unstruct::ReadGrid()
{
    string &cgnsFile = gridFileName;
    cout << "    CGNS file name: " << cgnsFile << "\n";
    cgsize_t iSize[3][1];
    int fileID, baseID;
    cgsize_t lowerRangID, upperRangID, iStart, iEnd;
    int isBoundary, iParentFlag;
    int cellDim, physDim;
    cgsize_t nTotalNode, nTotalCell;
    int nIFaceRegions;
    int nCoords, nSections, nBCRegions, nFamilies, FamBC, nGeo;
    cgsize_t nBCElem, normListFlag, nIFaceElem;
    int normalIndex, nDataSet;
    char zoneName[33], sectionName[33];
    char baseName[33], coordName[33], bocoName[33], iFaceName[33];
    cgsize_t               elementDataSize;
    ElementType_t          eType;
    ZoneType_t             zoneType;
    DataType_t             dataType;
    BCType_t               bocoType;
    PointSetType_t         pointSetType;
    DataType_t             normDataType;
    CGNSRawCoor            *rawCoor;
    DataType_t             interFaceDataType;
    GridConnectivityType_t iFaceType;

    //! Open the CGNS for reading and check if the file was found.
    if (cg_open(cgnsFile.c_str(), CG_MODE_READ, &fileID) != CG_OK)
    {
        TK_Exit::ExceptionExit(cg_get_error());
    }

    int hasVolumeCondition = 0;
    GlobalDataBase::UpdateData("hasVolumeCondition", &hasVolumeCondition, PHINT, 1);

    //! Determine the of bases in the grid.
    cg_nbases(fileID, &baseID);

    //! Determine the number of FamilyBC in the grid.
    cg_nfamilies(fileID, baseID, &nFamilies);

    typedef char char_33[33];
    char_33 *familyName=NULL;
    char_33 famBCName;
    BCType_t *familyBCType = new BCType_t [nFamilies];
    int bc=1;
    if (nFamilies > 0)
    {
        familyName = new char_33 [nFamilies];
        for (int iFam = 1; iFam <= nFamilies; ++ iFam)
        {
            cg_family_read(fileID, baseID, iFam, familyName[iFam-1], &FamBC, &nGeo);
            cg_fambc_read(fileID, baseID, iFam, bc, famBCName, &familyBCType[iFam-1]);
        }
    }

    //! Read the number of zones in the grid.
    cg_nzones(fileID, baseID, &nBlocks);

    //! Check the cell and physical dimensions of the bases.
    cg_base_read(fileID, baseID, baseName, &cellDim, &physDim);

    if (cellDim != PHSPACE::GetDim())
    {
        TK_Exit::ExceptionExit("Error: Dimension in CGNS file is not according with that in parameter file!");
    }

    factoryCGNS = new CGNSFactory(nBlocks);
    for (int iZone = 1; iZone <= nBlocks; ++ iZone)
    {
        CGNSBase *baseCGNS = factoryCGNS->GetCGNSBase(iZone-1);
        cg_goto(fileID, baseID, "Zone_t", iZone, "end");
        char_33 volumeName;
        cg_famname_read(volumeName);
        baseCGNS->SetVCName(volumeName);
        if ((strcmp(volumeName, "Unspecified") != 0) && (strcmp(volumeName, "") != 0))
        {
            hasVolumeCondition = 1;
            GlobalDataBase::UpdateData("hasVolumeCondition", &hasVolumeCondition, PHINT, 1);
        }

        int vcType = PHENGLEI::UNDEFINED;
        if (nFamilies > 0)
        {
            for (int iFam = 1; iFam <= nFamilies; ++ iFam)
            {
                if (0 == strcmp(volumeName, familyName[iFam-1]))
                {
                    char_33 name;
                    char *desc=NULL;
                    cg_goto(fileID, baseID, "Family_t", iFam, "end");
                    int nDesc;
                    cg_ndescriptors (&nDesc);
                    for (int iDesc = 1; iDesc <= nDesc; ++ iDesc)
                    {
                        cg_descriptor_read (iDesc, name, &desc);
                        if (0 == strcmp(name, "FamVC_TypeName"))
                        {
                            if (0 == strcmp(desc, "Fluid"))
                            {
                                vcType = PHENGLEI::FLUID;
                                break;
                            }
                            else if (0 == strcmp(desc, "Solid"))
                            {
                                vcType = PHENGLEI::SOLID;
                                break;
                            }
                        }
                    }
                    delete [] desc;    desc = nullptr;
                }
            }
        }
        baseCGNS->SetVCType(vcType);

        //! Check the zone type.
        cg_zone_type(fileID, baseID, iZone, &zoneType);

        if (Structured == zoneType)
        {
            TK_Exit::ExceptionExit("Error: the input CGNS mesh must be UNSTRUCTURED type!");
        }

        //! Determine the number of vertices and volume elements in this zone.
        cg_zone_read(fileID, baseID, iZone, zoneName, iSize[0]);
        nTotalNode = iSize[0][0];
        nTotalCell = iSize[1][0];
        baseCGNS->SetNTotalNode(nTotalNode);
        baseCGNS->SetNTotalCell(nTotalCell);
        cout << " nTotalNode = " << nTotalNode << " nTotalCell = " << nTotalCell << "\n";

        //! Determine the number and names of the coordinates.
        cg_ncoords(fileID, baseID, iZone, &nCoords);

        //! Lower range index.
        lowerRangID = 1;
        //! Upper range index of vertices.
        upperRangID = nTotalNode;

        //! Temporarily assume that grid is unstructured, it will be considered later.
        RDouble *x = new RDouble[nTotalNode]();
        RDouble *y = new RDouble[nTotalNode]();
        RDouble *z = new RDouble[nTotalNode]();

        RawGrid *rawGrid = factoryCGNS->GetRawGrid(iZone-1);
        rawGrid->SetNTotalNode(nTotalNode);
        rawGrid->SetX(x);
        rawGrid->SetY(y);
        rawGrid->SetZ(z);

        rawCoor = new CGNSRawCoor();
        for (int iCoor = 0; iCoor < nCoords; ++ iCoor)
        {
            cg_coord_info(fileID, baseID, iZone, iCoor+1, &dataType, coordName);
            rawCoor->AllocateData(iCoor, nTotalNode, dataType);

            //! Read the x-, y-, z-coordinates.
            cg_coord_read(fileID, baseID, iZone, coordName, dataType, &lowerRangID, &upperRangID, rawCoor->GetCoor(iCoor));
        }

        if (2 == nCoords)
        {
            for (int iNode = 0; iNode < nTotalNode; ++ iNode)
            {
                z[iNode] = 0.0;
            }
            rawCoor->AllocateData(2, nTotalNode, dataType);
        }
        rawCoor->SetAllData(x, y, z);
        if (2 == nCoords)
        {
            for (int iNodess = 0; iNodess < nTotalNode; ++ iNodess)
            {
                z[iNodess] = 0.0;
            }
        }
        delete rawCoor;

        //! Determine the number of sections for this zone. Note that
        //! surface elements can be stored in a volume zone, but they
        //! are NOT taken into account in the number obtained from 
        //! cg_zone_read.
        cg_nsections(fileID, baseID, iZone, &nSections);

        //! Loop over the number of sections and read the element
        //! connectivities. As CGNS starts the numbering at 1 the
        //! for-loop starts at 1 as well.
        baseCGNS->CreateElement(nSections);
        cgsize_t *baseCGNSIStart          = baseCGNS->GetIndexOfStart();
        cgsize_t *baseCGNSIEnd            = baseCGNS->GetIndexOfEnd();
        int      *baseCGNSElementType     = baseCGNS->GetElementType();
        cgsize_t **baseCGNSConnList       = baseCGNS->GetElementConnectionList();
        cgsize_t *baseCGNSElementDataSize = baseCGNS->GetElementDataSize();
        cgsize_t **baseCGNSMixedConnList  = baseCGNS->GetMixedConnectionIndex();
        cgsize_t *baseCGNSElementPoint    = factoryCGNS->GetBaseElement()->GetElementPoint();
#ifdef USE_VS2019
        cgsize_t **startOffset            = baseCGNS->GetConnectOffSet();
#endif
        cout << " nSections = " << nSections << "\n";
        for (int iSection = 1; iSection <= nSections; ++ iSection)
        {
            cg_section_read(fileID, baseID, iZone, iSection, sectionName, &eType,
                            &iStart, &iEnd, &isBoundary, &iParentFlag);
            cg_ElementDataSize(fileID, baseID, iZone, iSection, &elementDataSize);

            cout << "\nReading section " << iSection << " ...\n";
            cout << "   section name = " << sectionName << "\n";
            cout << "   section type = " << ElementTypeName[eType] << "\n";
            cout << "   istart,iend =  " << iStart << " " << iEnd << "\n";

            baseCGNSIStart         [iSection - 1] = iStart;
            baseCGNSIEnd           [iSection - 1] = iEnd;
            baseCGNSElementType    [iSection - 1] = eType;
            baseCGNSElementDataSize[iSection - 1] = elementDataSize;

            cgsize_t nElem = iEnd - iStart + 1;
            baseCGNSConnList[iSection - 1] = new cgsize_t [elementDataSize];
#ifdef USE_VS2019
            startOffset     [iSection - 1] = new cgsize_t [nElem + 1];
#endif

            //! Read the connectivity. Again, the node numbering of the 
            //! connectivities start at 1. If internally a starting index 
            //! of 0 is used (typical for C-codes) 1 must be subtracted 
            //! from the connectivities read.
#ifdef USE_VS2019
            if (NGON_n == eType || NFACE_n == eType)
            {
                cg_poly_elements_read(fileID, baseID, iZone, iSection, baseCGNSConnList[iSection - 1], startOffset[iSection - 1], NULL);
            }
            else
            {
                cg_elements_read(fileID, baseID, iZone, iSection, baseCGNSConnList[iSection - 1], NULL);
            }
#else
            cg_elements_read(fileID, baseID, iZone, iSection, baseCGNSConnList[iSection - 1], NULL);
#endif
            //! Determine whether the grid is Cartesian2D or Cartesian3D.
            if (NGON_n == eType)
            {
                if (!(baseCGNS->GetIsCartesian3D()))
                {
                    baseCGNS->SetIsCartesian2D(true);
                }
            }
            else if (NFACE_n == eType)
            {
                baseCGNS->SetIsCartesian2D(false);
                baseCGNS->SetIsCartesian3D(true);
            }

            if (MIXED == eType)
            {
                baseCGNSMixedConnList[iSection - 1] = new cgsize_t [nElem];
                int iNode = 0;
                int index = 0;
                for (; iNode < elementDataSize;)
                {
                    int pointNumber = static_cast<int>(baseCGNSElementPoint[baseCGNSConnList[iSection - 1][iNode]]);
                    baseCGNSMixedConnList[iSection - 1][index] = iNode;
                    iNode = iNode + pointNumber + 1;
                    ++ index;
                }
            }
        }

        //! Determine the number of boundary conditions for this zone.
        cg_nbocos(fileID, baseID, iZone, &nBCRegions);
        cout << "\nNumber of BC Regions: " << nBCRegions << endl;
        baseCGNS->CreateBC(nBCRegions);
        int      *baseCGNSnBCElem              = baseCGNS->GetNumberOfBCElements();
        int      *baseCGNSBoundaryElementType  = baseCGNS->GetBoundaryElementType();
        int      *baseCGNSBoundaryGridLocation = baseCGNS->GetBoundaryGridLocation();
        int      *baseCGNSBCType               = baseCGNS->GetBCType();
        cgsize_t **baseCGNSBCConnList          = baseCGNS->GetBoundaryElementConnectionList();
        GridLocation_t igr;

        //! Loop over the number of boundary conditions.
        string *boundaryName = new string [nBCRegions]();
        baseCGNS->SetBCName(boundaryName);
        for (int iBCRegion = 1; iBCRegion <= nBCRegions; ++ iBCRegion)
        {
            cout << "  Reading bcRegion " << iBCRegion << " ... " << endl;
            cg_goto(fileID, baseID, "Zone_t", iZone, "ZoneBC_t", 1, "BC_t", iBCRegion, "end");
            cg_gridlocation_read(&igr);
            cout << "  Grid Location Name = " << GridLocationName[igr] << "\n";
            if (FaceCenter == igr)
            {
                cout << "GridLocation = FaceCenter means BC data refers to elements, not nodes\n";
            }
            else if (Vertex == igr)
            {
                cout << "GridLocation = Vertex means BC data refers to nodes, not elements\n";
                cout << " ============================================================ " << endl;
                cout << "    Very Important Warning: " << endl;
                cout << "        if the cgns grid file is exported from Gridgen15.18,  " << endl; 
                cout << "        make sure the V16 Format button is toggle off, because" << endl;
                cout << "        the Face-based BCs are exported as Vertex-based type  " << endl;
                cout << "        wrongly by Gridgen15.18 !!!" << endl;
                cout << " =========================================================== "  << endl;
            }
            //! Read the info for this boundary condition.
            cg_boco_info(fileID, baseID, iZone, iBCRegion, bocoName, &bocoType, &pointSetType, &nBCElem, &normalIndex,
                         &normListFlag, &normDataType, &nDataSet);

            if (FamilySpecified == bocoType || nFamilies > 0)
            {
                for (int iFam = 1; iFam <= nFamilies; ++ iFam)
                {
                    if (0 == strcmp(bocoName, familyName[iFam-1]))
                    {
                        bocoType = familyBCType[iFam-1];
                        break;
                    }
                }
            }

            boundaryName[iBCRegion - 1] = bocoName;
            cout << "  Boundary Name = " << bocoName << "\n";
            cout << "  BCtype  = " << BCTypeName[bocoType] << "\n";
            cout << "  nBCElem = " << nBCElem << "\n";

            baseCGNSnBCElem             [iBCRegion-1] = static_cast<int>(nBCElem);
            baseCGNSBCConnList          [iBCRegion-1] = new cgsize_t [nBCElem];
            baseCGNSBCType              [iBCRegion-1] = bocoType;
            baseCGNSBoundaryElementType [iBCRegion-1] = pointSetType;
            baseCGNSBoundaryGridLocation[iBCRegion-1] = igr;

            cout << "  PointSetType_t = " << PointSetTypeName[pointSetType] << "\n";

            //! Read the element ID¡¯s.
            cg_boco_read(fileID, baseID, iZone, iBCRegion, baseCGNSBCConnList[iBCRegion-1], NULL);

            if (2 == nBCElem)
            {
                cout << "  Start & End point of Region " << ": " 
                     << baseCGNSBCConnList[iBCRegion-1][0] << " " << baseCGNSBCConnList[iBCRegion-1][1] << "\n";
            }
            else
            {
                cout << "  The first 20 points of Region " << ": ";
                for (int jPoint = 0; jPoint < MIN(20, (int)nBCElem); ++ jPoint)
                {
                    cout << baseCGNSBCConnList[iBCRegion-1][jPoint] << " ";
                }
                cout << "\n";
            }
            //! And much more to make it fit into the 
            //! internal data structures.

            cout << endl;
        }

        //! Read the interface information of Zones.
        cg_nconns(fileID, baseID, iZone, &nIFaceRegions);
        cout << "Number of InterFace Regions: " << nIFaceRegions << endl;
        baseCGNS->CreateInterFace(nIFaceRegions);
        int      *baseCGNSnIFaceElem        = baseCGNS->GetNumberOfInterFaceElements();
        int      *baseCGNSIFaceElementType  = baseCGNS->GetInterFaceElementType();
        int      *baseCGNSIFaceGridLocation = baseCGNS->GetInterFaceGridLocation();
        int      *baseCGNSIFaceType         = baseCGNS->GetInterFaceType();
        cgsize_t **baseCGNSIFaceConnList    = baseCGNS->GetInterFaceElementConnectionList();
        GridLocation_t igrIFace;

        string *interFaceName = new string [nIFaceRegions]();
        baseCGNS->SetIFaceName(interFaceName);
        for (int iIFaceRegion = 1; iIFaceRegion <= nIFaceRegions; ++ iIFaceRegion)
        {
            cout << "  Reading InterFaceRegion " << iIFaceRegion << " ... " << endl;
            cg_goto(fileID, baseID, "Zone_t", iZone, "ZoneBC_t", 1, "BC_t", iIFaceRegion, "end");
            cg_gridlocation_read(&igrIFace);
            cout << "  Grid Location Name = " << GridLocationName[igr] << "\n";
            if (FaceCenter == igrIFace)
            {
                cout << "GridLocation = FaceCenter means BC data refers to elements, not nodes\n";
            }
            else if (Vertex == igrIFace)
            {
                cout << "GridLocation = Vertex means BC data refers to nodes, not elements\n";
                cout << " ============================================================ " << endl;
                cout << "    Very Important Warning: " << endl;
                cout << "        if the cgns grid file is exported from Gridgen15.18,  " << endl;
                cout << "        make sure the V16 Format button is toggle off, because" << endl;
                cout << "        the Face-based BCs are exported as Vertex-based type  " << endl;
                cout << "        wrongly by Gridgen15.18 !!!" << endl;
                cout << " =========================================================== " << endl;
            }

            char           donorName[33];
            ZoneType_t     donorZoneType;
            PointSetType_t donorPtSetType;
            DataType_t     donorDataType;
            cgsize_t       nDataDonor;

            cg_conn_info(fileID, baseID, iZone, iIFaceRegion, iFaceName,&igrIFace,
                &iFaceType, &pointSetType, &nIFaceElem,donorName, &donorZoneType,
                &donorPtSetType, &donorDataType, &nDataDonor);

            interFaceName[iIFaceRegion - 1] = iFaceName;
            cout << "  InterFace Name = " << iFaceName << "\n";
            cout << "  InterFacetype  = " << GridConnectivityTypeName[iFaceType] << "\n";
            cout << "  nInterFaceElem = " << nIFaceElem << "\n";

            baseCGNSnIFaceElem       [iIFaceRegion - 1] = static_cast<int>(nIFaceElem);
            baseCGNSIFaceConnList    [iIFaceRegion - 1] = new cgsize_t[nIFaceElem];
            baseCGNSIFaceType        [iIFaceRegion - 1] = iFaceType;
            baseCGNSIFaceElementType [iIFaceRegion - 1] = pointSetType;
            baseCGNSIFaceGridLocation[iIFaceRegion - 1] = igrIFace;
            cout << "  PointSetType_t = " << PointSetTypeName[pointSetType] << "\n";
            cgsize_t *donorData = new cgsize_t [nDataDonor];
            cg_conn_read(fileID, baseID, iZone, iIFaceRegion,
                baseCGNSIFaceConnList[iIFaceRegion - 1], donorDataType, donorData);
            delete [] donorData;    donorData = nullptr;

            if (2 == nIFaceElem)
            {
                cout << "  Start & End point of Region " << ": "
                    << baseCGNSIFaceConnList[iIFaceRegion - 1][0] << " " << baseCGNSIFaceConnList[iIFaceRegion - 1][1] << "\n";
            }
            else
            {
                cout << "  The first 20 points of Region " << ": ";
                for (int jPoint = 0; jPoint < MIN(20, (int)nIFaceElem); ++ jPoint)
                {
                    cout << baseCGNSIFaceConnList[iIFaceRegion - 1][jPoint] << " ";
                }
                cout << "\n";
            }
            //! And much more to make it fit into the 
            //! internal data structures.
            cout << endl;
        }
    }
    delete [] familyBCType;    familyBCType = nullptr;

     //! Close CGNS file.
    cg_close(fileID);
}

void Pre_CGNSToPHengLEI_Unstruct::Conversion()
{
    factoryCGNS->ConvertGrid2Fantasy();
    grids = new Grid *[nBlocks];
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        GridID *index = new GridID(iZone);
        Grid *grid = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
        grids[iZone] = grid;
        CGNS2UnsGrid(factoryCGNS, iZone, UnstructGridCast(grid));
        grid->ComputeMinMaxBox();
    }
}

LIB_EXPORT void Pre_CGNSToPHengLEI_Unstruct::WriteAdditionalInformation(const string &targetGridFileName, Grid **grids_in)
{
    WriteFaceBoundaryName(targetGridFileName, grids_in);
    WriteVolumeName(targetGridFileName, grids_in);

    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
    if (dumpOldGrid)
    {
        WriteCellToNode(targetGridFileName, grids_in);

        WriteFaceBC(targetGridFileName, grids_in);
    }
}

}
