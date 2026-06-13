#include "Pre_FieldViewToPHengLEI_Unstruct.h"
#include "Math_BasisFunction.h"
#include "Pre_GridBase.h"
#include "Geo_UnstructBC.h"
#include "PHIO.h"
#include "Glb_Dimension.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
using namespace std;

namespace PHSPACE
{
LIB_EXPORT Pre_FieldViewToPHengLEI_Unstruct::Pre_FieldViewToPHengLEI_Unstruct(const string &gridFileName) : 
    Pre_GridConversion(gridFileName)
{
    factory_cgns = 0;
}

LIB_EXPORT Pre_FieldViewToPHengLEI_Unstruct::~Pre_FieldViewToPHengLEI_Unstruct()
{
    
}

void Pre_FieldViewToPHengLEI_Unstruct::ReadGrid()
{
    cout << "    Field view file name: " << gridFileName << "\n";

    fstream file;
    PHSPACE::OpenFile(file, gridFileName, ios_base::in);

    int nversion, nBCType;
    int nTotalNode, nTotalCell, nBoundFace;

    GetFieldViewVersion(file, nversion);
    GetFieldViewNBlocks(file, nversion);

    GetFieldViewNBCType(file, nversion, nBCType);
    string *bctype_name_list = new string[nBCType];
    int    *bctype_list      = new int   [nBCType];

    GetFieldViewBCNameTable(file, nversion, nBCType, bctype_name_list);
    GetFieldViewBCTypeTable(nBCType, nversion, bctype_name_list, bctype_list);
    Gridgen2CGNSBCTypeTable(nBCType, bctype_list);

    factory_cgns = new CGNSFactory(nBlocks);

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        factory_cgns->GetCGNSBase(iZone)->SetBCName(bctype_name_list);

        cout << " iZone = " << iZone << " nBlocks = " << nBlocks << "\n";
        GetFieldViewNTNode(file, nversion, nTotalNode);

        RDouble *x = new RDouble[nTotalNode];
        RDouble *y = new RDouble[nTotalNode];
        RDouble *z = new RDouble[nTotalNode];

        RawGrid *rawgrid = factory_cgns->GetRawGrid(iZone);

        rawgrid->SetX(x);
        rawgrid->SetY(y);
        rawgrid->SetZ(z);

        ReadFieldViewCoor(file, nTotalNode, x, y, z);
        GetFieldViewNBFace(file, nversion, nBoundFace);

        int nBC_Sections, nBCRegions;

        int *bc_elem_type;
        int *bc_istart, *bc_iend;
        cgsize_t **bc_elem_conn_list;
        cgsize_t **bc_conn_list;
        int *nBCElem;
        int *bc_type;

        GetFieldViewBCCell(file, nBoundFace, bctype_list, 
            nBC_Sections, bc_elem_type, bc_istart, bc_iend, bc_elem_conn_list,
            nBCRegions, bc_conn_list, nBCElem, bc_type);

        int nElem_Sections;

        int *elem_type;
        int *elem_istart, *elem_iend;
        cgsize_t **elem_conn_list;

        GetFieldViewMainCell(file, nTotalCell, nElem_Sections, elem_type, elem_istart, elem_iend, elem_conn_list);

        MergeFieldViewCell(factory_cgns, iZone, nBoundFace, nTotalNode, nTotalCell, nElem_Sections, elem_type, elem_istart, elem_iend, elem_conn_list,
            nBC_Sections, bc_elem_type, bc_istart, bc_iend, bc_elem_conn_list, nBCRegions, bc_conn_list, nBCElem, bc_type);
    }

    delete [] bctype_list;

    PHSPACE::CloseFile(file);
}

LIB_EXPORT void Pre_FieldViewToPHengLEI_Unstruct::WriteAdditionalInformation(const string &targetGridFileName, Grid **grids_in)
{
    //WriteCellToNode(targetGridFileName, grids_in);

    //WriteFaceBC(targetGridFileName, grids_in);

    WriteFaceBoundaryName(targetGridFileName, grids_in);
}

void Pre_FieldViewToPHengLEI_Unstruct::GetFieldViewVersion(fstream &file, int &nversion)
{
    string line,word;
    string separator = " =\t\r\n#$,;";
    string sep1 = "\r\n";
    string errmsg = "The file is not a Field-View file!";
    
    //! Read the format version.
    do
    { 
        getline(file, line);
    }while (line == "");

    nversion = FIELDVIEW_2_4;

    line = FindNextWord(line, word, sep1);
    transform(word.begin(), word.end(), word.begin(), ToLower());
    if (word.substr(0, 9) != "fieldview")
    {
        TK_Exit::ExceptionExit(errmsg);
    }
    else
    {
        string cs;
        string word1;
        int version1, version2;

        if (word.substr(0, 15) == "fieldview_grids")
        {
            //! Format of FIELDVIEW 3.0.
            cs = word.substr(15);
            cs = FindNextWord(cs, word1, separator);
            if (word1 == "")
            {
                TK_Exit::ExceptionExit("Error: main version number is not found!\n");
            }
            from_string< int >(version1, word1, std::dec);
            cs = FindNextWord(cs, word1, separator);
            if (word1 == "")
            {
                TK_Exit::ExceptionExit("Error: sub-version number is not found\n");
            }
            from_string< int >(version2, word1, std::dec);
            cout << "Format of FIELDVIEW: " << version1 << "." << version2 << "\n";
            nversion = FIELDVIEW_3_0;
        }
        else
        {
            //! Format of FIELDVIEW 2.4.
            cs = word.substr(9);
            cs = FindNextWord(cs, word1, separator);
            if (word1 == "")
            {
                TK_Exit::ExceptionExit("Error: main version number is not found!\n");
            }
            from_string< int >(version1, word1, std::dec);
            cs = FindNextWord(cs, word1, separator);
            if (word1 == "")
            {
                TK_Exit::ExceptionExit("Error: sub-version number is not found\n");
            }
            from_string< int >(version2, word1, std::dec);
            cout << "Format of FIELDVIEW: " << version1 << "." << version2 << "\n";
            nversion = FIELDVIEW_2_4;
        }
    }
}

void Pre_FieldViewToPHengLEI_Unstruct::GetFieldViewNBlocks(fstream &file, int nversion)
{
    string line, word;
    string separator = " =\t\r\n#$,;";
    string sep1 = "\r\n";
    do
    {
        getline(file, line);
        line = FindNextWord(line, word, sep1);
        transform(word.begin(), word.end(), word.begin(), ToLower());
    } while (word.substr(0,5) != "grids");

    if (nversion == FIELDVIEW_2_4)
    {
        getline(file, line);
    }
    else
    {
        line = word.substr(5);
    }

    line = FindNextWord(line, word, separator);
    from_string< int >(nBlocks, word, std::dec);
    cout << "Number of Block: " << nBlocks << "\n";
}

void Pre_FieldViewToPHengLEI_Unstruct::GetFieldViewNBCType(fstream &file, int nversion, int &nBCType)
{
    string line, word;
    string separator = " =\t\r\n#$,;";
    string sep1 = "\r\n";

    do
    {
        getline(file, line);
        line = FindNextWord(line, word, sep1);
        transform(word.begin(), word.end(), word.begin(), ToLower());
    } while (word.substr(0, 14) != "boundary table");

    //! Read in the number of BCs.
    if (nversion == FIELDVIEW_2_4)
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        from_string< int >(nBCType, word, std::dec);
    }
    else if (nversion == FIELDVIEW_3_0)
    {
        string word1 = word.substr(14);
        from_string< int >(nBCType, word1, std::dec);
    }
    cout << nBCType << " BCs have been found in this field view file.\n";
}

void Pre_FieldViewToPHengLEI_Unstruct::GetFieldViewBCNameTable(fstream &file, int nversion, int nBCType, string *bctype_name_list)
{
    string line, word;
    string separator = " =\t\r\n#$,;";
    string sep1 = "\r\n";

    int nskipwords = 1;

    if (nversion == FIELDVIEW_3_0) nskipwords = 3;

    //! Analyze the BCs.
    for (int i = 0; i < nBCType; ++ i)
    {
        getline(file, line);
        //! The first number is not necessary, so it was passed.
        for (int m = 0; m < nskipwords; ++ m)
        {
            line = FindNextWord(line, word, separator);
        }

        //! Analyze the remain string.
        line = FindNextWord(line, word, sep1);
        TrimBlanks(word);
        transform(word.begin(), word.end(), word.begin(), ToLower());

        bctype_name_list[i] = word;
        cout << "Boundary Condition: " << i << " = " << word << "\n";
    }
}

void Pre_FieldViewToPHengLEI_Unstruct::GetFieldViewNTNode(fstream &file, int nversion, int &nTotalNode)
{
    string line, word;
    string separator = " =\t\r\n#$,;";
    string sep1 = "\r\n";

    do
    {
        getline(file, line);
        line = FindNextWord(line, word, sep1);
        transform(word.begin(), word.end(), word.begin(), ToLower());
    } while (word.substr(0, 5) != "nodes");

    if (nversion == FIELDVIEW_2_4)
    {
        //! Read the node number.
        getline(file, line);
    }
    else
    {
        line = word.substr(5);
    }

    line = FindNextWord(line, word, separator);
    from_string< int >(nTotalNode, word, std::dec);
    cout << " Number of Nodes = " << nTotalNode << "\n";
}

void Pre_FieldViewToPHengLEI_Unstruct::GetFieldViewNBFace(fstream &file, int nversion, int &nBoundFace)
{
    string line, word;
    string separator = " =\t\r\n#$,;";
    string sep1 = "\r\n";
    string errmsg = "The file is not a Field-View file!";

    getline(file, line);
    line = FindNextWord(line, word, sep1);
    transform(word.begin(), word.end(), word.begin(), ToLower());
    if (word.substr(0, 14) != "boundary faces")
    {
        TK_Exit::ExceptionExit(errmsg);
    }

    if (nversion == FIELDVIEW_2_4)
    {
        //! Read the boundary face number.
        getline(file, line);
    }
    else
    {
        line = word.substr(14);
    }

    line = FindNextWord(line, word, separator);
    from_string< int >(nBoundFace, word, std::dec);
    cout << " Number of Boundary Faces = " << nBoundFace << "\n";
}

void Pre_FieldViewToPHengLEI_Unstruct::GetFieldViewBCTypeTable(int nBCType, int nversion, string *bctype_name_list, int *bctype_list)
{
    using namespace PHSPACE;
    map <string, int> fieldview_bcmap;
    if (nversion == FIELDVIEW_3_0)
    {
        fieldview_bcmap.insert(pair<string, int>("interface"            , GRIDGEN_SPACE::INTERFACE           ));
        fieldview_bcmap.insert(pair<string, int>("no boundary condition", GRIDGEN_SPACE::NO_BOUNDARY_CONDITION));
        fieldview_bcmap.insert(pair<string, int>("solid surface"        , GRIDGEN_SPACE::SOLID_SURFACE       ));
        fieldview_bcmap.insert(pair<string, int>("symmetry"             , GRIDGEN_SPACE::SYMMETRY            ));
        fieldview_bcmap.insert(pair<string, int>("farfield"             , GRIDGEN_SPACE::FARFIELD            ));
        fieldview_bcmap.insert(pair<string, int>("inflow"               , GRIDGEN_SPACE::INFLOW              ));
        fieldview_bcmap.insert(pair<string, int>("outflow"              , GRIDGEN_SPACE::OUTFLOW             ));
        fieldview_bcmap.insert(pair<string, int>("pole"                 , GRIDGEN_SPACE::POLE                ));
        fieldview_bcmap.insert(pair<string, int>("generic #1"           , GRIDGEN_SPACE::GENERIC_1           ));
        fieldview_bcmap.insert(pair<string, int>("generic #2"           , GRIDGEN_SPACE::GENERIC_2           ));
        fieldview_bcmap.insert(pair<string, int>("generic #3"           , GRIDGEN_SPACE::GENERIC_3           ));
    }
    else if (nversion == FIELDVIEW_2_4)
    {
        fieldview_bcmap.insert(pair<string, int>("interface"            , GRIDGEN_SPACE::INTERFACE           ));
        fieldview_bcmap.insert(pair<string, int>("no boundary condition", GRIDGEN_SPACE::NO_BOUNDARY_CONDITION));
        fieldview_bcmap.insert(pair<string, int>("wall"                 , GRIDGEN_SPACE::SOLID_SURFACE       ));
        fieldview_bcmap.insert(pair<string, int>("symmetry plane"       , GRIDGEN_SPACE::SYMMETRY            ));
        fieldview_bcmap.insert(pair<string, int>("farfield"             , GRIDGEN_SPACE::FARFIELD            ));
        fieldview_bcmap.insert(pair<string, int>("inflow"               , GRIDGEN_SPACE::INFLOW              ));
        fieldview_bcmap.insert(pair<string, int>("outflow"              , GRIDGEN_SPACE::OUTFLOW             ));
        fieldview_bcmap.insert(pair<string, int>("pole"                 , GRIDGEN_SPACE::POLE                ));
        fieldview_bcmap.insert(pair<string, int>("generic #1"           , GRIDGEN_SPACE::GENERIC_1           ));
        fieldview_bcmap.insert(pair<string, int>("generic #2"           , GRIDGEN_SPACE::GENERIC_2           ));
        fieldview_bcmap.insert(pair<string, int>("generic #3"           , GRIDGEN_SPACE::GENERIC_3           ));
    }

    int maxNumberOFBC = static_cast<int>(fieldview_bcmap.size());

    for (int i = 0; i < nBCType; ++ i)
    {
        map <string, int> :: const_iterator iter;
        iter = fieldview_bcmap.find(bctype_name_list[i]);
        if (iter == fieldview_bcmap.end())
        {
            // cout << "GetFieldViewBCTypeTable::没有找到对应的边界条件\n";
            // exit(0);
            cout << "Self defined Boundary Condition: " << bctype_name_list[ i ] << endl;
            bctype_list[i] = maxNumberOFBC;
            ++ maxNumberOFBC;
        }
        else
        {
            int bctype = iter->second;
            bctype_list[i] = bctype;
        }
    }
}

void Pre_FieldViewToPHengLEI_Unstruct::Gridgen2CGNSBCTypeTable(int nBCType, int *bctype_list)
{
    using namespace PHSPACE;

    map<int, int> bcmap;
    typedef pair <int, int> Int_Pair;
    map<int, int>::const_iterator iter;

    bcmap.insert(Int_Pair(PHENGLEI::NO_BOUNDARY_CONDITION , BCTypeNull      ));
    bcmap.insert(Int_Pair(PHENGLEI::INTERFACE             , BCTypeUserDefined));
    bcmap.insert(Int_Pair(PHENGLEI::EXTRAPOLATION         , BCExtrapolate   ));
    bcmap.insert(Int_Pair(PHENGLEI::SOLID_SURFACE         , BCWall          ));
    bcmap.insert(Int_Pair(PHENGLEI::SYMMETRY              , BCSymmetryPlane ));
    bcmap.insert(Int_Pair(PHENGLEI::INFLOW                , BCInflow        ));
    bcmap.insert(Int_Pair(PHENGLEI::OUTFLOW               , BCOutflow       ));
    bcmap.insert(Int_Pair(PHENGLEI::POLE                  , BCDegenerateLine));
    bcmap.insert(Int_Pair(PHENGLEI::FARFIELD              , BCFarfield      ));
    bcmap.insert(Int_Pair(PHENGLEI::GENERIC_1             , BCTypeUserDefined));
    bcmap.insert(Int_Pair(PHENGLEI::GENERIC_2             , BCTypeUserDefined));
    bcmap.insert(Int_Pair(PHENGLEI::GENERIC_3             , BCTypeUserDefined));
    bcmap.insert(Int_Pair(PHENGLEI::PERIODIC              , BCTypeUserDefined));

    const int maxNumberOfSelfDefineBC = 20;
    const int maxNumberOFFieldViewBC  = 11;
    const int maxNumberOFCGNSBC       = 26;
    for (int iBC = 0; iBC < maxNumberOfSelfDefineBC; iBC ++)
    {
        bcmap.insert(Int_Pair(maxNumberOFFieldViewBC + iBC, maxNumberOFCGNSBC + iBC));
    }

    for (int i = 0; i < nBCType; ++ i)
    {
        iter = bcmap.find(bctype_list[i]);

        if (iter == bcmap.end())
        {
            TK_Exit::ExceptionExit("Gridgen2CGNSBCTypeTable::No corresponding boundary conditions were found!\n");
        }
        else
        {
            bctype_list[i] = iter->second;
        }
    }
}

void Pre_FieldViewToPHengLEI_Unstruct::ReadFieldViewCoor(fstream &file, int nTotalNode, RDouble *x, RDouble *y, RDouble *z)
{
    string line,word;
    string separator = " =\t\r\n#$,;";

    RDouble xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin =   LARGE;
    xmax = ymax = zmax = - LARGE;

    cout << " Grid coordinates reading ... \n";

    int iskip = MAX(nTotalNode / 10, 1);
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (iNode > 0 && iNode % iskip == 0)
        {
            cout << "     " << (int)((iNode * 1.0) / nTotalNode * 100 + 1) << "% nodes read ..." << endl;
        }
        getline(file, line);
        line = FindNextWord(line, word, separator);
        from_string< RDouble >(x[iNode], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(y[iNode], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(z[iNode], word, std::dec);

        xmin = MIN(xmin, x[iNode]);
        ymin = MIN(ymin, y[iNode]);
        zmin = MIN(zmin, z[iNode]);

        xmax = MAX(xmax, x[iNode]);
        ymax = MAX(ymax, y[iNode]);
        zmax = MAX(zmax, z[iNode]);
    }

    cout << setiosflags(ios::right);
    cout << setprecision(8);
    cout << setiosflags(ios::scientific);
    cout << setiosflags(ios::showpoint);
    int wordwidth = 16;
    cout << " xmin = " << setw(wordwidth) << xmin << " xmax = " << setw(wordwidth) << xmax << "\n";
    cout << " ymin = " << setw(wordwidth) << ymin << " ymax = " << setw(wordwidth) << ymax << "\n";
    cout << " zmin = " << setw(wordwidth) << zmin << " zmax = " << setw(wordwidth) << zmax << "\n";
}

void Pre_FieldViewToPHengLEI_Unstruct::GetFieldViewBCCell(fstream &file, int nBoundFace, int *bctype_list, 
                        int &nBC_Sections, int *&elem_type, int *&istart, int *&iend, cgsize_t **&conn_list,
                        int &nBCRegions, cgsize_t **&bc_conn_list,
                        int *&nBCElem, int *&bc_type)
{
    string line, word;
    string separator = " =\t\r\n#$,;";

    cout << "Boundary elements reading ......\n";

    int *cell_type   = new int [nBoundFace];
    int *cell_bctype = new int [nBoundFace];

    const int max_face_node = 4;

    int *node2Type = new int[max_face_node+1];

    node2Type[3] = TRI_3;
    node2Type[4] = QUAD_4;

    cgsize_t *ref_conn_list = new cgsize_t[nBoundFace * max_face_node];

    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        getline(file, line);

        int bctype_index;
        line = FindNextWord(line, word, separator);
        from_string< int >(bctype_index, word, std::dec);

        cell_bctype[iFace] = bctype_list[bctype_index-1];

        //! The nodes number of each boundary face.
        int nNode;
        cgsize_t nodeIndex;
        line = FindNextWord(line, word, separator);
        from_string< int >(nNode, word, std::dec);

        cell_type[iFace] = node2Type[nNode];

        for (int iNode = 0; iNode < nNode; ++ iNode)
        {
            line = FindNextWord(line, word, separator);
            from_string< cgsize_t >(nodeIndex, word, std::dec);
            ref_conn_list[nodepos+iNode] = nodeIndex;
        }
        nodepos += max_face_node;
    }

    set<int> bctype_temp;
    vector<int> bctype_set;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int type = cell_bctype[iFace];
        set<int> :: const_iterator iter;
        iter = bctype_temp.find(type);
        if (iter == bctype_temp.end())
        {
            bctype_temp.insert(type);
            bctype_set.push_back(type);
        }
    }
    nBCRegions = static_cast<int>(bctype_set.size());

    bc_type = new int [nBCRegions];
    nBCElem = new int [nBCRegions];
    bc_conn_list = new cgsize_t *[nBCRegions];

    vector<int>::iterator iter;
    int icount = 0;
    for (iter = bctype_set.begin(); iter != bctype_set.end(); ++ iter)
    {
        bc_type[icount++] = *iter;
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        icount = 0;
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            if (cell_bctype[iFace] == bc_type[iBCRegion])
            {
                icount ++;
            }
        }
        nBCElem[iBCRegion] = icount;
        bc_conn_list[iBCRegion] = new cgsize_t [nBCElem[iBCRegion]];
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        icount = 0;
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            if (cell_bctype[iFace] == bc_type[iBCRegion])
            {
                bc_conn_list[iBCRegion][icount ++] = iFace;
            }
        }
    }
    int nTri, nQuad, nOthers;
    nTri    = 0;
    nQuad   = 0;
    nOthers = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (cell_type[iFace] == TRI_3)
        {
            nTri ++;
        }
        else if (cell_type[iFace] == QUAD_4)
        {
            nQuad ++;
        }
        else
        {
            nOthers ++;
        }
    }

    if (nOthers > 0)
    {
        TK_Exit::ExceptionExit("Unknown cell type!\n");
        
    }

    nBC_Sections = 0;
    if (nTri  ) nBC_Sections ++;
    if (nQuad ) nBC_Sections ++;
    if (nOthers) nBC_Sections ++;

    int *elem_pt = new int[NofValidElementTypes];

    elem_pt[NODE   ] = 1;
    elem_pt[BAR_2  ] = 2;
    elem_pt[TRI_3  ] = 3;
    elem_pt[QUAD_4 ] = 4;
    elem_pt[TETRA_4] = 4;
    elem_pt[PYRA_5 ] = 5;
    elem_pt[PENTA_6] = 6;
    elem_pt[HEXA_8 ] = 8;

    conn_list   = new cgsize_t * [nBC_Sections];
    elem_type   = new int [nBC_Sections];
    istart      = new int [nBC_Sections];
    iend        = new int [nBC_Sections];
    int *n_elem = new int [nBC_Sections];

    nBC_Sections = 0;
    int ptnum;
    if (nTri)
    {
        elem_type[nBC_Sections] = TRI_3;
        n_elem   [nBC_Sections] = nTri;
        nBC_Sections ++;
    }

    if (nQuad)
    {
        elem_type[nBC_Sections] = QUAD_4;
        n_elem   [nBC_Sections] = nQuad;
        nBC_Sections ++;
    }

    istart[0] = 1;
    for (int m = 0; m < nBC_Sections; ++ m)
    {
        iend[m] = istart[m] + n_elem[m] - 1;
        if (m < nBC_Sections - 1)
        {
            istart[m+1] = iend[m] + 1;
        }
    }
    delete [] n_elem;

    for (int m = 0; m < nBC_Sections; ++ m)
    {
        int etype = elem_type[m];
        ptnum = elem_pt[etype];
        int nelem = iend[m] - istart[m] + 1;
        conn_list[m] = new cgsize_t [nelem*ptnum];
    }

    for (int m = 0; m < nBC_Sections; ++ m)
    {
        nodepos = 0;
        int bc_cell_pos = 0;
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            if (cell_type[iFace] == elem_type[m])
            {
                ptnum = elem_pt[elem_type[m]];
                for (int j = 0; j < ptnum; ++ j)
                {
                    conn_list[m][bc_cell_pos+j] = ref_conn_list[nodepos+j];
                }
                bc_cell_pos += ptnum;
            }
            nodepos += max_face_node;
        }
    }

    delete [] ref_conn_list;

    int *sec_bc  = new int [nBoundFace];
    int *sec_pos = new int [nBoundFace];

    for (int m = 0; m < nBC_Sections; ++ m)
    {
        icount = 0;
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            if (cell_type[iFace] == elem_type[m])
            {
                sec_bc[iFace] = m;
                sec_pos[iFace] = icount;
                icount ++;
            }
        }
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        for (int i = 0; i < nBCElem[iBCRegion]; ++ i)
        {
            uint_t idx = bc_conn_list[iBCRegion][i];
            int m_sec = sec_bc [idx];
            int m_pos = sec_pos[idx];
            bc_conn_list[iBCRegion][i] = istart[m_sec] + m_pos;
        }
    }

    delete [] node2Type;
    delete [] sec_bc;
    delete [] sec_pos;
    delete [] cell_type;
    delete [] cell_bctype;
    delete [] elem_pt;
}

void Pre_FieldViewToPHengLEI_Unstruct::GetFieldViewMainCell(fstream &file, int &nTotalCell, int &nSections, int *&elem_type,
                          int *&istart, int *&iend, cgsize_t **&conn_list)
{
    string line, word;
    string separator = " =\t\r\n#$,;";
    string sep1 = "\r\n";
    string errmsg = "The file is not a Field-View file!";

    getline(file, line);
    line = FindNextWord(line, word, sep1);
    transform(word.begin(), word.end(), word.begin(), ToLower());
    if (word.substr(0, 8) != "elements")
    {
        TK_Exit::ExceptionExit(errmsg);
    }

    cout << "Interior cells start reading ......\n";

    //! Cell type.
    //! FieldView type：
    //! 1. Tetrahedron
    //! 2. Hexahedron
    //! 3. Prism
    //! 4. Pyramid

    int *elem_pt = new int[NofValidElementTypes];

    elem_pt[NODE   ] = 1;
    elem_pt[BAR_2  ] = 2;
    elem_pt[TRI_3  ] = 3;
    elem_pt[QUAD_4 ] = 4;
    elem_pt[TETRA_4] = 4;
    elem_pt[PYRA_5 ] = 5;
    elem_pt[PENTA_6] = 6;
    elem_pt[HEXA_8 ] = 8;

    int *TypeFieldView2CGNS = new int [NofValidElementTypes+1];

    using namespace PHSPACE;

    TypeFieldView2CGNS[FIELDVIEW_SPACE::TRIA] = TRI_3;
    TypeFieldView2CGNS[FIELDVIEW_SPACE::QUAD] = QUAD_4;
    TypeFieldView2CGNS[FIELDVIEW_SPACE::TETR] = TETRA_4;
    TypeFieldView2CGNS[FIELDVIEW_SPACE::HEXA] = HEXA_8;
    TypeFieldView2CGNS[FIELDVIEW_SPACE::PRIS] = PENTA_6;
    TypeFieldView2CGNS[FIELDVIEW_SPACE::PYRA] = PYRA_5;

    vector<int> elem_type_container;
    vector<int> elem_number_container;
    vector<cgsize_t> elem_index_container;

    cgsize_t node_swap[8], node_results[8];
    int total_elem = 0, iskip = 100000;
    while (true)
    {
        streampos file_pos = file.tellp();
        getline(file, line);
        if (file.eof()) break;
        string tmpstring = FindNextWord(line, word, separator);
        transform(word.begin(), word.end(), word.begin(), ToLower());

        if (!isdigit(word[0]))
        {
            file.seekp(file_pos);
            break;
        }

        if (total_elem % iskip == 0)
        {
            cout << total_elem << " elem read ...........................\n";
        }

        //! Read cell type.
        int fieldview_elem_type;
        line = FindNextWord(line, word, separator);
        from_string< int >(fieldview_elem_type, word, std::dec);

        //! Read the cell number of this cell type.
        int nelem_num;
        line = FindNextWord(line, word, separator);
        from_string< int >(nelem_num, word, std::dec);
        cgsize_t nodeindex;

        int etype = TypeFieldView2CGNS[fieldview_elem_type];
        int ptnum = elem_pt[etype];

        elem_type_container.push_back(etype);
        elem_number_container.push_back(nelem_num);

        for (int i = 0; i < nelem_num; ++ i)
        {
            for (int j = 0; j < ptnum; ++ j)
            {
                line = FindNextWord(line, word, separator);
                from_string< cgsize_t >(nodeindex, word, std::dec);
                node_swap[j] = nodeindex;
            }

            ElemReorderFieldView2CGNS(etype, node_swap, node_results);
            for (int j = 0; j < ptnum; ++ j)
            {
                elem_index_container.push_back(node_results[j]);
            }
        }
        total_elem += nelem_num;
    }

    cout << "Interior cells reading over.\n";
    cout << " Number of elements = " << total_elem << "\n";

    set<int> elem_type_set;

    for (std::size_t i = 0; i < elem_type_container.size(); ++ i)
    {
        int etype = elem_type_container[i];
        elem_type_set.insert(etype);
    }

    nSections = static_cast<int>(elem_type_set.size());
    conn_list   = new cgsize_t * [nSections];
    elem_type   = new int [nSections];
    istart      = new int [nSections];
    iend        = new int [nSections];
    int *n_elem = new int [nSections];

    Set2Array(elem_type_set, elem_type);

    for (int m = 0; m < nSections; ++ m)
    {
        int n_total_nelem = 0;
        for (std::size_t i = 0; i < elem_type_container.size(); ++ i)
        {
            int etype     = elem_type_container[i];
            int nelem_num = elem_number_container[i];
            if (etype == elem_type[m])
            {
                n_total_nelem += nelem_num;
            }
        }
        int pt_num   = elem_pt[elem_type[m]];
        conn_list[m] = new cgsize_t[n_total_nelem * pt_num];
        n_elem[m]    = n_total_nelem;
    }

    nTotalCell = 0;
    for (int m = 0; m < nSections; ++ m)
    {
        nTotalCell += n_elem[m];
    }

    istart[0] = 1;
    for (int m = 0; m < nSections; ++ m)
    {
        iend[m] = istart[m] + n_elem[m] - 1;
        if (m < nSections - 1)
        {
            istart[m+1] = iend[m] + 1;
        }
    }

    delete [] n_elem;

    for (int m = 0; m < nSections; ++ m)
    {
        int icount = 0;
        int jcount = 0;
        for (std::size_t i = 0; i < elem_type_container.size(); ++ i)
        {
            int etype     = elem_type_container[i];
            int nelem_num = elem_number_container[i];
            int ptnum     = elem_pt[etype];
            int ptmax     = nelem_num * ptnum;

            if (etype == elem_type[m])
            {
                for (int j = 0; j < ptmax; ++ j)
                {
                    conn_list[m][icount+j] = elem_index_container[jcount+j];
                }
                icount += ptmax;
            }
            jcount += ptmax;
        }
    }

    delete [] TypeFieldView2CGNS;
    delete [] elem_pt;
}

void Pre_FieldViewToPHengLEI_Unstruct::MergeFieldViewCell(CGNSFactory *factory_cgns, int iZone, int nBoundFace, int nTotalNode, int nTotalCell, 
                        int nElem_Sections, int *elem_type, int *elem_istart, int *elem_iend, cgsize_t **elem_conn_list,
                        int nBC_Sections, int *bc_elem_type, int *bc_istart, int *bc_iend, cgsize_t **bc_elem_conn_list,
                        int nBCRegions, cgsize_t **bc_conn_list, int *nBCElem, int *bc_type)
{
    CGNSBase *base_cgns = factory_cgns->GetCGNSBase(iZone);
    Base_Grid_Conn *gconn = factory_cgns->GetBaseGridConnection(iZone);

    cout << "Merging cell..............................\n";
    cout << "nElem_Sections = " << nElem_Sections << " nBC_Sections = " << nBC_Sections << "\n";

    int nSections = nElem_Sections + nBC_Sections;
    base_cgns->SetNBoundFace(nBoundFace);
    base_cgns->SetNTotalNode(nTotalNode);
    base_cgns->SetNTotalCell(nTotalCell);

    gconn->SetNBoundFace(nBoundFace);
    gconn->SetNTotalNode(nTotalNode);
    gconn->SetNTotalCell(nTotalCell);

    base_cgns->CreateElement(nSections);
    base_cgns->SetNumberOfBCRegions(nBCRegions);

    base_cgns->SetBoundaryElementConnectionList(new cgsize_t * [nBCRegions]);
    base_cgns->SetBoundaryElementType(new int [nBCRegions]);
    base_cgns->SetBoundaryGridLocation(new int [nBCRegions]);

    cgsize_t *base_cgns_istart = base_cgns->GetIndexOfStart();
    cgsize_t *base_cgns_iend   = base_cgns->GetIndexOfEnd();
    int *base_cgns_elem_type   = base_cgns->GetElementType();

    cgsize_t **base_cgns_conn_list    = base_cgns->GetElementConnectionList();
    cgsize_t **base_cgns_bc_conn_list = base_cgns->GetBoundaryElementConnectionList();
    cout << "    Warning: the cgsize_t type of conn_list has been convert to int type, when deal with FieldView grid. "
         << "If the number of grid nodes/cells/faces exceed 4-billion in single core, it would be wrong here!" << endl;

    for (int m = 0; m < nElem_Sections; ++ m)
    {
        base_cgns_istart[m] = elem_istart[m];
        base_cgns_iend  [m] = elem_iend  [m];
    }

    for (int m = 0; m < nBC_Sections; ++ m)
    {
        base_cgns_istart[nElem_Sections+m] = nTotalCell + bc_istart[m];
        base_cgns_iend  [nElem_Sections+m] = nTotalCell + bc_iend  [m];
    }

    for (int m = 0; m < nElem_Sections; ++ m)
    {
        base_cgns_conn_list[m] = elem_conn_list[m];
    }

    //! It seems to be something wrong in this place, and needs to be checked. Eric：201101030.
    //for (int m = 0; m < nElem_Sections; ++ m)
    for (int m = 0; m < nBC_Sections; ++ m)
    {
        base_cgns_conn_list[nElem_Sections+m] = bc_elem_conn_list[m];
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        base_cgns_bc_conn_list[iBCRegion] = bc_conn_list[iBCRegion];
        for (int j = 0; j < nBCElem[iBCRegion]; ++ j)
        {
            base_cgns_bc_conn_list[iBCRegion][j] += nTotalCell;
        }
    }

    base_cgns->SetNumberOfBCElements(nBCElem);
    int *base_cgns_bc_elem_set_type = base_cgns->GetBoundaryElementType();
    int *base_cgns_bc_grid_location = base_cgns->GetBoundaryGridLocation();

    base_cgns->SetBCType(bc_type);
    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        base_cgns_bc_elem_set_type[iBCRegion] = ElementList;
        base_cgns_bc_grid_location[iBCRegion] = CellCenter;
    }

    for (int m = 0; m < nElem_Sections; ++ m)
    {
        base_cgns_elem_type[m] = elem_type[m];
    }

    for (int m = 0; m < nBC_Sections; ++ m)
    {
        base_cgns_elem_type[nElem_Sections+m] = bc_elem_type[m];
    }

    delete [] elem_type;
    delete [] bc_elem_type;

    delete [] elem_istart;
    delete [] elem_iend;
    delete [] bc_istart;
    delete [] bc_iend;

    delete [] elem_conn_list;
    delete [] bc_elem_conn_list;
}

void Pre_FieldViewToPHengLEI_Unstruct::ElemReorderFieldView2CGNS(int etype, cgsize_t *node_swap, cgsize_t *node_results)
{
    if (etype == TETRA_4)
    {
        node_results[0] = node_swap[1];
        node_results[1] = node_swap[2];
        node_results[2] = node_swap[3];
        node_results[3] = node_swap[0];
    }
    else if (etype == PYRA_5)
    {
        //! This is actually the fixed here.
        node_results[0] = node_swap[0];
        node_results[1] = node_swap[1];
        node_results[2] = node_swap[2];
        node_results[3] = node_swap[3];
        node_results[4] = node_swap[4];
    }
    else if (etype == PENTA_6)
    {
        node_results[0] = node_swap[0];
        node_results[1] = node_swap[3];
        node_results[2] = node_swap[5];
        node_results[3] = node_swap[1];
        node_results[4] = node_swap[2];
        node_results[5] = node_swap[4];
    }
    else if (etype == HEXA_8)
    {
        node_results[0] = node_swap[0];
        node_results[1] = node_swap[1];
        node_results[2] = node_swap[3];
        node_results[3] = node_swap[2];
        node_results[4] = node_swap[4];
        node_results[5] = node_swap[5];
        node_results[6] = node_swap[7];
        node_results[7] = node_swap[6];
    }
}

void Pre_FieldViewToPHengLEI_Unstruct::Conversion()
{
    factory_cgns->ConvertGrid2Fantasy();

    grids = new Grid * [nBlocks];

    int nDim = PHSPACE::GetDim();
    if (nDim != THREE_D)
    {
        cout << "Error: Only 3D mesh is supported by FiledView !" << endl;
        exit(0);
    }

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        GridID *index = new GridID(iZone);
        Grid *grid = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
        grids[iZone] = grid;
        CGNS2UnsGrid(factory_cgns, iZone, UnstructGridCast(grid));
    }

    delete factory_cgns;
}

}
