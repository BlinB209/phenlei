#include "Pre_FluentToPHengLEI_Unstruct.h"
#include "Geo_UnstructBC.h"
#include "TK_Exit.h"
#include "IO_FileName.h"
#include "Glb_Dimension.h"
#include "Math_BasisFunction.h"
#include "PHIO.h"
#include "Constants.h"
#include "TK_Parse.h"
#include "Pre_GridBase.h"

namespace PHSPACE
{

LIB_EXPORT Pre_FluentToPHengLEI_Unstruct::Pre_FluentToPHengLEI_Unstruct(const string &gridFileName) :
    Pre_GridConversion(gridFileName)
{
    fluent_factory = 0;
}

LIB_EXPORT Pre_FluentToPHengLEI_Unstruct::~Pre_FluentToPHengLEI_Unstruct()
{

}

LIB_EXPORT void Pre_FluentToPHengLEI_Unstruct::WriteAdditionalInformation(const string &out_grid_file)
{
    //WriteFaceBC(out_grid_file, grids);
    WriteFaceBoundaryName(out_grid_file, grids);
}

void Pre_FluentToPHengLEI_Unstruct::ReadGrid()
{
    ReadFluent(fluent_factory, nBlocks, gridFileName);
    fluent_factory->DumpSelfDefinePartInformation();
}

void Pre_FluentToPHengLEI_Unstruct::Conversion()
{
    fluent_factory->DumpSelfDefinePartInformation();

    grids = new Grid * [nBlocks];

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        GridID *index = new GridID(iZone);
        Grid *grid = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
        grids[iZone] = grid;
        fluent_factory->Convert2UnsGrid(UnstructGridCast(grid));
    }
    delete fluent_factory;
}

void Pre_FluentToPHengLEI_Unstruct::ReadFluent(GridFormat_Fluent *&fluent_factory, int &nBlocks, const string &casefile)
{
    fstream file;
    file.open(casefile.c_str(), ios_base::in | ios_base::binary);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(casefile);
    }

    string line;

    nBlocks = 1;
    fluent_factory = new GridFormat_Fluent();

    while (true && !file.eof())
    {
        getline(file, line);
        cout << line << "\n";
        //! The following situation need to be judged.
        if (line.size() == 0) continue;
        if (line[0] == '(')
        {
            AnalysisFluentKeyWord(fluent_factory, file, line);
        }
    }

    file.close();
    file.clear();

    fluent_factory->Process();

    vector< FluentFacePatch * > &face_patchs = fluent_factory->GetFluentFacePatchS();

    cout << "Total number of boundaries: " << face_patchs.size() << "\n";
    for (std::size_t i = 0; i < face_patchs.size(); ++ i)
    {
        if (face_patchs[i]->GetBCType() != 100)
        {
            cout << face_patchs[i]->GetBCType() << " " << face_patchs[i]->GetPatchName() << "\n";
        }
    }
}

void Pre_FluentToPHengLEI_Unstruct::AnalysisFluentKeyWord(GridFormat_Fluent *fluent_factory, fstream &file, string &line)
{
    string word;
    string word0;
    string separator = " \r\n()";
    string sep0 = "\r\n()";

    //! Nodes 
    //! (10 (zone-id first-index last-index type ND))
    //! Cells 
    //! (12 (zone-id first-index last-index type element-type))
    //! Faces 
    //! (13 (zone-id first-index last-index bc-type face-type))
    //! Periodic Shadow Faces
    //! (18 (first-index last-index periodic-zone shadow-zone))
    //! Cell Tree 
    //! (58 (cell-id0 cell-id1 parent-zone-id child-zone-id))
    //! Face Tree
    //! (59 (face-id0 face-id1 parent-zone-id child-zone-id))
    //! Interface Face Parents
    //! (61 (face-id0 face-id1))

    int section_ID, dimension;
    int zone_id, first_index, last_index, type, element_type, bc_type, face_type;

    line = line.substr(1);

    line = FindNextWord(line, word, separator);
    from_string< int >(section_ID, word, std::dec);

    using namespace FLUENT_SPACE;

    if (section_ID == FLUENT_SPACE::FLUENT_COMMENT)
    {
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_HEADER)
    {
        line = FindNextWord(line, word, sep0);
        cout << word << "\n";
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_DIMENSIONS)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(dimension, word, std::hex);
        fluent_factory->SetDimension(dimension);

        if (GetDim() != dimension)
        {
            cout << "Error: Wrong dimensional information !" << endl;
            exit(0);
        }
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_NODES)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(zone_id, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(first_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(last_index, word, std::hex);

        if (zone_id == 0)
        {
            fluent_factory->SetNTotalNode(last_index - first_index + 1);

            int nTotalNode = last_index - first_index + 1;
            RDouble *x = new RDouble[nTotalNode];
            RDouble *y = new RDouble[nTotalNode];
            RDouble *z = new RDouble[nTotalNode];

            RawGrid *rawgrid = fluent_factory->GetRawGrid();
            rawgrid->SetNTotalNode(nTotalNode);
            rawgrid->SetX(x);
            rawgrid->SetY(y);
            rawgrid->SetZ(z);
        }
        else
        {
            int nodeType, nodeDimension;
            line = FindNextWord(line, word, separator);
            from_string< int >(nodeType, word, std::hex);

            line = FindNextWord(line, word, separator);
            from_string< int >(nodeDimension, word, std::hex);

            RawGrid *rawgrid = fluent_factory->GetRawGrid();
            ReadFluentNodes(rawgrid, first_index, last_index, fluent_factory->GetDimension(), file);
        }
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_NODES_BINARY)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(zone_id, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(first_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(last_index, word, std::hex);

        int nodeType, nodeDimension;
        line = FindNextWord(line, word, separator);
        from_string< int >(nodeType, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(nodeDimension, word, std::hex);

        RawGrid *rawgrid = fluent_factory->GetRawGrid();
        ReadFluentNodesBinary(rawgrid, first_index, last_index, fluent_factory->GetDimension(), file);
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_EDGES)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(zone_id, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(first_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(last_index, word, std::hex);

        if (zone_id != 0)
        {
            line = FindNextWord(line, word, separator);
            from_string< int >(type, word, std::hex);

            ReadFluentEdge(first_index, last_index, file);
        }

    }
    else if (section_ID == FLUENT_SPACE::FLUENT_CELLS)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(zone_id, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(first_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(last_index, word, std::hex);

        if (zone_id == 0)
        {
            fluent_factory->SetNTotalCell(last_index - first_index + 1);
        }
        else
        {
            line = FindNextWord(line, word, separator);
            from_string< int >(type, word, std::hex);

            line = FindNextWord(line, word, separator);
            from_string< int >(element_type, word, std::hex);

            if (element_type == FLUENT_SPACE::ELEMENTTYPE_MIXED)
            {
                ReadFluentCell(first_index, last_index, file);
            }
        }
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_CELLS_BINARY)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(zone_id, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(first_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(last_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(type, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(element_type, word, std::hex);

        if (!element_type)
        {
            ReadFluentCellBinary(first_index, last_index, file);
        }
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_FACES)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(zone_id, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(first_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(last_index, word, std::hex);

        if (zone_id == 0)
        {
            fluent_factory->SetNTotalFace(last_index - first_index + 1);
        }
        else
        {
            line = FindNextWord(line, word, separator);
            from_string< int >(bc_type, word, std::hex);

            line = FindNextWord(line, word, separator);
            from_string< int >(face_type, word, std::hex);

            ReadFluentFace(fluent_factory, zone_id, first_index, last_index, bc_type, face_type, file);
        }
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_FACES_BINARY)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(zone_id, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(first_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(last_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(bc_type, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(face_type, word, std::hex);

        ReadFluentFaceBinary(fluent_factory, zone_id, first_index, last_index, bc_type, face_type, file);
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_PERIODICFACES)
    {
        int periodicFaceZone, shadowFaceZone;

        line = FindNextWord(line, word, separator);
        from_string< int >(first_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(last_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(periodicFaceZone, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(shadowFaceZone, word, std::hex);

        ReadFluentPeriodicface(first_index, last_index, file);
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_C45)
    {
        line = FindNextWord(line, word0, separator);
        from_string< int >(zone_id, word0, std::dec);
        line = FindNextWord(line, word, separator);

        line = FindNextWord(line, word, separator);

        word0.insert(0, "-");
        string::size_type loc = word.find(word0);

        if (loc != string::npos)
        {
            word.erase(loc);
        }

        string selfDefineBCName = word;

        fluent_factory->SetBoundaryName(zone_id, selfDefineBCName);
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_C39)
    {
        line = FindNextWord(line, word0, separator);
        from_string< int >(zone_id, word0, std::dec);
        line = FindNextWord(line, word, separator);

        line = FindNextWord(line, word, separator);

        word0.insert(0, "-");
        string::size_type loc = word.find(word0);

        if (loc != string::npos)
        {
            word.erase(loc);
        }

        string selfDefineBCName = word;

        fluent_factory->SetBoundaryName(zone_id, selfDefineBCName);

        int frontBracket = count(line.begin(), line.end(), '(') + 2;
        int closedBracket = count(line.begin(), line.end(), ')');
        frontBracket = frontBracket - closedBracket;

        while (frontBracket)
        {
            getline(file, line);

            int frontBracket1 = count(line.begin(), line.end(), '(');
            int closedBracket1 = count(line.begin(), line.end(), ')');

            frontBracket = frontBracket + frontBracket1 - closedBracket1;
        }
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_CELLTREE)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(first_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(last_index, word, std::hex);

        int parentZoneId, childZoneId;
        line = FindNextWord(line, word, separator);
        from_string< int >(parentZoneId, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(childZoneId, word, std::hex);

        ReadFluentCellTree(first_index, last_index, file);
    }
    else if (section_ID == FLUENT_SPACE::FLUENT_FACETREE)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(first_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(last_index, word, std::hex);

        int parentZoneId, childZoneId;
        line = FindNextWord(line, word, separator);
        from_string< int >(parentZoneId, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(childZoneId, word, std::hex);

        ReadFluentFaceTree(first_index, last_index, file);
    }
    else if (section_ID == 3041)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(first_index, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(last_index, word, std::hex);

        getline(file, line);
    }
    else
    {
        int frontBracket = count(line.begin(), line.end(), '(') + 1;
        int closedBracket = count(line.begin(), line.end(), ')');
        frontBracket = frontBracket - closedBracket;

        while (frontBracket)
        {
            getline(file, line);

            int frontBracket1  = count(line.begin(), line.end(), '(');
            int closedBracket1 = count(line.begin(), line.end(), ')');

            frontBracket = frontBracket + frontBracket1 - closedBracket1;
        }
    }
}

void Pre_FluentToPHengLEI_Unstruct::ReadFluentNodes(RawGrid *rawgrid, int first_index, int last_index, int dimension, fstream &file)
{
    RDouble *x = rawgrid->GetX();
    RDouble *y = rawgrid->GetY();
    RDouble *z = rawgrid->GetZ();

    string line, word;
    string separator = " \r\n()";

    int i = first_index;
    do
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        if (word == "") continue;

        from_string< RDouble >(x[i-1], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(y[i-1], word, std::dec);

        if (dimension == TWO_D)
        {
            z[i-1] = 0.0;
        }
        else
        {
            line = FindNextWord(line, word, separator);
            from_string< RDouble >(z[i-1], word, std::dec);
        }

        ++ i;
    } while (i <= last_index);
}

void Pre_FluentToPHengLEI_Unstruct::ReadFluentNodesBinary(RawGrid *rawgrid, int first_index, int last_index, int dimension, fstream &file)
{
    RDouble *x = rawgrid->GetX();
    RDouble *y = rawgrid->GetY();
    RDouble *z = rawgrid->GetZ();

    char cbuf;
    file.read(&cbuf, 1);

    int iNode = first_index;
    do
    {
        double xx, yy, zz;
        file.read(reinterpret_cast<char *>(&xx), sizeof(double));
        file.read(reinterpret_cast<char *>(&yy), sizeof(double));
        if (dimension == TWO_D)
        {
            zz = 0.0;
        }
        else
        {
            file.read(reinterpret_cast<char *>(&zz), sizeof(double));
        }

        x[iNode - 1] = xx;
        y[iNode - 1] = yy;
        z[iNode - 1] = zz;

        ++ iNode;
    } while (iNode <= last_index);
}

void Pre_FluentToPHengLEI_Unstruct::ReadFluentFace(GridFormat_Fluent *fluent_factory, int zone_id, int first_index, int last_index, int bc_type, int face_type, fstream &file)
{
    int nTotalFace = last_index - first_index + 1;

    string line, word;
    string separator = " \r\n()";

    //! A non-zero zone-id indicates a regular face section, and will be followed by a body that contains information about the grid connectivity.
    //! Each line of the body will describe one face and will have the following format: 
    //! n0 n1 n2 c0 c1
    //! where, n* = defining nodes (vertices) of the face 
    //! c* = adjacent cells 
    //! This is the format for a 3D grid with a triangular face format.
    //! The actual number of nodes depends on the face-type.
    //! The order of the cell indices is important, and is determined by the right-hand rule: 
    //! if you curl the fingers of your right hand in the order of the nodes, your thumb will point toward c1.

    //! If the face zone is of mixed type (face-type = 0), 
    //! each line of the section body will begin with a reference to the number of nodes that make up that particular face,
    //! and has the following format: 
    //! x n0 n1 ... nf c0 c1
    //! where, 
    //! x = the number of nodes (vertices) of the face 
    //! nf = the final node of the face.

    vector< FluentFacePatch * > &face_patchs = fluent_factory->GetFluentFacePatchS();
    FluentFacePatch *face_patch = new FluentFacePatch();
    face_patchs.push_back(face_patch);

    face_patch->SetZoneID(zone_id);
    face_patch->SetBCType(bc_type);
    face_patch->SetFaceType(face_type);
    face_patch->SetNTotalFace(nTotalFace);

    vector<int> face_node;
    vector< vector<int> > &face2node = face_patch->GetFace2Node();
    vector<int> &left_cell_of_face  = face_patch->GetLeftCellOfFace();
    vector<int> &right_cell_of_face = face_patch->GetRightCellOfFace();

    int number, nodeid, le ,re;

    if (face_type == FLUENT_SPACE::FACETYPE_MIXED || face_type == FLUENT_SPACE::FACETYPE_POLYGONAL)
    {
        int i = first_index;
        do
        {
            getline(file, line);
            line = FindNextWord(line, word, separator);
            if (word == "") continue;

            from_string< int >(number, word, std::hex);
            face_node.resize(number);

            for (int j = 1; j <= number; ++ j)
            {
                line = FindNextWord(line, word, separator);
                from_string< int >(nodeid, word, std::hex);
                face_node[j-1] = nodeid;
            }
            face2node.push_back(face_node);

            line = FindNextWord(line, word, separator);
            from_string< int >(le, word, std::hex);
            line = FindNextWord(line, word, separator);
            from_string< int >(re, word, std::hex);

            if (bc_type == FLUENT_SPACE::INTERIOR  && fluent_factory->GetDimension() == THREE_D)
            {
                SWAP(le, re);
            }

            left_cell_of_face.push_back(le);
            right_cell_of_face.push_back(re);

            ++ i;
        } while (i <= last_index);
    }
    else
    {
        int i = first_index;
        do
        {
            getline(file, line);
            line = FindNextWord(line, word, separator);
            if (word == "") continue;

            face_node.resize(face_type);

            for (int j = 1; j <= face_type; ++ j)
            {
                from_string< int >(nodeid, word, std::hex);
                face_node[j-1] = nodeid;
                line = FindNextWord(line, word, separator);
            }
            face2node.push_back(face_node);

            from_string< int >(le, word, std::hex);
            line = FindNextWord(line, word, separator);
            from_string< int >(re, word, std::hex);

            if (bc_type == FLUENT_SPACE::INTERIOR  && fluent_factory->GetDimension() == THREE_D)
            {
                SWAP(le, re);
            }

            left_cell_of_face.push_back(le);
            right_cell_of_face.push_back(re);
            ++ i;
        } while (i <= last_index);
    }
    //delete face_patch;
}

void Pre_FluentToPHengLEI_Unstruct::ReadFluentFaceBinary(GridFormat_Fluent *fluent_factory, int zone_id, int first_index, int last_index, int bc_type, int face_type, fstream &file)
{
    char cbuf;
    file.read(&cbuf, 1);

    int nTotalFace = last_index - first_index + 1;

    FluentFacePatch *face_patch = new FluentFacePatch();
    vector< FluentFacePatch* > &face_patchs = fluent_factory->GetFluentFacePatchS();

    face_patchs.push_back(face_patch);
    face_patch->SetZoneID(zone_id);
    face_patch->SetBCType(bc_type);
    face_patch->SetFaceType(face_type);
    face_patch->SetNTotalFace(nTotalFace);

    vector< vector<int> > &face2node = face_patch->GetFace2Node();
    vector<int> &leftCellOfFace = face_patch->GetLeftCellOfFace();
    vector<int> &rightCellOfFace = face_patch->GetRightCellOfFace();

    vector<int> face_node;
    int nodeNumber, nodeId, leftCell, rightCell;
    if (face_type == 0 || face_type == 5)
    {
        int iFace = first_index;
        do
        {
            file.read(reinterpret_cast<char *>(&nodeNumber), sizeof(int));
            face_node.resize(nodeNumber);

            for (int iNode = 1; iNode <= nodeNumber; ++ iNode)
            {
                file.read(reinterpret_cast<char *>(&nodeId), sizeof(int));
                face_node[iNode - 1] = nodeId;
            }
            face2node.push_back(face_node);

            file.read(reinterpret_cast<char *>(&leftCell), sizeof(int));
            file.read(reinterpret_cast<char *>(&rightCell), sizeof(int));

            if (bc_type == FLUENT_SPACE::INTERIOR && fluent_factory->GetDimension() == THREE_D)
            {
                SWAP(leftCell, rightCell);
            }
            leftCellOfFace.push_back(leftCell);
            rightCellOfFace.push_back(rightCell);

            ++ iFace;
        } while (iFace <= last_index);
    }
    else
    {
        int iFace = first_index;
        do
        {
            face_node.resize(face_type);

            for (int iNode = 1; iNode <= face_type; ++ iNode)
            {
                file.read(reinterpret_cast<char *>(&nodeId), sizeof(int));
                face_node[iNode - 1] = nodeId;
            }
            face2node.push_back(face_node);

            file.read(reinterpret_cast<char *>(&leftCell), sizeof(int));
            file.read(reinterpret_cast<char *>(&rightCell), sizeof(int));

            if (bc_type == FLUENT_SPACE::INTERIOR && fluent_factory->GetDimension() == THREE_D)
            {
                SWAP(leftCell, rightCell);
            }
            leftCellOfFace.push_back(leftCell);
            rightCellOfFace.push_back(rightCell);

            ++ iFace;
        } while (iFace <= last_index);
    }
}

void Pre_FluentToPHengLEI_Unstruct::ReadFluentEdge(int first_index, int last_index, fstream &file)
{
    string line, word;
    string separator = " \r\n()";

    int node1, node2;
    int i = first_index;
    do
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        if (word == "") continue;

        from_string< int >(node1, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(node2, word, std::hex);

        ++ i;
    } while (i <= last_index);
}

void Pre_FluentToPHengLEI_Unstruct::ReadFluentPeriodicface(int first_index, int last_index, fstream &file)
{
    string line, word;
    string separator = " \r\n()";

    int periodicFaceIndex, shadowFaceIndex;
    int i = first_index;
    do
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        if (word == "") continue;

        from_string< int >(periodicFaceIndex, word, std::hex);

        line = FindNextWord(line, word, separator);
        from_string< int >(shadowFaceIndex, word, std::hex);

        ++ i;
    }while (i <= last_index);
}

void Pre_FluentToPHengLEI_Unstruct::ReadFluentCellTree(int first_index, int last_index, fstream &file)
{
    string line, word;
    string separator = " \r\n()";

    int numberOfKids, kidIndex;
    int i = first_index;
    do
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        if (word == "") continue;

        from_string< int >(numberOfKids, word, std::hex);

        for (int iKid = 0; iKid < numberOfKids; ++ iKid)
        {
            line = FindNextWord(line, word, separator);
            from_string< int >(kidIndex, word, std::hex);
        }

        ++ i;
    }while (i <= last_index);
}

void Pre_FluentToPHengLEI_Unstruct::ReadFluentFaceTree(int first_index, int last_index, fstream &file)
{
    string line, word;
    string separator = " \r\n()";

    int numberOfKids, kidIndex;
    int i = first_index;
    do
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        if (word == "") continue;

        from_string< int >(numberOfKids, word, std::hex);

        for (int iKid = 0; iKid < numberOfKids; ++ iKid)
        {
            line = FindNextWord(line, word, separator);
            from_string< int >(kidIndex, word, std::hex);
        }

        ++ i;
    }while (i <= last_index);
}

void Pre_FluentToPHengLEI_Unstruct::ReadFluentCell(int first_index, int last_index, fstream &file)
{
    string line, word;
    string separator = " \r\n()";

    int type;

    int i = first_index;
    do
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        if (word == "") continue;

        do
        {
            from_string< int >(type, word, std::hex);
            line = FindNextWord(line, word, separator);

            ++ i;
            if (word == "") break;
        } while (true);
    } while (i <= last_index);
}

void Pre_FluentToPHengLEI_Unstruct::ReadFluentCellBinary(int first_index, int last_index, fstream &file)
{
    char cbuf;
    file.read(&cbuf, 1);

    int cellType;
    int iCell = first_index;
    do
    {
        file.read(reinterpret_cast<char *>(&cellType), sizeof(int));
        ++ iCell;
    } while (iCell <= last_index);
}

GridFormat_Fluent::GridFormat_Fluent()
{
    rawgrid = new RawGrid();
    BuildBCMap();
}

GridFormat_Fluent::~GridFormat_Fluent()
{
    delete rawgrid;
}

void GridFormat_Fluent::BuildBCMap()
{
    typedef pair <int, int> Int_Pair;
    using namespace PHSPACE;
    bcmap.insert(Int_Pair(FLUENT_SPACE::INTERIOR          , PHENGLEI::INTERIOR            ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::WALL              , PHENGLEI::SOLID_SURFACE       ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::PRESSURE_INLET    , PHENGLEI::INFLOW              ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::INLET_VENT        , PHENGLEI::INFLOW              ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::INTAKE_FAN        , PHENGLEI::INFLOW              ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::VELOCITY_INLET    , PHENGLEI::INFLOW));
    bcmap.insert(Int_Pair(FLUENT_SPACE::PRESSURE_OUTLET   , PHENGLEI::OUTFLOW             ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::EXHAUST_FAN       , PHENGLEI::OUTFLOW             ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::OUTLET_VENT       , PHENGLEI::OUTFLOW             ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::SYMMETRY          , PHENGLEI::SYMMETRY            ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::PERIODIC_SHADOW   , PHENGLEI::PERIODIC            ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::PRESSURE_FAR_FIELD, PHENGLEI::FARFIELD            ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::PERIODIC          , PHENGLEI::PERIODIC            ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::FAN               , PHENGLEI::NO_BOUNDARY_CONDITION));
    bcmap.insert(Int_Pair(FLUENT_SPACE::POROUS_JUMP       , PHENGLEI::NO_BOUNDARY_CONDITION));
    bcmap.insert(Int_Pair(FLUENT_SPACE::RADIATOR          , PHENGLEI::NO_BOUNDARY_CONDITION));
    bcmap.insert(Int_Pair(FLUENT_SPACE::MASS_FLOW_INLET   , PHENGLEI::INFLOW              ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::INTERFACE         , PHENGLEI::OVERSET             ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::OUTFLOW           , PHENGLEI::OUTFLOW             ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::AXIS              , PHENGLEI::POLE                ));
    bcmap.insert(Int_Pair(FLUENT_SPACE::UNSPECIFIED       , PHENGLEI::GENERIC_1           ));

    int maxNumberOfBody = PHENGLEI::ENDNumberOfBodyOfHyperFLOW - PHENGLEI::StartNumberOfBodyOfHyperFLOW + 1;
    for (int iBody = 0; iBody < maxNumberOfBody; iBody ++)
    {
        bcmap.insert(Int_Pair(FLUENT_SPACE::MaxFluentBCType + iBody, PHENGLEI::StartNumberOfBodyOfHyperFLOW + iBody));
    }
}

int GridFormat_Fluent::Fluent2FantasyBC(int bctype)
{
    map<int, int>::const_iterator iter;

    iter = bcmap.find(bctype);
    if (iter == bcmap.end())
    {
        TK_Exit::ExceptionExit("Error: couldn't map fluent BC to HyperFLOW BC!\n");
    }
    return iter->second;
}

void GridFormat_Fluent::Process()
{
    ReorderFacePatchs();

    nBoundFace = 0;
    for (std::size_t m = 0; m < face_patchs.size(); ++ m)
    {
        int bc_type = face_patchs[m]->GetBCType();
        if (bc_type != FLUENT_SPACE::INTERIOR)
        {
            nBoundFace += face_patchs[m]->GetNTotalFace();
        }
    }

    for (std::size_t m = 0; m < face_patchs.size(); ++ m)
    {
        int bc_type = face_patchs[m]->GetBCType();
        int nPFace  = face_patchs[m]->GetNTotalFace();

        vector< vector<int> > &face2node = face_patchs[m]->GetFace2Node();

        vector<int> &left_cell_of_face  = face_patchs[m]->GetLeftCellOfFace();
        vector<int> &right_cell_of_face = face_patchs[m]->GetRightCellOfFace();

        for (int i = 0; i < nPFace; ++ i)
        {
            left_cell_of_face[i]  -= 1;
            right_cell_of_face[i] -= 1;

            uint_t npt = face2node[i].size();
            for (int j = 0; j < npt; ++ j)
            {
                face2node[i][j] -= 1;
            }
        }

        if (bc_type != FLUENT_SPACE::INTERIOR && this->dimension == THREE_D)
        {
            for (int i = 0; i < nPFace; ++ i)
            {
                std::reverse(face2node[i].begin(), face2node[i].end());
            }
        }
    }

    for (std::size_t m = 0; m < face_patchs.size(); ++ m)
    {
        face_patchs[m]->SetBCType(Fluent2FantasyBC(face_patchs[m]->GetBCType()));
    }
}

void GridFormat_Fluent::ReorderFacePatchs()
{
    nPatchs = static_cast<int>(face_patchs.size());

    int firstp = 0;
    int lastp = nPatchs - 1;

    while (firstp < lastp)
    {
        int bc_type;
        bool flag = false;
        for (int m = lastp; m >= 0; -- m)
        {
            bc_type = face_patchs[m]->GetBCType();
            if (bc_type != FLUENT_SPACE::INTERIOR)
            {
                lastp = m;
                flag  = true;
                break;
            }
        }

        if (!flag)
        {
            break;
        }

        flag = false;
        for (int m = firstp; m < lastp; ++ m)
        {
            bc_type = face_patchs[m]->GetBCType();
            if (bc_type == FLUENT_SPACE::INTERIOR)
            {
                firstp = m;
                flag = true;
                break;
            }
        }

        if (!flag)
        {
            break;
        }

        SWAP(face_patchs[firstp], face_patchs[lastp]);
    }
}

void GridFormat_Fluent::SetBoundaryName(int zoneIndex, const string &BCName)
{
    vector< FluentFacePatch * > &face_patchs = this->GetFluentFacePatchS();

    for (std::size_t iPatch = 0; iPatch < face_patchs.size(); ++ iPatch)
    {
        int zoneIndexOfFacePatch = face_patchs[iPatch]->GetZoneID();
        if (zoneIndexOfFacePatch == zoneIndex)
        {
            face_patchs[iPatch]->SetPatchName(BCName);
            return;
        }
    }
}

void GridFormat_Fluent::ResetBCTypeOfSelfDefinePart(int zoneIndex)
{
    int MaxBCType = -1;
    vector< FluentFacePatch * > &face_patchs = this->GetFluentFacePatchS();
    for (std::size_t iPatch = 0; iPatch < face_patchs.size(); ++ iPatch)
    {
        int bcType = face_patchs[iPatch]->GetBCType();
        MaxBCType = MAX(MaxBCType, bcType);
    }

    if (MaxBCType < FLUENT_SPACE::MaxFluentBCType)
    {
        MaxBCType = FLUENT_SPACE::MaxFluentBCType;
    }
    else
    {
        ++ MaxBCType;
    }

    for (std::size_t iPatch = 0; iPatch < face_patchs.size(); ++ iPatch)
    {
        int zoneIndexOfFacePatch = face_patchs[iPatch]->GetZoneID();
        if (zoneIndexOfFacePatch == zoneIndex)
        {
            face_patchs[iPatch]->SetBCType(MaxBCType);
        }
    }
}

void GridFormat_Fluent::DumpSelfDefinePartInformation()
{
    //! Collect all of the boundary condition.
    set< string > nameContainer;
    vector< int > bcType;
    vector< string > boundaryName;

    vector< FluentFacePatch * > &face_patchs = this->GetFluentFacePatchS();
    for (std::size_t iPatch = 0; iPatch < face_patchs.size(); ++ iPatch)
    {
        int fantasyBCType = face_patchs[iPatch]->GetBCType();

        string bcName = face_patchs[iPatch]->GetPatchName();

        if (fantasyBCType >= PHENGLEI::StartNumberOfBodyOfHyperFLOW && fantasyBCType <= PHENGLEI::ENDNumberOfBodyOfHyperFLOW)
        {
            set< string > :: iterator nameIter;
            nameIter = nameContainer.find(bcName);
            if (nameIter == nameContainer.end())
            {
                boundaryName.push_back(bcName);
                bcType.push_back(fantasyBCType);
            }
        }
    }

    uint_t nPart = bcType.size();
    if (!nPart)
    {
        return;
    }

    string fileName;
    GlobalDataBase::GetData("from_gfile", &fileName, PHSTRING, 1);
    fileName = ChangeExtensionOfFileName(fileName, "part");

    fstream outfile(fileName.c_str(), ios::out);
    if (!outfile)
    {
        ostringstream oss;
        oss << "Can not open file " << fileName << endl;
        TK_Exit::ExceptionExit(oss.str());
    }

    RDouble dm = 0.0, dl = 0.0, dn3 = 0.0;
    RFloat Lref = 1.0, Sref = 1.0, xref = 0.0, yref = 0.0, zref = 0.0;

    outfile << nPart << endl;
    outfile << dm << " " << dl << " " << dn3 << endl;

    vector< int > :: iterator typeIter = bcType.begin();
    vector< string > :: iterator nameIter = boundaryName.begin();

    for (int iPart = 0; iPart < nPart; ++ iPart)
    {
        int type = *typeIter;
        string name = *nameIter;

        outfile << name << "	" << type << endl;
        outfile << "		" << xref << "	" << yref << "	" << zref << "	" << Sref << "	" << Lref << endl;

        ++ typeIter;
        ++ nameIter;
    }

    outfile.close();
    outfile.clear();
}

void GridFormat_Fluent::Convert2UnsGrid(UnstructGrid *grid)
{
    int fcount, count;
    RawGrid *rawgrid = GetRawGrid();

    int nTotalCell = GetNTotalCell();
    int nTotalNode = GetNTotalNode();
    int nBoundFace = GetNBoundFace();
    grid->SetNTotalCell(nTotalCell);
    grid->SetNTotalNode(nTotalNode);
    grid->SetNBoundFace(nBoundFace);

    RDouble *x = new RDouble[nTotalNode];
    RDouble *y = new RDouble[nTotalNode];
    RDouble *z = new RDouble[nTotalNode];
    grid->SetX(x);
    grid->SetY(y);
    grid->SetZ(z);

    RDouble *xx = rawgrid->GetX();
    RDouble *yy = rawgrid->GetY();
    RDouble *zz = rawgrid->GetZ();

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        x[iNode] = xx[iNode];
        y[iNode] = yy[iNode];
        z[iNode] = zz[iNode];
    }

    vector <string> bcNameList;
    vector <int> bcTypeList;
    set<string> bcNameSet;
    count = 0;

    UnstructBCSet **bcrs = new UnstructBCSet *[nBoundFace];
    vector< FluentFacePatch * > &face_patchs = GetFluentFacePatchS();
    for (std::size_t m = 0; m < face_patchs.size(); ++ m)
    {
        int bc_type = face_patchs[m]->GetBCType();
        string bcNam = face_patchs[m]->GetPatchName();
        if (bc_type == PHENGLEI::INTERIOR)
        {
            continue;
        }

        int nfsum = face_patchs[m]->GetNTotalFace();
        for (int i = 0; i < nfsum; ++ i)
        {
            bcTypeList.push_back(bc_type);
            bcNameList.push_back(bcNam);
            bcNameSet.insert(bcNameList[count]);

            bcrs[count] = new UnstructBCSet();
            bcrs[count]->SetKey(bc_type);
            bcrs[count]->SetBCName(bcNam);

            ++ count;
        }
    }

    grid->SetBCRecord(bcrs);

    int nBCRegionUnstruct = static_cast<int>(bcNameSet.size());
    grid->CreateUnstructBCSet(nBCRegionUnstruct);

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = new int[nBoundFace];
    set<string>::iterator iter;
    
    count = 0;
    for (iter = bcNameSet.begin(); iter != bcNameSet.end(); iter++)
    {
        UnstructBC *unstructBC = new UnstructBC(count);
        unstructBCSet->SetBCRegion(count, unstructBC);
        unstructBC->SetBCName(*iter);

        for (int iFace = 0; iFace < nBoundFace; iFace++)
        {
            if (bcNameList[iFace] == unstructBC->GetBCName())
            {
                unstructBC->SetBCType(bcTypeList[iFace]);
                bcRegionIDofBCFace[iFace] = count;
                vector<int> *faceIndex = unstructBC->GetFaceIndex();
                faceIndex->push_back(iFace);
            }

        }
        count++;
    }
    unstructBCSet->SetBCFaceInfo(bcRegionIDofBCFace);

    int nTotalFace = GetNTotalFace();
    grid->SetNTotalFace(nTotalFace);

    int *node_number_of_each_face = new int[nTotalFace];
    count = 0;
    for (std::size_t m = 0; m < face_patchs.size(); ++ m)
    {
        int nfsum = face_patchs[m]->GetNTotalFace();
        vector< vector<int> > &face2nodetmp = face_patchs[m]->GetFace2Node();

        for (int i = 0; i < nfsum; ++ i)
        {
            node_number_of_each_face[count++] = static_cast<int>(face2nodetmp[i].size());
        }
    }

    grid->SetNodeNumberOfEachFace(node_number_of_each_face);

    int face2node_size = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        face2node_size += node_number_of_each_face[iFace];
    }
    int *face2node = new int[face2node_size];
    grid->SetFace2Node(face2node);

    fcount = 0;
    count  = 0;
    for (std::size_t m = 0; m < face_patchs.size(); ++ m)
    {
        int nfsum = face_patchs[m]->GetNTotalFace();
        vector< vector<int> > &face2nodetmp = face_patchs[m]->GetFace2Node();

        for (int i = 0; i < nfsum; ++ i)
        {
            for (int j = 0; j < node_number_of_each_face[fcount]; ++ j)
            {
                face2node[count++] = face2nodetmp[i][j];
            }
            fcount ++;
        }
    }

    int *left_cell_of_face  = new int[nTotalFace];
    int *right_cell_of_face = new int[nTotalFace];
    grid->SetLeftCellOfFace(left_cell_of_face);
    grid->SetRightCellOfFace(right_cell_of_face);

    fcount = 0;
    for (std::size_t m = 0; m < face_patchs.size(); ++ m)
    {
        int nfsum = face_patchs[m]->GetNTotalFace();
        vector<int> &left_cell_of_face_tmp  = face_patchs[m]->GetLeftCellOfFace();
        vector<int> &right_cell_of_face_tmp = face_patchs[m]->GetRightCellOfFace();

        for (int i = 0; i < nfsum; ++ i)
        {

            if (left_cell_of_face_tmp[i] < 0)
            {
                left_cell_of_face[fcount] = right_cell_of_face_tmp[i];
                right_cell_of_face[fcount] = left_cell_of_face_tmp[i];
            }
            else
            {
                left_cell_of_face[fcount] = left_cell_of_face_tmp[i];
                right_cell_of_face[fcount] = right_cell_of_face_tmp[i];
            }
            
            fcount ++;
        }
    }
}

}