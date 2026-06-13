#include "Gmesh.h"
#include "TK_Parse.h"
#include "Math_BasisFunction.h"
#include "Pre_GridBase.h"
#include "Glb_Dimension.h"
#include "TK_Exit.h"
using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{
void ElemReorderGmsh2CGNS(int etype, int *node, int *node_cgns, int nvertex)
{
    for (int i = 0; i < nvertex; ++ i)
    {
        node_cgns[i] = node[i];
    }

    if (etype == PYRA_14)
    {
        //! Pyramid
        node_cgns[6 ] = node[8 ];
        node_cgns[7 ] = node[10]; 
        node_cgns[8 ] = node[6 ];
        node_cgns[9 ] = node[7 ];
        node_cgns[10] = node[9 ];
    }
    else if (etype == TETRA_10)
    {
        //! Tetrahedron
        node_cgns[8 ] = node[9 ];
        node_cgns[9 ] = node[8 ]; 
    }
    else if (etype == PENTA_15)
    {
        //! TriangularPrism
        node_cgns[7 ] = node[9 ];
        node_cgns[8 ] = node[7 ];
        node_cgns[9 ] = node[8 ];
        node_cgns[13] = node[14];
        node_cgns[14] = node[13];
    }
    else if (etype == PENTA_18)
    {
        //! TriangularPrism
        node_cgns[7 ] = node[9 ];
        node_cgns[8 ] = node[7 ];
        node_cgns[9 ] = node[8 ];
        node_cgns[13] = node[14];
        node_cgns[14] = node[13];
        node_cgns[16] = node[17];
        node_cgns[17] = node[16];
    }
    else if (etype == HEXA_20)
    {
        //! Hexahedron
        node_cgns[9 ] = node[11];
        node_cgns[10] = node[13];
        node_cgns[11] = node[9 ];
        node_cgns[12] = node[10];
        node_cgns[13] = node[12];
        node_cgns[17] = node[18];
        node_cgns[18] = node[19];
        node_cgns[19] = node[17];
    }
    else if (etype == HEXA_27)
    {
        //! Hexahedron
        node_cgns[9 ] = node[11];
        node_cgns[10] = node[13];
        node_cgns[11] = node[9 ];
        node_cgns[12] = node[10];
        node_cgns[13] = node[12];
        node_cgns[17] = node[18];
        node_cgns[18] = node[19];
        node_cgns[19] = node[17];

        node_cgns[22] = node[23];
        node_cgns[23] = node[24];
        node_cgns[24] = node[22];
    }
    else
    {
        ;
    }
}

GmeshElem::GmeshElem()
{
    ntype      = 0;
    conn_list  = 0;
    type       = 0;
    istart     = 0;
    iend       = 0;
    number     = 0;
    nTotalCell = 0;
}

GmeshElem::~GmeshElem()
{
    //conn_list = new int * [ntype];
    for (int i = 0; i < ntype; ++ i)
    {
        delete [] conn_list[i];
    }
    delete [] conn_list;
    delete [] type;
    delete [] istart;
    delete [] iend;
    delete [] number;
}

void GmeshElem::Create(int ntype)
{
    this->ntype = ntype;
    conn_list = new int * [ntype];
    type      = new int [ntype];
    istart    = new int [ntype];
    iend      = new int [ntype];
    number    = new int [ntype];
}

void CreateValueSet(vector<int> &value_container, set<int> &value_set)
{
    for (std::size_t i = 0; i < value_container.size(); ++ i)
    { 
        value_set.insert(value_container[i]);
    }
}

void GmeshElem::TreatElements(CGNSFactory *factory_cgns, vector<int> &type_container, vector<int> &number_container, vector<int> &index_container, vector<int> &bctype_container, int *e_nvertex)
{
    //找出所有不同的单元类型 
    set<int> elem_type_set;
    for (std::size_t i = 0; i < type_container.size(); ++ i)
    {
        int etype = type_container[i];
        elem_type_set.insert(etype);
    }

    int nSections = static_cast<int>(elem_type_set.size());

    this->Create(nSections);

    Set2Array(elem_type_set, type);

    for (int m = 0; m < nSections; ++ m)
    {
        int nelem_this_type = 0;
        for (std::size_t i = 0; i < type_container.size(); ++ i)
        {
            int etype     = type_container[i];
            int nelem_num = number_container[i];
            if (etype == type[m])
            {
                nelem_this_type += nelem_num;
            }
        }
        int nvertex  = e_nvertex[type[m]];
        conn_list[m] = new int[nelem_this_type * nvertex];
        number[m]    = nelem_this_type;
    }

    nTotalCell = 0;
    for (int m = 0; m < nSections; ++ m)
    {
        nTotalCell += number[m];
    }

    istart[0] = 1;
    for (int m = 0; m < nSections; ++ m)
    {
        iend[m] = istart[m] + number[m] - 1;
        if (m < nSections - 1)
        {
            istart[m+1] = iend[m] + 1;
        }
    }

    for (int m = 0; m < nSections; ++ m)
    {
        int icount = 0;
        int jcount = 0;
        for (std::size_t i = 0; i < type_container.size(); ++ i)
        {
            int etype     = type_container[i];
            int nelem_num = number_container[i];
            int nvertex   = e_nvertex[ etype ];
            int ptmax     = nelem_num * nvertex;

            if (etype == type[m])
            {
                for (int j = 0; j < ptmax; ++ j)
                {
                    conn_list[m][icount + j] = index_container[jcount + j];
                }
                icount += ptmax;
            }
            jcount += ptmax;
        }
    }

    //找出所有不同的边界类型 
    set<int> bc_type_set;
    CreateValueSet(bctype_container, bc_type_set);

    //删除内边界 
    bc_type_set.erase(0);

    //! Determine the number of boundary conditions for this zone.
    int  nBCRegions   = static_cast<int>(bc_type_set.size());
    int *bctype_array = new int[nBCRegions];

    Set2Array(bc_type_set, bctype_array);

    CGNSBase *base_cgns = factory_cgns->GetCGNSBase(0);
    base_cgns->CreateBC(nBCRegions);

    int  *base_cgns_nBCElem          = base_cgns->GetNumberOfBCElements();
    int  *base_cgns_bc_elem_set_type = base_cgns->GetBoundaryElementType();
    int  *base_cgns_bc_grid_location = base_cgns->GetBoundaryGridLocation();
    int **base_cgns_bc_conn_list     = (int **)base_cgns->GetBoundaryElementConnectionList();
    cout << "    Warning: the cgsize_t type of conn_list has been convert to int type, when deal with GMSH grid!!!" << endl;
    PointSetType_t ptsetType = PointList;

    //! Loop over the number of boundary conditions.
    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        GridLocation_t igr = FaceCenter;
        cout << " Grid Location Name = " << GridLocationName[igr] << "\n";

        if (igr == FaceCenter)
        {
            cout << "\nGridLocation = FaceCenter means BC data refers to elements, not nodes\n";
        }
        else if (igr == Vertex)
        {
            cout << "\nGridLocation = Vertex means BC data refers to nodes, not elements\n";
        }

        int bctype  = bctype_array[iBCRegion];
        int nBCElem = 0;
        for (std::size_t i = 0; i < bctype_container.size(); ++ i)
        {
            if (bctype == bctype_container[i])
            {
                ++ nBCElem;
            }
        }

        base_cgns_nBCElem[iBCRegion]          = nBCElem;
        base_cgns_bc_conn_list[iBCRegion]     = new int[nBCElem];
        base_cgns_bc_elem_set_type[iBCRegion] = ptsetType;
        base_cgns_bc_grid_location[iBCRegion] = igr;

        cout << "PointSetType_t = " << ptsetType << "\n";

        if (nBCElem == 2)
        {
            cout << base_cgns_bc_conn_list[iBCRegion][0] << " " << base_cgns_bc_conn_list[iBCRegion][1] << "\n";
        }
        else
        {
            for (int j = 0; j < MIN(20, nBCElem); ++ j)
            {
                cout << base_cgns_bc_conn_list[iBCRegion][j] << " ";
            }
            cout << "\n";
        }
    }

    delete [] bctype_array;
}

Gmesh::Gmesh()
{
    //vector<string> keyword;
    keyword.push_back("MeshFormat");
    keyword.push_back("Nodes");
    keyword.push_back("$Elements");
    keyword.push_back("PhysicalNames");
    keyword.push_back("NodeData");
    keyword.push_back("ElementData");
    keyword.push_back("ElementNodeData");

    elem_pt = new int [NofValidElementTypes];

    elem_pt[NODE   ] = 1;
    elem_pt[BAR_2  ] = 2;
    elem_pt[TRI_3  ] = 3;
    elem_pt[QUAD_4 ] = 4;
    elem_pt[TETRA_4] = 4;
    elem_pt[PYRA_5 ] = 5;
    elem_pt[PENTA_6] = 6;
    elem_pt[HEXA_8 ] = 8;

    NofValidElementTypesGmsh = MSH_NUM_TYPE;
    TypeGmsh2CGNS = new int [NofValidElementTypesGmsh+1];

    InitGmsh2CGNSMap();

    gmesh_nvertex = new int [NofValidElementTypesGmsh];

    InitGmshNVertex();

    cgns_nvertex = new int [NofValidElementTypesGmsh];

    InitCGNSNVertex();
}

Gmesh::~Gmesh()
{
    delete [] elem_pt;
    delete [] TypeGmsh2CGNS;
    delete [] gmesh_nvertex;
    delete [] cgns_nvertex;
}

void Gmesh::InitGmshNVertex()
{
    for (int i = 0; i < NofValidElementTypesGmsh; ++ i)
    {
        gmesh_nvertex[i] = 1;
    }

    gmesh_nvertex[ PHSPACE::PH_MSH_LIN_2   ] = 2;
    gmesh_nvertex[ PHSPACE::PH_MSH_TRI_3   ] = 3;
    gmesh_nvertex[ PHSPACE::PH_MSH_QUA_4   ] = 4;
    gmesh_nvertex[ PHSPACE::PH_MSH_TET_4   ] = 4;
    gmesh_nvertex[ PHSPACE::PH_MSH_HEX_8   ] = 8;
    gmesh_nvertex[ PHSPACE::PH_MSH_PRI_6   ] = 6;
    gmesh_nvertex[ PHSPACE::PH_MSH_PYR_5   ] = 5;
    gmesh_nvertex[ PHSPACE::PH_MSH_LIN_3   ] = 3;
    gmesh_nvertex[ PHSPACE::PH_MSH_TRI_6   ] = 6;
    gmesh_nvertex[ PHSPACE::PH_MSH_QUA_9   ] = 9;
    gmesh_nvertex[ PHSPACE::PH_MSH_TET_10  ] = 10;
    gmesh_nvertex[ PHSPACE::PH_MSH_HEX_27  ] = 27;
    gmesh_nvertex[ PHSPACE::PH_MSH_PRI_18  ] = 18;
    gmesh_nvertex[ PHSPACE::PH_MSH_PYR_14  ] = 14;
    gmesh_nvertex[ PHSPACE::PH_MSH_PNT     ] = 1;
    gmesh_nvertex[ PHSPACE::PH_MSH_QUA_8   ] = 8;
    gmesh_nvertex[ PHSPACE::PH_MSH_HEX_20  ] = 20;
    gmesh_nvertex[ PHSPACE::PH_MSH_PRI_15  ] = 15;
    gmesh_nvertex[ PHSPACE::PH_MSH_PYR_13  ] = 13;
}

void Gmesh::InitGmsh2CGNSMap()
{
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_2   ] = BAR_2;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_3   ] = TRI_3;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_4   ] = QUAD_4;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_4   ] = TETRA_4;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_8   ] = HEXA_8;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PRI_6   ] = PENTA_6;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PYR_5   ] = PYRA_5;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_3   ] = BAR_3;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_6   ] = TRI_6;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_9   ] = QUAD_9;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_10  ] = TETRA_10;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_27  ] = HEXA_27;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PRI_18  ] = PENTA_18;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PYR_14  ] = PYRA_14;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PNT     ] = NODE;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_8   ] = QUAD_8;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_20  ] = HEXA_20;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PRI_15  ] = PENTA_15;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PYR_13  ] = PYRA_13;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_9   ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_10  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_12  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_15  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_15I ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_21  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_4   ] = BAR_3;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_5   ] = BAR_3;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_6   ] = BAR_3;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_20  ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_35  ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_56  ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_34  ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_52  ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_POLYG_  ] = NGON_n;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_POLYH_  ] = NGON_n;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_16  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_25  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_36  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_12  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_16I ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_20  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_28  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_36  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_45  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_55  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_66  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_49  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_64  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_81  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_100 ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_121 ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_18  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_21I ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_24  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_27  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_30  ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_24  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_28  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_32  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_36I ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_40  ] = QUAD_9;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_7   ] = BAR_3;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_8   ] = BAR_3;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_9   ] = BAR_3;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_10  ] = BAR_3;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_11  ] = BAR_3;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_B   ] = BAR_3;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_B   ] = TRI_6;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_POLYG_B ] = NGON_n;      //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_C   ] = BAR_3;       //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_84  ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_120 ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_165 ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_220 ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_286 ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_64  ] = HEXA_20;     //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_125 ] = HEXA_20;     //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_196 ] = HEXA_20;     //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_74  ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_100 ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_130 ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_164 ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_202 ] = TETRA_10;    //CGNS没有对应类型 
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_NUM_TYPE] = NGON_n;      //CGNS没有对应类型 
}

void Gmesh::InitCGNSNVertex()
{
    for (int i = 0; i < NofValidElementTypes; ++ i)
    {
        cgns_nvertex[i] = 1;
    }

    cgns_nvertex[ ElementTypeNull         ] = 0;
    cgns_nvertex[ ElementTypeUserDefined  ] = 0;
    cgns_nvertex[ NODE                    ] = NPE_NODE;
    cgns_nvertex[ BAR_2                   ] = NPE_BAR_2;
    cgns_nvertex[ BAR_3                   ] = NPE_BAR_3;
    cgns_nvertex[ TRI_3                   ] = NPE_TRI_3;
    cgns_nvertex[ TRI_6                   ] = NPE_TRI_6;
    cgns_nvertex[ QUAD_4                  ] = NPE_QUAD_4;
    cgns_nvertex[ QUAD_8                  ] = NPE_QUAD_8;
    cgns_nvertex[ QUAD_9                  ] = NPE_QUAD_9;
    cgns_nvertex[ TETRA_4                 ] = NPE_TETRA_4;
    cgns_nvertex[ TETRA_10                ] = NPE_TETRA_10;
    cgns_nvertex[ PYRA_5                  ] = NPE_PYRA_5;
    cgns_nvertex[ PYRA_14                 ] = NPE_PYRA_13;
    cgns_nvertex[ PENTA_6                 ] = NPE_PYRA_14;
    cgns_nvertex[ PENTA_15                ] = NPE_PENTA_6;
    cgns_nvertex[ PENTA_18                ] = NPE_PENTA_15; 
    cgns_nvertex[ HEXA_8                  ] = NPE_PENTA_18;
    cgns_nvertex[ HEXA_20                 ] = NPE_HEXA_8;
    cgns_nvertex[ HEXA_27                 ] = NPE_HEXA_20;
    cgns_nvertex[ MIXED                   ] = NPE_HEXA_27;
    cgns_nvertex[ PYRA_13                 ] = NPE_MIXED;
    cgns_nvertex[ NGON_n                  ] = NPE_NGON_n;
    cgns_nvertex[ NFACE_n                 ] = NPE_NFACE_n;
}

void Gmesh::Read(fstream &file, CGNSFactory *factory_cgns)
{
    ReadMeshFormat(file);
    ReadNodes(file, factory_cgns);
    ReadElements(file, factory_cgns);
}

void Gmesh::ReadMeshFormat(fstream &file)
{
    string line, word;
    string separator = " =\t\r\n#$,;";
    do
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        transform(word.begin(), word.end(), word.begin(), ToLower());
    } while (word != "meshformat");

    getline(file, line);

    line = FindNextWord(line, word, separator);
    from_string< RDouble >(version_number, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string< int >(file_type, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string< int >(data_size, word, std::dec);

    cout << " version-number = " << version_number << "\n";
    cout << " file-type      = " << file_type      << "\n";
    cout << " data-size      = " << data_size      << "\n";
}

void Gmesh::ReadNodes(fstream &file, CGNSFactory *factory_cgns)
{
    string line, word;
    string separator = " =\t\r\n#$,;";
    do
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        transform(word.begin(), word.end(), word.begin(), ToLower());
    } while (word != "nodes");

    getline(file, line);

    line = FindNextWord(line, word, separator);
    from_string< int >(nTotalNode, word, std::dec);

    cout << "\nRead grid coordinates ...\n";

    int nBlocks = 1;
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        cout << " iZone = " << iZone << " nBlocks = " << nBlocks << endl;
        cout << " nTotalNode = " << nTotalNode << endl;

        RDouble *x = new RDouble[nTotalNode];
        RDouble *y = new RDouble[nTotalNode];
        RDouble *z = new RDouble[nTotalNode];

        RawGrid *rawgrid = factory_cgns->GetRawGrid(iZone);

        rawgrid->SetX(x);
        rawgrid->SetY(y);
        rawgrid->SetZ(z);

        RDouble xmin, ymin, zmin, xmax, ymax, zmax;
        xmin = ymin = zmin =   LARGE;
        xmax = ymax = zmax = - LARGE;

        int index;
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            getline(file,line);
            line = FindNextWord(line, word, separator);
            from_string< int >(index, word, std::dec);

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
    cout << "Read grid coordinates end \n";
}

void Gmesh::ReadElements(fstream &file, CGNSFactory *factory_cgns)
{
    string line, word;
    string separator = " =\t\r\n#$,;";
    do
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        transform(word.begin(), word.end(), word.begin(), ToLower());
    } while (word != "elements");

    getline(file, line);

    //所有单元，包括边界单元 
    int elements = 0;

    line = FindNextWord(line, word, separator);
    from_string< int >(elements, word, std::dec);

    cout << "\nRead element topology ...\n";

    int elem_index, gmsh_elem_type;
    int number_of_tags;
    int elem_tag[3];

    vector <int>           elem_type_container;
    vector <int>           elem_physId_container;
    vector <int>           elem_start_container;
    vector <int>           elem_end_container;
    vector < vector<int> > elem_nodes_containers;

    int *node_gmsh = new int[NofValidElementTypesGmsh]();
    int *node_cgns = new int[NofValidElementTypesGmsh]();

    int iElement      = 0;
    int ori_elem_type = -1;
    int ori_PhysId    = -1;

    getline(file, line);

    line = FindNextWord(line, word, separator);
    from_string< int >(elem_index, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string< int >(gmsh_elem_type, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string< int >(number_of_tags, word, std::dec);

    for (int j = 0; j < number_of_tags; ++ j)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(elem_tag[j], word, std::dec);
    }

    do
    {
        ori_PhysId    = elem_tag[0];
        ori_elem_type = gmsh_elem_type;

        int cgns_elem_type = TypeGmsh2CGNS[gmsh_elem_type];
        int cgns_vnum = cgns_nvertex[cgns_elem_type];
        int gesh_vnum = gmesh_nvertex[gmsh_elem_type];

        elem_type_container.push_back(cgns_elem_type);
        elem_physId_container.push_back(elem_tag[0]);
        elem_start_container.push_back(elem_index);

        vector <int> elem_nodes_container;

        do
        {
            int nodeindex;
            for (int j = 0; j < gesh_vnum; ++j)
            {
                line = FindNextWord(line, word, separator);
                from_string< int >(nodeindex, word, std::dec);
                node_gmsh[j] = nodeindex;
            }

            ElemReorderGmsh2CGNS(cgns_elem_type, node_gmsh, node_cgns, gesh_vnum);
            for (int j = 0; j < cgns_vnum; ++j)
            {
                elem_nodes_container.push_back(node_cgns[j]);
            }

            iElement ++;
            if (iElement == elements)
            {
                elem_index ++;
                break;
            }

            getline(file, line);

            line = FindNextWord(line, word, separator);
            from_string< int >(elem_index, word, std::dec);

            line = FindNextWord(line, word, separator);
            from_string< int >(gmsh_elem_type, word, std::dec);

            line = FindNextWord(line, word, separator);
            from_string< int >(number_of_tags, word, std::dec);

            for (int j = 0; j < number_of_tags; ++j)
            {
                line = FindNextWord(line, word, separator);
                from_string< int >(elem_tag[j], word, std::dec);
            }
        } while (gmsh_elem_type == ori_elem_type && elem_tag[0] == ori_PhysId);

        elem_end_container.push_back(elem_index - 1);
        elem_nodes_containers.push_back(elem_nodes_container);

    } while (iElement < elements);

    vector <int>    PhysicalDim;
    vector <int>    PhysicalId;
    vector <string> PhysicalName;
    ReadBoundary(file, PhysicalDim, PhysicalId, PhysicalName);

    ConvertCGNS(factory_cgns, elem_type_container, elem_nodes_containers, elem_start_container, 
                elem_end_container, elem_physId_container, cgns_nvertex, PhysicalDim, PhysicalId, PhysicalName);

    delete [] node_gmsh;    node_gmsh = nullptr;
    delete [] node_cgns;    node_cgns = nullptr;

    cout << "Read element topology end\n\n";
}

void Gmesh::ReadBoundary(fstream &file, vector <int> &PhysicalDim, vector <int> &PhysicalId, vector <string> &PhysicalName)
{
    string line, word;
    string separator = " =\t\r\n#$,;\"";
    do
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        transform(word.begin(), word.end(), word.begin(), ToLower());
    } while (word != "physicalnames");

    int nNames, iDim, elem_tag;
    string physicalNames;

    getline(file, line);
    line = FindNextWord(line, word, separator);
    from_string< int >(nNames, word, std::dec);

    for (int iName = 0; iName < nNames; ++ iName)
    {
        getline(file, line);

        line = FindNextWord(line, word, separator);
        from_string< int >(iDim, word, std::dec);
        PhysicalDim.push_back(iDim);

        line = FindNextWord(line, word, separator);
        from_string< int >(elem_tag, word, std::dec);
        PhysicalId.push_back(elem_tag);

        line = FindNextWord(line, word, separator);
        from_string< string >(physicalNames, word, std::dec);
        PhysicalName.push_back(physicalNames);
    }
}

void Gmesh::ConvertCGNS(CGNSFactory *factory_cgns, vector<int> &type_container, vector < vector<int> > &elem_nodes_container, 
                        vector <int> &start_container, vector <int> &end_container, vector<int> &physId_container, 
                        int *e_nvertex, vector <int> &PhysicalDim, vector <int> &PhysicalId, vector <string> &PhysicalName)
{
    CGNSBase *base_cgns = factory_cgns->GetCGNSBase(0);

    int nSections = type_container.size();
    base_cgns->CreateElement(nSections);

    int       *base_cgns_elem_type       = base_cgns->GetElementType();
    cgsize_t  *base_cgns_istart          = base_cgns->GetIndexOfStart();
    cgsize_t  *base_cgns_iend            = base_cgns->GetIndexOfEnd();
    cgsize_t **base_cgns_conn_list       = base_cgns->GetElementConnectionList();
    cgsize_t  *base_cgns_elementdatasize = base_cgns->GetElementDataSize();

    for (int iSection = 0; iSection < nSections; ++ iSection)
    {
        cout << " Reading section "  << iSection + 1 << " ...\n";
        cout << "   section type = " << ElementTypeName[type_container[iSection]] << "\n";
        cout << "   istart,iend =  " << start_container[iSection] << " " << end_container[iSection] << "\n";

        base_cgns_istart[iSection]    = start_container[iSection];
        base_cgns_iend[iSection]      = end_container[iSection];
        base_cgns_elem_type[iSection] = type_container[iSection];

        int elementdatasize = elem_nodes_container[iSection].size();
        base_cgns_elementdatasize[iSection] = elementdatasize;

        base_cgns_conn_list[iSection] = new cgsize_t[elementdatasize];
        for (int iData = 0; iData < elementdatasize; ++ iData)
        {
            base_cgns_conn_list[iSection][iData] = elem_nodes_container[iSection][iData];
        }
    }

    int nBCRegions = 0;
    int globalDim  = GetDim();
    vector <int> BCRegionsIndex(nSections, -1);

    for (int jSection = 0; jSection < PhysicalDim.size(); ++ jSection)
    {
        if (PhysicalDim[jSection] == globalDim - 2 || PhysicalDim[jSection] > globalDim)
        {
            TK_Exit::ExceptionExit("Error: Dimension in Gmsh file is not according with that in parameter file!");
        }

        if (PhysicalDim[jSection] == globalDim)
        {
            continue;
        }

        for (int iSection = 0; iSection < nSections; ++ iSection)
        {
            if (physId_container[iSection] == PhysicalId[jSection])
            {
                BCRegionsIndex[iSection] = jSection;
                nBCRegions ++;
            }
        }
    }

    cout << "\n Number of BC Regions: " << nBCRegions << endl;
    base_cgns->CreateBC(nBCRegions);

    int       *base_cgns_nBCElem          = base_cgns->GetNumberOfBCElements();
    int       *base_cgns_bc_elem_set_type = base_cgns->GetBoundaryElementType();
    int       *base_cgns_bc_grid_location = base_cgns->GetBoundaryGridLocation();
    int       *base_cgns_bc_type          = base_cgns->GetBCType();
    cgsize_t **base_cgns_bc_conn_list     = base_cgns->GetBoundaryElementConnectionList();

    string *boundaryName = new string[nBCRegions]();
    base_cgns->SetBCName(boundaryName);

    PointSetType_t ptsetType = ElementRange;
    GridLocation_t igr       = FaceCenter;

    int iBCRegion = 0;
    for (int iSection = 0; iSection < nSections; ++ iSection)
    {
        if (BCRegionsIndex[iSection] < 0)
        {
            continue;
        }

        cout << " Reading bcRegion " << iBCRegion + 1 << " ...\n";
        cout << "   Boundary Name = " << PhysicalName[BCRegionsIndex[iSection]] << "\n";
        cout << "   BCtype  = " << BCTypeName[1] << "\n";

        boundaryName[iBCRegion]               = PhysicalName[BCRegionsIndex[iSection]];
        base_cgns_nBCElem[iBCRegion]          = 2;
        base_cgns_bc_conn_list[iBCRegion]     = new cgsize_t[2];
        base_cgns_bc_type[iBCRegion]          = 1;
        base_cgns_bc_elem_set_type[iBCRegion] = ptsetType;
        base_cgns_bc_grid_location[iBCRegion] = igr;

        base_cgns_bc_conn_list[iBCRegion][0] = start_container[iSection];
        base_cgns_bc_conn_list[iBCRegion][1] = end_container[iSection];

        cout << "  Start & End point of Region " << ": "
            << base_cgns_bc_conn_list[iBCRegion][0] << " " << base_cgns_bc_conn_list[iBCRegion][1] << "\n";

        iBCRegion ++;
    }

    nTotalCell = 0;
    for (int iSection = 0; iSection < nSections; ++ iSection)
    {
        if (BCRegionsIndex[iSection] < 0)
        {
            nTotalCell += end_container[iSection] - start_container[iSection] + 1;
        }
    }
    base_cgns->SetNTotalNode(nTotalNode);
    base_cgns->SetNTotalCell(nTotalCell);
}

}
