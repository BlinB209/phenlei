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
//! @file      Gmesh.h
//! @brief     Convert gmsh grid to fts grid.
//! @author    He Xin.

#pragma once
#include "Precision.h"
#include <vector>
#include <set>
#include <string>
using namespace std;

namespace PHSPACE
{
#include "GmshDefines.h"

typedef enum
{
    PH_MSH_LIN_2 = MSH_LIN_2,
    PH_MSH_TRI_3   ,
    PH_MSH_QUA_4   ,
    PH_MSH_TET_4   ,
    PH_MSH_HEX_8   ,
    PH_MSH_PRI_6   ,
    PH_MSH_PYR_5   ,
    PH_MSH_LIN_3   ,
    PH_MSH_TRI_6   ,
    PH_MSH_QUA_9   ,
    PH_MSH_TET_10  ,
    PH_MSH_HEX_27  ,
    PH_MSH_PRI_18  ,
    PH_MSH_PYR_14  ,
    PH_MSH_PNT     ,
    PH_MSH_QUA_8   ,
    PH_MSH_HEX_20  ,
    PH_MSH_PRI_15  ,
    PH_MSH_PYR_13  ,
    PH_MSH_TRI_9   ,
    PH_MSH_TRI_10  ,
    PH_MSH_TRI_12  ,
    PH_MSH_TRI_15  ,
    PH_MSH_TRI_15I ,
    PH_MSH_TRI_21  ,
    PH_MSH_LIN_4   ,
    PH_MSH_LIN_5   ,
    PH_MSH_LIN_6   ,
    PH_MSH_TET_20  ,
    PH_MSH_TET_35  ,
    PH_MSH_TET_56  ,
    PH_MSH_TET_34  ,
    PH_MSH_TET_52  ,
    PH_MSH_POLYG_  ,
    PH_MSH_POLYH_  ,
    PH_MSH_QUA_16  ,
    PH_MSH_QUA_25  ,
    PH_MSH_QUA_36  ,
    PH_MSH_QUA_12  ,
    PH_MSH_QUA_16I ,
    PH_MSH_QUA_20  ,
    PH_MSH_TRI_28  ,
    PH_MSH_TRI_36  ,
    PH_MSH_TRI_45  ,
    PH_MSH_TRI_55  ,
    PH_MSH_TRI_66  ,
    PH_MSH_QUA_49  ,
    PH_MSH_QUA_64  ,
    PH_MSH_QUA_81  ,
    PH_MSH_QUA_100 ,
    PH_MSH_QUA_121 ,
    PH_MSH_TRI_18  ,
    PH_MSH_TRI_21I ,
    PH_MSH_TRI_24  ,
    PH_MSH_TRI_27  ,
    PH_MSH_TRI_30  ,
    PH_MSH_QUA_24  ,
    PH_MSH_QUA_28  ,
    PH_MSH_QUA_32  ,
    PH_MSH_QUA_36I ,
    PH_MSH_QUA_40  ,
    PH_MSH_LIN_7   ,
    PH_MSH_LIN_8   ,
    PH_MSH_LIN_9   ,
    PH_MSH_LIN_10  ,
    PH_MSH_LIN_11  ,
    PH_MSH_LIN_B   ,
    PH_MSH_TRI_B   ,
    PH_MSH_POLYG_B ,
    PH_MSH_LIN_C   ,
    PH_MSH_TET_84  ,
    PH_MSH_TET_120 ,
    PH_MSH_TET_165 ,
    PH_MSH_TET_220 ,
    PH_MSH_TET_286 ,
    PH_MSH_HEX_64  ,
    PH_MSH_HEX_125 ,
    PH_MSH_HEX_196 ,
    PH_MSH_TET_74  ,
    PH_MSH_TET_100 ,
    PH_MSH_TET_130 ,
    PH_MSH_TET_164 ,
    PH_MSH_TET_202 ,
    PH_MSH_NUM_TYPE
} GmshElemEnum;

class CGNSFactory;

class GmeshElem
{
public:
    int ntype;       //有多少种不同的单元。
    int **conn_list; //连接表
    int *type;       //每种单元的类型
    int *istart;     //每种单元的起始标记
    int *iend;       //每种单元的终止标记
    int *number;     //每种单元的数量
    int nTotalCell;
public:
    GmeshElem();
    ~GmeshElem();
public:
    void Create(int ntype);
    int  GetNTotalCell() const { return nTotalCell; };
    void SetNTotalCell(int nTotalCell) { this->nTotalCell = nTotalCell; };
    void TreatElements(CGNSFactory *factory_cgns, vector<int> &type_container, vector<int> &number_container, vector<int> &index_container, vector<int> &bctype_container, int *e_nvertex);
};

class Gmesh
{
public:
    Gmesh();
    ~Gmesh();
public:
    void Read(fstream &file, CGNSFactory *factory_cgns);
private:
    void ReadMeshFormat(fstream &file);
    void ReadNodes(fstream &file, CGNSFactory *factory_cgns);
    void ReadElements(fstream &file, CGNSFactory *factory_cgns);
    void ReadBoundary(fstream &file, vector <int> &PhysicalDim, vector <int> &PhysicalId, vector <string> &PhysicalName);
    void ConvertCGNS(CGNSFactory *factory_cgns, vector<int> &type_container, vector < vector<int> > &elem_nodes_container, 
                     vector <int> &start_container, vector <int> &end_container, vector<int> &physId_container, 
                     int *e_nvertex, vector <int> &PhysicalDim, vector <int> &PhysicalId, vector <string> &PhysicalName);
public:
    int NofValidElementTypesGmsh;
    int *TypeGmsh2CGNS;
    int *gmesh_nvertex;
    int *cgns_nvertex;
    void InitGmshNVertex();
    void InitCGNSNVertex();
    void InitGmsh2CGNSMap();

private:
    int *elem_pt;
    //version-number is a real number equal to 2.2
    RDouble version_number;
    //
    //file-type is an integer equal to 0 in the ASCII file format.
    int file_type;

    //data-size is an integer equal to the size of the floating point numbers used in the file
    //(currently only data-size = sizeof(RDouble) is supported).
    int data_size;

    //number-of-nodes
    //is the number of nodes in the mesh.
    int nTotalNode, nTotalCell;

    vector<string> keyword;
};

void CreateValueSet(vector<int> &value_container, set<int> &value_set);

}

