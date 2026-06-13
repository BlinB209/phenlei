#include "Pre_GridBase.h"
#include "Math_BasisFunction.h"
#include "IO_FileName.h"
#include "HyList.h"
#include "TK_Exit.h"
#include "Geo_UnstructBC.h"
#include "Constants.h"
#include "TK_Log.h"
#include "Pre_Reorder.h"
#include <cmath>
#include <algorithm>


using namespace std;

namespace PHSPACE
{

CGNSRawCoor::CGNSRawCoor()
{
    type = new DataType_t[3];
    coor[0] = 0;
    coor[1] = 0;
    coor[2] = 0;
    nTotalNode = 0;
}

CGNSRawCoor::~CGNSRawCoor()
{
    DeAllocateData();
    delete [] type;
}

void CGNSRawCoor::AllocateData(int icoor, cgsize_t nTotalNode, DataType_t type)
{
    this->type[icoor] = type;
    this->nTotalNode = nTotalNode;
    if (type == RealSingle)
    {
        coor[icoor] = new RDouble[nTotalNode];
    }
    else
    {
        coor[icoor] = new RDouble[nTotalNode];
    }
}

void CGNSRawCoor::SetAllData(RDouble *x, RDouble *y, RDouble *z)
{
    RDouble *xyz[3];

    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;

    for (int icoor = 0; icoor < 3; ++ icoor)
    {
        DataType_t type = this->type[icoor];
        SetData(icoor, type, xyz[icoor]);
    }

    //! Compute the max && min box.
    RDouble xMin, yMin, zMin, xMax, yMax, zMax;
    xMin = yMin = zMin =   LARGE;
    xMax = yMax = zMax = - LARGE;

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        xMin = MIN(xMin, x[iNode]);
        yMin = MIN(yMin, y[iNode]);
        zMin = MIN(zMin, z[iNode]);

        xMax = MAX(xMax, x[iNode]);
        yMax = MAX(yMax, y[iNode]);
        zMax = MAX(zMax, z[iNode]);
    }

    ostringstream oss;
    oss << "Min && Max box of the CGNS grid: \n";
    oss << setiosflags(ios::right);
    oss << setprecision(8);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);
    int wordwidth = 16;
    oss << "  xMin = " << setw(wordwidth) << xMin << ", xMax = " << setw(wordwidth) << xMax << "\n";
    oss << "  yMin = " << setw(wordwidth) << yMin << ", yMax = " << setw(wordwidth) << yMax << "\n";
    oss << "  zMin = " << setw(wordwidth) << zMin << ", zMax = " << setw(wordwidth) << zMax << "\n";
    oss << "\n";
    PrintToWindow(oss);
}

void CGNSRawCoor::SetData(int icoor, DataType_t type, RDouble *x)
{
    if (type == RealSingle)
    {
        RDouble *data = static_cast<RDouble *>(coor[icoor]);
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            x[iNode] = data[iNode];
        }
    }
    else
    {
        RDouble *data = static_cast<RDouble *>(coor[icoor]);
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            x[iNode] = data[iNode];
        }
    }
}

void CGNSRawCoor::DeAllocateData()
{
    for (int icoor = 0; icoor < 3; ++ icoor)
    {
        int type = this->type[icoor];
        if (type == RealSingle)
        {
            RDouble *data = static_cast<RDouble *>(coor[icoor]);
            delete [] data;
        }
        else
        {
            RDouble *data = static_cast<RDouble *>(coor[icoor]);
            delete [] data;
        }
    }
}

MultiIndex::MultiIndex()
{
    this->size = 0;
    this->id   = -1;
    data       = 0;
}

MultiIndex::MultiIndex(int size, int id)
{
    this->size = size;
    this->id   = id;
    data       = new int [size];
}

MultiIndex::MultiIndex(const MultiIndex &rhs)
{
    this->size = rhs.size;
    this->id   = rhs.id;
    this->data = new int [size];

    for (int i = 0; i < size; ++ i)
    {
        this->data[i] = rhs.data[i];
    }
}

MultiIndex & MultiIndex::operator = (const MultiIndex &rhs)
{
    if (this == &rhs) return *this;

    if (this->size != rhs.size)
    {
        delete [] this->data;
        this->size = rhs.size;
        this->data = new int [this->size];
    }

    this->id = rhs.id;
    for (int i = 0; i < size; ++ i)
    {
        this->data[i] = rhs.data[i];
    }

    return *this;
}

MultiIndex::~MultiIndex()
{
    delete [] data;
}

bool MultiIndex::operator < (const MultiIndex &rhs) const
{
    if (this->size != rhs.size) return this->size < rhs.size;
    for (int i = 0; i < this->size; ++ i)
    {
        if (this->data[i] != rhs.data[i])
        {
            return this->data[i] < rhs.data[i];
        }
    }
    return false;
}

void MultiIndex::SetData(int *data)
{
    for (int i = 0; i < size; ++ i)
    {
        this->data[i] = data[i];
    }
}

BaseElement::BaseElement()
{
    size         = NofValidElementTypes;
    elementface  = new CElemFace[size];
    elementpoint = new cgsize_t[size];
    InitElementPoint();
    InitElementFace();
};

BaseElement::~BaseElement()
{
    delete [] elementface;
    delete [] elementpoint;
};

void BaseElement::InitElementPoint()
{
    elementpoint[NODE ]    = 1;
    elementpoint[BAR_2]    = 2;
    elementpoint[TRI_3]    = 3;
    elementpoint[TRI_6]    = 6;
    elementpoint[QUAD_4]   = 4;
    elementpoint[QUAD_9]   = 9;
    elementpoint[TETRA_4]  = 4;
    elementpoint[TETRA_10] = 10;
    elementpoint[PYRA_5 ]  = 5;
    elementpoint[PYRA_14]  = 14;
    elementpoint[PENTA_6]  = 6;
    elementpoint[PENTA_18] = 18;
    elementpoint[HEXA_8 ]  = 8;
    elementpoint[HEXA_27]  = 27;
};

void BaseElement::InitElementFace()
{
    elementface[NODE].Init(NODE);
    elementface[BAR_2].Init(BAR_2);
    elementface[TRI_3].Init(TRI_3);
    elementface[QUAD_4].Init(QUAD_4);
    elementface[TETRA_4].Init(TETRA_4);
    elementface[TETRA_10].Init(TETRA_10);
    elementface[PYRA_5].Init(PYRA_5);
    elementface[PYRA_14].Init(PYRA_14);
    elementface[PENTA_6].Init(PENTA_6);
    elementface[PENTA_18].Init(PENTA_18);
    elementface[HEXA_8].Init(HEXA_8);
    elementface[HEXA_27].Init(HEXA_27);
};

CElemFace::CElemFace()
{
    face_num     = 0;
    face_clkwise = 0;
    facept_num   = 0;
    face_id      = 0;
}

CElemFace::~CElemFace()
{
    if (face_id)
    {
        for (int i = 0; i < face_num; ++ i)
        {
            delete [] face_id[i];
        }
    }
    delete [] face_id;
    delete [] face_clkwise;
    delete [] facept_num;
}

void CElemFace::Init(int elem_type)
{
    if (elem_type == NODE)
    {
        face_num = 1;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 1;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
    }
    else if (elem_type == BAR_2)
    {
        face_num = 1;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 2;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 1;
    }
    else if (elem_type == TRI_3)
    {
        face_num = 3;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 2;
        facept_num[1] = 2;
        facept_num[2] = 2;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 1;

        face_id[1][0] = 1;
        face_id[1][1] = 2;

        face_id[2][0] = 2;
        face_id[2][1] = 0;
    }
    else if (elem_type == QUAD_4)
    {
        face_num = 4;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 2;
        facept_num[1] = 2;
        facept_num[2] = 2;
        facept_num[3] = 2;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 1;

        face_id[1][0] = 1;
        face_id[1][1] = 2;

        face_id[2][0] = 2;
        face_id[2][1] = 3;

        face_id[3][0] = 3;
        face_id[3][1] = 0;
    }
    else if (elem_type == TETRA_4)
    {
        face_num = 4;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 3;
        facept_num[1] = 3;
        facept_num[2] = 3;
        facept_num[3] = 3;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 2;
        face_id[0][2] = 1;

        face_id[1][0] = 0;
        face_id[1][1] = 1;
        face_id[1][2] = 3;

        face_id[2][0] = 1;
        face_id[2][1] = 2;
        face_id[2][2] = 3;

        face_id[3][0] = 2;
        face_id[3][1] = 0;
        face_id[3][2] = 3;
    }
    else if (elem_type == TETRA_10)
    {
        face_num = 4;
        face_id = new int * [face_num];
        facept_num = new int[face_num];

        facept_num[0] = 6;
        facept_num[1] = 6;
        facept_num[2] = 6;
        facept_num[3] = 6;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int[facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 1;
        face_id[0][2] = 2;
        face_id[0][3] = 4;
        face_id[0][4] = 5;
        face_id[0][5] = 6;

        face_id[1][0] = 0;
        face_id[1][1] = 1;
        face_id[1][2] = 3;
        face_id[1][3] = 4;
        face_id[1][4] = 8;
        face_id[1][5] = 7;

        face_id[2][0] = 1;
        face_id[2][1] = 2;
        face_id[2][2] = 3;
        face_id[2][3] = 5;
        face_id[2][4] = 9;
        face_id[2][5] = 8;

        face_id[3][0] = 2;
        face_id[3][1] = 3;
        face_id[3][2] = 0;
        face_id[3][3] = 9;
        face_id[3][4] = 7;
        face_id[3][5] = 6;
    }
    else if (elem_type == PYRA_5)
    {
        face_num = 5;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 4;
        facept_num[1] = 3;
        facept_num[2] = 3;
        facept_num[3] = 3;
        facept_num[4] = 3;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 3;
        face_id[0][2] = 2;
        face_id[0][3] = 1;

        face_id[1][0] = 0;
        face_id[1][1] = 1;
        face_id[1][2] = 4;

        face_id[2][0] = 1;
        face_id[2][1] = 2;
        face_id[2][2] = 4;

        face_id[3][0] = 2;
        face_id[3][1] = 3;
        face_id[3][2] = 4;

        face_id[4][0] = 3;
        face_id[4][1] = 0;
        face_id[4][2] = 4;
    }
    else if (elem_type == PYRA_14)
    {
        face_num = 5;
        face_id = new int * [face_num];
        facept_num = new int[face_num];

        facept_num[0] = 9;
        facept_num[1] = 6;
        facept_num[2] = 6;
        facept_num[3] = 6;
        facept_num[4] = 6;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int[facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 3;
        face_id[0][2] = 2;
        face_id[0][3] = 1;
        face_id[0][4] = 8;
        face_id[0][5] = 7;
        face_id[0][6] = 6;
        face_id[0][7] = 5;
        face_id[0][8] = 13;

        face_id[1][0] = 0;
        face_id[1][1] = 1;
        face_id[1][2] = 4;
        face_id[1][3] = 5;
        face_id[1][4] = 10;
        face_id[1][5] = 9;

        face_id[2][0] = 1;
        face_id[2][1] = 2;
        face_id[2][2] = 4;
        face_id[2][3] = 6;
        face_id[2][4] = 11;
        face_id[2][5] = 10;

        face_id[3][0] = 2;
        face_id[3][1] = 3;
        face_id[3][2] = 4;
        face_id[3][3] = 7;
        face_id[3][4] = 12;
        face_id[3][5] = 11;

        face_id[4][0] = 3;
        face_id[4][1] = 0;
        face_id[4][2] = 4;
        face_id[4][3] = 8;
        face_id[4][4] = 9;
        face_id[4][5] = 12;
    }
    else if (elem_type == PENTA_6)
    {
        face_num = 5;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 3;
        facept_num[1] = 3;
        facept_num[2] = 4;
        facept_num[3] = 4;
        facept_num[4] = 4;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 2;
        face_id[0][2] = 1;

        face_id[1][0] = 3;
        face_id[1][1] = 4;
        face_id[1][2] = 5;

        face_id[2][0] = 0;
        face_id[2][1] = 1;
        face_id[2][2] = 4;
        face_id[2][3] = 3;

        face_id[3][0] = 1;
        face_id[3][1] = 2;
        face_id[3][2] = 5;
        face_id[3][3] = 4;

        face_id[4][0] = 2;
        face_id[4][1] = 0;
        face_id[4][2] = 3;
        face_id[4][3] = 5;
    }
    else if (elem_type == PENTA_18)
    {
        face_num = 5;
        face_id = new int * [face_num];
        facept_num = new int[face_num];

        facept_num[0] = 9;
        facept_num[1] = 9;
        facept_num[2] = 9;
        facept_num[3] = 6;
        facept_num[4] = 6;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int[facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 1;
        face_id[0][2] = 4;
        face_id[0][3] = 3;
        face_id[0][4] = 6;
        face_id[0][5] = 10;
        face_id[0][6] = 12;
        face_id[0][7] = 9;
        face_id[0][8] = 15;

        face_id[1][0] = 1;
        face_id[1][1] = 2;
        face_id[1][2] = 5;
        face_id[1][3] = 4;
        face_id[1][4] = 7;
        face_id[1][5] = 11;
        face_id[1][6] = 13;
        face_id[1][7] = 10;
        face_id[1][8] = 16;

        face_id[2][0] = 2;
        face_id[2][1] = 0;
        face_id[2][2] = 3;
        face_id[2][3] = 5;
        face_id[2][4] = 8;
        face_id[2][5] = 9;
        face_id[2][6] = 14;
        face_id[2][7] = 11;
        face_id[2][8] = 17;

        face_id[3][0] = 0;
        face_id[3][1] = 2;
        face_id[3][2] = 1;
        face_id[3][3] = 8;
        face_id[3][4] = 7;
        face_id[3][5] = 6;

        face_id[4][0] = 3;
        face_id[4][1] = 4;
        face_id[4][2] = 5;
        face_id[4][3] = 12;
        face_id[4][4] = 13;
        face_id[4][5] = 14;
    }
    else if (elem_type == HEXA_8)
    {
        face_num = 6;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 4;
        facept_num[1] = 4;
        facept_num[2] = 4;
        facept_num[3] = 4;
        facept_num[4] = 4;
        facept_num[5] = 4;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 3;
        face_id[0][2] = 2;
        face_id[0][3] = 1;

        face_id[1][0] = 4;
        face_id[1][1] = 5;
        face_id[1][2] = 6;
        face_id[1][3] = 7;

        face_id[2][0] = 3;
        face_id[2][1] = 0;
        face_id[2][2] = 4;
        face_id[2][3] = 7;

        face_id[3][0] = 0;
        face_id[3][1] = 1;
        face_id[3][2] = 5;
        face_id[3][3] = 4;

        face_id[4][0] = 1;
        face_id[4][1] = 2;
        face_id[4][2] = 6;
        face_id[4][3] = 5;

        face_id[5][0] = 2;
        face_id[5][1] = 3;
        face_id[5][2] = 7;
        face_id[5][3] = 6;
    }
    else if (elem_type == HEXA_27)
    {
        face_num = 6;
        face_id = new int * [face_num];
        facept_num = new int[face_num];

        facept_num[0] = 9;
        facept_num[1] = 9;
        facept_num[2] = 9;
        facept_num[3] = 9;
        facept_num[4] = 9;
        facept_num[5] = 9;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int[facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 3;
        face_id[0][2] = 2;
        face_id[0][3] = 1;
        face_id[0][4] = 11;
        face_id[0][5] = 10;
        face_id[0][6] = 9;
        face_id[0][7] = 8;
        face_id[0][8] = 20;

        face_id[1][0] = 0;
        face_id[1][1] = 1;
        face_id[1][2] = 5;
        face_id[1][3] = 4;
        face_id[1][4] = 8;
        face_id[1][5] = 13;
        face_id[1][6] = 16;
        face_id[1][7] = 12;
        face_id[1][8] = 21;

        face_id[2][0] = 1;
        face_id[2][1] = 2;
        face_id[2][2] = 6;
        face_id[2][3] = 5;
        face_id[2][4] = 9;
        face_id[2][5] = 14;
        face_id[2][6] = 17;
        face_id[2][7] = 13;
        face_id[2][8] = 22;

        face_id[3][0] = 2;
        face_id[3][1] = 3;
        face_id[3][2] = 7;
        face_id[3][3] = 6;
        face_id[3][4] = 10;
        face_id[3][5] = 15;
        face_id[3][6] = 18;
        face_id[3][7] = 14;
        face_id[3][8] = 23;

        face_id[4][0] = 0;
        face_id[4][1] = 4;
        face_id[4][2] = 7;
        face_id[4][3] = 3;
        face_id[4][4] = 12;
        face_id[4][5] = 19;
        face_id[4][6] = 15;
        face_id[4][7] = 11;
        face_id[4][8] = 24;

        face_id[5][0] = 4;
        face_id[5][1] = 5;
        face_id[5][2] = 6;
        face_id[5][3] = 7;
        face_id[5][4] = 16;
        face_id[5][5] = 17;
        face_id[5][6] = 18;
        face_id[5][7] = 19;
        face_id[5][8] = 25;
    }
    else
    {
        TK_Exit::ExceptionExit(" CElemFace::Init, Unknown type\n");
    }
}

BaseBCType::BaseBCType()
{
    ibcmax = NofValidBCTypes;
    Init();
};

void BaseBCType::Init()
{
    using namespace PHSPACE;

    bctypemap.insert(pair<int, int>(BCTypeUserDefined  , PHENGLEI::USER_DEFINED         ));
    bctypemap.insert(pair<int, int>(BCSymmetryPlane    , PHENGLEI::SYMMETRY             ));
    bctypemap.insert(pair<int, int>(BCInflow           , PHENGLEI::INFLOW               ));
    bctypemap.insert(pair<int, int>(BCOutflow          , PHENGLEI::OUTFLOW              ));
    bctypemap.insert(pair<int, int>(BCTunnelOutflow    , PHENGLEI::OUTFLOW              ));
    bctypemap.insert(pair<int, int>(BCOutflowSupersonic, PHENGLEI::OUTFLOW              ));
    bctypemap.insert(pair<int, int>(BCWall             , PHENGLEI::SOLID_SURFACE        ));
    bctypemap.insert(pair<int, int>(BCWallInviscid     , PHENGLEI::SOLID_SURFACE        ));
    bctypemap.insert(pair<int, int>(BCWallViscous      , PHENGLEI::SOLID_SURFACE        ));
    bctypemap.insert(pair<int, int>(BCDegenerateLine   , PHENGLEI::POLE                 ));
    bctypemap.insert(pair<int, int>(BCTypeNull         , PHENGLEI::NO_BOUNDARY_CONDITION));
    bctypemap.insert(pair<int, int>(BCFarfield         , PHENGLEI::FARFIELD             ));
    bctypemap.insert(pair<int, int>(BCGeneral          , PHENGLEI::GENERIC_1            ));

    fts2cgns.insert(pair<int, int>(PHENGLEI::USER_DEFINED         , BCTypeUserDefined));
    fts2cgns.insert(pair<int, int>(PHENGLEI::SYMMETRY             , BCSymmetryPlane  ));
    fts2cgns.insert(pair<int, int>(PHENGLEI::INFLOW               , BCInflow         ));
    fts2cgns.insert(pair<int, int>(PHENGLEI::OUTFLOW              , BCOutflow        ));
    fts2cgns.insert(pair<int, int>(PHENGLEI::SOLID_SURFACE        , BCWall           ));
    fts2cgns.insert(pair<int, int>(PHENGLEI::POLE                 , BCDegenerateLine ));
    fts2cgns.insert(pair<int, int>(PHENGLEI::NO_BOUNDARY_CONDITION, BCTypeNull       ));
    fts2cgns.insert(pair<int, int>(PHENGLEI::FARFIELD             , BCFarfield       ));

    int maxNumberOfBody = 20;
    for (int iBody = 0; iBody < maxNumberOfBody; iBody ++)
    {
        bctypemap.insert(pair<int, int>(NofValidBCTypes + iBody, PHENGLEI::StartNumberOfBodyOfHyperFLOW + iBody));
        fts2cgns. insert(pair<int, int>(PHENGLEI::StartNumberOfBodyOfHyperFLOW + iBody, NofValidBCTypes + iBody));
    }
}

BaseBCType::~BaseBCType()
{
}

RawGrid::RawGrid()
{
    x          = 0;
    y          = 0;
    z          = 0;
    nTotalNode = 0;
}

RawGrid::~RawGrid()
{
    //delete [] x;
    //delete [] y;
    //delete [] z;
}

Base_Grid_Conn::Base_Grid_Conn()
{
    left_cell_of_face        = new vector<cgsize_t>;
    right_cell_of_face       = new vector<cgsize_t>;
    face_location            = new vector<int>;
    bcFace                   = new vector<CGNS_BCFace *>;
    node_number_of_each_face = new vector<int>;
    face2node                = new vector< vector<cgsize_t> >;
    faceset                  = new set<MultiIndex>;
    node_number_of_each_cell = new vector<int>;
    cell2node                = new vector< vector<cgsize_t> >;
    nTotalNode               = 0;
    nTotalCell               = 0;
    nBoundFace               = 0;
    nTotalFace               = 0;
}

Base_Grid_Conn::~Base_Grid_Conn()
{
    delete left_cell_of_face;
    delete right_cell_of_face;
    delete face_location;
    delete node_number_of_each_face;
    delete face2node ;
    delete faceset;
    delete node_number_of_each_cell;
    delete cell2node;

    vector<CGNS_BCFace *>::iterator iter;
    for (iter = bcFace->begin(); iter != bcFace->end(); ++ iter)
    {
        delete *iter;
    }
    delete bcFace;
}

void Base_Grid_Conn::FreeFaceSet()
{
    delete faceset;
    faceset = 0;
}

void Base_Grid_Conn::FreeFace2Node()
{
    delete face2node;
    face2node = 0;
}

void Base_Grid_Conn::FreeCell2Node()
{
    delete cell2node;
    cell2node = 0;
}

void Base_Grid_Conn::FreeCellOfFace()
{
    delete left_cell_of_face;
    left_cell_of_face = 0;

    delete right_cell_of_face;
    right_cell_of_face = 0;
}

void Base_Grid_Conn::FreeNodeNumberOfEachFace()
{
    delete node_number_of_each_face;
    node_number_of_each_face = 0;
}

void Base_Grid_Conn::FreeNodeNumberOfEachCell()
{
    delete node_number_of_each_cell;
    node_number_of_each_cell = 0;
}

CGNSBase::CGNSBase()
{
    elementType       = 0;
    iStart            = 0;
    iEnd              = 0;
    connList          = 0;
    nBCElem           = 0;
    vcType            = 0;
    bcType            = 0;
    bcConnList        = 0;
    boundaryName      = 0;
    elementDataSize   = 0;
    mixedConnID       = 0;
    connectOffSet     = 0;
    isCartesian2D     = 0;
    isCartesian3D     = 0;
    cellStartIndex    = -1;
    nIFaceElem        = 0;
    interFaceType     = 0;
    interFaceConnList = 0;
    interFaceName     = 0;
    nIFaceRegion      = 0;
}

CGNSBase::~CGNSBase()
{
    delete [] boundaryName;             boundaryName          = nullptr;
    delete [] elementType;              elementType           = nullptr;
    delete [] iStart;                   iStart                = nullptr;
    delete [] iEnd;                     iEnd                  = nullptr;
    delete [] nBCElem;                  nBCElem               = nullptr;
    delete [] bcType;                   bcType                = nullptr;
    delete [] elementDataSize;          elementDataSize       = nullptr;
    delete [] nIFaceElem;               nIFaceElem            = nullptr;
    delete [] interFaceType;            interFaceType         = nullptr;
    delete [] interFaceElementType;     interFaceElementType  = nullptr;
    delete [] interFaceGridLocation;    interFaceGridLocation = nullptr;
    delete [] interFaceName;            interFaceName         = nullptr;

    if (connList)
    {
        for (int iSections = 0; iSections < nSections; ++ iSections)
        {
            delete [] connList[iSections];    connList[iSections] = nullptr;
        }
    }
    delete [] connList;    connList = nullptr;

    if (mixedConnID)
    {
        for (int iSections = 0; iSections < nSections; ++ iSections)
        {
            delete [] mixedConnID[iSections];    mixedConnID[iSections] = nullptr;
        }
    }
    delete [] mixedConnID;    mixedConnID = nullptr;

    if (bcConnList)
    {
        for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
        {
            delete [] bcConnList[iBCRegion];    bcConnList[iBCRegion] = nullptr;
        }
    }
    delete [] bcConnList;    bcConnList = nullptr;

    if (connectOffSet)
    {
        for (int iSection = 0; iSection < nSections; ++ iSection)
        {
            delete [] connectOffSet[iSection];    connectOffSet[iSection] = nullptr;
        }
    }
    delete [] connectOffSet;    connectOffSet = nullptr;

    if (interFaceConnList)
    {
        for (int iIFaceRegion = 0; iIFaceRegion < nIFaceRegion; ++ iIFaceRegion)
        {
            delete [] interFaceConnList[iIFaceRegion];    interFaceConnList[iIFaceRegion] = nullptr;
        }
    }
    delete [] interFaceConnList;    interFaceConnList = nullptr;
}

void CGNSBase::CreateElement(int nSections)
{
    this->nSections = nSections;

    connList        = new cgsize_t *[nSections];
    iStart          = new cgsize_t  [nSections];
    iEnd            = new cgsize_t  [nSections];
    elementType     = new int       [nSections];
    elementDataSize = new cgsize_t  [nSections];
    mixedConnID     = new cgsize_t *[nSections];
    connectOffSet   = new cgsize_t *[nSections];;
    cout << " nSections = " << nSections << "\n";
}

void CGNSBase::CreateBC(int nBCRegions)
{
    this->nBCRegions = nBCRegions;
    bcConnList          = new cgsize_t *[nBCRegions];
    nBCElem             = new int       [nBCRegions];
    bcType              = new int       [nBCRegions];
    boundaryElementType = new int       [nBCRegions];
    bcGridLocation      = new int       [nBCRegions];
    boundaryName        = new string    [nBCRegions];
}

void CGNSBase::CreateInterFace(int nIFaceREgions)
{
    this->nIFaceRegion = nIFaceREgions;
    interFaceConnList     = new cgsize_t * [nIFaceREgions];
    nIFaceElem            = new int        [nIFaceREgions];
    interFaceType         = new int        [nIFaceREgions];
    interFaceElementType  = new int        [nIFaceREgions];
    interFaceGridLocation = new int        [nIFaceREgions];
    interFaceName         = new string     [nIFaceREgions];
}

void CGNSBase::PreProcess(Base_Grid_Conn *gConn, BaseElement *baseElem)
{
    cout << "PreProcess ...\n";
    gConn->SetNTotalCell(nTotalCell);
    gConn->SetNTotalNode(nTotalNode);
    gConn->SetNBoundFace(nBoundFace);

    for (int iSection = 0; iSection < nSections; ++ iSection)
    {
        iStart[iSection] -= 1;
        iEnd  [iSection] -= 1;
    }

    cgsize_t *elementPoint = baseElem->GetElementPoint();
    cgsize_t size;
    for (int iSection = 0; iSection < nSections; ++ iSection)
    {
        if (MIXED == elementType[iSection])
        {
            size = elementDataSize[iSection];
        }
        else if(NGON_n == elementType[iSection] || NFACE_n == elementType[iSection])
        {
            size = elementDataSize[iSection];
        }
        else
        {
            cgsize_t nElem = iEnd[iSection] - iStart[iSection] + 1;
            cgsize_t pt    = elementPoint[elementType[iSection]];
            size = nElem * pt;
        }

        if (NFACE_n != elementType[iSection])
        {
            for (cgsize_t jNode = 0; jNode < size; ++ jNode)
            {
                cgsize_t connListABS = abs(connList[iSection][jNode]);
                connList[iSection][jNode] = (connList[iSection][jNode] / connListABS) * (connListABS - 1);
            }
        }
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        //! It is true for all situation.
        for (int jBCElem = 0; jBCElem < nBCElem[iBCRegion]; ++ jBCElem)
        {
            cgsize_t bcConnListABS = abs(bcConnList[iBCRegion][jBCElem]);
            bcConnList[iBCRegion][jBCElem] = (bcConnList[iBCRegion][jBCElem] / bcConnListABS) * (bcConnListABS - 1);
        }
    }

    for (int iIFaceRegion = 0; iIFaceRegion < nIFaceRegion; ++ iIFaceRegion)
    {
        //! It is true for all situation.
        for (int jIFaceElem = 0; jIFaceElem < nIFaceElem[iIFaceRegion]; ++ jIFaceElem)
        {
            cgsize_t interFaceConnListABS = abs(interFaceConnList[iIFaceRegion][jIFaceElem]);
            interFaceConnList[iIFaceRegion][jIFaceElem] = (interFaceConnList[iIFaceRegion][jIFaceElem] / interFaceConnListABS) * (interFaceConnListABS - 1);
        }
    }
}

void CGNSBase::GetElem(cgsize_t *elem_source, cgsize_t iCell, int ptnum, cgsize_t *&elem)
{
    elem = &elem_source[ptnum * iCell];
}

void CGNSBase::GetFace(Base_Grid_Conn *gconn, cgsize_t *elem, cgsize_t iCell, CElemFace *elemface)
{
    int face_list[100], face_now[100];
    int maxface = elemface->GetFaceNumber();

    vector<cgsize_t> single_face(100);
    single_face.resize(0);

    vector< vector<cgsize_t> > *face2node                = gconn->GetFace2Node();
    vector<cgsize_t>           *left_cell_of_face        = gconn->GetLeftCellOfFace();
    vector<cgsize_t>           *right_cell_of_face       = gconn->GetRightCellOfFace();
    vector<int>                *node_number_of_each_face = gconn->GetNodeNumberOfEachFace();
    vector<CGNS_BCFace *>      *bcFace                   = gconn->GetBCFace();
    vector<int>                *face_location            = gconn->GetFaceLocation();
    set<MultiIndex>            *faceSet                  = gconn->GetFaceSet();

    int *elemface_facept_num = elemface->GetFacePointNumber();
    int **elemface_face_id   = elemface->GetFaceID();

    for (int iFace = 0; iFace < maxface; ++ iFace)
    {
        int facept_num = elemface_facept_num[iFace];
        int *face_id   = elemface_face_id[iFace];

        MultiIndex faceIndex(facept_num, static_cast<int>(left_cell_of_face->size()));
        for (int j = 0; j < facept_num; ++ j)
        {
            face_list[j] = static_cast<int>(elem[face_id[j]]);
            face_now[j]  = face_list[j];
        }
        sort(face_list, face_list + facept_num);
        faceIndex.SetData(face_list);

        set<MultiIndex>::iterator iter = faceSet->find(faceIndex);

        if (iter == faceSet->end())
        {
            faceSet->insert(faceIndex);
            left_cell_of_face->push_back(iCell);
            right_cell_of_face->push_back(-1);

            CGNS_BCFace *boundaryFace = new CGNS_BCFace(static_cast<int>(bcFace->size()), 0, "");
            bcFace->push_back(boundaryFace);

            face_location->push_back(iFace);

            single_face.resize(facept_num);
            for (int j = 0; j < facept_num; ++ j)
            {
                single_face[j] = face_now[j];
            }

            face2node->push_back(single_face);
            node_number_of_each_face->push_back(facept_num);
        }
        else
        {
            (*right_cell_of_face)[iter->id] = iCell;
        }
    }
}

void CGNSBase::GetFaceInformationofFace2Node(Base_Grid_Conn *gConn, int sectionID)
{
    vector< vector<cgsize_t> > *face2Node            = gConn->GetFace2Node();
    vector<int>                *nodeNumberofEachFace = gConn->GetNodeNumberOfEachFace();
    vector<CGNS_BCFace*>       *bcFace               = gConn->GetBCFace();

    for (int iElement = 0; iElement < (iEnd[sectionID] - iStart[sectionID] + 1); ++ iElement)
    {
        vector<cgsize_t> iElementFace2Node;
        nodeNumberofEachFace->push_back(connectOffSet[sectionID][iElement + 1] - connectOffSet[sectionID][iElement]);
        for (int iNode = 0; iNode < (connectOffSet[sectionID][iElement + 1] - connectOffSet[sectionID][iElement]);++ iNode)
        {
            cgsize_t nodeID;
            nodeID = connList[sectionID][iNode+ connectOffSet[sectionID][iElement]];
            iElementFace2Node.push_back(nodeID);
        }
        face2Node->push_back(iElementFace2Node);

        CGNS_BCFace *boundaryFace = new CGNS_BCFace(static_cast<int>(bcFace->size()), 0, "");
        bcFace->push_back(boundaryFace);
    }

    int nTotalFaces = (gConn->GetNTotalFace()) + (iEnd[sectionID] - iStart[sectionID] + 1);
    gConn->SetNTotalFace(nTotalFaces);
}

void CGNSBase::GetFaceInformationofFace2Cell(Base_Grid_Conn *gConn, int sectionID)
{
    vector<cgsize_t> *leftCellOfFace   = gConn->GetLeftCellOfFace();
    vector<cgsize_t> *rightCellOfFace  = gConn->GetRightCellOfFace();
    int              nTotalFaces       = gConn->GetNTotalFace();
    leftCellOfFace->resize(nTotalFaces);
    rightCellOfFace->resize(nTotalFaces);

    for (int iFace = 0; iFace < nTotalFaces; ++ iFace)
    {
        (*leftCellOfFace) [iFace] = -1;
        (*rightCellOfFace)[iFace] = -1;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int iFace = 0; iFace < (connectOffSet[sectionID][iCell + 1] - connectOffSet[sectionID][iCell]); ++ iFace)
        {
            cgsize_t faceID = connList[sectionID][connectOffSet[sectionID][iCell] + iFace];
            if (faceID > 0)
            {
                (*leftCellOfFace)[abs(faceID) - 1] = iCell;
            }
            else
            {
                (*rightCellOfFace)[abs(faceID) - 1] = iCell;
            }
        }
    }
}

void CGNSBase::GetCartesian2DFace(Base_Grid_Conn *gConn, cgsize_t *elem, cgsize_t iCell, int nPoint)
{
    int faceList[100], faceNow[100];
    int maxFace = nPoint;
    vector<cgsize_t> singleFace(100);
    singleFace.resize(0);

    vector< vector<cgsize_t> > *face2Node            = gConn->GetFace2Node();
    vector<cgsize_t>           *leftCellOfFace       = gConn->GetLeftCellOfFace();
    vector<cgsize_t>           *rightCellOfFace      = gConn->GetRightCellOfFace();
    vector<int>                *nodeNumberOfEachFace = gConn->GetNodeNumberOfEachFace();
    vector<CGNS_BCFace *>      *bcFace               = gConn->GetBCFace();
    vector<int>                *faceLocation         = gConn->GetFaceLocation();
    set<MultiIndex>            *faceSet              = gConn->GetFaceSet();

    int *elemFaceFacePtNum = new int [nPoint];
    for (int iFace = 0; iFace < nPoint; ++ iFace)
    {
        elemFaceFacePtNum[iFace] = 2;
    }
    int **elemFaceFaceID = new int * [nPoint];
    int count = 0;
    for (int iFace = 0; iFace < nPoint; ++ iFace)
    {
        elemFaceFaceID[iFace] = new int[2];
        elemFaceFaceID[iFace][0] = count;
        ++ count;
        if (count < nPoint)
        {
            elemFaceFaceID[iFace][1] = count;
        }
        else
        {
            elemFaceFaceID[iFace][1] = 0;
        }
    }

    for (int iFace = 0; iFace < maxFace; ++ iFace)
    {
        int facePtNum = elemFaceFacePtNum[iFace];
        int *faceID   = elemFaceFaceID[iFace];

        MultiIndex faceIndex(facePtNum, static_cast<int>(leftCellOfFace->size()));
        for (int j = 0; j < facePtNum; ++ j)
        {
            faceList[j] = static_cast<int>(elem[faceID[j]]);
            faceNow [j] = faceList[j];
        }
        sort(faceList, faceList + facePtNum);
        faceIndex.SetData(faceList);

        set<MultiIndex>::iterator iter = faceSet->find(faceIndex);

        if (iter == faceSet->end())
        {
            faceSet->insert(faceIndex);
            leftCellOfFace->push_back(iCell);
            rightCellOfFace->push_back(-1);

            CGNS_BCFace *boundaryFace = new CGNS_BCFace(static_cast<int>(bcFace->size()), 0, "");
            bcFace->push_back(boundaryFace);

            faceLocation->push_back(iFace);

            singleFace.resize(facePtNum);
            for (int j = 0; j < facePtNum; ++ j)
            {
                singleFace[j] = faceNow[j];
            }

            face2Node->push_back(singleFace);
            nodeNumberOfEachFace->push_back(facePtNum);
        }
        else
        {
            (*rightCellOfFace)[iter->id] = iCell;
        }
    }
}

void CGNSBase::ProcessBCFaceCartesian2D(Base_Grid_Conn *gConn, cgsize_t iElement, int sectionID,
                                           BCType_t bcType, const string& BCName)
{
    cgsize_t *bcFace = 0;
    int facePtNum;

    GetCartesian2DBCFace(iElement, sectionID, facePtNum, bcFace);

    if (!MarkBCFace(gConn, bcFace, facePtNum, bcType, BCName))
    {
        ostringstream oss;
        oss << " ERROR in MarkBCFace : iElement = " << iElement << "\n";
        TK_Exit::ExceptionExit(oss.str());
    }
}

void CGNSBase::GetCartesian2DBCFace(cgsize_t iCell, int sectionID, int &facePtNum, cgsize_t *&bcFace)
{
    facePtNum = 2;
    bcFace = &connList[sectionID][iCell * facePtNum];
}

void CGNSBase::PostProcess(Base_Grid_Conn *gConn)
{
    cout << "  Boundary Faces Reordering ..." << endl;
    gConn->FreeFaceSet();
    vector<cgsize_t> *leftCellOfFace  = gConn->GetLeftCellOfFace();
    vector<cgsize_t> *rightCellOfFace = gConn->GetRightCellOfFace();
    cgsize_t         nTotalFace       = rightCellOfFace->size();
    gConn->SetNTotalFace(nTotalFace);
    nBoundFace = 0;
    for (cgsize_t iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if (-1 == (*rightCellOfFace)[iFace] || -1 == (*leftCellOfFace)[iFace])
        {
            ++ nBoundFace;
        }
    }
    gConn->SetNBoundFace(nBoundFace);

    vector< vector<cgsize_t> > *face2Node            = gConn->GetFace2Node();
    vector<int>                *nodeNumberOfEachFace = gConn->GetNodeNumberOfEachFace();
    vector<CGNS_BCFace *>      *bcFace               = gConn->GetBCFace();

    //! Reorder faces, make sure the boundary faces lie in the first part.
    cgsize_t jStart = nBoundFace;
    for (cgsize_t iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if ((*rightCellOfFace)[iFace] != -1 && (*leftCellOfFace)[iFace] != -1)
        {
            //! It needs to swap.
            for (cgsize_t jFace = jStart; jFace < nTotalFace; ++ jFace)
            {
                if (-1 == (*rightCellOfFace)[jFace] || -1 == (*leftCellOfFace)[jFace])
                {
                    //! Swap i,j.
                    SWAP((*leftCellOfFace)[iFace], (*leftCellOfFace)[jFace]);
                    SWAP((*rightCellOfFace)[iFace], (*rightCellOfFace)[jFace]);

                    //! Swap bctype && bcname only, do not swap faceID.
                    //! So, did not use SWAP((*bcFace)[iFace], (*bcFace)[jFace]) here.
                    ASSERT((*bcFace)[iFace]->GetFaceID() == iFace);
                    int iFaceBCType = (*bcFace)[iFace]->GetBCType();
                    int jFaceBCType = (*bcFace)[jFace]->GetBCType();
                    SWAP(iFaceBCType, jFaceBCType);
                    (*bcFace)[iFace]->SetBCType(iFaceBCType);
                    (*bcFace)[jFace]->SetBCType(jFaceBCType);
                    string iFaceBCName = (*bcFace)[iFace]->GetBCName();
                    string jFaceBCName = (*bcFace)[jFace]->GetBCName();
                    SWAP(iFaceBCName, jFaceBCName);
                    (*bcFace)[iFace]->SetBCName(iFaceBCName);
                    (*bcFace)[jFace]->SetBCName(jFaceBCName);    
                    SWAP((*face2Node)[iFace], (*face2Node)[jFace]);
                    SWAP((*nodeNumberOfEachFace)[iFace], (*nodeNumberOfEachFace)[jFace]);
                    jStart = jFace + 1;
                    break;
                }
            }
        }
    }
}

void CGNSBase::GetISection(cgsize_t iCell, int &isection)
{
    cgsize_t iabscell = iCell;
    isection = 0;
    for (int i = 0; i < nSections; ++ i)
    {
        cgsize_t ist = iStart[i];
        cgsize_t ied = iEnd  [i];
        if (ist <= iabscell && iabscell <= ied)
        {
            isection = i;
            break;
        }
    }
}

void CGNSBase::GetBCFace(BaseElement *base_elem, cgsize_t iCell, int &facept_num, cgsize_t *&bcface)
{
    cgsize_t iabscell = iCell;
    int isection = 0;

    GetISection(iCell, isection);
    
    ElementType_t etype = static_cast<ElementType_t>(elementType[isection]);
    cgsize_t *base_elem_elem_pt = base_elem->GetElementPoint();

    if (etype != MIXED)
    {
        facept_num = static_cast<int>(base_elem_elem_pt[etype]);
        GetElem(connList[isection], iabscell - iStart[isection], facept_num, bcface);
    }
    else
    {
        int iNode;
        iNode = static_cast<int>(mixedConnID[isection][iabscell - iStart[isection]]);
        etype = static_cast<ElementType_t>(connList[isection][iNode] + 1);
        facept_num = static_cast<int>(base_elem_elem_pt[etype]);
        GetBCFaceMixed(connList[isection], iabscell - iStart[isection], isection, bcface);
    }
}

void CGNSBase::GetBCFaceMixed(cgsize_t *elem_source, cgsize_t iCell, int isection, cgsize_t *&elem)
{
    elem = &elem_source[mixedConnID[isection][iCell] + 1];
}

void CGNSBase::GetElem(BaseElement *base_elem, cgsize_t iCell, int &etype, cgsize_t *&elem)
{
    int isection = 0;
    int elempt_num;

    GetISection(iCell, isection);
    etype = elementType[isection];

    cgsize_t *base_elem_elem_pt = base_elem->GetElementPoint();

    if (etype == MIXED)
    {
        etype = static_cast<ElementType_t>(connList[isection][mixedConnID[isection][iCell]] + 1);
        elem = &connList[isection][mixedConnID[isection][iCell] + 1];
    }
    else
    {
        elempt_num = static_cast<int>(base_elem_elem_pt[etype]);
        GetElem(connList[isection], iCell - iStart[isection], elempt_num, elem);
    }
}

void CGNSBase::ProcessBCFace(Base_Grid_Conn *gconn, BaseElement *base_elem, 
                             cgsize_t iElement, BCType_t bctype, const string &BCName)
{
    cgsize_t *bcface = 0;
    int facept_num;

    GetBCFace(base_elem, iElement, facept_num, bcface);
    
    if (!MarkBCFace(gconn, bcface, facept_num, bctype, BCName))
    {
        ostringstream oss;
        oss << " ERROR in MarkBCFace : iElement = " << iElement << "\n";
        TK_Exit::ExceptionExit(oss.str());
    }
}

bool CGNSBase::MarkBCFace(Base_Grid_Conn *gconn, cgsize_t *bcface, int facept_num, BCType_t bctype, const string &BCName)
{
    int face_list[100];
    MultiIndex idx(facept_num, 0);

    for (int i = 0; i < facept_num; ++ i)
    {
        face_list[i] = static_cast<int>(bcface[i]);
    }

    sort(face_list, face_list + facept_num);

    idx.SetData(face_list);

    set<MultiIndex> *faceset = gconn->GetFaceSet();

    set<MultiIndex>::iterator iter = faceset->find(idx);
    vector<CGNS_BCFace *> *bcFace = gconn->GetBCFace();

    if (iter == faceset->end())
    {
        cout << "impossible, boundary face should be in the total face\n";
        for (int i = 0; i < facept_num; ++ i)
        {
            cout << bcface[ i ] << " ";
        }
        cout << "\n";
        return false;
    }
    else
    {
        (*bcFace)[iter->id]->SetBCType(bctype);
        (*bcFace)[iter->id]->SetBCName(BCName);
        return true;
    }
}

void CGNSBase::ComputeGridConn(Base_Grid_Conn *gConn, BaseElement *baseElem)
{
    PreProcess(gConn, baseElem);
    ComputeElemFace(gConn, baseElem);
    ComputeBCFace(gConn, baseElem);
    ComputeCellToNode(gConn, baseElem);
    PostProcess(gConn);

    int gridReorder = PHSPACE::GlobalDataBase::GetIntParaFromDB("gridReorder");
    if (1 == gridReorder)
    {
        using namespace Reorder;
        //! Compute bandwidth
        printf("Before ReOrder...\n");
        CallBandwidthEvaluate(gConn);
        CallMetricLeftRightCellOfFace(gConn);
        //! Cell Reorder
        CallReorderCellLabel(gConn);
        //! FaceReorder
        CallReorderFaceLabel(gConn);
        //! Metrics for left_cell_of_face and right_cell_of_face
        printf("After ReOrder...\n");
        CallBandwidthEvaluate(gConn);
        CallMetricLeftRightCellOfFace(gConn);
    }
}

void CGNSBase::ComputeCellToNode(Base_Grid_Conn *gConn, BaseElement *baseElem)
{
#ifdef USE_VS2019
    //! if the type of the grid is Cartesian and the dimension is THREE_D,
    //! we don't need to compute CellToNode.
    if (isCartesian3D)
    {
        return;
    }
#endif
    int                        nSections             = this->GetNumberOfSections();
    int                        *baseCGNSElementType  = this->GetElementType();
    int                        nTotalCell            = static_cast<int>(this->GetNTotalCell());
    cgsize_t                   **connList            = this->GetElementConnectionList();
    vector<int>                *nodeNumberOfEachCell = gConn->GetNodeNumberOfEachCell();
    vector< vector<cgsize_t> > *cell2Node            = gConn->GetCell2Node();

    if (!isCartesian3D)
    {
        cell2Node->resize(nTotalCell);
        nodeNumberOfEachCell->resize(nTotalCell);
        cgsize_t *elementPoint = baseElem->GetElementPoint();
        bool     flag          = false;
        for (int iSection = 0; iSection < nSections; ++ iSection)
        {
            cgsize_t      *baseCGNSIStart = this->GetIndexOfStart();
            cgsize_t      *baseCGNSIEnd   = this->GetIndexOfEnd();
            ElementType_t eType           = static_cast<ElementType_t>(baseCGNSElementType[iSection]);

            int whetherBC = 0;
            for (int iBC = 0; iBC < nBCRegions; ++ iBC)
            {
                if (iStart[iSection] >= bcConnList[iBC][0]
                 && iEnd  [iSection] <= bcConnList[iBC][nBCElem[iBC] - 1])
                {
                    whetherBC = 1;
                    break;
                }
            }

            if (!whetherBC)
            {
                for (int iIFace = 0; iIFace < nIFaceRegion; ++ iIFace)
                {
                    if (iStart[iSection] >= interFaceConnList[iIFace][0]
                     && iEnd  [iSection] <= interFaceConnList[iIFace][nIFaceElem[iIFace] - 1])
                    {
                        whetherBC = 1;
                        break;
                    }
                }
            }

            if (whetherBC)
            {
                continue;
            }

            if (eType != MIXED)
            {
                int count = 0;
                int nPoints;
                for (cgsize_t iCell = baseCGNSIStart[iSection]; iCell <= baseCGNSIEnd[iSection]; ++ iCell)
                {
                    if (NGON_n != eType)
                    {
                        nPoints = static_cast<int>(elementPoint[eType]);
                    }
                    else
                    {
                        nPoints = connectOffSet[iSection][iCell - baseCGNSIStart[iSection] + 1] - connectOffSet[iSection][iCell - baseCGNSIStart[iSection]];
                    }
                    (*nodeNumberOfEachCell)[iCell - cellStartIndex] = nPoints;
                    vector<cgsize_t> singleCell;
                    singleCell.resize(0);
                    singleCell.reserve(nPoints);
                    for (cgsize_t jPoint = 0; jPoint < nPoints; ++ jPoint)
                    {
                        singleCell.push_back(connList[iSection][count]);
                        count ++;
                    }
                    (*cell2Node)[iCell - cellStartIndex] = singleCell;
                }
            }
            else
            {
                int count = 0;
                for (cgsize_t iCell = baseCGNSIStart[iSection]; iCell <= baseCGNSIEnd[iSection]; ++ iCell)
                {
                    int nPoints = static_cast<int>(elementPoint[connList[iSection][count] + 1]);
                    count ++;
                    (*nodeNumberOfEachCell)[iCell - cellStartIndex] = nPoints;
                    vector<cgsize_t> singleCell;
                    singleCell.resize(0);
                    singleCell.reserve(nPoints);
                    for (cgsize_t jPoint = 0; jPoint < nPoints; ++ jPoint)
                    {
                        singleCell.push_back(connList[iSection][count]);
                        count ++;
                    }
                    (*cell2Node)[iCell - cellStartIndex] = singleCell;
                }
            }
        }
    }

    if (connList)
    {
        for (int iSection = 0; iSection < nSections; ++ iSection)
        {
            DelPointer(connList[iSection]);
        }
    }
    DelPointer2(connList);

    if (bcConnList)
    {
        for (int iBCRegion = 0; iBCRegion < nBCRegions; ++iBCRegion)
        {
            DelPointer(bcConnList[iBCRegion]);
        }
    }
    DelPointer(bcConnList);

    DelPointer(iStart);
    DelPointer(iEnd);
    DelPointer(nBCElem);
    DelPointer(elementDataSize);
}

void CGNSBase::ComputeElemFace(Base_Grid_Conn *gConn, BaseElement *baseElem)
{
    std::cout << "Reconstruction faces ..." << endl;
    bool      intoBCSection    = false;
    cgsize_t  *elem            = 0;
    int       count            = 0;
    CElemFace *faceList        = baseElem->GetElementFace();
    cgsize_t  *elementPoint    = baseElem->GetElementPoint();

    for (int iSection = 0; iSection < nSections; ++iSection)
    {
        int whetherBC = 0;
        for (int iBC = 0; iBC < nBCRegions; ++iBC)
        {
            if (iStart[iSection] >= bcConnList[iBC][0]
             && iEnd  [iSection] <= bcConnList[iBC][nBCElem[iBC] - 1])
            {
                whetherBC = 1;
                break;
            }
        }

        if (!whetherBC)
        {
            for (int iIFace = 0; iIFace < nIFaceRegion; ++iIFace)
            {
                if (iStart[iSection] >= interFaceConnList[iIFace][0]
                 && iEnd  [iSection] <= interFaceConnList[iIFace][nIFaceElem[iIFace] - 1])
                {
                    whetherBC = 1;
                    break;
                }
            }
        }

        if (!whetherBC)
        {
            if (cellStartIndex < 0)
            {
                cellStartIndex = iStart[iSection];
            }
            else
            {
                cellStartIndex = MIN(cellStartIndex, iStart[iSection]);
            }
        }
    }

    for (int iSection = 0; iSection < nSections; ++ iSection)
    {
        int whetherBC = 0;
        for (int iBC = 0; iBC < nBCRegions; ++ iBC)
        {
            if (iStart[iSection] >= bcConnList[iBC][0]
             && iEnd  [iSection] <= bcConnList[iBC][nBCElem[iBC]-1])
            {
                whetherBC = 1;
                break;
            }
        }

        if (!whetherBC)
        {
            for (int iIFace = 0; iIFace < nIFaceRegion; ++ iIFace)
            {
                if (iStart[iSection] >= interFaceConnList[iIFace][0]
                 && iEnd  [iSection] <= interFaceConnList[iIFace][nIFaceElem[iIFace] - 1])
                {
                    whetherBC = 1;
                    break;
                }
            }
        }

        if (!whetherBC)
        {
            //! Only dump information for interior computational cells.
            //! Do not dump for BC sections.
            std::cout << "  Process_Grid: isection = " << iSection << " nSections = " << nSections << "\n";
            std::cout << "    iStart = " << iStart[iSection] << " " << ", iEnd = " << iEnd[iSection] << "\n";
        }
        else
        {
            std::cout << "  Get into BC sections from section " << iSection << endl;
        }

#ifdef USE_VS2019
        //! Construct the face topological information of 3D Cartesian grid.
        if (NGON_n == elementType[iSection] && isCartesian3D)
        {
            GetFaceInformationofFace2Node(gConn, iSection);
        }
        else if (NFACE_n == elementType[iSection] && isCartesian3D)
        {
            GetFaceInformationofFace2Cell(gConn, iSection);
        }
        else
#endif
        {
            //! Reconstruct NON-BC faces ONLY.
            if (whetherBC)
            {
                continue;
            }
            for (cgsize_t iCell = iStart[iSection]; iCell <= iEnd[iSection]; ++ iCell)
            {
                ElementType_t eType = static_cast<ElementType_t>(elementType[iSection]);
                ++ count;
                if (count % 200000 == 0) std::cout << "  iCell = " << count << "\n";

                if (MIXED == eType)
                {
                    eType = static_cast<ElementType_t>(connList[iSection][mixedConnID[iSection][iCell - iStart[iSection]]] + 1);
                    elem  = &connList[iSection][mixedConnID[iSection][iCell - iStart[iSection]] + 1];
                    GetFace(gConn, elem, iCell - cellStartIndex, &faceList[eType]);
                }
                else if (NGON_n == eType)
                {
                    int      numberPoint = connectOffSet[iSection][iCell - iStart[iSection] + 1] - connectOffSet[iSection][iCell - iStart[iSection]];
                    cgsize_t *elem       = 0;
                    elem = &connList[iSection][connectOffSet[iSection][iCell - iStart[iSection]]];
                    GetCartesian2DFace(gConn, elem, iCell - cellStartIndex, numberPoint);
                }
                else
                {
                    int pointNumber = static_cast<int>(elementPoint[eType]);
                    GetElem(connList[iSection], iCell - iStart[iSection], pointNumber, elem);
                    GetFace(gConn, elem, iCell - cellStartIndex, &faceList[eType]);
                }
            }
        }
    }
}

void CGNSBase::ComputeBCFace(Base_Grid_Conn *gConn, BaseElement *baseElem)
{
    std::cout << "Boundary Faces Reconstruction ..." << endl;
    vector<cgsize_t> &leftCellOfFace  = *gConn->GetLeftCellOfFace();
    vector<cgsize_t> &rightCellOfFace = *gConn->GetRightCellOfFace();
    cgsize_t         nTotalFace       = rightCellOfFace.size();
    cgsize_t         nBoundFace       = 0;
    for (cgsize_t iFace = 0; iFace < static_cast<cgsize_t>(nTotalFace); ++ iFace)
    {
        if (-1 == rightCellOfFace[iFace] || -1 == leftCellOfFace[iFace])
        {
            nBoundFace ++;
        }
    }
    cgsize_t maxCell = this->nTotalCell + nBoundFace;

    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        BCType_t       boundryType = static_cast<BCType_t> (bcType[iBCRegion]);
        GridLocation_t igr         = static_cast<GridLocation_t>(bcGridLocation[iBCRegion]);
        string         bcName      = "";
        if (boundaryName)
        {
            bcName = boundaryName[iBCRegion];
            if (boundryType < NofValidBCTypes)
            {
                std::cout << "  Boundary Condition : the " << iBCRegion + 1 
                    << " BCName = "   << bcName
                    << ", BCType = "  << BCTypeName[boundryType]
                    << "\n";
            }
            else
            {
                std::cout << "  Boundary Condition : the " << iBCRegion + 1 
                    << " BCName = " << bcName << "(SelfDefined)\n";
            }
        }

        if (PointRange == boundaryElementType[iBCRegion])
        {
            if (EdgeCenter == igr || FaceCenter == igr)
            {
                boundaryElementType[iBCRegion] = ElementRange;
            }
        }
        else if (PointList == boundaryElementType[iBCRegion])
        {
            if (EdgeCenter == igr || FaceCenter == igr)
            {
                boundaryElementType[iBCRegion] = ElementList;
            }
        }

        int sectionID;
        vector<CGNS_BCFace*> *boundaryFace = gConn->GetBCFace();
        if (ElementRange == boundaryElementType[iBCRegion])
        {
            for (int iSection = 0; iSection < nSections; ++ iSection)
            {
                if ((bcConnList[iBCRegion][0] <= iStart[iSection])
                 && (bcConnList[iBCRegion][1] >= iEnd  [iSection]))
                {
                    sectionID = iSection;
                    break;
                }
            }

#ifdef USE_VS2019
            if (isCartesian3D)
            {
                for (int iBCElement = 0; iBCElement < (iEnd[sectionID] - iStart[sectionID] + 1); ++ iBCElement)
                {
                    (*boundaryFace)[iBCElement + iStart[sectionID]]->SetBCType(boundryType);
                    (*boundaryFace)[iBCElement + iStart[sectionID]]->SetBCName(bcName);
                }
                continue;
            }
#endif
            if (bcConnList[iBCRegion][0] >= maxCell)
            {
                continue;
            }
            cgsize_t startElement = bcConnList[iBCRegion][0];
            cgsize_t endElement   = bcConnList[iBCRegion][1];
            for (cgsize_t iElement = startElement; iElement <= endElement; ++ iElement)
            {
                if (isCartesian2D)
                {
                    ProcessBCFaceCartesian2D(gConn, iElement - startElement, sectionID, boundryType, bcName);
                }
                else
                {
                    ProcessBCFace(gConn, baseElem, iElement, boundryType, bcName);
                }
            }
        }
        else if (ElementList == boundaryElementType[iBCRegion])
        {
            if (bcConnList[iBCRegion][0] >= maxCell)
            {
                continue;
            }
            for (int iBCElem = 0; iBCElem < nBCElem[iBCRegion]; ++ iBCElem)
            {
                ProcessBCFace(gConn, baseElem, bcConnList[iBCRegion][iBCElem], boundryType, bcName);
            }
        }
        //! This situation has not been encountered, we can not ensure the code is right! 
        else if (PointRange == boundaryElementType[iBCRegion]
              || PointList  == boundaryElementType[iBCRegion])
        {
            ProcessBCVertex(gConn, baseElem, boundaryElementType[iBCRegion],
                nBCElem[iBCRegion], bcConnList[iBCRegion], boundryType, bcName);
        }
        else
        {
            ostringstream oss;
            oss << "Error: this situation has not been considered, for points set type is " << PointSetTypeName[boundaryElementType[iBCRegion]];
            oss << " in BCRegions " << iBCRegion << endl;
            TK_Exit::ExceptionExit(oss);
        }
    }

    MarkBCInterface(gConn);
}

void CGNSBase::MarkBCInterface(Base_Grid_Conn *gconn)
{
    vector<cgsize_t> &right_cell_of_face = *gconn->GetRightCellOfFace();
    vector<cgsize_t> &left_cell_of_face  = *gconn->GetLeftCellOfFace();

    uint_t nTotalFace = right_cell_of_face.size();
    int nBoundFace = 0;
    for (cgsize_t iFace = 0; iFace <static_cast<cgsize_t>(nTotalFace); ++ iFace)
    {
        if (right_cell_of_face[iFace] == -1 || left_cell_of_face[iFace] == -1)
        {
            nBoundFace ++;
        }
    }
}

void CGNSBase::ProcessBCVertex(Base_Grid_Conn *gconn, BaseElement *base_elem, 
                               int bc_elem_set_type, int nBCElem, cgsize_t *bc_conn_list, int bctype, const string &BCName)
{
    set<cgsize_t> bc_vertex;
    if (bc_elem_set_type == PointRange)
    {
        //ProcessBCVertex(gconn, base_elem, bc_elem_set_type, nBCElem, bc_conn_list, bctype);
        for (int i = static_cast<int>(bc_conn_list[0]); i <= static_cast<int>(bc_conn_list[1]); ++ i)
        {
            bc_vertex.insert(i);
        }
    }
    else
    {
        for (int i = 0; i < nBCElem; ++ i)
        {
            bc_vertex.insert(bc_conn_list[i]);
        }
    }

    ProcessBCVertex(gconn, base_elem, bc_vertex, bctype, BCName);
}

void CGNSBase::GetFace(BaseElement *base_elem, cgsize_t iCell, int &facept_num, cgsize_t *&bcface)
{
    cgsize_t iabscell = iCell;
    int isection = 0;

    GetISection(iCell, isection);
    ElementType_t etype = static_cast<ElementType_t>(elementType[isection]);
    cgsize_t *base_elem_elem_pt = base_elem->GetElementPoint();
    facept_num = static_cast<int>(base_elem_elem_pt[etype]);
    GetElem(connList[isection], iabscell - iStart[isection], facept_num, bcface);
}

void CGNSBase::GetFace(cgsize_t *elem, int loc, int *bc_face, int &nfacept, CElemFace *elemface)
{
    int * elemface_facept_num = elemface->GetFacePointNumber();
    int **elemface_face_id    = elemface->GetFaceID();

    nfacept      = elemface_facept_num[loc];
    int *face_id = elemface_face_id[loc];

    for (int i = 0; i < nfacept; ++ i)
    {
        bc_face[i] = static_cast<int>(elem[face_id[i]]);
    }
}

bool CGNSBase::CheckBCFace(int *bc_face, int nfacept, set<cgsize_t> &bc_vertex)
{
    set<cgsize_t>::iterator iter;
    for (int i = 0; i < nfacept; ++ i)
    {
        iter = bc_vertex.find(bc_face[i]);
        if (iter == bc_vertex.end())
        {
            return false;
        }
    }
    return true;
}

void CGNSBase::ProcessBCVertex(Base_Grid_Conn *gconn, BaseElement *base_elem, 
                               set<cgsize_t> &bc_vertex, int bctype, const string &BCName)
{
    vector<cgsize_t> &left_cell_of_face  = *gconn->GetLeftCellOfFace();
    vector<cgsize_t> &right_cell_of_face = *gconn->GetRightCellOfFace();

    uint_t nTotalFace = right_cell_of_face.size();
    int nBoundFace = 0;
    for (cgsize_t iFace = 0; iFace < static_cast<cgsize_t>(nTotalFace); ++ iFace)
    {
        if (right_cell_of_face[iFace] == -1)
        {
            nBoundFace ++;
        }
    }

    cout << "ProcessBCVertex BCtype = " << bctype << " ........................................\n";
    cout << "All Candidate bc num is : " << nBoundFace << "\n";

    cgsize_t *elem;
    cgsize_t iCell;
    int etype, loc, nfacept;
    const int maxfacept = 100;
    int bc_face[maxfacept];

    vector<CGNS_BCFace *> *bcFace = gconn->GetBCFace();
    vector<int> *face_location = gconn->GetFaceLocation();

    CElemFace *face_list = base_elem->GetElementFace();

    for (cgsize_t iFace = 0; iFace < static_cast<cgsize_t>(nTotalFace); ++ iFace)
    {
        int old_bc_type = (*bcFace)[iFace]->GetBCType();
        if (right_cell_of_face[iFace] == - 1 && old_bc_type == 0)
        {
            iCell = left_cell_of_face[iFace];
            loc   = (*face_location)[iFace];

            GetElem(base_elem, iCell, etype, elem);
            GetFace(elem, loc, bc_face, nfacept, &face_list[etype]);

            if (CheckBCFace(bc_face, nfacept, bc_vertex))
            {
                (*bcFace)[iFace]->SetBCType(bctype);
                (*bcFace)[iFace]->SetBCName(BCName);
                if (bctype == 0)
                {
                    cout << "impossible bctype == 0 in ProcessBCVertex\n";
                }
            }
        }
    }
}

CGNSFactory::CGNSFactory(int nBlocks)
{
    this->nBlocks = nBlocks;
    baseConn = new Base_Grid_Conn *[nBlocks];
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        baseConn[iZone] = new Base_Grid_Conn();
    }

    rawGrid = new RawGrid *[nBlocks];
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        rawGrid[iZone] = new RawGrid();
    }

    baseCGNS = new CGNSBase *[nBlocks];
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        baseCGNS[iZone] = new CGNSBase();
    }

    baseElem  = new BaseElement();
    fantasyBC = new BaseBCType();
}

CGNSFactory::~CGNSFactory()
{
    delete baseElem;
    delete fantasyBC;

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        delete baseConn[iZone];
    }
    delete [] baseConn;    baseConn = nullptr;

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        delete rawGrid[iZone];
    }
    delete [] rawGrid;    rawGrid = nullptr;

    //for (int iZone = 0; iZone < nBlocks; ++iZone)
    //{
    //    delete baseCGNS[iZone];
    //}
    delete [] baseCGNS;    baseCGNS = nullptr;
}

void CGNSFactory::ConvertGrid2Fantasy()
{
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        std::cout << "Converting Grid to PHengLEI for Zone " << iZone << " of Total Blocks " << nBlocks << "\n";
        ComputeGridConn(this->GetCGNSBase(iZone), baseConn[iZone]);
    }

    DumpSelfDefinePartInformation();
}

void CGNSFactory::ComputeGridConn(CGNSBase *cgnsGrid, Base_Grid_Conn *gConn)
{
    cgnsGrid->ComputeGridConn(gConn, this->baseElem);
}

void CGNSFactory::ReLableSelfDefineBC()
{
    int countOfSelfDefine = 0;
    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        CGNSBase *cgnsGrid = this->GetCGNSBase(iBlock);

        int *bc_type = cgnsGrid->GetBCType();
        int nBCRegions = cgnsGrid->GetNumberOfBCRegions();

        for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
        {
            BCType_t bctype = static_cast<BCType_t>(bc_type[iBCRegion]);
            //! Since the number of custom boundary condition is 1, all of them are forced to add the maximum number NofValidBCTypes.
            if (bctype == BCTypeUserDefined)
            {
                bctype = static_cast<BCType_t> (NofValidBCTypes + countOfSelfDefine);
                ++ countOfSelfDefine;
            }

            if (bctype >= NofValidBCTypes)
            {
                bc_type[iBCRegion] = bctype;
            }
        }
    }
}

void CGNSFactory::DumpSelfDefinePartInformation()
{
    //! Collect all of the boundary condition.
    set<string> nameContainer;
    vector<int> bcType;
    vector<string> boundaryName;

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        CGNSBase *cgnsgrid = this->GetCGNSBase(iZone);

        int *bc_type = cgnsgrid->GetBCType();
        string *bc_name = cgnsgrid->GetBCName();
        int nBCRegions = cgnsgrid->GetNumberOfBCRegions();

        for (int iBC = 0; iBC < nBCRegions; ++ iBC)
        {
            int cgnsBCType = bc_type[iBC];
            int fantasyBCType = this->GetBaseBCType()->GetBCType(cgnsBCType);

            string name = bc_name[iBC];

            set<string> :: iterator nameIter;

            nameIter = nameContainer.find(name);
            //if (nameIter == nameContainer.end())
            if (nameIter == nameContainer.end() && fantasyBCType >= PHENGLEI::StartNumberOfBodyOfHyperFLOW && fantasyBCType <= PHENGLEI::ENDNumberOfBodyOfHyperFLOW)
            {
                boundaryName.push_back(name);
                bcType.push_back(fantasyBCType);
            }
        }
    }

    uint_t nPart = bcType.size();
    if (!nPart)
    {
        return;
    }

    //! Dump out.
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
    RDouble Lref = 1.0, Sref = 1.0, xref = 0.0, yref = 0.0, zref = 0.0;

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

void CGNS2UnsGrid(CGNSFactory *factoryCGNS, int iZone, UnstructGrid *grid)
{
    Base_Grid_Conn *gConn    = factoryCGNS->GetBaseGridConnection(iZone);
    RawGrid        *rawGrid  = factoryCGNS->GetRawGrid(iZone);
    CGNSBase       *cgnsGrid = factoryCGNS->GetCGNSBase(iZone);

    int nTotalCell = static_cast<int>(gConn->GetNTotalCell());
    int nTotalNode = static_cast<int>(gConn->GetNTotalNode());
    int nBoundFace = static_cast<int>(gConn->GetNBoundFace());

    grid->SetNTotalCell(nTotalCell);
    grid->SetNTotalNode(nTotalNode);
    grid->SetNBoundFace(nBoundFace);

    RDouble *x = new RDouble[nTotalNode];
    RDouble *y = new RDouble[nTotalNode];
    RDouble *z = new RDouble[nTotalNode];

    grid->SetX(x);
    grid->SetY(y);
    grid->SetZ(z);

    RDouble *xx = rawGrid->GetX();
    RDouble *yy = rawGrid->GetY();
    RDouble *zz = rawGrid->GetZ();

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        x[iNode] = xx[iNode];
        y[iNode] = yy[iNode];
        z[iNode] = zz[iNode];
    }
    delete [] xx;    xx = nullptr;
    delete [] yy;    yy = nullptr;
    delete [] zz;    zz = nullptr;

    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
    RDouble *gridTranslationVector = new RDouble[3]();
    if (GlobalDataBase::IsExist("gridTranslationVector", PHDOUBLE, 3))
    {
        GlobalDataBase::GetData("gridTranslationVector", gridTranslationVector, PHDOUBLE, 3);
    }

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        x[iNode] *= gridScaleFactor;
        x[iNode] += gridTranslationVector[0];
        y[iNode] *= gridScaleFactor;
        y[iNode] += gridTranslationVector[1];
        z[iNode] *= gridScaleFactor;
        z[iNode] += gridTranslationVector[2];
    }

    cout << "  Boundary condition setting ......\n";
    //! Set the boundary conditions.
    UnstructBCSet **bcrs = new UnstructBCSet *[nBoundFace];
    HyList<int> cgnsBCList;
    HyList<int> fantasyBCList;

    map<int, string>      fantasyBCTypeMap  = CreatFantasybctype();
    vector<CGNS_BCFace *> *bcFace           = gConn->GetBCFace();
    int                   nBCRegionUnstruct = cgnsGrid->GetNumberOfBCRegions();
    string                *bcRegionName     = cgnsGrid->GetBCName();

    for (cgsize_t iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        const string bcName        = (*bcFace)[iFace]->GetBCName();
        int          cgnsBCType    = (*bcFace)[iFace]->GetBCType();
        int          fantasyBCType = factoryCGNS->GetBaseBCType()->GetBCType(cgnsBCType);

        //!Collect faces which is NO_BOUNDARY_CONDITION into a BCRegion
        if (0 == fantasyBCType)
        {
            nBCRegionUnstruct ++;
            string *boundaryName = new string [nBCRegionUnstruct]();
            for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; ++ iBCRegionUnstruct)
            {
                if (iBCRegionUnstruct < nBCRegionUnstruct - 1)
                {
                    boundaryName[iBCRegionUnstruct] = bcRegionName[iBCRegionUnstruct];
                }
                else
                {
                    boundaryName[iBCRegionUnstruct] = bcName;
                }
            }
            cgnsGrid->SetBCName(boundaryName);
            break;
        }
    }
    string *boundaryName = cgnsGrid->GetBCName();

    grid->CreateUnstructBCSet(nBCRegionUnstruct);
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = new int[nBoundFace];
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; ++ iBCRegionUnstruct)
    {
        UnstructBC *unstructBC = new UnstructBC(iBCRegionUnstruct);
        unstructBCSet->SetBCRegion(iBCRegionUnstruct, unstructBC);
        unstructBC->SetBCName(boundaryName[iBCRegionUnstruct]);

        for (cgsize_t iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            const string bcName        = (*bcFace)[iFace]->GetBCName();
            int          cgnsBCType    = (*bcFace)[iFace]->GetBCType();
            int          fantasyBCType = factoryCGNS->GetBaseBCType()->GetBCType(cgnsBCType);
            cgnsBCList.insert(cgnsBCType);
            fantasyBCList.insert(fantasyBCType);
            if (bcName == unstructBC->GetBCName())
            {
                unstructBC->SetBCType(fantasyBCType);
                bcRegionIDofBCFace[iFace] = iBCRegionUnstruct;
                vector<int> *faceIndex = unstructBC->GetFaceIndex();
                faceIndex->push_back(static_cast<int>(iFace));
            }
        }
    }
    unstructBCSet->SetBCFaceInfo(bcRegionIDofBCFace);

    for (cgsize_t iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcrs[iFace] = new UnstructBCSet();
        const string &bcName       = (*bcFace)[iFace]->GetBCName();
        int          cgnsBCType    = (*bcFace)[iFace]->GetBCType();
        int          fantasyBCType = factoryCGNS->GetBaseBCType()->GetBCType(cgnsBCType);
        bcrs[iFace]->SetKey(fantasyBCType);
        bcrs[iFace]->SetBCName(bcName);
        cgnsBCList.insert(cgnsBCType);
        fantasyBCList.insert(fantasyBCType);
        delete (*bcFace)[iFace];
    }
    (*bcFace).clear();
    grid->SetBCRecord(bcrs);

    cout << "After converting CGNS grid,\n";
    cout << "Number of Original Boundary Condition: " << cgnsBCList.size() << "\n";
    for (int iBC = 0; iBC < cgnsBCList.size(); ++ iBC)
    {
        cout << "  The " << iBC + 1 << " CGNS BC Name is: " << BCTypeName[cgnsBCList[ iBC ]] << endl;
    }

    cout << "Number of Final Boundary Condition: " << fantasyBCList.size() << "\n";
    for (int iBC = 0; iBC < fantasyBCList.size(); ++ iBC)
    {
        cout << "  The " << iBC + 1 << " PHengLEI BC Name is: ";
        cout << fantasyBCTypeMap[fantasyBCList[iBC]] << "\n";
    }

    vector<cgsize_t> *gConnLeftCellOfFace       = gConn->GetLeftCellOfFace();
    vector<cgsize_t> *gConnRightCellOfFace      = gConn->GetRightCellOfFace();
    vector<int>      *gConnNodeNumberOfEachFace = gConn->GetNodeNumberOfEachFace();
    int              nTotalFace                 = static_cast<int>(gConnLeftCellOfFace->size());
    grid->SetNTotalFace(nTotalFace);
    cout << "Number of total faces: " << nTotalFace << "\n";
    int *nodeNumberOfEachFace = new int [nTotalFace];
    for (cgsize_t iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        nodeNumberOfEachFace[iFace] = (*gConnNodeNumberOfEachFace)[iFace];
    }
    grid->SetNodeNumberOfEachFace(nodeNumberOfEachFace);
    gConn->FreeNodeNumberOfEachFace();

    long long int face2NodeSize = 0;
    for (cgsize_t iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        face2NodeSize += nodeNumberOfEachFace[iFace];
    }
    int *face2Node = new int [face2NodeSize];
    grid->SetFace2Node(face2Node);
    vector< vector<cgsize_t> > *gFace2Node = gConn->GetFace2Node();
    long long int              count       = 0;
    for (cgsize_t iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        for (std::size_t jNode = 0; jNode < (*gFace2Node)[iFace].size(); ++ jNode)
        {
            face2Node[count] = static_cast<int>((*gFace2Node)[iFace][jNode]);
            ++ count;
        }
    }
    gConn->FreeFace2Node();

    int *leftCellOfFace  = new int[nTotalFace];
    int *rightCellOfFace = new int[nTotalFace];
    grid->SetLeftCellOfFace (leftCellOfFace);
    grid->SetRightCellOfFace(rightCellOfFace);
    for (cgsize_t iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        leftCellOfFace [iFace] = static_cast<int>((*gConnLeftCellOfFace)[iFace]);
        //! Set the label of ghostcell.
        rightCellOfFace[iFace] = static_cast<int>((*gConnRightCellOfFace)[iFace]);
    }
    gConn->FreeCellOfFace();

    //! if the type of the grid is Cartesian,don't need to compute Cell2Node
#ifdef USE_VS2019
    int isCartesian3D    = factoryCGNS->GetCGNSBase()->GetIsCartesian3D();
    if (!isCartesian3D)
    {
#endif
        int *nodeNumberOfEachCell = new int[nTotalCell];
        vector<int> *gConnNodeNumberOfEachCell = gConn->GetNodeNumberOfEachCell();
        for (cgsize_t iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            nodeNumberOfEachCell[iCell] = (*gConnNodeNumberOfEachCell)[iCell];
        }
        grid->SetNodeNumberOfEachCell(nodeNumberOfEachCell);
        gConn->FreeNodeNumberOfEachCell();

        long long int cell2NodeSize = 0;
        for (cgsize_t iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            cell2NodeSize += nodeNumberOfEachCell[iCell];
        }
        int *cell2Node = new int[cell2NodeSize];
        grid->SetCell2Node(cell2Node);
        vector< vector<cgsize_t> > *gCell2Node = gConn->GetCell2Node();
        count = 0;
        for (cgsize_t iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            for (int iNode = 0; iNode < nodeNumberOfEachCell[iCell]; ++ iNode)
            {
                cell2Node[count] = static_cast<int>((*gCell2Node)[iCell][iNode]);
                ++ count;
            }
        }
        gConn->FreeCell2Node();
#ifdef USE_VS2019
    }
#endif
    SimpleVC *volumeCondition = new SimpleVC();
    string vcName = cgnsGrid->GetVCName();
    volumeCondition->SetVCName(vcName);
    int vcType = cgnsGrid->GetVCType();
    volumeCondition->SetVCType(vcType);
    grid->SetVolumeCondition(volumeCondition);
}

}
