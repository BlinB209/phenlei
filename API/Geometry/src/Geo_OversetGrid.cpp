#include "Geo_OversetGrid.h"
#include "Geo_UnstructGrid.h"
#include "PHMatrix.h"
#include "PHIO.h"
#include "TK_Exit.h"
using namespace std;

namespace PHSPACE
{
Geo_OversetGrid::Geo_OversetGrid(Grid *grid_in, int rank, int nBlocks)
{
    grid                  = grid_in;
    this->rank            = rank;
    this->nBlocks         = nBlocks;
    xb                    = 0;
    yb                    = 0;
    zb                    = 0;
    numberOfBoundaryPoints = 0;

    oversetInfoProxy  = new OversetInfoProxy();
    //! Transfer oversetInfoProxy release right to grid.
    grid->SetOversetInfoProxy(oversetInfoProxy);
}

Geo_OversetGrid::~Geo_OversetGrid()
{
    delete [] xb;
    delete [] yb;
    delete [] zb;
}

bool Geo_OversetGrid::IsExteriorPoint(int flagOfPoint)
{
    return flagOfPoint == ExteriorPoint;
}

void Geo_OversetGrid::InitReceivingInformation(vector<int> &cellIndexRef1, vector<int> &zoneIndexRef2, vector<int> &cellIndexRef2)
{
    oversetInfoProxy->InitReceivingInformation(cellIndexRef1, zoneIndexRef2, cellIndexRef2);
}

void Geo_OversetGrid::PrepareSendingInformation(Geo_OversetGrid **oversetGrids)
{
    oversetInfoProxy->PrepareSendingInformation(oversetGrids);
}

void Geo_OversetGrid::PreProcess()
{
    Grid *grid = this->GetGrid();
    grid->ComputeMetrics();
}

void Geo_OversetGrid::PostProcess(Geo_OversetGrid **oversetGrids)
{
    oversetInfoProxy->PostProcess(oversetGrids);
}

void Geo_OversetGrid::DumpOversetGrid(fstream &file)
{
    DataContainer *cdata = new DataContainer();

    grid->Encode(cdata, 1);
    WriteFile(file, cdata);

    delete cdata;    cdata = nullptr;
}

bool InBox(RDouble &xb, RDouble &yb, RDouble &zb, RDouble *pmin, RDouble *pmax)
{
    if (xb < pmin[0] || xb > pmax[0]) return false;
    if (yb < pmin[1] || yb > pmax[1]) return false;
    if (zb < pmin[2] || zb > pmax[2]) return false;
    return true;
}

bool p_On_RHSide(RDouble &x, RDouble &y, RDouble &z, RDouble *p1, RDouble *p2)
{
    RDouble dx1,dy1,dz1,dx2,dy2,dz2;
    dx1 = x - p1[0];
    dy1 = y - p1[1];
    dz1 = z - p1[2];

    dx2 = p2[0] - p1[0];
    dy2 = p2[1] - p1[1];
    dz2 = p2[2] - p1[2];

    //RDouble sx = dy1 * dz2 - dy2 * dz1;
    //RDouble sy = dz1 * dx2 - dz2 * dx1;
    RDouble sz = dx1 * dy2 - dx2 * dy1;

    if (sz <= 0.0) return false;
    return true;
}

bool RightSide(RDouble &xm, RDouble &ym, RDouble &zm, RDouble &xfc, RDouble &yfc, RDouble &zfc, RDouble &xfn, RDouble &yfn, RDouble &zfn)
{
    RDouble dx = xm - xfc;
    RDouble dy = ym - yfc;
    RDouble dz = zm - zfc;

    RDouble result = dx * xfn + dy * yfn + dz * zfn;

    if (result > 0) return true;
    return false;
}

vector<string> GetTecplotTitle(const string &field_name)
{
    vector<string> title_tecplot;
    title_tecplot.push_back("title=\"Flow Fields of PHengLEI\"");
    title_tecplot.push_back("variables = ");
    title_tecplot.push_back("\"x\"");
    title_tecplot.push_back("\"y\"");
    title_tecplot.push_back("\"z\"");

    string str = "\"" + field_name + "\"";
    title_tecplot.push_back(str);

    return title_tecplot;
}

void DumpField(UnstructGrid *grid, fstream &file, RDouble *field, const string &name)
{
    int nTotalNode = grid->GetNTotalNode();

    int nplot = 1;
    std::ostringstream oss;

    RDouble ** qv = NewPointer2<RDouble>(nplot, nTotalNode);

    CompNodeVar(grid, qv[0], field);

    vector<string> title_tecplot = GetTecplotTitle(name);
    FieldVisualization(grid, oss, title_tecplot, qv, nplot);

    DelPointer2(qv);

    file << oss.str();
}

void Visualization(UnstructGrid *grid, RDouble *field, const string &filename, const string &fieldname)
{
    fstream file;
    file.open(filename.c_str(),ios_base::out|ios_base::trunc);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(filename);
    }
    DumpField(grid, file, field, fieldname);
    file.close();
    file.clear();
}


}
