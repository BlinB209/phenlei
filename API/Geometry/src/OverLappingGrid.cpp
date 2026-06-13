#include "OverLappingGrid.h"
#include "Geo_OversetGrid.h"
#include "TK_Exit.h"
using namespace std;

namespace PHSPACE
{
void DumpFlags(Geo_OversetGrid **oversetGrids, int nBlocks, const string &filename)
{
    fstream file;
    file.open(filename.c_str(),ios_base::out);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(filename);
    }

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        file << " iZone = " << iZone << " nBlocks = " << nBlocks << "\n";
        oversetGrids[iZone]->DumpFlags(file);
        //file << "\n";
    }
    file.close();
    file.clear();
}

void DumpOversetGrid(Geo_OversetGrid **oversetGrids, int nBlocks, const string &filename)
{
    fstream file;
    file.open(filename.c_str(), ios_base::out|ios_base::binary);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(filename);
    }

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        oversetGrids[iZone]->DumpOversetGrid(file);
    }
    file.close();
    file.clear();
}

}
