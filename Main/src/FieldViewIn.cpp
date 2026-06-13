#include "FieldViewIn.h"
#include "Pre_FieldViewToPHengLEI_Unstruct.h"
#include "OverLappingGrid.h"
using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{
void FieldView2PHengLEI()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridFile  = GlobalDataBase::GetStrParaFromDB("out_gfile");

    //! Read field view grid and convert it's format.
    Pre_FieldViewToPHengLEI_Unstruct fieldViewGridConversion(fromGridFile);

    fieldViewGridConversion.Run();

    int numberOfBlocks = fieldViewGridConversion.GetNumberofBlocks();
    Grid **grids = fieldViewGridConversion.GetGrids();

    //! Construct block connection and dump out to file.
    ConstructGlobalInterfaceLink(grids, numberOfBlocks);
    DumpGrid(outGridFile, grids, numberOfBlocks);

    fieldViewGridConversion.WriteAdditionalInformation(outGridFile, grids);
    DumpCharacteristicBoundary(outGridFile, grids, numberOfBlocks);

    //! Clear the memory.
    fieldViewGridConversion.Clear();
}

}
