#include "Geo_LineTopo_Unstruct.h"

namespace PHSPACE
{
    LIB_EXPORT Geo_LineTopo_Unstruct::Geo_LineTopo_Unstruct()
    {
        nLine = 0;
        starIndex     = 0;
        endIndex      = 0;
        nCellsInLines = NULL;
        CellOfLine = NULL;
        FaceOfLine = NULL;
        LineOfCell = NULL;
        Line2Node     = NULL;
        Face2Line     = NULL;
    }

    LIB_EXPORT Geo_LineTopo_Unstruct::~Geo_LineTopo_Unstruct()
    {
        if (nCellsInLines != NULL)
        {
            DelPointer(nCellsInLines);
        }

        if (CellOfLine != NULL)
        {
            DelPointer2(CellOfLine);
        }
        
        if (FaceOfLine != NULL)
        {
            DelPointer2(FaceOfLine);
        }

        if (LineOfCell != NULL)
        {
            DelPointer(LineOfCell);
        }

        if (Line2Node != NULL)
        {
            DelPointer(Line2Node);
        }

        if (Face2Line != NULL)
        {
            DelPointer2(Face2Line);
        }
    }

}