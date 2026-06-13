inline void Geo_LineTopo_Unstruct::SetLI_nLine(int nLine_in)
{
    this->nLine = nLine_in;
}

inline void Geo_LineTopo_Unstruct::SetLI_nCellsInLine(int *nCellsInLines_in)
{
    this->nCellsInLines = nCellsInLines_in;
}

inline void Geo_LineTopo_Unstruct::SetLI_CellOfLine(int **CellOfLine_in)
{
    this->CellOfLine = CellOfLine_in;
}

inline void Geo_LineTopo_Unstruct::SetLI_FaceOfLine(int **FaceOfLine_in)
{
    this->FaceOfLine = FaceOfLine_in;
}

inline void Geo_LineTopo_Unstruct::SetLine2Node(int *Line2NodeIn)
{
    this->Line2Node = Line2NodeIn;
}

inline void Geo_LineTopo_Unstruct::SetFace2Line(int** Face2LineIn)
{
    this->Face2Line = Face2LineIn;
}

inline void Geo_LineTopo_Unstruct::SetLI_LineOfCell(int *LineOfCell_in)
{
    this->LineOfCell = LineOfCell_in;
}

inline int Geo_LineTopo_Unstruct::GetLI_nLine() const
{
    return nLine;
}

inline int * Geo_LineTopo_Unstruct::GetLI_nCellsInLine() const
{
    return nCellsInLines;
}

inline int ** Geo_LineTopo_Unstruct::GetLI_CellOfLine() const
{
    return CellOfLine;
}

inline int ** Geo_LineTopo_Unstruct::GetLI_FaceOfLine() const
{
    return FaceOfLine;
}

inline int * Geo_LineTopo_Unstruct::GetLI_LineOfCell() const
{
    return LineOfCell;
}

