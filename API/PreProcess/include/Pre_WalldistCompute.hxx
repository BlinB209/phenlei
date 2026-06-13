inline void Pre_WalldistCompute::RotateAxis(const RDouble *oldR, RDouble *newR)
{
    if (walldistComputeMethod != FAST_METHOD)
    {
        PHSPACE::SetField(newR, oldR, 3);
        return;
    }

    switch (rotateAxis)
    {
    case 0:
        PHSPACE::SetField(newR, oldR, 3);
        return;
        break;
    case 1:
        RotateAxisX(oldR, newR, rotateDegree);
        break;
    case 2:
        RotateAxisY(oldR, newR, rotateDegree);
        break;
    case 3:
        RotateAxisZ(oldR, newR, rotateDegree);
        break;
    default:
        PHSPACE::SetField(newR, oldR, 3);
        break;
    }
}

inline std::size_t WallStructure::GetNumberOfWallFaces()
{
    return xfc.size();
}

inline std::size_t WallStructure::GetNumberOfWallPoints()
{
    return x.size();
}

inline WallStructure::value_type & WallStructure::GetX()
{
    return x;
}

inline WallStructure::value_type & WallStructure::GetY()
{
    return y;
}

inline WallStructure::value_type & WallStructure::GetZ()
{
    return z;
}

inline WallStructure::value_type & WallStructure::GetXFaceCenter()
{
    return xfc;
}

inline WallStructure::value_type & WallStructure::GetYFaceCenter()
{
    return yfc;
}

inline WallStructure::value_type & WallStructure::GetZFaceCenter()
{
    return zfc;
}

inline WallStructure::value_type & WallStructure::GetXFaceNormal()
{
    return xfn;
}

inline WallStructure::value_type & WallStructure::GetYFaceNormal()
{
    return yfn;
}

inline WallStructure::value_type & WallStructure::GetZFaceNormal()
{
    return zfn;
}

inline int * WallStructure::GetWallFace2Node()
{
    return wallFace2Node;
}

inline int * WallStructure::GetnPointPerFace()
{
    return nPointPerFace;
}

inline int * WallStructure::GetNodePosition()
{
    return nodePosition;
}

inline void WallStructure::SetWallFace2Node(int *wallFace2Node)
{
    this->wallFace2Node = wallFace2Node;
}

inline void WallStructure::SetnPointPerFace(int *nPointPerFace)
{
    this->nPointPerFace = nPointPerFace;
}

inline void WallStructure::SetNodePosition(int *nodePosition)
{
    this->nodePosition = nodePosition;
}

inline RDouble * WallStructure::GetBox()
{
    return this->box;
}

inline RDouble WallStructure::GetDistance()
{
    return distance;
}

inline void WallStructure::SetDistance(RDouble dist)
{
    distance = dist;
}

inline int WallStructure::GetIBlock()
{
    return iBlock;
}

inline void WallStructure::SetIBlock(int iBlockIn)
{
    this->iBlock = iBlockIn;
}