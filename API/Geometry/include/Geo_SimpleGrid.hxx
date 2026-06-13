inline void	SimpleGrid::SetX(RDouble *x)
{
    this->x = x;
}

inline void	SimpleGrid::SetY(RDouble *y)
{
    this->y = y;
}

inline void	SimpleGrid::SetZ(RDouble *z)
{
    this->z = z;
}

inline RDouble * SimpleGrid::GetX() const
{
    return x;
}

inline RDouble * SimpleGrid::GetY() const
{
    return y;
}

inline RDouble * SimpleGrid::GetZ() const
{
    return z;
}

inline void SimpleGrid::SetNTotalNode(int nTotalNode)
{
    this->nTotalNode = nTotalNode;
}

inline int SimpleGrid::GetNTotalNode() const
{
    return nTotalNode;
}