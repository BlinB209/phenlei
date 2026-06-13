inline void Pre_BcFace_Struct::SetRegion(vector<int> &nodeStart, vector<int> &nodeEnd)
{
    this->nodeStart = nodeStart;
    this->nodeEnd = nodeEnd; 
}

inline void Pre_BcFace_Struct::SetAxisLabel(vector<int> &axisLabel)
{
    this->axisLabel = axisLabel;
}

inline void Pre_BcFace_Struct::SetSimpleBlock(Pre_Block_Struct *simpleBlock)
{
    this->simpleBlock = simpleBlock;
}

inline void Pre_BcFace_Struct::SetNext(Pre_Patch_Struct *next)
{
    this->next = next;
}

inline void Pre_BcFace_Struct::SetChild(vector<Pre_BcFace_Struct *> *child)
{
    this->child = child;
}

inline void Pre_BcFace_Struct::SetBoundaryType(int boundaryType)
{
    this->boundaryType = boundaryType;
}

inline void Pre_BcFace_Struct::SetBoundaryName(string bcName)
{
    this->bcName = bcName;
}

inline void Pre_BcFace_Struct::SetFaceMatchingTargetDirIndex(int *dir3dIn)
{
    for (int i = 0; i < 3; i ++)
    {
        s_dir3d[i] = dir3dIn[i];
    }
}

inline vector<int> & Pre_BcFace_Struct::GetNodeStart()
{
    return this->nodeStart;
}

inline vector<int> & Pre_BcFace_Struct::GetNodeEnd()
{
    return this->nodeEnd;
}

inline vector<int> & Pre_BcFace_Struct::GetAxisLabel()
{
    return this->axisLabel;
}

inline Pre_Block_Struct * Pre_BcFace_Struct::GetSimpleBlock()
{
    return this->simpleBlock;
}

inline Pre_Patch_Struct * Pre_BcFace_Struct::GetNext()
{
    return this->next;
}

inline int Pre_BcFace_Struct::GetBoundaryType()
{
    return this->boundaryType;
}

inline string Pre_BcFace_Struct::GetBoundaryName()
{
    return this->bcName;
}

inline int * Pre_BcFace_Struct::GetFaceMatchingTargetDirIndex()
{
    return s_dir3d;
}

inline vector<Pre_BcFace_Struct *> * Pre_BcFace_Struct::GetChild()
{
    return this->child;
}

inline void Pre_Patch_Struct::SetRegion(vector<int> &nodeStart, vector<int> &nodeEnd)
{
    this->nodeStart = nodeStart;
    this->nodeEnd   = nodeEnd;
}

inline void Pre_Patch_Struct::SetNodeMapping(vector<vector<int> > &frameOfAxes)
{
    this->frameOfAxes = frameOfAxes;
}

inline void Pre_Patch_Struct::SetAxisLabel(vector<int> &axisLabel)
{
    this->axisLabel = axisLabel;
}

inline void Pre_Patch_Struct::SetSimpleBlock(Pre_Block_Struct *simpleBlock)
{
    this->simpleBlock = simpleBlock;
}

inline void Pre_Patch_Struct::SetTargetBlockLabel(int i)
{
    this->targetBlockLabel = i;
}

inline vector<int> & Pre_Patch_Struct::GetNodeStart()
{
    return this->nodeStart;
}

inline vector<int> & Pre_Patch_Struct::GetNodeEnd()
{
    return this->nodeEnd;
}

inline vector<int> & Pre_Patch_Struct::GetAxisLabel()
{
    return this->axisLabel;
}

inline Pre_Block_Struct * Pre_Patch_Struct::GetSimpleBlock()
{
    return this->simpleBlock;
}

inline vector<vector<int> > & Pre_Patch_Struct::GetFrameOfAxes()
{
    return this->frameOfAxes;
}

inline int Pre_Block_Struct::GetNumberOfUnitCells()
{
    return Pre_Block_Struct::numberOfUnitCells;
}

inline int Pre_Block_Struct::GetNumberOfNodes()
{
    return nodeDimension[0] * nodeDimension[1] * nodeDimension[2];
}

inline int Pre_Block_Struct::GetNumberOfCells()
{
    return (nodeDimension[0] - 1) * (nodeDimension[1] - 1) * max((nodeDimension[2] - 1), 1);
}

inline void Pre_Block_Struct::AddBoundaryFace(Pre_BcFace_Struct *boundaryFace)
{
    boundaryFaceList.push_back(boundaryFace);
}

inline void Pre_Block_Struct::SetZoneIndex(int iZone)
{
    this->zoneIndex = iZone;
}

inline void Pre_Block_Struct::SetProcessorIndex(int iProcessor)
{
    this->processorIndex = iProcessor;
}

inline void Pre_Block_Struct::SetOriginalZoneIndex(int iZone)
{
    this->originalZoneIndex = iZone;
}

inline void Pre_Block_Struct::SetParentBlock(Pre_Block_Struct *parent)
{
    this->parent = parent;
}

inline void Pre_Block_Struct::SetBoundaryFace(int i, Pre_BcFace_Struct *boundaryFace)
{
    boundaryFaceList[i] = boundaryFace;
}

inline int Pre_Block_Struct::GetZoneIndex()
{
    return zoneIndex;
}

inline int Pre_Block_Struct::GetOriginalZoneIndex()
{
    return originalZoneIndex;
}

inline Pre_BcFace_Struct * Pre_Block_Struct::GetBoundaryFace(int i)
{
    return boundaryFaceList[i];
}

inline void Pre_Block_Struct::SetOriginalIndex(vector<int> &od)
{
    originalIndex = od;
}
inline void Pre_Block_Struct::SetNodeDimension(vector<int> &nd)
{
    nodeDimension = nd;
}

inline vector<int> & Pre_Block_Struct::GetOriginalIndex()
{
    return originalIndex;
}

inline vector<int> & Pre_Block_Struct::GetNodeDimension()
{
    return nodeDimension;
}

inline int Pre_Block_Struct::GetNI() const
{
    return ni;
}

inline int Pre_Block_Struct::GetNJ() const
{
    return nj;
}

inline int Pre_Block_Struct::GetNK() const
{
    return nk;
}

inline void Pre_Block_Struct::SetNI(int ni)
{
    this->ni = ni;
}

inline void Pre_Block_Struct::SetNJ(int nj)
{
    this->nj = nj;
}

inline void Pre_Block_Struct::SetNK(int nk)
{
    this->nk = nk;
}

inline RDouble * Pre_Block_Struct::GetX()
{
    return &x[0];
}

inline RDouble * Pre_Block_Struct::GetY()
{
    return &y[0];
}

inline RDouble * Pre_Block_Struct::GetZ()
{
    return &z[0];
}