#include "Pre_Block_Struct.h"
#include "PHIO.h"
#include "PHMatrix.h"
#include "Math_BasisFunction.h"
#include "Glb_Dimension.h"

namespace PHSPACE
{
LIB_EXPORT Pre_BcFace_Struct::Pre_BcFace_Struct()
{
    simpleBlock  = 0;
    next         = 0;
    child        = 0;
    bcDirection  = -1;

    nodeStart.resize(3);
    nodeEnd.resize(3);
    axisLabel.resize(3);

    InitUnitMatrix();
}

LIB_EXPORT Pre_BcFace_Struct::~Pre_BcFace_Struct()
{
    Free();
}

LIB_EXPORT void Pre_BcFace_Struct::Free()
{
    if (child)
    {
        uint_t numberOfChildren = child->size();
        for (int i = 0; i < numberOfChildren; ++ i)
        {
            (* child)[i]->Free();
        }
        delete child;
    }

    delete next;
}

LIB_EXPORT void Pre_BcFace_Struct::ErectVertexMapping(int geometricDimension)
{
    vector< int > &s_st = this->nodeStart;
    vector< int > &s_ed = this->nodeEnd;
    vector< int > &s_nd = this->axisLabel;
    int nsDriection = this->GetDirection();

    GetSurfaceDirection(s_st, s_ed, s_nd, nsDriection, geometricDimension);

    for (int m = 0; m < 3; ++ m)
    {
        s_st[m] = abs(s_st[m]);
        s_ed[m] = abs(s_ed[m]);
    }
    
    vector< int > &t_st = this->next->GetNodeStart();
    vector< int > &t_ed = this->next->GetNodeEnd();
    vector< int > &t_nd = this->next->GetAxisLabel();
    int ntDriection = this->next->GetDirection();

    GetSurfaceDirection(t_st, t_ed, t_nd, ntDriection, geometricDimension);

    for (int m = 0; m < 3; ++ m)
    {
        t_st[m] = abs(t_st[m]);
        t_ed[m] = abs(t_ed[m]);
    }

    vector< int > sourceToTargetAxisMapping(3), targetToSourceAxisMapping(3), rate(3);
    for (int m = 0; m < 3; ++ m)
    {
        sourceToTargetAxisMapping[s_nd[m]] = t_nd[m];
        targetToSourceAxisMapping[t_nd[m]] = s_nd[m];
    }
    
    for (int m = 0; m < 3; ++ m)
    {
        if (s_st[m] > s_ed[m])
        {
            int n = sourceToTargetAxisMapping[m];
            PHSPACE::SWAP(s_st[m], s_ed[m]);
            PHSPACE::SWAP(t_st[n], t_ed[n]);
        }
    }

    for (int m = 0; m < 3; ++ m)
    {
        int n = sourceToTargetAxisMapping[m];
        rate[m] = SignFunction(t_ed[n] - t_st[n]);
    }
    
    vector< vector< int > > &frameOfAxes = this->next->GetFrameOfAxes();

    for (int n = 0; n < 3; ++ n)
    {
        int ms = targetToSourceAxisMapping[n];
        for (int m = 0; m < 3; ++ m)
        {
            frameOfAxes[n][m] = rate[m] * unitMatrix[m][ms];
        }
    }

    for (int n = 0; n < 3; ++ n)
    {
        frameOfAxes[n][3] = t_st[n];
        for (int m = 0; m < 3; ++ m)
        {
            frameOfAxes[n][3] -= frameOfAxes[n][m] * s_st[m];
        }
    }

    for (int m = 0; m < 3; ++ m)
    {
        if (t_st[m] > t_ed[m])
        {
            PHSPACE::SWAP(t_st[m], t_ed[m]);
        }
    }

    return;
}

LIB_EXPORT void Pre_BcFace_Struct::Normalize(int geometricDimension)
{
    if ((IsInterface(boundaryType)))
    {
        ErectVertexMapping(geometricDimension);
    }
    else
    {
        for (int m = 0; m < 3; ++ m)
        {
            if (nodeStart[m] > nodeEnd[m]) PHSPACE::SWAP(nodeStart[m], nodeEnd[m]);
        }
        //axisLabel[0] = GetSurfaceLabel(nodeStart, nodeEnd);
        axisLabel[0] = bcDirection;
    }
}

LIB_EXPORT void Pre_BcFace_Struct::ComputeLocalCoordinate()
{
    vector< int > &originalIndex = simpleBlock->GetOriginalIndex();
    for (int m = 0; m < 3; ++ m)
    {
        nodeStart[m] -= originalIndex[m];
        nodeEnd[m]   -= originalIndex[m];
    }
}

LIB_EXPORT void Pre_BcFace_Struct::OutputBoundaryInformation(StructBC *bcregion, int geometricDimension)
{
    int imin, imax, jmin, jmax, kmin, kmax, nbt;
    vector< int > s_st(nodeStart), s_ed(nodeEnd);
    
    if (IsInterface(boundaryType))
    {
        int patch = axisLabel[1];
        s_st[patch] = - abs(s_st[patch]);
        s_ed[patch] = - abs(s_ed[patch]);
    }

    imin = s_st[0];
    imax = s_ed[0];

    jmin = s_st[1];
    jmax = s_ed[1];

    kmin = s_st[2];
    kmax = s_ed[2];

    bcregion->SetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
    bcregion->SetBCType(boundaryType);
    bcregion->SetBCName(bcName);

    if (IsInterface(boundaryType))
    {
        bcregion->SetBCName("Connection");

        vector< int > t_st(3), t_ed(3);
        CalculateVertexIndexOnTargetBlock(next->GetFrameOfAxes(), nodeStart, t_st);
        CalculateVertexIndexOnTargetBlock(next->GetFrameOfAxes(), nodeEnd, t_ed);

        int patch = (next->GetAxisLabel())[1];
        t_st[patch] = - abs(t_st[patch]);
        t_ed[patch] = - abs(t_ed[patch]);

        imin = t_st[0];
        imax = t_ed[0];

        jmin = t_st[1];
        jmax = t_ed[1];

        kmin = t_st[2];
        kmax = t_ed[2];

        nbt  = next->GetSimpleBlock()->GetZoneIndex();

        bcregion->SetTargetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
        bcregion->SetTargetRegionBlock(nbt);
        if (GetDim() == THREE_D)
        {
            bcregion->ComputeRelativeParameters();
        }
    }
}

LIB_EXPORT void Pre_BcFace_Struct::OutputBoundaryInformation(fstream &file, int geometricDimension)
{
    int width = 8;
    int imin, imax, jmin, jmax, kmin, kmax, nbt;

    vector< int > s_st(nodeStart), s_ed(nodeEnd);

    if (IsInterface(boundaryType))
    {
        int patch = axisLabel[1];
        s_st[patch] = - abs(s_st[patch]);
        s_ed[patch] = - abs(s_ed[patch]);
    }

    imin = s_st[0];
    imax = s_ed[0];

    jmin = s_st[1];
    jmax = s_ed[1];

    kmin = s_st[2];
    kmax = s_ed[2];

    file << setw(width) << imin;
    file << setw(width) << imax;
    file << setw(width) << jmin;
    file << setw(width) << jmax;

    if (geometricDimension == 3)
    {
        file << setw(width) << kmin;
        file << setw(width) << kmax;
    }
    file << setw(width) << boundaryType << "\n";

    if (IsInterface(boundaryType))
    {
        vector< int > t_st(3), t_ed(3);
        CalculateVertexIndexOnTargetBlock(next->GetFrameOfAxes(), nodeStart, t_st);
        CalculateVertexIndexOnTargetBlock(next->GetFrameOfAxes(), nodeEnd, t_ed);

        int patch = (next->GetAxisLabel())[1];
        t_st[patch] = - abs(t_st[patch]);
        t_ed[patch] = - abs(t_ed[patch]);

        imin = t_st[0];
        imax = t_ed[0];

        jmin = t_st[1];
        jmax = t_ed[1];

        kmin = t_st[2];
        kmax = t_ed[2];

        nbt = next->GetSimpleBlock()->GetZoneIndex() + 1;

        file << setw(width) << imin;
        file << setw(width) << imax;
        file << setw(width) << jmin;
        file << setw(width) << jmax;
        if (geometricDimension == 3)
        {
            file << setw(width) << kmin;
            file << setw(width) << kmax;
        }
        file << setw(width) << nbt << "\n";
    }
}

LIB_EXPORT vector< Pre_BcFace_Struct * > Pre_BcFace_Struct::GetLeaves()
{
    vector< Pre_BcFace_Struct * > leaves;
    if (!child)
    {
        leaves.push_back(this);
        return leaves;
    }

    uint_t numberOfChildren = child->size();
    for (int i = 0; i < numberOfChildren; ++ i)
    {
        vector< Pre_BcFace_Struct * > cl = (*child)[i]->GetLeaves();

        uint_t numberOfLeaves = cl.size();
        for (int j = 0; j < numberOfLeaves; ++ j)
        {
            leaves.push_back(cl[j]);
        }
    }

    return leaves;
}

void Pre_BcFace_Struct::GetSurfaceDirection(vector< int > &st, vector< int > &ed, vector< int > &nd, int nodeDirection, int geometricDimension)
{
    vector< int >(3, -1).swap(nd);

    if (geometricDimension == 3)
    {
        for (int m = 0; m < geometricDimension; ++ m)
        {
            if (st[m] == ed[m])
            {
                nd[0] = m;
            }
            else if (st[m] < 0)
            {
                nd[1] = m;
            }
        }
    }
    else
    {
        for (int m = 0; m < geometricDimension; ++ m)
        {
            if (st[m] == ed[m])
            {
                nd[0] = m;
            }
            else
            {
                nd[1] = m;
            }
        }
    }

    nd[2] = 3 - abs(nd[0]) - abs(nd[1]);

    for (int m = 0; m < 3; ++ m)
    {
        if (nd[m] == -1)
        {
            cout << "Gridgen boundary file is incorrect!" << endl;
        }
    }

    return;
}

void Pre_BcFace_Struct::InitUnitMatrix()
{
    PHSPACE::AllocateVector(unitMatrix, 3, 3);
    for (int i = 0; i < 3; ++ i)
    {
        for (int j = 0; j < 3; ++ j)
        {
            if (i == j)
            {
                unitMatrix[i][j] = 1;
            }
            else
            {
                unitMatrix[i][j] = 0;
            }
        }
    }
    return;
}

LIB_EXPORT Pre_Patch_Struct::Pre_Patch_Struct()
{
    nodeStart.resize(3);
    nodeEnd.resize(3);
    axisLabel.resize(3);

    simpleBlock = 0;
    AllocateVector(frameOfAxes, 3, 4);
    targetBlockLabel = 0;
    bcDirection = 0;
}

LIB_EXPORT Pre_Patch_Struct::~Pre_Patch_Struct()
{

}

int Pre_Block_Struct::numberOfUnitCells = 1;

LIB_EXPORT Pre_Block_Struct::Pre_Block_Struct()
{
    parent = 0;
    child  = 0;

    originalIndex.resize(3, 0);
    nodeDimension.resize(3, 1);
}

LIB_EXPORT Pre_Block_Struct::~Pre_Block_Struct()
{
    Free();
}

LIB_EXPORT void Pre_Block_Struct::Split(RDouble numberOfCells, vector< int > &leftDimension, vector < int > &rightDimension, vector< int > &leftOriginalIndex, vector< int > &rightOriginalIndex)
{
    const int geometricDimension = 3;
    vector< int > faceCellNumber(geometricDimension);

    leftDimension  = this->nodeDimension;
    rightDimension = this->nodeDimension;

    leftOriginalIndex.resize(geometricDimension, 0);
    rightOriginalIndex.resize(geometricDimension, 0);
    
    for (int iDimension = 0; iDimension < geometricDimension; ++ iDimension)
    {
        faceCellNumber[iDimension] = max(nodeDimension[(iDimension+1)%geometricDimension] - 1, 1) * max(nodeDimension[(iDimension+2)%geometricDimension] - 1, 1);
    }

    int axisSelect = 0, maxNodeNumber = nodeDimension[0];
    for (int iDimension = 0; iDimension < geometricDimension; ++ iDimension)
    {
        if (IsPoleBoundaryFace(iDimension))
        {
             continue;
        }
  
        if (maxNodeNumber < nodeDimension[iDimension])
        {
            maxNodeNumber = nodeDimension[iDimension];
            axisSelect    = iDimension;
        }
    }

    RDouble dsplit = static_cast< RDouble >(numberOfCells) / faceCellNumber[axisSelect];
    int usplit = static_cast< int >(floor(dsplit + 0.5)) / numberOfUnitCells;
    int isplit = 1 + max(usplit, 1) * numberOfUnitCells;

    //For multigrid splitting, it's should be avoid to adjust the split position in the following way.
    CorrectSplitDimension(axisSelect, numberOfCells, geometricDimension, isplit);
    
    leftDimension[axisSelect]      = isplit;
    rightDimension[axisSelect]     = nodeDimension[axisSelect] - leftDimension[axisSelect] + 1;
    rightOriginalIndex[axisSelect] = isplit - 1;
    if (leftDimension[axisSelect] < 5 || rightDimension[axisSelect] < 5)
    {
        cout << "Too close splitting point from extreme point!\n";
    }
}

void Pre_Block_Struct::CorrectSplitDimension(int axisSelect, RDouble numberOfCells, int geometricDimension, int &isplit)
{

    //vector< int > leftDimension (geometricDimension);
    //vector< int > rightDimension(geometricDimension);
    //vector< int > faceCellNumber(geometricDimension);

    uint_t numberOfBoundaryFaces = boundaryFaceList.size();
    for (int iBoundaryFace = 0; iBoundaryFace < numberOfBoundaryFaces; ++ iBoundaryFace)
    {
        if (axisSelect != boundaryFaceList[iBoundaryFace]->GetDirection())
        {
            vector< int > s_st = boundaryFaceList[iBoundaryFace]->GetNodeStart();
            vector< int > s_ed = boundaryFaceList[iBoundaryFace]->GetNodeEnd();
            if (s_st[axisSelect] < isplit && isplit - s_st[axisSelect] < 3)
            {
                if (isplit == s_ed[axisSelect])
                {
                    continue;
                }
                cout << "Too close splitting point in surface boundary point!\n";
                isplit = s_st[axisSelect];
            }

            if (isplit < s_ed[axisSelect] && s_ed[axisSelect] - isplit < 3)
            {
                if (isplit == s_st[axisSelect])
                {
                    continue;
                }
                cout << "Too close splitting point in surface boundary point!\n";
                isplit = s_ed[axisSelect];
            }
        }
    }
}

bool Pre_Block_Struct::IsOverlapBoundaryFace(vector< int > &st, vector< int > &ed, int geometricDimension)
{
    int count = 0;
    for (int m = 0; m < 3; ++ m)
    {
        if (st[m] > ed[m])
        {
            return false;
        }

        if (st[m] == ed[m])
        {
            ++ count;
        }
    }

    if (geometricDimension == 3)
    {
        if (count < 2)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        if (count < 3)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

bool Pre_Block_Struct::IsPoleBoundaryFace(int iDimension)
{
    using namespace GRIDGEN_SPACE;
    uint_t numberOfBoundaryFaces = this->boundaryFaceList.size();
    for (int iBoundaryFace = 0; iBoundaryFace < numberOfBoundaryFaces; ++ iBoundaryFace)
    {
        int boundaryType = boundaryFaceList[iBoundaryFace]->GetBoundaryType();
        if (boundaryType / 10 == POLE && boundaryType - 71 == iDimension)
        {
            return true;
        }
    }

    return false;
}

LIB_EXPORT void Pre_Block_Struct::CreateBoundaryFace(int geometricDimension)
{
    uint_t numberOfBoundaryFaces = parent->boundaryFaceList.size();
    for (int i = 0; i < numberOfBoundaryFaces; ++ i)
    {
        Pre_BcFace_Struct *boundaryFace = parent->boundaryFaceList[i];
        CreateBoundaryFace(boundaryFace, geometricDimension);
    }
    return;
}

LIB_EXPORT void Pre_Block_Struct::NormalizeBoundaryFace(int geometricDimension)
{
    uint_t numberOfBoundaryFaces = this->boundaryFaceList.size();
    for (int i = 0; i < numberOfBoundaryFaces; ++ i)
    {
        boundaryFaceList[i]->Normalize(geometricDimension);
    }
    return;
}

LIB_EXPORT void Pre_Block_Struct::CalculateNumberOfUnitCells(int numberOfMultigrid)
{
    RDouble units = floor(pow(2.0, numberOfMultigrid - 1) + 0.5);
    Pre_Block_Struct::numberOfUnitCells = static_cast< int >(units);
    return;
}

LIB_EXPORT void Pre_Block_Struct::GenerateCoordinateSpace()
{
    int numberOfNodes = ni * nj * nk;

    x.resize(numberOfNodes);
    y.resize(numberOfNodes);
    z.resize(numberOfNodes);

    return;
}

Pre_BcFace_Struct * Pre_Block_Struct::GetOverlapBoundaryFace(vector< int > &st, vector< int > &ed, vector< int > &st_ref, vector< int > &ed_ref, int geometricDimension)
{
    int sourceSurfaceLabel = GetSurfaceLabel(st    , ed   );
    int targetSurfaceLabel = GetSurfaceLabel(st_ref, ed_ref);

    if (sourceSurfaceLabel != targetSurfaceLabel) return 0;

    if (st[sourceSurfaceLabel] != st_ref[sourceSurfaceLabel]) return 0;

    for (int m = 0; m < 3; ++ m)
    {
        int smin = min(st[m], ed[m]);
        int smax = max(st[m], ed[m]);

        if (smin > ed_ref[m] || smax < st_ref[m]) return 0;
    }

    vector< int > st_cross(3), ed_cross(3);

    for (int m = 0; m < 3; ++ m)
    {
        st_cross[m] = max(st[m], st_ref[m]);
        ed_cross[m] = min(ed[m], ed_ref[m]);
    }

    if (IsOverlapBoundaryFace(st_cross, ed_cross, geometricDimension))
    {
        Pre_BcFace_Struct *boundaryFace = new Pre_BcFace_Struct();
        boundaryFace->SetRegion(st_cross, ed_cross);
        return boundaryFace;
    }

    return 0;
}

void Pre_Block_Struct::PatchChild(vector< int > &t_st, vector< int > &t_ed, Pre_Block_Struct *sourceBlock, int geometricDimension)
{
    uint_t numberOfBoundaryFaces = boundaryFaceList.size();
    for (int i = 0; i < numberOfBoundaryFaces; ++ i)
    {
        if (boundaryFaceList[i]->GetBoundaryType() >= 0) continue;

        Pre_BcFace_Struct *boundaryFace = GetOverlapBoundaryFace(boundaryFaceList[i]->GetNodeStart(), boundaryFaceList[i]->GetNodeEnd(), t_st, t_ed, geometricDimension);

        boundaryFace = NULL;
        
        if (PatchChild(boundaryFaceList[i], t_st, t_ed, sourceBlock, geometricDimension)) return;
    }
}

bool Pre_Block_Struct::PatchChild(Pre_BcFace_Struct *boundaryFace, vector< int > &t_st, vector< int > &t_ed, Pre_Block_Struct *sourceBlock, int geometricDimension)
{
    Pre_BcFace_Struct *overlapBoundaryFace = GetOverlapBoundaryFace(boundaryFace->GetNodeStart(), boundaryFace->GetNodeEnd(), t_st, t_ed, geometricDimension);

    if (!overlapBoundaryFace) return false;

    if (boundaryFace->GetChild())
    {
        uint_t numberOfChildren = boundaryFace->GetChild()->size();
        for (int i = 0; i < numberOfChildren; ++ i)
        {
            Pre_BcFace_Struct *cbc = (* boundaryFace->GetChild())[i];
            if (PatchChild(cbc, t_st, t_ed, sourceBlock, geometricDimension)) return true;
        }
    }
    else
    {
        boundaryFace->SetChild(new vector< Pre_BcFace_Struct * >());
    }

    boundaryFace->GetChild()->push_back(overlapBoundaryFace);

    overlapBoundaryFace->SetBoundaryType(-1);
    overlapBoundaryFace->SetSimpleBlock(boundaryFace->GetSimpleBlock());
    overlapBoundaryFace->SetAxisLabel(boundaryFace->GetAxisLabel());

    overlapBoundaryFace->SetNext(new Pre_Patch_Struct());
    Pre_Patch_Struct *simplePatch = overlapBoundaryFace->GetNext();

    simplePatch->SetSimpleBlock(sourceBlock);
    simplePatch->SetNodeMapping(boundaryFace->GetNext()->GetFrameOfAxes());

    simplePatch->SetAxisLabel(boundaryFace->GetNext()->GetAxisLabel());

    vector< vector< int > > &frameOfAxes = simplePatch->GetFrameOfAxes();
    for (int m = 0; m < 3; ++ m)
    {
        frameOfAxes[m][3] -= sourceBlock->originalIndex[m];
    }

    return true;
}

void Pre_Block_Struct::PatchSplit(Pre_BcFace_Struct *boundaryFace, Pre_BcFace_Struct *overlapBoundaryFace, int geometricDimension)
{
    if (!IsInterface(overlapBoundaryFace->GetBoundaryType())) return;

    Pre_Patch_Struct *simplePatch = new Pre_Patch_Struct();
    overlapBoundaryFace->SetNext(simplePatch);

    simplePatch->SetSimpleBlock(boundaryFace->GetNext()->GetSimpleBlock());
    simplePatch->SetAxisLabel(boundaryFace->GetNext()->GetAxisLabel());
    simplePatch->SetNodeMapping(boundaryFace->GetNext()->GetFrameOfAxes());

    vector< vector< int > > &frameOfAxes = simplePatch->GetFrameOfAxes();
    for (int m = 0; m < 3; ++ m)
    {
        for (int n = 0; n < 3; ++ n)
        {
            frameOfAxes[m][3] += frameOfAxes[m][n] * this->originalIndex[n];
        }
    }

    vector< int > t_st(3), t_ed(3);
    PHSPACE::CalculateVertexIndexOnTargetBlock(simplePatch->GetFrameOfAxes(), overlapBoundaryFace->GetNodeStart(), t_st);
    PHSPACE::CalculateVertexIndexOnTargetBlock(simplePatch->GetFrameOfAxes(), overlapBoundaryFace->GetNodeEnd()  , t_ed);

    for (int m = 0; m < 3; ++ m)
    {
        if (t_st[m] > t_ed[m]) PHSPACE::SWAP(t_st[m], t_ed[m]);
    }

    Pre_Block_Struct *targetBlock = simplePatch->GetSimpleBlock();

    targetBlock->PatchChild(t_st, t_ed, this, geometricDimension);

    return;
}

void Pre_Block_Struct::GetAbsoluteCoordinate(vector< int > &st, vector< int > &ed)
{
    for (int m = 0; m < 3; ++ m)
    {
        st[m] = originalIndex[m] + 1;
        ed[m] = originalIndex[m] + nodeDimension[m];
    }
    return;
}

void Pre_Block_Struct::GetLocalCoordinate(vector< int > &st, vector< int > &ed, vector< int > &local_st, vector< int > &local_ed)
{
    for (int m = 0; m < 3; ++ m)
    {
        local_st[m] = st[m] - originalIndex[m];
        local_ed[m] = ed[m] - originalIndex[m];
    }
    return;
}

void Pre_Block_Struct::OutputBoundaryInformation(fstream &file, int geometricDimension)
{
    int width = 8;
    vector< Pre_BcFace_Struct * > leafBoundaryFace;

    uint_t numberOfBoundaryFaces = boundaryFaceList.size();
    for (int i = 0; i < numberOfBoundaryFaces; ++ i)
    {
        vector< Pre_BcFace_Struct * > cl = boundaryFaceList[i]->GetLeaves();
        uint_t numberOfLeaves = cl.size();
        for (int j = 0; j < numberOfLeaves; ++ j)
        {
            leafBoundaryFace.push_back(cl[j]);
        }
    }

    file << setw(width) << leafBoundaryFace.size() << endl;

    uint_t numberOfLeafBoundaryFaces = leafBoundaryFace.size();
    for (int i = 0; i < numberOfLeafBoundaryFaces; ++ i)
    {
        leafBoundaryFace[i]->OutputBoundaryInformation(file, geometricDimension);
    }

    return;
}

void Pre_Block_Struct::GetRootInformation(Pre_Block_Struct *&root, vector< int > &ref)
{
    for (int m = 0; m < 3; ++ m)
    {
        ref[m] += this->originalIndex[m];
    }

    if (this->parent)
    {
        parent->GetRootInformation(root, ref);
    }
    else
    {
        root = this;
    }

    return;
}

void Pre_Block_Struct::GenerateBoundaryFaceList(int numberOfBoundaryFaces)
{
    boundaryFaceList.resize(numberOfBoundaryFaces);

    return;
}

void Pre_Block_Struct::WriteGridCoordinate(fstream &file, vector< int > &nst, vector< int > &ned)
{
    vector< RDouble * > s(3, NULL);

    s[0] = this->GetX();
    s[1] = this->GetY();
    s[2] = this->GetZ();

    for (int m = 0; m < 3; ++ m)
    {
        for (int k = nst[2]-1; k < ned[2]; ++ k)
        {
            for (int j = nst[1]-1; j < ned[1]; ++ j)
            {
                for (int i = nst[0]-1; i < ned[0]; ++ i)
                {
                    int id = ni * nj * k + ni * j + i;

                    PHSPACE::PHWrite(file, s[m][id]);
                }
            }
        }
    }
    return;
}

void Pre_Block_Struct::Free()
{
    if (child)
    {
        uint_t numberOfChildren = child->size();
        for (int i = 0; i < numberOfChildren; ++ i)
        {
            (*child)[i]->Free();
        }
        delete child;
    }

    uint_t numberOfBoundaryFaces = boundaryFaceList.size();
    for (int i = 0; i < numberOfBoundaryFaces; ++ i)
    {
        delete boundaryFaceList[i];
    }

    return;
}

void Pre_Block_Struct::CreateBoundaryFace(Pre_BcFace_Struct *boundaryFace, int geometricDimension)
{
    if (!boundaryFace->GetChild())
    {
        vector< int > bst(3), bed(3);
        GetAbsoluteCoordinate(bst, bed);

        vector< int > ref_st(3), ref_ed(3);
        for (int m = 0; m < 6; ++ m)
        {
            PHSPACE::GetIJKRegion(bst, bed, ref_st, ref_ed, m);
            Pre_BcFace_Struct *overlapBoundaryFace = GetOverlapBoundaryFace(boundaryFace->GetNodeStart(), boundaryFace->GetNodeEnd(), ref_st, ref_ed, geometricDimension);
            GerneralSplit(boundaryFace, overlapBoundaryFace, geometricDimension);
        }
    }
    else
    {
        uint_t numberOfChildren = boundaryFace->GetChild()->size();
        for (int i = 0; i < numberOfChildren; ++ i)
        {
            CreateBoundaryFace((*boundaryFace->GetChild())[i], geometricDimension);
        }
    }
    return;
}

void Pre_Block_Struct::GerneralSplit(Pre_BcFace_Struct *boundaryFace, Pre_BcFace_Struct *overlapBoundaryFace, int geometricDimension)
{
    if (!overlapBoundaryFace) return;

    boundaryFaceList.push_back(overlapBoundaryFace);
    overlapBoundaryFace->SetBoundaryType(boundaryFace->GetBoundaryType());
    overlapBoundaryFace->SetBoundaryName(boundaryFace->GetBoundaryName());
    overlapBoundaryFace->SetSimpleBlock(this);

    overlapBoundaryFace->ComputeLocalCoordinate();
    overlapBoundaryFace->SetAxisLabel(boundaryFace->GetAxisLabel());

    if (overlapBoundaryFace->GetBoundaryType() < 0)
    {
        PatchSplit(boundaryFace, overlapBoundaryFace, geometricDimension);
    }

    return;
}

int GetSurfaceLabel(vector< int > &st, vector< int > &ed)
{
    for (int m = 0; m < 3; ++ m)
    {
        if (st[m] == ed[m]) return m;
    }
    return -1;
}

int SignFunction(int value)
{
    if (value < 0)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

void CalculateVertexIndexOnTargetBlock(vector< vector< int > > &frameOfAxes, vector< int > &p1, vector< int > &p2)
{
    for (int m = 0; m < 3; ++ m)
    {
        p2[m] = frameOfAxes[m][3];
    }

    for (int m = 0; m < 3; ++ m)
    {
        for (int n = 0; n < 3; ++ n)
        {
            p2[m] += frameOfAxes[m][n] * p1[n];
        }
    }
}

void GetIJKRegion(vector< int > &blockNodeStart, vector< int > &blockNodeEnd, vector< int > &faceNodeStart, vector< int > &faceNodeEnd, int localFaceLabel)
{
    faceNodeStart = blockNodeStart;
    faceNodeEnd   = blockNodeEnd;

    int axisLabel = localFaceLabel / 2, priority = localFaceLabel % 2;
    if (priority == 0)
    {
        faceNodeEnd[axisLabel] = faceNodeStart[axisLabel];
    }
    else
    {
        faceNodeStart[axisLabel] = faceNodeEnd[axisLabel];
    }

    return;
}

}