inline int StructBC::GetRegionBlock() const
{
    return nbs;
}

inline void StructBC::SetRegionBlock(int currentSourceZoneID)
{
    this->nbs = currentSourceZoneID;
}

inline void StructBC::SetFaceMatchingTargetDirIndex(int *dir3dIn)
{
    for (int i = 0; i < 3; i ++)
    {
        s_dir3d[i] = dir3dIn[i];
    }
}

inline int * StructBC::GetStam()
{
    return stam;
}

inline int * StructBC::GetRate()
{
    return rate;
}

inline int * StructBC::GetFaceDirectionIndex()
{
    return s_lr3d;
}

inline int * StructBC::GetFaceMatchingTargetDirIndex()
{
    return s_dir3d;
}

inline int * StructBC::GetFaceDirectionIndexOfTargetBlock()
{
    return t_lr3d;
}

inline void StructBC::SetIJKRegion(int ist, int ied, int jst, int jed, int kst, int ked)
{
    s_st[0] = ist;
    s_ed[0] = ied;
    s_st[1] = jst;
    s_ed[1] = jed;
    s_st[2] = kst;
    s_ed[2] = ked;
}

inline void StructBC::GetIJKRegion(int &ist, int &ied, int &jst, int &jed, int &kst, int &ked)
{
    ist = s_st[0];
    ied = s_ed[0];
    jst = s_st[1];
    jed = s_ed[1];
    kst = s_st[2];
    ked = s_ed[2];
}

inline void StructBC::SetTargetIJKRegion(int ist, int ied, int jst, int jed, int kst, int ked)
{
    t_st[0] = ist;
    t_ed[0] = ied;
    t_st[1] = jst;
    t_ed[1] = jed;
    t_st[2] = kst;
    t_ed[2] = ked;
}

inline int * StructBC::GetStartPoint()
{
    return s_st;
}

inline int StructBC::GetStartPoint(int m)
{
    return s_st[m];
}

inline int StructBC::GetTargetStart(int m)
{
    return t_st[m];
}

inline int * StructBC::GetTargetStart()
{
    return t_st;
}

inline int StructBC::GetTargetEnd(int m)
{
    return t_ed[m];
}

inline int * StructBC::GetTargetEnd()
{
    return t_ed;
}

inline int * StructBC::GetEndPoint(void)
{
    return s_ed;
}

inline int * StructBC::GetInkInfo()
{
    return lnkInfo;
}

inline int StructBC::GetEndPoint(int m)
{
    return s_ed[m];
}

inline int StructBC::GetTargetRegionBlock() const
{
    return nbt;
}

inline void StructBC::SetTargetRegionBlock(int targetZoneID)
{
    this->nbt = targetZoneID;
}

inline int StructBC::GetFaceLeftOrRightIndex() const
{
    return s_lr;
}

inline void StructBC::SetFaceLeftOrRightIndex(int s_lr)
{
    this->s_lr = s_lr;
}

inline int StructBC::GetFaceLeftOrRightIndexOfTargetBlock() const
{
    return t_lr;
}

inline void StructBC::SetFaceLeftOrRightIndexOfTargetBlock(int t_lr)
{
    this->t_lr = t_lr;
}

inline int StructBC::GetFaceDirection() const
{
    return s_nd;
}

inline void StructBC::SetFaceDirection(int s_nd)
{
    this->s_nd = s_nd;
}

inline int StructBC::GetFaceDirectionOfTargetBlock() const
{
    return t_nd;
}

inline void StructBC::SetFaceDirectionOfTargetBlock(int t_nd)
{
    this->t_nd = t_nd;
}

inline void StructBC::GetInsideCellIndex(int i, int j, int k, int &is, int &js, int &ks, int nss)
{
    is = i - s_lr3d[0] * (nss - 1);
    js = j - s_lr3d[1] * (nss - 1);
    ks = k - s_lr3d[2] * (nss - 1);
}

inline void StructBC::GetGhostCellIndex(int i, int j, int k, int &it, int &jt, int &kt, int ntt)
{
    it = i + s_lr3d[0] * ntt;
    jt = j + s_lr3d[1] * ntt;
    kt = k + s_lr3d[2] * ntt;
}

inline void StructBC::GetBoundaryFaceIndex(int i, int j, int k, int &ib, int &jb, int &kb)
{
    ib = i;
    jb = j;
    kb = k;
    
    if (s_lr3d[0] == 1)
    {
        ib = i + 1;
    }
    else if (s_lr3d[1] == 1)
    {
        jb = j + 1;
    }
    else if (s_lr3d[2] == 1)
    {
        kb = k + 1;
    }
}

inline void StructBCSet::CreateBCRegion(uint_t nBCRegion)
{
    bcRegions = new vector<StructBC *>(nBCRegion);
}

inline void StructBCSet::SetBCRegion(uint_t iBCRegion, StructBC *bcregion)
{
    (*bcRegions)[iBCRegion] = bcregion;
}

inline int StructBCSet::GetnBCRegion() const
{
    return static_cast<int>(bcRegions->size());
}