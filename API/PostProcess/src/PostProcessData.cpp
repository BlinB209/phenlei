#include "PostProcessData.h"

using namespace std;
namespace PHSPACE
{

PostProcessData::PostProcessData()
{
    partBCID = 0;

    cfx = 0;
    cfy = 0;
    cfz = 0;
    cpx = 0;
    cpy = 0;
    cpz = 0;
    cmx = 0;
    cmy = 0;
    cmz = 0;

    clTot    = 0;
    cdTot    = 0;
    cd_pr    = 0;
    cd_sf    = 0;
    cd_cl2pa = 0;
    xcp      = 0;
    side     = 0;

    hingeMoment = 0;
    numberofSolidWallPart = 0;
    numberofTotalPart = 0;
}

PostProcessData::~PostProcessData()
{
    delete [] partBCID;

    delete [] cfx;
    delete [] cfy;
    delete [] cfz;
    delete [] cpx;
    delete [] cpy;
    delete [] cpz;
    delete [] cmx;
    delete [] cmy;
    delete [] cmz;

    delete [] clTot;
    delete [] cdTot;
    delete [] cd_pr;
    delete [] cd_sf;
    delete [] cd_cl2pa;
    delete [] xcp;
    delete [] side;

    delete [] hingeMoment;
}

void PostProcessData::SetNWallPart(int numberofSolidWallPartIn)
{
    this->numberofSolidWallPart = numberofSolidWallPartIn;
}

void PostProcessData::SetNPart(int numberofTotalPartIn)
{
    this->numberofTotalPart = numberofTotalPartIn;
}

void PostProcessData::SetPartBCID(int *partBCIDIn)
{
    if (this->partBCID)
    {
        delete this->partBCID;
    }
    this->partBCID = partBCIDIn;
}

void PostProcessData::SetCfx(RDouble *cfxIn)
{
    if (this->cfx)
    {
        delete this->cfx;
    }
    this->cfx = cfxIn;
}

void PostProcessData::SetCfy(RDouble *cfyIn)
{
    if (this->cfy)
    {
        delete this->cfy;
    }
    this->cfy = cfyIn;
}

void PostProcessData::SetCfz(RDouble *cfzIn)
{
    if (this->cfz)
    {
        delete this->cfz;
    }
    this->cfz = cfzIn;
}

void PostProcessData::SetCpx(RDouble *cpxIn)
{
    if (this->cpx)
    {
        delete this->cpx;
    }
    this->cpx = cpxIn;
}

void PostProcessData::SetCpy(RDouble *cpyIn)
{
    if (this->cpy)
    {
        delete this->cpy;
    }
    this->cpy = cpyIn;
}

void PostProcessData::SetCpz(RDouble *cpzIn)
{
    if (this->cpz)
    {
        delete this->cpz;
    }
    this->cpz = cpzIn;
}

void PostProcessData::SetCmx(RDouble *cmxIn)
{
    if (this->cmx)
    {
        delete this->cmx;
    }
    this->cmx = cmxIn;
}

void PostProcessData::SetCmy(RDouble *cmyIn)
{
    if (this->cmy)
    {
        delete this->cmy;
    }
    this->cmy = cmyIn;
}

void PostProcessData::SetCmz(RDouble *cmzIn)
{
    if (this->cmz)
    {
        delete this->cmz;
    }
    this->cmz = cmzIn;
}

void PostProcessData::SetClTot(RDouble *clTotIn)
{
    if (this->clTot)
    {
        delete this->clTot;
    }
    this->clTot = clTotIn;
}

void PostProcessData::SetCdTot(RDouble *cdTotIn)
{
    if (this->cdTot)
    {
        delete this->cdTot;
    }
    this->cdTot = cdTotIn;
}

void PostProcessData::SetCdPr(RDouble *cd_prIn)
{
    if (this->cd_pr)
    {
        delete this->cd_pr;
    }
    this->cd_pr = cd_prIn;
}

void PostProcessData::SetCdSf(RDouble *cd_sfIn)
{
    if (this->cd_sf)
    {
        delete this->cd_sf;
    }
    this->cd_sf = cd_sfIn;
}

void PostProcessData::SetCdCl2Pa(RDouble *cd_cl2paIn)
{
    if (this->cd_cl2pa)
    {
        delete this->cd_cl2pa;
    }
    this->cd_cl2pa = cd_cl2paIn;
}

void PostProcessData::SetXcp(RDouble *xcpIn)
{
    if (this->xcp)
    {
        delete this->xcp;
    }
    this->xcp = xcpIn;
}

void PostProcessData::SetSide(RDouble *sideIn)
{
    if (this->side)
    {
        delete this->side;
    }
    this->side = sideIn;
}

void PostProcessData::SetHingeMoment(RDouble *hingeMomentIn)
{
    if (this->hingeMoment)
    {
        delete this->hingeMoment;
    }
    this->hingeMoment = hingeMomentIn;
}

int PostProcessData::GetNWallPart()
{
    return this->numberofSolidWallPart;
}

int PostProcessData::GetNPart()
{
    return this->numberofTotalPart;
}

int * PostProcessData::GetPartBCID()
{
    return this->partBCID;
};

RDouble * PostProcessData::GetCfx()
{
    return this->cfx; 
};

RDouble * PostProcessData::GetCfy() 
{
    return this->cfy; 
};

RDouble * PostProcessData::GetCfz() 
{
    return this->cfz; 
};

RDouble * PostProcessData::GetCpx() 
{
    return this->cpx; 
};

RDouble * PostProcessData::GetCpy() 
{
    return this->cpy; 
};

RDouble * PostProcessData::GetCpz() 
{
    return this->cpz; 
};

RDouble * PostProcessData::GetCmx() 
{
    return this->cmx; 
};

RDouble * PostProcessData::GetCmy() 
{
    return this->cmy; 
};

RDouble * PostProcessData::GetCmz() 
{
    return this->cmz; 
};

RDouble * PostProcessData::GetClTot() 
{
    return this->clTot; 
};

RDouble * PostProcessData::GetCdTot() 
{
    return this->cdTot; 
};

RDouble * PostProcessData::GetCdPr() 
{
    return this->cd_pr; 
};

RDouble * PostProcessData::GetCdSf() 
{
    return this->cd_sf; 
};

RDouble * PostProcessData::GetCdCl2Pa() 
{
    return this->cd_cl2pa; 
};

RDouble * PostProcessData::GetXcp() 
{
    return this->xcp; 
};

RDouble * PostProcessData::GetSide() 
{
    return this->side; 
};

RDouble * PostProcessData::GetHingeMoment() 
{
    return this->hingeMoment; 
};

}