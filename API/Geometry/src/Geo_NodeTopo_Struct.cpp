#include "Geo_NodeTopo_Struct.h"
#include "GlobalDataBase.h"
#include "TypeDefine.h"
#include "Constants.h"
namespace PHSPACE
{

Geo_NodeTopo_Struct::Geo_NodeTopo_Struct()
{
    ni = 0;
    nj = 0;
    nk = 0;
    structx = NULL;
    structy = NULL;
    structz = NULL;
}

Geo_NodeTopo_Struct::~Geo_NodeTopo_Struct()
{
    DeallocateAll();
}

void Geo_NodeTopo_Struct::SetArrayLayout(RDouble *x, RDouble *y, RDouble *z)
{
    int nsimutask = -1;
    nsimutask = GlobalDataBase::GetIntParaFromDB("nsimutask");
    if (nsimutask == 0)
    {
        int isFVMOrFDM = PHSPACE::GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
        if (isFVMOrFDM == FD_METHOD)
        {
            //! struct hight order, xyz ghost more layers.
            SetArrayLayoutStructHighOrder(x, y, z);
            return;
        }

        int compressible = COMPRESSIBLE;
        compressible = GlobalDataBase::GetIntParaFromDB("compressible");
        if (compressible == INCOMPRESSIBLE)
        {
            Range I(1, ni);
            Range J(1, nj);
            Range K(1, nk);
            DeallocateAll();
            RDouble3D *strx = new RDouble3D(x, I, J, K, neverDeleteData, fortranArray);
            RDouble3D *stry = new RDouble3D(y, I, J, K, neverDeleteData, fortranArray);
            RDouble3D *strz = new RDouble3D(z, I, J, K, neverDeleteData, fortranArray);

            Range II(0, ni+1);
            Range JJ(0, nj+1);
            Range KK(0, nk+1);
            if (nk == 1) KK.setRange(1, 1);
            structx = new RDouble3D(II, JJ, KK, fortranArray);
            structy = new RDouble3D(II, JJ, KK, fortranArray);
            structz = new RDouble3D(II, JJ, KK, fortranArray);

            (*structx)(I, J, K) = (*strx)(I, J, K);
            (*structy)(I, J, K) = (*stry)(I, J, K);
            (*structz)(I, J, K) = (*strz)(I, J, K);
            return;
        }

        bool Rwithghost = true;
        
        int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
        //int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
        int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
        //int nMGLevel = PHSPACE::GlobalDataBase::GetIntParaFromDB("nMGLevel");
        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
        int ifLowSpeedPrecon = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");
        int wallFunctionType = GlobalDataBase::GetIntParaFromDB("wallFunctionType");

        if (tscheme != 4) Rwithghost = false;
        //if (systemGridType != 1) Rwithghost = false;
        if (isOverset) Rwithghost = false;
        //if (nMGLevel != 1) Rwithghost = false;
        if (compressible != COMPRESSIBLE) Rwithghost = false;
        if (isUnsteady) Rwithghost = false;
        if (nChemical != 0) Rwithghost = false;
        if (ifLowSpeedPrecon != 0) Rwithghost = false;
        if (wallFunctionType != WALLFUNCTION::NONE) Rwithghost = false;

        if (Rwithghost)
        {
            Range I(1, ni);
            Range J(1, nj);
            Range K(1, nk);
            DeallocateAll();
            RDouble3D *strx = new RDouble3D(x, I, J, K, neverDeleteData, fortranArray);
            RDouble3D *stry = new RDouble3D(y, I, J, K, neverDeleteData, fortranArray);
            RDouble3D *strz = new RDouble3D(z, I, J, K, neverDeleteData, fortranArray);
            
            Range II(0, ni+1);
            Range JJ(0, nj+1);
            Range KK(0, nk+1);
            if (nk == 1) KK.setRange(1, 1);
            structx = new RDouble3D(II, JJ, KK, fortranArray);
            structy = new RDouble3D(II, JJ, KK, fortranArray);
            structz = new RDouble3D(II, JJ, KK, fortranArray);
            
            (*structx)(I, J, K) = (*strx)(I, J, K);
            (*structy)(I, J, K) = (*stry)(I, J, K);
            (*structz)(I, J, K) = (*strz)(I, J, K);

            for (int k = 1; k <= nk; ++ k)
            {
                for (int j = 1; j <= nj; ++ j)
                {
                    (*structx)(0, j, k) = 2.0 * (*structx)(1, j, k) - (*structx)(2, j, k);
                    (*structy)(0, j, k) = 2.0 * (*structy)(1, j, k) - (*structy)(2, j, k);
                    (*structz)(0, j, k) = 2.0 * (*structz)(1, j, k) - (*structz)(2, j, k);
                    
                    (*structx)(ni+1, j, k) = 2.0 * (*structx)(ni, j, k) - (*structx)(ni-1, j, k);
                    (*structy)(ni+1, j, k) = 2.0 * (*structy)(ni, j, k) - (*structy)(ni-1, j, k);
                    (*structz)(ni+1, j, k) = 2.0 * (*structz)(ni, j, k) - (*structz)(ni-1, j, k);
                }
            }
            for (int k = 1; k <= nk; ++ k)
            {
                for (int i = 1; i <= ni; ++ i)
                {
                    (*structx)(i, 0, k) = 2.0 * (*structx)(i, 1, k) - (*structx)(i, 2, k);
                    (*structy)(i, 0, k) = 2.0 * (*structy)(i, 1, k) - (*structy)(i, 2, k);
                    (*structz)(i, 0, k) = 2.0 * (*structz)(i, 1, k) - (*structz)(i, 2, k);
                    
                    (*structx)(i, nj+1, k) = 2.0 * (*structx)(i, nj, k) - (*structx)(i, nj-1, k);
                    (*structy)(i, nj+1, k) = 2.0 * (*structy)(i, nj, k) - (*structy)(i, nj-1, k);
                    (*structz)(i, nj+1, k) = 2.0 * (*structz)(i, nj, k) - (*structz)(i, nj-1, k);
                }
            }
            if (nk != 1)
            {
                for (int j = 1; j <= nj; ++ j)
                {
                    for (int i = 1; i <= ni; ++ i)
                    {
                        (*structx)(i, j, 0) = 2.0 * (*structx)(i, j, 1) - (*structx)(i, j, 2);
                        (*structy)(i, j, 0) = 2.0 * (*structy)(i, j, 1) - (*structy)(i, j, 2);
                        (*structz)(i, j, 0) = 2.0 * (*structz)(i, j, 1) - (*structz)(i, j, 2);
                        
                        (*structx)(i, j, nk+1) = 2.0 * (*structx)(i, j, nk) - (*structx)(i, j, nk-1);
                        (*structy)(i, j, nk+1) = 2.0 * (*structy)(i, j, nk) - (*structy)(i, j, nk-1);
                        (*structz)(i, j, nk+1) = 2.0 * (*structz)(i, j, nk) - (*structz)(i, j, nk-1);
                    }
                }
            }
            delete strx;    strx = NULL;
            delete stry;    stry = NULL;
            delete strz;    strz = NULL;
            return;
        }
    }

    Range I(1, ni);
    Range J(1, nj);
    Range K(1, nk);
    DeallocateAll();
    structx = new RDouble3D(x, I, J, K, neverDeleteData, fortranArray);
    structy = new RDouble3D(y, I, J, K, neverDeleteData, fortranArray);
    structz = new RDouble3D(z, I, J, K, neverDeleteData, fortranArray);
}

void Geo_NodeTopo_Struct::SetArrayLayoutStructHighOrder(RDouble *x, RDouble *y, RDouble *z)
{
    string highordersolvername = GlobalDataBase::GetStrParaFromDB("str_highorder_solver");
    if (highordersolvername.substr(0, 4) == "WCNS")
    {
        Range I(1, ni);
        Range J(1, nj);
        Range K(1, nk);
        DeallocateAll();

        RDouble3D *strx = new RDouble3D(x, I, J, K, neverDeleteData, fortranArray);
        RDouble3D *stry = new RDouble3D(y, I, J, K, neverDeleteData, fortranArray);
        RDouble3D *strz = new RDouble3D(z, I, J, K, neverDeleteData, fortranArray);

        Range II(0, ni+1);
        Range JJ(0, nj+1);
        Range KK(0, nk+1);
        if (nk == 1) KK.setRange(1, 1);
        structx = new RDouble3D(II, JJ, KK, fortranArray);
        structy = new RDouble3D(II, JJ, KK, fortranArray);
        structz = new RDouble3D(II, JJ, KK, fortranArray);

        (*structx)(I, J, K) = (*strx)(I, J, K);
        (*structy)(I, J, K) = (*stry)(I, J, K);
        (*structz)(I, J, K) = (*strz)(I, J, K);

        delete strx;
        delete stry;
        delete strz;

        for (int k = 1; k <= nk; ++ k)
        {
            for (int j = 1; j <= nj; ++ j)
            {
                (*structx)(0, j, k) = 2.0 * (*structx)(1, j, k) - (*structx)(2, j, k);
                (*structy)(0, j, k) = 2.0 * (*structy)(1, j, k) - (*structy)(2, j, k);
                (*structz)(0, j, k) = 2.0 * (*structz)(1, j, k) - (*structz)(2, j, k);
                
                (*structx)(ni + 1, j, k) = 2.0 * (*structx)(ni, j, k) - (*structx)(ni - 1, j, k);
                (*structy)(ni + 1, j, k) = 2.0 * (*structy)(ni, j, k) - (*structy)(ni - 1, j, k);
                (*structz)(ni + 1, j, k) = 2.0 * (*structz)(ni, j, k) - (*structz)(ni - 1, j, k);
            }
        }

        for (int k = 1; k <= nk; ++ k)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                (*structx)(i, 0, k) = 2.0 * (*structx)(i, 1, k) - (*structx)(i, 2, k);
                (*structy)(i, 0, k) = 2.0 * (*structy)(i, 1, k) - (*structy)(i, 2, k);
                (*structz)(i, 0, k) = 2.0 * (*structz)(i, 1, k) - (*structz)(i, 2, k);
                
                (*structx)(i, nj + 1, k) = 2.0 * (*structx)(i, nj, k) - (*structx)(i, nj - 1, k);
                (*structy)(i, nj + 1, k) = 2.0 * (*structy)(i, nj, k) - (*structy)(i, nj - 1, k);
                (*structz)(i, nj + 1, k) = 2.0 * (*structz)(i, nj, k) - (*structz)(i, nj - 1, k);
            }
        }

        if (nk == 1)
        {
            int k = 1;
            (*structx)(0, 0, k) = 2.0 * (*structx)(1, 0, k) - (*structx)(2, 0, k);
            (*structy)(0, 0, k) = 2.0 * (*structy)(1, 0, k) - (*structy)(2, 0, k);
            (*structz)(0, 0, k) = 2.0 * (*structz)(1, 0, k) - (*structz)(2, 0, k);

            (*structx)(ni + 1, 0, k) = 2.0 * (*structx)(ni, 0, k) - (*structx)(ni - 1, 0, k);
            (*structy)(ni + 1, 0, k) = 2.0 * (*structy)(ni, 0, k) - (*structy)(ni - 1, 0, k);
            (*structz)(ni + 1, 0, k) = 2.0 * (*structz)(ni, 0, k) - (*structz)(ni - 1, 0, k);

            (*structx)(0, nj + 1, k) = 2.0 * (*structx)(1, nj + 1, k) - (*structx)(2, nj + 1, k);
            (*structy)(0, nj + 1, k) = 2.0 * (*structy)(1, nj + 1, k) - (*structy)(2, nj + 1, k);
            (*structz)(0, nj + 1, k) = 2.0 * (*structz)(1, nj + 1, k) - (*structz)(2, nj + 1, k);

            (*structx)(ni + 1, nj + 1, k) = 2.0 * (*structx)(ni, nj + 1, k) - (*structx)(ni - 1, nj + 1, k);
            (*structy)(ni + 1, nj + 1, k) = 2.0 * (*structy)(ni, nj + 1, k) - (*structy)(ni - 1, nj + 1, k);
            (*structz)(ni + 1, nj + 1, k) = 2.0 * (*structz)(ni, nj + 1, k) - (*structz)(ni - 1, nj + 1, k);

            return;
        }

        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                (*structx)(i, j, 0) = 2.0 * (*structx)(i, j, 1) - (*structx)(i, j, 2);
                (*structy)(i, j, 0) = 2.0 * (*structy)(i, j, 1) - (*structy)(i, j, 2);
                (*structz)(i, j, 0) = 2.0 * (*structz)(i, j, 1) - (*structz)(i, j, 2);
                
                (*structx)(i, j, nk + 1) = 2.0 * (*structx)(i, j, nk) - (*structx)(i, j, nk - 1);
                (*structy)(i, j, nk + 1) = 2.0 * (*structy)(i, j, nk) - (*structy)(i, j, nk - 1);
                (*structz)(i, j, nk + 1) = 2.0 * (*structz)(i, j, nk) - (*structz)(i, j, nk - 1);
            }
        }

        for (int j = 1; j <= nj; ++ j)
        {
            int i = 0;
            (*structx)(i, j, 0) = 2.0 * (*structx)(i, j, 1) - (*structx)(i, j, 2);
            (*structy)(i, j, 0) = 2.0 * (*structy)(i, j, 1) - (*structy)(i, j, 2);
            (*structz)(i, j, 0) = 2.0 * (*structz)(i, j, 1) - (*structz)(i, j, 2);

            (*structx)(i, j, nk + 1) = 2.0 * (*structx)(i, j, nk) - (*structx)(i, j, nk - 1);
            (*structy)(i, j, nk + 1) = 2.0 * (*structy)(i, j, nk) - (*structy)(i, j, nk - 1);
            (*structz)(i, j, nk + 1) = 2.0 * (*structz)(i, j, nk) - (*structz)(i, j, nk - 1);

            i = ni + 1;
            (*structx)(i, j, 0) = 2.0 * (*structx)(i, j, 1) - (*structx)(i, j, 2);
            (*structy)(i, j, 0) = 2.0 * (*structy)(i, j, 1) - (*structy)(i, j, 2);
            (*structz)(i, j, 0) = 2.0 * (*structz)(i, j, 1) - (*structz)(i, j, 2);

            (*structx)(i, j, nk + 1) = 2.0 * (*structx)(i, j, nk) - (*structx)(i, j, nk - 1);
            (*structy)(i, j, nk + 1) = 2.0 * (*structy)(i, j, nk) - (*structy)(i, j, nk - 1);
            (*structz)(i, j, nk + 1) = 2.0 * (*structz)(i, j, nk) - (*structz)(i, j, nk - 1);
        }

        for (int k = 1; k <= nk; ++ k)
        {
            int i = 0;
            (*structx)(i, 0, k) = 2.0 * (*structx)(i, 1, k) - (*structx)(i, 2, k);
            (*structy)(i, 0, k) = 2.0 * (*structy)(i, 1, k) - (*structy)(i, 2, k);
            (*structz)(i, 0, k) = 2.0 * (*structz)(i, 1, k) - (*structz)(i, 2, k);

            (*structx)(i, nj + 1, k) = 2.0 * (*structx)(i, nj, k) - (*structx)(i, nj - 1, k);
            (*structy)(i, nj + 1, k) = 2.0 * (*structy)(i, nj, k) - (*structy)(i, nj - 1, k);
            (*structz)(i, nj + 1, k) = 2.0 * (*structz)(i, nj, k) - (*structz)(i, nj - 1, k);

            i = ni + 1;
            (*structx)(i, 0, k) = 2.0 * (*structx)(i, 1, k) - (*structx)(i, 2, k);
            (*structy)(i, 0, k) = 2.0 * (*structy)(i, 1, k) - (*structy)(i, 2, k);
            (*structz)(i, 0, k) = 2.0 * (*structz)(i, 1, k) - (*structz)(i, 2, k);

            (*structx)(i, nj + 1, k) = 2.0 * (*structx)(i, nj, k) - (*structx)(i, nj - 1, k);
            (*structy)(i, nj + 1, k) = 2.0 * (*structy)(i, nj, k) - (*structy)(i, nj - 1, k);
            (*structz)(i, nj + 1, k) = 2.0 * (*structz)(i, nj, k) - (*structz)(i, nj - 1, k);
        }

        for (int i = 1; i <= ni; ++ i)
        {
            int j = 0;
            (*structx)(i, j, 0) = 2.0 * (*structx)(i, j, 1) - (*structx)(i, j, 2);
            (*structy)(i, j, 0) = 2.0 * (*structy)(i, j, 1) - (*structy)(i, j, 2);
            (*structz)(i, j, 0) = 2.0 * (*structz)(i, j, 1) - (*structz)(i, j, 2);

            (*structx)(i, j, nk + 1) = 2.0 * (*structx)(i, j, nk) - (*structx)(i, j, nk - 1);
            (*structy)(i, j, nk + 1) = 2.0 * (*structy)(i, j, nk) - (*structy)(i, j, nk - 1);
            (*structz)(i, j, nk + 1) = 2.0 * (*structz)(i, j, nk) - (*structz)(i, j, nk - 1);

            j = nj + 1;
            (*structx)(i, j, 0) = 2.0 * (*structx)(i, j, 1) - (*structx)(i, j, 2);
            (*structy)(i, j, 0) = 2.0 * (*structy)(i, j, 1) - (*structy)(i, j, 2);
            (*structz)(i, j, 0) = 2.0 * (*structz)(i, j, 1) - (*structz)(i, j, 2);

            (*structx)(i, j, nk + 1) = 2.0 * (*structx)(i, j, nk) - (*structx)(i, j, nk - 1);
            (*structy)(i, j, nk + 1) = 2.0 * (*structy)(i, j, nk) - (*structy)(i, j, nk - 1);
            (*structz)(i, j, nk + 1) = 2.0 * (*structz)(i, j, nk) - (*structz)(i, j, nk - 1);
        }

        for (int k = 1; k <= nk; ++ k)
        {
            int j = 0;
            (*structx)(0, j, k) = 2.0 * (*structx)(1, j, k) - (*structx)(2, j, k);
            (*structy)(0, j, k) = 2.0 * (*structy)(1, j, k) - (*structy)(2, j, k);
            (*structz)(0, j, k) = 2.0 * (*structz)(1, j, k) - (*structz)(2, j, k);

            (*structx)(ni + 1, j, k) = 2.0 * (*structx)(ni, j, k) - (*structx)(ni - 1, j, k);
            (*structy)(ni + 1, j, k) = 2.0 * (*structy)(ni, j, k) - (*structy)(ni - 1, j, k);
            (*structz)(ni + 1, j, k) = 2.0 * (*structz)(ni, j, k) - (*structz)(ni - 1, j, k);

            j = nj + 1;
            (*structx)(0, j, k) = 2.0 * (*structx)(1, j, k) - (*structx)(2, j, k);
            (*structy)(0, j, k) = 2.0 * (*structy)(1, j, k) - (*structy)(2, j, k);
            (*structz)(0, j, k) = 2.0 * (*structz)(1, j, k) - (*structz)(2, j, k);

            (*structx)(ni + 1, j, k) = 2.0 * (*structx)(ni, j, k) - (*structx)(ni - 1, j, k);
            (*structy)(ni + 1, j, k) = 2.0 * (*structy)(ni, j, k) - (*structy)(ni - 1, j, k);
            (*structz)(ni + 1, j, k) = 2.0 * (*structz)(ni, j, k) - (*structz)(ni - 1, j, k);
        }

        for (int i = 1; i <= ni; ++ i)
        {
            int k = 0;
            (*structx)(i, 0, k) = 2.0 * (*structx)(i, 1, k) - (*structx)(i, 2, k);
            (*structy)(i, 0, k) = 2.0 * (*structy)(i, 1, k) - (*structy)(i, 2, k);
            (*structz)(i, 0, k) = 2.0 * (*structz)(i, 1, k) - (*structz)(i, 2, k);

            (*structx)(i, nj + 1, k) = 2.0 * (*structx)(i, nj, k) - (*structx)(i, nj - 1, k);
            (*structy)(i, nj + 1, k) = 2.0 * (*structy)(i, nj, k) - (*structy)(i, nj - 1, k);
            (*structz)(i, nj + 1, k) = 2.0 * (*structz)(i, nj, k) - (*structz)(i, nj - 1, k);

            k = nk + 1;
            (*structx)(i, 0, k) = 2.0 * (*structx)(i, 1, k) - (*structx)(i, 2, k);
            (*structy)(i, 0, k) = 2.0 * (*structy)(i, 1, k) - (*structy)(i, 2, k);
            (*structz)(i, 0, k) = 2.0 * (*structz)(i, 1, k) - (*structz)(i, 2, k);
            
            (*structx)(i, nj + 1, k) = 2.0 * (*structx)(i, nj, k) - (*structx)(i, nj - 1, k);
            (*structy)(i, nj + 1, k) = 2.0 * (*structy)(i, nj, k) - (*structy)(i, nj - 1, k);
            (*structz)(i, nj + 1, k) = 2.0 * (*structz)(i, nj, k) - (*structz)(i, nj - 1, k);
        }

        for (int j = 1; j <= nj; ++ j)
        {
            int k = 0;
            (*structx)(0, j, k) = 2.0 * (*structx)(1, j, k) - (*structx)(2, j, k);
            (*structy)(0, j, k) = 2.0 * (*structy)(1, j, k) - (*structy)(2, j, k);
            (*structz)(0, j, k) = 2.0 * (*structz)(1, j, k) - (*structz)(2, j, k);

            (*structx)(ni + 1, j, k) = 2.0 * (*structx)(ni, j, k) - (*structx)(ni - 1, j, k);
            (*structy)(ni + 1, j, k) = 2.0 * (*structy)(ni, j, k) - (*structy)(ni - 1, j, k);
            (*structz)(ni + 1, j, k) = 2.0 * (*structz)(ni, j, k) - (*structz)(ni - 1, j, k);

            k = nk + 1;
            (*structx)(0, j, k) = 2.0 * (*structx)(1, j, k) - (*structx)(2, j, k);
            (*structy)(0, j, k) = 2.0 * (*structy)(1, j, k) - (*structy)(2, j, k);
            (*structz)(0, j, k) = 2.0 * (*structz)(1, j, k) - (*structz)(2, j, k);

            (*structx)(ni + 1, j, k) = 2.0 * (*structx)(ni, j, k) - (*structx)(ni - 1, j, k);
            (*structy)(ni + 1, j, k) = 2.0 * (*structy)(ni, j, k) - (*structy)(ni - 1, j, k);
            (*structz)(ni + 1, j, k) = 2.0 * (*structz)(ni, j, k) - (*structz)(ni - 1, j, k);
        }
    }
    else if (highordersolvername.substr(0, 4) == "HDCS")
    {
        Range I(1, ni);
        Range J(1, nj);
        Range K(1, nk);
        DeallocateAll();

        RDouble3D *strx = new RDouble3D(x, I, J, K, neverDeleteData, fortranArray);
        RDouble3D *stry = new RDouble3D(y, I, J, K, neverDeleteData, fortranArray);
        RDouble3D *strz = new RDouble3D(z, I, J, K, neverDeleteData, fortranArray);

        Range II(-1, ni+2);
        Range JJ(-1, nj+2);
        Range KK(-1, nk+2);
        if (nk == 1) KK.setRange(1, 1);
        structx = new RDouble3D(II, JJ, KK, fortranArray);
        structy = new RDouble3D(II, JJ, KK, fortranArray);
        structz = new RDouble3D(II, JJ, KK, fortranArray);

        (*structx)(I, J, K) = (*strx)(I, J, K);
        (*structy)(I, J, K) = (*stry)(I, J, K);
        (*structz)(I, J, K) = (*strz)(I, J, K);

        delete strx;
        delete stry;
        delete strz;

        for (int k = 1; k <= nk; ++ k)
        {
            for (int j = 1; j <= nj; ++ j)
            {
                (*structx)(0, j, k) = 2.0 * (*structx)(1, j, k) - (*structx)(2, j, k);
                (*structy)(0, j, k) = 2.0 * (*structy)(1, j, k) - (*structy)(2, j, k);
                (*structz)(0, j, k) = 2.0 * (*structz)(1, j, k) - (*structz)(2, j, k);

                (*structx)(-1, j, k) = 2.0 * (*structx)(0, j, k) - (*structx)(1, j, k);
                (*structy)(-1, j, k) = 2.0 * (*structy)(0, j, k) - (*structy)(1, j, k);
                (*structz)(-1, j, k) = 2.0 * (*structz)(0, j, k) - (*structz)(1, j, k);
                
                (*structx)(ni + 1, j, k) = 2.0 * (*structx)(ni, j, k) - (*structx)(ni - 1, j, k);
                (*structy)(ni + 1, j, k) = 2.0 * (*structy)(ni, j, k) - (*structy)(ni - 1, j, k);
                (*structz)(ni + 1, j, k) = 2.0 * (*structz)(ni, j, k) - (*structz)(ni - 1, j, k);

                (*structx)(ni + 2, j, k) = 2.0 * (*structx)(ni + 1, j, k) - (*structx)(ni, j, k);
                (*structy)(ni + 2, j, k) = 2.0 * (*structy)(ni + 1, j, k) - (*structy)(ni, j, k);
                (*structz)(ni + 2, j, k) = 2.0 * (*structz)(ni + 1, j, k) - (*structz)(ni, j, k);
            }
        }

        for (int k = 1; k <= nk; ++ k)
        {
            for (int i = -1; i <= ni+2; ++ i)
            {
                (*structx)(i, 0, k) = 2.0 * (*structx)(i, 1, k) - (*structx)(i, 2, k);
                (*structy)(i, 0, k) = 2.0 * (*structy)(i, 1, k) - (*structy)(i, 2, k);
                (*structz)(i, 0, k) = 2.0 * (*structz)(i, 1, k) - (*structz)(i, 2, k);

                (*structx)(i, -1, k) = 2.0 * (*structx)(i, 0, k) - (*structx)(i, 1, k);
                (*structy)(i, -1, k) = 2.0 * (*structy)(i, 0, k) - (*structy)(i, 1, k);
                (*structz)(i, -1, k) = 2.0 * (*structz)(i, 0, k) - (*structz)(i, 1, k);
                
                (*structx)(i, nj + 1, k) = 2.0 * (*structx)(i, nj, k) - (*structx)(i, nj - 1, k);
                (*structy)(i, nj + 1, k) = 2.0 * (*structy)(i, nj, k) - (*structy)(i, nj - 1, k);
                (*structz)(i, nj + 1, k) = 2.0 * (*structz)(i, nj, k) - (*structz)(i, nj - 1, k);

                (*structx)(i, nj + 2, k) = 2.0 * (*structx)(i, nj + 1, k) - (*structx)(i, nj, k);
                (*structy)(i, nj + 2, k) = 2.0 * (*structy)(i, nj + 1, k) - (*structy)(i, nj, k);
                (*structz)(i, nj + 2, k) = 2.0 * (*structz)(i, nj + 1, k) - (*structz)(i, nj, k);
            }
        }

        if (nk == 1)
        {
            return;
        }

        for (int j = -1; j <= nj+2; ++ j)
        {
            for (int i = -1; i <= ni+2; ++ i)
            {
                (*structx)(i, j, 0) = 2.0 * (*structx)(i, j, 1) - (*structx)(i, j, 2);
                (*structy)(i, j, 0) = 2.0 * (*structy)(i, j, 1) - (*structy)(i, j, 2);
                (*structz)(i, j, 0) = 2.0 * (*structz)(i, j, 1) - (*structz)(i, j, 2);

                (*structx)(i, j, -1) = 2.0 * (*structx)(i, j, 0) - (*structx)(i, j, 1);
                (*structy)(i, j, -1) = 2.0 * (*structy)(i, j, 0) - (*structy)(i, j, 1);
                (*structz)(i, j, -1) = 2.0 * (*structz)(i, j, 0) - (*structz)(i, j, 1);

                (*structx)(i, j, nk + 1) = 2.0 * (*structx)(i, j, nk) - (*structx)(i, j, nk - 1);
                (*structy)(i, j, nk + 1) = 2.0 * (*structy)(i, j, nk) - (*structy)(i, j, nk - 1);
                (*structz)(i, j, nk + 1) = 2.0 * (*structz)(i, j, nk) - (*structz)(i, j, nk - 1);

                (*structx)(i, j, nk + 2) = 2.0 * (*structx)(i, j, nk + 1) - (*structx)(i, j, nk);
                (*structy)(i, j, nk + 2) = 2.0 * (*structy)(i, j, nk + 1) - (*structy)(i, j, nk);
                (*structz)(i, j, nk + 2) = 2.0 * (*structz)(i, j, nk + 1) - (*structz)(i, j, nk);
            }
        }
    }
}

}