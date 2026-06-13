#include "Geo_StructGrid.h"
#include "Geo_UnstructGrid.h"
#include "Constants.h"
#include "Geo_MultiGridInfo_Struct.h"
#include "Geo_SimpleBC.h"
#include "Geo_StructBC.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "Geo_NodeTopo_Struct.h"
#include "Geo_FaceMetrics_Struct.h"
#include "Geo_CellMetrics_Struct.h"
#include "Geo_DynamicGridMetrics_Struct.h"
#include "Geo_OversetGridTopo_Struct.h"
#include "LinkStruct.h"
#include "PHMpi.h"

namespace PHSPACE
{
void StructGrid::AllocateMetricsStructHighOrder(ActionKey *actkey)
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();
    
    Range I, J, K;
    GetRange(ni, nj, nk, -4, 3, I, J, K);
    Range D(1, 3);
    
    faceMetrics->AllocateMetrics(I, J, K, D);
    cellMetrics->AllocateMetrics(I, J, K);
    
    RDouble4D *vgn = dynamicGridMetrics->GetFaceNormalVelocity();
    if (vgn == 0) vgn = new RDouble4D(I, J, K, Range(1, 3), fortranArray);
    
    AllocateMetricsALE(actkey);
    
    (*vgn) = zero;
    dynamicGridMetrics->SetFaceNormalVelocity(vgn);
}

void StructGrid::ComputeMetricsStructHighOrder2D(ActionKey *actkey)
{
    string highordersolvername = GlobalDataBase::GetStrParaFromDB("str_highorder_solver");
    if (highordersolvername.substr(0, 4) == "HDCS")
    {
        ComputeMetricsStructHighOrder2D_HDCS();
        return;
    }
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    int ist =  1;
    int ied = ni;
    int jst =  1;
    int jed = nj;

    int ni2 = 2 * ni - 1;
    int nj2 = 2 * nj - 1;
    int nk2 = 2 * nk - 1;

    int i2st =   1;
    int i2ed = ni2;
    int j2st =   1;
    int j2ed = nj2;

    RDouble3D &x2 = *new RDouble3D(Range(1, ni2), Range(1, nj2), Range(1, nk2), fortranArray);
    RDouble3D &y2 = *new RDouble3D(Range(1, ni2), Range(1, nj2), Range(1, nk2), fortranArray);
    RDouble3D &z2 = *new RDouble3D(Range(1, ni2), Range(1, nj2), Range(1, nk2), fortranArray);
    
    RDouble3D &xx = *this->GetStructX();
    RDouble3D &yy = *this->GetStructY();
    RDouble3D &zz = *this->GetStructZ();

    //! grid refined block by block

    int i, j, k, i2, j2, k2;
    
    k  = 1;
    k2 = 1;

    //! grid refined : copy primitive grid

    for (j = jst; j <= jed; ++ j)
    {
        j2 = 2 * j - 1;
        for (i = ist; i <= ied; ++ i)
        {
            i2 = 2 * i - 1;
            x2(i2, j2, k2) = xx(i, j, k);
            y2(i2, j2, k2) = yy(i, j, k);
            z2(i2, j2, k2) = zz(i, j, k);
        }
    }

    //! grid refined : I direction

    for (j = jst; j <= jed; ++ j)
    {
        j2 = 2 * j - 1;
        for (i = ist; i <= ied - 1; ++ i)
        {
            i2 = 2 * i;
            x2(i2, j2, k2) = half * (x2(i2 - 1, j2, k2) + x2(i2 + 1, j2, k2));
            y2(i2, j2, k2) = half * (y2(i2 - 1, j2, k2) + y2(i2 + 1, j2, k2));
            z2(i2, j2, k2) = half * (z2(i2 - 1, j2, k2) + z2(i2 + 1, j2, k2));
        }
    }

    //! grid refined : J direction

    for (j = jst; j <= jed - 1; ++ j)
    {
        j2 = 2 * j;
        for (i2 = i2st; i2 <= i2ed; ++ i2)
        {
            x2(i2, j2, k2) = half * (x2(i2, j2 - 1, k2) + x2(i2, j2 + 1, k2));
            y2(i2, j2, k2) = half * (y2(i2, j2 - 1, k2) + y2(i2, j2 + 1, k2));
            z2(i2, j2, k2) = half * (z2(i2, j2 - 1, k2) + z2(i2, j2 + 1, k2));
        }
    }

    //! calculate mesh transform metrics

    RDouble2D * lineData1,* lineData2;

    //! calculate mesh transform metrics : first grid delta

    RDouble4D &gridDeltaFirst = * new RDouble4D(Range(1, 4), Range(1, ni2), Range(1, nj2), Range(1, nk2), fortranArray);

    //! calculate mesh transform metrics : first grid delta (I direction)

    lineData1 = new RDouble2D(Range(1, 2), Range(1, ni2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 2), Range(1, ni2), fortranArray);

    for (j2 = j2st; j2 <= j2ed; ++ j2)
    {
            
        for (i2 = i2st; i2 <= i2ed; ++ i2)
        {
            (*lineData1)(1, i2) = x2(i2, j2, k2);
            (*lineData1)(2, i2) = y2(i2, j2, k2);
        }
            
        GridDelta(Range(1, 2), Range(1, ni2), lineData1, lineData2);

        for (i2 = i2st; i2 <= i2ed; ++ i2)
        {
            gridDeltaFirst(1, i2, j2, k2) = (*lineData2)(1, i2);
            gridDeltaFirst(2, i2, j2, k2) = (*lineData2)(2, i2);
        }

    }

    delete lineData1;
    delete lineData2;

    //! calculate mesh transform metrics : first grid delta (J direction)

    lineData1 = new RDouble2D(Range(1, 2), Range(1, nj2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 2), Range(1, nj2), fortranArray);

    for (i2 = i2st; i2 <= i2ed; ++ i2)
    {
            
        for (j2 = j2st; j2 <= j2ed; ++ j2)
        {
            (*lineData1)(1, j2) = x2(i2, j2, k2);
            (*lineData1)(2, j2) = y2(i2, j2, k2);
        }
            
        GridDelta(Range(1, 2), Range(1, nj2), lineData1, lineData2);

        for (j2 = j2st; j2 <= j2ed; ++ j2)
        {
            gridDeltaFirst(3, i2, j2, k2) = (*lineData2)(1, j2);
            gridDeltaFirst(4, i2, j2, k2) = (*lineData2)(2, j2);
        }

    }

    delete lineData1;
    delete lineData2;

    //! calculate mesh transform metrics by grid delta

    RDouble4D &xfn  = *(this->GetFaceNormalX());
    RDouble4D &yfn  = *(this->GetFaceNormalY());
    RDouble4D &zfn  = *(this->GetFaceNormalZ());
    RDouble4D &area = *(this->GetFaceArea());
    RDouble4D &xfv  = *(this->GetFaceVectorX());
    RDouble4D &yfv  = *(this->GetFaceVectorY());
    RDouble4D &zfv  = *(this->GetFaceVectorZ());

    xfn = 0.0;
    yfn = 0.0;
    zfn = 0.0;
    area= 0.0;
    xfv = 0.0;
    yfv = 0.0;
    zfv = 0.0;

    RDouble vkc[2],vet[2],vkcn[3],vetn[3];

    //! calculate mesh transform metrics by grid delta (I direction)

    for (j = 1; j <= nj - 1; ++ j)
    {
        j2 = 2 * j;

        for (i = 1; i <= ni; ++ i)
        {

            i2 = 2 * i - 1;

            vkc[0] =   gridDeltaFirst(4, i2, j2, k2);
            vkc[1] = - gridDeltaFirst(3, i2, j2, k2);

            xfv(i, j, k, IDX::IKES) = vkc[0];
            yfv(i, j, k, IDX::IKES) = vkc[1];
            zfv(i, j, k, IDX::IKES) = zero  ;
            
            NormalizeVector(vkc,vkcn, 2);

            area(i, j, k, IDX::IKES) = vkcn[2];
            xfn (i, j, k, IDX::IKES) = vkcn[0];
            yfn (i, j, k, IDX::IKES) = vkcn[1];
            zfn (i, j, k, IDX::IKES) = zero   ;

        }
    }

    //! calculate mesh transform metrics by grid delta (J direction)

    for (i = 1; i <= ni - 1; ++ i)
    {
        i2 = 2 * i;

        for (j = 1; j <= nj; ++ j)
        {

            j2 = 2 * j - 1;

            vet[0] = - gridDeltaFirst(2, i2, j2, k2);
            vet[1] =   gridDeltaFirst(1, i2, j2, k2);

            xfv(i, j, k, IDX::IETA) = vet[0];
            yfv(i, j, k, IDX::IETA) = vet[1];
            zfv(i, j, k, IDX::IETA) = zero  ;

            NormalizeVector(vet,vetn, 2);

            area(i, j, k, IDX::IETA) = vetn[2];
            xfn (i, j, k, IDX::IETA) = vetn[0];
            yfn (i, j, k, IDX::IETA) = vetn[1];
            zfn (i, j, k, IDX::IETA) = zero   ;

        }
    }

    //! calculate mesh transform Jacobian : second grid delta

    RDouble4D &gridDeltaSecond = * new RDouble4D(Range(1, 2),Range(1, ni2),Range(1, nj2),Range(1, nk2), fortranArray);

    //! calculate mesh transform Jacobian : second grid delta (I direction)

    lineData1 = new RDouble2D(Range(1, 1), Range(1, ni2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 1), Range(1, ni2), fortranArray);

    for (j2 = j2st; j2 <= j2ed; ++ j2)
    {
            
        for (i2 = i2st; i2 <= i2ed; ++ i2)
        {
            (*lineData1)(1, i2) =   x2(i2, j2, k2) * gridDeltaFirst(4, i2, j2, k2)
                                 - y2(i2, j2, k2) * gridDeltaFirst(3, i2, j2, k2);
        }
            
        GridDelta(Range(1, 1), Range(1, ni2), lineData1, lineData2);

        for (i2 = i2st; i2 <= i2ed; ++ i2)
        {
            gridDeltaSecond(1, i2, j2, k2) = (*lineData2)(1, i2);
        }

    }

    delete lineData1;
    delete lineData2;

    //! calculate mesh transform Jacobian : second grid delta (J direction)

    lineData1 = new RDouble2D(Range(1, 1), Range(1, nj2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 1), Range(1, nj2), fortranArray);

    for (i2 = i2st; i2 <= i2ed; ++ i2)
    {
            
        for (j2 = j2st; j2 <= j2ed; ++ j2)
        {
            (*lineData1)(1, j2) = - x2(i2, j2, k2) * gridDeltaFirst(2, i2, j2, k2)
                                 + y2(i2, j2, k2) * gridDeltaFirst(1, i2, j2, k2);
        }
            
        GridDelta(Range(1, 1), Range(1, nj2), lineData1, lineData2);

        for (j2 = j2st; j2 <= j2ed; ++ j2)
        {
            gridDeltaSecond(2, i2, j2, k2) = (*lineData2)(1, j2);
        }

    }

    delete lineData1;
    delete lineData2;

    delete &x2;
    delete &y2;
    delete &z2;

    delete &gridDeltaFirst;

    //! calculate mesh transform Jacobian by second grid delta

    RDouble3D &vol  = *(this->GetCellVolume());
    int ZoneID = this->GetZoneID();
    int imin, jmin, kmin, imax, jmax, kmax;
    RDouble vv,volmin,volmax;

    volmin =   LARGE;
    volmax = - LARGE;

    std::ostringstream oss;

    for (j = 1; j <= nj - 1; ++ j)
    {
        j2 = 2 * j;
        for (i = 1; i <= ni - 1; ++ i)
        {
            i2 = 2 * i;
            vv = half * (gridDeltaSecond(1, i2, j2, k2) + gridDeltaSecond(2, i2, j2, k2));
            vol(i, j, k) = vv;

            if (vv <= zero)
            {
                oss << " vol < 0 " << " ZoneID = " << ZoneID << " i = " << i << " j = " << j << " vol = " << vv << "\n";
            }           

            if (vv < volmin)
            {
                volmin = vv;
                imin = i;
                jmin = j;
                kmin = k;
            }

            if (vv > volmax)
            {
                volmax = vv;
                imax = i;
                jmax = j;
                kmax = k;
            }

        }
    }

    delete &gridDeltaSecond;
    GhostCell3D(vol, ni, nj, 1);
    GhostMetrics2D();
    ComputeCellCenter();
    StreamToActionKey(actkey, oss);
}

void StructGrid::ComputeMetricsStructHighOrder3D(ActionKey *actkey)
{
    string highordersolvername = GlobalDataBase::GetStrParaFromDB("str_highorder_solver");
    if (highordersolvername.substr(0, 4) == "HDCS")
    {
        ComputeMetricsStructHighOrder3D_HDCS();
        return;
    }
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    int ist =  1;
    int ied = ni;
    int jst =  1;
    int jed = nj;
    int kst =  1;
    int ked = nk;

    int ni2 = 2 * ni - 1;
    int nj2 = 2 * nj - 1;
    int nk2 = 2 * nk - 1;

    int i2st =   1;
    int i2ed = ni2;
    int j2st =   1;
    int j2ed = nj2;
    int k2st =   1;
    int k2ed = nk2;
    
    RDouble3D &x2 = * new RDouble3D(Range(1, ni2), Range(1, nj2), Range(1, nk2), fortranArray);
    RDouble3D &y2 = * new RDouble3D(Range(1, ni2), Range(1, nj2), Range(1, nk2), fortranArray);
    RDouble3D &z2 = * new RDouble3D(Range(1, ni2), Range(1, nj2), Range(1, nk2), fortranArray);
    
    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    //! grid refined block by block

    int i, j, k, i2, j2, k2;

    //! grid refined : copy primitive grid

    for (k = kst; k <= ked; ++ k)
    {
        k2 = 2 * k - 1;
        for (j = jst; j <= jed; ++ j)
        {
            j2 = 2 * j - 1;
            for (i = ist; i <= ied; ++ i)
            {
                i2 = 2 * i - 1;
                x2(i2, j2, k2) = xx(i, j, k);
                y2(i2, j2, k2) = yy(i, j, k);
                z2(i2, j2, k2) = zz(i, j, k);
            }
        }
    }

    //! grid refined : I direction

    for (k = kst; k <= ked; ++ k)
    {
        k2 = 2 * k - 1;
        for (j = jst; j <= jed; ++ j)
        {
            j2 = 2 * j - 1;
            for (i = ist; i <= ied - 1; ++ i)
            {
                i2 = 2 * i;
                x2(i2, j2, k2) = half * (x2(i2 - 1, j2, k2) + x2(i2 + 1, j2, k2));
                y2(i2, j2, k2) = half * (y2(i2 - 1, j2, k2) + y2(i2 + 1, j2, k2));
                z2(i2, j2, k2) = half * (z2(i2 - 1, j2, k2) + z2(i2 + 1, j2, k2));
            }
        }
    }

    //! grid refined : J direction

    for (k = kst; k <= ked; ++ k)
    {
        k2 = 2 * k - 1;
        for (j = jst; j <= jed - 1; ++ j)
        {
            j2 = 2 * j;
            for (i2 = i2st; i2 <= i2ed; ++ i2)
            {
                x2(i2, j2, k2) = half * (x2(i2, j2 - 1, k2) + x2(i2, j2 + 1, k2));
                y2(i2, j2, k2) = half * (y2(i2, j2 - 1, k2) + y2(i2, j2 + 1, k2));
                z2(i2, j2, k2) = half * (z2(i2, j2 - 1, k2) + z2(i2, j2 + 1, k2));
            }
        }
    }

    //! grid refined : K direction

    for (k = kst; k <= ked - 1; ++ k)
    {
        k2 = 2 * k;
        for (j2 = j2st; j2 <= j2ed; ++ j2)
        {
            for (i2 = i2st; i2 <= i2ed; ++ i2)
            {
                x2(i2, j2, k2) = half * (x2(i2, j2, k2-1) + x2(i2, j2, k2+1));
                y2(i2, j2, k2) = half * (y2(i2, j2, k2-1) + y2(i2, j2, k2+1));
                z2(i2, j2, k2) = half * (z2(i2, j2, k2-1) + z2(i2, j2, k2+1));
            }
        }
    }

    //! calculate mesh transform metrics

    RDouble2D * lineData1,* lineData2;

    //! calculate mesh transform metrics : first grid delta

    RDouble4D &gridDeltaFirst = * new RDouble4D(Range(1, 9), Range(1, ni2), Range(1, nj2), Range(1, nk2), fortranArray);

    //! calculate mesh transform metrics : first grid delta (I direction)

    lineData1 = new RDouble2D(Range(1, 3), Range(1, ni2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 3), Range(1, ni2), fortranArray);

    for (k2 = k2st; k2 <= k2ed; ++ k2)
    {
        for (j2 = j2st; j2 <= j2ed; ++ j2)
        {
            
            for (i2 = i2st; i2 <= i2ed; ++ i2)
            {
                (*lineData1)(1, i2) = x2(i2, j2, k2);
                (*lineData1)(2, i2) = y2(i2, j2, k2);
                (*lineData1)(3, i2) = z2(i2, j2, k2);
            }
            
            GridDelta(Range(1, 3), Range(1, ni2), lineData1, lineData2);

            for (i2 = i2st; i2 <= i2ed; ++ i2)
            {
                gridDeltaFirst(1, i2, j2, k2) = (*lineData2)(1, i2);
                gridDeltaFirst(2, i2, j2, k2) = (*lineData2)(2, i2);
                gridDeltaFirst(3, i2, j2, k2) = (*lineData2)(3, i2);
            }

        }
    }

    delete lineData1;
    delete lineData2;

    //! calculate mesh transform metrics : first grid delta (J direction)

    lineData1 = new RDouble2D(Range(1, 3), Range(1, nj2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 3), Range(1, nj2), fortranArray);

    for (k2 = k2st; k2 <= k2ed; ++ k2)
    {
        for (i2 = i2st; i2 <= i2ed; ++ i2)
        {
            
            for (j2 = j2st; j2 <= j2ed; ++ j2)
            {
                (*lineData1)(1, j2) = x2(i2, j2, k2);
                (*lineData1)(2, j2) = y2(i2, j2, k2);
                (*lineData1)(3, j2) = z2(i2, j2, k2);
            }
            
            GridDelta(Range(1, 3), Range(1, nj2), lineData1, lineData2);

            for (j2 = j2st; j2 <= j2ed; ++ j2)
            {
                gridDeltaFirst(4, i2, j2, k2) = (*lineData2)(1, j2);
                gridDeltaFirst(5, i2, j2, k2) = (*lineData2)(2, j2);
                gridDeltaFirst(6, i2, j2, k2) = (*lineData2)(3, j2);
            }

        }
    }

    delete lineData1;
    delete lineData2;

    //! calculate mesh transform metrics : first grid delta (K direction)

    lineData1 = new RDouble2D(Range(1, 3), Range(1, nk2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 3), Range(1, nk2), fortranArray);

    for (j2 = j2st; j2 <= j2ed; ++ j2)
    {
        for (i2 = i2st; i2 <= i2ed; ++ i2)
        {
            
            for (k2 = k2st; k2 <= k2ed; ++ k2)
            {
                (*lineData1)(1, k2) = x2(i2, j2, k2);
                (*lineData1)(2, k2) = y2(i2, j2, k2);
                (*lineData1)(3, k2) = z2(i2, j2, k2);
            }
            
            GridDelta(Range(1, 3), Range(1, nk2), lineData1, lineData2);

            for (k2 = k2st; k2 <= k2ed; ++ k2)
            {
                gridDeltaFirst(7, i2, j2, k2) = (*lineData2)(1, k2);
                gridDeltaFirst(8, i2, j2, k2) = (*lineData2)(2, k2);
                gridDeltaFirst(9, i2, j2, k2) = (*lineData2)(3, k2);
            }

        }
    }

    delete lineData1;
    delete lineData2;

    //! calculate mesh transform metrics : second grid delta

    RDouble4D &gridDeltaSecond = * new RDouble4D(Range(1, 18), Range(1, ni2), Range(1, nj2), Range(1, nk2), fortranArray);
    
    //! calculate mesh transform metrics : second grid delta (I direction)

    lineData1 = new RDouble2D(Range(1, 6), Range(1, ni2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 6), Range(1, ni2), fortranArray);

    for (k2 = k2st; k2 <= k2ed; ++ k2)
    {
        for (j2 = j2st; j2 <= j2ed; ++ j2)
        {
            
            for (i2 = i2st; i2 <= i2ed; ++ i2)
            {
                (*lineData1)(1, i2) = half * (z2(i2, j2, k2) * gridDeltaFirst(8, i2, j2, k2) - y2(i2, j2, k2) * gridDeltaFirst(9, i2, j2, k2));
                (*lineData1)(2, i2) = half * (x2(i2, j2, k2) * gridDeltaFirst(9, i2, j2, k2) - z2(i2, j2, k2) * gridDeltaFirst(7, i2, j2, k2));
                (*lineData1)(3, i2) = half * (y2(i2, j2, k2) * gridDeltaFirst(7, i2, j2, k2) - x2(i2, j2, k2) * gridDeltaFirst(8, i2, j2, k2));
                (*lineData1)(4, i2) = half * (y2(i2, j2, k2) * gridDeltaFirst(6, i2, j2, k2) - z2(i2, j2, k2) * gridDeltaFirst(5, i2, j2, k2));
                (*lineData1)(5, i2) = half * (z2(i2, j2, k2) * gridDeltaFirst(4, i2, j2, k2) - x2(i2, j2, k2) * gridDeltaFirst(6, i2, j2, k2));
                (*lineData1)(6, i2) = half * (x2(i2, j2, k2) * gridDeltaFirst(5, i2, j2, k2) - y2(i2, j2, k2) * gridDeltaFirst(4, i2, j2, k2));
            }
            
            GridDelta(Range(1, 6), Range(1, ni2), lineData1, lineData2);

            for (i2 = i2st; i2 <= i2ed; ++ i2)
            {
                for (int n = 1; n <= 6; ++ n)
                {
                    gridDeltaSecond(n, i2, j2, k2) = (*lineData2)(n, i2);
                }
            }

        }
    }

    delete lineData1;
    delete lineData2;

    //! calculate mesh transform metrics : second grid delta (J direction)

    lineData1 = new RDouble2D(Range(1, 6), Range(1, nj2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 6), Range(1, nj2), fortranArray);

    for (k2 = k2st; k2 <= k2ed; ++ k2)
    {
        for (i2 = i2st; i2 <= i2ed; ++ i2)
        {
            
            for (j2 = j2st; j2 <= j2ed; ++ j2)
            {
                (*lineData1)(1, j2) = half * (y2(i2, j2, k2) * gridDeltaFirst(9, i2, j2, k2) - z2(i2, j2, k2) * gridDeltaFirst(8, i2, j2, k2));
                (*lineData1)(2, j2) = half * (z2(i2, j2, k2) * gridDeltaFirst(7, i2, j2, k2) - x2(i2, j2, k2) * gridDeltaFirst(9, i2, j2, k2));
                (*lineData1)(3, j2) = half * (x2(i2, j2, k2) * gridDeltaFirst(8, i2, j2, k2) - y2(i2, j2, k2) * gridDeltaFirst(7, i2, j2, k2));
                (*lineData1)(4, j2) = half * (z2(i2, j2, k2) * gridDeltaFirst(2, i2, j2, k2) - y2(i2, j2, k2) * gridDeltaFirst(3, i2, j2, k2));
                (*lineData1)(5, j2) = half * (x2(i2, j2, k2) * gridDeltaFirst(3, i2, j2, k2) - z2(i2, j2, k2) * gridDeltaFirst(1, i2, j2, k2));
                (*lineData1)(6, j2) = half * (y2(i2, j2, k2) * gridDeltaFirst(1, i2, j2, k2) - x2(i2, j2, k2) * gridDeltaFirst(2, i2, j2, k2));
            }
            
            GridDelta(Range(1, 6), Range(1, nj2), lineData1, lineData2);

            for (j2 = j2st; j2 <= j2ed; ++ j2)
            {
                for (int n = 1; n <= 6; ++ n)
                {
                    gridDeltaSecond(n + 6, i2, j2, k2) = (*lineData2)(n, j2);
                }
            }

        }
    }

    delete lineData1;
    delete lineData2;

    //! calculate mesh transform metrics : second grid delta (K direction)

    lineData1 = new RDouble2D(Range(1, 6), Range(1, nk2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 6), Range(1, nk2), fortranArray);

    for (j2 = j2st; j2 <= j2ed; ++ j2)
    {
        for (i2 = i2st; i2 <= i2ed; ++ i2)
        {
            for (k2 = k2st; k2 <= k2ed; ++ k2)
            {
                (*lineData1)(1, k2) = half * (z2(i2, j2, k2) * gridDeltaFirst(5, i2, j2, k2) - y2(i2, j2, k2) * gridDeltaFirst(6, i2, j2, k2));
                (*lineData1)(2, k2) = half * (x2(i2, j2, k2) * gridDeltaFirst(6, i2, j2, k2) - z2(i2, j2, k2) * gridDeltaFirst(4, i2, j2, k2));
                (*lineData1)(3, k2) = half * (y2(i2, j2, k2) * gridDeltaFirst(4, i2, j2, k2) - x2(i2, j2, k2) * gridDeltaFirst(5, i2, j2, k2));
                (*lineData1)(4, k2) = half * (y2(i2, j2, k2) * gridDeltaFirst(3, i2, j2, k2) - z2(i2, j2, k2) * gridDeltaFirst(2, i2, j2, k2));
                (*lineData1)(5, k2) = half * (z2(i2, j2, k2) * gridDeltaFirst(1, i2, j2, k2) - x2(i2, j2, k2) * gridDeltaFirst(3, i2, j2, k2));
                (*lineData1)(6, k2) = half * (x2(i2, j2, k2) * gridDeltaFirst(2, i2, j2, k2) - y2(i2, j2, k2) * gridDeltaFirst(1, i2, j2, k2));
            }
            
            GridDelta(Range(1, 6), Range(1, nk2), lineData1, lineData2);

            for (k2 = k2st; k2 <= k2ed; ++ k2)
            {
                for (int n = 1; n <= 6; ++ n)
                {
                    gridDeltaSecond(n + 12, i2, j2, k2) = (*lineData2)(n, k2);
                }
            }

        }
    }

    delete lineData1;
    delete lineData2;

    delete &gridDeltaFirst;

    //! calculate mesh transform metrics by second grid delta

    RDouble4D &xfn  = *(this->GetFaceNormalX());
    RDouble4D &yfn  = *(this->GetFaceNormalY());
    RDouble4D &zfn  = *(this->GetFaceNormalZ());
    RDouble4D &area = *(this->GetFaceArea());
    RDouble4D &xfv  = *(this->GetFaceVectorX());
    RDouble4D &yfv  = *(this->GetFaceVectorY());
    RDouble4D &zfv  = *(this->GetFaceVectorZ());

    xfn = 0.0;
    yfn = 0.0;
    zfn = 0.0;
    area= 0.0;
    xfv = 0.0;
    yfv = 0.0;
    zfv = 0.0;

    RDouble vkc[3],vet[3],vct[3],vkcn[4],vetn[4],vctn[4];

    //! calculate mesh transform metrics by second grid delta (I direction)

    for (k = 1; k <= nk - 1; ++ k)
    {
        k2 = 2 * k;

        for (j = 1; j <= nj - 1; ++ j)
        {
            j2 = 2 * j;

            for (i = 1; i <= ni; ++ i)
            {

                i2 = 2 * i - 1;

                vkc[0] = gridDeltaSecond(7, i2, j2, k2) + gridDeltaSecond(13, i2, j2, k2);
                vkc[1] = gridDeltaSecond(8, i2, j2, k2) + gridDeltaSecond(14, i2, j2, k2);
                vkc[2] = gridDeltaSecond(9, i2, j2, k2) + gridDeltaSecond(15, i2, j2, k2);

                xfv(i, j, k, IDX::IKES) = vkc[0];
                yfv(i, j, k, IDX::IKES) = vkc[1];
                zfv(i, j, k, IDX::IKES) = vkc[2];

                NormalizeVector(vkc,vkcn, 3);

                area(i, j, k, IDX::IKES) = vkcn[3];
                xfn (i, j, k, IDX::IKES) = vkcn[0];
                yfn (i, j, k, IDX::IKES) = vkcn[1];
                zfn (i, j, k, IDX::IKES) = vkcn[2];

            }
        }
    }

    //! calculate mesh transform metrics by second grid delta (J direction)

    for (k = 1; k <= nk - 1; ++ k)
    {
        k2 = 2 * k;

        for (i = 1; i <= ni - 1; ++ i)
        {
            i2 = 2 * i;

            for (j = 1; j <= nj; ++ j)
            {

                j2 = 2 * j - 1;

                vet[0] = gridDeltaSecond(1, i2, j2, k2) + gridDeltaSecond(16, i2, j2, k2);
                vet[1] = gridDeltaSecond(2, i2, j2, k2) + gridDeltaSecond(17, i2, j2, k2);
                vet[2] = gridDeltaSecond(3, i2, j2, k2) + gridDeltaSecond(18, i2, j2, k2);

                xfv(i, j, k, IDX::IETA) = vet[0];
                yfv(i, j, k, IDX::IETA) = vet[1];
                zfv(i, j, k, IDX::IETA) = vet[2];

                NormalizeVector(vet,vetn, 3);

                area(i, j, k, IDX::IETA) = vetn[3];
                xfn (i, j, k, IDX::IETA) = vetn[0];
                yfn (i, j, k, IDX::IETA) = vetn[1];
                zfn (i, j, k, IDX::IETA) = vetn[2];

            }
        }
    }

    //! calculate mesh transform metrics by second grid delta (K direction)

    for (j = 1; j <= nj - 1; ++ j)
    {
        j2 = 2 * j;

        for (i = 1; i <= ni - 1; ++ i)
        {
            i2 = 2 * i;

            for (k = 1; k <= nk; ++ k)
            {

                k2 = 2 * k - 1;

                vct[0] = gridDeltaSecond(4, i2, j2, k2) + gridDeltaSecond(10, i2, j2, k2);
                vct[1] = gridDeltaSecond(5, i2, j2, k2) + gridDeltaSecond(11, i2, j2, k2);
                vct[2] = gridDeltaSecond(6, i2, j2, k2) + gridDeltaSecond(12, i2, j2, k2);

                xfv(i, j, k, IDX::IZTA) = vct[0];
                yfv(i, j, k, IDX::IZTA) = vct[1];
                zfv(i, j, k, IDX::IZTA) = vct[2];

                NormalizeVector(vct,vctn, 3);

                area(i, j, k, IDX::IZTA) = vctn[3];
                xfn (i, j, k, IDX::IZTA) = vctn[0];
                yfn (i, j, k, IDX::IZTA) = vctn[1];
                zfn (i, j, k, IDX::IZTA) = vctn[2];

            }
        }
    }

    //! calculate mesh transform Jacobian : third grid delta

    RDouble4D &gridDeltaThird = * new RDouble4D(Range(1, 3), Range(1, ni2), Range(1, nj2), Range(1, nk2), fortranArray);

    //! calculate mesh transform Jacobian : third grid delta (I direction)

    lineData1 = new RDouble2D(Range(1, 1), Range(1, ni2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 1), Range(1, ni2), fortranArray);

    for (k2 = k2st; k2 <= k2ed; ++ k2)
    {
        for (j2 = j2st; j2 <= j2ed; ++ j2)
        {
            
            for (i2 = i2st; i2 <= i2ed; ++ i2)
            {
                (*lineData1)(1, i2) =   x2(i2, j2, k2) * (gridDeltaSecond(7, i2, j2, k2) + gridDeltaSecond(13, i2, j2, k2))
                                     + y2(i2, j2, k2) * (gridDeltaSecond(8, i2, j2, k2) + gridDeltaSecond(14, i2, j2, k2))
                                     + z2(i2, j2, k2) * (gridDeltaSecond(9, i2, j2, k2) + gridDeltaSecond(15, i2, j2, k2));
            }
            
            GridDelta(Range(1, 1), Range(1, ni2), lineData1, lineData2);

            for (i2 = i2st; i2 <= i2ed; ++ i2)
            {
                gridDeltaThird(1, i2, j2, k2) = (*lineData2)(1, i2);
            }

        }
    }

    delete lineData1;
    delete lineData2;

    //! calculate mesh transform Jacobian : third grid delta (J direction)

    lineData1 = new RDouble2D(Range(1, 1), Range(1, nj2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 1), Range(1, nj2), fortranArray);

    for (k2 = k2st; k2 <= k2ed; ++ k2)
    {
        for (i2 = i2st; i2 <= i2ed; ++ i2)
        {
            
            for (j2 = j2st; j2 <= j2ed; ++ j2)
            {
                (*lineData1)(1, j2) =   x2(i2, j2, k2) * (gridDeltaSecond(1, i2, j2, k2) + gridDeltaSecond(16, i2, j2, k2))
                                     + y2(i2, j2, k2) * (gridDeltaSecond(2, i2, j2, k2) + gridDeltaSecond(17, i2, j2, k2))
                                     + z2(i2, j2, k2) * (gridDeltaSecond(3, i2, j2, k2) + gridDeltaSecond(18, i2, j2, k2));
            }
            
            GridDelta(Range(1, 1), Range(1, nj2), lineData1, lineData2);

            for (j2 = j2st; j2 <= j2ed; ++ j2)
            {
                gridDeltaThird(2, i2, j2, k2) = (*lineData2)(1, j2);
            }

        }
    }

    delete lineData1;
    delete lineData2;

    //! calculate mesh transform Jacobian : third grid delta (K direction)

    lineData1 = new RDouble2D(Range(1, 1), Range(1, nk2), fortranArray);
    lineData2 = new RDouble2D(Range(1, 1), Range(1, nk2), fortranArray);

    for (j2 = j2st; j2 <= j2ed; ++ j2)
    {
        for (i2 = i2st; i2 <= i2ed; ++ i2)
        {
            
            for (k2 = k2st; k2 <= k2ed; ++ k2)
            {
                (*lineData1)(1, k2) =   x2(i2, j2, k2) * (gridDeltaSecond(4, i2, j2, k2) + gridDeltaSecond(10, i2, j2, k2))
                                     + y2(i2, j2, k2) * (gridDeltaSecond(5, i2, j2, k2) + gridDeltaSecond(11, i2, j2, k2))
                                     + z2(i2, j2, k2) * (gridDeltaSecond(6, i2, j2, k2) + gridDeltaSecond(12, i2, j2, k2));
            }
            
            GridDelta(Range(1, 1), Range(1, nk2), lineData1, lineData2);

            for (k2 = k2st; k2 <= k2ed; ++ k2)
            {
                gridDeltaThird(3, i2, j2, k2) = (*lineData2)(1, k2);
            }

        }
    }

    delete lineData1;
    delete lineData2;
    
    delete &x2;
    delete &y2;
    delete &z2;

    delete &gridDeltaSecond;

    //! calculate mesh transform Jacobian by third grid delta

    RDouble3D &vol = * this->GetCellVolume();
    int ZoneID = this->GetZoneID();
    int imin, jmin, kmin, imax, jmax, kmax;
    RDouble vv,volmin,volmax;

    volmin =   LARGE;
    volmax = - LARGE;

    std::ostringstream oss;

    for (k = 1; k <= nk - 1; ++ k)
    {
        k2 = 2 * k;
        for (j = 1; j <= nj - 1; ++ j)
        {
            j2 = 2 * j;
            for (i = 1; i <= ni - 1; ++ i)
            {
                i2 = 2 * i;
                vv = third * (gridDeltaThird(1, i2, j2, k2) + gridDeltaThird(2, i2, j2, k2) + gridDeltaThird(3, i2, j2, k2));
                vol(i, j, k) = vv;

                if (vv <= zero)
                {
                    //oss << " vol < 0 " << ZoneID << " " << i << " " << j << " " << k << vv << "\n";
                    oss << " vol < 0 " << " ZoneID = " << ZoneID << " i = " << i << " j = " << j << " k = " << k << " vol = " << vv << "\n";
                }

                if (vv < volmin)
                {
                    volmin = vv;
                    imin = i;
                    jmin = j;
                    kmin = k;
                }

                if (vv > volmax)
                {
                    volmax = vv;
                    imax = i;
                    jmax = j;
                    kmax = k;
                }

            }
        }
    }

    delete &gridDeltaThird;

    GhostCell3D(vol, ni, nj, nk);
    GhostMetrics3D();
    ComputeCellCenter();
    StreamToActionKey(actkey, oss);
}

void StructGrid::GridDelta(const Range &NP, const Range &NI, RDouble2D *lineData1, RDouble2D *lineData2)
{
    int nst = NP.first();
    int ned = NP.last();
    int ist = NI.first();
    int ied = NI.last();
    
    int i, n;
    
    i = ist;
    for (n = nst; n <= ned; ++ n)
    {
        (*lineData2)(n, i) = - 3.0 * (*lineData1)(n, i) + 4.0 * (*lineData1)(n, i + 1) - (*lineData1)(n, i + 2);
    }

    i = ist + 1;
    for (n = nst; n <= ned; ++ n)
    {
        (*lineData2)(n, i) = (*lineData1)(n, i + 1) - (*lineData1)(n, i - 1);
    }

    i = ist + 2;
    for (n = nst; n <= ned; ++ n)
    {
        (*lineData2)(n, i) = (*lineData1)(n, i + 1) - (*lineData1)(n, i - 1);
    }

    for (i = ist + 3; i <= ied - 3; ++ i)
    {
        for (n = nst; n <= ned; ++ n)
        {
            
            //! 2nd order
            //(*lineData2)(n, i) = (*lineData1)(n, i + 1) - (*lineData1)(n, i - 1);
            
            //! 4th order
            (*lineData2)(n, i) =   9.0 * ((*lineData1)(n, i + 1) - (*lineData1)(n, i - 1)) /  8.0
                               -         ((*lineData1)(n, i + 3) - (*lineData1)(n, i - 3)) / 24.0;
        
        }
    }

    i = ied - 2;
    for (n = nst; n <= ned; ++ n)
    {
        (*lineData2)(n, i) = (*lineData1)(n, i + 1) - (*lineData1)(n, i - 1);
    }

    i = ied - 1;
    for (n = nst; n <= ned; ++ n)
    {
        (*lineData2)(n, i) = (*lineData1)(n, i + 1) - (*lineData1)(n, i - 1);
    }

    i = ied;
    for (n = nst; n <= ned; ++ n)
    {
        (*lineData2)(n, i) = 3.0 * (*lineData1)(n, i) - 4.0 * (*lineData1)(n, i - 1) + (*lineData1)(n, i - 2);
    }

}

void StructGrid::ComputeMetricsStructHighOrder2D_HDCS()
{
    RDouble4D &xfn  = *(this->GetFaceNormalX());
    RDouble4D &yfn  = *(this->GetFaceNormalY());
    RDouble4D &zfn  = *(this->GetFaceNormalZ());
    RDouble4D &area = *(this->GetFaceArea());
    RDouble4D &xfv  = *(this->GetFaceVectorX());
    RDouble4D &yfv  = *(this->GetFaceVectorY());
    RDouble4D &zfv  = *(this->GetFaceVectorZ());

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();

    zfn = 0.0;
    zfv = 0.0;

    int ni = this->GetNI();
    int nj = this->GetNJ();
    const int k = 1;
    
    const RDouble alfa_sixth =  64.0 / 45.0;
    const RDouble    a_sixth = - 2.0 / 9.0;
    const RDouble    b_sixth =   1.0 / 180.0;

    for (int i = -1; i <= ni + 2; ++ i)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            RDouble xeta, yeta;

            RDouble x_CR2 = half * (xx(i, j + 2, k) + xx(i, j + 3, k));
            RDouble x_CR1 = half * (xx(i, j + 1, k) + xx(i, j + 2, k));
            RDouble x_CL1 = half * (xx(i, j    , k) + xx(i, j - 1, k));
            RDouble x_CL2 = half * (xx(i, j - 1, k) + xx(i, j - 2, k));

            RDouble x_fR = xx(i, j + 1, k);
            RDouble x_fL = xx(i, j    , k);

            RDouble y_CR2 = half * (yy(i, j + 2, k) + yy(i, j + 3, k));
            RDouble y_CR1 = half * (yy(i, j + 1, k) + yy(i, j + 2, k));
            RDouble y_CL1 = half * (yy(i, j    , k) + yy(i, j - 1, k));
            RDouble y_CL2 = half * (yy(i, j - 1, k) + yy(i, j - 2, k));

            RDouble y_fR = yy(i, j + 1, k);
            RDouble y_fL = yy(i, j    , k);

            xeta = alfa_sixth * (x_fR - x_fL) + a_sixth * (x_CR1 - x_CL1) + b_sixth * (x_CR2 - x_CL2);
            yeta = alfa_sixth * (y_fR - y_fL) + a_sixth * (y_CR1 - y_CL1) + b_sixth * (y_CR2 - y_CL2);

            RDouble ds = sqrt(xeta * xeta + yeta * yeta);

            xfv(i, j, k, 1) =   yeta;
            yfv(i, j, k, 1) = - xeta;

            area(i, j, k, 1) = ds;

            xfn(i, j, k, 1) =   yeta / (ds + TINY);
            yfn(i, j, k, 1) = - xeta / (ds + TINY);
        }
    }

    for (int j = -1; j <= nj + 2; ++ j)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            RDouble xkxi, ykxi;

            RDouble x_CR2 = half * (xx(i + 2, j, k) + xx(i + 3, j, k));
            RDouble x_CR1 = half * (xx(i + 1, j, k) + xx(i + 2, j, k));
            RDouble x_CL1 = half * (xx(i    , j, k) + xx(i - 1, j, k));
            RDouble x_CL2 = half * (xx(i - 1, j, k) + xx(i - 2, j, k));

            RDouble x_fR = xx(i + 1, j, k);
            RDouble x_fL = xx(i    , j, k);

            RDouble y_CR2 = half * (yy(i + 2, j, k) + yy(i + 3, j, k));
            RDouble y_CR1 = half * (yy(i + 1, j, k) + yy(i + 2, j, k));
            RDouble y_CL1 = half * (yy(i    , j, k) + yy(i - 1, j, k));
            RDouble y_CL2 = half * (yy(i - 1, j, k) + yy(i - 2, j, k));

            RDouble y_fR = yy(i + 1, j, k);
            RDouble y_fL = yy(i    , j, k);

            xkxi = alfa_sixth * (x_fR - x_fL) + a_sixth * (x_CR1 - x_CL1) + b_sixth * (x_CR2 - x_CL2);
            ykxi = alfa_sixth * (y_fR - y_fL) + a_sixth * (y_CR1 - y_CL1) + b_sixth * (y_CR2 - y_CL2);

            RDouble ds = sqrt(xkxi * xkxi + ykxi * ykxi);

            xfv(i, j, k, 2) = - ykxi;
            yfv(i, j, k, 2) =   xkxi;

            area(i, j, k, 2) = ds;

            xfn(i, j, k, 2) = - ykxi / (ds + TINY);
            yfn(i, j, k, 2) =   xkxi / (ds + TINY);
        }
    }
    //!=
    RDouble3D &vol = *this->GetCellVolume();

    RDouble3D &xfaceI = *new RDouble3D(Range(-1, ni + 2), Range(1, nj - 1), Range(1, 1), fortranArray);
    RDouble3D &yfaceI = *new RDouble3D(Range(-1, ni + 2), Range(1, nj - 1), Range(1, 1), fortranArray);

    RDouble3D &xfaceJ = *new RDouble3D(Range(1, ni - 1), Range(-1, nj + 2), Range(1, 1), fortranArray);
    RDouble3D &yfaceJ = *new RDouble3D(Range(1, ni - 1), Range(-1, nj + 2), Range(1, 1), fortranArray);

    for (int j = 1; j <= nj - 1; ++ j)
    {
        for (int i = -1; i <= ni + 2; ++ i)
        {
            RDouble x1 = xx(i, j,  k);
            RDouble x2 = xx(i, j + 1, k);

            RDouble y1 = yy(i, j,  k);
            RDouble y2 = yy(i, j + 1, k);
            
            xfaceI(i, j, k) = half * (x1 + x2);
            yfaceI(i, j, k) = half * (y1 + y2);
        }
    }
    for (int i = 1; i <= ni - 1; ++ i)
    {
        for (int j = -1; j <= nj + 2; ++ j)
        {
            RDouble x1 = xx(i,  j, k);
            RDouble x2 = xx(i + 1, j, k);

            RDouble y1 = yy(i,  j, k);
            RDouble y2 = yy(i + 1, j, k);
            
            xfaceJ(i, j, k) = half * (x1 + x2);
            yfaceJ(i, j, k) = half * (y1 + y2);
        }
    }

    int ZoneID = this->GetZoneID();
    std::ostringstream oss;
    for (int i = 1; i <= ni - 1; ++ i)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            RDouble vol_kxi;
            RDouble vol_eta;

            RDouble kxi_fL = xfaceI(i    , j, k) * xfv(i    , j, k, 1) + yfaceI(i    , j, k) * yfv(i    , j, k, 1);
            RDouble kxi_fR = xfaceI(i + 1, j, k) * xfv(i + 1, j, k, 1) + yfaceI(i + 1, j, k) * yfv(i + 1, j, k, 1);

            RDouble kxi_CL2 = 0.25 * (xfaceI(i - 2, j, k) + xfaceI(i - 1, j, k)) * (xfv(i - 2, j, k, 1) + xfv(i - 1, j, k, 1)) + 0.25 * (yfaceI(i - 2, j, k) + yfaceI(i - 1, j, k)) * (yfv(i - 2, j, k, 1) + yfv(i - 1, j, k, 1));
            RDouble kxi_CL1 = 0.25 * (xfaceI(i - 1, j, k) + xfaceI(i    , j, k)) * (xfv(i - 1, j, k, 1) + xfv(i    , j, k, 1)) + 0.25 * (yfaceI(i - 1, j, k) + yfaceI(i    , j, k)) * (yfv(i - 1, j, k, 1) + yfv(i    , j, k, 1));
            RDouble kxi_CR1 = 0.25 * (xfaceI(i + 1, j, k) + xfaceI(i + 2, j, k)) * (xfv(i + 1, j, k, 1) + xfv(i + 2, j, k, 1)) + 0.25 * (yfaceI(i + 1, j, k) + yfaceI(i + 2, j, k)) * (yfv(i + 1, j, k, 1) + yfv(i + 2, j, k, 1));
            RDouble kxi_CR2 = 0.25 * (xfaceI(i + 2, j, k) + xfaceI(i + 3, j, k)) * (xfv(i + 2, j, k, 1) + xfv(i + 3, j, k, 1)) + 0.25 * (yfaceI(i + 2, j, k) + yfaceI(i + 3, j, k)) * (yfv(i + 2, j, k, 1) + yfv(i + 3, j, k, 1));

            RDouble eta_fL = xfaceJ(i, j    , k) * xfv(i, j    , k, 2) + yfaceJ(i, j    , k) * yfv(i, j    , k, 2);
            RDouble eta_fR = xfaceJ(i, j + 1, k) * xfv(i, j + 1, k, 2) + yfaceJ(i, j + 1, k) * yfv(i, j + 1, k, 2);

            RDouble eta_CL2 = 0.25 * (xfaceJ(i, j - 2, k) + xfaceJ(i, j - 1, k)) * (xfv(i, j - 2, k, 2) + xfv(i, j - 1, k, 2)) + 0.25 * (yfaceJ(i, j - 2, k) + yfaceJ(i, j - 1, k)) * (yfv(i, j - 2, k, 2) + yfv(i, j - 1, k, 2));
            RDouble eta_CL1 = 0.25 * (xfaceJ(i, j - 1, k) + xfaceJ(i, j    , k)) * (xfv(i, j - 1, k, 2) + xfv(i, j    , k, 2)) + 0.25 * (yfaceJ(i, j - 1, k) + yfaceJ(i, j    , k)) * (yfv(i, j - 1, k, 2) + yfv(i, j    , k, 2));
            RDouble eta_CR1 = 0.25 * (xfaceJ(i, j + 1, k) + xfaceJ(i, j + 2, k)) * (xfv(i, j + 1, k, 2) + xfv(i, j + 2, k, 2)) + 0.25 * (yfaceJ(i, j + 1, k) + yfaceJ(i, j + 2, k)) * (yfv(i, j + 1, k, 2) + yfv(i, j + 2, k, 2));
            RDouble eta_CR2 = 0.25 * (xfaceJ(i, j + 2, k) + xfaceJ(i, j + 3, k)) * (xfv(i, j + 2, k, 2) + xfv(i, j + 3, k, 2)) + 0.25 * (yfaceJ(i, j + 2, k) + yfaceJ(i, j + 3, k)) * (yfv(i, j + 2, k, 2) + yfv(i, j + 3, k, 2));


            vol_kxi= alfa_sixth * (kxi_fR - kxi_fL) + a_sixth * (kxi_CR1 - kxi_CL1) + b_sixth * (kxi_CR2 - kxi_CL2);
            vol_eta= alfa_sixth * (eta_fR - eta_fL) + a_sixth * (eta_CR1 - eta_CL1) + b_sixth * (eta_CR2 - eta_CL2);

            RDouble vol_FD = half * (vol_kxi + vol_eta);

            //!RDouble vol_FV = vol(i, j, k);
            //!
            //!if (vol_FD <= zero || vol_FD > 1.2 * vol_FV || vol_FD < 0.8 * vol_FV)
            //!{
            //!    vol_FD = vol_FV;
            //!}
            if (vol_FD < SMALL)
            {
                oss << " vol < 0 " << " ZoneID = " << ZoneID << " i = " << i << " j = " << j << " vol = " << vol_FD << "\n";
            }
            vol(i, j, k) = vol_FD;
        }
    }
    ComputeCellCenter();

    delete &xfaceI;
    delete &yfaceI;
    delete &xfaceJ;
    delete &yfaceJ;
}

void StructGrid::ComputeMetricsStructHighOrder3D_HDCS()
{
    RDouble4D &xfn  = *(this->GetFaceNormalX());
    RDouble4D &yfn  = *(this->GetFaceNormalY());
    RDouble4D &zfn  = *(this->GetFaceNormalZ());
    RDouble4D &area = *(this->GetFaceArea());
    RDouble4D &xfv  = *(this->GetFaceVectorX());
    RDouble4D &yfv  = *(this->GetFaceVectorY());
    RDouble4D &zfv  = *(this->GetFaceVectorZ());

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    const RDouble alfa_sixth =  64.0 / 45.0;
    const RDouble    a_sixth = - 2.0 / 9.0;
    const RDouble    b_sixth =   1.0 / 180.0;

    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();  
    RDouble3D &xA    = *new RDouble3D(Range(1, ni - 1), Range(-1, nj + 2), Range(-1, nk + 2), fortranArray);
    RDouble3D &yA    = *new RDouble3D(Range(1, ni - 1), Range(-1, nj + 2), Range(-1, nk + 2), fortranArray);
    RDouble3D &zA    = *new RDouble3D(Range(1, ni - 1), Range(-1, nj + 2), Range(-1, nk + 2), fortranArray);
    RDouble3D &xkxi  = *new RDouble3D(Range(1, ni - 1), Range(-1, nj + 2), Range(-1, nk + 2), fortranArray);
    RDouble3D &ykxi  = *new RDouble3D(Range(1, ni - 1), Range(-1, nj + 2), Range(-1, nk + 2), fortranArray);
    RDouble3D &zkxi  = *new RDouble3D(Range(1, ni - 1), Range(-1, nj + 2), Range(-1, nk + 2), fortranArray);

    RDouble3D &xB    = *new RDouble3D(Range(-1, ni + 2), Range(1, nj - 1), Range(-1, nk + 2), fortranArray);
    RDouble3D &yB    = *new RDouble3D(Range(-1, ni + 2), Range(1, nj - 1), Range(-1, nk + 2), fortranArray);
    RDouble3D &zB    = *new RDouble3D(Range(-1, ni + 2), Range(1, nj - 1), Range(-1, nk + 2), fortranArray);
    RDouble3D &xeta  = *new RDouble3D(Range(-1, ni + 2), Range(1, nj - 1), Range(-1, nk + 2), fortranArray);
    RDouble3D &yeta  = *new RDouble3D(Range(-1, ni + 2), Range(1, nj - 1), Range(-1, nk + 2), fortranArray);
    RDouble3D &zeta  = *new RDouble3D(Range(-1, ni + 2), Range(1, nj - 1), Range(-1, nk + 2), fortranArray);

    RDouble3D &xC    = *new RDouble3D(Range(-1, ni + 2), Range(-1, nj + 2), Range(1, nk - 1), fortranArray);
    RDouble3D &yC    = *new RDouble3D(Range(-1, ni + 2), Range(-1, nj + 2), Range(1, nk - 1), fortranArray);
    RDouble3D &zC    = *new RDouble3D(Range(-1, ni + 2), Range(-1, nj + 2), Range(1, nk - 1), fortranArray);
    RDouble3D &xcta  = *new RDouble3D(Range(-1, ni + 2), Range(-1, nj + 2), Range(1, nk - 1), fortranArray);
    RDouble3D &ycta  = *new RDouble3D(Range(-1, ni + 2), Range(-1, nj + 2), Range(1, nk - 1), fortranArray);
    RDouble3D &zcta  = *new RDouble3D(Range(-1, ni + 2), Range(-1, nj + 2), Range(1, nk - 1), fortranArray);

    for (int i = 1; i <= ni - 1; ++ i)
    {
        for (int j = -1; j <= nj + 2; ++ j)
        {
            for (int k = -1; k <= nk + 2; ++ k)
            {
                xA(i, j, k) = half * (xx(i + 1, j, k) + xx(i, j, k));
                yA(i, j, k) = half * (yy(i + 1, j, k) + yy(i, j, k));
                zA(i, j, k) = half * (zz(i + 1, j, k) + zz(i, j, k));

                RDouble x_CL2 = half * (xx(i - 1, j, k) + xx(i - 2, j, k));
                RDouble x_CL1 = half * (xx(i    , j, k) + xx(i - 1, j, k));
                RDouble x_fL  =         xx(i    , j, k);
                RDouble x_fR  =         xx(i + 1, j, k);
                RDouble x_CR1 = half * (xx(i + 1, j, k) + xx(i + 2, j, k));
                RDouble x_CR2 = half * (xx(i + 2, j, k) + xx(i + 3, j, k));

                RDouble y_CL2 = half * (yy(i - 1, j, k) + yy(i - 2, j, k));
                RDouble y_CL1 = half * (yy(i    , j, k) + yy(i - 1, j, k));
                RDouble y_fL  =         yy(i    , j, k);
                RDouble y_fR  =         yy(i + 1, j, k);
                RDouble y_CR1 = half * (yy(i + 1, j, k) + yy(i + 2, j, k));
                RDouble y_CR2 = half * (yy(i + 2, j, k) + yy(i + 3, j, k));

                RDouble z_CL2 = half * (zz(i - 1, j, k) + zz(i - 2, j, k));
                RDouble z_CL1 = half * (zz(i    , j, k) + zz(i - 1, j, k));
                RDouble z_fL  =         zz(i    , j, k);
                RDouble z_fR  =         zz(i + 1, j, k);
                RDouble z_CR1 = half * (zz(i + 1, j, k) + zz(i + 2, j, k));
                RDouble z_CR2 = half * (zz(i + 2, j, k) + zz(i + 3, j, k));

                xkxi(i, j, k) = alfa_sixth * (x_fR - x_fL) + a_sixth * (x_CR1 - x_CL1) + b_sixth * (x_CR2 - x_CL2);
                ykxi(i, j, k) = alfa_sixth * (y_fR - y_fL) + a_sixth * (y_CR1 - y_CL1) + b_sixth * (y_CR2 - y_CL2);
                zkxi(i, j, k) = alfa_sixth * (z_fR - z_fL) + a_sixth * (z_CR1 - z_CL1) + b_sixth * (z_CR2 - z_CL2);
            }
        }
    }

    for (int j = 1; j <= nj - 1; ++ j)
    {
        for (int k = -1; k <= nk + 2; ++ k)
        {
            for (int i = -1; i <= ni + 2; ++ i)
            {
                xB(i, j, k) = half * (xx(i, j + 1, k) + xx(i, j, k));
                yB(i, j, k) = half * (yy(i, j + 1, k) + yy(i, j, k));
                zB(i, j, k) = half * (zz(i, j + 1, k) + zz(i, j, k));

                RDouble x_CL2 = half * (xx(i, j - 1, k) + xx(i, j - 2, k));
                RDouble x_CL1 = half * (xx(i, j    , k) + xx(i, j - 1, k));
                RDouble x_fL  =         xx(i, j    , k);
                RDouble x_fR  =         xx(i, j + 1, k);
                RDouble x_CR1 = half * (xx(i, j + 1, k) + xx(i, j + 2, k));
                RDouble x_CR2 = half * (xx(i, j + 2, k) + xx(i, j + 3, k));

                RDouble y_CL2 = half * (yy(i, j - 1, k) + yy(i, j - 2, k));
                RDouble y_CL1 = half * (yy(i, j    , k) + yy(i, j - 1, k));
                RDouble y_fL  =         yy(i, j    , k);
                RDouble y_fR  =         yy(i, j + 1, k);
                RDouble y_CR1 = half * (yy(i, j + 1, k) + yy(i, j + 2, k));
                RDouble y_CR2 = half * (yy(i, j + 2, k) + yy(i, j + 3, k));

                RDouble z_CL2 = half * (zz(i, j - 1, k) + zz(i, j - 2, k));
                RDouble z_CL1 = half * (zz(i, j    , k) + zz(i, j - 1, k));
                RDouble z_fL  =         zz(i, j    , k);
                RDouble z_fR  =         zz(i, j + 1, k);
                RDouble z_CR1 = half * (zz(i, j + 1, k) + zz(i, j + 2, k));
                RDouble z_CR2 = half * (zz(i, j + 2, k) + zz(i, j + 3, k));

                xeta(i, j, k) = alfa_sixth * (x_fR - x_fL) + a_sixth * (x_CR1 - x_CL1) + b_sixth * (x_CR2 - x_CL2);
                yeta(i, j, k) = alfa_sixth * (y_fR - y_fL) + a_sixth * (y_CR1 - y_CL1) + b_sixth * (y_CR2 - y_CL2);
                zeta(i, j, k) = alfa_sixth * (z_fR - z_fL) + a_sixth * (z_CR1 - z_CL1) + b_sixth * (z_CR2 - z_CL2);
            }
        }
    }
    
    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int i = -1; i <= ni + 2; ++ i)
        {
            for (int j = -1; j <= nj + 2; ++ j)
            {
                xC(i, j, k) = half * (xx(i, j, k + 1) + xx(i, j, k));
                yC(i, j, k) = half * (yy(i, j, k + 1) + yy(i, j, k));
                zC(i, j, k) = half * (zz(i, j, k + 1) + zz(i, j, k));

                RDouble x_CL2 = half * (xx(i, j, k - 1) + xx(i, j, k - 2));
                RDouble x_CL1 = half * (xx(i, j, k   ) + xx(i, j, k - 1));
                RDouble x_fL  =         xx(i, j, k   );
                RDouble x_fR  =         xx(i, j, k + 1);
                RDouble x_CR1 = half * (xx(i, j, k + 1) + xx(i, j, k + 2));
                RDouble x_CR2 = half * (xx(i, j, k + 2) + xx(i, j, k + 3));

                RDouble y_CL2 = half * (yy(i, j, k - 1) + yy(i, j, k - 2));
                RDouble y_CL1 = half * (yy(i, j, k   ) + yy(i, j, k - 1));
                RDouble y_fL  =         yy(i, j, k   );
                RDouble y_fR  =         yy(i, j, k + 1);
                RDouble y_CR1 = half * (yy(i, j, k + 1) + yy(i, j, k + 2));
                RDouble y_CR2 = half * (yy(i, j, k + 2) + yy(i, j, k + 3));

                RDouble z_CL2 = half * (zz(i, j, k - 1) + zz(i, j, k - 2));
                RDouble z_CL1 = half * (zz(i, j, k   ) + zz(i, j, k - 1));
                RDouble z_fL  =         zz(i, j, k   );
                RDouble z_fR  =         zz(i, j, k + 1);
                RDouble z_CR1 = half * (zz(i, j, k + 1) + zz(i, j, k + 2));
                RDouble z_CR2 = half * (zz(i, j, k + 2) + zz(i, j, k + 3));

                xcta(i, j, k) = alfa_sixth * (x_fR - x_fL) + a_sixth * (x_CR1 - x_CL1) + b_sixth * (x_CR2 - x_CL2);
                ycta(i, j, k) = alfa_sixth * (y_fR - y_fL) + a_sixth * (y_CR1 - y_CL1) + b_sixth * (y_CR2 - y_CL2);
                zcta(i, j, k) = alfa_sixth * (z_fR - z_fL) + a_sixth * (z_CR1 - z_CL1) + b_sixth * (z_CR2 - z_CL2);
            }
        }
    }
    //!=
    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = -1; i <= ni + 2; ++ i)
            {
                RDouble KX = 0.0;
                RDouble KY = 0.0;
                RDouble KZ = 0.0;

                RDouble sx_CL2 = 0.25 * (yeta(i, j, k - 1) + yeta(i, j, k - 2)) * (zB(i, j, k - 1) + zB(i, j, k - 2)) - 0.25 * (yB(i, j, k - 1) + yB(i, j, k - 2)) * (zeta(i, j, k - 1) + zeta(i, j, k - 2));
                RDouble sx_CL1 = 0.25 * (yeta(i, j, k   ) + yeta(i, j, k - 1)) * (zB(i, j, k   ) + zB(i, j, k - 1)) - 0.25 * (yB(i, j, k   ) + yB(i, j, k - 1)) * (zeta(i, j, k   ) + zeta(i, j, k - 1));
                RDouble sx_fL  =         yeta(i, j, k   )                      *  zB(i, j, k   )                    -         yB(i, j, k   )                    *  zeta(i, j, k   )                     ;
                RDouble sx_fR  =         yeta(i, j, k + 1)                      *  zB(i, j, k + 1)                    -         yB(i, j, k + 1)                    *  zeta(i, j, k + 1)                     ;
                RDouble sx_CR1 = 0.25 * (yeta(i, j, k + 1) + yeta(i, j, k + 2)) * (zB(i, j, k + 1) + zB(i, j, k + 2)) - 0.25 * (yB(i, j, k + 1) + yB(i, j, k + 2)) * (zeta(i, j, k + 1) + zeta(i, j, k + 2));
                RDouble sx_CR2 = 0.25 * (yeta(i, j, k + 2) + yeta(i, j, k + 3)) * (zB(i, j, k + 2) + zB(i, j, k + 3)) - 0.25 * (yB(i, j, k + 2) + yB(i, j, k + 3)) * (zeta(i, j, k + 2) + zeta(i, j, k + 3));

                RDouble sy_CL2 = 0.25 * (zeta(i, j, k - 1) + zeta(i, j, k - 2)) * (xB(i, j, k - 1) + xB(i, j, k - 2)) - 0.25 * (zB(i, j, k - 1) + zB(i, j, k - 2)) * (xeta(i, j, k - 1) + xeta(i, j, k - 2));
                RDouble sy_CL1 = 0.25 * (zeta(i, j, k   ) + zeta(i, j, k - 1)) * (xB(i, j, k   ) + xB(i, j, k - 1)) - 0.25 * (zB(i, j, k   ) + zB(i, j, k - 1)) * (xeta(i, j, k   ) + xeta(i, j, k - 1));
                RDouble sy_fL  =         zeta(i, j, k   )                      *  xB(i, j, k   )                    -         zB(i, j, k   )                    *  xeta(i, j, k   )                     ;
                RDouble sy_fR  =         zeta(i, j, k + 1)                      *  xB(i, j, k + 1)                    -         zB(i, j, k + 1)                    *  xeta(i, j, k + 1)                     ;
                RDouble sy_CR1 = 0.25 * (zeta(i, j, k + 1) + zeta(i, j, k + 2)) * (xB(i, j, k + 1) + xB(i, j, k + 2)) - 0.25 * (zB(i, j, k + 1) + zB(i, j, k + 2)) * (xeta(i, j, k + 1) + xeta(i, j, k + 2));
                RDouble sy_CR2 = 0.25 * (zeta(i, j, k + 2) + zeta(i, j, k + 3)) * (xB(i, j, k + 2) + xB(i, j, k + 3)) - 0.25 * (zB(i, j, k + 2) + zB(i, j, k + 3)) * (xeta(i, j, k + 2) + xeta(i, j, k + 3));

                RDouble sz_CL2 = 0.25 * (xeta(i, j, k - 1) + xeta(i, j, k - 2)) * (yB(i, j, k - 1) + yB(i, j, k - 2)) - 0.25 * (xB(i, j, k - 1) + xB(i, j, k - 2)) * (yeta(i, j, k - 1) + yeta(i, j, k - 2));
                RDouble sz_CL1 = 0.25 * (xeta(i, j, k   ) + xeta(i, j, k - 1)) * (yB(i, j, k   ) + yB(i, j, k - 1)) - 0.25 * (xB(i, j, k   ) + xB(i, j, k - 1)) * (yeta(i, j, k   ) + yeta(i, j, k - 1));
                RDouble sz_fL  =         xeta(i, j, k   )                      *  yB(i, j, k   )                    -         xB(i, j, k   )                    *  yeta(i, j, k   )                     ;
                RDouble sz_fR  =         xeta(i, j, k + 1)                      *  yB(i, j, k + 1)                    -         xB(i, j, k + 1)                    *  yeta(i, j, k + 1)                     ;
                RDouble sz_CR1 = 0.25 * (xeta(i, j, k + 1) + xeta(i, j, k + 2)) * (yB(i, j, k + 1) + yB(i, j, k + 2)) - 0.25 * (xB(i, j, k + 1) + xB(i, j, k + 2)) * (yeta(i, j, k + 1) + yeta(i, j, k + 2));
                RDouble sz_CR2 = 0.25 * (xeta(i, j, k + 2) + xeta(i, j, k + 3)) * (yB(i, j, k + 2) + yB(i, j, k + 3)) - 0.25 * (xB(i, j, k + 2) + xB(i, j, k + 3)) * (yeta(i, j, k + 2) + yeta(i, j, k + 3));
                
                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));

                sx_CL2 = 0.25 * (zcta(i, j - 1, k) + zcta(i, j - 2, k)) * (yC(i, j - 1, k) + yC(i, j - 2, k)) - 0.25 * (zC(i, j - 1, k) + zC(i, j - 2, k)) * (ycta(i, j - 1, k) + ycta(i, j - 2, k));
                sx_CL1 = 0.25 * (zcta(i, j    , k) + zcta(i, j - 1, k)) * (yC(i, j    , k) + yC(i, j - 1, k)) - 0.25 * (zC(i, j    , k) + zC(i, j - 1, k)) * (ycta(i, j    , k) + ycta(i, j - 1, k));
                sx_fL  =         zcta(i, j    , k)                      *  yC(i, j    , k)                    -         zC(i, j    , k)                    *  ycta(i, j    , k)                     ;
                sx_fR  =         zcta(i, j + 1, k)                      *  yC(i, j + 1, k)                    -         zC(i, j + 1, k)                    *  ycta(i, j + 1, k)                     ;
                sx_CR1 = 0.25 * (zcta(i, j + 1, k) + zcta(i, j + 2, k)) * (yC(i, j + 1, k) + yC(i, j + 2, k)) - 0.25 * (zC(i, j + 1, k) + zC(i, j + 2, k)) * (ycta(i, j + 1, k) + ycta(i, j + 2, k));
                sx_CR2 = 0.25 * (zcta(i, j + 2, k) + zcta(i, j + 3, k)) * (yC(i, j + 2, k) + yC(i, j + 3, k)) - 0.25 * (zC(i, j + 2, k) + zC(i, j + 3, k)) * (ycta(i, j + 2, k) + ycta(i, j + 3, k));

                sy_CL2 = 0.25 * (xcta(i, j - 1, k) + xcta(i, j - 2, k)) * (zC(i, j - 1, k) + zC(i, j - 2, k)) - 0.25 * (xC(i, j - 1, k) + xC(i, j - 2, k)) * (zcta(i, j - 1, k) + zcta(i, j - 2, k));
                sy_CL1 = 0.25 * (xcta(i, j    , k) + xcta(i, j - 1, k)) * (zC(i, j    , k) + zC(i, j - 1, k)) - 0.25 * (xC(i, j    , k) + xC(i, j - 1, k)) * (zcta(i, j    , k) + zcta(i, j - 1, k));
                sy_fL  =         xcta(i, j    , k)                      *  zC(i, j    , k)                    -         xC(i, j    , k)                    *  zcta(i, j    , k)                     ;
                sy_fR  =         xcta(i, j + 1, k)                      *  zC(i, j + 1, k)                    -         xC(i, j + 1, k)                    *  zcta(i, j + 1, k)                     ;
                sy_CR1 = 0.25 * (xcta(i, j + 1, k) + xcta(i, j + 2, k)) * (zC(i, j + 1, k) + zC(i, j + 2, k)) - 0.25 * (xC(i, j + 1, k) + xC(i, j + 2, k)) * (zcta(i, j + 1, k) + zcta(i, j + 2, k));
                sy_CR2 = 0.25 * (xcta(i, j + 2, k) + xcta(i, j + 3, k)) * (zC(i, j + 2, k) + zC(i, j + 3, k)) - 0.25 * (xC(i, j + 2, k) + xC(i, j + 3, k)) * (zcta(i, j + 2, k) + zcta(i, j + 3, k));

                sz_CL2 = 0.25 * (ycta(i, j - 1, k) + ycta(i, j - 2, k)) * (xC(i, j - 1, k) + xC(i, j - 2, k)) - 0.25 * (yC(i, j - 1, k) + yC(i, j - 2, k)) * (xcta(i, j - 1, k) + xcta(i, j - 2, k));
                sz_CL1 = 0.25 * (ycta(i, j    , k) + ycta(i, j - 1, k)) * (xC(i, j    , k) + xC(i, j - 1, k)) - 0.25 * (yC(i, j    , k) + yC(i, j - 1, k)) * (xcta(i, j    , k) + xcta(i, j - 1, k));
                sz_fL  =         ycta(i, j    , k)                      *  xC(i, j    , k)                    -         yC(i, j    , k)                    *  xcta(i, j    , k)                     ;
                sz_fR  =         ycta(i, j + 1, k)                      *  xC(i, j + 1, k)                    -         yC(i, j + 1, k)                    *  xcta(i, j + 1, k)                     ;
                sz_CR1 = 0.25 * (ycta(i, j + 1, k) + ycta(i, j + 2, k)) * (xC(i, j + 1, k) + xC(i, j + 2, k)) - 0.25 * (yC(i, j + 1, k) + yC(i, j + 2, k)) * (xcta(i, j + 1, k) + xcta(i, j + 2, k));
                sz_CR2 = 0.25 * (ycta(i, j + 2, k) + ycta(i, j + 3, k)) * (xC(i, j + 2, k) + xC(i, j + 3, k)) - 0.25 * (yC(i, j + 2, k) + yC(i, j + 3, k)) * (xcta(i, j + 2, k) + xcta(i, j + 3, k));
                
                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));

                RDouble ds = sqrt(KX * KX + KY * KY + KZ * KZ);

                xfv (i, j, k, 1) = KX;
                yfv (i, j, k, 1) = KY;
                zfv (i, j, k, 1) = KZ;
                area(i, j, k, 1) = ds;
                xfn (i, j, k, 1) = KX / (ds + TINY);
                yfn (i, j, k, 1) = KY / (ds + TINY);
                zfn (i, j, k, 1) = KZ / (ds + TINY);
            }
        }
    }

    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int j = -1; j <= nj + 2; ++ j)
            {
                RDouble KX = 0.0;
                RDouble KY = 0.0;
                RDouble KZ = 0.0;

                RDouble sx_CL2 = 0.25 * (ycta(i - 1, j, k) + ycta(i - 2, j, k)) * (zC(i - 1, j, k) + zC(i - 2, j, k)) - 0.25 * (yC(i - 1, j, k) + yC(i - 2, j, k)) * (zcta(i - 1, j, k) + zcta(i - 2, j, k));
                RDouble sx_CL1 = 0.25 * (ycta(i    , j, k) + ycta(i - 1, j, k)) * (zC(i    , j, k) + zC(i - 1, j, k)) - 0.25 * (yC(i    , j, k) + yC(i - 1, j, k)) * (zcta(i    , j, k) + zcta(i - 1, j, k));
                RDouble sx_fL  =         ycta(i    , j, k)                      *  zC(i    , j, k)                    -         yC(i    , j, k)                    *  zcta(i    , j, k)                     ;
                RDouble sx_fR  =         ycta(i + 1, j, k)                      *  zC(i + 1, j, k)                    -         yC(i + 1, j, k)                    *  zcta(i + 1, j, k)                     ;
                RDouble sx_CR1 = 0.25 * (ycta(i + 1, j, k) + ycta(i + 2, j, k)) * (zC(i + 1, j, k) + zC(i + 2, j, k)) - 0.25 * (yC(i + 1, j, k) + yC(i + 2, j, k)) * (zcta(i + 1, j, k) + zcta(i + 2, j, k));
                RDouble sx_CR2 = 0.25 * (ycta(i + 2, j, k) + ycta(i + 3, j, k)) * (zC(i + 2, j, k) + zC(i + 3, j, k)) - 0.25 * (yC(i + 2, j, k) + yC(i + 3, j, k)) * (zcta(i + 2, j, k) + zcta(i + 3, j, k));

                RDouble sy_CL2 = 0.25 * (zcta(i - 1, j, k) + zcta(i - 2, j, k)) * (xC(i - 1, j, k) + xC(i - 2, j, k)) - 0.25 * (zC(i - 1, j, k) + zC(i - 2, j, k)) * (xcta(i - 1, j, k) + xcta(i - 2, j, k));
                RDouble sy_CL1 = 0.25 * (zcta(i    , j, k) + zcta(i - 1, j, k)) * (xC(i    , j, k) + xC(i - 1, j, k)) - 0.25 * (zC(i    , j, k) + zC(i - 1, j, k)) * (xcta(i    , j, k) + xcta(i - 1, j, k));
                RDouble sy_fL  =         zcta(i    , j, k)                      *  xC(i    , j, k)                    -         zC(i    , j, k)                    *  xcta(i    , j, k)                     ;
                RDouble sy_fR  =         zcta(i + 1, j, k)                      *  xC(i + 1, j, k)                    -         zC(i + 1, j, k)                    *  xcta(i + 1, j, k)                     ;
                RDouble sy_CR1 = 0.25 * (zcta(i + 1, j, k) + zcta(i + 2, j, k)) * (xC(i + 1, j, k) + xC(i + 2, j, k)) - 0.25 * (zC(i + 1, j, k) + zC(i + 2, j, k)) * (xcta(i + 1, j, k) + xcta(i + 2, j, k));
                RDouble sy_CR2 = 0.25 * (zcta(i + 2, j, k) + zcta(i + 3, j, k)) * (xC(i + 2, j, k) + xC(i + 3, j, k)) - 0.25 * (zC(i + 2, j, k) + zC(i + 3, j, k)) * (xcta(i + 2, j, k) + xcta(i + 3, j, k));

                RDouble sz_CL2 = 0.25 * (xcta(i - 1, j, k) + xcta(i - 2, j, k)) * (yC(i - 1, j, k) + yC(i - 2, j, k)) - 0.25 * (xC(i - 1, j, k) + xC(i - 2, j, k)) * (ycta(i - 1, j, k) + ycta(i - 2, j, k));
                RDouble sz_CL1 = 0.25 * (xcta(i    , j, k) + xcta(i - 1, j, k)) * (yC(i    , j, k) + yC(i - 1, j, k)) - 0.25 * (xC(i    , j, k) + xC(i - 1, j, k)) * (ycta(i    , j, k) + ycta(i - 1, j, k));
                RDouble sz_fL  =         xcta(i    , j, k)                      *  yC(i    , j, k)                    -         xC(i    , j, k)                    *  ycta(i    , j, k)                     ;
                RDouble sz_fR  =         xcta(i + 1, j, k)                      *  yC(i + 1, j, k)                    -         xC(i + 1, j, k)                    *  ycta(i + 1, j, k)                     ;
                RDouble sz_CR1 = 0.25 * (xcta(i + 1, j, k) + xcta(i + 2, j, k)) * (yC(i + 1, j, k) + yC(i + 2, j, k)) - 0.25 * (xC(i + 1, j, k) + xC(i + 2, j, k)) * (ycta(i + 1, j, k) + ycta(i + 2, j, k));
                RDouble sz_CR2 = 0.25 * (xcta(i + 2, j, k) + xcta(i + 3, j, k)) * (yC(i + 2, j, k) + yC(i + 3, j, k)) - 0.25 * (xC(i + 2, j, k) + xC(i + 3, j, k)) * (ycta(i + 2, j, k) + ycta(i + 3, j, k));
                
                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));

                sx_CL2 = 0.25 * (zkxi(i, j, k - 1) + zkxi(i, j, k - 2)) * (yA(i, j, k - 1) + yA(i, j, k - 2)) - 0.25 * (zA(i, j, k - 1) + zA(i, j, k - 2)) * (ykxi(i, j, k - 1) + ykxi(i, j, k - 2));
                sx_CL1 = 0.25 * (zkxi(i, j, k   ) + zkxi(i, j, k - 1)) * (yA(i, j, k   ) + yA(i, j, k - 1)) - 0.25 * (zA(i, j, k   ) + zA(i, j, k - 1)) * (ykxi(i, j, k   ) + ykxi(i, j, k - 1));
                sx_fL  =         zkxi(i, j, k   )                      *  yA(i, j, k   )                    -         zA(i, j, k   )                    *  ykxi(i, j, k   )                     ;
                sx_fR  =         zkxi(i, j, k + 1)                      *  yA(i, j, k + 1)                    -         zA(i, j, k + 1)                    *  ykxi(i, j, k + 1)                     ;
                sx_CR1 = 0.25 * (zkxi(i, j, k + 1) + zkxi(i, j, k + 2)) * (yA(i, j, k + 1) + yA(i, j, k + 2)) - 0.25 * (zA(i, j, k + 1) + zA(i, j, k + 2)) * (ykxi(i, j, k + 1) + ykxi(i, j, k + 2));
                sx_CR2 = 0.25 * (zkxi(i, j, k + 2) + zkxi(i, j, k + 3)) * (yA(i, j, k + 2) + yA(i, j, k + 3)) - 0.25 * (zA(i, j, k + 2) + zA(i, j, k + 3)) * (ykxi(i, j, k + 2) + ykxi(i, j, k + 3));

                sy_CL2 = 0.25 * (xkxi(i, j, k - 1) + xkxi(i, j, k - 2)) * (zA(i, j, k - 1) + zA(i, j, k - 2)) - 0.25 * (xA(i, j, k - 1) + xA(i, j, k - 2)) * (zkxi(i, j, k - 1) + zkxi(i, j, k - 2));
                sy_CL1 = 0.25 * (xkxi(i, j, k   ) + xkxi(i, j, k - 1)) * (zA(i, j, k   ) + zA(i, j, k - 1)) - 0.25 * (xA(i, j, k   ) + xA(i, j, k - 1)) * (zkxi(i, j, k   ) + zkxi(i, j, k - 1));
                sy_fL  =         xkxi(i, j, k   )                      *  zA(i, j, k   )                    -         xA(i, j, k   )                    *  zkxi(i, j, k   )                     ;
                sy_fR  =         xkxi(i, j, k + 1)                      *  zA(i, j, k + 1)                    -         xA(i, j, k + 1)                    *  zkxi(i, j, k + 1)                     ;
                sy_CR1 = 0.25 * (xkxi(i, j, k + 1) + xkxi(i, j, k + 2)) * (zA(i, j, k + 1) + zA(i, j, k + 2)) - 0.25 * (xA(i, j, k + 1) + xA(i, j, k + 2)) * (zkxi(i, j, k + 1) + zkxi(i, j, k + 2));
                sy_CR2 = 0.25 * (xkxi(i, j, k + 2) + xkxi(i, j, k + 3)) * (zA(i, j, k + 2) + zA(i, j, k + 3)) - 0.25 * (xA(i, j, k + 2) + xA(i, j, k + 3)) * (zkxi(i, j, k + 2) + zkxi(i, j, k + 3));

                sz_CL2 = 0.25 * (ykxi(i, j, k - 1) + ykxi(i, j, k - 2)) * (xA(i, j, k - 1) + xA(i, j, k - 2)) - 0.25 * (yA(i, j, k - 1) + yA(i, j, k - 2)) * (xkxi(i, j, k - 1) + xkxi(i, j, k - 2));
                sz_CL1 = 0.25 * (ykxi(i, j, k   ) + ykxi(i, j, k - 1)) * (xA(i, j, k   ) + xA(i, j, k - 1)) - 0.25 * (yA(i, j, k   ) + yA(i, j, k - 1)) * (xkxi(i, j, k   ) + xkxi(i, j, k - 1));
                sz_fL  =         ykxi(i, j, k   )                      *  xA(i, j, k   )                    -         yA(i, j, k   )                    *  xkxi(i, j, k   )                     ;
                sz_fR  =         ykxi(i, j, k + 1)                      *  xA(i, j, k + 1)                    -         yA(i, j, k + 1)                    *  xkxi(i, j, k + 1)                     ;
                sz_CR1 = 0.25 * (ykxi(i, j, k + 1) + ykxi(i, j, k + 2)) * (xA(i, j, k + 1) + xA(i, j, k + 2)) - 0.25 * (yA(i, j, k + 1) + yA(i, j, k + 2)) * (xkxi(i, j, k + 1) + xkxi(i, j, k + 2));
                sz_CR2 = 0.25 * (ykxi(i, j, k + 2) + ykxi(i, j, k + 3)) * (xA(i, j, k + 2) + xA(i, j, k + 3)) - 0.25 * (yA(i, j, k + 2) + yA(i, j, k + 3)) * (xkxi(i, j, k + 2) + xkxi(i, j, k + 3));
                
                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));

                RDouble ds = sqrt(KX * KX + KY * KY + KZ * KZ);

                xfv (i, j, k, 2) = KX;
                yfv (i, j, k, 2) = KY;
                zfv (i, j, k, 2) = KZ;
                area(i, j, k, 2) = ds;
                xfn (i, j, k, 2) = KX / (ds + TINY);
                yfn (i, j, k, 2) = KY / (ds + TINY);
                zfn (i, j, k, 2) = KZ / (ds + TINY);
            }
        }
    }

    for (int j = 1; j <= nj - 1; ++ j)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int k = -1; k <= nk + 2; ++ k)
            {
                RDouble KX = 0.0;
                RDouble KY = 0.0;
                RDouble KZ = 0.0;

                RDouble sx_CL2 = 0.25 * (ykxi(i, j - 1, k) + ykxi(i, j - 2, k)) * (zA(i, j - 1, k) + zA(i, j - 2, k)) - 0.25 * (yA(i, j - 1, k) + yA(i, j - 2, k)) * (zkxi(i, j - 1, k) + zkxi(i, j - 2, k));
                RDouble sx_CL1 = 0.25 * (ykxi(i, j    , k) + ykxi(i, j - 1, k)) * (zA(i, j    , k) + zA(i, j - 1, k)) - 0.25 * (yA(i, j    , k) + yA(i, j - 1, k)) * (zkxi(i, j    , k) + zkxi(i, j - 1, k));
                RDouble sx_fL  =         ykxi(i, j    , k)                      *  zA(i, j    , k)                    -         yA(i, j    , k)                    *  zkxi(i, j    , k)                     ;
                RDouble sx_fR  =         ykxi(i, j + 1, k)                      *  zA(i, j + 1, k)                    -         yA(i, j + 1, k)                    *  zkxi(i, j + 1, k)                     ;
                RDouble sx_CR1 = 0.25 * (ykxi(i, j + 1, k) + ykxi(i, j + 2, k)) * (zA(i, j + 1, k) + zA(i, j + 2, k)) - 0.25 * (yA(i, j + 1, k) + yA(i, j + 2, k)) * (zkxi(i, j + 1, k) + zkxi(i, j + 2, k));
                RDouble sx_CR2 = 0.25 * (ykxi(i, j + 2, k) + ykxi(i, j + 3, k)) * (zA(i, j + 2, k) + zA(i, j + 3, k)) - 0.25 * (yA(i, j + 2, k) + yA(i, j + 3, k)) * (zkxi(i, j + 2, k) + zkxi(i, j + 3, k));

                RDouble sy_CL2 = 0.25 * (zkxi(i, j - 1, k) + zkxi(i, j - 2, k)) * (xA(i, j - 1, k) + xA(i, j - 2, k)) - 0.25 * (zA(i, j - 1, k) + zA(i, j - 2, k)) * (xkxi(i, j - 1, k) + xkxi(i, j - 2, k));
                RDouble sy_CL1 = 0.25 * (zkxi(i, j    , k) + zkxi(i, j - 1, k)) * (xA(i, j    , k) + xA(i, j - 1, k)) - 0.25 * (zA(i, j    , k) + zA(i, j - 1, k)) * (xkxi(i, j    , k) + xkxi(i, j - 1, k));
                RDouble sy_fL  =         zkxi(i, j    , k)                      *  xA(i, j    , k)                    -         zA(i, j    , k)                    *  xkxi(i, j    , k)                     ;
                RDouble sy_fR  =         zkxi(i, j + 1, k)                      *  xA(i, j + 1, k)                    -         zA(i, j + 1, k)                    *  xkxi(i, j + 1, k)                     ;
                RDouble sy_CR1 = 0.25 * (zkxi(i, j + 1, k) + zkxi(i, j + 2, k)) * (xA(i, j + 1, k) + xA(i, j + 2, k)) - 0.25 * (zA(i, j + 1, k) + zA(i, j + 2, k)) * (xkxi(i, j + 1, k) + xkxi(i, j + 2, k));
                RDouble sy_CR2 = 0.25 * (zkxi(i, j + 2, k) + zkxi(i, j + 3, k)) * (xA(i, j + 2, k) + xA(i, j + 3, k)) - 0.25 * (zA(i, j + 2, k) + zA(i, j + 3, k)) * (xkxi(i, j + 2, k) + xkxi(i, j + 3, k));

                RDouble sz_CL2 = 0.25 * (xkxi(i, j - 1, k) + xkxi(i, j - 2, k)) * (yA(i, j - 1, k) + yA(i, j - 2, k)) - 0.25 * (xA(i, j - 1, k) + xA(i, j - 2, k)) * (ykxi(i, j - 1, k) + ykxi(i, j - 2, k));
                RDouble sz_CL1 = 0.25 * (xkxi(i, j    , k) + xkxi(i, j - 1, k)) * (yA(i, j    , k) + yA(i, j - 1, k)) - 0.25 * (xA(i, j    , k) + xA(i, j - 1, k)) * (ykxi(i, j    , k) + ykxi(i, j - 1, k));
                RDouble sz_fL  =         xkxi(i, j    , k)                      *  yA(i, j    , k)                    -         xA(i, j    , k)                    *  ykxi(i, j    , k)                     ;
                RDouble sz_fR  =         xkxi(i, j + 1, k)                      *  yA(i, j + 1, k)                    -         xA(i, j + 1, k)                    *  ykxi(i, j + 1, k)                     ;
                RDouble sz_CR1 = 0.25 * (xkxi(i, j + 1, k) + xkxi(i, j + 2, k)) * (yA(i, j + 1, k) + yA(i, j + 2, k)) - 0.25 * (xA(i, j + 1, k) + xA(i, j + 2, k)) * (ykxi(i, j + 1, k) + ykxi(i, j + 2, k));
                RDouble sz_CR2 = 0.25 * (xkxi(i, j + 2, k) + xkxi(i, j + 3, k)) * (yA(i, j + 2, k) + yA(i, j + 3, k)) - 0.25 * (xA(i, j + 2, k) + xA(i, j + 3, k)) * (ykxi(i, j + 2, k) + ykxi(i, j + 3, k));
                
                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));

                sx_CL2 = 0.25 * (zeta(i - 1, j, k) + zeta(i - 2, j, k)) * (yB(i - 1, j, k) + yB(i - 2, j, k)) - 0.25 * (zB(i - 1, j, k) + zB(i - 2, j, k)) * (yeta(i - 1, j, k) + yeta(i - 2, j, k));
                sx_CL1 = 0.25 * (zeta(i    , j, k) + zeta(i - 1, j, k)) * (yB(i    , j, k) + yB(i - 1, j, k)) - 0.25 * (zB(i    , j, k) + zB(i - 1, j, k)) * (yeta(i    , j, k) + yeta(i - 1, j, k));
                sx_fL  =         zeta(i    , j, k)                      *  yB(i    , j, k)                    -         zB(i    , j, k)                    *  yeta(i    , j, k)                     ;
                sx_fR  =         zeta(i + 1, j, k)                      *  yB(i + 1, j, k)                    -         zB(i + 1, j, k)                    *  yeta(i + 1, j, k)                     ;
                sx_CR1 = 0.25 * (zeta(i + 1, j, k) + zeta(i + 2, j, k)) * (yB(i + 1, j, k) + yB(i + 2, j, k)) - 0.25 * (zB(i + 1, j, k) + zB(i + 2, j, k)) * (yeta(i + 1, j, k) + yeta(i + 2, j, k));
                sx_CR2 = 0.25 * (zeta(i + 2, j, k) + zeta(i + 3, j, k)) * (yB(i + 2, j, k) + yB(i + 3, j, k)) - 0.25 * (zB(i + 2, j, k) + zB(i + 3, j, k)) * (yeta(i + 2, j, k) + yeta(i + 3, j, k));

                sy_CL2 = 0.25 * (xeta(i - 1, j, k) + xeta(i - 2, j, k)) * (zB(i - 1, j, k) + zB(i - 2, j, k)) - 0.25 * (xB(i - 1, j, k) + xB(i - 2, j, k)) * (zeta(i - 1, j, k) + zeta(i - 2, j, k));
                sy_CL1 = 0.25 * (xeta(i    , j, k) + xeta(i - 1, j, k)) * (zB(i    , j, k) + zB(i - 1, j, k)) - 0.25 * (xB(i    , j, k) + xB(i - 1, j, k)) * (zeta(i    , j, k) + zeta(i - 1, j, k));
                sy_fL  =         xeta(i    , j, k)                      *  zB(i    , j, k)                    -         xB(i    , j, k)                    *  zeta(i    , j, k)                     ;
                sy_fR  =         xeta(i + 1, j, k)                      *  zB(i + 1, j, k)                    -         xB(i + 1, j, k)                    *  zeta(i + 1, j, k)                     ;
                sy_CR1 = 0.25 * (xeta(i + 1, j, k) + xeta(i + 2, j, k)) * (zB(i + 1, j, k) + zB(i + 2, j, k)) - 0.25 * (xB(i + 1, j, k) + xB(i + 2, j, k)) * (zeta(i + 1, j, k) + zeta(i + 2, j, k));
                sy_CR2 = 0.25 * (xeta(i + 2, j, k) + xeta(i + 3, j, k)) * (zB(i + 2, j, k) + zB(i + 3, j, k)) - 0.25 * (xB(i + 2, j, k) + xB(i + 3, j, k)) * (zeta(i + 2, j, k) + zeta(i + 3, j, k));

                sz_CL2 = 0.25 * (yeta(i - 1, j, k) + yeta(i - 2, j, k)) * (xB(i - 1, j, k) + xB(i - 2, j, k)) - 0.25 * (yB(i - 1, j, k) + yB(i - 2, j, k)) * (xeta(i - 1, j, k) + xeta(i - 2, j, k));
                sz_CL1 = 0.25 * (yeta(i    , j, k) + yeta(i - 1, j, k)) * (xB(i    , j, k) + xB(i - 1, j, k)) - 0.25 * (yB(i    , j, k) + yB(i - 1, j, k)) * (xeta(i    , j, k) + xeta(i - 1, j, k));
                sz_fL  =         yeta(i    , j, k)                      *  xB(i    , j, k)                    -         yB(i    , j, k)                    *  xeta(i    , j, k)                     ;
                sz_fR  =         yeta(i + 1, j, k)                      *  xB(i + 1, j, k)                    -         yB(i + 1, j, k)                    *  xeta(i + 1, j, k)                     ;
                sz_CR1 = 0.25 * (yeta(i + 1, j, k) + yeta(i + 2, j, k)) * (xB(i + 1, j, k) + xB(i + 2, j, k)) - 0.25 * (yB(i + 1, j, k) + yB(i + 2, j, k)) * (xeta(i + 1, j, k) + xeta(i + 2, j, k));
                sz_CR2 = 0.25 * (yeta(i + 2, j, k) + yeta(i + 3, j, k)) * (xB(i + 2, j, k) + xB(i + 3, j, k)) - 0.25 * (yB(i + 2, j, k) + yB(i + 3, j, k)) * (xeta(i + 2, j, k) + xeta(i + 3, j, k));
                
                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));

                RDouble ds = sqrt(KX * KX + KY * KY + KZ * KZ);

                xfv (i, j, k, 3) = KX;
                yfv (i, j, k, 3) = KY;
                zfv (i, j, k, 3) = KZ;
                area(i, j, k, 3) = ds;
                xfn (i, j, k, 3) = KX / (ds + TINY);
                yfn (i, j, k, 3) = KY / (ds + TINY);
                zfn (i, j, k, 3) = KZ / (ds + TINY);
            }
        }
    }
    //!=
    RDouble4D &Cellxfv = *new RDouble4D(Range(-1, ni + 1), Range(-1, nj + 1), Range(-1, nk + 1), Range(1, 3), fortranArray);
    RDouble4D &Cellyfv = *new RDouble4D(Range(-1, ni + 1), Range(-1, nj + 1), Range(-1, nk + 1), Range(1, 3), fortranArray);
    RDouble4D &Cellzfv = *new RDouble4D(Range(-1, ni + 1), Range(-1, nj + 1), Range(-1, nk + 1), Range(1, 3), fortranArray);

    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = -1; i <= ni + 1; ++ i)
            {
                RDouble KX = 0.0;
                RDouble KY = 0.0;
                RDouble KZ = 0.0;

                RDouble sx_CL2 = 0.0625 * (yeta(i, j, k - 1) + yeta(i, j, k - 2) + yeta(i + 1, j, k - 1) + yeta(i + 1, j, k - 2)) * (zB(i, j, k - 1) + zB(i, j, k - 2) + zB(i + 1, j, k - 1) + zB(i + 1, j, k - 2)) - 0.0625 * (yB(i, j, k - 1) + yB(i, j, k - 2) + yB(i + 1, j, k - 1) + yB(i + 1, j, k - 2)) * (zeta(i, j, k - 1) + zeta(i, j, k - 2) + zeta(i + 1, j, k - 1) + zeta(i + 1, j, k - 2));
                RDouble sx_CL1 = 0.0625 * (yeta(i, j, k   ) + yeta(i, j, k - 1) + yeta(i + 1, j, k   ) + yeta(i + 1, j, k - 1)) * (zB(i, j, k   ) + zB(i, j, k - 1) + zB(i + 1, j, k   ) + zB(i + 1, j, k - 1)) - 0.0625 * (yB(i, j, k   ) + yB(i, j, k - 1) + yB(i + 1, j, k   ) + yB(i + 1, j, k - 1)) * (zeta(i, j, k   ) + zeta(i, j, k - 1) + zeta(i + 1, j, k   ) + zeta(i + 1, j, k - 1));
                RDouble sx_fL  = 0.0625 * (yeta(i, j, k   ) + yeta(i, j, k   ) + yeta(i + 1, j, k   ) + yeta(i + 1, j, k   )) * (zB(i, j, k   ) + zB(i, j, k   ) + zB(i + 1, j, k   ) + zB(i + 1, j, k   )) - 0.0625 * (yB(i, j, k   ) + yB(i, j, k   ) + yB(i + 1, j, k   ) + yB(i + 1, j, k   )) * (zeta(i, j, k   ) + zeta(i, j, k   ) + zeta(i + 1, j, k   ) + zeta(i + 1, j, k   ));
                RDouble sx_fR  = 0.0625 * (yeta(i, j, k + 1) + yeta(i, j, k + 1) + yeta(i + 1, j, k + 1) + yeta(i + 1, j, k + 1)) * (zB(i, j, k + 1) + zB(i, j, k + 1) + zB(i + 1, j, k + 1) + zB(i + 1, j, k + 1)) - 0.0625 * (yB(i, j, k + 1) + yB(i, j, k + 1) + yB(i + 1, j, k + 1) + yB(i + 1, j, k + 1)) * (zeta(i, j, k + 1) + zeta(i, j, k + 1) + zeta(i + 1, j, k + 1) + zeta(i + 1, j, k + 1));
                RDouble sx_CR1 = 0.0625 * (yeta(i, j, k + 1) + yeta(i, j, k + 2) + yeta(i + 1, j, k + 1) + yeta(i + 1, j, k + 2)) * (zB(i, j, k + 1) + zB(i, j, k + 2) + zB(i + 1, j, k + 1) + zB(i + 1, j, k + 2)) - 0.0625 * (yB(i, j, k + 1) + yB(i, j, k + 2) + yB(i + 1, j, k + 1) + yB(i + 1, j, k + 2)) * (zeta(i, j, k + 1) + zeta(i, j, k + 2) + zeta(i + 1, j, k + 1) + zeta(i + 1, j, k + 2));
                RDouble sx_CR2 = 0.0625 * (yeta(i, j, k + 2) + yeta(i, j, k + 3) + yeta(i + 1, j, k + 2) + yeta(i + 1, j, k + 3)) * (zB(i, j, k + 2) + zB(i, j, k + 3) + zB(i + 1, j, k + 2) + zB(i + 1, j, k + 3)) - 0.0625 * (yB(i, j, k + 2) + yB(i, j, k + 3) + yB(i + 1, j, k + 2) + yB(i + 1, j, k + 3)) * (zeta(i, j, k + 2) + zeta(i, j, k + 3) + zeta(i + 1, j, k + 2) + zeta(i + 1, j, k + 3));
                RDouble sy_CL2 = 0.0625 * (zeta(i, j, k - 1) + zeta(i, j, k - 2) + zeta(i + 1, j, k - 1) + zeta(i + 1, j, k - 2)) * (xB(i, j, k - 1) + xB(i, j, k - 2) + xB(i + 1, j, k - 1) + xB(i + 1, j, k - 2)) - 0.0625 * (zB(i, j, k - 1) + zB(i, j, k - 2) + zB(i + 1, j, k - 1) + zB(i + 1, j, k - 2)) * (xeta(i, j, k - 1) + xeta(i, j, k - 2) + xeta(i + 1, j, k - 1) + xeta(i + 1, j, k - 2));
                RDouble sy_CL1 = 0.0625 * (zeta(i, j, k   ) + zeta(i, j, k - 1) + zeta(i + 1, j, k   ) + zeta(i + 1, j, k - 1)) * (xB(i, j, k   ) + xB(i, j, k - 1) + xB(i + 1, j, k   ) + xB(i + 1, j, k - 1)) - 0.0625 * (zB(i, j, k   ) + zB(i, j, k - 1) + zB(i + 1, j, k   ) + zB(i + 1, j, k - 1)) * (xeta(i, j, k   ) + xeta(i, j, k - 1) + xeta(i + 1, j, k   ) + xeta(i + 1, j, k - 1));
                RDouble sy_fL  = 0.0625 * (zeta(i, j, k   ) + zeta(i, j, k   ) + zeta(i + 1, j, k   ) + zeta(i + 1, j, k   )) * (xB(i, j, k   ) + xB(i, j, k   ) + xB(i + 1, j, k   ) + xB(i + 1, j, k   )) - 0.0625 * (zB(i, j, k   ) + zB(i, j, k   ) + zB(i + 1, j, k   ) + zB(i + 1, j, k   )) * (xeta(i, j, k   ) + xeta(i, j, k   ) + xeta(i + 1, j, k   ) + xeta(i + 1, j, k   ));
                RDouble sy_fR  = 0.0625 * (zeta(i, j, k + 1) + zeta(i, j, k + 1) + zeta(i + 1, j, k + 1) + zeta(i + 1, j, k + 1)) * (xB(i, j, k + 1) + xB(i, j, k + 1) + xB(i + 1, j, k + 1) + xB(i + 1, j, k + 1)) - 0.0625 * (zB(i, j, k + 1) + zB(i, j, k + 1) + zB(i + 1, j, k + 1) + zB(i + 1, j, k + 1)) * (xeta(i, j, k + 1) + xeta(i, j, k + 1) + xeta(i + 1, j, k + 1) + xeta(i + 1, j, k + 1));
                RDouble sy_CR1 = 0.0625 * (zeta(i, j, k + 1) + zeta(i, j, k + 2) + zeta(i + 1, j, k + 1) + zeta(i + 1, j, k + 2)) * (xB(i, j, k + 1) + xB(i, j, k + 2) + xB(i + 1, j, k + 1) + xB(i + 1, j, k + 2)) - 0.0625 * (zB(i, j, k + 1) + zB(i, j, k + 2) + zB(i + 1, j, k + 1) + zB(i + 1, j, k + 2)) * (xeta(i, j, k + 1) + xeta(i, j, k + 2) + xeta(i + 1, j, k + 1) + xeta(i + 1, j, k + 2));
                RDouble sy_CR2 = 0.0625 * (zeta(i, j, k + 2) + zeta(i, j, k + 3) + zeta(i + 1, j, k + 2) + zeta(i + 1, j, k + 3)) * (xB(i, j, k + 2) + xB(i, j, k + 3) + xB(i + 1, j, k + 2) + xB(i + 1, j, k + 3)) - 0.0625 * (zB(i, j, k + 2) + zB(i, j, k + 3) + zB(i + 1, j, k + 2) + zB(i + 1, j, k + 3)) * (xeta(i, j, k + 2) + xeta(i, j, k + 3) + xeta(i + 1, j, k + 2) + xeta(i + 1, j, k + 3));
                RDouble sz_CL2 = 0.0625 * (xeta(i, j, k - 1) + xeta(i, j, k - 2) + xeta(i + 1, j, k - 1) + xeta(i + 1, j, k - 2)) * (yB(i, j, k - 1) + yB(i, j, k - 2) + yB(i + 1, j, k - 1) + yB(i + 1, j, k - 2)) - 0.0625 * (xB(i, j, k - 1) + xB(i, j, k - 2) + xB(i + 1, j, k - 1) + xB(i + 1, j, k - 2)) * (yeta(i, j, k - 1) + yeta(i, j, k - 2) + yeta(i + 1, j, k - 1) + yeta(i + 1, j, k - 2));
                RDouble sz_CL1 = 0.0625 * (xeta(i, j, k   ) + xeta(i, j, k - 1) + xeta(i + 1, j, k   ) + xeta(i + 1, j, k - 1)) * (yB(i, j, k   ) + yB(i, j, k - 1) + yB(i + 1, j, k   ) + yB(i + 1, j, k - 1)) - 0.0625 * (xB(i, j, k   ) + xB(i, j, k - 1) + xB(i + 1, j, k   ) + xB(i + 1, j, k - 1)) * (yeta(i, j, k   ) + yeta(i, j, k - 1) + yeta(i + 1, j, k   ) + yeta(i + 1, j, k - 1));
                RDouble sz_fL  = 0.0625 * (xeta(i, j, k   ) + xeta(i, j, k   ) + xeta(i + 1, j, k   ) + xeta(i + 1, j, k   )) * (yB(i, j, k   ) + yB(i, j, k   ) + yB(i + 1, j, k   ) + yB(i + 1, j, k   )) - 0.0625 * (xB(i, j, k   ) + xB(i, j, k   ) + xB(i + 1, j, k   ) + xB(i + 1, j, k   )) * (yeta(i, j, k   ) + yeta(i, j, k   ) + yeta(i + 1, j, k   ) + yeta(i + 1, j, k   ));
                RDouble sz_fR  = 0.0625 * (xeta(i, j, k + 1) + xeta(i, j, k + 1) + xeta(i + 1, j, k + 1) + xeta(i + 1, j, k + 1)) * (yB(i, j, k + 1) + yB(i, j, k + 1) + yB(i + 1, j, k + 1) + yB(i + 1, j, k + 1)) - 0.0625 * (xB(i, j, k + 1) + xB(i, j, k + 1) + xB(i + 1, j, k + 1) + xB(i + 1, j, k + 1)) * (yeta(i, j, k + 1) + yeta(i, j, k + 1) + yeta(i + 1, j, k + 1) + yeta(i + 1, j, k + 1));
                RDouble sz_CR1 = 0.0625 * (xeta(i, j, k + 1) + xeta(i, j, k + 2) + xeta(i + 1, j, k + 1) + xeta(i + 1, j, k + 2)) * (yB(i, j, k + 1) + yB(i, j, k + 2) + yB(i + 1, j, k + 1) + yB(i + 1, j, k + 2)) - 0.0625 * (xB(i, j, k + 1) + xB(i, j, k + 2) + xB(i + 1, j, k + 1) + xB(i + 1, j, k + 2)) * (yeta(i, j, k + 1) + yeta(i, j, k + 2) + yeta(i + 1, j, k + 1) + yeta(i + 1, j, k + 2));
                RDouble sz_CR2 = 0.0625 * (xeta(i, j, k + 2) + xeta(i, j, k + 3) + xeta(i + 1, j, k + 2) + xeta(i + 1, j, k + 3)) * (yB(i, j, k + 2) + yB(i, j, k + 3) + yB(i + 1, j, k + 2) + yB(i + 1, j, k + 3)) - 0.0625 * (xB(i, j, k + 2) + xB(i, j, k + 3) + xB(i + 1, j, k + 2) + xB(i + 1, j, k + 3)) * (yeta(i, j, k + 2) + yeta(i, j, k + 3) + yeta(i + 1, j, k + 2) + yeta(i + 1, j, k + 3));

                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));

                sx_CL2 =  0.0625 * (yC(i, j - 1, k) + yC(i, j - 2, k) + yC(i + 1, j - 1, k) + yC(i + 1, j - 2, k)) * (zcta(i, j - 1, k) + zcta(i, j - 2, k) + zcta(i + 1, j - 1, k) + zcta(i + 1, j - 2, k)) - 0.0625 * (ycta(i, j - 1, k) + ycta(i, j - 2, k) + ycta(i + 1, j - 1, k) + ycta(i + 1, j - 2, k)) * (zC(i, j - 1, k) + zC(i, j - 2, k) + zC(i + 1, j - 1, k) + zC(i + 1, j - 2, k));
                sx_CL1 =  0.0625 * (yC(i, j    , k) + yC(i, j - 1, k) + yC(i + 1, j    , k) + yC(i + 1, j - 1, k)) * (zcta(i, j    , k) + zcta(i, j - 1, k) + zcta(i + 1, j    , k) + zcta(i + 1, j - 1, k)) - 0.0625 * (ycta(i, j    , k) + ycta(i, j - 1, k) + ycta(i + 1, j    , k) + ycta(i + 1, j - 1, k)) * (zC(i, j    , k) + zC(i, j - 1, k) + zC(i + 1, j    , k) + zC(i + 1, j - 1, k));
                sx_fL  =  0.0625 * (yC(i, j    , k) + yC(i, j    , k) + yC(i + 1, j    , k) + yC(i + 1, j    , k)) * (zcta(i, j    , k) + zcta(i, j    , k) + zcta(i + 1, j    , k) + zcta(i + 1, j    , k)) - 0.0625 * (ycta(i, j    , k) + ycta(i, j    , k) + ycta(i + 1, j    , k) + ycta(i + 1, j    , k)) * (zC(i, j    , k) + zC(i, j    , k) + zC(i + 1, j    , k) + zC(i + 1, j    , k));
                sx_fR  =  0.0625 * (yC(i, j + 1, k) + yC(i, j + 1, k) + yC(i + 1, j + 1, k) + yC(i + 1, j + 1, k)) * (zcta(i, j + 1, k) + zcta(i, j + 1, k) + zcta(i + 1, j + 1, k) + zcta(i + 1, j + 1, k)) - 0.0625 * (ycta(i, j + 1, k) + ycta(i, j + 1, k) + ycta(i + 1, j + 1, k) + ycta(i + 1, j + 1, k)) * (zC(i, j + 1, k) + zC(i, j + 1, k) + zC(i + 1, j + 1, k) + zC(i + 1, j + 1, k));
                sx_CR1 =  0.0625 * (yC(i, j + 1, k) + yC(i, j + 2, k) + yC(i + 1, j + 1, k) + yC(i + 1, j + 2, k)) * (zcta(i, j + 1, k) + zcta(i, j + 2, k) + zcta(i + 1, j + 1, k) + zcta(i + 1, j + 2, k)) - 0.0625 * (ycta(i, j + 1, k) + ycta(i, j + 2, k) + ycta(i + 1, j + 1, k) + ycta(i + 1, j + 2, k)) * (zC(i, j + 1, k) + zC(i, j + 2, k) + zC(i + 1, j + 1, k) + zC(i + 1, j + 2, k));
                sx_CR2 =  0.0625 * (yC(i, j + 2, k) + yC(i, j + 3, k) + yC(i + 1, j + 2, k) + yC(i + 1, j + 3, k)) * (zcta(i, j + 2, k) + zcta(i, j + 3, k) + zcta(i + 1, j + 2, k) + zcta(i + 1, j + 3, k)) - 0.0625 * (ycta(i, j + 2, k) + ycta(i, j + 3, k) + ycta(i + 1, j + 2, k) + ycta(i + 1, j + 3, k)) * (zC(i, j + 2, k) + zC(i, j + 3, k) + zC(i + 1, j + 2, k) + zC(i + 1, j + 3, k));
                sy_CL2 =  0.0625 * (zC(i, j - 1, k) + zC(i, j - 2, k) + zC(i + 1, j - 1, k) + zC(i + 1, j - 2, k)) * (xcta(i, j - 1, k) + xcta(i, j - 2, k) + xcta(i + 1, j - 1, k) + xcta(i + 1, j - 2, k)) - 0.0625 * (zcta(i, j - 1, k) + zcta(i, j - 2, k) + zcta(i + 1, j - 1, k) + zcta(i + 1, j - 2, k)) * (xC(i, j - 1, k) + xC(i, j - 2, k) + xC(i + 1, j - 1, k) + xC(i + 1, j - 2, k));
                sy_CL1 =  0.0625 * (zC(i, j    , k) + zC(i, j - 1, k) + zC(i + 1, j    , k) + zC(i + 1, j - 1, k)) * (xcta(i, j    , k) + xcta(i, j - 1, k) + xcta(i + 1, j    , k) + xcta(i + 1, j - 1, k)) - 0.0625 * (zcta(i, j    , k) + zcta(i, j - 1, k) + zcta(i + 1, j    , k) + zcta(i + 1, j - 1, k)) * (xC(i, j    , k) + xC(i, j - 1, k) + xC(i + 1, j    , k) + xC(i + 1, j - 1, k));
                sy_fL  =  0.0625 * (zC(i, j    , k) + zC(i, j    , k) + zC(i + 1, j    , k) + zC(i + 1, j    , k)) * (xcta(i, j    , k) + xcta(i, j    , k) + xcta(i + 1, j    , k) + xcta(i + 1, j    , k)) - 0.0625 * (zcta(i, j    , k) + zcta(i, j    , k) + zcta(i + 1, j    , k) + zcta(i + 1, j    , k)) * (xC(i, j    , k) + xC(i, j    , k) + xC(i + 1, j    , k) + xC(i + 1, j    , k));
                sy_fR  =  0.0625 * (zC(i, j + 1, k) + zC(i, j + 1, k) + zC(i + 1, j + 1, k) + zC(i + 1, j + 1, k)) * (xcta(i, j + 1, k) + xcta(i, j + 1, k) + xcta(i + 1, j + 1, k) + xcta(i + 1, j + 1, k)) - 0.0625 * (zcta(i, j + 1, k) + zcta(i, j + 1, k) + zcta(i + 1, j + 1, k) + zcta(i + 1, j + 1, k)) * (xC(i, j + 1, k) + xC(i, j + 1, k) + xC(i + 1, j + 1, k) + xC(i + 1, j + 1, k));
                sy_CR1 =  0.0625 * (zC(i, j + 1, k) + zC(i, j + 2, k) + zC(i + 1, j + 1, k) + zC(i + 1, j + 2, k)) * (xcta(i, j + 1, k) + xcta(i, j + 2, k) + xcta(i + 1, j + 1, k) + xcta(i + 1, j + 2, k)) - 0.0625 * (zcta(i, j + 1, k) + zcta(i, j + 2, k) + zcta(i + 1, j + 1, k) + zcta(i + 1, j + 2, k)) * (xC(i, j + 1, k) + xC(i, j + 2, k) + xC(i + 1, j + 1, k) + xC(i + 1, j + 2, k));
                sy_CR2 =  0.0625 * (zC(i, j + 2, k) + zC(i, j + 3, k) + zC(i + 1, j + 2, k) + zC(i + 1, j + 3, k)) * (xcta(i, j + 2, k) + xcta(i, j + 3, k) + xcta(i + 1, j + 2, k) + xcta(i + 1, j + 3, k)) - 0.0625 * (zcta(i, j + 2, k) + zcta(i, j + 3, k) + zcta(i + 1, j + 2, k) + zcta(i + 1, j + 3, k)) * (xC(i, j + 2, k) + xC(i, j + 3, k) + xC(i + 1, j + 2, k) + xC(i + 1, j + 3, k));
                sz_CL2 =  0.0625 * (xC(i, j - 1, k) + xC(i, j - 2, k) + xC(i + 1, j - 1, k) + xC(i + 1, j - 2, k)) * (ycta(i, j - 1, k) + ycta(i, j - 2, k) + ycta(i + 1, j - 1, k) + ycta(i + 1, j - 2, k)) - 0.0625 * (xcta(i, j - 1, k) + xcta(i, j - 2, k) + xcta(i + 1, j - 1, k) + xcta(i + 1, j - 2, k)) * (yC(i, j - 1, k) + yC(i, j - 2, k) + yC(i + 1, j - 1, k) + yC(i + 1, j - 2, k));
                sz_CL1 =  0.0625 * (xC(i, j    , k) + xC(i, j - 1, k) + xC(i + 1, j    , k) + xC(i + 1, j - 1, k)) * (ycta(i, j    , k) + ycta(i, j - 1, k) + ycta(i + 1, j    , k) + ycta(i + 1, j - 1, k)) - 0.0625 * (xcta(i, j    , k) + xcta(i, j - 1, k) + xcta(i + 1, j    , k) + xcta(i + 1, j - 1, k)) * (yC(i, j    , k) + yC(i, j - 1, k) + yC(i + 1, j    , k) + yC(i + 1, j - 1, k));
                sz_fL  =  0.0625 * (xC(i, j    , k) + xC(i, j    , k) + xC(i + 1, j    , k) + xC(i + 1, j    , k)) * (ycta(i, j    , k) + ycta(i, j    , k) + ycta(i + 1, j    , k) + ycta(i + 1, j    , k)) - 0.0625 * (xcta(i, j    , k) + xcta(i, j    , k) + xcta(i + 1, j    , k) + xcta(i + 1, j    , k)) * (yC(i, j    , k) + yC(i, j    , k) + yC(i + 1, j    , k) + yC(i + 1, j    , k));
                sz_fR  =  0.0625 * (xC(i, j + 1, k) + xC(i, j + 1, k) + xC(i + 1, j + 1, k) + xC(i + 1, j + 1, k)) * (ycta(i, j + 1, k) + ycta(i, j + 1, k) + ycta(i + 1, j + 1, k) + ycta(i + 1, j + 1, k)) - 0.0625 * (xcta(i, j + 1, k) + xcta(i, j + 1, k) + xcta(i + 1, j + 1, k) + xcta(i + 1, j + 1, k)) * (yC(i, j + 1, k) + yC(i, j + 1, k) + yC(i + 1, j + 1, k) + yC(i + 1, j + 1, k));
                sz_CR1 =  0.0625 * (xC(i, j + 1, k) + xC(i, j + 2, k) + xC(i + 1, j + 1, k) + xC(i + 1, j + 2, k)) * (ycta(i, j + 1, k) + ycta(i, j + 2, k) + ycta(i + 1, j + 1, k) + ycta(i + 1, j + 2, k)) - 0.0625 * (xcta(i, j + 1, k) + xcta(i, j + 2, k) + xcta(i + 1, j + 1, k) + xcta(i + 1, j + 2, k)) * (yC(i, j + 1, k) + yC(i, j + 2, k) + yC(i + 1, j + 1, k) + yC(i + 1, j + 2, k));
                sz_CR2 =  0.0625 * (xC(i, j + 2, k) + xC(i, j + 3, k) + xC(i + 1, j + 2, k) + xC(i + 1, j + 3, k)) * (ycta(i, j + 2, k) + ycta(i, j + 3, k) + ycta(i + 1, j + 2, k) + ycta(i + 1, j + 3, k)) - 0.0625 * (xcta(i, j + 2, k) + xcta(i, j + 3, k) + xcta(i + 1, j + 2, k) + xcta(i + 1, j + 3, k)) * (yC(i, j + 2, k) + yC(i, j + 3, k) + yC(i + 1, j + 2, k) + yC(i + 1, j + 3, k));

                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));

                Cellxfv(i, j, k, 1) = KX;
                Cellyfv(i, j, k, 1) = KY;
                Cellzfv(i, j, k, 1) = KZ;
            }
        }
    }

    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int j = -1; j <= nj + 1; ++ j)
            {
                RDouble KX = 0.0;
                RDouble KY = 0.0;
                RDouble KZ = 0.0;

                RDouble sx_CL2 = 0.0625 * (ycta(i - 1, j, k) + ycta(i - 2, j, k) + ycta(i - 1, j + 1, k) + ycta(i - 2, j + 1, k)) * (zC(i - 1, j, k) + zC(i - 2, j, k) + zC(i - 1, j + 1, k) + zC(i - 2, j + 1, k)) - 0.0625 * (yC(i - 1, j, k) + yC(i - 2, j, k) + yC(i - 1, j + 1, k) + yC(i - 2, j + 1, k)) * (zcta(i - 1, j, k) + zcta(i - 2, j, k) + zcta(i - 1, j + 1, k) + zcta(i - 2, j + 1, k));
                RDouble sx_CL1 = 0.0625 * (ycta(i    , j, k) + ycta(i - 1, j, k) + ycta(i    , j + 1, k) + ycta(i - 1, j + 1, k)) * (zC(i    , j, k) + zC(i - 1, j, k) + zC(i    , j + 1, k) + zC(i - 1, j + 1, k)) - 0.0625 * (yC(i    , j, k) + yC(i - 1, j, k) + yC(i    , j + 1, k) + yC(i - 1, j + 1, k)) * (zcta(i    , j, k) + zcta(i - 1, j, k) + zcta(i    , j + 1, k) + zcta(i - 1, j + 1, k));
                RDouble sx_fL  = 0.0625 * (ycta(i    , j, k) + ycta(i    , j, k) + ycta(i    , j + 1, k) + ycta(i    , j + 1, k)) * (zC(i    , j, k) + zC(i    , j, k) + zC(i    , j + 1, k) + zC(i    , j + 1, k)) - 0.0625 * (yC(i    , j, k) + yC(i    , j, k) + yC(i    , j + 1, k) + yC(i    , j + 1, k)) * (zcta(i    , j, k) + zcta(i    , j, k) + zcta(i    , j + 1, k) + zcta(i    , j + 1, k));
                RDouble sx_fR  = 0.0625 * (ycta(i + 1, j, k) + ycta(i + 1, j, k) + ycta(i + 1, j + 1, k) + ycta(i + 1, j + 1, k)) * (zC(i + 1, j, k) + zC(i + 1, j, k) + zC(i + 1, j + 1, k) + zC(i + 1, j + 1, k)) - 0.0625 * (yC(i + 1, j, k) + yC(i + 1, j, k) + yC(i + 1, j + 1, k) + yC(i + 1, j + 1, k)) * (zcta(i + 1, j, k) + zcta(i + 1, j, k) + zcta(i + 1, j + 1, k) + zcta(i + 1, j + 1, k));
                RDouble sx_CR1 = 0.0625 * (ycta(i + 1, j, k) + ycta(i + 2, j, k) + ycta(i + 1, j + 1, k) + ycta(i + 2, j + 1, k)) * (zC(i + 1, j, k) + zC(i + 2, j, k) + zC(i + 1, j + 1, k) + zC(i + 2, j + 1, k)) - 0.0625 * (yC(i + 1, j, k) + yC(i + 2, j, k) + yC(i + 1, j + 1, k) + yC(i + 2, j + 1, k)) * (zcta(i + 1, j, k) + zcta(i + 2, j, k) + zcta(i + 1, j + 1, k) + zcta(i + 2, j + 1, k));
                RDouble sx_CR2 = 0.0625 * (ycta(i + 2, j, k) + ycta(i + 3, j, k) + ycta(i + 2, j + 1, k) + ycta(i + 3, j + 1, k)) * (zC(i + 2, j, k) + zC(i + 3, j, k) + zC(i + 2, j + 1, k) + zC(i + 3, j + 1, k)) - 0.0625 * (yC(i + 2, j, k) + yC(i + 3, j, k) + yC(i + 2, j + 1, k) + yC(i + 3, j + 1, k)) * (zcta(i + 2, j, k) + zcta(i + 3, j, k) + zcta(i + 2, j + 1, k) + zcta(i + 3, j + 1, k));
                RDouble sy_CL2 = 0.0625 * (zcta(i - 1, j, k) + zcta(i - 2, j, k) + zcta(i - 1, j + 1, k) + zcta(i - 2, j + 1, k)) * (xC(i - 1, j, k) + xC(i - 2, j, k) + xC(i - 1, j + 1, k) + xC(i - 2, j + 1, k)) - 0.0625 * (zC(i - 1, j, k) + zC(i - 2, j, k) + zC(i - 1, j + 1, k) + zC(i - 2, j + 1, k)) * (xcta(i - 1, j, k) + xcta(i - 2, j, k) + xcta(i - 1, j + 1, k) + xcta(i - 2, j + 1, k));
                RDouble sy_CL1 = 0.0625 * (zcta(i    , j, k) + zcta(i - 1, j, k) + zcta(i    , j + 1, k) + zcta(i - 1, j + 1, k)) * (xC(i    , j, k) + xC(i - 1, j, k) + xC(i    , j + 1, k) + xC(i - 1, j + 1, k)) - 0.0625 * (zC(i    , j, k) + zC(i - 1, j, k) + zC(i    , j + 1, k) + zC(i - 1, j + 1, k)) * (xcta(i    , j, k) + xcta(i - 1, j, k) + xcta(i    , j + 1, k) + xcta(i - 1, j + 1, k));
                RDouble sy_fL  = 0.0625 * (zcta(i    , j, k) + zcta(i    , j, k) + zcta(i    , j + 1, k) + zcta(i    , j + 1, k)) * (xC(i    , j, k) + xC(i    , j, k) + xC(i    , j + 1, k) + xC(i    , j + 1, k)) - 0.0625 * (zC(i    , j, k) + zC(i    , j, k) + zC(i    , j + 1, k) + zC(i    , j + 1, k)) * (xcta(i    , j, k) + xcta(i    , j, k) + xcta(i    , j + 1, k) + xcta(i    , j + 1, k));
                RDouble sy_fR  = 0.0625 * (zcta(i + 1, j, k) + zcta(i + 1, j, k) + zcta(i + 1, j + 1, k) + zcta(i + 1, j + 1, k)) * (xC(i + 1, j, k) + xC(i + 1, j, k) + xC(i + 1, j + 1, k) + xC(i + 1, j + 1, k)) - 0.0625 * (zC(i + 1, j, k) + zC(i + 1, j, k) + zC(i + 1, j + 1, k) + zC(i + 1, j + 1, k)) * (xcta(i + 1, j, k) + xcta(i + 1, j, k) + xcta(i + 1, j + 1, k) + xcta(i + 1, j + 1, k));
                RDouble sy_CR1 = 0.0625 * (zcta(i + 1, j, k) + zcta(i + 2, j, k) + zcta(i + 1, j + 1, k) + zcta(i + 2, j + 1, k)) * (xC(i + 1, j, k) + xC(i + 2, j, k) + xC(i + 1, j + 1, k) + xC(i + 2, j + 1, k)) - 0.0625 * (zC(i + 1, j, k) + zC(i + 2, j, k) + zC(i + 1, j + 1, k) + zC(i + 2, j + 1, k)) * (xcta(i + 1, j, k) + xcta(i + 2, j, k) + xcta(i + 1, j + 1, k) + xcta(i + 2, j + 1, k));
                RDouble sy_CR2 = 0.0625 * (zcta(i + 2, j, k) + zcta(i + 3, j, k) + zcta(i + 2, j + 1, k) + zcta(i + 3, j + 1, k)) * (xC(i + 2, j, k) + xC(i + 3, j, k) + xC(i + 2, j + 1, k) + xC(i + 3, j + 1, k)) - 0.0625 * (zC(i + 2, j, k) + zC(i + 3, j, k) + zC(i + 2, j + 1, k) + zC(i + 3, j + 1, k)) * (xcta(i + 2, j, k) + xcta(i + 3, j, k) + xcta(i + 2, j + 1, k) + xcta(i + 3, j + 1, k));
                RDouble sz_CL2 = 0.0625 * (xcta(i - 1, j, k) + xcta(i - 2, j, k) + xcta(i - 1, j + 1, k) + xcta(i - 2, j + 1, k)) * (yC(i - 1, j, k) + yC(i - 2, j, k) + yC(i - 1, j + 1, k) + yC(i - 2, j + 1, k)) - 0.0625 * (xC(i - 1, j, k) + xC(i - 2, j, k) + xC(i - 1, j + 1, k) + xC(i - 2, j + 1, k)) * (ycta(i - 1, j, k) + ycta(i - 2, j, k) + ycta(i - 1, j + 1, k) + ycta(i - 2, j + 1, k));
                RDouble sz_CL1 = 0.0625 * (xcta(i    , j, k) + xcta(i - 1, j, k) + xcta(i    , j + 1, k) + xcta(i - 1, j + 1, k)) * (yC(i    , j, k) + yC(i - 1, j, k) + yC(i    , j + 1, k) + yC(i - 1, j + 1, k)) - 0.0625 * (xC(i    , j, k) + xC(i - 1, j, k) + xC(i    , j + 1, k) + xC(i - 1, j + 1, k)) * (ycta(i    , j, k) + ycta(i - 1, j, k) + ycta(i    , j + 1, k) + ycta(i - 1, j + 1, k));
                RDouble sz_fL  = 0.0625 * (xcta(i    , j, k) + xcta(i    , j, k) + xcta(i    , j + 1, k) + xcta(i    , j + 1, k)) * (yC(i    , j, k) + yC(i    , j, k) + yC(i    , j + 1, k) + yC(i    , j + 1, k)) - 0.0625 * (xC(i    , j, k) + xC(i    , j, k) + xC(i    , j + 1, k) + xC(i    , j + 1, k)) * (ycta(i    , j, k) + ycta(i    , j, k) + ycta(i    , j + 1, k) + ycta(i    , j + 1, k));
                RDouble sz_fR  = 0.0625 * (xcta(i + 1, j, k) + xcta(i + 1, j, k) + xcta(i + 1, j + 1, k) + xcta(i + 1, j + 1, k)) * (yC(i + 1, j, k) + yC(i + 1, j, k) + yC(i + 1, j + 1, k) + yC(i + 1, j + 1, k)) - 0.0625 * (xC(i + 1, j, k) + xC(i + 1, j, k) + xC(i + 1, j + 1, k) + xC(i + 1, j + 1, k)) * (ycta(i + 1, j, k) + ycta(i + 1, j, k) + ycta(i + 1, j + 1, k) + ycta(i + 1, j + 1, k));
                RDouble sz_CR1 = 0.0625 * (xcta(i + 1, j, k) + xcta(i + 2, j, k) + xcta(i + 1, j + 1, k) + xcta(i + 2, j + 1, k)) * (yC(i + 1, j, k) + yC(i + 2, j, k) + yC(i + 1, j + 1, k) + yC(i + 2, j + 1, k)) - 0.0625 * (xC(i + 1, j, k) + xC(i + 2, j, k) + xC(i + 1, j + 1, k) + xC(i + 2, j + 1, k)) * (ycta(i + 1, j, k) + ycta(i + 2, j, k) + ycta(i + 1, j + 1, k) + ycta(i + 2, j + 1, k));
                RDouble sz_CR2 = 0.0625 * (xcta(i + 2, j, k) + xcta(i + 3, j, k) + xcta(i + 2, j + 1, k) + xcta(i + 3, j + 1, k)) * (yC(i + 2, j, k) + yC(i + 3, j, k) + yC(i + 2, j + 1, k) + yC(i + 3, j + 1, k)) - 0.0625 * (xC(i + 2, j, k) + xC(i + 3, j, k) + xC(i + 2, j + 1, k) + xC(i + 3, j + 1, k)) * (ycta(i + 2, j, k) + ycta(i + 3, j, k) + ycta(i + 2, j + 1, k) + ycta(i + 3, j + 1, k));

                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));

                sx_CL2 =  0.0625 * (yA(i, j, k - 1) + yA(i, j, k - 2) + yA(i, j + 1, k - 1) + yA(i, j + 1, k - 2)) * (zkxi(i, j, k - 1) + zkxi(i, j, k - 2) + zkxi(i, j + 1, k - 1) + zkxi(i, j + 1, k - 2)) - 0.0625 * (ykxi(i, j, k - 1) + ykxi(i, j, k - 2) + ykxi(i, j + 1, k - 1) + ykxi(i, j + 1, k - 2)) * (zA(i, j, k - 1) + zA(i, j, k - 2) + zA(i, j + 1, k - 1) + zA(i, j + 1, k - 2));
                sx_CL1 =  0.0625 * (yA(i, j, k   ) + yA(i, j, k - 1) + yA(i, j + 1, k   ) + yA(i, j + 1, k - 1)) * (zkxi(i, j, k   ) + zkxi(i, j, k - 1) + zkxi(i, j + 1, k   ) + zkxi(i, j + 1, k - 1)) - 0.0625 * (ykxi(i, j, k   ) + ykxi(i, j, k - 1) + ykxi(i, j + 1, k   ) + ykxi(i, j + 1, k - 1)) * (zA(i, j, k   ) + zA(i, j, k - 1) + zA(i, j + 1, k   ) + zA(i, j + 1, k - 1));
                sx_fL  =  0.0625 * (yA(i, j, k   ) + yA(i, j, k   ) + yA(i, j + 1, k   ) + yA(i, j + 1, k   )) * (zkxi(i, j, k   ) + zkxi(i, j, k   ) + zkxi(i, j + 1, k   ) + zkxi(i, j + 1, k   )) - 0.0625 * (ykxi(i, j, k   ) + ykxi(i, j, k   ) + ykxi(i, j + 1, k   ) + ykxi(i, j + 1, k   )) * (zA(i, j, k   ) + zA(i, j, k   ) + zA(i, j + 1, k   ) + zA(i, j + 1, k   ));
                sx_fR  =  0.0625 * (yA(i, j, k + 1) + yA(i, j, k + 1) + yA(i, j + 1, k + 1) + yA(i, j + 1, k + 1)) * (zkxi(i, j, k + 1) + zkxi(i, j, k + 1) + zkxi(i, j + 1, k + 1) + zkxi(i, j + 1, k + 1)) - 0.0625 * (ykxi(i, j, k + 1) + ykxi(i, j, k + 1) + ykxi(i, j + 1, k + 1) + ykxi(i, j + 1, k + 1)) * (zA(i, j, k + 1) + zA(i, j, k + 1) + zA(i, j + 1, k + 1) + zA(i, j + 1, k + 1));
                sx_CR1 =  0.0625 * (yA(i, j, k + 1) + yA(i, j, k + 2) + yA(i, j + 1, k + 1) + yA(i, j + 1, k + 2)) * (zkxi(i, j, k + 1) + zkxi(i, j, k + 2) + zkxi(i, j + 1, k + 1) + zkxi(i, j + 1, k + 2)) - 0.0625 * (ykxi(i, j, k + 1) + ykxi(i, j, k + 2) + ykxi(i, j + 1, k + 1) + ykxi(i, j + 1, k + 2)) * (zA(i, j, k + 1) + zA(i, j, k + 2) + zA(i, j + 1, k + 1) + zA(i, j + 1, k + 2));
                sx_CR2 =  0.0625 * (yA(i, j, k + 2) + yA(i, j, k + 3) + yA(i, j + 1, k + 2) + yA(i, j + 1, k + 3)) * (zkxi(i, j, k + 2) + zkxi(i, j, k + 3) + zkxi(i, j + 1, k + 2) + zkxi(i, j + 1, k + 3)) - 0.0625 * (ykxi(i, j, k + 2) + ykxi(i, j, k + 3) + ykxi(i, j + 1, k + 2) + ykxi(i, j + 1, k + 3)) * (zA(i, j, k + 2) + zA(i, j, k + 3) + zA(i, j + 1, k + 2) + zA(i, j + 1, k + 3));
                sy_CL2 =  0.0625 * (zA(i, j, k - 1) + zA(i, j, k - 2) + zA(i, j + 1, k - 1) + zA(i, j + 1, k - 2)) * (xkxi(i, j, k - 1) + xkxi(i, j, k - 2) + xkxi(i, j + 1, k - 1) + xkxi(i, j + 1, k - 2)) - 0.0625 * (zkxi(i, j, k - 1) + zkxi(i, j, k - 2) + zkxi(i, j + 1, k - 1) + zkxi(i, j + 1, k - 2)) * (xA(i, j, k - 1) + xA(i, j, k - 2) + xA(i, j + 1, k - 1) + xA(i, j + 1, k - 2));
                sy_CL1 =  0.0625 * (zA(i, j, k   ) + zA(i, j, k - 1) + zA(i, j + 1, k   ) + zA(i, j + 1, k - 1)) * (xkxi(i, j, k   ) + xkxi(i, j, k - 1) + xkxi(i, j + 1, k   ) + xkxi(i, j + 1, k - 1)) - 0.0625 * (zkxi(i, j, k   ) + zkxi(i, j, k - 1) + zkxi(i, j + 1, k   ) + zkxi(i, j + 1, k - 1)) * (xA(i, j, k   ) + xA(i, j, k - 1) + xA(i, j + 1, k   ) + xA(i, j + 1, k - 1));
                sy_fL  =  0.0625 * (zA(i, j, k   ) + zA(i, j, k   ) + zA(i, j + 1, k   ) + zA(i, j + 1, k   )) * (xkxi(i, j, k   ) + xkxi(i, j, k   ) + xkxi(i, j + 1, k   ) + xkxi(i, j + 1, k   )) - 0.0625 * (zkxi(i, j, k   ) + zkxi(i, j, k   ) + zkxi(i, j + 1, k   ) + zkxi(i, j + 1, k   )) * (xA(i, j, k   ) + xA(i, j, k   ) + xA(i, j + 1, k   ) + xA(i, j + 1, k   ));
                sy_fR  =  0.0625 * (zA(i, j, k + 1) + zA(i, j, k + 1) + zA(i, j + 1, k + 1) + zA(i, j + 1, k + 1)) * (xkxi(i, j, k + 1) + xkxi(i, j, k + 1) + xkxi(i, j + 1, k + 1) + xkxi(i, j + 1, k + 1)) - 0.0625 * (zkxi(i, j, k + 1) + zkxi(i, j, k + 1) + zkxi(i, j + 1, k + 1) + zkxi(i, j + 1, k + 1)) * (xA(i, j, k + 1) + xA(i, j, k + 1) + xA(i, j + 1, k + 1) + xA(i, j + 1, k + 1));
                sy_CR1 =  0.0625 * (zA(i, j, k + 1) + zA(i, j, k + 2) + zA(i, j + 1, k + 1) + zA(i, j + 1, k + 2)) * (xkxi(i, j, k + 1) + xkxi(i, j, k + 2) + xkxi(i, j + 1, k + 1) + xkxi(i, j + 1, k + 2)) - 0.0625 * (zkxi(i, j, k + 1) + zkxi(i, j, k + 2) + zkxi(i, j + 1, k + 1) + zkxi(i, j + 1, k + 2)) * (xA(i, j, k + 1) + xA(i, j, k + 2) + xA(i, j + 1, k + 1) + xA(i, j + 1, k + 2));
                sy_CR2 =  0.0625 * (zA(i, j, k + 2) + zA(i, j, k + 3) + zA(i, j + 1, k + 2) + zA(i, j + 1, k + 3)) * (xkxi(i, j, k + 2) + xkxi(i, j, k + 3) + xkxi(i, j + 1, k + 2) + xkxi(i, j + 1, k + 3)) - 0.0625 * (zkxi(i, j, k + 2) + zkxi(i, j, k + 3) + zkxi(i, j + 1, k + 2) + zkxi(i, j + 1, k + 3)) * (xA(i, j, k + 2) + xA(i, j, k + 3) + xA(i, j + 1, k + 2) + xA(i, j + 1, k + 3));
                sz_CL2 =  0.0625 * (xA(i, j, k - 1) + xA(i, j, k - 2) + xA(i, j + 1, k - 1) + xA(i, j + 1, k - 2)) * (ykxi(i, j, k - 1) + ykxi(i, j, k - 2) + ykxi(i, j + 1, k - 1) + ykxi(i, j + 1, k - 2)) - 0.0625 * (xkxi(i, j, k - 1) + xkxi(i, j, k - 2) + xkxi(i, j + 1, k - 1) + xkxi(i, j + 1, k - 2)) * (yA(i, j, k - 1) + yA(i, j, k - 2) + yA(i, j + 1, k - 1) + yA(i, j + 1, k - 2));
                sz_CL1 =  0.0625 * (xA(i, j, k   ) + xA(i, j, k - 1) + xA(i, j + 1, k   ) + xA(i, j + 1, k - 1)) * (ykxi(i, j, k   ) + ykxi(i, j, k - 1) + ykxi(i, j + 1, k   ) + ykxi(i, j + 1, k - 1)) - 0.0625 * (xkxi(i, j, k   ) + xkxi(i, j, k - 1) + xkxi(i, j + 1, k   ) + xkxi(i, j + 1, k - 1)) * (yA(i, j, k   ) + yA(i, j, k - 1) + yA(i, j + 1, k   ) + yA(i, j + 1, k - 1));
                sz_fL  =  0.0625 * (xA(i, j, k   ) + xA(i, j, k   ) + xA(i, j + 1, k   ) + xA(i, j + 1, k   )) * (ykxi(i, j, k   ) + ykxi(i, j, k   ) + ykxi(i, j + 1, k   ) + ykxi(i, j + 1, k   )) - 0.0625 * (xkxi(i, j, k   ) + xkxi(i, j, k   ) + xkxi(i, j + 1, k   ) + xkxi(i, j + 1, k   )) * (yA(i, j, k   ) + yA(i, j, k   ) + yA(i, j + 1, k   ) + yA(i, j + 1, k   ));
                sz_fR  =  0.0625 * (xA(i, j, k + 1) + xA(i, j, k + 1) + xA(i, j + 1, k + 1) + xA(i, j + 1, k + 1)) * (ykxi(i, j, k + 1) + ykxi(i, j, k + 1) + ykxi(i, j + 1, k + 1) + ykxi(i, j + 1, k + 1)) - 0.0625 * (xkxi(i, j, k + 1) + xkxi(i, j, k + 1) + xkxi(i, j + 1, k + 1) + xkxi(i, j + 1, k + 1)) * (yA(i, j, k + 1) + yA(i, j, k + 1) + yA(i, j + 1, k + 1) + yA(i, j + 1, k + 1));
                sz_CR1 =  0.0625 * (xA(i, j, k + 1) + xA(i, j, k + 2) + xA(i, j + 1, k + 1) + xA(i, j + 1, k + 2)) * (ykxi(i, j, k + 1) + ykxi(i, j, k + 2) + ykxi(i, j + 1, k + 1) + ykxi(i, j + 1, k + 2)) - 0.0625 * (xkxi(i, j, k + 1) + xkxi(i, j, k + 2) + xkxi(i, j + 1, k + 1) + xkxi(i, j + 1, k + 2)) * (yA(i, j, k + 1) + yA(i, j, k + 2) + yA(i, j + 1, k + 1) + yA(i, j + 1, k + 2));
                sz_CR2 =  0.0625 * (xA(i, j, k + 2) + xA(i, j, k + 3) + xA(i, j + 1, k + 2) + xA(i, j + 1, k + 3)) * (ykxi(i, j, k + 2) + ykxi(i, j, k + 3) + ykxi(i, j + 1, k + 2) + ykxi(i, j + 1, k + 3)) - 0.0625 * (xkxi(i, j, k + 2) + xkxi(i, j, k + 3) + xkxi(i, j + 1, k + 2) + xkxi(i, j + 1, k + 3)) * (yA(i, j, k + 2) + yA(i, j, k + 3) + yA(i, j + 1, k + 2) + yA(i, j + 1, k + 3));

                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));

                Cellxfv(i, j, k, 2) = KX;
                Cellyfv(i, j, k, 2) = KY;
                Cellzfv(i, j, k, 2) = KZ;
            }
        }
    }

    for (int j = 1; j <= nj - 1; ++ j)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int k = -1; k <= nk + 1; ++ k)
            {
                RDouble KX = 0.0;
                RDouble KY = 0.0;
                RDouble KZ = 0.0;

                RDouble sx_CL2 = 0.0625 * (ykxi(i, j - 1, k) + ykxi(i, j - 2, k) + ykxi(i, j - 1, k + 1) + ykxi(i, j - 2, k + 1)) * (zA(i, j - 1, k) + zA(i, j - 2, k) + zA(i, j - 1, k + 1) + zA(i, j - 2, k + 1)) - 0.0625 * (yA(i, j - 1, k) + yA(i, j - 2, k) + yA(i, j - 1, k + 1) + yA(i, j - 2, k + 1)) * (zkxi(i, j - 1, k) + zkxi(i, j - 2, k) + zkxi(i, j - 1, k + 1) + zkxi(i, j - 2, k + 1));
                RDouble sx_CL1 = 0.0625 * (ykxi(i, j    , k) + ykxi(i, j - 1, k) + ykxi(i, j    , k + 1) + ykxi(i, j - 1, k + 1)) * (zA(i, j    , k) + zA(i, j - 1, k) + zA(i, j    , k + 1) + zA(i, j - 1, k + 1)) - 0.0625 * (yA(i, j    , k) + yA(i, j - 1, k) + yA(i, j    , k + 1) + yA(i, j - 1, k + 1)) * (zkxi(i, j    , k) + zkxi(i, j - 1, k) + zkxi(i, j    , k + 1) + zkxi(i, j - 1, k + 1));
                RDouble sx_fL  = 0.0625 * (ykxi(i, j    , k) + ykxi(i, j    , k) + ykxi(i, j    , k + 1) + ykxi(i, j    , k + 1)) * (zA(i, j    , k) + zA(i, j    , k) + zA(i, j    , k + 1) + zA(i, j    , k + 1)) - 0.0625 * (yA(i, j    , k) + yA(i, j    , k) + yA(i, j    , k + 1) + yA(i, j    , k + 1)) * (zkxi(i, j    , k) + zkxi(i, j    , k) + zkxi(i, j    , k + 1) + zkxi(i, j    , k + 1));
                RDouble sx_fR  = 0.0625 * (ykxi(i, j + 1, k) + ykxi(i, j + 1, k) + ykxi(i, j + 1, k + 1) + ykxi(i, j + 1, k + 1)) * (zA(i, j + 1, k) + zA(i, j + 1, k) + zA(i, j + 1, k + 1) + zA(i, j + 1, k + 1)) - 0.0625 * (yA(i, j + 1, k) + yA(i, j + 1, k) + yA(i, j + 1, k + 1) + yA(i, j + 1, k + 1)) * (zkxi(i, j + 1, k) + zkxi(i, j + 1, k) + zkxi(i, j + 1, k + 1) + zkxi(i, j + 1, k + 1));
                RDouble sx_CR1 = 0.0625 * (ykxi(i, j + 1, k) + ykxi(i, j + 2, k) + ykxi(i, j + 1, k + 1) + ykxi(i, j + 2, k + 1)) * (zA(i, j + 1, k) + zA(i, j + 2, k) + zA(i, j + 1, k + 1) + zA(i, j + 2, k + 1)) - 0.0625 * (yA(i, j + 1, k) + yA(i, j + 2, k) + yA(i, j + 1, k + 1) + yA(i, j + 2, k + 1)) * (zkxi(i, j + 1, k) + zkxi(i, j + 2, k) + zkxi(i, j + 1, k + 1) + zkxi(i, j + 2, k + 1));
                RDouble sx_CR2 = 0.0625 * (ykxi(i, j + 2, k) + ykxi(i, j + 3, k) + ykxi(i, j + 2, k + 1) + ykxi(i, j + 3, k + 1)) * (zA(i, j + 2, k) + zA(i, j + 3, k) + zA(i, j + 2, k + 1) + zA(i, j + 3, k + 1)) - 0.0625 * (yA(i, j + 2, k) + yA(i, j + 3, k) + yA(i, j + 2, k + 1) + yA(i, j + 3, k + 1)) * (zkxi(i, j + 2, k) + zkxi(i, j + 3, k) + zkxi(i, j + 2, k + 1) + zkxi(i, j + 3, k + 1));
                RDouble sy_CL2 = 0.0625 * (zkxi(i, j - 1, k) + zkxi(i, j - 2, k) + zkxi(i, j - 1, k + 1) + zkxi(i, j - 2, k + 1)) * (xA(i, j - 1, k) + xA(i, j - 2, k) + xA(i, j - 1, k + 1) + xA(i, j - 2, k + 1)) - 0.0625 * (zA(i, j - 1, k) + zA(i, j - 2, k) + zA(i, j - 1, k + 1) + zA(i, j - 2, k + 1)) * (xkxi(i, j - 1, k) + xkxi(i, j - 2, k) + xkxi(i, j - 1, k + 1) + xkxi(i, j - 2, k + 1));
                RDouble sy_CL1 = 0.0625 * (zkxi(i, j    , k) + zkxi(i, j - 1, k) + zkxi(i, j    , k + 1) + zkxi(i, j - 1, k + 1)) * (xA(i, j    , k) + xA(i, j - 1, k) + xA(i, j    , k + 1) + xA(i, j - 1, k + 1)) - 0.0625 * (zA(i, j    , k) + zA(i, j - 1, k) + zA(i, j    , k + 1) + zA(i, j - 1, k + 1)) * (xkxi(i, j    , k) + xkxi(i, j - 1, k) + xkxi(i, j    , k + 1) + xkxi(i, j - 1, k + 1));
                RDouble sy_fL  = 0.0625 * (zkxi(i, j    , k) + zkxi(i, j    , k) + zkxi(i, j    , k + 1) + zkxi(i, j    , k + 1)) * (xA(i, j    , k) + xA(i, j    , k) + xA(i, j    , k + 1) + xA(i, j    , k + 1)) - 0.0625 * (zA(i, j    , k) + zA(i, j    , k) + zA(i, j    , k + 1) + zA(i, j    , k + 1)) * (xkxi(i, j    , k) + xkxi(i, j    , k) + xkxi(i, j    , k + 1) + xkxi(i, j    , k + 1));
                RDouble sy_fR  = 0.0625 * (zkxi(i, j + 1, k) + zkxi(i, j + 1, k) + zkxi(i, j + 1, k + 1) + zkxi(i, j + 1, k + 1)) * (xA(i, j + 1, k) + xA(i, j + 1, k) + xA(i, j + 1, k + 1) + xA(i, j + 1, k + 1)) - 0.0625 * (zA(i, j + 1, k) + zA(i, j + 1, k) + zA(i, j + 1, k + 1) + zA(i, j + 1, k + 1)) * (xkxi(i, j + 1, k) + xkxi(i, j + 1, k) + xkxi(i, j + 1, k + 1) + xkxi(i, j + 1, k + 1));
                RDouble sy_CR1 = 0.0625 * (zkxi(i, j + 1, k) + zkxi(i, j + 2, k) + zkxi(i, j + 1, k + 1) + zkxi(i, j + 2, k + 1)) * (xA(i, j + 1, k) + xA(i, j + 2, k) + xA(i, j + 1, k + 1) + xA(i, j + 2, k + 1)) - 0.0625 * (zA(i, j + 1, k) + zA(i, j + 2, k) + zA(i, j + 1, k + 1) + zA(i, j + 2, k + 1)) * (xkxi(i, j + 1, k) + xkxi(i, j + 2, k) + xkxi(i, j + 1, k + 1) + xkxi(i, j + 2, k + 1));
                RDouble sy_CR2 = 0.0625 * (zkxi(i, j + 2, k) + zkxi(i, j + 3, k) + zkxi(i, j + 2, k + 1) + zkxi(i, j + 3, k + 1)) * (xA(i, j + 2, k) + xA(i, j + 3, k) + xA(i, j + 2, k + 1) + xA(i, j + 3, k + 1)) - 0.0625 * (zA(i, j + 2, k) + zA(i, j + 3, k) + zA(i, j + 2, k + 1) + zA(i, j + 3, k + 1)) * (xkxi(i, j + 2, k) + xkxi(i, j + 3, k) + xkxi(i, j + 2, k + 1) + xkxi(i, j + 3, k + 1));
                RDouble sz_CL2 = 0.0625 * (xkxi(i, j - 1, k) + xkxi(i, j - 2, k) + xkxi(i, j - 1, k + 1) + xkxi(i, j - 2, k + 1)) * (yA(i, j - 1, k) + yA(i, j - 2, k) + yA(i, j - 1, k + 1) + yA(i, j - 2, k + 1)) - 0.0625 * (xA(i, j - 1, k) + xA(i, j - 2, k) + xA(i, j - 1, k + 1) + xA(i, j - 2, k + 1)) * (ykxi(i, j - 1, k) + ykxi(i, j - 2, k) + ykxi(i, j - 1, k + 1) + ykxi(i, j - 2, k + 1));
                RDouble sz_CL1 = 0.0625 * (xkxi(i, j    , k) + xkxi(i, j - 1, k) + xkxi(i, j    , k + 1) + xkxi(i, j - 1, k + 1)) * (yA(i, j    , k) + yA(i, j - 1, k) + yA(i, j    , k + 1) + yA(i, j - 1, k + 1)) - 0.0625 * (xA(i, j    , k) + xA(i, j - 1, k) + xA(i, j    , k + 1) + xA(i, j - 1, k + 1)) * (ykxi(i, j    , k) + ykxi(i, j - 1, k) + ykxi(i, j    , k + 1) + ykxi(i, j - 1, k + 1));
                RDouble sz_fL  = 0.0625 * (xkxi(i, j    , k) + xkxi(i, j    , k) + xkxi(i, j    , k + 1) + xkxi(i, j    , k + 1)) * (yA(i, j    , k) + yA(i, j    , k) + yA(i, j    , k + 1) + yA(i, j    , k + 1)) - 0.0625 * (xA(i, j    , k) + xA(i, j    , k) + xA(i, j    , k + 1) + xA(i, j    , k + 1)) * (ykxi(i, j    , k) + ykxi(i, j    , k) + ykxi(i, j    , k + 1) + ykxi(i, j    , k + 1));
                RDouble sz_fR  = 0.0625 * (xkxi(i, j + 1, k) + xkxi(i, j + 1, k) + xkxi(i, j + 1, k + 1) + xkxi(i, j + 1, k + 1)) * (yA(i, j + 1, k) + yA(i, j + 1, k) + yA(i, j + 1, k + 1) + yA(i, j + 1, k + 1)) - 0.0625 * (xA(i, j + 1, k) + xA(i, j + 1, k) + xA(i, j + 1, k + 1) + xA(i, j + 1, k + 1)) * (ykxi(i, j + 1, k) + ykxi(i, j + 1, k) + ykxi(i, j + 1, k + 1) + ykxi(i, j + 1, k + 1));
                RDouble sz_CR1 = 0.0625 * (xkxi(i, j + 1, k) + xkxi(i, j + 2, k) + xkxi(i, j + 1, k + 1) + xkxi(i, j + 2, k + 1)) * (yA(i, j + 1, k) + yA(i, j + 2, k) + yA(i, j + 1, k + 1) + yA(i, j + 2, k + 1)) - 0.0625 * (xA(i, j + 1, k) + xA(i, j + 2, k) + xA(i, j + 1, k + 1) + xA(i, j + 2, k + 1)) * (ykxi(i, j + 1, k) + ykxi(i, j + 2, k) + ykxi(i, j + 1, k + 1) + ykxi(i, j + 2, k + 1));
                RDouble sz_CR2 = 0.0625 * (xkxi(i, j + 2, k) + xkxi(i, j + 3, k) + xkxi(i, j + 2, k + 1) + xkxi(i, j + 3, k + 1)) * (yA(i, j + 2, k) + yA(i, j + 3, k) + yA(i, j + 2, k + 1) + yA(i, j + 3, k + 1)) - 0.0625 * (xA(i, j + 2, k) + xA(i, j + 3, k) + xA(i, j + 2, k + 1) + xA(i, j + 3, k + 1)) * (ykxi(i, j + 2, k) + ykxi(i, j + 3, k) + ykxi(i, j + 2, k + 1) + ykxi(i, j + 3, k + 1));

                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));

                sx_CL2 =  0.0625 * (yB(i - 1, j, k) + yB(i - 2, j, k) + yB(i - 1, j, k + 1) + yB(i - 2, j, k + 1)) * (zeta(i - 1, j, k) + zeta(i - 2, j, k) + zeta(i - 1, j, k + 1) + zeta(i - 2, j, k + 1)) - 0.0625 * (yeta(i - 1, j, k) + yeta(i - 2, j, k) + yeta(i - 1, j, k + 1) + yeta(i - 2, j, k + 1)) * (zB(i - 1, j, k) + zB(i - 2, j, k) + zB(i - 1, j, k + 1) + zB(i - 2, j, k + 1));
                sx_CL1 =  0.0625 * (yB(i    , j, k) + yB(i - 1, j, k) + yB(i    , j, k + 1) + yB(i - 1, j, k + 1)) * (zeta(i    , j, k) + zeta(i - 1, j, k) + zeta(i    , j, k + 1) + zeta(i - 1, j, k + 1)) - 0.0625 * (yeta(i    , j, k) + yeta(i - 1, j, k) + yeta(i    , j, k + 1) + yeta(i - 1, j, k + 1)) * (zB(i    , j, k) + zB(i - 1, j, k) + zB(i    , j, k + 1) + zB(i - 1, j, k + 1));
                sx_fL  =  0.0625 * (yB(i    , j, k) + yB(i    , j, k) + yB(i    , j, k + 1) + yB(i    , j, k + 1)) * (zeta(i    , j, k) + zeta(i    , j, k) + zeta(i    , j, k + 1) + zeta(i    , j, k + 1)) - 0.0625 * (yeta(i    , j, k) + yeta(i    , j, k) + yeta(i    , j, k + 1) + yeta(i    , j, k + 1)) * (zB(i    , j, k) + zB(i    , j, k) + zB(i    , j, k + 1) + zB(i    , j, k + 1));
                sx_fR  =  0.0625 * (yB(i + 1, j, k) + yB(i + 1, j, k) + yB(i + 1, j, k + 1) + yB(i + 1, j, k + 1)) * (zeta(i + 1, j, k) + zeta(i + 1, j, k) + zeta(i + 1, j, k + 1) + zeta(i + 1, j, k + 1)) - 0.0625 * (yeta(i + 1, j, k) + yeta(i + 1, j, k) + yeta(i + 1, j, k + 1) + yeta(i + 1, j, k + 1)) * (zB(i + 1, j, k) + zB(i + 1, j, k) + zB(i + 1, j, k + 1) + zB(i + 1, j, k + 1));
                sx_CR1 =  0.0625 * (yB(i + 1, j, k) + yB(i + 2, j, k) + yB(i + 1, j, k + 1) + yB(i + 2, j, k + 1)) * (zeta(i + 1, j, k) + zeta(i + 2, j, k) + zeta(i + 1, j, k + 1) + zeta(i + 2, j, k + 1)) - 0.0625 * (yeta(i + 1, j, k) + yeta(i + 2, j, k) + yeta(i + 1, j, k + 1) + yeta(i + 2, j, k + 1)) * (zB(i + 1, j, k) + zB(i + 2, j, k) + zB(i + 1, j, k + 1) + zB(i + 2, j, k + 1));
                sx_CR2 =  0.0625 * (yB(i + 2, j, k) + yB(i + 3, j, k) + yB(i + 2, j, k + 1) + yB(i + 3, j, k + 1)) * (zeta(i + 2, j, k) + zeta(i + 3, j, k) + zeta(i + 2, j, k + 1) + zeta(i + 3, j, k + 1)) - 0.0625 * (yeta(i + 2, j, k) + yeta(i + 3, j, k) + yeta(i + 2, j, k + 1) + yeta(i + 3, j, k + 1)) * (zB(i + 2, j, k) + zB(i + 3, j, k) + zB(i + 2, j, k + 1) + zB(i + 3, j, k + 1));
                sy_CL2 =  0.0625 * (zB(i - 1, j, k) + zB(i - 2, j, k) + zB(i - 1, j, k + 1) + zB(i - 2, j, k + 1)) * (xeta(i - 1, j, k) + xeta(i - 2, j, k) + xeta(i - 1, j, k + 1) + xeta(i - 2, j, k + 1)) - 0.0625 * (zeta(i - 1, j, k) + zeta(i - 2, j, k) + zeta(i - 1, j, k + 1) + zeta(i - 2, j, k + 1)) * (xB(i - 1, j, k) + xB(i - 2, j, k) + xB(i - 1, j, k + 1) + xB(i - 2, j, k + 1));
                sy_CL1 =  0.0625 * (zB(i    , j, k) + zB(i - 1, j, k) + zB(i    , j, k + 1) + zB(i - 1, j, k + 1)) * (xeta(i    , j, k) + xeta(i - 1, j, k) + xeta(i    , j, k + 1) + xeta(i - 1, j, k + 1)) - 0.0625 * (zeta(i    , j, k) + zeta(i - 1, j, k) + zeta(i    , j, k + 1) + zeta(i - 1, j, k + 1)) * (xB(i    , j, k) + xB(i - 1, j, k) + xB(i    , j, k + 1) + xB(i - 1, j, k + 1));
                sy_fL  =  0.0625 * (zB(i    , j, k) + zB(i    , j, k) + zB(i    , j, k + 1) + zB(i    , j, k + 1)) * (xeta(i    , j, k) + xeta(i    , j, k) + xeta(i    , j, k + 1) + xeta(i    , j, k + 1)) - 0.0625 * (zeta(i    , j, k) + zeta(i    , j, k) + zeta(i    , j, k + 1) + zeta(i    , j, k + 1)) * (xB(i    , j, k) + xB(i    , j, k) + xB(i    , j, k + 1) + xB(i    , j, k + 1));
                sy_fR  =  0.0625 * (zB(i + 1, j, k) + zB(i + 1, j, k) + zB(i + 1, j, k + 1) + zB(i + 1, j, k + 1)) * (xeta(i + 1, j, k) + xeta(i + 1, j, k) + xeta(i + 1, j, k + 1) + xeta(i + 1, j, k + 1)) - 0.0625 * (zeta(i + 1, j, k) + zeta(i + 1, j, k) + zeta(i + 1, j, k + 1) + zeta(i + 1, j, k + 1)) * (xB(i + 1, j, k) + xB(i + 1, j, k) + xB(i + 1, j, k + 1) + xB(i + 1, j, k + 1));
                sy_CR1 =  0.0625 * (zB(i + 1, j, k) + zB(i + 2, j, k) + zB(i + 1, j, k + 1) + zB(i + 2, j, k + 1)) * (xeta(i + 1, j, k) + xeta(i + 2, j, k) + xeta(i + 1, j, k + 1) + xeta(i + 2, j, k + 1)) - 0.0625 * (zeta(i + 1, j, k) + zeta(i + 2, j, k) + zeta(i + 1, j, k + 1) + zeta(i + 2, j, k + 1)) * (xB(i + 1, j, k) + xB(i + 2, j, k) + xB(i + 1, j, k + 1) + xB(i + 2, j, k + 1));
                sy_CR2 =  0.0625 * (zB(i + 2, j, k) + zB(i + 3, j, k) + zB(i + 2, j, k + 1) + zB(i + 3, j, k + 1)) * (xeta(i + 2, j, k) + xeta(i + 3, j, k) + xeta(i + 2, j, k + 1) + xeta(i + 3, j, k + 1)) - 0.0625 * (zeta(i + 2, j, k) + zeta(i + 3, j, k) + zeta(i + 2, j, k + 1) + zeta(i + 3, j, k + 1)) * (xB(i + 2, j, k) + xB(i + 3, j, k) + xB(i + 2, j, k + 1) + xB(i + 3, j, k + 1));
                sz_CL2 =  0.0625 * (xB(i - 1, j, k) + xB(i - 2, j, k) + xB(i - 1, j, k + 1) + xB(i - 2, j, k + 1)) * (yeta(i - 1, j, k) + yeta(i - 2, j, k) + yeta(i - 1, j, k + 1) + yeta(i - 2, j, k + 1)) - 0.0625 * (xeta(i - 1, j, k) + xeta(i - 2, j, k) + xeta(i - 1, j, k + 1) + xeta(i - 2, j, k + 1)) * (yB(i - 1, j, k) + yB(i - 2, j, k) + yB(i - 1, j, k + 1) + yB(i - 2, j, k + 1));
                sz_CL1 =  0.0625 * (xB(i    , j, k) + xB(i - 1, j, k) + xB(i    , j, k + 1) + xB(i - 1, j, k + 1)) * (yeta(i    , j, k) + yeta(i - 1, j, k) + yeta(i    , j, k + 1) + yeta(i - 1, j, k + 1)) - 0.0625 * (xeta(i    , j, k) + xeta(i - 1, j, k) + xeta(i    , j, k + 1) + xeta(i - 1, j, k + 1)) * (yB(i    , j, k) + yB(i - 1, j, k) + yB(i    , j, k + 1) + yB(i - 1, j, k + 1));
                sz_fL  =  0.0625 * (xB(i    , j, k) + xB(i    , j, k) + xB(i    , j, k + 1) + xB(i    , j, k + 1)) * (yeta(i    , j, k) + yeta(i    , j, k) + yeta(i    , j, k + 1) + yeta(i    , j, k + 1)) - 0.0625 * (xeta(i    , j, k) + xeta(i    , j, k) + xeta(i    , j, k + 1) + xeta(i    , j, k + 1)) * (yB(i    , j, k) + yB(i    , j, k) + yB(i    , j, k + 1) + yB(i    , j, k + 1));
                sz_fR  =  0.0625 * (xB(i + 1, j, k) + xB(i + 1, j, k) + xB(i + 1, j, k + 1) + xB(i + 1, j, k + 1)) * (yeta(i + 1, j, k) + yeta(i + 1, j, k) + yeta(i + 1, j, k + 1) + yeta(i + 1, j, k + 1)) - 0.0625 * (xeta(i + 1, j, k) + xeta(i + 1, j, k) + xeta(i + 1, j, k + 1) + xeta(i + 1, j, k + 1)) * (yB(i + 1, j, k) + yB(i + 1, j, k) + yB(i + 1, j, k + 1) + yB(i + 1, j, k + 1));
                sz_CR1 =  0.0625 * (xB(i + 1, j, k) + xB(i + 2, j, k) + xB(i + 1, j, k + 1) + xB(i + 2, j, k + 1)) * (yeta(i + 1, j, k) + yeta(i + 2, j, k) + yeta(i + 1, j, k + 1) + yeta(i + 2, j, k + 1)) - 0.0625 * (xeta(i + 1, j, k) + xeta(i + 2, j, k) + xeta(i + 1, j, k + 1) + xeta(i + 2, j, k + 1)) * (yB(i + 1, j, k) + yB(i + 2, j, k) + yB(i + 1, j, k + 1) + yB(i + 2, j, k + 1));
                sz_CR2 =  0.0625 * (xB(i + 2, j, k) + xB(i + 3, j, k) + xB(i + 2, j, k + 1) + xB(i + 3, j, k + 1)) * (yeta(i + 2, j, k) + yeta(i + 3, j, k) + yeta(i + 2, j, k + 1) + yeta(i + 3, j, k + 1)) - 0.0625 * (xeta(i + 2, j, k) + xeta(i + 3, j, k) + xeta(i + 2, j, k + 1) + xeta(i + 3, j, k + 1)) * (yB(i + 2, j, k) + yB(i + 3, j, k) + yB(i + 2, j, k + 1) + yB(i + 3, j, k + 1));

                KX += half * (alfa_sixth * (sx_fR - sx_fL) + a_sixth * (sx_CR1 - sx_CL1) + b_sixth * (sx_CR2 - sx_CL2));
                KY += half * (alfa_sixth * (sy_fR - sy_fL) + a_sixth * (sy_CR1 - sy_CL1) + b_sixth * (sy_CR2 - sy_CL2));
                KZ += half * (alfa_sixth * (sz_fR - sz_fL) + a_sixth * (sz_CR1 - sz_CL1) + b_sixth * (sz_CR2 - sz_CL2));
                
                Cellxfv(i, j, k, 3) = KX;
                Cellyfv(i, j, k, 3) = KY;
                Cellzfv(i, j, k, 3) = KZ;
            }
        }
    }


    delete &xA   ;
    delete &yA   ;
    delete &zA   ;
    delete &xkxi;
    delete &ykxi;
    delete &zkxi;

    delete &xB   ;
    delete &yB   ;
    delete &zB   ;
    delete &xeta;
    delete &yeta;
    delete &zeta;

    delete &xC   ;
    delete &yC   ;
    delete &zC   ;
    delete &xcta;
    delete &ycta;
    delete &zcta;
    //!=
    RDouble3D &vol = *this->GetCellVolume();

    RDouble3D &xfaceI = *new RDouble3D(Range(-1, ni + 2), Range(1, nj - 1), Range(1, nk - 1), fortranArray);
    RDouble3D &yfaceI = *new RDouble3D(Range(-1, ni + 2), Range(1, nj - 1), Range(1, nk - 1), fortranArray);
    RDouble3D &zfaceI = *new RDouble3D(Range(-1, ni + 2), Range(1, nj - 1), Range(1, nk - 1), fortranArray);

    RDouble3D &xfaceJ = *new RDouble3D(Range(1, ni - 1), Range(-1, nj + 2), Range(1, nk - 1), fortranArray);
    RDouble3D &yfaceJ = *new RDouble3D(Range(1, ni - 1), Range(-1, nj + 2), Range(1, nk - 1), fortranArray);
    RDouble3D &zfaceJ = *new RDouble3D(Range(1, ni - 1), Range(-1, nj + 2), Range(1, nk - 1), fortranArray);

    RDouble3D &xfaceK = *new RDouble3D(Range(1, ni - 1), Range(1, nj - 1), Range(-1, nk + 2), fortranArray);
    RDouble3D &yfaceK = *new RDouble3D(Range(1, ni - 1), Range(1, nj - 1), Range(-1, nk + 2), fortranArray);
    RDouble3D &zfaceK = *new RDouble3D(Range(1, ni - 1), Range(1, nj - 1), Range(-1, nk + 2), fortranArray);

    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = -1; i <= ni + 2; ++ i)
            {
                RDouble x1 = xx(i, j,     k);
                RDouble x2 = xx(i, j + 1, k);
                RDouble x3 = xx(i, j,     k + 1);
                RDouble x4 = xx(i, j + 1, k + 1);

                RDouble y1 = yy(i, j,     k);
                RDouble y2 = yy(i, j + 1, k);
                RDouble y3 = yy(i, j,     k + 1);
                RDouble y4 = yy(i, j + 1, k + 1);

                RDouble z1 = zz(i, j,     k);
                RDouble z2 = zz(i, j + 1, k);
                RDouble z3 = zz(i, j,     k + 1);
                RDouble z4 = zz(i, j + 1, k + 1);
                
                xfaceI(i, j, k) = 0.25 * (x1 + x2 + x3 + x4);
                yfaceI(i, j, k) = 0.25 * (y1 + y2 + y3 + y4);
                zfaceI(i, j, k) = 0.25 * (z1 + z2 + z3 + z4);
            }
        }
    }

    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int j = -1; j <= nj + 2; ++ j)
            {
                RDouble x1 = xx(i,     j, k);
                RDouble x2 = xx(i + 1, j, k);
                RDouble x3 = xx(i,     j, k + 1);
                RDouble x4 = xx(i + 1, j, k + 1);

                RDouble y1 = yy(i,     j, k);
                RDouble y2 = yy(i + 1, j, k);
                RDouble y3 = yy(i,     j, k + 1);
                RDouble y4 = yy(i + 1, j, k + 1);

                RDouble z1 = zz(i,     j, k);
                RDouble z2 = zz(i + 1, j, k);
                RDouble z3 = zz(i,     j, k + 1);
                RDouble z4 = zz(i + 1, j, k + 1);
                
                xfaceJ(i, j, k) = 0.25 * (x1 + x2 + x3 + x4);
                yfaceJ(i, j, k) = 0.25 * (y1 + y2 + y3 + y4);
                zfaceJ(i, j, k) = 0.25 * (z1 + z2 + z3 + z4);
            }
        }
    }

    for (int j = 1; j <= nj - 1; ++ j)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int k = -1; k <= nk + 2; ++ k)
            {
                RDouble x1 = xx(i,     j,     k);
                RDouble x2 = xx(i + 1, j,     k);
                RDouble x3 = xx(i,     j + 1, k);
                RDouble x4 = xx(i + 1, j + 1, k);

                RDouble y1 = yy(i,     j,     k);
                RDouble y2 = yy(i + 1, j,     k);
                RDouble y3 = yy(i,     j + 1, k);
                RDouble y4 = yy(i + 1, j + 1, k);

                RDouble z1 = zz(i,     j,     k);
                RDouble z2 = zz(i + 1, j,     k);
                RDouble z3 = zz(i,     j + 1, k);
                RDouble z4 = zz(i + 1, j + 1, k);
                
                xfaceK(i, j, k) = 0.25 * (x1 + x2 + x3 + x4);
                yfaceK(i, j, k) = 0.25 * (y1 + y2 + y3 + y4);
                zfaceK(i, j, k) = 0.25 * (z1 + z2 + z3 + z4);
            }
        }
    }

    int ZoneID = this->GetZoneID();
    std::ostringstream oss;
    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = 1; i <= ni - 1; ++ i)
            {
                RDouble vol_kxi, vol_eta, vol_cta;

                RDouble kxi_CL2 = half * (xfaceI(i - 1, j, k) + xfaceI(i - 2, j, k)) * Cellxfv(i - 2, j, k, 1) + half * (yfaceI(i - 1, j, k) + yfaceI(i - 2, j, k)) * Cellyfv(i - 2, j, k, 1) + half * (zfaceI(i - 1, j, k) + zfaceI(i - 2, j, k)) * Cellzfv(i - 2, j, k, 1);
                RDouble kxi_CL1 = half * (xfaceI(i    , j, k) + xfaceI(i - 1, j, k)) * Cellxfv(i - 1, j, k, 1) + half * (yfaceI(i    , j, k) + yfaceI(i - 1, j, k)) * Cellyfv(i - 1, j, k, 1) + half * (zfaceI(i    , j, k) + zfaceI(i - 1, j, k)) * Cellzfv(i - 1, j, k, 1);
                RDouble kxi_fL  = half * (xfaceI(i    , j, k) + xfaceI(i    , j, k)) *     xfv(i    , j, k, 1) + half * (yfaceI(i    , j, k) + yfaceI(i    , j, k)) *     yfv(i    , j, k, 1) + half * (zfaceI(i    , j, k) + zfaceI(i    , j, k)) *     zfv(i    , j, k, 1);
                RDouble kxi_fR  = half * (xfaceI(i + 1, j, k) + xfaceI(i + 1, j, k)) *     xfv(i + 1, j, k, 1) + half * (yfaceI(i + 1, j, k) + yfaceI(i + 1, j, k)) *     yfv(i + 1, j, k, 1) + half * (zfaceI(i + 1, j, k) + zfaceI(i + 1, j, k)) *     zfv(i + 1, j, k, 1);
                RDouble kxi_CR1 = half * (xfaceI(i + 1, j, k) + xfaceI(i + 2, j, k)) * Cellxfv(i + 1, j, k, 1) + half * (yfaceI(i + 1, j, k) + yfaceI(i + 2, j, k)) * Cellyfv(i + 1, j, k, 1) + half * (zfaceI(i + 1, j, k) + zfaceI(i + 2, j, k)) * Cellzfv(i + 1, j, k, 1);
                RDouble kxi_CR2 = half * (xfaceI(i + 2, j, k) + xfaceI(i + 3, j, k)) * Cellxfv(i + 2, j, k, 1) + half * (yfaceI(i + 2, j, k) + yfaceI(i + 3, j, k)) * Cellyfv(i + 2, j, k, 1) + half * (zfaceI(i + 2, j, k) + zfaceI(i + 3, j, k)) * Cellzfv(i + 2, j, k, 1);
                
                vol_kxi = alfa_sixth * (kxi_fR - kxi_fL) + a_sixth * (kxi_CR1 - kxi_CL1) + b_sixth * (kxi_CR2 - kxi_CL2);

                RDouble eta_CL2 = half * (xfaceJ(i, j - 1, k) + xfaceJ(i, j - 2, k)) * Cellxfv(i, j - 2, k, 2) + half * (yfaceJ(i, j - 1, k) + yfaceJ(i, j - 2, k)) * Cellyfv(i, j - 2, k, 2) + half * (zfaceJ(i, j - 1, k) + zfaceJ(i, j - 2, k)) * Cellzfv(i, j - 2, k, 2);
                RDouble eta_CL1 = half * (xfaceJ(i, j    , k) + xfaceJ(i, j - 1, k)) * Cellxfv(i, j - 1, k, 2) + half * (yfaceJ(i, j    , k) + yfaceJ(i, j - 1, k)) * Cellyfv(i, j - 1, k, 2) + half * (zfaceJ(i, j    , k) + zfaceJ(i, j - 1, k)) * Cellzfv(i, j - 1, k, 2);
                RDouble eta_fL  = half * (xfaceJ(i, j    , k) + xfaceJ(i, j    , k)) *     xfv(i, j    , k, 2) + half * (yfaceJ(i, j    , k) + yfaceJ(i, j    , k)) *     yfv(i, j    , k, 2) + half * (zfaceJ(i, j    , k) + zfaceJ(i, j    , k)) *     zfv(i, j    , k, 2);
                RDouble eta_fR  = half * (xfaceJ(i, j + 1, k) + xfaceJ(i, j + 1, k)) *     xfv(i, j + 1, k, 2) + half * (yfaceJ(i, j + 1, k) + yfaceJ(i, j + 1, k)) *     yfv(i, j + 1, k, 2) + half * (zfaceJ(i, j + 1, k) + zfaceJ(i, j + 1, k)) *     zfv(i, j + 1, k, 2);
                RDouble eta_CR1 = half * (xfaceJ(i, j + 1, k) + xfaceJ(i, j + 2, k)) * Cellxfv(i, j + 1, k, 2) + half * (yfaceJ(i, j + 1, k) + yfaceJ(i, j + 2, k)) * Cellyfv(i, j + 1, k, 2) + half * (zfaceJ(i, j + 1, k) + zfaceJ(i, j + 2, k)) * Cellzfv(i, j + 1, k, 2);
                RDouble eta_CR2 = half * (xfaceJ(i, j + 2, k) + xfaceJ(i, j + 3, k)) * Cellxfv(i, j + 2, k, 2) + half * (yfaceJ(i, j + 2, k) + yfaceJ(i, j + 3, k)) * Cellyfv(i, j + 2, k, 2) + half * (zfaceJ(i, j + 2, k) + zfaceJ(i, j + 3, k)) * Cellzfv(i, j + 2, k, 2);
                
                vol_eta = alfa_sixth * (eta_fR - eta_fL) + a_sixth * (eta_CR1 - eta_CL1) + b_sixth * (eta_CR2 - eta_CL2);
                
                RDouble cta_CL2 = half * (xfaceK(i, j, k - 1) + xfaceK(i, j, k - 2)) * Cellxfv(i, j, k - 2, 3) + half * (yfaceK(i, j, k - 1) + yfaceK(i, j, k - 2)) * Cellyfv(i, j, k - 2, 3) + half * (zfaceK(i, j, k - 1) + zfaceK(i, j, k - 2)) * Cellzfv(i, j, k - 2, 3);
                RDouble cta_CL1 = half * (xfaceK(i, j, k   ) + xfaceK(i, j, k - 1)) * Cellxfv(i, j, k - 1, 3) + half * (yfaceK(i, j, k   ) + yfaceK(i, j, k - 1)) * Cellyfv(i, j, k - 1, 3) + half * (zfaceK(i, j, k   ) + zfaceK(i, j, k - 1)) * Cellzfv(i, j, k - 1, 3);
                RDouble cta_fL  = half * (xfaceK(i, j, k   ) + xfaceK(i, j, k   )) *     xfv(i, j, k    , 3) + half * (yfaceK(i, j, k   ) + yfaceK(i, j, k   )) *     yfv(i, j, k    , 3) + half * (zfaceK(i, j, k   ) + zfaceK(i, j, k   )) *     zfv(i, j, k    , 3);
                RDouble cta_fR  = half * (xfaceK(i, j, k + 1) + xfaceK(i, j, k + 1)) *     xfv(i, j, k + 1, 3) + half * (yfaceK(i, j, k + 1) + yfaceK(i, j, k + 1)) *     yfv(i, j, k + 1, 3) + half * (zfaceK(i, j, k + 1) + zfaceK(i, j, k + 1)) *     zfv(i, j, k + 1, 3);
                RDouble cta_CR1 = half * (xfaceK(i, j, k + 1) + xfaceK(i, j, k + 2)) * Cellxfv(i, j, k + 1, 3) + half * (yfaceK(i, j, k + 1) + yfaceK(i, j, k + 2)) * Cellyfv(i, j, k + 1, 3) + half * (zfaceK(i, j, k + 1) + zfaceK(i, j, k + 2)) * Cellzfv(i, j, k + 1, 3);
                RDouble cta_CR2 = half * (xfaceK(i, j, k + 2) + xfaceK(i, j, k + 3)) * Cellxfv(i, j, k + 2, 3) + half * (yfaceK(i, j, k + 2) + yfaceK(i, j, k + 3)) * Cellyfv(i, j, k + 2, 3) + half * (zfaceK(i, j, k + 2) + zfaceK(i, j, k + 3)) * Cellzfv(i, j, k + 2, 3);
                
                vol_cta = alfa_sixth * (cta_fR - cta_fL) + a_sixth * (cta_CR1 - cta_CL1) + b_sixth * (cta_CR2 - cta_CL2);
                
                RDouble vol_FD = (vol_kxi + vol_eta + vol_cta) / 3.0;

                if (vol_FD < SMALL)
                {
                    oss << " vol < 0 " << " ZoneID = " << ZoneID << " i = " << i << " j = " << j << " vol = " << vol_FD << "\n";
                }
                //!RDouble vol_FV = vol(i, j, k);
                //!
                //!if (vol_FD <= zero || vol_FD > 1.2 * vol_FV || vol_FD < 0.8 * vol_FV)
                //!{
                //!    vol_FD = vol_FV;
                //!}
                vol(i ,j ,k) = vol_FD;
            }
        }
    }

    delete &xfaceI;
    delete &yfaceI;
    delete &zfaceI;

    delete &xfaceJ;
    delete &yfaceJ;
    delete &zfaceJ;

    delete &xfaceK;
    delete &yfaceK;
    delete &zfaceK;

    //!=
    ComputeCellCenter();
}

void StructGrid::GhostMetricsStructHighOrder2D()
{
    int ni = this->GetNI();
    int nj = this->GetNJ();

    RDouble3D &vol  = *(this->GetCellVolume());
    RDouble3D &xcc = *(this->GetCellCenterX());
    RDouble3D &ycc = *(this->GetCellCenterY());
    RDouble3D &zcc = *(this->GetCellCenterZ());

    const int k = 1;

    for (int layer = 1; layer <= 4; ++layer)
    {
        for (int j = 1; j <= nj - 1; ++j)
        {
            vol(    1 - layer, j, k) = vol(    layer, j, k);
            vol(ni - 1 + layer, j, k) = vol(ni - layer, j, k);

            int nsurf = 1;
            RDouble xfc, yfc, zfc;

            FaceCoor(1, j, k, nsurf, xfc, yfc, zfc);            
            xcc(1-layer,j,k) = two * xfc - xcc(layer,j,k);
            ycc(1-layer,j,k) = two * yfc - ycc(layer,j,k);
            zcc(1-layer,j,k) = two * zfc - zcc(layer,j,k);

            
            FaceCoor(ni, j, k, nsurf, xfc, yfc, zfc);            
            xcc(ni-1+layer,j,k) = two * xfc - xcc(ni-layer,j,k);
            ycc(ni-1+layer,j,k) = two * yfc - ycc(ni-layer,j,k);
            zcc(ni-1+layer,j,k) = two * zfc - zcc(ni-layer,j,k);
        }

        for (int i = 1; i <= ni - 1; ++i)
        {
            vol(i,      1 - layer, k) = vol(i,      layer, k);
            vol(i, nj - 1 + layer, k) = vol(i, nj - layer, k);

            int nsurf = 2;
            RDouble xfc, yfc, zfc;

            FaceCoor(i, 1, k, nsurf, xfc, yfc, zfc);            
            xcc(i,1-layer,k) = two * xfc - xcc(i,layer,k);
            ycc(i,1-layer,k) = two * yfc - ycc(i,layer,k);
            zcc(i,1-layer,k) = two * zfc - zcc(i,layer,k);

            
            FaceCoor(i, nj, k, nsurf, xfc, yfc, zfc);            
            xcc(i,nj-1+layer,k) = two * xfc - xcc(i,nj-layer,k);
            ycc(i,nj-1+layer,k) = two * yfc - ycc(i,nj-layer,k);
            zcc(i,nj-1+layer,k) = two * zfc - zcc(i,nj-layer,k);
        }
    }    
}

void StructGrid::GhostMetricsStructHighOrder3D()
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    RDouble3D &vol  = *(this->GetCellVolume());
    RDouble3D &xcc = *(this->GetCellCenterX());
    RDouble3D &ycc = *(this->GetCellCenterY());
    RDouble3D &zcc = *(this->GetCellCenterZ());
    
    for (int layer = 1; layer <= 4; ++layer)
    {
        for (int k = 1; k <= nk - 1; ++k)
        {
            for (int j = 1; j <= nj - 1; ++j)
            {
                vol(    1 - layer, j, k) = vol(    layer, j, k);
                vol(ni - 1 + layer, j, k) = vol(ni - layer, j, k);
                
                int nsurf = 1;
                RDouble xfc, yfc, zfc;
                
                FaceCoor(1, j, k, nsurf, xfc, yfc, zfc); 
                xcc(1-layer,j,k) = two * xfc - xcc(layer,j,k);
                ycc(1-layer,j,k) = two * yfc - ycc(layer,j,k);
                zcc(1-layer,j,k) = two * zfc - zcc(layer,j,k);
                
                FaceCoor(ni, j, k, nsurf, xfc, yfc, zfc);
                xcc(ni-1+layer,j,k) = two * xfc - xcc(ni-layer,j,k);
                ycc(ni-1+layer,j,k) = two * yfc - ycc(ni-layer,j,k);
                zcc(ni-1+layer,j,k) = two * zfc - zcc(ni-layer,j,k);
            }
        }
        
        for (int k = 1; k <= nk - 1; ++k)
        {
            for (int i = 1; i <= ni - 1; ++i)
            {
                vol(i,      1 - layer, k) = vol(i,      layer, k);
                vol(i, nj - 1 + layer, k) = vol(i, nj - layer, k);
                
                int nsurf = 2;
                RDouble xfc, yfc, zfc;
                
                FaceCoor(i, 1, k, nsurf, xfc, yfc, zfc);
                xcc(i,1-layer,k) = two * xfc - xcc(i,layer,k);
                ycc(i,1-layer,k) = two * yfc - ycc(i,layer,k);
                zcc(i,1-layer,k) = two * zfc - zcc(i,layer,k);
                
                FaceCoor(i, nj, k, nsurf, xfc, yfc, zfc);
                xcc(i,nj-1+layer,k) = two * xfc - xcc(i,nj-layer,k);
                ycc(i,nj-1+layer,k) = two * yfc - ycc(i,nj-layer,k);
                zcc(i,nj-1+layer,k) = two * zfc - zcc(i,nj-layer,k);
            }
        }

        for (int j = 1; j <= nj - 1; ++j)
        {
            for (int i = 1; i <= ni - 1; ++i)
            {
                vol(i, j,      1 - layer) = vol(i, j,      layer);
                vol(i, j, nk - 1 + layer) = vol(i, j, nk - layer);
                
                int nsurf = 3;
                RDouble xfc, yfc, zfc;
                
                FaceCoor(i, j, 1, nsurf, xfc, yfc, zfc);
                xcc(i,j,1-layer) = two * xfc - xcc(i,j,layer);
                ycc(i,j,1-layer) = two * yfc - ycc(i,j,layer);
                zcc(i,j,1-layer) = two * zfc - zcc(i,j,layer);
                
                FaceCoor(i, j, nk, nsurf, xfc, yfc, zfc);
                xcc(i,j,nk-1+layer) = two * xfc - xcc(i,j,nk-layer);
                ycc(i,j,nk-1+layer) = two * yfc - ycc(i,j,nk-layer);
                zcc(i,j,nk-1+layer) = two * zfc - zcc(i,j,nk-layer);
            }
        }
        
    }    
}

void StructGrid::InitialDiffOrderforStructHighOrder()
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    DiscretePrecisionKXI = new Int1D(Range(-1, ni + 1), fortranArray);
    DiscretePrecisionETA = new Int1D(Range(-1, nj + 1), fortranArray);
    DiscretePrecisionCTA = new Int1D(Range(-1, nk + 1), fortranArray);
    
    for (int i = -1; i <= ni + 1; ++ i)
    {
        (*DiscretePrecisionKXI)(i) = 4;
    }

    for (int j = -1; j <= nj + 1; ++ j)
    {
        (*DiscretePrecisionETA)(j) = 4;
    }

    for (int k = -1; k <= nk + 1; ++ k)
    {
        (*DiscretePrecisionCTA)(k) = 4;
    }

    const RDouble LowerValue = cos(15.0 / 180.0 * PI);

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    StructBCSet *structBCSet = this->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int BCType = bcregion->GetBCType();
        int nsurf_bc = bcregion->GetFaceDirection() + 1;
        if (IsInterface(BCType))
        {
            continue;
        }
        int ist, ied, jst, jed, kst, ked;
        bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);
        
        int *s_lr3d = bcregion->GetFaceDirectionIndex();  
        int id, jd, kd;
        GetBCFaceIDX(s_lr3d, id, jd, kd);
        if (nsurf_bc == 1)
        {
            if (id + jd + kd == 0)
            {
                (*DiscretePrecisionKXI)(-1) = 2;
                (*DiscretePrecisionKXI)(0) = 2;
                (*DiscretePrecisionKXI)(1) = 2;
            }
            else
            {
                (*DiscretePrecisionKXI)(ni - 1) = 2;
                (*DiscretePrecisionKXI)(ni   ) = 2;
                (*DiscretePrecisionKXI)(ni + 1) = 2;
            } 
        }
        else if (nsurf_bc == 2)
        {
            if (id + jd + kd == 0)
            {
                (*DiscretePrecisionETA)(-1) = 2;
                (*DiscretePrecisionETA)(0) = 2;
                (*DiscretePrecisionETA)(1) = 2;
            }
            else
            {
                (*DiscretePrecisionETA)(nj - 1) = 2;
                (*DiscretePrecisionETA)(nj   ) = 2;
                (*DiscretePrecisionETA)(nj + 1) = 2;
            } 
        }
        else
        {
            if (id + jd + kd == 0)
            {
                (*DiscretePrecisionCTA)(-1) = 2;
                (*DiscretePrecisionCTA)(0) = 2;
                (*DiscretePrecisionCTA)(1) = 2;
            }
            else
            {
                (*DiscretePrecisionCTA)(nk - 1) = 2;
                (*DiscretePrecisionCTA)(nk   ) = 2;
                (*DiscretePrecisionCTA)(nk + 1) = 2;
            } 
        }          
    }    

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni;
    int jed = nj;
    int ked = nk;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble ax = - xx(i - 1, j, k) + xx(i    , j, k);
                RDouble ay = - yy(i - 1, j, k) + yy(i    , j, k);
                RDouble az = - zz(i - 1, j, k) + zz(i    , j, k);
                RDouble bx = - xx(i    , j, k) + xx(i + 1, j, k);
                RDouble by = - yy(i    , j, k) + yy(i + 1, j, k);
                RDouble bz = - zz(i    , j, k) + zz(i + 1, j, k);

                RDouble la = sqrt(ax * ax + ay * ay + az * az);
                RDouble lb = sqrt(bx * bx + by * by + bz * bz);                
                RDouble ab = ax * bx + ay * by + az * bz;
                
                if (ab < la * lb * LowerValue)
                {
                    (*DiscretePrecisionKXI)(i - 1) = 2;
                    (*DiscretePrecisionKXI)(i   ) = 2;
                }
            }            
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble ax = - xx(i, j - 1, k) + xx(i, j    , k);
                RDouble ay = - yy(i, j - 1, k) + yy(i, j    , k);
                RDouble az = - zz(i, j - 1, k) + zz(i, j    , k);
                RDouble bx = - xx(i, j    , k) + xx(i, j + 1, k);
                RDouble by = - yy(i, j    , k) + yy(i, j + 1, k);
                RDouble bz = - zz(i, j    , k) + zz(i, j + 1, k);

                RDouble la = sqrt(ax * ax + ay * ay + az * az);
                RDouble lb = sqrt(bx * bx + by * by + bz * bz);                
                RDouble ab = ax * bx + ay * by + az * bz;
                
                if (ab < la * lb * LowerValue)
                {
                    (*DiscretePrecisionETA)(j - 1) = 2;
                    (*DiscretePrecisionETA)(j   ) = 2;
                }
            }            
        }
    }
    
    if (nk == 1)
    {
        return;
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble ax = - xx(i, j, k - 1) + xx(i, j, k   );
                RDouble ay = - yy(i, j, k - 1) + yy(i, j, k   );
                RDouble az = - zz(i, j, k - 1) + zz(i, j, k   );
                RDouble bx = - xx(i, j, k   ) + xx(i, j, k + 1);
                RDouble by = - yy(i, j, k   ) + yy(i, j, k + 1);
                RDouble bz = - zz(i, j, k   ) + zz(i, j, k + 1);

                RDouble la = sqrt(ax * ax + ay * ay + az * az);
                RDouble lb = sqrt(bx * bx + by * by + bz * bz);                
                RDouble ab = ax * bx + ay * by + az * bz;
                
                if (ab < la * lb * LowerValue)
                {
                    (*DiscretePrecisionCTA)(k - 1) = 2;
                    (*DiscretePrecisionCTA)(k   ) = 2;
                }
            }            
        }
    }    
}

void StructGrid::CorrectDiffOrderNearUnstructBlockforStructHighOrder(const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);
        
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];

        for (int layer = 1; layer <= 1; ++ layer)
        {            
            int id, jd, kd;
            int i_lr, j_lr, k_lr;
            finestGrid->GetTargetIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);
            if (i_lr != 0)
            {
                //! I_face
                if (i_lr == -1)
                {
                    //! left boundary     
                    (*DiscretePrecisionKXI)(-1) = 2;
                    (*DiscretePrecisionKXI)(0) = 2;
                    (*DiscretePrecisionKXI)(1) = 2;
                }
                else
                {
                    //! right boundary
                    (*DiscretePrecisionKXI)(ni - 1) = 2;
                    (*DiscretePrecisionKXI)(ni   ) = 2;
                    (*DiscretePrecisionKXI)(ni + 1) = 2;
                }
            }
            else if (j_lr != 0)
            {
                //! J_face
                if (j_lr == -1)
                {
                    //! left boundary                    
                    (*DiscretePrecisionETA)(-1) = 2;
                    (*DiscretePrecisionETA)(0) = 2;
                    (*DiscretePrecisionETA)(1) = 2;
                }
                else
                {
                    //! right boundary
                    (*DiscretePrecisionETA)(nj - 1) = 2;
                    (*DiscretePrecisionETA)(nj   ) = 2;
                    (*DiscretePrecisionETA)(nj + 1) = 2;
                }
            }
            else
            {
                //! K_face
                if (k_lr == -1)
                {
                    //! left boundary                    
                    (*DiscretePrecisionCTA)(-1) = 2;
                    (*DiscretePrecisionCTA)(0) = 2;
                    (*DiscretePrecisionCTA)(1) = 2;
                }
                else
                {
                    //! right boundary
                    (*DiscretePrecisionCTA)(nk - 1) = 2;
                    (*DiscretePrecisionCTA)(nk   ) = 2;
                    (*DiscretePrecisionCTA)(nk + 1) = 2;
                }
            }  
        }
    }  
}

void StructGrid::ComputeMetricsWithBCforStructHighOrder()
{
    if (GetDim() == 2)
    {
        ComputeMetricsWithBCforStructHighOrder2D();
    }
    else
    {
        ComputeMetricsWithBCforStructHighOrder3D();
    }
}

void StructGrid::ComputeMetricsWithBCforStructHighOrder2D()
{
    RDouble4D &xfn  = *(this->GetFaceNormalX());
    RDouble4D &yfn  = *(this->GetFaceNormalY());
    RDouble4D &zfn  = *(this->GetFaceNormalZ());
    RDouble4D &area = *(this->GetFaceArea());
    RDouble4D &xfv  = *(this->GetFaceVectorX());
    RDouble4D &yfv  = *(this->GetFaceVectorY());
    RDouble4D &zfv  = *(this->GetFaceVectorZ());

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();

    zfn = 0.0;
    zfv = 0.0;

    int ni = this->GetNI();
    int nj = this->GetNJ();
    const int k = 1;
    
    for (int i = 1; i <= ni; ++ i)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            RDouble xeta, yeta;

            RDouble x0 = xx(i, j - 1, k);
            RDouble x1 = xx(i, j,     k);
            RDouble x2 = xx(i, j + 1, k);
            RDouble x3 = xx(i, j + 2, k);

            RDouble y0 = yy(i, j - 1, k);
            RDouble y1 = yy(i, j,     k);
            RDouble y2 = yy(i, j + 1, k);
            RDouble y3 = yy(i, j + 2, k);

            if ((*DiscretePrecisionETA)(j) == 2)
            {
                xeta = (x2 - x1);
                yeta = (y2 - y1);
            }
            else
            {
                xeta = (27.0 * (x2 - x1) - (x3 - x0)) / 24.0;
                yeta = (27.0 * (y2 - y1) - (y3 - y0)) / 24.0;
            }
            
            RDouble ds = sqrt(xeta * xeta + yeta * yeta);

            xfv(i, j, k, 1) =   yeta;
            yfv(i, j, k, 1) = - xeta;

            area(i, j, k, 1) = ds;

            xfn(i, j, k, 1) =   yeta / (ds + TINY);
            yfn(i, j, k, 1) = - xeta / (ds + TINY);
        }
    }

    for (int j = 1; j <= nj; ++ j)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            RDouble xkxi, ykxi;

            RDouble x0 = xx(i - 1, j, k);
            RDouble x1 = xx(i,     j, k);
            RDouble x2 = xx(i + 1, j, k);
            RDouble x3 = xx(i + 2, j, k);

            RDouble y0 = yy(i - 1, j, k);
            RDouble y1 = yy(i,     j, k);
            RDouble y2 = yy(i + 1, j, k);
            RDouble y3 = yy(i + 2, j, k);

            if ((*DiscretePrecisionKXI)(i) == 2)
            {
                xkxi = (x2 - x1);
                ykxi = (y2 - y1);
            }
            else
            {
                xkxi = (27.0 * (x2 - x1) - (x3 - x0)) / 24.0;
                ykxi = (27.0 * (y2 - y1) - (y3 - y0)) / 24.0;
            }
            
            RDouble ds = sqrt(xkxi * xkxi + ykxi * ykxi);

            xfv(i, j, k, 2) = - ykxi;
            yfv(i, j, k, 2) =   xkxi;

            area(i, j, k, 2) = ds;

            xfn(i, j, k, 2) = - ykxi / (ds + TINY);
            yfn(i, j, k, 2) =   xkxi / (ds + TINY);
        }        
    }

    //! BC metrics
    //! interior BC would be replaced in the next function CommPartofMetrics_FD.
    GhostMetrics2D();
}

void StructGrid::ComputeMetricsWithBCforStructHighOrder3D()
{
    RDouble4D &xfn  = *(this->GetFaceNormalX());
    RDouble4D &yfn  = *(this->GetFaceNormalY());
    RDouble4D &zfn  = *(this->GetFaceNormalZ());
    RDouble4D &area = *(this->GetFaceArea());
    RDouble4D &xfv  = *(this->GetFaceVectorX());
    RDouble4D &yfv  = *(this->GetFaceVectorY());
    RDouble4D &zfv  = *(this->GetFaceVectorZ());

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();  
    RDouble3D &xA    = *new RDouble3D(Range(1, ni - 1), Range(0, nj + 1), Range(0, nk + 1), fortranArray);
    RDouble3D &yA    = *new RDouble3D(Range(1, ni - 1), Range(0, nj + 1), Range(0, nk + 1), fortranArray);
    RDouble3D &zA    = *new RDouble3D(Range(1, ni - 1), Range(0, nj + 1), Range(0, nk + 1), fortranArray);
    RDouble3D &xkxi = *new RDouble3D(Range(1, ni - 1), Range(0, nj + 1), Range(0, nk + 1), fortranArray);
    RDouble3D &ykxi = *new RDouble3D(Range(1, ni - 1), Range(0, nj + 1), Range(0, nk + 1), fortranArray);
    RDouble3D &zkxi = *new RDouble3D(Range(1, ni - 1), Range(0, nj + 1), Range(0, nk + 1), fortranArray);

    RDouble3D &xB    = *new RDouble3D(Range(0, ni + 1), Range(1, nj - 1), Range(0, nk + 1), fortranArray);
    RDouble3D &yB    = *new RDouble3D(Range(0, ni + 1), Range(1, nj - 1), Range(0, nk + 1), fortranArray);
    RDouble3D &zB    = *new RDouble3D(Range(0, ni + 1), Range(1, nj - 1), Range(0, nk + 1), fortranArray);
    RDouble3D &xeta = *new RDouble3D(Range(0, ni + 1), Range(1, nj - 1), Range(0, nk + 1), fortranArray);
    RDouble3D &yeta = *new RDouble3D(Range(0, ni + 1), Range(1, nj - 1), Range(0, nk + 1), fortranArray);
    RDouble3D &zeta = *new RDouble3D(Range(0, ni + 1), Range(1, nj - 1), Range(0, nk + 1), fortranArray);

    RDouble3D &xC    = *new RDouble3D(Range(0, ni + 1), Range(0, nj + 1), Range(1, nk - 1), fortranArray);
    RDouble3D &yC    = *new RDouble3D(Range(0, ni + 1), Range(0, nj + 1), Range(1, nk - 1), fortranArray);
    RDouble3D &zC    = *new RDouble3D(Range(0, ni + 1), Range(0, nj + 1), Range(1, nk - 1), fortranArray);
    RDouble3D &xcta = *new RDouble3D(Range(0, ni + 1), Range(0, nj + 1), Range(1, nk - 1), fortranArray);
    RDouble3D &ycta = *new RDouble3D(Range(0, ni + 1), Range(0, nj + 1), Range(1, nk - 1), fortranArray);
    RDouble3D &zcta = *new RDouble3D(Range(0, ni + 1), Range(0, nj + 1), Range(1, nk - 1), fortranArray);
    
    int ibc[2];
    int jbc[2];
    int kbc[2];

    ibc[0] = 0;
    jbc[0] = 0;
    kbc[0] = 0;
    ibc[1] = ni + 1;
    jbc[1] = nj + 1;
    kbc[1] = nk + 1;   

    for (int i = 1; i <= ni - 1; ++ i)
    {
        for (int j = 1; j <= nj; ++ j)
        {
            for (int k = 1; k <= nk; ++ k)
            {
                if ((*DiscretePrecisionKXI)(i) == 2)
                {
                    xA(i, j, k) = half * (xx(i + 1, j, k) + xx(i, j, k));
                    yA(i, j, k) = half * (yy(i + 1, j, k) + yy(i, j, k));
                    zA(i, j, k) = half * (zz(i + 1, j, k) + zz(i, j, k));
                    
                    xkxi(i, j, k) = xx(i + 1, j, k) - xx(i, j, k);
                    ykxi(i, j, k) = yy(i + 1, j, k) - yy(i, j, k);
                    zkxi(i, j, k) = zz(i + 1, j, k) - zz(i, j, k);
                }
                else
                {
                    xA(i, j, k) = half * (xx(i + 1, j, k) + xx(i, j, k));
                    yA(i, j, k) = half * (yy(i + 1, j, k) + yy(i, j, k));
                    zA(i, j, k) = half * (zz(i + 1, j, k) + zz(i, j, k));
                    
                    xkxi(i, j, k) = (27.0 * (xx(i + 1, j, k) - xx(i, j, k)) - (xx(i + 2, j, k) - xx(i - 1, j, k))) / 24.0;
                    ykxi(i, j, k) = (27.0 * (yy(i + 1, j, k) - yy(i, j, k)) - (yy(i + 2, j, k) - yy(i - 1, j, k))) / 24.0;
                    zkxi(i, j, k) = (27.0 * (zz(i + 1, j, k) - zz(i, j, k)) - (zz(i + 2, j, k) - zz(i - 1, j, k))) / 24.0;
                }                
            }
        }
        for (int bcindex = 0; bcindex <= 1; ++ bcindex)
        {
            int j = jbc[bcindex];
            for (int k = 1; k <= nk; ++ k)
            {
                if ((*DiscretePrecisionKXI)(i) == 2)
                {
                    xA(i, j, k) = half * (xx(i + 1, j, k) + xx(i, j, k));
                    yA(i, j, k) = half * (yy(i + 1, j, k) + yy(i, j, k));
                    zA(i, j, k) = half * (zz(i + 1, j, k) + zz(i, j, k));
                    
                    xkxi(i, j, k) = xx(i + 1, j, k) - xx(i, j, k);
                    ykxi(i, j, k) = yy(i + 1, j, k) - yy(i, j, k);
                    zkxi(i, j, k) = zz(i + 1, j, k) - zz(i, j, k);
                }
                else
                {
                    xA(i, j, k) = half * (xx(i + 1, j, k) + xx(i, j, k));
                    yA(i, j, k) = half * (yy(i + 1, j, k) + yy(i, j, k));
                    zA(i, j, k) = half * (zz(i + 1, j, k) + zz(i, j, k));
                    
                    xkxi(i, j, k) = (27.0 * (xx(i + 1, j, k) - xx(i, j, k)) - (xx(i + 2, j, k) - xx(i - 1, j, k))) / 24.0;
                    ykxi(i, j, k) = (27.0 * (yy(i + 1, j, k) - yy(i, j, k)) - (yy(i + 2, j, k) - yy(i - 1, j, k))) / 24.0;
                    zkxi(i, j, k) = (27.0 * (zz(i + 1, j, k) - zz(i, j, k)) - (zz(i + 2, j, k) - zz(i - 1, j, k))) / 24.0;
                }  
            }
        }
        for (int bcindex = 0; bcindex <= 1; ++bcindex)
        {
            int k = kbc[bcindex];
            for (int j = 1; j <= nj; ++ j)
            {
                if ((*DiscretePrecisionKXI)(i) == 2)
                {
                    xA(i, j, k) = half * (xx(i + 1, j, k) + xx(i, j, k));
                    yA(i, j, k) = half * (yy(i + 1, j, k) + yy(i, j, k));
                    zA(i, j, k) = half * (zz(i + 1, j, k) + zz(i, j, k));
                    
                    xkxi(i, j, k) = xx(i + 1, j, k) - xx(i, j, k);
                    ykxi(i, j, k) = yy(i + 1, j, k) - yy(i, j, k);
                    zkxi(i, j, k) = zz(i + 1, j, k) - zz(i, j, k);
                }
                else
                {
                    xA(i, j, k) = half * (xx(i + 1, j, k) + xx(i, j, k));
                    yA(i, j, k) = half * (yy(i + 1, j, k) + yy(i, j, k));
                    zA(i, j, k) = half * (zz(i + 1, j, k) + zz(i, j, k));
                    
                    xkxi(i, j, k) = (27.0 * (xx(i + 1, j, k) - xx(i, j, k)) - (xx(i + 2, j, k) - xx(i - 1, j, k))) / 24.0;
                    ykxi(i, j, k) = (27.0 * (yy(i + 1, j, k) - yy(i, j, k)) - (yy(i + 2, j, k) - yy(i - 1, j, k))) / 24.0;
                    zkxi(i, j, k) = (27.0 * (zz(i + 1, j, k) - zz(i, j, k)) - (zz(i + 2, j, k) - zz(i - 1, j, k))) / 24.0;
                }  
            }
        }
    }

    for (int j = 1; j <= nj - 1; ++ j)
    {
        for (int k = 1; k <= nk; ++ k)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                if ((*DiscretePrecisionETA)(j) == 2)
                {
                    xB(i, j, k) = half * (xx(i, j + 1, k) + xx(i, j, k));
                    yB(i, j, k) = half * (yy(i, j + 1, k) + yy(i, j, k));
                    zB(i, j, k) = half * (zz(i, j + 1, k) + zz(i, j, k));
                    
                    xeta(i, j, k) = xx(i, j + 1, k) - xx(i, j, k);
                    yeta(i, j, k) = yy(i, j + 1, k) - yy(i, j, k);
                    zeta(i, j, k) = zz(i, j + 1, k) - zz(i, j, k);
                }
                else
                {
                    xB(i, j, k) = half * (xx(i, j + 1, k) + xx(i, j, k));
                    yB(i, j, k) = half * (yy(i, j + 1, k) + yy(i, j, k));
                    zB(i, j, k) = half * (zz(i, j + 1, k) + zz(i, j, k));
                    
                    xeta(i, j, k) = (27.0 * (xx(i, j + 1, k) - xx(i, j, k)) - (xx(i, j + 2, k) - xx(i, j - 1, k))) / 24.0;
                    yeta(i, j, k) = (27.0 * (yy(i, j + 1, k) - yy(i, j, k)) - (yy(i, j + 2, k) - yy(i, j - 1, k))) / 24.0;
                    zeta(i, j, k) = (27.0 * (zz(i, j + 1, k) - zz(i, j, k)) - (zz(i, j + 2, k) - zz(i, j - 1, k))) / 24.0;
                }                
            }
        }
        for (int bcindex = 0; bcindex <= 1; ++bcindex)
        {
            int k = kbc[bcindex];
            for (int i = 1; i <= ni; ++ i)
            {
                if ((*DiscretePrecisionETA)(j) == 2)
                {
                    xB(i, j, k) = half * (xx(i, j + 1, k) + xx(i, j, k));
                    yB(i, j, k) = half * (yy(i, j + 1, k) + yy(i, j, k));
                    zB(i, j, k) = half * (zz(i, j + 1, k) + zz(i, j, k));
                    
                    xeta(i, j, k) = xx(i, j + 1, k) - xx(i, j, k);
                    yeta(i, j, k) = yy(i, j + 1, k) - yy(i, j, k);
                    zeta(i, j, k) = zz(i, j + 1, k) - zz(i, j, k);
                }
                else
                {
                    xB(i, j, k) = half * (xx(i, j + 1, k) + xx(i, j, k));
                    yB(i, j, k) = half * (yy(i, j + 1, k) + yy(i, j, k));
                    zB(i, j, k) = half * (zz(i, j + 1, k) + zz(i, j, k));
                    
                    xeta(i, j, k) = (27.0 * (xx(i, j + 1, k) - xx(i, j, k)) - (xx(i, j + 2, k) - xx(i, j - 1, k))) / 24.0;
                    yeta(i, j, k) = (27.0 * (yy(i, j + 1, k) - yy(i, j, k)) - (yy(i, j + 2, k) - yy(i, j - 1, k))) / 24.0;
                    zeta(i, j, k) = (27.0 * (zz(i, j + 1, k) - zz(i, j, k)) - (zz(i, j + 2, k) - zz(i, j - 1, k))) / 24.0;
                } 
            }
        }        
        for (int bcindex = 0; bcindex <= 1; ++bcindex)
        {
            int i = ibc[bcindex];
            for (int k = 1; k <= nk; ++ k)
            {
                if ((*DiscretePrecisionETA)(j) == 2)
                {
                    xB(i, j, k) = half * (xx(i, j + 1, k) + xx(i, j, k));
                    yB(i, j, k) = half * (yy(i, j + 1, k) + yy(i, j, k));
                    zB(i, j, k) = half * (zz(i, j + 1, k) + zz(i, j, k));
                    
                    xeta(i, j, k) = xx(i, j + 1, k) - xx(i, j, k);
                    yeta(i, j, k) = yy(i, j + 1, k) - yy(i, j, k);
                    zeta(i, j, k) = zz(i, j + 1, k) - zz(i, j, k);
                }
                else
                {
                    xB(i, j, k) = half * (xx(i, j + 1, k) + xx(i, j, k));
                    yB(i, j, k) = half * (yy(i, j + 1, k) + yy(i, j, k));
                    zB(i, j, k) = half * (zz(i, j + 1, k) + zz(i, j, k));
                    
                    xeta(i, j, k) = (27.0 * (xx(i, j + 1, k) - xx(i, j, k)) - (xx(i, j + 2, k) - xx(i, j - 1, k))) / 24.0;
                    yeta(i, j, k) = (27.0 * (yy(i, j + 1, k) - yy(i, j, k)) - (yy(i, j + 2, k) - yy(i, j - 1, k))) / 24.0;
                    zeta(i, j, k) = (27.0 * (zz(i, j + 1, k) - zz(i, j, k)) - (zz(i, j + 2, k) - zz(i, j - 1, k))) / 24.0;
                } 
            }
        }
    }
    
    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int i = 1; i <= ni; ++ i)
        {
            for (int j = 1; j <= nj; ++ j)
            {
                if ((*DiscretePrecisionCTA)(k) == 2)
                {
                    xC(i, j, k) = half * (xx(i, j, k + 1) + xx(i, j, k));
                    yC(i, j, k) = half * (yy(i, j, k + 1) + yy(i, j, k));
                    zC(i, j, k) = half * (zz(i, j, k + 1) + zz(i, j, k));
                    
                    xcta(i, j, k) = xx(i, j, k + 1) - xx(i, j, k);
                    ycta(i, j, k) = yy(i, j, k + 1) - yy(i, j, k);
                    zcta(i, j, k) = zz(i, j, k + 1) - zz(i, j, k);
                }
                else
                {
                    xC(i, j, k) = half * (xx(i, j, k + 1) + xx(i, j, k));
                    yC(i, j, k) = half * (yy(i, j, k + 1) + yy(i, j, k));
                    zC(i, j, k) = half * (zz(i, j, k + 1) + zz(i, j, k));
                    
                    xcta(i, j, k) = (27.0 * (xx(i, j, k + 1) - xx(i, j, k)) - (xx(i, j, k + 2) - xx(i, j, k - 1))) / 24.0;
                    ycta(i, j, k) = (27.0 * (yy(i, j, k + 1) - yy(i, j, k)) - (yy(i, j, k + 2) - yy(i, j, k - 1))) / 24.0;
                    zcta(i, j, k) = (27.0 * (zz(i, j, k + 1) - zz(i, j, k)) - (zz(i, j, k + 2) - zz(i, j, k - 1))) / 24.0;
                }                
            }
        }
        for (int bcindex = 0; bcindex <= 1; ++bcindex)
        {
            int i = ibc[bcindex];
            for (int j = 1; j <= nj; ++ j)
            {
                if ((*DiscretePrecisionCTA)(k) == 2)
                {
                    xC(i, j, k) = half * (xx(i, j, k + 1) + xx(i, j, k));
                    yC(i, j, k) = half * (yy(i, j, k + 1) + yy(i, j, k));
                    zC(i, j, k) = half * (zz(i, j, k + 1) + zz(i, j, k));
                    
                    xcta(i, j, k) = xx(i, j, k + 1) - xx(i, j, k);
                    ycta(i, j, k) = yy(i, j, k + 1) - yy(i, j, k);
                    zcta(i, j, k) = zz(i, j, k + 1) - zz(i, j, k);
                }
                else
                {
                    xC(i, j, k) = half * (xx(i, j, k + 1) + xx(i, j, k));
                    yC(i, j, k) = half * (yy(i, j, k + 1) + yy(i, j, k));
                    zC(i, j, k) = half * (zz(i, j, k + 1) + zz(i, j, k));
                    
                    xcta(i, j, k) = (27.0 * (xx(i, j, k + 1) - xx(i, j, k)) - (xx(i, j, k + 2) - xx(i, j, k - 1))) / 24.0;
                    ycta(i, j, k) = (27.0 * (yy(i, j, k + 1) - yy(i, j, k)) - (yy(i, j, k + 2) - yy(i, j, k - 1))) / 24.0;
                    zcta(i, j, k) = (27.0 * (zz(i, j, k + 1) - zz(i, j, k)) - (zz(i, j, k + 2) - zz(i, j, k - 1))) / 24.0;
                } 
            }
        }
        for (int bcindex = 0; bcindex <= 1; ++bcindex)
        {
            int j = jbc[bcindex];
            for (int i = 1; i <= ni; ++ i)
            {
                if ((*DiscretePrecisionCTA)(k) == 2)
                {
                    xC(i, j, k) = half * (xx(i, j, k + 1) + xx(i, j, k));
                    yC(i, j, k) = half * (yy(i, j, k + 1) + yy(i, j, k));
                    zC(i, j, k) = half * (zz(i, j, k + 1) + zz(i, j, k));
                    
                    xcta(i, j, k) = xx(i, j, k + 1) - xx(i, j, k);
                    ycta(i, j, k) = yy(i, j, k + 1) - yy(i, j, k);
                    zcta(i, j, k) = zz(i, j, k + 1) - zz(i, j, k);
                }
                else
                {
                    xC(i, j, k) = half * (xx(i, j, k + 1) + xx(i, j, k));
                    yC(i, j, k) = half * (yy(i, j, k + 1) + yy(i, j, k));
                    zC(i, j, k) = half * (zz(i, j, k + 1) + zz(i, j, k));
                    
                    xcta(i, j, k) = (27.0 * (xx(i, j, k + 1) - xx(i, j, k)) - (xx(i, j, k + 2) - xx(i, j, k - 1))) / 24.0;
                    ycta(i, j, k) = (27.0 * (yy(i, j, k + 1) - yy(i, j, k)) - (yy(i, j, k + 2) - yy(i, j, k - 1))) / 24.0;
                    zcta(i, j, k) = (27.0 * (zz(i, j, k + 1) - zz(i, j, k)) - (zz(i, j, k + 2) - zz(i, j, k - 1))) / 24.0;
                } 
            }
        }
    }

    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                RDouble sx0, sy0, sz0;
                RDouble sx1, sy1, sz1;
                RDouble sx2, sy2, sz2;
                RDouble sx3, sy3, sz3;

                RDouble KX = 0.0;
                RDouble KY = 0.0;
                RDouble KZ = 0.0;

                sx0 = half * (yeta(i, j, k - 1) * zB(i, j, k - 1) - yB(i, j, k - 1) * zeta(i, j, k - 1));
                sy0 = half * (zeta(i, j, k - 1) * xB(i, j, k - 1) - zB(i, j, k - 1) * xeta(i, j, k - 1));
                sz0 = half * (xeta(i, j, k - 1) * yB(i, j, k - 1) - xB(i, j, k - 1) * yeta(i, j, k - 1));
                sx1 = half * (yeta(i, j, k) * zB(i, j, k) - yB(i, j, k) * zeta(i, j, k));
                sy1 = half * (zeta(i, j, k) * xB(i, j, k) - zB(i, j, k) * xeta(i, j, k));
                sz1 = half * (xeta(i, j, k) * yB(i, j, k) - xB(i, j, k) * yeta(i, j, k));
                sx2 = half * (yeta(i, j, k + 1) * zB(i, j, k + 1) - yB(i, j, k + 1) * zeta(i, j, k + 1));
                sy2 = half * (zeta(i, j, k + 1) * xB(i, j, k + 1) - zB(i, j, k + 1) * xeta(i, j, k + 1));
                sz2 = half * (xeta(i, j, k + 1) * yB(i, j, k + 1) - xB(i, j, k + 1) * yeta(i, j, k + 1));
                sx3 = half * (yeta(i, j, k + 2) * zB(i, j, k + 2) - yB(i, j, k + 2) * zeta(i, j, k + 2));
                sy3 = half * (zeta(i, j, k + 2) * xB(i, j, k + 2) - zB(i, j, k + 2) * xeta(i, j, k + 2));
                sz3 = half * (xeta(i, j, k + 2) * yB(i, j, k + 2) - xB(i, j, k + 2) * yeta(i, j, k + 2));                

                if ((*DiscretePrecisionCTA)(k) == 2)
                {
                    KX += (sx2 - sx1);
                    KY += (sy2 - sy1);
                    KZ += (sz2 - sz1);
                }
                else
                {
                    KX += (27.0 * (sx2 - sx1) - (sx3 - sx0)) / 24.0;
                    KY += (27.0 * (sy2 - sy1) - (sy3 - sy0)) / 24.0;
                    KZ += (27.0 * (sz2 - sz1) - (sz3 - sz0)) / 24.0;
                } 
                
                sx0 = half * (yC(i, j - 1, k) * zcta(i, j - 1, k) - ycta(i, j - 1, k) * zC(i, j - 1, k));
                sy0 = half * (zC(i, j - 1, k) * xcta(i, j - 1, k) - zcta(i, j - 1, k) * xC(i, j - 1, k));
                sz0 = half * (xC(i, j - 1, k) * ycta(i, j - 1, k) - xcta(i, j - 1, k) * yC(i, j - 1, k));
                sx1 = half * (yC(i, j,  k) * zcta(i, j,  k) - ycta(i, j, k) * zC(i, j, k));
                sy1 = half * (zC(i, j,  k) * xcta(i, j,  k) - zcta(i, j, k) * xC(i, j, k));
                sz1 = half * (xC(i, j,  k) * ycta(i, j,  k) - xcta(i, j, k) * yC(i, j, k));
                sx2 = half * (yC(i, j + 1, k) * zcta(i, j + 1, k) - ycta(i, j + 1, k) * zC(i, j + 1, k));
                sy2 = half * (zC(i, j + 1, k) * xcta(i, j + 1, k) - zcta(i, j + 1, k) * xC(i, j + 1, k));
                sz2 = half * (xC(i, j + 1, k) * ycta(i, j + 1, k) - xcta(i, j + 1, k) * yC(i, j + 1, k));
                sx3 = half * (yC(i, j + 2, k) * zcta(i, j + 2, k) - ycta(i, j + 2, k) * zC(i, j + 2, k));
                sy3 = half * (zC(i, j + 2, k) * xcta(i, j + 2, k) - zcta(i, j + 2, k) * xC(i, j + 2, k));
                sz3 = half * (xC(i, j + 2, k) * ycta(i, j + 2, k) - xcta(i, j + 2, k) * yC(i, j + 2, k));                

                if ((*DiscretePrecisionETA)(j) == 2)
                {
                    KX += (sx2 - sx1);
                    KY += (sy2 - sy1);
                    KZ += (sz2 - sz1);
                }
                else
                {
                    KX += (27.0 * (sx2 - sx1) - (sx3 - sx0)) / 24.0;
                    KY += (27.0 * (sy2 - sy1) - (sy3 - sy0)) / 24.0;
                    KZ += (27.0 * (sz2 - sz1) - (sz3 - sz0)) / 24.0;
                }
                
                RDouble ds = sqrt(KX * KX + KY * KY + KZ * KZ);

                xfv (i, j, k, 1) = KX;
                yfv (i, j, k, 1) = KY;
                zfv (i, j, k, 1) = KZ;
                area(i, j, k, 1) = ds;
                xfn (i, j, k, 1) = KX / (ds + TINY);
                yfn (i, j, k, 1) = KY / (ds + TINY);
                zfn (i, j, k, 1) = KZ / (ds + TINY);
            }
        }
    }

    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int j = 1; j <= nj; ++ j)
            {
                RDouble sx0, sy0, sz0;
                RDouble sx1, sy1, sz1;
                RDouble sx2, sy2, sz2;
                RDouble sx3, sy3, sz3;

                RDouble KX = 0.0;
                RDouble KY = 0.0;
                RDouble KZ = 0.0;

                sx0 = half * (ycta(i - 1, j, k) * zC(i - 1, j, k) - yC(i - 1, j, k) * zcta(i - 1, j, k));
                sy0 = half * (zcta(i - 1, j, k) * xC(i - 1, j, k) - zC(i - 1, j, k) * xcta(i - 1, j, k));
                sz0 = half * (xcta(i - 1, j, k) * yC(i - 1, j, k) - xC(i - 1, j, k) * ycta(i - 1, j, k));
                sx1 = half * (ycta(i  , j, k) * zC(i  , j, k) - yC(i  , j, k) * zcta(i  , j, k));
                sy1 = half * (zcta(i  , j, k) * xC(i  , j, k) - zC(i  , j, k) * xcta(i  , j, k));
                sz1 = half * (xcta(i  , j, k) * yC(i  , j, k) - xC(i  , j, k) * ycta(i  , j, k));
                sx2 = half * (ycta(i + 1, j, k) * zC(i + 1, j, k) - yC(i + 1, j, k) * zcta(i + 1, j, k));
                sy2 = half * (zcta(i + 1, j, k) * xC(i + 1, j, k) - zC(i + 1, j, k) * xcta(i + 1, j, k));
                sz2 = half * (xcta(i + 1, j, k) * yC(i + 1, j, k) - xC(i + 1, j, k) * ycta(i + 1, j, k));
                sx3 = half * (ycta(i + 2, j, k) * zC(i + 2, j, k) - yC(i + 2, j, k) * zcta(i + 2, j, k));
                sy3 = half * (zcta(i + 2, j, k) * xC(i + 2, j, k) - zC(i + 2, j, k) * xcta(i + 2, j, k));
                sz3 = half * (xcta(i + 2, j, k) * yC(i + 2, j, k) - xC(i + 2, j, k) * ycta(i + 2, j, k));                

                if ((*DiscretePrecisionKXI)(i) == 2)
                {
                    KX += (sx2 - sx1);
                    KY += (sy2 - sy1);
                    KZ += (sz2 - sz1);
                }
                else
                {
                    KX += (27.0 * (sx2 - sx1) - (sx3 - sx0)) / 24.0;
                    KY += (27.0 * (sy2 - sy1) - (sy3 - sy0)) / 24.0;
                    KZ += (27.0 * (sz2 - sz1) - (sz3 - sz0)) / 24.0;
                }                

                sx0 = half * (yA(i, j, k - 1) * zkxi(i, j, k - 1) - ykxi(i, j, k - 1) * zA(i, j, k - 1));
                sy0 = half * (zA(i, j, k - 1) * xkxi(i, j, k - 1) - zkxi(i, j, k - 1) * xA(i, j, k - 1));
                sz0 = half * (xA(i, j, k - 1) * ykxi(i, j, k - 1) - xkxi(i, j, k - 1) * yA(i, j, k - 1));
                sx1 = half * (yA(i, j, k) * zkxi(i, j, k) - ykxi(i, j, k) * zA(i, j, k));
                sy1 = half * (zA(i, j, k) * xkxi(i, j, k) - zkxi(i, j, k) * xA(i, j, k));
                sz1 = half * (xA(i, j, k) * ykxi(i, j, k) - xkxi(i, j, k) * yA(i, j, k));
                sx2 = half * (yA(i, j, k + 1) * zkxi(i, j, k + 1) - ykxi(i, j, k + 1) * zA(i, j, k + 1));
                sy2 = half * (zA(i, j, k + 1) * xkxi(i, j, k + 1) - zkxi(i, j, k + 1) * xA(i, j, k + 1));
                sz2 = half * (xA(i, j, k + 1) * ykxi(i, j, k + 1) - xkxi(i, j, k + 1) * yA(i, j, k + 1));
                sx3 = half * (yA(i, j, k + 2) * zkxi(i, j, k + 2) - ykxi(i, j, k + 2) * zA(i, j, k + 2));
                sy3 = half * (zA(i, j, k + 2) * xkxi(i, j, k + 2) - zkxi(i, j, k + 2) * xA(i, j, k + 2));
                sz3 = half * (xA(i, j, k + 2) * ykxi(i, j, k + 2) - xkxi(i, j, k + 2) * yA(i, j, k + 2));                

                if ((*DiscretePrecisionCTA)(k) == 2)
                {
                    KX += (sx2 - sx1);
                    KY += (sy2 - sy1);
                    KZ += (sz2 - sz1);
                }
                else
                {
                    KX += (27.0 * (sx2 - sx1) - (sx3 - sx0)) / 24.0;
                    KY += (27.0 * (sy2 - sy1) - (sy3 - sy0)) / 24.0;
                    KZ += (27.0 * (sz2 - sz1) - (sz3 - sz0)) / 24.0;
                }

                RDouble ds = sqrt(KX * KX + KY * KY + KZ * KZ);

                xfv (i, j, k, 2) = KX;
                yfv (i, j, k, 2) = KY;
                zfv (i, j, k, 2) = KZ;
                area(i, j, k, 2) = ds;
                xfn (i, j, k, 2) = KX / (ds + TINY);
                yfn (i, j, k, 2) = KY / (ds + TINY);
                zfn (i, j, k, 2) = KZ / (ds + TINY);
            }
        }
    }

    for (int j = 1; j <= nj - 1; ++ j)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int k = 1; k <= nk; ++ k)
            {
                RDouble sx0, sy0, sz0;
                RDouble sx1, sy1, sz1;
                RDouble sx2, sy2, sz2;
                RDouble sx3, sy3, sz3;

                RDouble KX = 0.0;
                RDouble KY = 0.0;
                RDouble KZ = 0.0;

                sx0 = half * (ykxi(i, j - 1, k) * zA(i, j - 1, k) - yA(i, j - 1, k) * zkxi(i, j - 1, k));
                sy0 = half * (zkxi(i, j - 1, k) * xA(i, j - 1, k) - zA(i, j - 1, k) * xkxi(i, j - 1, k));
                sz0 = half * (xkxi(i, j - 1, k) * yA(i, j - 1, k) - xA(i, j - 1, k) * ykxi(i, j - 1, k));
                sx1 = half * (ykxi(i, j,  k) * zA(i, j,  k) - yA(i, j,  k) * zkxi(i, j,  k));
                sy1 = half * (zkxi(i, j,  k) * xA(i, j,  k) - zA(i, j,  k) * xkxi(i, j,  k));
                sz1 = half * (xkxi(i, j,  k) * yA(i, j,  k) - xA(i, j,  k) * ykxi(i, j,  k));
                sx2 = half * (ykxi(i, j + 1, k) * zA(i, j + 1, k) - yA(i, j + 1, k) * zkxi(i, j + 1, k));
                sy2 = half * (zkxi(i, j + 1, k) * xA(i, j + 1, k) - zA(i, j + 1, k) * xkxi(i, j + 1, k));
                sz2 = half * (xkxi(i, j + 1, k) * yA(i, j + 1, k) - xA(i, j + 1, k) * ykxi(i, j + 1, k));
                sx3 = half * (ykxi(i, j + 2, k) * zA(i, j + 2, k) - yA(i, j + 2, k) * zkxi(i, j + 2, k));
                sy3 = half * (zkxi(i, j + 2, k) * xA(i, j + 2, k) - zA(i, j + 2, k) * xkxi(i, j + 2, k));
                sz3 = half * (xkxi(i, j + 2, k) * yA(i, j + 2, k) - xA(i, j + 2, k) * ykxi(i, j + 2, k));                

                if ((*DiscretePrecisionETA)(j) == 2)
                {
                    KX += (sx2 - sx1);
                    KY += (sy2 - sy1);
                    KZ += (sz2 - sz1);
                }
                else
                {
                    KX += (27.0 * (sx2 - sx1) - (sx3 - sx0)) / 24.0;
                    KY += (27.0 * (sy2 - sy1) - (sy3 - sy0)) / 24.0;
                    KZ += (27.0 * (sz2 - sz1) - (sz3 - sz0)) / 24.0;
                }


                sx0 = half * (yB(i - 1, j, k) * zeta(i - 1, j, k) - yeta(i - 1, j, k) * zB(i - 1, j, k));
                sy0 = half * (zB(i - 1, j, k) * xeta(i - 1, j, k) - zeta(i - 1, j, k) * xB(i - 1, j, k));
                sz0 = half * (xB(i - 1, j, k) * yeta(i - 1, j, k) - xeta(i - 1, j, k) * yB(i - 1, j, k));
                sx1 = half * (yB(i  , j, k) * zeta(i  , j, k) - yeta(i  , j, k) * zB(i  , j, k));
                sy1 = half * (zB(i  , j, k) * xeta(i  , j, k) - zeta(i  , j, k) * xB(i  , j, k));
                sz1 = half * (xB(i  , j, k) * yeta(i  , j, k) - xeta(i  , j, k) * yB(i  , j, k));
                sx2 = half * (yB(i + 1, j, k) * zeta(i + 1, j, k) - yeta(i + 1, j, k) * zB(i + 1, j, k));
                sy2 = half * (zB(i + 1, j, k) * xeta(i + 1, j, k) - zeta(i + 1, j, k) * xB(i + 1, j, k));
                sz2 = half * (xB(i + 1, j, k) * yeta(i + 1, j, k) - xeta(i + 1, j, k) * yB(i + 1, j, k));
                sx3 = half * (yB(i + 2, j, k) * zeta(i + 2, j, k) - yeta(i + 2, j, k) * zB(i + 2, j, k));
                sy3 = half * (zB(i + 2, j, k) * xeta(i + 2, j, k) - zeta(i + 2, j, k) * xB(i + 2, j, k));
                sz3 = half * (xB(i + 2, j, k) * yeta(i + 2, j, k) - xeta(i + 2, j, k) * yB(i + 2, j, k));                

                if ((*DiscretePrecisionKXI)(i) == 2)
                {
                    KX += (sx2 - sx1);
                    KY += (sy2 - sy1);
                    KZ += (sz2 - sz1);
                }
                else
                {
                    KX += (27.0 * (sx2 - sx1) - (sx3 - sx0)) / 24.0;
                    KY += (27.0 * (sy2 - sy1) - (sy3 - sy0)) / 24.0;
                    KZ += (27.0 * (sz2 - sz1) - (sz3 - sz0)) / 24.0;
                }
                
                RDouble ds = sqrt(KX * KX + KY * KY + KZ * KZ);

                xfv (i, j, k, 3) = KX;
                yfv (i, j, k, 3) = KY;
                zfv (i, j, k, 3) = KZ;
                area(i, j, k, 3) = ds;
                xfn (i, j, k, 3) = KX / (ds + TINY);
                yfn (i, j, k, 3) = KY / (ds + TINY);
                zfn (i, j, k, 3) = KZ / (ds + TINY);
            }
        }
    }
    
    delete &xA   ;
    delete &yA   ;
    delete &zA   ;
    delete &xkxi;
    delete &ykxi;
    delete &zkxi;

    delete &xB   ;
    delete &yB   ;
    delete &zB   ;
    delete &xeta;
    delete &yeta;
    delete &zeta;

    delete &xC   ;
    delete &yC   ;
    delete &zC   ;
    delete &xcta;
    delete &ycta;
    delete &zcta;

    //! BC metrics
    //! interior BC would be replaced in the next function CommPartofMetrics_FD.
    GhostMetrics3D();
}

void StructGrid::ComputeJacobianforStructHighOrder()
{
    if (GetDim() == 2)
    {
        ComputeJacobianforStructHighOrder2D();
    }
    else
    {
        ComputeJacobianforStructHighOrder3D();
    }
}

void StructGrid::ComputeJacobianforStructHighOrder2D()
{
    RDouble4D &xfv  = *(this->GetFaceVectorX());
    RDouble4D &yfv  = *(this->GetFaceVectorY());

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();

    RDouble3D &vol = *this->GetCellVolume();    

    int ni = this->GetNI();
    int nj = this->GetNJ();
    const int k = 1;

    RDouble3D &xfaceI = *new RDouble3D(Range(0, ni + 1), Range(1, nj - 1), Range(1, 1), fortranArray);
    RDouble3D &yfaceI = *new RDouble3D(Range(0, ni + 1), Range(1, nj - 1), Range(1, 1), fortranArray);

    RDouble3D &xfaceJ = *new RDouble3D(Range(1, ni - 1), Range(0, nj + 1), Range(1, 1), fortranArray);
    RDouble3D &yfaceJ = *new RDouble3D(Range(1, ni - 1), Range(0, nj + 1), Range(1, 1), fortranArray);

    for (int j = 1; j <= nj - 1; ++ j)
    {
        for (int i = 0; i <= ni + 1; ++ i)
        {
            RDouble x1 = xx(i, j,  k);
            RDouble x2 = xx(i, j + 1, k);

            RDouble y1 = yy(i, j,  k);
            RDouble y2 = yy(i, j + 1, k);
            
            xfaceI(i, j, k) = half * (x1 + x2);
            yfaceI(i, j, k) = half * (y1 + y2);
        }
    }
    for (int i = 1; i <= ni - 1; ++ i)
    {
        for (int j = 0; j <= nj + 1; ++ j)
        {
            RDouble x1 = xx(i,  j, k);
            RDouble x2 = xx(i + 1, j, k);

            RDouble y1 = yy(i,  j, k);
            RDouble y2 = yy(i + 1, j, k);
            
            xfaceJ(i, j, k) = half * (x1 + x2);
            yfaceJ(i, j, k) = half * (y1 + y2);
        }
    }

    for (int i = 1; i <= ni - 1; ++ i)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            RDouble vol_kxi;
            RDouble vol_eta;

            RDouble kxi0 = xfaceI(i - 1, j, k) * xfv(i - 1, j, k, 1) + yfaceI(i - 1, j, k) * yfv(i - 1, j, k, 1);
            RDouble kxi1 = xfaceI(i,     j, k) * xfv(i,     j, k, 1) + yfaceI(i,     j, k) * yfv(i,     j, k, 1);
            RDouble kxi2 = xfaceI(i + 1, j, k) * xfv(i + 1, j, k, 1) + yfaceI(i + 1, j, k) * yfv(i + 1, j, k, 1);
            RDouble kxi3 = xfaceI(i + 2, j, k) * xfv(i + 2, j, k, 1) + yfaceI(i + 2, j, k) * yfv(i + 2, j, k, 1);

            RDouble eta0 = xfaceJ(i, j - 1, k) * xfv(i, j - 1, k, 2) + yfaceJ(i, j - 1, k) * yfv(i, j - 1, k, 2);
            RDouble eta1 = xfaceJ(i, j,     k) * xfv(i, j,     k, 2) + yfaceJ(i, j,     k) * yfv(i, j,     k, 2);
            RDouble eta2 = xfaceJ(i, j + 1, k) * xfv(i, j + 1, k, 2) + yfaceJ(i, j + 1, k) * yfv(i, j + 1, k, 2);
            RDouble eta3 = xfaceJ(i, j + 2, k) * xfv(i, j + 2, k, 2) + yfaceJ(i, j + 2, k) * yfv(i, j + 2, k, 2);

            if ((*DiscretePrecisionKXI)(i) == 2)
            {
                vol_kxi= (kxi2 - kxi1);
            }
            else
            {
                vol_kxi= (27.0 * (kxi2 - kxi1) - (kxi3 - kxi0)) / 24.0;
            }
            
            if ((*DiscretePrecisionETA)(j) == 2)
            {
                vol_eta= (eta2 - eta1);
            }
            else
            {
                vol_eta= (27.0 * (eta2 - eta1) - (eta3 - eta0)) / 24.0;
            }

            RDouble vol_FD = half * (vol_kxi + vol_eta);
            RDouble vol_FV = vol(i, j, k);

            if (vol_FD <= zero || vol_FD > 1.2 * vol_FV || vol_FD < 0.8 * vol_FV)
            {
                vol_FD = vol_FV;
            }

            vol(i, j, k) = vol_FD;
        }        
    }
    
    GhostCell3D(vol, ni, nj, 1);
    //! ghost vol
    //! interior BC would be replaced in the next function CommJacobian_FD.

    delete &xfaceI;
    delete &yfaceI;
    delete &xfaceJ;
    delete &yfaceJ;
}

void StructGrid::ComputeJacobianforStructHighOrder3D()
{
    RDouble4D &xfv  = *(this->GetFaceVectorX());
    RDouble4D &yfv  = *(this->GetFaceVectorY());
    RDouble4D &zfv  = *(this->GetFaceVectorZ());

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    RDouble3D &vol = *this->GetCellVolume();

    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    RDouble3D &xfaceI = *new RDouble3D(Range(0, ni + 1), Range(1, nj - 1), Range(1, nk - 1), fortranArray);
    RDouble3D &yfaceI = *new RDouble3D(Range(0, ni + 1), Range(1, nj - 1), Range(1, nk - 1), fortranArray);
    RDouble3D &zfaceI = *new RDouble3D(Range(0, ni + 1), Range(1, nj - 1), Range(1, nk - 1), fortranArray);

    RDouble3D &xfaceJ = *new RDouble3D(Range(1, ni - 1), Range(0, nj + 1), Range(1, nk - 1), fortranArray);
    RDouble3D &yfaceJ = *new RDouble3D(Range(1, ni - 1), Range(0, nj + 1), Range(1, nk - 1), fortranArray);
    RDouble3D &zfaceJ = *new RDouble3D(Range(1, ni - 1), Range(0, nj + 1), Range(1, nk - 1), fortranArray);

    RDouble3D &xfaceK = *new RDouble3D(Range(1, ni - 1), Range(1, nj - 1), Range(0, nk + 1), fortranArray);
    RDouble3D &yfaceK = *new RDouble3D(Range(1, ni - 1), Range(1, nj - 1), Range(0, nk + 1), fortranArray);
    RDouble3D &zfaceK = *new RDouble3D(Range(1, ni - 1), Range(1, nj - 1), Range(0, nk + 1), fortranArray);

    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = 0; i <= ni + 1; ++ i)
            {
                RDouble x1 = xx(i, j,     k);
                RDouble x2 = xx(i, j + 1, k);
                RDouble x3 = xx(i, j,     k + 1);
                RDouble x4 = xx(i, j + 1, k + 1);

                RDouble y1 = yy(i, j,     k);
                RDouble y2 = yy(i, j + 1, k);
                RDouble y3 = yy(i, j,     k + 1);
                RDouble y4 = yy(i, j + 1, k + 1);

                RDouble z1 = zz(i, j,     k);
                RDouble z2 = zz(i, j + 1, k);
                RDouble z3 = zz(i, j,     k + 1);
                RDouble z4 = zz(i, j + 1, k + 1);
                
                xfaceI(i, j, k) = 0.25 * (x1 + x2 + x3 + x4);
                yfaceI(i, j, k) = 0.25 * (y1 + y2 + y3 + y4);
                zfaceI(i, j, k) = 0.25 * (z1 + z2 + z3 + z4);
            }
        }
    }

    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int j = 0; j <= nj + 1; ++ j)
            {
                RDouble x1 = xx(i,     j, k);
                RDouble x2 = xx(i + 1, j, k);
                RDouble x3 = xx(i,     j, k + 1);
                RDouble x4 = xx(i + 1, j, k + 1);

                RDouble y1 = yy(i,     j, k);
                RDouble y2 = yy(i + 1, j, k);
                RDouble y3 = yy(i,     j, k + 1);
                RDouble y4 = yy(i + 1, j, k + 1);

                RDouble z1 = zz(i,     j, k);
                RDouble z2 = zz(i + 1, j, k);
                RDouble z3 = zz(i,     j, k + 1);
                RDouble z4 = zz(i + 1, j, k + 1);
                
                xfaceJ(i, j, k) = 0.25 * (x1 + x2 + x3 + x4);
                yfaceJ(i, j, k) = 0.25 * (y1 + y2 + y3 + y4);
                zfaceJ(i, j, k) = 0.25 * (z1 + z2 + z3 + z4);
            }
        }
    }

    for (int j = 1; j <= nj - 1; ++ j)
    {
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int k = 0; k <= nk + 1; ++ k)
            {
                RDouble x1 = xx(i,     j,     k);
                RDouble x2 = xx(i + 1, j,     k);
                RDouble x3 = xx(i,     j + 1, k);
                RDouble x4 = xx(i + 1, j + 1, k);

                RDouble y1 = yy(i,     j,     k);
                RDouble y2 = yy(i + 1, j,     k);
                RDouble y3 = yy(i,     j + 1, k);
                RDouble y4 = yy(i + 1, j + 1, k);

                RDouble z1 = zz(i,     j,     k);
                RDouble z2 = zz(i + 1, j,     k);
                RDouble z3 = zz(i,     j + 1, k);
                RDouble z4 = zz(i + 1, j + 1, k);
                
                xfaceK(i, j, k) = 0.25 * (x1 + x2 + x3 + x4);
                yfaceK(i, j, k) = 0.25 * (y1 + y2 + y3 + y4);
                zfaceK(i, j, k) = 0.25 * (z1 + z2 + z3 + z4);
            }
        }
    }

    for (int k = 1; k <= nk - 1; ++ k)
    {
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = 1; i <= ni - 1; ++ i)
            {
                RDouble vol_kxi, vol_eta, vol_cta;

                RDouble kxi0 = xfaceI(i - 1, j, k) * xfv(i - 1, j, k, 1) + yfaceI(i - 1, j, k) * yfv(i - 1, j, k, 1) + zfaceI(i - 1, j, k) * zfv(i - 1, j, k, 1);
                RDouble kxi1 = xfaceI(i,     j, k) * xfv(i,     j, k, 1) + yfaceI(i,     j, k) * yfv(i,     j, k, 1) + zfaceI(i,     j, k) * zfv(i,     j, k, 1);
                RDouble kxi2 = xfaceI(i + 1, j, k) * xfv(i + 1, j, k, 1) + yfaceI(i + 1, j, k) * yfv(i + 1, j, k, 1) + zfaceI(i + 1, j, k) * zfv(i + 1, j, k, 1);
                RDouble kxi3 = xfaceI(i + 2, j, k) * xfv(i + 2, j, k, 1) + yfaceI(i + 2, j, k) * yfv(i + 2, j, k, 1) + zfaceI(i + 2, j, k) * zfv(i + 2, j, k, 1);
                
                if ((*DiscretePrecisionKXI)(i) == 2)
                {
                    vol_kxi= (kxi2 - kxi1);
                }
                else
                {
                    vol_kxi= (27.0 * (kxi2 - kxi1) - (kxi3 - kxi0)) / 24.0;
                }                

                RDouble eta0 = xfaceJ(i, j - 1, k) * xfv(i, j - 1, k, 2) + yfaceJ(i, j - 1, k) * yfv(i, j - 1, k, 2) + zfaceJ(i, j - 1, k) * zfv(i, j - 1, k, 2);
                RDouble eta1 = xfaceJ(i, j,     k) * xfv(i, j,     k, 2) + yfaceJ(i, j,     k) * yfv(i, j,     k, 2) + zfaceJ(i, j,     k) * zfv(i, j,     k, 2);
                RDouble eta2 = xfaceJ(i, j + 1, k) * xfv(i, j + 1, k, 2) + yfaceJ(i, j + 1, k) * yfv(i, j + 1, k, 2) + zfaceJ(i, j + 1, k) * zfv(i, j + 1, k, 2);
                RDouble eta3 = xfaceJ(i, j + 2, k) * xfv(i, j + 2, k, 2) + yfaceJ(i, j + 2, k) * yfv(i, j + 2, k, 2) + zfaceJ(i, j + 2, k) * zfv(i, j + 2, k, 2);
                
                if ((*DiscretePrecisionETA)(j) == 2)
                {
                    vol_eta= (eta2 - eta1);
                }
                else
                {
                    vol_eta= (27.0 * (eta2 - eta1) - (eta3 - eta0)) / 24.0;
                }
                
                RDouble cta0 = xfaceK(i, j, k - 1) * xfv(i, j, k - 1, 3) + yfaceK(i, j, k - 1) * yfv(i, j, k - 1, 3) + zfaceK(i, j, k - 1) * zfv(i, j, k - 1, 3);
                RDouble cta1 = xfaceK(i, j, k   ) * xfv(i, j, k,     3) + yfaceK(i, j, k   ) * yfv(i, j, k,     3) + zfaceK(i, j, k   ) * zfv(i, j, k,     3);
                RDouble cta2 = xfaceK(i, j, k + 1) * xfv(i, j, k + 1, 3) + yfaceK(i, j, k + 1) * yfv(i, j, k + 1, 3) + zfaceK(i, j, k + 1) * zfv(i, j, k + 1, 3);
                RDouble cta3 = xfaceK(i, j, k + 2) * xfv(i, j, k + 2, 3) + yfaceK(i, j, k + 2) * yfv(i, j, k + 2, 3) + zfaceK(i, j, k + 2) * zfv(i, j, k + 2, 3);
                
                if ((*DiscretePrecisionCTA)(k) == 2)
                {
                    vol_cta= (cta2 - cta1);
                }
                else
                {
                    vol_cta= (27.0 * (cta2 - cta1) - (cta3 - cta0)) / 24.0;
                } 
                
                RDouble vol_FD = (vol_kxi + vol_eta + vol_cta) / 3.0;
                RDouble vol_FV = vol(i, j, k);
                
                if (vol_FD <= zero || vol_FD > 1.2 * vol_FV || vol_FD < 0.8 * vol_FV)
                {
                    vol_FD = vol_FV;
                }
                vol(i ,j ,k) = vol_FD;
            }
        }
    }

    GhostCell3D(vol, ni, nj, nk);
    //! ghost vol
    //! interior BC would be replaced in the next function CommJacobian_FD.

    delete &xfaceI;
    delete &yfaceI;
    delete &zfaceI;

    delete &xfaceJ;
    delete &yfaceJ;
    delete &zfaceJ;

    delete &xfaceK;
    delete &yfaceK;
    delete &zfaceK;
}

void StructGrid::GetNodeIJKfromSourceIndexghostofghost(StructGrid *fgrid, int iFace, int layer, int IJKofNodes[12][3])
{
    int is, js, ks;
    fgrid->GetSourceIndexIJK(iFace, layer, is, js, ks);

    int id, jd, kd;
    int i_lr, j_lr, k_lr;
    fgrid->GetSourceIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);

    int a[3];
    int b[3];
    int c[3];
    int d[3];
    int ab[3];
    int ad[3];
    int bc[3];
    int ba[3];
    int cd[3];
    int cb[3];
    int da[3];
    int dc[3];
    
    if (i_lr != 0)
    {
        //! I_face
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;

        b[0] = 0;
        b[1] = 1;
        b[2] = 0;

        c[0] = 0;
        c[1] = 1;
        c[2] = 1;

        d[0] = 0;
        d[1] = 0;
        d[2] = 1;
    }
    else if (j_lr != 0)
    {
        //! J_face
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;

        b[0] = 0;
        b[1] = 0;
        b[2] = 1;

        c[0] = 1;
        c[1] = 0;
        c[2] = 1;

        d[0] = 1;
        d[1] = 0;
        d[2] = 0;
    }
    else
    {
        //! K_face
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;

        b[0] = 1;
        b[1] = 0;
        b[2] = 0;

        c[0] = 1;
        c[1] = 1;
        c[2] = 0;

        d[0] = 0;
        d[1] = 1;
        d[2] = 0;
    }

    for (int m = 0; m <= 2; ++ m)
    {
        ab[m] = 2 * a[m] - b[m];
        ad[m] = 2 * a[m] - d[m];
        bc[m] = 2 * b[m] - c[m];
        ba[m] = 2 * b[m] - a[m];
        cd[m] = 2 * c[m] - d[m];
        cb[m] = 2 * c[m] - b[m];
        da[m] = 2 * d[m] - a[m];
        dc[m] = 2 * d[m] - c[m];
    }

    int basepoint[3];
    basepoint[0] = is + id;
    basepoint[1] = js + jd;
    basepoint[2] = ks + kd;

    for (int m = 0; m <= 2; ++ m)
    {
        IJKofNodes[ 0][m] = basepoint[m] + a[m];
        IJKofNodes[ 1][m] = basepoint[m] + b[m];
        IJKofNodes[ 2][m] = basepoint[m] + c[m];
        IJKofNodes[ 3][m] = basepoint[m] + d[m];

        IJKofNodes[ 4][m] = basepoint[m] + ab[m];
        IJKofNodes[ 5][m] = basepoint[m] + ad[m];
        IJKofNodes[ 6][m] = basepoint[m] + bc[m];
        IJKofNodes[ 7][m] = basepoint[m] + ba[m];
        IJKofNodes[ 8][m] = basepoint[m] + cd[m];
        IJKofNodes[ 9][m] = basepoint[m] + cb[m];
        IJKofNodes[10][m] = basepoint[m] + da[m];
        IJKofNodes[11][m] = basepoint[m] + dc[m];
    }
}

void StructGrid::GetNodeIJKfromTargetIndexghostofghost(StructGrid *fgrid, int iFace, int layer, int IJKofNodes[12][3])
{
    int it, jt, kt;
    fgrid->GetTargetIndexIJK(iFace, layer, it, jt, kt);

    int id, jd, kd;
    int i_lr, j_lr, k_lr;
    fgrid->GetTargetIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);

    int a[3];
    int b[3];
    int c[3];
    int d[3];
    int ab[3];
    int ad[3];
    int bc[3];
    int ba[3];
    int cd[3];
    int cb[3];
    int da[3];
    int dc[3];
    
    if (i_lr != 0)
    {
        //! I_face
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;

        b[0] = 0;
        b[1] = 1;
        b[2] = 0;

        c[0] = 0;
        c[1] = 1;
        c[2] = 1;

        d[0] = 0;
        d[1] = 0;
        d[2] = 1;
    }
    else if (j_lr != 0)
    {
        //! J_face
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;

        b[0] = 0;
        b[1] = 0;
        b[2] = 1;

        c[0] = 1;
        c[1] = 0;
        c[2] = 1;

        d[0] = 1;
        d[1] = 0;
        d[2] = 0;
    }
    else
    {
        //! K_face
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;

        b[0] = 1;
        b[1] = 0;
        b[2] = 0;

        c[0] = 1;
        c[1] = 1;
        c[2] = 0;

        d[0] = 0;
        d[1] = 1;
        d[2] = 0;
    }

    for (int m = 0; m <= 2; ++ m)
    {
        ab[m] = 2 * a[m] - b[m];
        ad[m] = 2 * a[m] - d[m];
        bc[m] = 2 * b[m] - c[m];
        ba[m] = 2 * b[m] - a[m];
        cd[m] = 2 * c[m] - d[m];
        cb[m] = 2 * c[m] - b[m];
        da[m] = 2 * d[m] - a[m];
        dc[m] = 2 * d[m] - c[m];
    }

    int basepoint[3];
    basepoint[0] = it + id;
    basepoint[1] = jt + jd;
    basepoint[2] = kt + kd;

    for (int m = 0; m <= 2; ++ m)
    {
        IJKofNodes[ 0][m] = basepoint[m] + a[m];
        IJKofNodes[ 1][m] = basepoint[m] + b[m];
        IJKofNodes[ 2][m] = basepoint[m] + c[m];
        IJKofNodes[ 3][m] = basepoint[m] + d[m];

        IJKofNodes[ 4][m] = basepoint[m] + ab[m];
        IJKofNodes[ 5][m] = basepoint[m] + ad[m];
        IJKofNodes[ 6][m] = basepoint[m] + bc[m];
        IJKofNodes[ 7][m] = basepoint[m] + ba[m];
        IJKofNodes[ 8][m] = basepoint[m] + cd[m];
        IJKofNodes[ 9][m] = basepoint[m] + cb[m];
        IJKofNodes[10][m] = basepoint[m] + da[m];
        IJKofNodes[11][m] = basepoint[m] + dc[m];
    }
}

int StructGrid::GetRelationshipBetweenTwoPairsofNodes(RDouble NodesofABCD[4][3], RDouble Nodesofabcd[4][3])
{
    RDouble length2ofAa = 0.0;
    RDouble length2ofAb = 0.0;
    RDouble length2ofAc = 0.0;
    RDouble length2ofAd = 0.0;
    
    RDouble length2ofBa = 0.0;
    RDouble length2ofBb = 0.0;
    RDouble length2ofBc = 0.0;
    RDouble length2ofBd = 0.0;
    
    RDouble length2ofCa = 0.0;
    RDouble length2ofCb = 0.0;
    RDouble length2ofCc = 0.0;
    RDouble length2ofCd = 0.0;
    
    RDouble length2ofDa = 0.0;
    RDouble length2ofDb = 0.0;
    RDouble length2ofDc = 0.0;
    RDouble length2ofDd = 0.0;

    for (int m = 0; m <= 2; ++ m)
    {
        length2ofAa += (NodesofABCD[0][m] - Nodesofabcd[0][m]) * (NodesofABCD[0][m] - Nodesofabcd[0][m]);
        length2ofAb += (NodesofABCD[0][m] - Nodesofabcd[1][m]) * (NodesofABCD[0][m] - Nodesofabcd[1][m]);
        length2ofAc += (NodesofABCD[0][m] - Nodesofabcd[2][m]) * (NodesofABCD[0][m] - Nodesofabcd[2][m]);
        length2ofAd += (NodesofABCD[0][m] - Nodesofabcd[3][m]) * (NodesofABCD[0][m] - Nodesofabcd[3][m]);
        
        length2ofBa += (NodesofABCD[1][m] - Nodesofabcd[0][m]) * (NodesofABCD[1][m] - Nodesofabcd[0][m]);
        length2ofBb += (NodesofABCD[1][m] - Nodesofabcd[1][m]) * (NodesofABCD[1][m] - Nodesofabcd[1][m]);
        length2ofBc += (NodesofABCD[1][m] - Nodesofabcd[2][m]) * (NodesofABCD[1][m] - Nodesofabcd[2][m]);
        length2ofBd += (NodesofABCD[1][m] - Nodesofabcd[3][m]) * (NodesofABCD[1][m] - Nodesofabcd[3][m]);
        
        length2ofCa += (NodesofABCD[2][m] - Nodesofabcd[0][m]) * (NodesofABCD[2][m] - Nodesofabcd[0][m]);
        length2ofCb += (NodesofABCD[2][m] - Nodesofabcd[1][m]) * (NodesofABCD[2][m] - Nodesofabcd[1][m]);
        length2ofCc += (NodesofABCD[2][m] - Nodesofabcd[2][m]) * (NodesofABCD[2][m] - Nodesofabcd[2][m]);
        length2ofCd += (NodesofABCD[2][m] - Nodesofabcd[3][m]) * (NodesofABCD[2][m] - Nodesofabcd[3][m]);
        
        length2ofDa += (NodesofABCD[3][m] - Nodesofabcd[0][m]) * (NodesofABCD[3][m] - Nodesofabcd[0][m]);
        length2ofDb += (NodesofABCD[3][m] - Nodesofabcd[1][m]) * (NodesofABCD[3][m] - Nodesofabcd[1][m]);
        length2ofDc += (NodesofABCD[3][m] - Nodesofabcd[2][m]) * (NodesofABCD[3][m] - Nodesofabcd[2][m]);
        length2ofDd += (NodesofABCD[3][m] - Nodesofabcd[3][m]) * (NodesofABCD[3][m] - Nodesofabcd[3][m]);
    }

    RDouble errorofcase1 = length2ofAa + length2ofBb + length2ofCc + length2ofDd;
    RDouble errorofcase2 = length2ofAa + length2ofBd + length2ofCc + length2ofDb;
    RDouble errorofcase3 = length2ofAb + length2ofBa + length2ofCd + length2ofDc;
    RDouble errorofcase4 = length2ofAb + length2ofBc + length2ofCd + length2ofDa;
    RDouble errorofcase5 = length2ofAc + length2ofBb + length2ofCa + length2ofDd;
    RDouble errorofcase6 = length2ofAc + length2ofBd + length2ofCa + length2ofDb;
    RDouble errorofcase7 = length2ofAd + length2ofBa + length2ofCb + length2ofDc;
    RDouble errorofcase8 = length2ofAd + length2ofBc + length2ofCb + length2ofDa;

    RDouble errorMin = errorofcase1;
    int mapcase = 1;
    
    if (errorofcase2 < errorMin)
    {
        errorMin = errorofcase2;
        mapcase = 2;
    }

    if (errorofcase3 < errorMin)
    {
        errorMin = errorofcase3;
        mapcase = 3;
    }

    if (errorofcase4 < errorMin)
    {
        errorMin = errorofcase4;
        mapcase = 4;
    }

    if (errorofcase5 < errorMin)
    {
        errorMin = errorofcase5;
        mapcase = 5;
    }

    if (errorofcase6 < errorMin)
    {
        errorMin = errorofcase6;
        mapcase = 6;
    }

    if (errorofcase7 < errorMin)
    {
        errorMin = errorofcase7;
        mapcase = 7;
    }

    if (errorofcase8 < errorMin)
    {
        errorMin = errorofcase8;
        mapcase = 8;
    }

    return mapcase;
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder(DataContainer *&dataContainer, const int &neighborZoneIndex, int mission)
{
    //int mission_Comm(XYZghost                                ) = 1;
    //int mission_Comm(XYZghostofghos                          ) = 2;
    //int mission_Comm(Metrics                                 ) = 3;
    //int mission_Comm(Jacobian                                ) = 4;
    //int mission_Comm(q_FD                                    ) = 5;
    //int mission_Comm(rtem                                    ) = 6;
    //int mission_Comm(rtem, dqdxyznokxi,dqdxyznoeta,dqdxyznocta) = 7;
    //int mission_Comm(rtem, gradUVWTCellCenterXYZ             ) = 8;
    //int mission_Comm(vist, gradqturb                         ) = 9;
    //int mission_Comm(vist, gradqturb, blending               ) = 10;
    //int mission_Comm(vist                                    ) = 11 (for LES);
    
    if (mission == 1)
    {
        UploadGridInfoOnlyforStructHighOrder_mission1(dataContainer, neighborZoneIndex);
    }
    else if (mission == 2)
    {
        UploadGridInfoOnlyforStructHighOrder_mission2(dataContainer, neighborZoneIndex);
    }
    else if (mission == 3)
    {
        UploadGridInfoOnlyforStructHighOrder_mission3(dataContainer, neighborZoneIndex);
    }
    else if (mission == 4)
    {
        UploadGridInfoOnlyforStructHighOrder_mission4(dataContainer, neighborZoneIndex);
    }
    else if (mission == 5)
    {
        UploadGridInfoOnlyforStructHighOrder_mission5(dataContainer, neighborZoneIndex);
    }
    else if (mission == 6)
    {
        UploadGridInfoOnlyforStructHighOrder_mission6(dataContainer, neighborZoneIndex);
    }
    else if (mission == 7)
    {
        UploadGridInfoOnlyforStructHighOrder_mission7(dataContainer, neighborZoneIndex);
    }
    else if (mission == 8)
    {
        UploadGridInfoOnlyforStructHighOrder_mission8(dataContainer, neighborZoneIndex);
    }
    else if (mission == 9)
    {
        UploadGridInfoOnlyforStructHighOrder_mission9(dataContainer, neighborZoneIndex);
    }
    else if (mission == 10)
    {
        UploadGridInfoOnlyforStructHighOrder_mission10(dataContainer, neighborZoneIndex);
    }
    else if (mission == 11)
    {
        UploadGridInfoOnlyforStructHighOrder_mission11(dataContainer, neighborZoneIndex);
    }
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder(DataContainer *&dataContainer, const int &neighborZoneIndex, int mission)
{
    //int mission_Comm(XYZghost                                ) = 1;
    //int mission_Comm(XYZghostofghos                          ) = 2;
    //int mission_Comm(Metrics                                 ) = 3;
    //int mission_Comm(Jacobian                                ) = 4;
    //int mission_Comm(q_FD                                    ) = 5;
    //int mission_Comm(rtem                                    ) = 6;
    //int mission_Comm(rtem, dqdxyznokxi,dqdxyznoeta,dqdxyznocta) = 7;
    //int mission_Comm(rtem, gradUVWTCellCenterXYZ             ) = 8;
    //int mission_Comm(vist, gradqturb                         ) = 9;
    //int mission_Comm(vist, gradqturb, blending               ) = 10;
    //int mission_Comm(vist                                    ) = 11 (for LES);
    
    if (mission == 1)
    {
        DownloadGridInfoOnlyforStructHighOrder_mission1(dataContainer, neighborZoneIndex);
    }
    else if (mission == 2)
    {
        DownloadGridInfoOnlyforStructHighOrder_mission2(dataContainer, neighborZoneIndex);
    }
    else if (mission == 3)
    {
        DownloadGridInfoOnlyforStructHighOrder_mission3(dataContainer, neighborZoneIndex);
    }
    else if (mission == 4)
    {
        DownloadGridInfoOnlyforStructHighOrder_mission4(dataContainer, neighborZoneIndex);
    }
    else if (mission == 5)
    {
        DownloadGridInfoOnlyforStructHighOrder_mission5(dataContainer, neighborZoneIndex);
    }
    else if (mission == 6)
    {
        DownloadGridInfoOnlyforStructHighOrder_mission6(dataContainer, neighborZoneIndex);
    }
    else if (mission == 7)
    {
        DownloadGridInfoOnlyforStructHighOrder_mission7(dataContainer, neighborZoneIndex);
    }
    else if (mission == 8)
    {
        DownloadGridInfoOnlyforStructHighOrder_mission8(dataContainer, neighborZoneIndex);
    }
    else if (mission == 9)
    {
        DownloadGridInfoOnlyforStructHighOrder_mission9(dataContainer, neighborZoneIndex);
    }
    else if (mission == 10)
    {
        DownloadGridInfoOnlyforStructHighOrder_mission10(dataContainer, neighborZoneIndex);
    }
    else if (mission == 11)
    {
        DownloadGridInfoOnlyforStructHighOrder_mission11(dataContainer, neighborZoneIndex);
    }
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder_mission1(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &vol = *this->GetCellVolume();
    RDouble3D &xcc = *this->GetCellCenterX();
    RDouble3D &ycc = *this->GetCellCenterY();
    RDouble3D &zcc = *this->GetCellCenterZ();

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();
    
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];

        for (int layer = 1; layer <= 4; ++layer)
        {
            int is, js, ks;
            finestGrid->GetSourceIndexIJK(iFace,layer, is, js, ks);
            finestGrid->RemapMultigridIJK(level, is, js, ks);
            
            PHWrite(dataContainer, xcc(is, js, ks));
            PHWrite(dataContainer, ycc(is, js, ks));
            PHWrite(dataContainer, zcc(is, js, ks));
            PHWrite(dataContainer, vol(is, js, ks));
        }

        int IJKofNodes1stLayer[4][3];
        int IJKofNodes2ndLayer[4][3];

        GetNodeIJKfromSourceIndex(finestGrid, iFace, 1, IJKofNodes1stLayer);
        GetNodeIJKfromSourceIndex(finestGrid, iFace, 2, IJKofNodes2ndLayer);
        
        for (int iNode = 0; iNode < 4; ++ iNode)
        {
            PHWrite(dataContainer, xx(IJKofNodes1stLayer[iNode][0], IJKofNodes1stLayer[iNode][1], IJKofNodes1stLayer[iNode][2]));
            PHWrite(dataContainer, yy(IJKofNodes1stLayer[iNode][0], IJKofNodes1stLayer[iNode][1], IJKofNodes1stLayer[iNode][2]));
            PHWrite(dataContainer, zz(IJKofNodes1stLayer[iNode][0], IJKofNodes1stLayer[iNode][1], IJKofNodes1stLayer[iNode][2]));
            
            PHWrite(dataContainer, xx(IJKofNodes2ndLayer[iNode][0], IJKofNodes2ndLayer[iNode][1], IJKofNodes2ndLayer[iNode][2]));
            PHWrite(dataContainer, yy(IJKofNodes2ndLayer[iNode][0], IJKofNodes2ndLayer[iNode][1], IJKofNodes2ndLayer[iNode][2]));
            PHWrite(dataContainer, zz(IJKofNodes2ndLayer[iNode][0], IJKofNodes2ndLayer[iNode][1], IJKofNodes2ndLayer[iNode][2]));
        }
    }
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder_mission2(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();
    
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];

        int IJKofNodes1stLayer[12][3];
        int IJKofNodes2ndLayer[12][3];

        GetNodeIJKfromSourceIndexghostofghost(finestGrid, iFace, 1, IJKofNodes1stLayer);
        GetNodeIJKfromSourceIndexghostofghost(finestGrid, iFace, 2, IJKofNodes2ndLayer);

        for (int nodeID = 0; nodeID <= 3; ++nodeID)
        {
            PHWrite(dataContainer, xx(IJKofNodes1stLayer[nodeID][0], IJKofNodes1stLayer[nodeID][1], IJKofNodes1stLayer[nodeID][2]));
            PHWrite(dataContainer, yy(IJKofNodes1stLayer[nodeID][0], IJKofNodes1stLayer[nodeID][1], IJKofNodes1stLayer[nodeID][2]));
            PHWrite(dataContainer, zz(IJKofNodes1stLayer[nodeID][0], IJKofNodes1stLayer[nodeID][1], IJKofNodes1stLayer[nodeID][2]));
        }  

        for (int nodeID = 4; nodeID <= 11; ++nodeID)
        {
            PHWrite(dataContainer, xx(IJKofNodes2ndLayer[nodeID][0], IJKofNodes2ndLayer[nodeID][1], IJKofNodes2ndLayer[nodeID][2]));
            PHWrite(dataContainer, yy(IJKofNodes2ndLayer[nodeID][0], IJKofNodes2ndLayer[nodeID][1], IJKofNodes2ndLayer[nodeID][2]));
            PHWrite(dataContainer, zz(IJKofNodes2ndLayer[nodeID][0], IJKofNodes2ndLayer[nodeID][1], IJKofNodes2ndLayer[nodeID][2]));
        }
    }
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder_mission3(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();
    
    RDouble4D &xfv  = *(this->GetFaceVectorX());
    RDouble4D &yfv  = *(this->GetFaceVectorY());
    RDouble4D &zfv  = *(this->GetFaceVectorZ());
    
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];

        for (int layer = 1; layer <= 2; ++ layer)
        {
            int is, js, ks;
            finestGrid->GetSourceIndexIJK(iFace,layer, is, js, ks);
            finestGrid->RemapMultigridIJK(level, is, js, ks);

            RDouble KX, KY, KZ;

            int id, jd, kd;
            int i_lr, j_lr, k_lr;
            finestGrid->GetSourceIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);
            
            if (i_lr != 0)
            {
                //! I_face
                if (i_lr == -1)
                {
                    //! left boundary
                    KX = xfv(is+1, js, ks, 1);
                    KY = yfv(is+1, js, ks, 1);
                    KZ = zfv(is+1, js, ks, 1);
                }
                else
                {
                    //! right boundary
                    KX = - xfv(is, js, ks, 1);
                    KY = - yfv(is, js, ks, 1);
                    KZ = - zfv(is, js, ks, 1);
                }
            }
            else if (j_lr != 0)
            {
                //! J_face
                if (j_lr == -1)
                {
                    //! left boundary
                    KX = xfv(is, js+1, ks, 2);
                    KY = yfv(is, js+1, ks, 2);
                    KZ = zfv(is, js+1, ks, 2);
                }
                else
                {
                    //! right boundary
                    KX = - xfv(is, js, ks, 2);
                    KY = - yfv(is, js, ks, 2);
                    KZ = - zfv(is, js, ks, 2);
                }
            }
            else
            {
                //! K_face
                if (k_lr == -1)
                {
                    //! left boundary
                    KX = xfv(is, js, ks+1, 3);
                    KY = yfv(is, js, ks+1, 3);
                    KZ = zfv(is, js, ks+1, 3);
                }
                else
                {
                    //! right boundary
                    KX = - xfv(is, js, ks, 3);
                    KY = - yfv(is, js, ks, 3);
                    KZ = - zfv(is, js, ks, 3);
                }
            }
            PHWrite(dataContainer, KX);
            PHWrite(dataContainer, KY);
            PHWrite(dataContainer, KZ);
        }
    }
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder_mission4(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();
    
    RDouble3D &vol = *this->GetCellVolume();
    
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];

        for (int layer = 1; layer <= 4; ++layer)
        {
            int is, js, ks;
            finestGrid->GetSourceIndexIJK(iFace,layer, is, js, ks);
            finestGrid->RemapMultigridIJK(level, is, js, ks);            
            
            PHWrite(dataContainer, vol(is, js, ks));
        }
    }
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder_mission5(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble4D &q_FD = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("q_FD"));

    int nEquation = 5;

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        for (int layer = 3; layer <= 4; ++layer)
        {
            int is, js, ks;
            finestGrid->GetSourceIndexIJK(iFace, layer, is, js, ks);
            finestGrid->RemapMultigridIJK(0, is, js, ks);
            for (int m = 0; m < nEquation; ++ m)
            {
                PHWrite(dataContainer, q_FD(is, js, ks, m));
            }
        }
    }
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder_mission6(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (this->GetDataPtr("rtem"));

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        for (int layer = 1; layer <= 2; ++layer)
        {
            int is, js, ks;
            finestGrid->GetSourceIndexIJK(iFace, layer, is, js, ks);
            finestGrid->RemapMultigridIJK(0, is, js, ks);
            PHWrite(dataContainer, rtem(is, js, ks));
        }
    }
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder_mission7(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (this->GetDataPtr("rtem"));

    RDouble4D &dqdxnokxi = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdxnokxi"));
    RDouble4D &dqdynokxi = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdynokxi"));
    RDouble4D &dqdznokxi = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdznokxi"));

    RDouble4D &dqdxnoeta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdxnoeta"));
    RDouble4D &dqdynoeta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdynoeta"));
    RDouble4D &dqdznoeta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdznoeta"));

    RDouble4D &dqdxnocta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdxnocta"));
    RDouble4D &dqdynocta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdynocta"));
    RDouble4D &dqdznocta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdznocta"));

    RDouble dqdx_layer1[4];
    RDouble dqdy_layer1[4];
    RDouble dqdz_layer1[4];

    RDouble dqdx_layer2[4];
    RDouble dqdy_layer2[4];
    RDouble dqdz_layer2[4];

    RDouble dqdx_layer3[4];
    RDouble dqdy_layer3[4];
    RDouble dqdz_layer3[4];

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        
        int is1, js1, ks1; //! layer = 1
        int is2, js2, ks2; //! layer = 2
        int is3, js3, ks3; //! layer = 3

        finestGrid->GetSourceIndexIJK(iFace, 1, is1, js1, ks1);
        finestGrid->GetSourceIndexIJK(iFace, 2, is2, js2, ks2);
        finestGrid->GetSourceIndexIJK(iFace, 3, is3, js3, ks3);

        finestGrid->RemapMultigridIJK(0, is1, js1, ks1);
        finestGrid->RemapMultigridIJK(0, is2, js2, ks2);
        finestGrid->RemapMultigridIJK(0, is3, js3, ks3);
        
        int id, jd, kd;
        int i_lr, j_lr, k_lr;
        finestGrid->GetSourceIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);
        
        if (i_lr != 0)
        {
            //! I_face
            for (int m = 0; m < 4; ++ m)
            {
                dqdx_layer1[m] = dqdxnokxi(is1, js1, ks1, m);
                dqdy_layer1[m] = dqdynokxi(is1, js1, ks1, m);
                dqdz_layer1[m] = dqdznokxi(is1, js1, ks1, m);

                dqdx_layer2[m] = dqdxnokxi(is2, js2, ks2, m);
                dqdy_layer2[m] = dqdynokxi(is2, js2, ks2, m);
                dqdz_layer2[m] = dqdznokxi(is2, js2, ks2, m);

                dqdx_layer3[m] = dqdxnokxi(is3, js3, ks3, m);
                dqdy_layer3[m] = dqdynokxi(is3, js3, ks3, m);
                dqdz_layer3[m] = dqdznokxi(is3, js3, ks3, m);
            }
        }
        else if (j_lr != 0)
        {
            //! J_face
            for (int m = 0; m < 4; ++ m)
            {
                dqdx_layer1[m] = dqdxnoeta(is1, js1, ks1, m);
                dqdy_layer1[m] = dqdynoeta(is1, js1, ks1, m);
                dqdz_layer1[m] = dqdznoeta(is1, js1, ks1, m);

                dqdx_layer2[m] = dqdxnoeta(is2, js2, ks2, m);
                dqdy_layer2[m] = dqdynoeta(is2, js2, ks2, m);
                dqdz_layer2[m] = dqdznoeta(is2, js2, ks2, m);

                dqdx_layer3[m] = dqdxnoeta(is3, js3, ks3, m);
                dqdy_layer3[m] = dqdynoeta(is3, js3, ks3, m);
                dqdz_layer3[m] = dqdznoeta(is3, js3, ks3, m);
            }
        }
        else
        {
            //! K_face
            for (int m = 0; m < 4; ++ m)
            {
                dqdx_layer1[m] = dqdxnocta(is1, js1, ks1, m);
                dqdy_layer1[m] = dqdynocta(is1, js1, ks1, m);
                dqdz_layer1[m] = dqdznocta(is1, js1, ks1, m);

                dqdx_layer2[m] = dqdxnocta(is2, js2, ks2, m);
                dqdy_layer2[m] = dqdynocta(is2, js2, ks2, m);
                dqdz_layer2[m] = dqdznocta(is2, js2, ks2, m);

                dqdx_layer3[m] = dqdxnocta(is3, js3, ks3, m);
                dqdy_layer3[m] = dqdynocta(is3, js3, ks3, m);
                dqdz_layer3[m] = dqdznocta(is3, js3, ks3, m);
            }
        }
        
        PHWrite(dataContainer, rtem(is1, js1, ks1));
        PHWrite(dataContainer, rtem(is2, js2, ks2));

        for (int m = 0; m < 4; ++ m)
        {
            PHWrite(dataContainer, dqdx_layer1[m]);
            PHWrite(dataContainer, dqdy_layer1[m]);
            PHWrite(dataContainer, dqdz_layer1[m]);

            PHWrite(dataContainer, dqdx_layer2[m]);
            PHWrite(dataContainer, dqdy_layer2[m]);
            PHWrite(dataContainer, dqdz_layer2[m]);

            PHWrite(dataContainer, dqdx_layer3[m]);
            PHWrite(dataContainer, dqdy_layer3[m]);
            PHWrite(dataContainer, dqdz_layer3[m]);
        }
    }
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder_mission8(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (this->GetDataPtr("rtem"));
    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("gradUVWTCellCenterZ"));

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        
        int is1, js1, ks1; //! layer = 1
        int is2, js2, ks2; //! layer = 2
        int is3, js3, ks3; //! layer = 3

        finestGrid->GetSourceIndexIJK(iFace, 1, is1, js1, ks1);
        finestGrid->GetSourceIndexIJK(iFace, 2, is2, js2, ks2);
        finestGrid->GetSourceIndexIJK(iFace, 3, is3, js3, ks3);

        finestGrid->RemapMultigridIJK(0, is1, js1, ks1);
        finestGrid->RemapMultigridIJK(0, is2, js2, ks2);
        finestGrid->RemapMultigridIJK(0, is3, js3, ks3);

        PHWrite(dataContainer, rtem(is1, js1, ks1));
        PHWrite(dataContainer, rtem(is2, js2, ks2));

        PHWrite(dataContainer, gradUVWTCellCenterX(is1, js1, ks1,0));
        PHWrite(dataContainer, gradUVWTCellCenterX(is1, js1, ks1,1));
        PHWrite(dataContainer, gradUVWTCellCenterX(is1, js1, ks1,2));
        PHWrite(dataContainer, gradUVWTCellCenterX(is1, js1, ks1,3));

        PHWrite(dataContainer, gradUVWTCellCenterX(is2, js2, ks2,0));
        PHWrite(dataContainer, gradUVWTCellCenterX(is2, js2, ks2,1));
        PHWrite(dataContainer, gradUVWTCellCenterX(is2, js2, ks2,2));
        PHWrite(dataContainer, gradUVWTCellCenterX(is2, js2, ks2,3));

        PHWrite(dataContainer, gradUVWTCellCenterX(is3, js3, ks3,0));
        PHWrite(dataContainer, gradUVWTCellCenterX(is3, js3, ks3,1));
        PHWrite(dataContainer, gradUVWTCellCenterX(is3, js3, ks3,2));
        PHWrite(dataContainer, gradUVWTCellCenterX(is3, js3, ks3,3));

        PHWrite(dataContainer, gradUVWTCellCenterY(is1, js1, ks1,0));
        PHWrite(dataContainer, gradUVWTCellCenterY(is1, js1, ks1,1));
        PHWrite(dataContainer, gradUVWTCellCenterY(is1, js1, ks1,2));
        PHWrite(dataContainer, gradUVWTCellCenterY(is1, js1, ks1,3));

        PHWrite(dataContainer, gradUVWTCellCenterY(is2, js2, ks2,0));
        PHWrite(dataContainer, gradUVWTCellCenterY(is2, js2, ks2,1));
        PHWrite(dataContainer, gradUVWTCellCenterY(is2, js2, ks2,2));
        PHWrite(dataContainer, gradUVWTCellCenterY(is2, js2, ks2,3));

        PHWrite(dataContainer, gradUVWTCellCenterY(is3, js3, ks3,0));
        PHWrite(dataContainer, gradUVWTCellCenterY(is3, js3, ks3,1));
        PHWrite(dataContainer, gradUVWTCellCenterY(is3, js3, ks3,2));
        PHWrite(dataContainer, gradUVWTCellCenterY(is3, js3, ks3,3));

        PHWrite(dataContainer, gradUVWTCellCenterZ(is1, js1, ks1,0));
        PHWrite(dataContainer, gradUVWTCellCenterZ(is1, js1, ks1,1));
        PHWrite(dataContainer, gradUVWTCellCenterZ(is1, js1, ks1,2));
        PHWrite(dataContainer, gradUVWTCellCenterZ(is1, js1, ks1,3));

        PHWrite(dataContainer, gradUVWTCellCenterZ(is2, js2, ks2,0));
        PHWrite(dataContainer, gradUVWTCellCenterZ(is2, js2, ks2,1));
        PHWrite(dataContainer, gradUVWTCellCenterZ(is2, js2, ks2,2));
        PHWrite(dataContainer, gradUVWTCellCenterZ(is2, js2, ks2,3));

        PHWrite(dataContainer, gradUVWTCellCenterZ(is3, js3, ks3,0));
        PHWrite(dataContainer, gradUVWTCellCenterZ(is3, js3, ks3,1));
        PHWrite(dataContainer, gradUVWTCellCenterZ(is3, js3, ks3,2));
        PHWrite(dataContainer, gradUVWTCellCenterZ(is3, js3, ks3,3));
    }
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder_mission9(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &vist = *reinterpret_cast<RDouble3D *> (this->GetDataPtr("vist"));
    RDouble4D &gradTurbulenceCellCenterX = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterX"));
    RDouble4D &gradTurbulenceCellCenterY = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterY"));
    RDouble4D &gradTurbulenceCellCenterZ = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterZ"));

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        
        int is1, js1, ks1; //! layer = 1

        finestGrid->GetSourceIndexIJK(iFace, 1, is1, js1, ks1);

        finestGrid->RemapMultigridIJK(0, is1, js1, ks1);

        PHWrite(dataContainer, vist(is1, js1, ks1));

        PHWrite(dataContainer, gradTurbulenceCellCenterX(is1, js1, ks1,ISA));
        PHWrite(dataContainer, gradTurbulenceCellCenterY(is1, js1, ks1,ISA));
        PHWrite(dataContainer, gradTurbulenceCellCenterZ(is1, js1, ks1,ISA));
    }
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder_mission10(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &blend = *reinterpret_cast<RDouble3D *> (this->GetDataPtr("blend"));
    RDouble3D &vist  = *reinterpret_cast<RDouble3D *> (this->GetDataPtr("vist"));
    RDouble4D &gradTurbulenceCellCenterX = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterX"));
    RDouble4D &gradTurbulenceCellCenterY = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterY"));
    RDouble4D &gradTurbulenceCellCenterZ = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterZ"));

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        
        int is1, js1, ks1; //! layer = 1

        finestGrid->GetSourceIndexIJK(iFace, 1, is1, js1, ks1);

        finestGrid->RemapMultigridIJK(0, is1, js1, ks1);

        PHWrite(dataContainer, blend(is1, js1, ks1));

        PHWrite(dataContainer, vist(is1, js1, ks1));

        PHWrite(dataContainer, gradTurbulenceCellCenterX(is1, js1, ks1,IKE));
        PHWrite(dataContainer, gradTurbulenceCellCenterX(is1, js1, ks1,IKW));

        PHWrite(dataContainer, gradTurbulenceCellCenterY(is1, js1, ks1,IKE));
        PHWrite(dataContainer, gradTurbulenceCellCenterY(is1, js1, ks1,IKW));

        PHWrite(dataContainer, gradTurbulenceCellCenterZ(is1, js1, ks1,IKE));
        PHWrite(dataContainer, gradTurbulenceCellCenterZ(is1, js1, ks1,IKW));
    }
}

void StructGrid::UploadGridInfoOnlyforStructHighOrder_mission11(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &vist  = *reinterpret_cast<RDouble3D *> (this->GetDataPtr("vist"));
    RDouble3D &subgridScaleEnergy  = *reinterpret_cast<RDouble3D *> (this->GetDataPtr("subgridScaleEnergy"));
    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        
        int is1, js1, ks1; //! layer = 1

        finestGrid->GetSourceIndexIJK(iFace, 1, is1, js1, ks1);

        finestGrid->RemapMultigridIJK(0, is1, js1, ks1);

        PHWrite(dataContainer, vist(is1, js1, ks1));

        PHWrite(dataContainer, subgridScaleEnergy(is1, js1, ks1));

        PHWrite(dataContainer, turbulentPrandtlNumber(is1, js1, ks1));
    }
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder_mission1(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();
    
    RDouble3D &vol = *this->GetCellVolume();
    RDouble3D &xcc = *this->GetCellCenterX();
    RDouble3D &ycc = *this->GetCellCenterY();
    RDouble3D &zcc = *this->GetCellCenterZ();
    
    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();
    
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];

        for (int layer = 1; layer <= 4; ++ layer)
        {
            int it, jt, kt;
            finestGrid->GetTargetIndexIJK(iFace,layer, it, jt, kt);
            finestGrid->RemapMultigridIJK(level, it, jt, kt);
            
            PHRead(dataContainer, xcc(it, jt, kt));
            PHRead(dataContainer, ycc(it, jt, kt));
            PHRead(dataContainer, zcc(it, jt, kt));
            PHRead(dataContainer, vol(it, jt, kt));
        }

        RDouble Data1stLayer[4][3];
        RDouble Data2ndLayer[4][3];
        for (int nodeID = 0; nodeID < 4; ++nodeID)
        {
            PHRead(dataContainer, Data1stLayer[nodeID][0]);
            PHRead(dataContainer, Data1stLayer[nodeID][1]);
            PHRead(dataContainer, Data1stLayer[nodeID][2]);

            PHRead(dataContainer, Data2ndLayer[nodeID][0]);
            PHRead(dataContainer, Data2ndLayer[nodeID][1]);
            PHRead(dataContainer, Data2ndLayer[nodeID][2]);
        }  

        int IJKofNodes1stLayer[4][3];
        int IJKofNodes2ndLayer[4][3];
        GetNodeIJKfromTargetIndex(finestGrid, iFace, 1, IJKofNodes1stLayer);
        GetNodeIJKfromTargetIndex(finestGrid, iFace, 2, IJKofNodes2ndLayer);

        RDouble XYZofNodes1stLayer[4][3];
        for (int nodeID = 0; nodeID < 4; ++nodeID)
        {
            XYZofNodes1stLayer[nodeID][0] = xx(IJKofNodes1stLayer[nodeID][0], IJKofNodes1stLayer[nodeID][1], IJKofNodes1stLayer[nodeID][2]);
            XYZofNodes1stLayer[nodeID][1] = yy(IJKofNodes1stLayer[nodeID][0], IJKofNodes1stLayer[nodeID][1], IJKofNodes1stLayer[nodeID][2]);
            XYZofNodes1stLayer[nodeID][2] = zz(IJKofNodes1stLayer[nodeID][0], IJKofNodes1stLayer[nodeID][1], IJKofNodes1stLayer[nodeID][2]);
        }
        
        int RelationBetweenTwoPairsofNodes[4];
        GetRelationBetweenTwoPairsofNodes(RelationBetweenTwoPairsofNodes, XYZofNodes1stLayer, Data1stLayer);
        
        for (int nodeID = 0; nodeID < 4; ++nodeID)
        {
            int abcdID = RelationBetweenTwoPairsofNodes[nodeID];            
            xx(IJKofNodes2ndLayer[nodeID][0], IJKofNodes2ndLayer[nodeID][1], IJKofNodes2ndLayer[nodeID][2]) = Data2ndLayer[abcdID][0];
            yy(IJKofNodes2ndLayer[nodeID][0], IJKofNodes2ndLayer[nodeID][1], IJKofNodes2ndLayer[nodeID][2]) = Data2ndLayer[abcdID][1];
            zz(IJKofNodes2ndLayer[nodeID][0], IJKofNodes2ndLayer[nodeID][1], IJKofNodes2ndLayer[nodeID][2]) = Data2ndLayer[abcdID][2];
        }        
    }  
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder_mission2(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();
    
    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();
    
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];

        RDouble Data1stLayer[4][3];
        for (int nodeID = 0; nodeID <= 3; ++nodeID)
        {
            PHRead(dataContainer, Data1stLayer[nodeID][0]);
            PHRead(dataContainer, Data1stLayer[nodeID][1]);
            PHRead(dataContainer, Data1stLayer[nodeID][2]);
        }

        RDouble ab[3];
        RDouble ad[3];
        RDouble bc[3];
        RDouble ba[3];
        RDouble cd[3];
        RDouble cb[3];
        RDouble da[3];
        RDouble dc[3];

        PHRead(dataContainer, ab[0]);
        PHRead(dataContainer, ab[1]);
        PHRead(dataContainer, ab[2]);

        PHRead(dataContainer, ad[0]);
        PHRead(dataContainer, ad[1]);
        PHRead(dataContainer, ad[2]);

        PHRead(dataContainer, bc[0]);
        PHRead(dataContainer, bc[1]);
        PHRead(dataContainer, bc[2]);

        PHRead(dataContainer, ba[0]);
        PHRead(dataContainer, ba[1]);
        PHRead(dataContainer, ba[2]);

        PHRead(dataContainer, cd[0]);
        PHRead(dataContainer, cd[1]);
        PHRead(dataContainer, cd[2]);

        PHRead(dataContainer, cb[0]);
        PHRead(dataContainer, cb[1]);
        PHRead(dataContainer, cb[2]);

        PHRead(dataContainer, da[0]);
        PHRead(dataContainer, da[1]);
        PHRead(dataContainer, da[2]);

        PHRead(dataContainer, dc[0]);
        PHRead(dataContainer, dc[1]);
        PHRead(dataContainer, dc[2]);

        int IJKofNodes1stLayer[12][3];
        int IJKofNodes2ndLayer[12][3];
        GetNodeIJKfromTargetIndexghostofghost(finestGrid, iFace, 1, IJKofNodes1stLayer);
        GetNodeIJKfromTargetIndexghostofghost(finestGrid, iFace, 2, IJKofNodes2ndLayer);

        RDouble XYZofNodes1stLayer[4][3];
        for (int nodeID = 0; nodeID <= 3; ++nodeID)
        {
            XYZofNodes1stLayer[nodeID][0] = xx(IJKofNodes1stLayer[nodeID][0], IJKofNodes1stLayer[nodeID][1], IJKofNodes1stLayer[nodeID][2]);
            XYZofNodes1stLayer[nodeID][1] = yy(IJKofNodes1stLayer[nodeID][0], IJKofNodes1stLayer[nodeID][1], IJKofNodes1stLayer[nodeID][2]);
            XYZofNodes1stLayer[nodeID][2] = zz(IJKofNodes1stLayer[nodeID][0], IJKofNodes1stLayer[nodeID][1], IJKofNodes1stLayer[nodeID][2]);
        }
        
        int mapcase = GetRelationshipBetweenTwoPairsofNodes(XYZofNodes1stLayer, Data1stLayer);

        if (mapcase == 1)
        {
            xx(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = ab[0];
            yy(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = ab[1];
            zz(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = ab[2];

            xx(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = ad[0];
            yy(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = ad[1];
            zz(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = ad[2];

            xx(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = bc[0];
            yy(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = bc[1];
            zz(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = bc[2];

            xx(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = ba[0];
            yy(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = ba[1];
            zz(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = ba[2];

            xx(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = cd[0];
            yy(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = cd[1];
            zz(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = cd[2];

            xx(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = cb[0];
            yy(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = cb[1];
            zz(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = cb[2];

            xx(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = da[0];
            yy(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = da[1];
            zz(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = da[2];

            xx(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = dc[0];
            yy(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = dc[1];
            zz(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = dc[2];
        }
        else if (mapcase == 2)
        {
            xx(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = ad[0];
            yy(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = ad[1];
            zz(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = ad[2];

            xx(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = ab[0];
            yy(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = ab[1];
            zz(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = ab[2];

            xx(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = dc[0];
            yy(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = dc[1];
            zz(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = dc[2];

            xx(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = da[0];
            yy(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = da[1];
            zz(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = da[2];

            xx(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = cb[0];
            yy(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = cb[1];
            zz(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = cb[2];

            xx(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = cd[0];
            yy(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = cd[1];
            zz(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = cd[2];

            xx(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = ba[0];
            yy(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = ba[1];
            zz(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = ba[2];

            xx(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = bc[0];
            yy(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = bc[1];
            zz(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = bc[2];
        }
        else if (mapcase == 3)
        {
            xx(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = ba[0];
            yy(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = ba[1];
            zz(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = ba[2];

            xx(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = bc[0];
            yy(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = bc[1];
            zz(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = bc[2];

            xx(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = ad[0];
            yy(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = ad[1];
            zz(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = ad[2];

            xx(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = ab[0];
            yy(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = ab[1];
            zz(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = ab[2];

            xx(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = dc[0];
            yy(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = dc[1];
            zz(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = dc[2];

            xx(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = da[0];
            yy(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = da[1];
            zz(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = da[2];

            xx(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = cb[0];
            yy(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = cb[1];
            zz(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = cb[2];

            xx(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = cd[0];
            yy(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = cd[1];
            zz(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = cd[2];
        }
        else if (mapcase == 4)
        {
            xx(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = bc[0];
            yy(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = bc[1];
            zz(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = bc[2];

            xx(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = ba[0];
            yy(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = ba[1];
            zz(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = ba[2];

            xx(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = cd[0];
            yy(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = cd[1];
            zz(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = cd[2];

            xx(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = cb[0];
            yy(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = cb[1];
            zz(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = cb[2];

            xx(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = da[0];
            yy(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = da[1];
            zz(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = da[2];

            xx(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = dc[0];
            yy(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = dc[1];
            zz(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = dc[2];

            xx(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = ab[0];
            yy(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = ab[1];
            zz(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = ab[2];

            xx(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = ad[0];
            yy(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = ad[1];
            zz(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = ad[2];
        }
        else if (mapcase == 5)
        {
            xx(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = cb[0];
            yy(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = cb[1];
            zz(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = cb[2];

            xx(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = cd[0];
            yy(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = cd[1];
            zz(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = cd[2];

            xx(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = ba[0];
            yy(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = ba[1];
            zz(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = ba[2];

            xx(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = bc[0];
            yy(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = bc[1];
            zz(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = bc[2];

            xx(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = ad[0];
            yy(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = ad[1];
            zz(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = ad[2];

            xx(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = ab[0];
            yy(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = ab[1];
            zz(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = ab[2];

            xx(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = dc[0];
            yy(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = dc[1];
            zz(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = dc[2];

            xx(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = da[0];
            yy(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = da[1];
            zz(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = da[2];
        }
        else if (mapcase == 6)
        {
            xx(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = cd[0];
            yy(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = cd[1];
            zz(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = cd[2];

            xx(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = cb[0];
            yy(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = cb[1];
            zz(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = cb[2];

            xx(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = da[0];
            yy(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = da[1];
            zz(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = da[2];

            xx(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = dc[0];
            yy(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = dc[1];
            zz(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = dc[2];

            xx(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = ab[0];
            yy(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = ab[1];
            zz(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = ab[2];

            xx(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = ad[0];
            yy(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = ad[1];
            zz(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = ad[2];

            xx(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = bc[0];
            yy(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = bc[1];
            zz(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = bc[2];

            xx(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = ba[0];
            yy(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = ba[1];
            zz(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = ba[2];
        }
        else if (mapcase == 7)
        {
            xx(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = da[0];
            yy(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = da[1];
            zz(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = da[2];

            xx(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = dc[0];
            yy(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = dc[1];
            zz(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = dc[2];

            xx(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = ab[0];
            yy(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = ab[1];
            zz(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = ab[2];

            xx(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = ad[0];
            yy(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = ad[1];
            zz(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = ad[2];

            xx(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = bc[0];
            yy(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = bc[1];
            zz(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = bc[2];

            xx(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = ba[0];
            yy(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = ba[1];
            zz(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = ba[2];

            xx(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = cd[0];
            yy(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = cd[1];
            zz(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = cd[2];

            xx(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = cb[0];
            yy(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = cb[1];
            zz(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = cb[2];
        }
        else if (mapcase == 8)
        {
            xx(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = dc[0];
            yy(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = dc[1];
            zz(IJKofNodes2ndLayer[ 4][0], IJKofNodes2ndLayer[ 4][1], IJKofNodes2ndLayer[ 4][2]) = dc[2];

            xx(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = da[0];
            yy(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = da[1];
            zz(IJKofNodes2ndLayer[ 5][0], IJKofNodes2ndLayer[ 5][1], IJKofNodes2ndLayer[ 5][2]) = da[2];

            xx(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = cb[0];
            yy(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = cb[1];
            zz(IJKofNodes2ndLayer[ 6][0], IJKofNodes2ndLayer[ 6][1], IJKofNodes2ndLayer[ 6][2]) = cb[2];

            xx(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = cd[0];
            yy(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = cd[1];
            zz(IJKofNodes2ndLayer[ 7][0], IJKofNodes2ndLayer[ 7][1], IJKofNodes2ndLayer[ 7][2]) = cd[2];

            xx(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = ba[0];
            yy(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = ba[1];
            zz(IJKofNodes2ndLayer[ 8][0], IJKofNodes2ndLayer[ 8][1], IJKofNodes2ndLayer[ 8][2]) = ba[2];

            xx(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = bc[0];
            yy(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = bc[1];
            zz(IJKofNodes2ndLayer[ 9][0], IJKofNodes2ndLayer[ 9][1], IJKofNodes2ndLayer[ 9][2]) = bc[2];

            xx(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = ad[0];
            yy(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = ad[1];
            zz(IJKofNodes2ndLayer[10][0], IJKofNodes2ndLayer[10][1], IJKofNodes2ndLayer[10][2]) = ad[2];

            xx(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = ab[0];
            yy(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = ab[1];
            zz(IJKofNodes2ndLayer[11][0], IJKofNodes2ndLayer[11][1], IJKofNodes2ndLayer[11][2]) = ab[2];
        }
    }  
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder_mission3(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();
    
    RDouble4D &xfn  = *(this->GetFaceNormalX());
    RDouble4D &yfn  = *(this->GetFaceNormalY());
    RDouble4D &zfn  = *(this->GetFaceNormalZ());
    RDouble4D &area = *(this->GetFaceArea());
    RDouble4D &xfv  = *(this->GetFaceVectorX());
    RDouble4D &yfv  = *(this->GetFaceVectorY());
    RDouble4D &zfv  = *(this->GetFaceVectorZ());
    
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];

        for (int layer = 1; layer <= 2; ++ layer)
        {
            int it, jt, kt;
            finestGrid->GetTargetIndexIJK(iFace,layer, it, jt, kt);
            finestGrid->RemapMultigridIJK(level, it, jt, kt);

            RDouble KX, KY, KZ;

            PHRead(dataContainer, KX);
            PHRead(dataContainer, KY);
            PHRead(dataContainer, KZ);
            
            int id, jd, kd;
            int i_lr, j_lr, k_lr;
            finestGrid->GetTargetIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);
            if (i_lr != 0)
            {
                //! I_face
                if (i_lr == -1)
                {
                    //! left boundary                    
                    RDouble nx = - KX;
                    RDouble ny = - KY;
                    RDouble nz = - KZ;
                    
                    RDouble ds = sqrt(nx * nx + ny * ny + nz * nz);
                    
                    xfv (it, jt, kt, 1) = nx;
                    yfv (it, jt, kt, 1) = ny;
                    zfv (it, jt, kt, 1) = nz;
                    area(it, jt, kt, 1) = ds;
                    xfn (it, jt, kt, 1) = nx / (ds + TINY);
                    yfn (it, jt, kt, 1) = ny / (ds + TINY);
                    zfn (it, jt, kt, 1) = nz / (ds + TINY);
                }
                else
                {
                    //! right boundary
                    RDouble nx = KX;
                    RDouble ny = KY;
                    RDouble nz = KZ;
                    
                    RDouble ds = sqrt(nx * nx + ny * ny + nz * nz);
                    
                    xfv (it+1, jt, kt, 1) = nx;
                    yfv (it+1, jt, kt, 1) = ny;
                    zfv (it+1, jt, kt, 1) = nz;
                    area(it+1, jt, kt, 1) = ds;
                    xfn (it+1, jt, kt, 1) = nx / (ds + TINY);
                    yfn (it+1, jt, kt, 1) = ny / (ds + TINY);
                    zfn (it+1, jt, kt, 1) = nz / (ds + TINY);
                }
            }
            else if (j_lr != 0)
            {
                //! J_face
                if (j_lr == -1)
                {
                    //! left boundary                    
                    RDouble nx = - KX;
                    RDouble ny = - KY;
                    RDouble nz = - KZ;
                    
                    RDouble ds = sqrt(nx * nx + ny * ny + nz * nz);
                    
                    xfv (it, jt, kt, 2) = nx;
                    yfv (it, jt, kt, 2) = ny;
                    zfv (it, jt, kt, 2) = nz;
                    area(it, jt, kt, 2) = ds;
                    xfn (it, jt, kt, 2) = nx / (ds + TINY);
                    yfn (it, jt, kt, 2) = ny / (ds + TINY);
                    zfn (it, jt, kt, 2) = nz / (ds + TINY);
                }
                else
                {
                    //! right boundary
                    RDouble nx = KX;
                    RDouble ny = KY;
                    RDouble nz = KZ;
                    
                    RDouble ds = sqrt(nx * nx + ny * ny + nz * nz);
                    
                    xfv (it, jt+1, kt, 2) = nx;
                    yfv (it, jt+1, kt, 2) = ny;
                    zfv (it, jt+1, kt, 2) = nz;
                    area(it, jt+1, kt, 2) = ds;
                    xfn (it, jt+1, kt, 2) = nx / (ds + TINY);
                    yfn (it, jt+1, kt, 2) = ny / (ds + TINY);
                    zfn (it, jt+1, kt, 2) = nz / (ds + TINY);
                }
            }
            else
            {
                //! K_face
                if (k_lr == -1)
                {
                    //! left boundary                    
                    RDouble nx = - KX;
                    RDouble ny = - KY;
                    RDouble nz = - KZ;
                    
                    RDouble ds = sqrt(nx * nx + ny * ny + nz * nz);
                    
                    xfv (it, jt, kt, 3) = nx;
                    yfv (it, jt, kt, 3) = ny;
                    zfv (it, jt, kt, 3) = nz;
                    area(it, jt, kt, 3) = ds;
                    xfn (it, jt, kt, 3) = nx / (ds + TINY);
                    yfn (it, jt, kt, 3) = ny / (ds + TINY);
                    zfn (it, jt, kt, 3) = nz / (ds + TINY);
                }
                else
                {
                    //! right boundary
                    RDouble nx = KX;
                    RDouble ny = KY;
                    RDouble nz = KZ;
                    
                    RDouble ds = sqrt(nx * nx + ny * ny + nz * nz);
                    
                    xfv (it, jt, kt+1, 3) = nx;
                    yfv (it, jt, kt+1, 3) = ny;
                    zfv (it, jt, kt+1, 3) = nz;
                    area(it, jt, kt+1, 3) = ds;
                    xfn (it, jt, kt+1, 3) = nx / (ds + TINY);
                    yfn (it, jt, kt+1, 3) = ny / (ds + TINY);
                    zfn (it, jt, kt+1, 3) = nz / (ds + TINY);
                }
            }  
        }
    }  
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder_mission4(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();
    
    RDouble3D &vol = *this->GetCellVolume();
    
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];

        for (int layer = 1; layer <= 4; ++ layer)
        {
            int it, jt, kt;
            finestGrid->GetTargetIndexIJK(iFace,layer, it, jt, kt);
            finestGrid->RemapMultigridIJK(level, it, jt, kt);

            PHRead(dataContainer, vol(it, jt, kt));
        }
    }  
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder_mission5(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble4D &q_FD = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("q_FD"));

    int nEquation = 5;

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        for (int layer = 3; layer <= 4; ++layer)
        {
            int it, jt, kt;
            finestGrid->GetTargetIndexIJK(iFace, layer, it, jt, kt);
            finestGrid->RemapMultigridIJK(0, it, jt, kt);
            for (int m = 0; m < nEquation; ++ m)
            {
                PHRead(dataContainer, q_FD(it, jt, kt, m));
            }
        }
    }
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder_mission6(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (this->GetDataPtr("rtem"));

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        for (int layer = 1; layer <= 2; ++layer)
        {
            int it, jt, kt;
            finestGrid->GetTargetIndexIJK(iFace, layer, it, jt, kt);
            finestGrid->RemapMultigridIJK(0, it, jt, kt);
            PHRead(dataContainer, rtem(it, jt, kt));
        }
    }
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder_mission7(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (this->GetDataPtr("rtem"));
    RDouble4D &dqdxnokxi = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdxnokxi"));
    RDouble4D &dqdynokxi = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdynokxi"));
    RDouble4D &dqdznokxi = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdznokxi"));

    RDouble4D &dqdxnoeta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdxnoeta"));
    RDouble4D &dqdynoeta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdynoeta"));
    RDouble4D &dqdznoeta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdznoeta"));

    RDouble4D &dqdxnocta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdxnocta"));
    RDouble4D &dqdynocta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdynocta"));
    RDouble4D &dqdznocta = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("dqdznocta"));

    RDouble dqdx_layer1[4];
    RDouble dqdy_layer1[4];
    RDouble dqdz_layer1[4];

    RDouble dqdx_layer2[4];
    RDouble dqdy_layer2[4];
    RDouble dqdz_layer2[4];

    RDouble dqdx_layer3[4];
    RDouble dqdy_layer3[4];
    RDouble dqdz_layer3[4];

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        
        int it1, jt1, kt1; //! layer = 1;
        int it2, jt2, kt2; //! layer = 2;
        int it3, jt3, kt3; //! layer = 3;

        finestGrid->GetTargetIndexIJK(iFace, 1, it1, jt1, kt1);
        finestGrid->GetTargetIndexIJK(iFace, 2, it2, jt2, kt2);
        finestGrid->GetTargetIndexIJK(iFace, 3, it3, jt3, kt3);

        finestGrid->RemapMultigridIJK(0, it1, jt1, kt1);
        finestGrid->RemapMultigridIJK(0, it2, jt2, kt2);
        finestGrid->RemapMultigridIJK(0, it3, jt3, kt3);

        PHRead(dataContainer, rtem(it1, jt1, kt1));
        PHRead(dataContainer, rtem(it2, jt2, kt2));

        for (int m = 0; m < 4; ++ m)
        {
            PHRead(dataContainer, dqdx_layer1[m]);
            PHRead(dataContainer, dqdy_layer1[m]);
            PHRead(dataContainer, dqdz_layer1[m]);

            PHRead(dataContainer, dqdx_layer2[m]);
            PHRead(dataContainer, dqdy_layer2[m]);
            PHRead(dataContainer, dqdz_layer2[m]);

            PHRead(dataContainer, dqdx_layer3[m]);
            PHRead(dataContainer, dqdy_layer3[m]);
            PHRead(dataContainer, dqdz_layer3[m]);
        }
        
        int id, jd, kd;
        int i_lr, j_lr, k_lr;
        finestGrid->GetTargetIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);
        if (i_lr != 0)
        {
            //! I_face.
            for (int m = 0; m < 4; ++ m)
            {
                dqdxnokxi(it1, jt1, kt1, m) = dqdx_layer1[m];
                dqdynokxi(it1, jt1, kt1, m) = dqdy_layer1[m];
                dqdznokxi(it1, jt1, kt1, m) = dqdz_layer1[m];

                dqdxnokxi(it2, jt2, kt2, m) = dqdx_layer2[m];
                dqdynokxi(it2, jt2, kt2, m) = dqdy_layer2[m];
                dqdznokxi(it2, jt2, kt2, m) = dqdz_layer2[m];

                dqdxnokxi(it3, jt3, kt3, m) = dqdx_layer3[m];
                dqdynokxi(it3, jt3, kt3, m) = dqdy_layer3[m];
                dqdznokxi(it3, jt3, kt3, m) = dqdz_layer3[m];
            }
        }
        else if (j_lr != 0)
        {
            //! J_face.
            for (int m = 0; m < 4; ++ m)
            {
                dqdxnoeta(it1, jt1, kt1, m) = dqdx_layer1[m];
                dqdynoeta(it1, jt1, kt1, m) = dqdy_layer1[m];
                dqdznoeta(it1, jt1, kt1, m) = dqdz_layer1[m];

                dqdxnoeta(it2, jt2, kt2, m) = dqdx_layer2[m];
                dqdynoeta(it2, jt2, kt2, m) = dqdy_layer2[m];
                dqdznoeta(it2, jt2, kt2, m) = dqdz_layer2[m];

                dqdxnoeta(it3, jt3, kt3, m) = dqdx_layer3[m];
                dqdynoeta(it3, jt3, kt3, m) = dqdy_layer3[m];
                dqdznoeta(it3, jt3, kt3, m) = dqdz_layer3[m];
            }
        }
        else
        {
            //! K_face.
            for (int m = 0; m < 4; ++ m)
            {
                dqdxnocta(it1, jt1, kt1, m) = dqdx_layer1[m];
                dqdynocta(it1, jt1, kt1, m) = dqdy_layer1[m];
                dqdznocta(it1, jt1, kt1, m) = dqdz_layer1[m];

                dqdxnocta(it2, jt2, kt2, m) = dqdx_layer2[m];
                dqdynocta(it2, jt2, kt2, m) = dqdy_layer2[m];
                dqdznocta(it2, jt2, kt2, m) = dqdz_layer2[m];

                dqdxnocta(it3, jt3, kt3, m) = dqdx_layer3[m];
                dqdynocta(it3, jt3, kt3, m) = dqdy_layer3[m];
                dqdznocta(it3, jt3, kt3, m) = dqdz_layer3[m];
            }
        }
    }
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder_mission8(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (this->GetDataPtr("rtem"));
    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast <RDouble4D *> (this->GetDataPtr("gradUVWTCellCenterZ"));

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        
        int it1, jt1, kt1; //! layer = 1;
        int it2, jt2, kt2; //! layer = 2;
        int it3, jt3, kt3; //! layer = 3;

        finestGrid->GetTargetIndexIJK(iFace, 1, it1, jt1, kt1);
        finestGrid->GetTargetIndexIJK(iFace, 2, it2, jt2, kt2);
        finestGrid->GetTargetIndexIJK(iFace, 3, it3, jt3, kt3);

        finestGrid->RemapMultigridIJK(0, it1, jt1, kt1);
        finestGrid->RemapMultigridIJK(0, it2, jt2, kt2);
        finestGrid->RemapMultigridIJK(0, it3, jt3, kt3);

        PHRead(dataContainer, rtem(it1, jt1, kt1));
        PHRead(dataContainer, rtem(it2, jt2, kt2));

        PHRead(dataContainer, gradUVWTCellCenterX(it1, jt1, kt1,0));
        PHRead(dataContainer, gradUVWTCellCenterX(it1, jt1, kt1,1));
        PHRead(dataContainer, gradUVWTCellCenterX(it1, jt1, kt1,2));
        PHRead(dataContainer, gradUVWTCellCenterX(it1, jt1, kt1,3));

        PHRead(dataContainer, gradUVWTCellCenterX(it2, jt2, kt2,0));
        PHRead(dataContainer, gradUVWTCellCenterX(it2, jt2, kt2,1));
        PHRead(dataContainer, gradUVWTCellCenterX(it2, jt2, kt2,2));
        PHRead(dataContainer, gradUVWTCellCenterX(it2, jt2, kt2,3));

        PHRead(dataContainer, gradUVWTCellCenterX(it3, jt3, kt3,0));
        PHRead(dataContainer, gradUVWTCellCenterX(it3, jt3, kt3,1));
        PHRead(dataContainer, gradUVWTCellCenterX(it3, jt3, kt3,2));
        PHRead(dataContainer, gradUVWTCellCenterX(it3, jt3, kt3,3));

        PHRead(dataContainer, gradUVWTCellCenterY(it1, jt1, kt1,0));
        PHRead(dataContainer, gradUVWTCellCenterY(it1, jt1, kt1,1));
        PHRead(dataContainer, gradUVWTCellCenterY(it1, jt1, kt1,2));
        PHRead(dataContainer, gradUVWTCellCenterY(it1, jt1, kt1,3));

        PHRead(dataContainer, gradUVWTCellCenterY(it2, jt2, kt2,0));
        PHRead(dataContainer, gradUVWTCellCenterY(it2, jt2, kt2,1));
        PHRead(dataContainer, gradUVWTCellCenterY(it2, jt2, kt2,2));
        PHRead(dataContainer, gradUVWTCellCenterY(it2, jt2, kt2,3));

        PHRead(dataContainer, gradUVWTCellCenterY(it3, jt3, kt3,0));
        PHRead(dataContainer, gradUVWTCellCenterY(it3, jt3, kt3,1));
        PHRead(dataContainer, gradUVWTCellCenterY(it3, jt3, kt3,2));
        PHRead(dataContainer, gradUVWTCellCenterY(it3, jt3, kt3,3));

        PHRead(dataContainer, gradUVWTCellCenterZ(it1, jt1, kt1,0));
        PHRead(dataContainer, gradUVWTCellCenterZ(it1, jt1, kt1,1));
        PHRead(dataContainer, gradUVWTCellCenterZ(it1, jt1, kt1,2));
        PHRead(dataContainer, gradUVWTCellCenterZ(it1, jt1, kt1,3));

        PHRead(dataContainer, gradUVWTCellCenterZ(it2, jt2, kt2,0));
        PHRead(dataContainer, gradUVWTCellCenterZ(it2, jt2, kt2,1));
        PHRead(dataContainer, gradUVWTCellCenterZ(it2, jt2, kt2,2));
        PHRead(dataContainer, gradUVWTCellCenterZ(it2, jt2, kt2,3));

        PHRead(dataContainer, gradUVWTCellCenterZ(it3, jt3, kt3,0));
        PHRead(dataContainer, gradUVWTCellCenterZ(it3, jt3, kt3,1));
        PHRead(dataContainer, gradUVWTCellCenterZ(it3, jt3, kt3,2));
        PHRead(dataContainer, gradUVWTCellCenterZ(it3, jt3, kt3,3));
    }
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder_mission9(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &vist = *reinterpret_cast<RDouble3D *> (this->GetDataPtr("vist"));
    RDouble4D &gradTurbulenceCellCenterX = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterX"));
    RDouble4D &gradTurbulenceCellCenterY = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterY"));
    RDouble4D &gradTurbulenceCellCenterZ = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterZ"));

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        
        int it1, jt1, kt1; //! layer = 1;

        finestGrid->GetTargetIndexIJK(iFace, 1, it1, jt1, kt1);

        finestGrid->RemapMultigridIJK(0, it1, jt1, kt1);

        PHRead(dataContainer, vist(it1, jt1, kt1));

        PHRead(dataContainer, gradTurbulenceCellCenterX(it1, jt1, kt1,ISA));
        PHRead(dataContainer, gradTurbulenceCellCenterY(it1, jt1, kt1,ISA));
        PHRead(dataContainer, gradTurbulenceCellCenterZ(it1, jt1, kt1,ISA));
    }
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder_mission10(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &blend = *reinterpret_cast<RDouble3D *> (this->GetDataPtr("blend"));
    RDouble3D &vist = *reinterpret_cast<RDouble3D *> (this->GetDataPtr("vist"));
    RDouble4D &gradTurbulenceCellCenterX = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterX"));
    RDouble4D &gradTurbulenceCellCenterY = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterY"));
    RDouble4D &gradTurbulenceCellCenterZ = *reinterpret_cast<RDouble4D *> (this->GetDataPtr("gradTurbulenceCellCenterZ"));

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        
        int it1, jt1, kt1; //! layer = 1;

        finestGrid->GetTargetIndexIJK(iFace, 1, it1, jt1, kt1);

        finestGrid->RemapMultigridIJK(0, it1, jt1, kt1);

        PHRead(dataContainer, blend(it1, jt1, kt1));

        PHRead(dataContainer, vist(it1, jt1, kt1));

        PHRead(dataContainer, gradTurbulenceCellCenterX(it1, jt1, kt1,IKE));
        PHRead(dataContainer, gradTurbulenceCellCenterX(it1, jt1, kt1,IKW));

        PHRead(dataContainer, gradTurbulenceCellCenterY(it1, jt1, kt1,IKE));
        PHRead(dataContainer, gradTurbulenceCellCenterY(it1, jt1, kt1,IKW));

        PHRead(dataContainer, gradTurbulenceCellCenterZ(it1, jt1, kt1,IKE));
        PHRead(dataContainer, gradTurbulenceCellCenterZ(it1, jt1, kt1,IKW));
    }
}

void StructGrid::DownloadGridInfoOnlyforStructHighOrder_mission11(DataContainer *&dataContainer, const int &neighborZoneIndex)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(this->GetFinestGrid());
    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->MoveToBegin();

    RDouble3D &vist = *reinterpret_cast<RDouble3D *> (this->GetDataPtr("vist"));

    RDouble3D &subgridScaleEnergy = *reinterpret_cast<RDouble3D *> (this->GetDataPtr("subgridScaleEnergy"));

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        
        int it1, jt1, kt1; //! layer = 1;

        finestGrid->GetTargetIndexIJK(iFace, 1, it1, jt1, kt1);

        finestGrid->RemapMultigridIJK(0, it1, jt1, kt1);

        PHRead(dataContainer, vist(it1, jt1, kt1));

        PHRead(dataContainer, subgridScaleEnergy(it1, jt1, kt1));

        PHRead(dataContainer, turbulentPrandtlNumber(it1, jt1, kt1));
    }
}

}
