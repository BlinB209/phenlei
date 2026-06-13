#include "GlobalDataBase.h"
#include "Geo_FaceMetrics_Struct.h"

namespace PHSPACE
{

    Geo_FaceMetrics_Struct::Geo_FaceMetrics_Struct()
    {
        xfn = NULL;
        yfn = NULL;
        zfn = NULL;
        area = NULL;
        xFaceVector = NULL;
        yFaceVector = NULL;
        zFaceVector = NULL;

        xyzFaceVector_FD = NULL;
        xyzFaceNormal_FD = NULL;
    }

    Geo_FaceMetrics_Struct::~Geo_FaceMetrics_Struct()
    {
        DeallocateAll();
    }

    void Geo_FaceMetrics_Struct::AllocateMetrics(Range &I, Range &J, Range &K, Range &D)
    {
        DeallocateAll();

        xfn  = new RDouble4D(I, J, K, D, fortranArray);
        yfn  = new RDouble4D(I, J, K, D, fortranArray);
        zfn  = new RDouble4D(I, J, K, D, fortranArray);
        area = new RDouble4D(I, J, K, D, fortranArray);

        xFaceVector = new RDouble4D(I, J, K, D, fortranArray);
        yFaceVector = new RDouble4D(I, J, K, D, fortranArray);
        zFaceVector = new RDouble4D(I, J, K, D, fortranArray);

        int SolverStructOrder = PHSPACE::GlobalDataBase::GetIntParaFromDB("SolverStructOrder");
        if (SolverStructOrder > 2)
        {
            Range DA(D.first() - 1, D.last());

            xyzFaceVector_FD = new RDouble5D(D , I, J, K, D, fortranArray);
            xyzFaceNormal_FD = new RDouble5D(DA, I, J, K, D, fortranArray);
        }
    }
}