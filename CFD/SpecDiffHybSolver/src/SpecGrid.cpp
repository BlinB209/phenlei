#include "TK_Time.h"
#include "TK_Exit.h"
#include "PHMpi.h"
#include "PHHeader.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "SpecGrid.h"
#include "p3dfft.h"

using namespace std;

namespace PHSPACE
{
SpecGrid::SpecGrid()
{
}

SpecGrid::~SpecGrid()
{
    delete realKX; realKX = NULL;
    delete complexKX; complexKX = NULL;
    delete realKY; realKY = NULL;
    delete complexKY; complexKY = NULL;
    delete physicalStartIndex; physicalStartIndex = NULL;
    delete physicalEndIndex; physicalEndIndex = NULL;
    delete physicalIndexSize; physicalIndexSize = NULL;
    delete fourierStartIndex; fourierStartIndex = NULL;
    delete fourierEndIndex; fourierEndIndex = NULL; 
    delete fourierIndexSize; fourierIndexSize = NULL;
    delete dimensionOfP3DFFT; dimensionOfP3DFFT = NULL;
    delete startAndEndFourierIndexOfAll; startAndEndFourierIndexOfAll = NULL;
    delete startAndEndPhysicalIndexOfAll; startAndEndPhysicalIndexOfAll = NULL;
}

void SpecGrid::InitGridData()
{
    InitGridSize();

    InitGridFileName();

    AllocGridData();

    InitWavenumber();
}

void SpecGrid::InitGridSize()
{
    procIDWithKx0Ky0 = 0;
    nX = GlobalDataBase::GetIntParaFromDB("nx");
    nY = GlobalDataBase::GetIntParaFromDB("ny");
    nZ = GlobalDataBase::GetIntParaFromDB("nz");
    nX2 = 3 * nX / 2;
    nY2 = 3 * nY / 2;
    nZ2 = 3 * nZ / 2;
    nXYAll = nX2 * nY2;
    nXYZAll = nX2 * nY2 * nZ2;

    //		int numberofProcessor = PHMPI::GetNumberOfProcessor();
    Range I(1, 3);
    physicalStartIndex = new Int1D(I, fortranArray);
    physicalEndIndex   = new Int1D(I, fortranArray);
    physicalIndexSize  = new Int1D(I, fortranArray);
    fourierStartIndex  = new Int1D(I, fortranArray);
    fourierEndIndex    = new Int1D(I, fortranArray);
    fourierIndexSize   = new Int1D(I, fortranArray);
    dimensionOfP3DFFT  = new Int1D(Range(1, 2), fortranArray);
    GlobalDataBase::GetData("dimensionOfP3DFFT", &(*dimensionOfP3DFFT)(1), PHINT, 2);

    InitP3DFFT();

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int nTProcessor = GetNumberOfProcessor();

    // !To get startAndEndFourierIndexOfAll, write to griddistribute.dat
    startAndEndFourierIndexOfAll = new Int2D(Range(1, 6), Range(0, nTProcessor - 1), fortranArray);
    Int1D *startAndEndFourierIndex = new Int1D(Range(1, 6), fortranArray);
    (*startAndEndFourierIndex)(I) = (*fourierStartIndex)(I);
    (*startAndEndFourierIndex)(Range(4, 6)) = (*fourierEndIndex)(I);
#ifdef PH_PARALLEL
    MPI_Status status; 
    if(myid == server)
    {
#endif
        (*startAndEndFourierIndexOfAll)(Range(1, 6), myid) = (*startAndEndFourierIndex)(Range(1, 6));
#ifdef PH_PARALLEL
        for(int ip = 1; ip <= nTProcessor - 1; ++ip)
        {
            MPI_Recv(&(*startAndEndFourierIndexOfAll)(1, ip), 6, MPI_INTEGER, ip, ip, MPI_COMM_WORLD, &status);
        }
    }
    else
    {
        MPI_Send(&(*startAndEndFourierIndex)(1), 6, MPI_INTEGER, 0, myid, MPI_COMM_WORLD);
    }

    MPI_Bcast(&(*startAndEndFourierIndexOfAll)(1, 0), 6 * nTProcessor, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    int wordwidth = 20;
    string filename = "griddistribute.dat";
    fstream file;
    if( myid == server )
    {
        file.open(filename.c_str(), ios_base::out);
        if ( !file )
        {
            TK_Exit::ExceptionExit("could not open griddistribute.dat\n");
        }
        file << setw(wordwidth) << nTProcessor
            << setw(wordwidth) << (*dimensionOfP3DFFT)(1)
            << setw(wordwidth) << (*dimensionOfP3DFFT)(2) 
            << "\n";
        for(int ip = 0; ip <= nTProcessor - 1; ++ip)
        {
            file << setw(wordwidth) << ip 
                << setw(wordwidth) << (*startAndEndFourierIndexOfAll)(1, ip)
                << setw(wordwidth) << (*startAndEndFourierIndexOfAll)(2, ip)
                << setw(wordwidth) << (*startAndEndFourierIndexOfAll)(3, ip)
                << setw(wordwidth) << (*startAndEndFourierIndexOfAll)(4, ip)
                << setw(wordwidth) << (*startAndEndFourierIndexOfAll)(5, ip)
                << setw(wordwidth) << (*startAndEndFourierIndexOfAll)(6, ip)
                << "\n";
        }
    }
    file.close();
    file.clear();

    // !To get startAndEndPhysicalIndexOfAll, write to griddistribute_phy.dat
    startAndEndPhysicalIndexOfAll = new Int2D(Range(1, 6), Range(0, nTProcessor - 1), fortranArray);
    Int1D *startAndEndPhysicalIndex = new Int1D(Range(1, 6), fortranArray);
    (*startAndEndPhysicalIndex)(I) = (*physicalStartIndex)(I);
    (*startAndEndPhysicalIndex)(Range(4, 6)) = (*physicalEndIndex)(I);
#ifdef PH_PARALLEL
    if(myid == server)
    {
#endif
        (*startAndEndPhysicalIndexOfAll)(Range(1, 6), myid) = (*startAndEndPhysicalIndex)(Range(1, 6));
#ifdef PH_PARALLEL
        for(int ip = 1; ip <= nTProcessor - 1; ++ip)
        {
            MPI_Recv(&(*startAndEndPhysicalIndexOfAll)(1, ip), 6, MPI_INTEGER, ip, ip, MPI_COMM_WORLD, &status);
        }
    }
    else
    {
        MPI_Send(&(*startAndEndPhysicalIndex)(1), 6, MPI_INTEGER, 0, myid, MPI_COMM_WORLD);
    }

    MPI_Bcast(&(*startAndEndPhysicalIndexOfAll)(1, 0), 6 * nTProcessor, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    filename = "griddistribute_phy.dat";
    if( myid == server )
    {
        file.open( filename.c_str(), ios_base::out );
        if ( !file )
        {
            TK_Exit::ExceptionExit("could not open griddistribute_phy.dat\n");
        }
        file << setw(wordwidth) << nTProcessor 
            << setw(wordwidth) << (*dimensionOfP3DFFT)(1) 
            << setw(wordwidth) << (*dimensionOfP3DFFT)(2) 
            << "\n";
        for(int ip = 0; ip <= nTProcessor - 1; ++ip)
        {
            file << setw(wordwidth) << ip 
                << setw(wordwidth) << (*startAndEndPhysicalIndexOfAll)(1, ip)
                << setw(wordwidth) << (*startAndEndPhysicalIndexOfAll)(2, ip)
                << setw(wordwidth) << (*startAndEndPhysicalIndexOfAll)(3, ip)
                << setw(wordwidth) << (*startAndEndPhysicalIndexOfAll)(4, ip)
                << setw(wordwidth) << (*startAndEndPhysicalIndexOfAll)(5, ip)
                << setw(wordwidth) << (*startAndEndPhysicalIndexOfAll)(6, ip)
                << "\n";
        }
    }
    file.close();
    file.clear();

    nXYZFourier  = (*fourierIndexSize )(1) * (*fourierIndexSize )(2) * (*fourierIndexSize )(3);
    nXYZPhysical = (*physicalIndexSize)(1) * (*physicalIndexSize)(2) * (*physicalIndexSize)(3);

    iStartFourierIndex = (*fourierStartIndex)(1);
    jStartFourierIndex = (*fourierStartIndex)(2);
    kStartFourierIndex = (*fourierStartIndex)(3);
    iEndFourierIndex   = (*fourierEndIndex  )(1);
    jEndFourierIndex   = (*fourierEndIndex  )(2);
    kEndFourierIndex   = (*fourierEndIndex  )(3);

    iStartPhysicalIndex = (*physicalStartIndex)(1);
    jStartPhysicalIndex = (*physicalStartIndex)(2);
    kStartPhysicalIndex = (*physicalStartIndex)(3);
    iEndPhysicalIndex   = (*physicalEndIndex  )(1);
    jEndPhysicalIndex   = (*physicalEndIndex  )(2);
    kEndPhysicalIndex   = (*physicalEndIndex  )(3);

    Int1D *fourierStartIndexGlobal = new Int1D(Range(1, 3), fortranArray);
    Int1D *fourierEndIndexGlobal = new Int1D(Range(1, 3), fortranArray);

    (*fourierStartIndexGlobal)(1) = iStartFourierIndex;
    (*fourierStartIndexGlobal)(2) = jStartFourierIndex;
    (*fourierStartIndexGlobal)(3) = kStartFourierIndex;
    (*fourierEndIndexGlobal)(1) = iEndFourierIndex;
    (*fourierEndIndexGlobal)(2) = jEndFourierIndex;
    (*fourierEndIndexGlobal)(3) = kEndFourierIndex;

#ifdef PH_PARALLEL
    using namespace PHMPI;
    MPI_Reduce(&iStartFourierIndex, &(*fourierStartIndexGlobal)(1), 1, MPI_INTEGER, MPI_MIN, procIDWithKx0Ky0, MPI_COMM_WORLD);
    MPI_Reduce(&jStartFourierIndex, &(*fourierStartIndexGlobal)(2), 1, MPI_INTEGER, MPI_MIN, procIDWithKx0Ky0, MPI_COMM_WORLD);
    MPI_Reduce(&kStartFourierIndex, &(*fourierStartIndexGlobal)(3), 1, MPI_INTEGER, MPI_MIN, procIDWithKx0Ky0, MPI_COMM_WORLD);

    MPI_Reduce(&iEndFourierIndex, &(*fourierEndIndexGlobal)(1), 1, MPI_INTEGER, MPI_MAX, procIDWithKx0Ky0, MPI_COMM_WORLD);
    MPI_Reduce(&jEndFourierIndex, &(*fourierEndIndexGlobal)(2), 1, MPI_INTEGER, MPI_MAX, procIDWithKx0Ky0, MPI_COMM_WORLD);
    MPI_Reduce(&kEndFourierIndex, &(*fourierEndIndexGlobal)(3), 1, MPI_INTEGER, MPI_MAX, procIDWithKx0Ky0, MPI_COMM_WORLD);

    MPI_Bcast(&(*fourierStartIndexGlobal)(1) , 3, MPI_INTEGER, procIDWithKx0Ky0, MPI_COMM_WORLD);
    MPI_Bcast(&(*fourierEndIndexGlobal)(1), 3, MPI_INTEGER, procIDWithKx0Ky0, MPI_COMM_WORLD);
#endif

    iStartFourierIndexGlobal = (*fourierStartIndexGlobal)(1);
    jStartFourierIndexGlobal = (*fourierStartIndexGlobal)(2);
    kStartFourierIndexGlobal = (*fourierStartIndexGlobal)(3);

    iEndFourierIndexGlobal = (*fourierEndIndexGlobal)(1);
    jEndFourierIndexGlobal = (*fourierEndIndexGlobal)(2);
    kEndFourierIndexGlobal = (*fourierEndIndexGlobal)(3);

    nyquistX = kStartFourierIndexGlobal + nX / 2;
    nyquistY = jStartFourierIndexGlobal + nY / 2;

    //if ( (jStartFourierIndex == jStartFourierIndexGlobal) && (kStartFourierIndex == kStartFourierIndexGlobal) )
    //{
    //    procIDWithKx0Ky0 = myid;
    //}

    delete startAndEndFourierIndex; startAndEndFourierIndex = NULL;
    delete startAndEndPhysicalIndex; startAndEndPhysicalIndex = NULL;
    delete fourierStartIndexGlobal; fourierStartIndexGlobal = NULL;
    delete fourierEndIndexGlobal; fourierEndIndexGlobal = NULL;
}

void SpecGrid::InitGridFileName()
{
    string outputdir = "./results";
    outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");

    gridFileName = outputdir + "/grid.x"; 
    zgridFileName = outputdir + "/rXYZ.dat";
}

void SpecGrid::InitP3DFFT()
{
    Int1D *memsize = new Int1D(Range(1, 3), fortranArray);

    nX = GlobalDataBase::GetIntParaFromDB("nx");
    nY = GlobalDataBase::GetIntParaFromDB("ny");
    nZ = GlobalDataBase::GetIntParaFromDB("nz");

    (*fourierStartIndex)(1) = 1;
    (*fourierStartIndex)(2) = 1;
    (*fourierStartIndex)(3) = 1;
    
    (*fourierEndIndex)(1) = nZ;
    (*fourierEndIndex)(2) = nY;
    (*fourierEndIndex)(3) = nX / 2 + 1;

    (*physicalStartIndex)(1) = 1;
    (*physicalStartIndex)(2) = 1;
    (*physicalStartIndex)(3) = 1;

    (*physicalEndIndex)(1) = 3 * nX / 2;
    (*physicalEndIndex)(2) = 3 * nY / 2;
    (*physicalEndIndex)(3) = 3 * nZ / 2;

    Cp3dfft_setup(&(*dimensionOfP3DFFT)(1), nX2, nY2, nZ2, MPI_COMM_WORLD, nX, nY, nZ, 1, &(*memsize)(1));

    Cp3dfft_get_dims(&(*physicalStartIndex)(1), &(*physicalEndIndex)(1), &(*physicalIndexSize)(1), 1);

    Cp3dfft_get_dims(&(*fourierStartIndex)(1), &(*fourierEndIndex)(1), &(*fourierIndexSize)(1), 2);

    delete memsize; memsize = NULL;
}

void SpecGrid::AllocGridData()
{
    Range K(kStartFourierIndex, kEndFourierIndex);
    Range J(jStartFourierIndex, jEndFourierIndex);
    Range I(iStartFourierIndex, iEndFourierIndex);

    realKX = new RDouble1D(K, fortranArray);
    realKY = new RDouble1D(J, fortranArray);
    realKZ = new RDouble1D(I, fortranArray);

    complexKX = new Complex1D(K, fortranArray);
    complexKY = new Complex1D(J, fortranArray);
    complexKZ = new Complex1D(I, fortranArray);
}

void SpecGrid::InitWavenumber()
{
    RDouble realXL = GlobalDataBase::GetDoubleParaFromDB("XL");
    RDouble realAlphaX = 1.0 / realXL;
    for (int i3 = kStartFourierIndex; i3 <= kEndFourierIndex; ++i3)
    {
        int ix = i3 - kStartFourierIndexGlobal;
        (*realKX)(i3) = static_cast<double>(ix) * realAlphaX;
        (*complexKX)(i3) = PHComplex(0.0, 1.0) * (*realKX)(i3);
    }
    realKNyquistX = static_cast<double>(nX) * realAlphaX / 2.0;

    RDouble realYL = GlobalDataBase::GetDoubleParaFromDB("YL");
    RDouble realAlphaY = 1.0 / realYL;
    for ( int i2 = jStartFourierIndex; i2 <= jEndFourierIndex; ++i2 )
    {
        int iy = i2 - jStartFourierIndexGlobal;
        if ( iy <= nY/2 )
        {
            (*realKY)(i2) = static_cast<double>(iy) * realAlphaY;
        }
        else
        {
            (*realKY)(i2) = static_cast<double>(iy - nY) * realAlphaY;
        }
        (*complexKY)(i2) = PHComplex(0.0, 1.0) * (*realKY)(i2);
    }
    realKNyquistY = static_cast<double>(nY) * realAlphaY / 2.0;

    RDouble realZL = GlobalDataBase::GetDoubleParaFromDB("ZL");
    RDouble realAlphaZ = 1.0 / realZL;
    for ( int i1 = iStartFourierIndex; i1 <= iEndFourierIndex; ++i1 )
    {
        int iz = i1 - iStartFourierIndexGlobal;
        if (iz <= nZ/2)
        {
            (*realKZ)(i1) = static_cast<double>(iz) * realAlphaZ;
        }
        else
        {
            (*realKZ)(i1) = static_cast<double>(iz - nZ) * realAlphaZ;
        }
        (*complexKZ)(i1) = PHComplex(0.0, 1.0) * (*realKZ)(i1);
    }
    realKNyquistZ = static_cast<double>(nZ) * realAlphaZ / 2.0;
}

int SpecGrid::GetProcIDWithKx0Ky0()
{
    return procIDWithKx0Ky0;
}

}

