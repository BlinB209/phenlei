#include "SpecDiffHybFieldView.h"
#include "TK_Exit.h"
#include "PHHeader.h"
#include "Complex.h"
#include "fftw3.h"
#include "GlobalDataBase.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include "TK_Parse.h"

namespace PHSPACE
{

	SpecDiffHybFieldView::SpecDiffHybFieldView()
	{
	}

	SpecDiffHybFieldView::~SpecDiffHybFieldView()
	{
        delete realX; realX = NULL;
        delete realY; realY = NULL;
        delete realZ; realZ = NULL;

        delete realU; realU = NULL;
        delete realV; realV = NULL;
        delete realW; realW = NULL;
        delete realP; realP = NULL;

        delete complexU; complexU = NULL;
        delete complexV; complexV = NULL;
        delete complexW; complexW = NULL;
        delete complexP; complexP = NULL;

        delete realOut; realOut = NULL;
        delete complexIn; complexIn = NULL;

        delete fourierStartIndex; fourierStartIndex = NULL;
        delete fourierEndIndex; fourierEndIndex = NULL;
	}

	void SpecDiffHybFieldView::Run()
	{
        Init();

        InputGridDistribute();

        InputVelocity();

        InputPressure();

        FieldView();
	}

	void SpecDiffHybFieldView::Init()
	{
        nx = GlobalDataBase::GetIntParaFromDB("nx");
        ny = GlobalDataBase::GetIntParaFromDB("ny");
        nz = GlobalDataBase::GetIntParaFromDB("nz");

        AllocFieldData();
	}

    void SpecDiffHybFieldView::AllocFieldData()
    {
        Range IP(1, nx);
        Range JP(1, ny);
        Range KP(1, nz+1);

        Range IF(1, nz+1);
        Range JF(1, ny);
        Range KF(1, nx/2+1);

        realX = new RDouble1D(IP, fortranArray);
        realY = new RDouble1D(JP, fortranArray);
        realZ = new RDouble1D(KP, fortranArray);

        realU = new RDouble3D(IP, JP, KP, fortranArray);
        realV = new RDouble3D(IP, JP, KP, fortranArray);
        realW = new RDouble3D(IP, JP, KP, fortranArray);
        realP = new RDouble3D(IP, JP, KP, fortranArray);

        complexU = new Complex3D(IF, JF, KF, fortranArray);
        complexV = new Complex3D(IF, JF, KF, fortranArray);
        complexW = new Complex3D(IF, JF, KF, fortranArray);
        complexP = new Complex3D(IF, JF, KF, fortranArray);

        realOut = new RDouble2D(IP, JP, fortranArray);
        complexIn = new Complex2D(KF, JF, fortranArray);
    }

    void SpecDiffHybFieldView::InputGridDistribute()
    {
        string filename = "griddistribute.dat";
        fstream file;
        string line,word;
        string separator  = " =\t\r\n#$,;\"'";

        file.open(filename.c_str(), ios_base::in);
        if ( !file )
        {
            TK_Exit::ExceptionExit("could not open griddistribute.dat\n");
        }

        ReadNewLine(file, line);
        line = FindNextWord(line, word, separator);
        from_string<int>(nTProcessor, word, std::dec);

        fourierStartIndex = new Int2D(Range(1, 3), Range(0, nTProcessor-1), fortranArray);
        fourierEndIndex = new Int2D(Range(1, 3), Range(0, nTProcessor-1), fortranArray);
        
        int inttmp;
        for(int ip = 0; ip <= nTProcessor - 1; ++ip)
        {
            ReadNewLine(file, line);

            line = FindNextWord(line, word, separator);
            from_string<int>(inttmp, word, std::dec);
           
            line = FindNextWord(line, word, separator);
            from_string<int>((*fourierStartIndex)(1, ip), word, std::dec);
            
            line = FindNextWord(line, word, separator);
            from_string<int>((*fourierStartIndex)(2, ip), word, std::dec);
            
            line = FindNextWord(line, word, separator);
            from_string<int>((*fourierStartIndex)(3, ip), word, std::dec);
            
            line = FindNextWord(line, word, separator);
            from_string<int>((*fourierEndIndex)(1, ip), word, std::dec);
            
            line = FindNextWord(line, word, separator);
            from_string<int>((*fourierEndIndex)(2, ip), word, std::dec);
            
            line = FindNextWord(line, word, separator);
            from_string<int>((*fourierEndIndex)(3, ip), word, std::dec);
        }

        file.close();
        file.clear();
    }

    void SpecDiffHybFieldView::InputVelocity()
    {
        string outputdir = "./results";
        outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");

        for(int ip = 0; ip <= nTProcessor - 1; ++ip)
        {
            int istf = (*fourierStartIndex)(1, ip);
            int jstf = (*fourierStartIndex)(2, ip);
            int kstf = (*fourierStartIndex)(3, ip);
            int iedf = (*fourierEndIndex)(1, ip);
            int jedf = (*fourierEndIndex)(2, ip);
            int kedf = (*fourierEndIndex)(3, ip);

            stringstream strmip;
            strmip << ip;
            string filename = outputdir + "/Velocity_" + strmip.str() + ".rsta";
            fstream file;
            file.open(filename.c_str(), ios_base::in|ios_base::binary);
            if ( !file )
            {
                TK_Exit::FileOpenErrorExit(filename);
            }

            int inttmp[5];
            file.read(reinterpret_cast<char*>(inttmp), 5*sizeof(int));

            for(int ix = kstf; ix <= kedf; ++ix)
            {
                for(int iy = jstf; iy <= jedf; ++iy)
                {
                    for(int iz = istf; iz <= iedf; ++iz)
                    {
                        file.read(reinterpret_cast<char*>(&(*complexU)(iz, iy, ix)), sizeof(PHComplex));
                    }
                }
            }

            for(int ix = kstf; ix <= kedf; ++ix)
            {
                for(int iy = jstf; iy <= jedf; ++iy)
                {
                    for(int iz = istf; iz <= iedf; ++iz)
                    {
                        file.read(reinterpret_cast<char*>(&(*complexV)(iz, iy, ix)), sizeof(PHComplex));
                    }
                }
            }

            for(int ix = kstf; ix <= kedf; ++ix)
            {
                for(int iy = jstf; iy <= jedf; ++iy)
                {
                    for(int iz = istf; iz <= iedf; ++iz)
                    {
                        file.read(reinterpret_cast<char*>(&(*complexW)(iz, iy, ix)), sizeof(PHComplex));
                    }
                }
            }

            file.close();
            file.clear();
        }
    }

    void SpecDiffHybFieldView::InputPressure()
    {
        string outputdir = "./results";
        outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");

        for(int ip = 0; ip <= nTProcessor - 1; ++ip)
        {
            int istf = (*fourierStartIndex)(1, ip);
            int jstf = (*fourierStartIndex)(2, ip);
            int kstf = (*fourierStartIndex)(3, ip);
            int iedf = (*fourierEndIndex)(1, ip);
            int jedf = (*fourierEndIndex)(2, ip);
            int kedf = (*fourierEndIndex)(3, ip);

            stringstream strmip;
            strmip << ip;
            string filename = outputdir + "/Pressure_" + strmip.str() + ".rsta";
            fstream file;
            file.open(filename.c_str(), ios_base::in|ios_base::binary);
            if ( !file )
            {
                TK_Exit::FileOpenErrorExit(filename);
            }

            int inttmp[5];
            file.read(reinterpret_cast<char*>(inttmp), 5*sizeof(int));

            for(int ix = kstf; ix <= kedf; ++ix)
            {
                for(int iy = jstf; iy <= jedf; ++iy)
                {
                    for(int iz = istf; iz <= iedf; ++iz)
                    {
                        file.read(reinterpret_cast<char*>(&(*complexP)(iz, iy, ix)), sizeof(PHComplex));
                    }
                }
            }

            file.close();
            file.clear();
        }
    }

    void SpecDiffHybFieldView::FieldView()
    {
        fftw_plan iplan = fftw_plan_dft_c2r_2d(ny, nx, reinterpret_cast<fftw_complex*>(&(*complexIn)(1, 1)), &(*realOut)(1, 1), FFTW_ESTIMATE);

        for(int iz = 1; iz <= (nz+1); ++iz)
        {
            for(int ix = 1; ix <= (nx/2+1); ++ix)
            {
                for(int iy = 1; iy <= ny; ++iy)
                {
                    (*complexIn)(ix, iy) = (*complexU)(iz, iy, ix);
                }
            }

            /*test fieldview
            stringstream strmiz;
            strmiz << iz;
            string filename = "UIn_" + strmiz.str() + ".dat";
            fstream file;
            file.open(filename.c_str(), ios_base::out);
            
            for(int ix = 1; ix <= (nx/2+1); ++ix)
            {
                for(int iy = 1; iy <= ny; ++iy)
                {
                    file << setw(20) << ix
                         << setw(20) << iy
                         << setw(20) << (*complexIn)(ix, iy)
                         << "\n";
                }
            }
            file.close();
            file.clear();
            //test fieldview*/

            fftw_execute_dft_c2r(iplan, reinterpret_cast<fftw_complex*>(&(*complexIn)(1, 1)), &(*realOut)(1, 1));

            /*test fieldview
            filename = "UOut_" + strmiz.str() + ".dat";
            file.open(filename.c_str(), ios_base::out);
            
            for(int ix = 1; ix <= nx; ++ix)
            {
                for(int iy = 1; iy <= ny; ++iy)
                {
                    file << setw(20) << ix
                         << setw(20) << iy
                         << setw(20) << (*realOut)(ix, iy)
                         << "\n";
                }
            }
            file.close();
            file.clear();
            //test fieldview*/

            for(int ix = 1; ix <= nx; ++ix)
            {
                for(int iy = 1; iy <= ny; ++iy)
                {
                    (*realU)(ix, iy, iz) = (*realOut)(ix, iy);
                }
            }

            for(int ix = 1; ix <= (nx/2+1); ++ix)
            {
                for(int iy = 1; iy <= ny; ++iy)
                {
                    (*complexIn)(ix, iy) = (*complexV)(iz, iy, ix);
                }
            }
            fftw_execute_dft_c2r(iplan, reinterpret_cast<fftw_complex*>(&(*complexIn)(1, 1)), &(*realOut)(1, 1));
            for(int ix = 1; ix <= nx; ++ix)
            {
                for(int iy = 1; iy <= ny; ++iy)
                {
                    (*realV)(ix, iy, iz) = (*realOut)(ix, iy);
                }
            }

            for(int ix = 1; ix <= (nx/2+1); ++ix)
            {
                for(int iy = 1; iy <= ny; ++iy)
                {
                    (*complexIn)(ix, iy) = (*complexW)(iz, iy, ix);
                }
            }
            fftw_execute_dft_c2r(iplan, reinterpret_cast<fftw_complex*>(&(*complexIn)(1, 1)), &(*realOut)(1, 1));
            for(int ix = 1; ix <= nx; ++ix)
            {
                for(int iy = 1; iy <= ny; ++iy)
                {
                    (*realW)(ix, iy, iz) = (*realOut)(ix, iy);
                }
            }

            for(int ix = 1; ix <= (nx/2+1); ++ix)
            {
                for(int iy = 1; iy <= ny; ++iy)
                {
                    (*complexIn)(ix, iy) = (*complexP)(iz, iy, ix);
                }
            }
            fftw_execute_dft_c2r(iplan, reinterpret_cast<fftw_complex*>(&(*complexIn)(1, 1)), &(*realOut)(1, 1));
            for(int ix = 1; ix <= nx; ++ix)
            {
                for(int iy = 1; iy <= ny; ++iy)
                {
                    (*realP)(ix, iy, iz) = (*realOut)(ix, iy);
                }
            }
        } //!for iz
 
        /*test fieldview
        string testgridFileName = "testUVWP.dat";
        fstream file;
        file.open(testgridFileName.c_str(), ios_base::out);
        if (!file)
        {
            TK_Exit::FileOpenErrorExit(testgridFileName);
        }
        file << "P, U, V, W = \n";
        for(int iz = 1; iz <= (nz+1); ++iz)
        {
            for(int iy = 1; iy <= ny; ++iy)
            {
                for(int ix = 1; ix <= nx; ++ix)
                {
                    file << setw(20) << iz
                         << setw(20) << iy
                         << setw(20) << ix
                         << setw(20) << (*realP)(ix, iy, iz)
                         << setw(20) << (*realU)(ix, iy, iz)
                         << setw(20) << (*realV)(ix, iy, iz)
                         << setw(20) << (*realW)(ix, iy, iz)
                         << "\n";
                }
            }
        }
        file.close();
        file.clear();
        //test fieldview*/

        OutputGrid();

        OutputField();

        fftw_destroy_plan(iplan);
    }

    void SpecDiffHybFieldView::OutputGrid()
    {
        const double pi = 2.0 * acos(0.0);

        double XL = GlobalDataBase::GetDoubleParaFromDB("XL");
        double YL = GlobalDataBase::GetDoubleParaFromDB("YL");

        for(int ix = 1; ix <= nx; ++ix)
        {
            (*realX)(ix) = XL * 2.0 * pi / static_cast<double>(nx) * static_cast<double>(ix-1);
        }

        for(int iy = 1; iy <= ny; ++iy)
        {
            (*realY)(iy) = YL * 2.0 * pi / static_cast<double>(ny) * static_cast<double>(iy-1);
        }

        string outputdir = "./results";
        outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");
        string zgridFileName = outputdir + "/rZ.dat";
        fstream file;
        file.open(zgridFileName.c_str(), ios_base::in);
        if (!file)
        {
            TK_Exit::FileOpenErrorExit(zgridFileName);
        }

        int inttmp;
        string line,word;
        string separator  = " =\t\r\n#$,;\"'";
        for(int iz = 1; iz <= nz+1; ++iz)
        {
            ReadNewLine(file, line);
            line = FindNextWord(line, word, separator);
            from_string<int>(inttmp, word, std::dec);
            line = FindNextWord(line, word, separator);
            from_string<double>((*realZ)(iz), word, std::dec);
        }
        file.close();
        file.clear();

        string gridFileName = "grid.x";
        file.open(gridFileName.c_str(), ios_base::out|ios_base::binary);
        if (!file)
        {
            TK_Exit::FileOpenErrorExit(zgridFileName);
        }
        inttmp = 1;
        file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
        file.write(reinterpret_cast<char*>(&nx), sizeof(int));
        file.write(reinterpret_cast<char*>(&ny), sizeof(int));
        inttmp = nz + 1;
        file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));

        for(int iz = 1; iz <= nz+1; ++iz)
        {
            for(int iy = 1; iy <= ny; ++iy)
            {
                for(int ix = 1; ix <= nx; ++ix)
                {
                    file.write(reinterpret_cast<char*>(&(*realX)(ix)), sizeof(double));
                }
            }
        }
        
        for(int iz = 1; iz <= nz+1; ++iz)
        {
            for(int iy = 1; iy <= ny; ++iy)
            {
                for(int ix = 1; ix <= nx; ++ix)
                {
                    file.write(reinterpret_cast<char*>(&(*realY)(iy)), sizeof(double));
                }
            }
        }

        for(int iz = 1; iz <= nz+1; ++iz)
        {
            for(int iy = 1; iy <= ny; ++iy)
            {
                for(int ix = 1; ix <= nx; ++ix)
                {
                    file.write(reinterpret_cast<char*>(&(*realZ)(iz)), sizeof(double));
                }
            }
        }
        file.close();
        file.clear();

        /*test fieldview
        string testgridFileName = "testgrid.x";
        file.open(testgridFileName.c_str(), ios_base::out);
        if (!file)
        {
            TK_Exit::FileOpenErrorExit(zgridFileName);
        }
        file << setw(20) << 1 << "\n";
        file << setw(20) << nx 
             << setw(20) << ny
             << setw(20) << nz+1
             << "\n";
        for(int ix = 1; ix <= nx; ++ix)
        {
            file << setw(20) << ix << setw(20) << (*realX)(ix) << "\n";
        }
        for(int iy = 1; iy <= ny; ++iy)
        {
            file << setw(20) << iy << setw(20) << (*realY)(iy) << "\n";
        }
        for(int iz = 1; iz <= nz+1; ++iz)
        {
            file << setw(20) << iz << setw(20) << (*realZ)(iz) << "\n";
        }
        file.close();
        file.clear();
        //test fieldview*/
    }

    void SpecDiffHybFieldView::OutputField()
    {
        string fieldFileName = "fields.f";
        fstream file;
        file.open(fieldFileName.c_str(), ios_base::out|ios_base::binary);
        if (!file)
        {
            TK_Exit::FileOpenErrorExit(fieldFileName);
        }

        int inttmp = 1;
        file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
        file.write(reinterpret_cast<char*>(&nx), sizeof(int));
        file.write(reinterpret_cast<char*>(&ny), sizeof(int));
        inttmp = nz + 1;
        file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
        inttmp = 4;
        file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
        int nout = nx * ny * (nz + 1);
        file.write(reinterpret_cast<char*>(&(*realP)(1, 1, 1)), nout * sizeof(double));
        file.write(reinterpret_cast<char*>(&(*realU)(1, 1, 1)), nout * sizeof(double));
        file.write(reinterpret_cast<char*>(&(*realV)(1, 1, 1)), nout * sizeof(double));
        file.write(reinterpret_cast<char*>(&(*realW)(1, 1, 1)), nout * sizeof(double));
        file.close();
        file.clear();
    }
}