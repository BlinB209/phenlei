#include "GridGroup.h"
#include "TK_Exit.h"
#include "GridType.h"
#include "PHIO.h"
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{

GridGroup::GridGroup(int zoneStart)
{
    block_proc_dump = 0; // Bell 20131124 add
    nzones     = 0;
    block_proc = 0;
    block_proc_grid = 0;
    block_type = 0;
    block_idx  = 0;
    block_fileProc = 0;
    this->zoneStart = zoneStart;
    region     = 0;
}

GridGroup::~GridGroup()
{
    delete [] block_proc_dump; // Bell 20131124 add
    delete [] block_proc;
    delete [] block_proc_grid;
    delete [] block_type;
    delete [] block_idx ;
    delete [] block_fileProc;
}

void GridGroup::InitZoneLayoutInformation(const string &filename)
{
    fstream file;
    OpenSerialFile(file, filename, ios_base::in|ios_base::binary);

    InitZoneLayout(file);
    SetGlobalZoneLayout();

    ParallelCloseFile(file);
}

void GridGroup::SetGlobalZoneLayout()
{
    using namespace PHMPI;

    if (!GetZoneGridType())
    {
        SetNumberOfGlobalZones(nzones);
        CreateZoneProcessorIDDump(nzones);
        CreateZoneProcessorID(nzones);
        CreateZoneGridID     (nzones);
        CreateZoneGridType   (nzones);
        CreateZoneFileID     (nzones);
    }

    int * global_block_proc      = GetZoneProcessorID();
    int * global_block_proc_grid = GetZoneProcessorIDForGrid();
    int * global_block_type      = GetZoneGridType();
    int * global_block_idx       = GetZoneGridID();
    int * global_block_file      = GetZoneFileID();
    int * global_block_proc_dump = GetZoneProcessorIDDump();

    set< int > zoneIDinCurrentProc;
    int myid = GetCurrentProcessorID();
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        int multi_izone = zoneStart + iZone;

        global_block_proc     [multi_izone] = block_proc     [iZone];
        global_block_proc_grid[multi_izone] = block_proc_grid[iZone];
        global_block_type     [multi_izone] = block_type     [iZone];
        global_block_idx      [multi_izone] = block_idx      [iZone];
        global_block_file     [multi_izone] = block_fileProc [iZone];

        if (block_proc[iZone] == myid || block_proc_grid[iZone] == myid)
        {
            zoneIDinCurrentProc.insert(iZone);
        }
    }

    ostringstream oss;
    oss << "Zone index in current processor : ";
    for (set<int>::iterator iter = zoneIDinCurrentProc.begin(); iter != zoneIDinCurrentProc.end(); ++ iter)
    {
        int partID = * iter;
        oss << partID << " ";
    }
    oss << endl;
    WriteLogFile(oss.str());

    //! global_block_proc_dump would be deleted?
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        int multi_izone = zoneStart + iZone;
        global_block_proc_dump[multi_izone] = block_proc_dump[iZone];
    }
}

void GridGroup::InitZoneLayout(fstream &file)
{
    using namespace PHMPI;

    int file_proc;
    GetFileProc(file_proc);

    PH_Read_Bcast(file, &nzones, sizeof(int), file_proc);

    block_proc = new int [nzones];
    block_proc_grid = new int [nzones];
    block_fileProc = new int [nzones];

    //! Initialization.
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        block_proc[iZone] = 0;
        block_proc_grid[iZone] = 0;
        block_fileProc[iZone] = 0;
    }

    PH_Read_Bcast(file, block_proc, nzones*sizeof(int), file_proc);
    PH_Read_Bcast(file, block_proc_grid, nzones*sizeof(int), file_proc);

    //! block_proc_dump is deleted in ~GridGroup() function.
    block_proc_dump = new int [nzones];
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        block_proc_dump[iZone] = block_proc[iZone];
    }

    block_idx  = new int [nzones];
    block_type = new int [nzones];

    //! Read the index of each block,actually nb defaults to iZone.
    PH_Read_Bcast(file, block_idx, nzones*sizeof(int), file_proc);

    //! Read the grid type of each block,it can be unstructured and structured.
    PH_Read_Bcast(file, block_type, nzones*sizeof(int), file_proc);

    int m_block_proc = PHMPI::GetZoneDistributionMethod();

    int number_of_processor = GetNumberOfProcessor();

    if (m_block_proc == REDISTRIBUTION)
    {
        if (number_of_processor <= nzones)
        {
            for (int iZone = 0; iZone < nzones; ++ iZone)
            {
                block_proc[iZone] = (zoneStart + iZone) % number_of_processor;
            }
        }
        else
        {
            string cs;
            std::ostringstream oss;
            oss << "When the parallel policy is determined by the system (m_block_proc = 1), the number of processes started must be less than the number of partitions!\n" 
                << "   ---- Number of processes = " << number_of_processor << ", Number of partitions = " << nzones << " ----\n" 
                << "Abnormally Exit Program!!!\n";
            cs = oss.str();

            cout << cs << endl;

            TK_Exit::ExceptionExit(cs);
        }
    }
    else
    {
        if (number_of_processor == 1)
        {
            if (IsConvertGridToMixGrid())
            {
                //! For structured + unstructured = mix grid,move the specified process of  unstructured grid partition backward.
                if (block_type[0] == PHSPACE::STRUCTGRID)
                {
                    int maxproc = block_proc[0];
                    for (int iZone = 0; iZone < nzones; ++ iZone)
                    {
                        if (block_proc[iZone] > maxproc)
                        {
                            maxproc = block_proc[iZone];
                        }
                    }
                    maxproc += 1;
                    SetNumberOfProcStructUsed(maxproc);
                }
                else
                {
                    for (int iZone = 0; iZone < nzones; ++ iZone)
                    {
                        block_proc_dump[iZone] += GetNumberOfProcStructUsed();
                    }
                }
            }
            for (int iZone = 0; iZone < nzones; ++ iZone)
            {
                block_proc[iZone] = 0;
            }
        }
        else
        {
            int maxproc = block_proc[0];
            for (int iZone = 0; iZone < nzones; ++ iZone)
            {
                if (block_proc[iZone] > maxproc)
                {
                    maxproc = block_proc[iZone];
                }
            }

            maxproc += 1;
            if (number_of_processor != maxproc)
            {
                string cs;
                std::ostringstream oss;
                oss << "When a parallel policy is run with a specified partition (m_block_proc = 0), the number of processes started must be less than the number of partitions!\n" 
                    << "   ---- Number of processes = " << number_of_processor << ", Number of partitions = " << maxproc << " ----\n" 
                    << "Abnormally Exit Program!!!\n";
                cs = oss.str();

                cout << cs << endl;

                TK_Exit::ExceptionExit(cs);
            }
        }
    }
}

void GridGroup::AddGrid(Grid *grid)
{
    grids.push_back(grid);
}

bool GridGroup::IsZoneLayoutInitialized()
{
    if (this->block_idx || this->block_proc || this->block_type)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int GridGroup::GetNumberofGrid() const
{
    return static_cast<int>(grids.size());
}

}


