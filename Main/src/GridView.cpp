#include "GridView.h"
#include "Geo_SimpleBC.h"
#include "PHMpi.h"
#include "TK_Parse.h"
#include "Constants.h"
#include "GlobalDataBase.h"

using namespace std;

namespace PHSPACE
{
PHSimpleFace::PHSimpleFace()
{
    simpleFaceNext = 0;
}

PHSimpleFace::~PHSimpleFace()
{
    
}

void PHSimpleFace::SetKeyParameter(int *keyNumber)
{
    for (int i = 0; i < 7; ++ i)
    {
        this->keyNumber[i] = keyNumber[i];
    }

    return;
}

int PHSimpleFace::GetKeyParameter()
{
    return this->keyNumber[6];
}

int * PHSimpleFace::GetLocalParameter()
{
    return keyNumber;
}

void PHSimpleFace::CreateNewSimpleFace()
{
    this->simpleFaceNext = new PHSimpleFace();

    return;
}

PHSimpleFace * PHSimpleFace::GetNextPointer()
{
    return simpleFaceNext;
}

void PHSimpleFace::ShowSimpleFace()
{
    for (int i = 0; i < 7; ++ i)
    {
        cout << setw(8) << keyNumber[i];
    }
    cout << endl;

    return;
}


BoundaryCluster::BoundaryCluster()
{
    simpleFaceHead = 0;
    simpleFaceGuard = 0;
    slaveBlockList = 0;
}
BoundaryCluster::~BoundaryCluster()
{
    DeleteLinkedTable();

    if (slaveBlockList)
    {
        delete [] slaveBlockList;
    }
}

void BoundaryCluster::SetBlockDimension(int iDimension, int jDimension, int kDimension)
{
    this->iDimension = iDimension;
    this->jDimension = jDimension;
    this->kDimension = kDimension;

    return;
}

void BoundaryCluster::SetBlockLabel(int blockLabel)
{
    this->blockLabel = blockLabel;

    return;
}

void BoundaryCluster::SetProcessLabel(int processLabel)
{
    this->processLabel = processLabel;

    return;
}

void BoundaryCluster::SetOverallFaceNumber(int overallFaceNumber)
{
    this->overallFaceNumber = overallFaceNumber;

    return;
}

void BoundaryCluster::SetExternalFaceNumber(int externalFaceNumber)
{
    this->externalFaceNumber = externalFaceNumber;

    return;
}

int BoundaryCluster::GetExternalFaceNumber()
{
    return this->externalFaceNumber;
}

vector <vector <int> > BoundaryCluster::ExtractPhysicalBoundaryIndex()
{
    vector <vector<int> > pbi;
    vector <int> trans(6);

    simpleFaceGuard = simpleFaceHead;
    while(simpleFaceGuard)
    {
        const int *localParameter = simpleFaceGuard->GetLocalParameter();

        if (localParameter[6] == GRIDGEN_SPACE::INTERFACE)
        {
            simpleFaceGuard = simpleFaceGuard->GetNextPointer();
        }
        else
        {
            for (int i = 0; i < 6; ++ i)
            {
                trans[i] = localParameter[i] - 1;
            }

            pbi.push_back(trans);
        }

        simpleFaceGuard = simpleFaceGuard->GetNextPointer();
    }

    return pbi;
}

int BoundaryCluster::GetOverallFaceNumber()
{
    return this->overallFaceNumber;
}

int BoundaryCluster::GetProcessLabel()
{
    return this->processLabel;
}

void BoundaryCluster::SetKeyNumberPointer(int *keyNumber)
{
    this->keyNumber = keyNumber;

    return;
}

void BoundaryCluster::CreateNewSimpleFace(int iFace)
{
    if (iFace == 0)
    {
        simpleFaceHead = new PHSimpleFace();
        simpleFaceGuard = simpleFaceHead;
    }
    else
    {
        simpleFaceGuard->CreateNewSimpleFace();
        simpleFaceGuard = simpleFaceGuard->GetNextPointer();
    }
    simpleFaceGuard->SetKeyParameter(keyNumber);

    return;
}

void BoundaryCluster::DeleteLinkedTable()
{
    if (simpleFaceHead)
    {
        simpleFaceGuard = simpleFaceHead;
        simpleFaceHead = simpleFaceHead->GetNextPointer();
        delete simpleFaceGuard;
    }

    return;
}

void BoundaryCluster::ShowBoundaryCluster()
{
    simpleFaceGuard = simpleFaceHead;
    while(simpleFaceGuard)
    {
        simpleFaceGuard->ShowSimpleFace();
        simpleFaceGuard = simpleFaceGuard->GetNextPointer();
    }

    return;
}

PHSimpleFace * BoundaryCluster::GetSimpleFaceHead()
{
    return simpleFaceHead;
}

void BoundaryCluster::OutputHoleSurfaceGrid(ofstream &outFile, BasalGrid *basalGrid)
{
    int ni = basalGrid->GetIDimension();
    int nj = basalGrid->GetJDimension();

    RDouble *x = basalGrid->GetX();
    RDouble *y = basalGrid->GetY();
    RDouble *z = basalGrid->GetZ();

    int is, js, ks;
    for (int m = 0; m < numberOfSolidFaces; ++ m)
    {
        for (int k = kst[m]; k <= ked[m]; ++ k)
        {
            ks = k-1;
            for (int j = jst[m]; j <= jed[m]; ++ j)
            {
                js = j-1;
                for (int i = ist[m]; i <= ied[m]; ++ i)
                {
                    is = i-1;
                    int mm = ni*nj*ks + ni*js + is;
                    outFile.write((char *)&x[mm], sizeof(RDouble));
                }
            }
        }

        for (int k = kst[m]; k <= ked[m]; ++ k)
        {
            ks = k-1;
            for (int j = jst[m]; j <= jed[m]; ++ j)
            {
                js = j-1;
                for (int i = ist[m]; i <= ied[m]; ++ i)
                {
                    is = i-1;
                    int mm = ni*nj*ks + ni*js + is;
                    outFile.write((char *)&y[mm], sizeof(RDouble));
                }
            }
        }

        for (int k = kst[m]; k <= ked[m]; ++ k)
        {
            ks = k-1;
            for (int j = jst[m]; j <= jed[m]; ++ j)
            {
                js = j-1;
                for (int i = ist[m]; i <= ied[m]; ++ i)
                {
                    is = i-1;
                    int mm = ni*nj*ks + ni*js + is;
                    outFile.write((char *)&z[mm], sizeof(RDouble));
                }
            }
        }
    }

    return;
}

void BoundaryCluster::OutputHoleSurfaceDimension(ofstream &outFile)
{
    int ni, nj, nk;
    for (int iFace = 0; iFace < numberOfSolidFaces; ++ iFace)
    {
        ni = ied[iFace]-ist[iFace]+1;
        nj = jed[iFace]-jst[iFace]+1;
        nk = ked[iFace]-kst[iFace]+1;
        if (nk == 1)
        {
            outFile.write((char *)&ni, sizeof(int));
            outFile.write((char *)&nj, sizeof(int));
            outFile.write((char *)&nk, sizeof(int));
        }
        else if (nj == 1)
        {
            outFile.write((char *)&ni, sizeof(int));
            outFile.write((char *)&nk, sizeof(int));
            outFile.write((char *)&nj, sizeof(int));
        }
        else
        {
            outFile.write((char *)&nj, sizeof(int));
            outFile.write((char *)&nk, sizeof(int));
            outFile.write((char *)&ni, sizeof(int));
        }
    }

    return;
}

void BoundaryCluster::RecordSolidFaceInformation()
{
    simpleFaceGuard = simpleFaceHead;
    int counter = 0;

    while (simpleFaceGuard)
    {
        const int *localParameter = simpleFaceGuard->GetLocalParameter();

        if (localParameter[6] == -1)
        {
            simpleFaceGuard = simpleFaceGuard->GetNextPointer();
        }
        else if (localParameter[6] == GRIDGEN_SPACE::SOLID_SURFACE)
        {
            ist.push_back(localParameter[0]);
            ied.push_back(localParameter[1]);

            jst.push_back(localParameter[2]);
            jed.push_back(localParameter[3]);

            kst.push_back(localParameter[4]);
            ked.push_back(localParameter[5]);

            counter += 1;
        }

        simpleFaceGuard = simpleFaceGuard->GetNextPointer();
    }

    numberOfSolidFaces = counter;

    return;
}

void BoundaryCluster::SetSlaveBlockList()
{
    if (externalFaceNumber > 0)
    {
        int id = 0;
        slaveBlockList = new int[externalFaceNumber];
        simpleFaceGuard = simpleFaceHead;
        while (simpleFaceGuard)
        {
            int sign = simpleFaceGuard->GetKeyParameter();
            simpleFaceGuard = simpleFaceGuard->GetNextPointer();
            if (sign == -1)
            {
                int slaveBlockLabel = simpleFaceGuard->GetKeyParameter();
                slaveBlockList[id] = slaveBlockLabel;
                id += 1;
            }
        }
    }

    return;
}

int * BoundaryCluster::GetSlaveBlockList()
{
    return slaveBlockList;
}

LinkedBlocks::LinkedBlocks()
{
    linkedBlocksNext = 0;
    blockLabel       = 0;
}

LinkedBlocks::~LinkedBlocks()
{

}

void LinkedBlocks::CreateNewLinkedBlock()
{
    this->linkedBlocksNext = new LinkedBlocks();

    return;
}

void LinkedBlocks::SetBlockLabel(int blockLabel)
{
    this->blockLabel = blockLabel;

    return;
}

int LinkedBlocks::GetBlockLabel()
{
    return this->blockLabel;
}

LinkedBlocks * LinkedBlocks::GetNextPointer()
{
    return linkedBlocksNext;
}

BlockCluster::BlockCluster()
{
    groupLabel = 0;
    newBlockLabel = 0;
    linkedBlocksHead = 0;
    linkedBlocksGuard = 0;
    listColor = WHITE;
}

BlockCluster::~BlockCluster()
{
    DeleteLinkedTable();
}

void BlockCluster::CreatNewLinkedBlock()
{
    if (!linkedBlocksHead)
    {
        linkedBlocksHead = new LinkedBlocks();
        linkedBlocksGuard = linkedBlocksHead;
    }
    else
    {
        linkedBlocksGuard->CreateNewLinkedBlock();
        linkedBlocksGuard = linkedBlocksGuard->GetNextPointer();
    }

    return;
}

void BlockCluster::DeleteLinkedTable()
{
    if (linkedBlocksHead)
    {
        linkedBlocksGuard = linkedBlocksHead;
        linkedBlocksHead = linkedBlocksHead->GetNextPointer();
        delete linkedBlocksGuard;
    }

    return;
}

void BlockCluster::SetBlockLabel(int blockLabel)
{
    linkedBlocksGuard->SetBlockLabel(blockLabel);

    return;
}

void BlockCluster::ShowLinkedBlocks()
{
    linkedBlocksGuard = linkedBlocksHead;
    while (linkedBlocksGuard)
    {
        int blockLabel = linkedBlocksGuard->GetBlockLabel();
        linkedBlocksGuard = linkedBlocksGuard->GetNextPointer();
        cout << '\t' << blockLabel;
    }
    cout << endl;

    return;
}

void BlockCluster::SetListColor(int listColor)
{
    this->listColor = listColor;

    return;
}

void BlockCluster::SetGroupLabel(int groupLabel)
{
    this->groupLabel = groupLabel;

    return;
}

void BlockCluster::SetNewBlockLabel(int newBlockLabel)
{
    this->newBlockLabel = newBlockLabel;

    return;
}

int BlockCluster::GetNewBlockLabel()
{
    return this->newBlockLabel;
}

int BlockCluster::GetListColor()
{
    return this->listColor;
}

int BlockCluster::GetGroupLabel()
{
    return this->groupLabel;
}

LinkedBlocks * BlockCluster::GetLinkedBlocksHead()
{
    return linkedBlocksHead;
}

BasalGrid::BasalGrid()
{
    xCoordinate = 0;
    yCoordinate = 0;
    zCoordinate = 0;
}

BasalGrid::~BasalGrid()
{
    delete [] xCoordinate;
    delete [] yCoordinate;
    delete [] zCoordinate;
}

void BasalGrid::GenerateBasalSpace(int iDimension, int jDimension, int kDimension)
{
    this->iDimension = iDimension;
    this->jDimension = jDimension;
    this->kDimension = kDimension;
    pointNumber = iDimension * jDimension * kDimension;
    xCoordinate = new RDouble[pointNumber];
    yCoordinate = new RDouble[pointNumber];
    zCoordinate = new RDouble[pointNumber];

    return;
}

void BasalGrid::SetCoordinate(int aspect, int i, int j, int k, RDouble coordinate)
{
    int pointLabel = iDimension*jDimension*k+iDimension*j+i;
    if (aspect == ASPECTX)
    {
        xCoordinate[pointLabel] = coordinate;
    }
    else if (aspect == ASPECTY)
    {
        yCoordinate[pointLabel] = coordinate;
    }
    else
    {
        zCoordinate[pointLabel] = coordinate;
    }

    return;
}

int BasalGrid::GetIDimension()
{
    return this->iDimension;
}

int BasalGrid::GetJDimension()
{
    return this->jDimension;
}

int BasalGrid::GetKDimension()
{
    return this->kDimension;
}

RDouble BasalGrid::GetCoordinate(int aspect, int i, int j, int k)
{
    int pointLabel = iDimension*jDimension*k + iDimension*j + i;
    RDouble coordinate = 0;
    if (aspect == ASPECTX)
    {
        coordinate = xCoordinate[pointLabel];
    }
    else if (aspect == ASPECTY)
    {
        coordinate = yCoordinate[pointLabel];
    }
    else
    {
        coordinate = zCoordinate[pointLabel];
    }

    return coordinate;
}

BlockGroupManager::BlockGroupManager()
{
    basalGrid = 0;
    boundaryCluster = 0;
    blockCluster = 0;
    InverseMapping = 0;
    blockStart = 0;
    blockEnd = 0;
}

BlockGroupManager::~BlockGroupManager()
{
    delete [] basalGrid;
    delete [] boundaryCluster;
    delete [] blockCluster;
    delete [] InverseMapping;
    delete [] blockStart;
    delete [] blockEnd;
}

void BlockGroupManager::Run()
{
    InputOriginalBoundary();

    InputOriginalGrid();

    ShowOriginalBoundary();

    OutputStandardBoundary();

    GatherLinkedBlocks();

    ShowLinkedBlocks();

    MarkerBlocksRelation();

    PermuteBlockLabel();

    SeparateGrid();

    //OutputNewGrid();
    OutputNewSurfaceGrid();

    OutputNewBoundary();

    OutputZoneInverseMapping();

    return;
}

void BlockGroupManager::InputOriginalGrid()
{
    int iDimension, jDimension, kDimension;
    RDouble coordinate;
    string originalGridFileName;
    GlobalDataBase::GetData("originalGridFileName", &originalGridFileName, PHSTRING,1);

    ifstream infile(originalGridFileName.c_str(), ios_base::binary);
    infile.read((char *)&globalBlockNumber, sizeof(int));
    basalGrid = new BasalGrid[globalBlockNumber];
    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        infile.read((char *)&iDimension, sizeof(int));
        infile.read((char *)&jDimension, sizeof(int));
        infile.read((char *)&kDimension, sizeof(int));
        basalGrid[iBlock].GenerateBasalSpace(iDimension, jDimension, kDimension);
    }

    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        iDimension = basalGrid[iBlock].GetIDimension();
        jDimension = basalGrid[iBlock].GetJDimension();
        kDimension = basalGrid[iBlock].GetKDimension();
        for (int aspect = ASPECTX; aspect <= ASPECTZ; ++ aspect)
        {
            for (int k = 0; k < kDimension; ++ k)
            {
                for (int j = 0; j < jDimension; ++ j)
                {
                    for (int i = 0; i < iDimension; ++ i)
                    {
                        infile.read((char *)&coordinate, sizeof(RDouble));
                        basalGrid[iBlock].SetCoordinate(aspect, i, j, k, coordinate);
                    }
                }
            }
        }
    }

    infile.close();

    return;
}

void BlockGroupManager::InputOriginalBoundary()
{
    string line, word;
    string separator = " =\t\r\n#$,;";

    int iSolver, processLabel, iDimension, jDimension, kDimension;
    processLabel = 0;

    string temp, originalBoundaryFileName;
    fstream file;
    GlobalDataBase::GetData("originalBoundaryFileName", &originalBoundaryFileName, PHSTRING, 1);
    file.open(originalBoundaryFileName.c_str(), ios_base::in);
    file >> iSolver;
    file >> globalBlockNumber;
    boundaryCluster = new BoundaryCluster[globalBlockNumber];
    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        ReadNewLine(file, line);
        line = FindNextWord(line, word, separator);
        from_string< int >(iDimension, word, std::dec);
        line = FindNextWord(line, word, separator);
        from_string< int >(jDimension, word, std::dec);
        line = FindNextWord(line, word, separator);
        from_string< int >(kDimension, word, std::dec);

        line = FindNextWord(line, word, separator);
        if (word == "proc")
        {
            readProc = true;
            line = FindNextWord(line, word, separator);
            from_string< int >(processLabel, word, std::dec);
        }

        boundaryCluster[iBlock].SetBlockDimension(iDimension, jDimension, kDimension);
        boundaryCluster[iBlock].SetProcessLabel(processLabel);
        file >> temp;

        int blockLabel = iBlock + 1;
        boundaryCluster[iBlock].SetBlockLabel(blockLabel);
        int overallFaceNumber;
        file >> overallFaceNumber;
        boundaryCluster[iBlock].SetOverallFaceNumber(overallFaceNumber);
        int iFace = 0, externalFaceNumber = 0;
        for (int overallFace = 0; overallFace < overallFaceNumber; ++ overallFace)
        {
            file >> keyNumber[0] >> keyNumber[1] >> keyNumber[2] >> keyNumber[3] >> keyNumber[4] >> keyNumber[5] >> keyNumber[6];
            boundaryCluster[iBlock].SetKeyNumberPointer(keyNumber);
            boundaryCluster[iBlock].CreateNewSimpleFace(iFace);
            iFace += 1;
            if (keyNumber[6] == -1)
            {
                file >> keyNumber[0] >> keyNumber[1] >> keyNumber[2] >> keyNumber[3] >> keyNumber[4] >> keyNumber[5] >> keyNumber[6];
                boundaryCluster[iBlock].SetKeyNumberPointer(keyNumber);
                boundaryCluster[iBlock].CreateNewSimpleFace(iFace);
                iFace += 1;
                externalFaceNumber += 1;
            }
        }
        boundaryCluster[iBlock].SetExternalFaceNumber(externalFaceNumber);
        boundaryCluster[iBlock].SetSlaveBlockList();
    }

    file.close();

    return;
}


void BlockGroupManager::ShowOriginalBoundary()
{
    cout << "ShowOriginalData" << endl;
    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        cout << "iBlock = " << setw(8) << iBlock << endl;
        boundaryCluster[iBlock].ShowBoundaryCluster();
    }

    return;
}

void BlockGroupManager::GatherLinkedBlocks()
{
    blockCluster = new BlockCluster[globalBlockNumber];
    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        int externalFaceNumber = boundaryCluster[iBlock].GetExternalFaceNumber();
        if (externalFaceNumber > 0)
        {
            const int *slaveBlockList = boundaryCluster[iBlock].GetSlaveBlockList();
            for (int iFace = 0; iFace < externalFaceNumber; ++ iFace)
            {
                blockCluster[iBlock].CreatNewLinkedBlock();
                blockCluster[iBlock].SetBlockLabel(slaveBlockList[iFace]);
            }
        }
    }

    return;
}

void BlockGroupManager::ShowLinkedBlocks()
{
    cout << "ShowLinkedBlocks" << endl;
    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        cout << "iBlock = " << iBlock << endl;
        blockCluster[iBlock].ShowLinkedBlocks();
    }

    return;
}

void BlockGroupManager::MarkerBlocksRelation()
{
    int groupCounter = 0;
    while (true)
    {
        int whiteLabel = SearchColorfulList(WHITE);
        if (whiteLabel < 0)
        {
            break;
        }
        else
        {
            if (!blockCluster[whiteLabel].GetLinkedBlocksHead())
            {
                blockCluster[whiteLabel].SetListColor(RED);
                blockCluster[whiteLabel].SetGroupLabel(groupCounter);
                groupCounter += 1;
            }
            else
            {
                blockCluster[whiteLabel].SetListColor(BLUE);
                while (true)
                {
                    int blueLabel = SearchColorfulList(BLUE);
                    LinkedBlocks *listGuide = blockCluster[blueLabel].GetLinkedBlocksHead();
                    if (blueLabel >= 0)
                    {
                        while (listGuide)
                        {
                            int iSlaveBlock = listGuide->GetBlockLabel();
                            blockCluster[iSlaveBlock-1].SetListColor(BLUE);
                            listGuide = listGuide->GetNextPointer();
                        }
                        blockCluster[blueLabel].DeleteLinkedTable();
                        blockCluster[blueLabel].SetListColor(PURPLE);
                    }
                    else
                    {
                        break;
                    }
                }

                for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
                {
                    int color = blockCluster[iBlock].GetListColor();
                    if (color == PURPLE)
                    {
                        blockCluster[iBlock].SetListColor(RED);
                        blockCluster[iBlock].SetGroupLabel(groupCounter);
                    }
                }
                groupCounter += 1;
            }
        }
    }

    this->groupNumber = groupCounter;

    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        int i = blockCluster[iBlock].GetGroupLabel();
        cout << "iBlock = " << iBlock+1 << '\t' << "group = " << i << endl;
    }

    return;
}

int BlockGroupManager::SearchColorfulList(int spectrum)
{
    int blockLabel = -1;
    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        int color = blockCluster[iBlock].GetListColor();
        if (color == spectrum)
        {
            blockLabel = iBlock;
            break;
        }
    }

    return blockLabel;
}

void BlockGroupManager::PermuteBlockLabel()
{
    int newBlockLabel = 0;
    for (int iGroup = 0; iGroup < groupNumber; ++ iGroup)
    {
        for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
        {
            int jGroup = blockCluster[iBlock].GetGroupLabel();
            if (iGroup == jGroup)
            {
                blockCluster[iBlock].SetNewBlockLabel(newBlockLabel);
                newBlockLabel += 1;
            }
        }
    }

    InverseMapping = new int [globalBlockNumber];

    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        newBlockLabel = blockCluster[iBlock].GetNewBlockLabel();
        InverseMapping[newBlockLabel] = iBlock;
    }

    blockStart = new int [groupNumber];
    blockStart[0] = 0;
    blockEnd = new int [groupNumber];
    blockEnd[groupNumber-1] = globalBlockNumber - 1;

    int groupLabel = 0;
    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        int oldBlockLabel = InverseMapping[iBlock];
        int groupGuard = blockCluster[oldBlockLabel].GetGroupLabel();
        if (groupGuard != groupLabel)
        {
            blockEnd[groupLabel] = iBlock - 1;
            groupLabel += 1;
            blockStart[groupLabel] = iBlock;
        }
    }

    return;
}

void BlockGroupManager::SeparateGrid()
{
    string gridGroupFileName;

    GlobalDataBase::GetData("gridGroupFileName", &gridGroupFileName, PHSTRING, 1);

    ofstream myFile(gridGroupFileName.c_str());
    
    for (int iGroup = 0; iGroup < groupNumber; ++ iGroup)
    {
        OutputNewBlockGroup(iGroup, myFile);

        OutputHoleGroup(iGroup);
    }

    myFile.close();

    return;
}

void BlockGroupManager::OutputNewBlockGroup(int iGroup, ofstream &myFile)
{
    ostringstream oss;
    string gridGroupHeadName;
    GlobalDataBase::GetData("gridGroupHeadName", &gridGroupHeadName, PHSTRING, 1);
    oss << gridGroupHeadName << iGroup << ".grd";
    string file = oss.str();
    ofstream outFile(file.c_str(), ios_base::binary);

    int localBlockNumber = blockEnd[iGroup] - blockStart[iGroup] + 1;
    myFile << "iGroup = " << setw(8) << iGroup << "  blockStart = " << setw(8) << blockStart[iGroup] << "  blockEnd = " << setw(8) << blockEnd[iGroup] << endl;

    vector< vector < vector < int > > > bcIndexTable;
    bcIndexTable.resize(localBlockNumber);

    int iCounter = 0;
    for (int iBlock = blockStart[iGroup]; iBlock <= blockEnd[iGroup]; ++ iBlock)
    {
        int oldBlockLabel = InverseMapping[iBlock];
        BoundaryCluster &bc = boundaryCluster[oldBlockLabel];
        bcIndexTable[iCounter] = bc.ExtractPhysicalBoundaryIndex();

        iCounter += 1;
    }

    uint_t numberOfPhysicsFaces = 0;
    for (int m = 0; m < localBlockNumber; ++ m)
    {
        numberOfPhysicsFaces += bcIndexTable[m].size();
    }
    outFile.write((char *)&numberOfPhysicsFaces, sizeof(int));

    int si, sj, sk;
    for (int m = 0; m < localBlockNumber; ++ m)
    {
        if (bcIndexTable[m].size() == 0) continue;

        for (std::size_t n = 0; n < bcIndexTable[m].size(); ++ n)
        {
            si = bcIndexTable[m][n][1] - bcIndexTable[m][n][0] + 1;
            sj = bcIndexTable[m][n][3] - bcIndexTable[m][n][2] + 1;
            sk = bcIndexTable[m][n][5] - bcIndexTable[m][n][4] + 1;
            if (sk == 1)
            {
                outFile.write((char *)&si,  sizeof(int));
                outFile.write((char *)&sj,  sizeof(int));
                outFile.write((char *)&sk,  sizeof(int));
            }
            else if (sj == 1)
            {
                outFile.write((char *)&si,  sizeof(int));
                outFile.write((char *)&sk,  sizeof(int));
                outFile.write((char *)&sj,  sizeof(int));
            }
            else
            {
                outFile.write((char *)&sj,  sizeof(int));
                outFile.write((char *)&sk,  sizeof(int));
                outFile.write((char *)&si,  sizeof(int));
            }
        }
    }

    for (int m = 0; m < localBlockNumber; ++ m)
    {
        if (bcIndexTable[m].size() == 0) continue;

        int oldBlockLabel = InverseMapping[blockStart[iGroup] + m];
        for (std::size_t n = 0; n < bcIndexTable[m].size(); ++ n)
        {
            int ist = bcIndexTable[m][n][0];
            int ied = bcIndexTable[m][n][1];

            int jst = bcIndexTable[m][n][2];
            int jed = bcIndexTable[m][n][3];

            int kst = bcIndexTable[m][n][4];
            int ked = bcIndexTable[m][n][5];

            for (int aspect = ASPECTX; aspect <= ASPECTZ; ++ aspect)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            RDouble coordinate = basalGrid[oldBlockLabel].GetCoordinate(aspect, i, j, k);
                            outFile.write((char *)&coordinate, sizeof(RDouble));
                        }
                    }
                }
            }
        }
    }

    outFile.close();
    return;
}

void BlockGroupManager::OutputNewSurfaceGrid()
{
    vector <int> localNumberOfSurfaces, localNumberOfPoints;
    vector <int> ni, nj;
    vector <RDouble> coordinate;

    for (int iGroup = 0; iGroup < groupNumber; ++ iGroup)
    {
        ostringstream oss;
        string gridGroupHeadName;
        GlobalDataBase::GetData("gridGroupHeadName", &gridGroupHeadName, PHSTRING, 1);
        oss << gridGroupHeadName << iGroup << ".grd";
        string file = oss.str();
        ifstream infile(file.c_str(), ios_base::binary);

        int ns;
        infile.read((char *)&ns, sizeof(int));
        localNumberOfSurfaces.push_back(ns);

        int is, js, ks, np = 0;
        for (int id = 0; id < ns; ++ id)
        {
            infile.read((char *)&is, sizeof(int));
            ni.push_back(is);

            infile.read((char *)&js, sizeof(int));
            nj.push_back(js);

            infile.read((char *)&ks, sizeof(int));
            np += is * js;
        }

        localNumberOfPoints.push_back(np);

        RDouble ss;
        for (int ip = 0; ip < 3 * np; ++ ip)
        {
            infile.read((char *)&ss, sizeof(RDouble));
            coordinate.push_back(ss);
        }
        infile.close();
    }

    string newGridFileName;
    GlobalDataBase::GetData("newGridFileName", &newGridFileName, PHSTRING,1);
    ofstream outFile(newGridFileName.c_str(), ios_base::binary);

    int nst = 0, npt = 0;
    for (int iGroup = 0; iGroup < groupNumber; ++ iGroup)
    {
        nst += localNumberOfSurfaces[iGroup];
        npt += localNumberOfPoints[iGroup];
    }
    outFile.write((char *)&nst, sizeof(int));

    int kk = 1;
    for (int iSurface = 0; iSurface < nst; ++ iSurface)
    {
        outFile.write((char *)&ni[iSurface], sizeof(int));
        outFile.write((char *)&nj[iSurface], sizeof(int));
        outFile.write((char *)&kk          , sizeof(int));
    }

    outFile.write((char *)&coordinate[0], 3 * npt * sizeof(RDouble));
    outFile.close();

    string surfaceGroupFileName;
    GlobalDataBase::GetData("surfaceGroupFileName", &surfaceGroupFileName, PHSTRING,1);
    ofstream surfile(surfaceGroupFileName.c_str());
    int istart = 0;
    for (int iGroup = 0; iGroup < groupNumber; ++ iGroup)
    {
        int iSurfaceStart = istart, iSurfaceEnd = istart + localNumberOfSurfaces[iGroup] - 1;
        surfile << "iGroup = " << setw(8) << iGroup << "  surfaceStart = " << setw(8) << iSurfaceStart  << "  surfaceEnd = " << setw(8) << iSurfaceEnd << endl;
        istart = iSurfaceEnd + 1;
    }
    surfile.close();

    return;
}

void BlockGroupManager::OutputNewGrid()
{
    int iDimension, jDimension, kDimension;
    string newGridFileName;
    GlobalDataBase::GetData("newGridFileName", &newGridFileName, PHSTRING,1);
    ofstream outFile(newGridFileName.c_str(), ios::binary);
    outFile.write((char *)&globalBlockNumber, sizeof(int));

    for (int iGroup = 0; iGroup < groupNumber; ++ iGroup)
    {
        for (int iBlock = blockStart[iGroup]; iBlock <= blockEnd[iGroup]; ++ iBlock)
        {
            int oldBlockLabel = InverseMapping[iBlock];
            iDimension = basalGrid[oldBlockLabel].GetIDimension();
            jDimension = basalGrid[oldBlockLabel].GetJDimension();
            kDimension = basalGrid[oldBlockLabel].GetKDimension();

            outFile.write((char *)&iDimension, sizeof(int));
            outFile.write((char *)&jDimension, sizeof(int));
            outFile.write((char *)&kDimension, sizeof(int));
        }
    }

    for (int iGroup = 0; iGroup < groupNumber; ++ iGroup)
    {
        for (int iBlock = blockStart[iGroup]; iBlock <= blockEnd[iGroup]; ++ iBlock)
        {
            int oldBlockLabel = InverseMapping[iBlock];
            iDimension = basalGrid[oldBlockLabel].GetIDimension();
            jDimension = basalGrid[oldBlockLabel].GetJDimension();
            kDimension = basalGrid[oldBlockLabel].GetKDimension();

            for (int aspect = ASPECTX; aspect <= ASPECTZ; ++ aspect)
            {
                for (int k = 0; k < kDimension; ++ k)
                {
                    for (int j = 0; j < jDimension; ++ j)
                    {
                        for (int i = 0; i < iDimension; ++ i)
                        {
                            RDouble coordinate = basalGrid[oldBlockLabel].GetCoordinate(aspect, i, j, k);
                            outFile.write((char *)&coordinate, sizeof(RDouble));
                        }
                    }
                }
            }
        }
    }
    outFile.close();

    return;
}

void BlockGroupManager::OutputNewBoundary()
{
    int iDimension, jDimension, kDimension, lastPara;
    int iSolver = 1;
    PHSimpleFace *simpleFaceGuard = 0;
    string newBoundaryFileName;
    GlobalDataBase::GetData("newBoundaryFileName", &newBoundaryFileName, PHSTRING,1);
    ofstream outFile(newBoundaryFileName.c_str());
    outFile << iSolver << endl;
    outFile << globalBlockNumber << endl;
    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        int oldBlockLabel = InverseMapping[iBlock];
        iDimension = basalGrid[oldBlockLabel].GetIDimension();
        jDimension = basalGrid[oldBlockLabel].GetJDimension();
        kDimension = basalGrid[oldBlockLabel].GetKDimension();

        outFile << iDimension << setw(4) << jDimension << setw(4) << kDimension << setw(8) << "proc " << iBlock << endl;
        
        outFile << "block" << iBlock+1 << endl;
        int overallFaceNumber = boundaryCluster[oldBlockLabel].GetOverallFaceNumber();
        outFile << setw(8) << overallFaceNumber << endl;
        int externalSign = 0;
        simpleFaceGuard = boundaryCluster[oldBlockLabel].GetSimpleFaceHead();

        while (simpleFaceGuard)
        {
            const int *localParameter = simpleFaceGuard->GetLocalParameter();
            if (externalSign == 0)
            {
                lastPara = localParameter[6];
                if (lastPara == -1)
                {
                    externalSign = 1;
                }
            }
            else
            {
                int oldSlaveLabel = localParameter[6]-1;
                lastPara = blockCluster[oldSlaveLabel].GetNewBlockLabel()+1;
                externalSign = 0;
            }
            
            for (int i = 0; i < 6; ++ i)
            {
                outFile << setw(8) << localParameter[i];
            }
            outFile << setw(8) << lastPara << endl;
            simpleFaceGuard = simpleFaceGuard->GetNextPointer();
        }
    }
    outFile.close();

    return;
}

void BlockGroupManager::OutputStandardBoundary()
{
    using namespace PHMPI;

    GlobalDataBase::GetData("originalExternalLabel", &originalExternalLabel, PHINT, 1);

    int iDimension, jDimension, kDimension, lastPara;
    int iSolver = 1;
    PHSimpleFace *simpleFaceGuard = 0;

    string standardBoundaryFileName;
    GlobalDataBase::GetData("standardBoundaryFileName", &standardBoundaryFileName, PHSTRING,1);
    ofstream outFile(standardBoundaryFileName.c_str());
    outFile << iSolver << endl;
    outFile << globalBlockNumber << endl;
    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        int oldBlockLabel = iBlock;
        iDimension = basalGrid[oldBlockLabel].GetIDimension();
        jDimension = basalGrid[oldBlockLabel].GetJDimension();
        kDimension = basalGrid[oldBlockLabel].GetKDimension();

        int proc = boundaryCluster[iBlock].GetProcessLabel();

        if (!readProc)
        {
            outFile << iDimension << setw(4) << jDimension << setw(4) << kDimension << endl;
        }
        else
        {
            outFile << iDimension << setw(4) << jDimension << setw(4) << kDimension << setw(8) << "proc " << proc << endl;
        }

        outFile << "block" << iBlock+1 << endl;
        int internalFaceNumber = boundaryCluster[oldBlockLabel].GetOverallFaceNumber();
        outFile << setw(8) << internalFaceNumber << endl;
        int externalSign = 0;
        simpleFaceGuard = boundaryCluster[oldBlockLabel].GetSimpleFaceHead();

        while (simpleFaceGuard)
        {
            const int *localParameter = simpleFaceGuard->GetLocalParameter();
            if (externalSign == 0)
            {
                lastPara = localParameter[6];

                if (lastPara == originalExternalLabel)
                {
                    lastPara = EXTERNAL_BC;
                }

                if (lastPara == -1)
                {
                    externalSign = 1;
                }
            }
            else
            {
                int oldSlaveLabel = localParameter[6]-1;
                lastPara = oldSlaveLabel+1;
                externalSign = 0;
            }
            
            for (int i = 0; i < 6; ++ i)
            {
                outFile << setw(8) << localParameter[i];
            }
            outFile << setw(8) << lastPara << endl;
            simpleFaceGuard = simpleFaceGuard->GetNextPointer();
        }
    }
    outFile.close();

    return;
}

void BlockGroupManager::OutputZoneInverseMapping()
{
    string inverseMappingFileName;
    GlobalDataBase::GetData("inverseMappingFileName", &inverseMappingFileName, PHSTRING,1);
    ofstream outFile(inverseMappingFileName.c_str(), ios::binary);
    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        int oldBlockLabel = InverseMapping[iBlock];
        cout << "iBlock = " << iBlock << "  oldBlockLabel = " << oldBlockLabel << endl;
        outFile.write((char *)&oldBlockLabel, sizeof(int));
    }
    outFile.close();

    return;
}

void BlockGroupManager::OutputHoleGroup(int iGroup)
{
    ostringstream oss;
    string holeGroupHeadName;
    GlobalDataBase::GetData("holeGroupHeadName", &holeGroupHeadName, PHSTRING, 1);
    oss << holeGroupHeadName << iGroup << ".grd";
    string name = oss.str();
    ofstream outFile(name.c_str(), ios_base::binary);

    int solidFaceCounter = 0;

    for (int iBlock = blockStart[iGroup]; iBlock <= blockEnd[iGroup]; ++ iBlock)
    {
        int oldBlockLabel = InverseMapping[iBlock];
        boundaryCluster[oldBlockLabel].RecordSolidFaceInformation();
        solidFaceCounter += boundaryCluster[oldBlockLabel].GetNumberOfSolidFaces();
    }

    outFile.write((char *)&solidFaceCounter, sizeof(int));

    for (int iBlock = blockStart[iGroup]; iBlock <= blockEnd[iGroup]; ++ iBlock)
    {
        int oldBlockLabel = InverseMapping[iBlock];
        boundaryCluster[oldBlockLabel].OutputHoleSurfaceDimension(outFile);
    }

    for (int iBlock = blockStart[iGroup]; iBlock <= blockEnd[iGroup]; ++ iBlock)
    {
        int oldBlockLabel = InverseMapping[iBlock];
        boundaryCluster[oldBlockLabel].OutputHoleSurfaceGrid(outFile, &basalGrid[oldBlockLabel]);
    }

    outFile.close();

    return;
}

}