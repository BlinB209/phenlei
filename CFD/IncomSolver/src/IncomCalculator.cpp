#include "IncomCalculator.h"




namespace PHSPACE
{

    int IncomCalculator::GetSolverIndex(int iEquation)
    {
        int equationIndex;
        if (iEquation < IDX::S_IP)
        {
            equationIndex = iEquation;
        }
        else
        {
            equationIndex = GetSolverIndex();
        }
        return equationIndex;
    }

    void IncomCalculator::calculateLinearEquation(Grid *grid, int iEquation, RDouble *x, RDouble *diagCoef, RDouble *b, RDouble *upACOFFace, RDouble *downACOFFace)
    {
        int *linearSolverID = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("linearSolverID"));
        int *constVarNumToSolveVarNum = reinterpret_cast<int *>(grid->GetDataPtr("constVarNumToSolveVarNum"));

        AxEqualbCalculator->mathLibSolve(grid, iEquation, x, diagCoef, b, upACOFFace, downACOFFace);
        CalcRes(grid, iEquation, x, diagCoef, b, upACOFFace, downACOFFace);
    }

    void IncomCalculator::SetAxEqualbCalculator(LinearEquationSystemCalculator *AxEqualbCalculatorIn)
    {
        this->AxEqualbCalculator = AxEqualbCalculatorIn;
    }

    void IncomCalculator::CalcRes(Grid *gridIn, int iEquation,
        RDouble *x, RDouble *diagCoef, RDouble *b, RDouble *upACOFFace, RDouble *downACOFFace)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nTotalFace = grid->GetNTotalFace();
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();
        RDouble *vol = grid->GetCellVolume();

        RDouble *residualError = new RDouble[nTotalCell];
        string *varNameIncom = reinterpret_cast <string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
        string varName = varNameIncom[iEquation];

        for (int iCell = 0; iCell < nTotalCell; iCell++ )
        {
            residualError[iCell] = b[iCell] - diagCoef[iCell] * x[iCell];
        }

        for (int iFace = nBoundFace; iFace < nTotalFace; iFace++)
        {
            int re = rightCellOfFace[iFace];
            int le = leftCellOfFace[iFace];

            residualError[le] -= upACOFFace[iFace] * x[re];
            residualError[re] -= downACOFFace[iFace] * x[le];
        }

        UnstructBCSet **bcr = grid->GetBCRecord();
        for (int iFace = 0; iFace < nBoundFace; iFace++)
        {
            int bcType = bcr[iFace]->GetKey();
            if (bcType == PHENGLEI::INTERFACE)
            {
                int le = leftCellOfFace[iFace];
                int re = rightCellOfFace[iFace];
                residualError[le] -= upACOFFace[iFace] * x[re];
            }
        }

        if (iEquation == IDX::S_IP)
        {
            RDouble globalTotalCells = GlobalDataBase::GetDoubleParaFromDB("GlobalTotalCells");
            for (int iCell = 0; iCell < nTotalCell; iCell++ )
            {
                residualError[iCell] = x[iCell] / globalTotalCells;
            }
        }

        RDouble varTotalRes = 0.0;        
        RDouble globalTotalCells = 1.0;
        int ResType = L1Res;
        if (GlobalDataBase::IsExist("ResType", PHINT, 1))
        {
            ResType = GlobalDataBase::GetIntParaFromDB("ResType");
        }

        if (ResType == L1Res)
        {
            for (int iCell = 0; iCell < nTotalCell; ++iCell)
            {
                varTotalRes += ABS(residualError[iCell]) / globalTotalCells;
            }
            RDouble temp1 = varTotalRes ;
            PH_AllreduceSepMode(&temp1, &varTotalRes, 1, MPI_DOUBLE, MPI_SUM);
        }
        else if (ResType == L2Res)
        {
            for (int iCell = 0; iCell < nTotalCell; ++iCell)
            {
                varTotalRes += pow(residualError[iCell], 2) / globalTotalCells;
            }
            RDouble temp1 = varTotalRes;
            PH_AllreduceSepMode(&temp1, &varTotalRes, 1, MPI_DOUBLE, MPI_SUM);
        }
        else if (ResType == LinfRes)
        {
            for (int iCell = 0; iCell < nTotalCell; ++iCell)
            {
                if (varTotalRes < ABS(residualError[iCell]))
                {
                    varTotalRes = ABS(residualError[iCell]);
                }
            }
            RDouble temp1 = varTotalRes;
            PH_AllreduceSepMode(&temp1, &varTotalRes, 1, MPI_DOUBLE, MPI_MAX);
        }

        grid->UpdateData(varName + "ResNow", &varTotalRes, PHDOUBLE, 1);

        delete [] residualError;
    }

    void IncomCalculator::SkewnessCorrection(Grid *grid, RDouble *phi, RDouble *dphidx, RDouble *dphidy, RDouble *dphidz)
    {
        int nInteriorFaces = grid->GetNTotalFace() - grid->GetNBoundFace() + grid->GetNIFace();
        RDouble *dphidxf = new RDouble[nInteriorFaces];
        RDouble *dphidyf = new RDouble[nInteriorFaces];
        RDouble *dphidzf = new RDouble[nInteriorFaces];
        RDouble *phif = new RDouble[nInteriorFaces];

        for (int iter = 0; iter < 3; ++iter)
        {
            FaceGradCorrection(grid, dphidxf, dphidyf, dphidzf, phi, dphidx, dphidy, dphidz);
            FaceValueCorretion(grid, dphidxf, dphidyf, dphidzf, phif, phi);
            GradCorrection(grid, phif, phi, dphidx, dphidy, dphidz);
        }

        delete [] dphidxf;
        delete [] dphidyf;
        delete [] dphidzf;
        delete [] phif;
    }

    void IncomCalculator::FaceGradCorrection(Grid *gridIn, RDouble *dphidxf, RDouble *dphidyf, RDouble *dphidzf, RDouble *phi, RDouble *dphidx, RDouble *dphidy, RDouble *dphidz)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nbfaces = grid->GetNBoundFace();
        int nIFaces = grid->GetNIFace();
        int nInteriorFaces = grid->GetNTotalFace() - nbfaces + nIFaces;
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();
        RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));

        RDouble *xcc = grid->GetCellCenterX();
        RDouble *ycc = grid->GetCellCenterY();
        RDouble *zcc = grid->GetCellCenterZ();

        for (int iface = nIFaces; iface < nInteriorFaces; ++iface)
        {
            int faceIndex = iface + nbfaces - nIFaces;
            int le = leftCellOfFace[faceIndex];
            int re = rightCellOfFace[faceIndex];

            RDouble adphidxf = faceWeightOfLeftCell[faceIndex] * dphidx[le] + (1.0 - faceWeightOfLeftCell[faceIndex]) * dphidx[re];
            RDouble adphidyf = faceWeightOfLeftCell[faceIndex] * dphidy[le] + (1.0 - faceWeightOfLeftCell[faceIndex]) * dphidy[re];
            RDouble adphidzf = faceWeightOfLeftCell[faceIndex] * dphidz[le] + (1.0 - faceWeightOfLeftCell[faceIndex]) * dphidz[re];

            RDouble dCFx = xcc[re] - xcc[le];
            RDouble dCFy = ycc[re] - ycc[le];
            RDouble dCFz = zcc[re] - zcc[le];
            RDouble dCFNorm = sqrt(dCFx * dCFx + dCFy * dCFy + dCFz * dCFz);
            RDouble correction1 = (phi[re] - phi[le]) / dCFNorm;

            RDouble eCFx = dCFx / dCFNorm;
            RDouble eCFy = dCFy / dCFNorm;
            RDouble eCFz = dCFz / dCFNorm;
            RDouble a2 = correction1 - adphidxf * eCFx + adphidyf * eCFy + adphidzf * eCFz;
            RDouble correction2x = a2 * eCFx;
            RDouble correction2y = a2 * eCFy;
            RDouble correction2z = a2 * eCFz;

            dphidxf[iface] = adphidxf + correction2x;
            dphidyf[iface] = adphidyf + correction2y;
            dphidzf[iface] = adphidzf + correction2z;
        }

        for (int iface = 0; iface < nIFaces; ++iface)
        {
            int faceIndex = iface + nbfaces - nIFaces;
            int le = leftCellOfFace[faceIndex];
            int re = rightCellOfFace[faceIndex];

            RDouble adphidxf = dphidx[le];
            RDouble adphidyf = dphidy[le];
            RDouble adphidzf = dphidz[le];

            RDouble dCFx = xcc[re] - xcc[le];
            RDouble dCFy = ycc[re] - ycc[le];
            RDouble dCFz = zcc[re] - zcc[le];
            RDouble dCFNorm = sqrt(dCFx * dCFx + dCFy * dCFy + dCFz * dCFz);
            RDouble correction1 = (phi[re] - phi[le]) / dCFNorm;

            RDouble eCFx = dCFx / dCFNorm;
            RDouble eCFy = dCFy / dCFNorm;
            RDouble eCFz = dCFz / dCFNorm;
            RDouble a2 = correction1 - adphidxf * eCFx + adphidyf * eCFy + adphidzf * eCFz;
            RDouble correction2x = a2 * eCFx;
            RDouble correction2y = a2 * eCFy;
            RDouble correction2z = a2 * eCFz;

            dphidxf[iface] = adphidxf + correction2x;
            dphidyf[iface] = adphidyf + correction2y;
            dphidzf[iface] = adphidzf + correction2z;
        }
    }

    void IncomCalculator::FaceValueCorretion(Grid *gridIn, RDouble *dphidxf, RDouble *dphidyf, RDouble *dphidzf, RDouble *phif, RDouble *phi)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nbfaces = grid->GetNBoundFace();
        int nIFaces = grid->GetNIFace();
        int nInteriorFaces = grid->GetNTotalFace() - nbfaces + nIFaces;
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();
        RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));
        RDouble *area = grid->GetFaceArea();
        RDouble *xfn = grid->GetFaceNormalX();
        RDouble *yfn = grid->GetFaceNormalY();
        RDouble *zfn = grid->GetFaceNormalZ();
        RDouble *xfc = grid->GetFaceCenterX();
        RDouble *yfc = grid->GetFaceCenterY();
        RDouble *zfc = grid->GetFaceCenterZ();
        RDouble *xcc = grid->GetCellCenterX();
        RDouble *ycc = grid->GetCellCenterY();
        RDouble *zcc = grid->GetCellCenterZ();

        for (int iface = 0; iface < nInteriorFaces; ++iface)
        {
            int ifaceIndex = iface + nbfaces - nIFaces;
            int le = leftCellOfFace[ifaceIndex];
            int re = rightCellOfFace[ifaceIndex];

            RDouble dCfx = xfc[ifaceIndex] - xcc[le];
            RDouble dCfy = yfc[ifaceIndex] - ycc[le];
            RDouble dCfz = zfc[ifaceIndex] - zcc[le];
            RDouble dCfNorm = sqrt(dCfx * dCfx + dCfy * dCfy + dCfz * dCfz);

            RDouble dCFx = xcc[re] - xcc[le];
            RDouble dCFy = ycc[re] - ycc[le];
            RDouble dCFz = zcc[re] - zcc[le];
            RDouble dCFNorm = sqrt(dCFx * dCFx + dCFy * dCFy + dCFz * dCFz);

            RDouble ax = area[ifaceIndex] * xfn[ifaceIndex];
            RDouble ay = area[ifaceIndex] * yfn[ifaceIndex];
            RDouble az = area[ifaceIndex] * zfn[ifaceIndex];

            RDouble proj1 = dCfx * ax + dCfy * ay + dCfz * az;
            RDouble proj2 = dCFx * ax + dCFy * ay + dCFz * az;
            RDouble ratio = proj1 / proj2;

            // f': CF intersects the surface
            RDouble fQuotex = xcc[le] + dCFx * ratio;
            RDouble fQuotey = ycc[le] + dCFy * ratio;
            RDouble fQuotez = zcc[le] + dCFz * ratio;

            // ffQuotex: vec f'f
            RDouble ffQuotex = xfc[ifaceIndex] - fQuotex;
            RDouble ffQuotey = yfc[ifaceIndex] - fQuotey;
            RDouble ffQuotez = zfc[ifaceIndex] - fQuotez;

            RDouble correction = dphidxf[iface] * ffQuotex + dphidyf[iface] * ffQuotey + dphidzf[iface] * ffQuotez;
            RDouble phifQuote = faceWeightOfLeftCell[ifaceIndex] * phi[le] + (1.0 - faceWeightOfLeftCell[ifaceIndex]) * phi[re];
            phif[iface] = phifQuote + correction;
        }
    }

    void IncomCalculator::GradCorrection(Grid *gridIn, RDouble *phif, RDouble *phi, RDouble *dphidx, RDouble *dphidy, RDouble *dphidz)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);

        RDouble *xfn = grid->GetFaceNormalX();
        RDouble *yfn = grid->GetFaceNormalY();
        RDouble *zfn = grid->GetFaceNormalZ();
        RDouble *area = grid->GetFaceArea();
        RDouble *vol = grid->GetCellVolume();

        int *leftCellofFace = grid->GetLeftCellOfFace();
        int *rightCellofFace = grid->GetRightCellOfFace();

        int nTotalCell = grid->GetNTotalCell();
        int nTotalFace = grid->GetNTotalFace();
        int nBoundFace = grid->GetNBoundFace();
        int nInterFace = grid->GetNIFace();
        int nBoundFaces = nBoundFace - nInterFace;          //! except interface
        int nTotal = nTotalCell + nBoundFace;

        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            dphidx[iCell] = 0.0;
            dphidy[iCell] = 0.0;
            dphidz[iCell] = 0.0;
        }

        for (int iFace = 0; iFace < nBoundFaces; ++iFace)
        {
            int le = leftCellofFace[iFace];
            int re = iFace + nTotalCell;

            RDouble qfc = phi[re];
            RDouble nx = xfn[iFace] * area[iFace];
            RDouble ny = yfn[iFace] * area[iFace];
            RDouble nz = zfn[iFace] * area[iFace];
            dphidx[le] += qfc * nx;
            dphidy[le] += qfc * ny;
            dphidz[le] += qfc * nz;
        }

        for (int iFace = nBoundFaces; iFace < nBoundFace; ++iFace)
        {
            int le = leftCellofFace[iFace];
            int re = rightCellofFace[iFace];
            int faceIndex = iFace - nBoundFaces;

            RDouble nx = xfn[iFace] * area[iFace];
            RDouble ny = yfn[iFace] * area[iFace];
            RDouble nz = zfn[iFace] * area[iFace];

            RDouble qfc = phif[faceIndex];
            dphidx[le] += qfc * nx;
            dphidy[le] += qfc * ny;
            dphidz[le] += qfc * nz;
        }

        for (int iFace = nBoundFace; iFace < nTotalFace; ++iFace)
        {
            int le = leftCellofFace[iFace];
            int re = rightCellofFace[iFace];
            int faceIndex = iFace - nBoundFaces;

            RDouble nx = xfn[iFace] * area[iFace];
            RDouble ny = yfn[iFace] * area[iFace];
            RDouble nz = zfn[iFace] * area[iFace];

            RDouble qfc = phif[faceIndex];
            dphidx[le] += qfc * nx;
            dphidy[le] += qfc * ny;
            dphidz[le] += qfc * nz;
            dphidx[re] -= qfc * nx;
            dphidy[re] -= qfc * ny;
            dphidz[re] -= qfc * nz;
        }

        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            dphidx[iCell] *= one / vol[iCell];
            dphidy[iCell] *= one / vol[iCell];
            dphidz[iCell] *= one / vol[iCell];
        }
    }

    void IncomCalculator::GradientCalculation(Grid *gridIn, string phiName, string dphidxName, string dphidyName, string dphidzName, bool isComm)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);

        RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
        RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr(dphidxName));
        RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr(dphidyName));
        RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr(dphidzName));

        string GradCalcMethod = GlobalDataBase::GetStrParaFromDB("GradCalcMethod");
        if (GradCalcMethod == "GAUSS")
        {
            RDouble *faceWeightOfLeftCell= reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));
            grid->CompGradientGGCell2(phi, faceWeightOfLeftCell, dphidx, dphidy, dphidz);

            int isSkewness = GlobalDataBase::GetIntParaFromDB("isSkewness");
            if (isSkewness)
            {
                SkewnessCorrection(grid, phi, dphidx, dphidy, dphidz);
            }
        }
        else if (GradCalcMethod == "LSQ")
        {
            grid->ComputeWeight();
            grid->CompGradientLSQ(phi, dphidx, dphidy, dphidz);
        }

        CommunicateAnInterfaceVar(dphidx);
        CommunicateAnInterfaceVar(dphidy);
        CommunicateAnInterfaceVar(dphidz);
    }

    void IncomCalculator::calcFaceWeight(Grid *gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));
        int nBoundFace = grid->GetNBoundFace();
        int nTotalFace = grid->GetNTotalFace();
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();

        RDouble *xfc = grid->GetFaceCenterX();
        RDouble *yfc = grid->GetFaceCenterY();
        RDouble *zfc = grid->GetFaceCenterZ();
        RDouble *xcc = grid->GetCellCenterX();
        RDouble *ycc = grid->GetCellCenterY();
        RDouble *zcc = grid->GetCellCenterZ();

        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            int le = leftCellOfFace [ iFace ];
            int re = rightCellOfFace[ iFace ];

            RDouble dxL = xfc[iFace] - xcc[le];
            RDouble dyL = yfc[iFace] - ycc[le];
            RDouble dzL = zfc[iFace] - zcc[le];

            RDouble dxR = xfc[iFace] - xcc[re];
            RDouble dyR = yfc[iFace] - ycc[re];
            RDouble dzR = zfc[iFace] - zcc[re];

            // Left
            RDouble delta1 = DISTANCE(dxL, dyL, dzL);
            RDouble delta2 = DISTANCE(dxR, dyR, dzR);
            RDouble delta = one / (delta1 + delta2 + SMALL);

            faceWeightOfLeftCell[iFace] = delta2 * delta;
        }

        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
            int bcType = bcRegion->GetBCType();

            if (!IsInterface(bcType))
            {
                faceWeightOfLeftCell[iFace] = 1.0;
            }
        }
    }

    void IncomCalculator::CompressAnInterfaceVar(RDouble *field)
    {
        int currentProcessor = PHMPI::GetCurrentProcessorID();
        int numberOfSimulationZones = PHMPI::GetNumberofGlobalZones();
        receivedDataBuffer.resize(numberOfSimulationZones);
        sendDataBuffer.resize(numberOfSimulationZones);

        //! Step 0: Compressing data firstly.
        for (int iZone = 0; iZone < numberOfSimulationZones; ++ iZone)
        {
            int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

            ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
            int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
            if (currentProcessor == sendProcessor)
            {
                sendDataBuffer[iZone].resize(numberOfNeighbor);
            }

            for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
            {
                int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

                int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

                if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                {
                    continue;
                }

                //! Allocate and compress the buffers for sending.
                DataContainer *sendBuffer = 0;
                if (currentProcessor == sendProcessor)
                {
                    sendBuffer = new DataContainer;
                    sendDataBuffer[iZone][iNeighbor] = sendBuffer;

                    int localZoneID = PHMPI::GetGlobalZoneIDToLocalZoneID(iZone);

                    //! Compress the send information into the actkey.
                    Grid *grid = PHSPACE::GetGrid(iZone, 0);
                    CompressAnInterfaceVar(sendBuffer, field, grid, iZone, neighborZone);
                }
            }
        }
    }

    void IncomCalculator::CommunicateAnInterfaceVar()
    {
        int currentProcessor = PHMPI::GetCurrentProcessorID();
        int numberOfSimulationZones = PHMPI::GetNumberofGlobalZones();
        //! Step 1: Communication.
        for (int iZone = 0; iZone < numberOfSimulationZones; ++ iZone)
        {
            ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
            int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

            //! Communicating.
            for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
            {
                int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

                int sendProcessor    = PHMPI::GetZoneProcessorID(iZone);
                int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

                if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                {
                    continue;
                }

                //! Allocate the buffers for sending and receiving.
                DataContainer *receivedBuffer = 0;
                DataContainer *sendBuffer = 0;
                if (currentProcessor == receiveProcessor)
                {
                    receivedBuffer = new DataContainer;
                    receivedDataBuffer[iZone].push_back(receivedBuffer);
                }
                if (currentProcessor == sendProcessor)
                {
                    sendBuffer = sendDataBuffer[iZone][iNeighbor];
                }

                int tag = iZone;

                //! Communication.
                if (currentProcessor == sendProcessor)
                {
                    if (sendProcessor != receiveProcessor)
                    {
                        streamsize nlen = sendBuffer->Size();

                        // Send the data to neighbors.
                        send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                    }
                }
                if (currentProcessor == receiveProcessor)
                {
                    //! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                    int neighborOrer = -1;
                    ZoneNeighbor *globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                    int numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                    for (int neighborID = 0; neighborID < numberOfNeighborTemp; ++ neighborID)
                    {
                        int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(neighborID);
                        if (neighborZoneTemp == iZone)
                        {
                            neighborOrer = neighborID;
                            break;
                        }
                    }

                    ASSERT(neighborOrer != -1);

                    if (sendProcessor != receiveProcessor)
                    {
                        //! The data length of the received data is same to the length that send to the neighbor.
                        CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->Size();

                        //! Receive data from neighbors.
                        receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                    }
                    else
                    {
                        SWAP(receivedDataBuffer[iZone].back(), sendDataBuffer[iZone][iNeighbor]);
                    }
                }
            }
        }

        //! Step2: MPI waiting for Non-blocking communications.
        if (PHMPI::GetNumberOfProcessor() > 1)
        {
            PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
        }
    }

    void IncomCalculator::DecompressAnInterfaceVar(RDouble *field)
    {
        int currentProcessor = PHMPI::GetCurrentProcessorID();
        int numberOfSimulationZones = PHMPI::GetNumberofGlobalZones();
        //! Step3: Translate the data container.
        for (int iZone = 0; iZone < numberOfSimulationZones; ++ iZone)
        {
            ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
            int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

            int count = 0;
            for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
            {
                int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

                int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

                if (currentProcessor == receiveProcessor)
                {
                    DataContainer *receiveData = receivedDataBuffer[iZone][count];

                    //! Because get into here when currentProcessor==receiveProcessor,
                    //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                    //! Decompress the interface data from data-container.
                    int localZoneID = PHMPI::GetGlobalZoneIDToLocalZoneID(neighborZone);

                    Grid *grid = PHSPACE::GetGrid(neighborZone, 0);
                    DecompressAnInterfaceVar(receiveData, field, grid, neighborZone, iZone);

                    ++ count;
                }
            }
        }

        //! Step4: Free the buffers.
        for (int iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
        {
            for (int jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
            {
                delete receivedDataBuffer[iDim][jDim];
            }
        }
        receivedDataBuffer.clear();
        for (int iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
        {
            for (int jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
            {
                delete sendDataBuffer[iDim][jDim];
            }
        }
        sendDataBuffer.clear();
        requestContainer.clear();
    }

    void IncomCalculator::CommunicateAnInterfaceVar(RDouble *field)
    {
        CompressAnInterfaceVar(field);
        CommunicateAnInterfaceVar();
        DecompressAnInterfaceVar(field);
    }

    void IncomCalculator::CompressAnInterfaceVar(DataContainer *&dataContainer, RDouble *field, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex)
    {
        UnstructGrid * grid = UnstructGridCast(gridIn);

        InterfaceInfo * interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation) return;

        int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
        int * interfaceIndexContainerForSend         = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
        dataContainer->MoveToBegin();

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int sourceCell;
                int iFace = interfaceIndexContainerForSend[ iLocalFace ];
                grid->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
                PHWrite(dataContainer, field[sourceCell]);
            }
        }
    }

    void IncomCalculator::DecompressAnInterfaceVar(DataContainer *&dataContainer, RDouble *field, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex)
    {
        UnstructGrid * grid = UnstructGridCast(gridIn);

        InterfaceInfo * interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation) return;

        int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
        int * interfaceIndexContainerForReceive      = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);

        dataContainer->MoveToBegin();

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[ iLocalFace ];
                int targetCell;
                grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);

                PHRead(dataContainer, field[targetCell]);
            }
        }
    }

    void IncomCalculator::InitializeMatrixACoeff(Grid *gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn); 
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalFace = grid->GetNTotalFace();
        RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
        RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("upperMatrixCoeff"));
        RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("lowerMatrixCoeff"));

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            diagMatrixCoeff[iCell] = 0.0;
        }

        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            upperMatrixCoeff[iFace] = 0.0;
            lowerMatrixCoeff[iFace] = 0.0;
        }
    }

    void IncomCalculator::InitializeBCoeff(Grid *gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn); 
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalFace = grid->GetNTotalFace();
        RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            bCoeff[iCell] = 0.0;
        }
    }



}