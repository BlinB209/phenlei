#include "HODefine.h"
#include "HOGeometryStructure.h"
#include "HOStandardElement.h"
#include <vector>
#include "Geo_UnstructBC.h"
#include "HOSolverUnstructParam.h"

using namespace std;
using namespace HOUnstruct;
using namespace PHENGLEI;

namespace HOUnstruct
{

HighOrderCellSol::HighOrderCellSol(void) 
{
    nPolySol                = 0;
    nFVPolySol              = 0;
    nDOFsSol                = 0;

    integNumberOfGaussPoint = 0;
    ID                      = 0;
    offsetDOFsSolGlobal     = 0;
    offsetDOFsSolLocal      = 0;
}

void HighOrderCellSol::ComputeCoorOfIntegrationPoint(int nodenumber, const unsigned long *nodeID, 
    const RDouble *pointX, const RDouble *pointY, const RDouble *pointZ, const vector < vector < RDouble > > & shapeFunctionVolumPoint)
{
    //The dof of shape function is equal to the number of nodes;
    const int dofOfShapeFunction = nodenumber;
    for (int iGaussPoint = 0; iGaussPoint < integNumberOfGaussPoint; ++iGaussPoint)
    {
        coorXIntegrationPoints[iGaussPoint] = 0.0;
        coorYIntegrationPoints[iGaussPoint] = 0.0;
        coorZIntegrationPoints[iGaussPoint] = 0.0;

        for (int iDofShapeFun = 0; iDofShapeFun < dofOfShapeFunction; ++iDofShapeFun)
        {
            coorXIntegrationPoints[iGaussPoint] += pointX[nodeID[iDofShapeFun]] * shapeFunctionVolumPoint[iGaussPoint][iDofShapeFun];
            coorYIntegrationPoints[iGaussPoint] += pointY[nodeID[iDofShapeFun]] * shapeFunctionVolumPoint[iGaussPoint][iDofShapeFun];
            coorZIntegrationPoints[iGaussPoint] += pointZ[nodeID[iDofShapeFun]] * shapeFunctionVolumPoint[iGaussPoint][iDofShapeFun];
        }
    }
}

void HighOrderCellSol::ComputeOriginalMassmatrix(HighOrderMatrix & taylorBasisMassMatrix, const RDouble *weightCoeff)
{
    //! compute upper mass matrix
    RDouble volumJacobian = 0.0;
    for (int iGaussPoint = 0; iGaussPoint < integNumberOfGaussPoint; ++iGaussPoint)
    {
        volumJacobian += weightCoeff[iGaussPoint] * JacDetCellIntegration[iGaussPoint];
    }

    for (int idof = 0; idof < nDOFsSol; ++idof)
    {
        for (int jdof = idof; jdof < nDOFsSol; ++jdof)
        {
            for (int iGaussPoint = 0; iGaussPoint < integNumberOfGaussPoint; ++iGaussPoint)
            {
                taylorBasisMassMatrix(idof, jdof) += weightCoeff[iGaussPoint] * JacDetCellIntegration[iGaussPoint]
                * basisFuncCellIntegration[iGaussPoint][idof] * basisFuncCellIntegration[iGaussPoint][jdof];
            }
            taylorBasisMassMatrix(idof, jdof) /= volumJacobian;
        }
    }

    //! Get lower mass matrix
    for (int idof = 1; idof < nDOFsSol; ++idof)
    {
        for (int jdof = 0; jdof < idof; ++jdof)
        {
            taylorBasisMassMatrix(idof, jdof) = taylorBasisMassMatrix(jdof, idof);
        }
    }
}

void HighOrderCellSol::ComputeOriginalBasisFunction(int currentOrder, RDouble lenScaleOfCell, const RDouble Xcc, const RDouble Ycc, const RDouble Zcc)
{
    //! The dof of shape function is equal to the number of nodes;
    RDouble refCoord[3];
    RDouble oneDLenScal = 1.0/lenScaleOfCell;
    //! currentOrder is the accurry of polynomial 
    if (currentOrder == 0)
    {
        for (int iGaussPoint = 0; iGaussPoint < integNumberOfGaussPoint; ++iGaussPoint)
        {
            basisFuncCellIntegration[iGaussPoint][0] = 1.0;  
        }
    }
    else if (currentOrder == 1)
    {
        for (int iGaussPoint = 0; iGaussPoint < integNumberOfGaussPoint; ++iGaussPoint)
        {
            refCoord[0] = (coorXIntegrationPoints[iGaussPoint] - Xcc) * oneDLenScal;
            refCoord[1] = (coorYIntegrationPoints[iGaussPoint] - Ycc) * oneDLenScal;
            refCoord[2] = (coorZIntegrationPoints[iGaussPoint] - Zcc) * oneDLenScal;

            basisFuncCellIntegration[iGaussPoint][0] = 1.0;  
            basisFuncCellIntegration[iGaussPoint][1] = refCoord[0];
            basisFuncCellIntegration[iGaussPoint][2] = refCoord[1];
            basisFuncCellIntegration[iGaussPoint][3] = refCoord[2];
        }
    }
    else if (currentOrder == 2)
    {
        for (int iGaussPoint = 0; iGaussPoint < integNumberOfGaussPoint; ++iGaussPoint)
        {
            refCoord[0] = (coorXIntegrationPoints[iGaussPoint] - Xcc) * oneDLenScal;
            refCoord[1] = (coorYIntegrationPoints[iGaussPoint] - Ycc) * oneDLenScal;
            refCoord[2] = (coorZIntegrationPoints[iGaussPoint] - Zcc) * oneDLenScal;

            basisFuncCellIntegration[iGaussPoint][0] = 1.0;  
            basisFuncCellIntegration[iGaussPoint][1] = refCoord[0];
            basisFuncCellIntegration[iGaussPoint][2] = refCoord[1];
            basisFuncCellIntegration[iGaussPoint][3] = refCoord[2];

            basisFuncCellIntegration[iGaussPoint][4] = refCoord[0] * refCoord[0];
            basisFuncCellIntegration[iGaussPoint][5] = refCoord[1] * refCoord[1];
            basisFuncCellIntegration[iGaussPoint][6] = refCoord[2] * refCoord[2];
            basisFuncCellIntegration[iGaussPoint][7] = refCoord[0] * refCoord[1];
            basisFuncCellIntegration[iGaussPoint][8] = refCoord[0] * refCoord[2];
            basisFuncCellIntegration[iGaussPoint][9] = refCoord[1] * refCoord[2];
        }
    }
}

void HighOrderCellSol::ComputeJacobianDeterminant(int iGaussPoint, int iCellVTKType, int nodeNumber, const unsigned long *nodeID, 
    const RDouble *pointX, const RDouble *pointY, const RDouble *pointZ, RDouble *gradXiShapeFunction, RDouble *gradEtShapeFunction, RDouble *gradZtShapeFunction)
{
    RDouble dxDxi = 0.0;
    RDouble dxDet = 0.0;
    RDouble dxDzt = 0.0;
    RDouble dyDxi = 0.0;
    RDouble dyDet = 0.0;
    RDouble dyDzt = 0.0;
    RDouble dzDxi = 0.0;
    RDouble dzDet = 0.0;
    RDouble dzDzt = 0.0;

    //The dof of shape function is equal to the number of nodes;
    JacDetCellIntegration[iGaussPoint] = 0.0;
    for (int localNode = 0; localNode < nodeNumber; ++localNode)
    {
        dxDxi += pointX[nodeID[localNode]] * gradXiShapeFunction[localNode];
        dxDet += pointX[nodeID[localNode]] * gradEtShapeFunction[localNode];
        dxDzt += pointX[nodeID[localNode]] * gradZtShapeFunction[localNode];

        dyDxi += pointY[nodeID[localNode]] * gradXiShapeFunction[localNode];
        dyDet += pointY[nodeID[localNode]] * gradEtShapeFunction[localNode];
        dyDzt += pointY[nodeID[localNode]] * gradZtShapeFunction[localNode];

        dzDxi += pointZ[nodeID[localNode]] * gradXiShapeFunction[localNode];
        dzDet += pointZ[nodeID[localNode]] * gradEtShapeFunction[localNode];
        dzDzt += pointZ[nodeID[localNode]] * gradZtShapeFunction[localNode];
    }

    JacDetCellIntegration[iGaussPoint] = dxDxi * dyDet * dzDzt + dxDzt * dyDxi * dzDet + dxDet * dyDzt * dzDxi - 
        dxDzt * dyDet * dzDxi - dxDxi * dyDzt * dzDet - dxDet * dyDxi * dzDzt;

    switch(iCellVTKType)
    {
        case HOUnstruct::TETRAHEDRON:
            JacDetCellIntegration[iGaussPoint] /= 6.0;
            break;
        case HOUnstruct::PYRAMID:
            JacDetCellIntegration[iGaussPoint] /= (9.0/5.0);
            break;
        case HOUnstruct::PRISM:
            JacDetCellIntegration[iGaussPoint] /= 2.0;
            break;
        case HOUnstruct::HEXAHEDRON:
            JacDetCellIntegration[iGaussPoint] *= 1.0;
            break;
        default:
            break;
    }
}

void HighOrderCellSol::ComputeOrthogonalBasisFunction(int degreeOfFreedom, HighOrderMatrix & orthogonalizationMatrix)
{
    //20 is the dof for P3,if
    if (degreeOfFreedom > 20)
    {
        cout << "Currect code only suports fourth-order accuracy,if the accuracy is higher, the number of this function needs to be modified" << endl;
    }

    RDouble tmpBasisFunction[20];
    for (int iGaussPoint = 0; iGaussPoint < integNumberOfGaussPoint; ++iGaussPoint)
    {
        for (int idof = 0; idof < nDOFsSol; ++idof)
        {
            tmpBasisFunction[idof] = 0.0;
            for (int jdof = 0; jdof < nDOFsSol; ++jdof)
            {
                tmpBasisFunction[idof] += orthogonalizationMatrix(idof,jdof) * basisFuncCellIntegration[iGaussPoint][jdof]; 
            }
        }
        for (int idof = 0; idof < nDOFsSol; ++idof)
        {
            basisFuncCellIntegration[iGaussPoint][idof] = tmpBasisFunction[idof];
        }
    }
}

void HighOrderCellSol::TestOrthogonalBasisFunction(const RDouble *weight)
{
    //! 20 is the dof for P3,if
    RDouble tmpMatrix[10][10];
    for (int n=0;n<10;n++)
    {
        for (int m=0;m<10;m++)
        {
            tmpMatrix[n][m] = 0.0;
            for (int iGaussPoint = 0; iGaussPoint < integNumberOfGaussPoint; ++iGaussPoint)
            {
                tmpMatrix[n][m] += weight[iGaussPoint] *JacDetCellIntegration[iGaussPoint] * basisFuncCellIntegration[iGaussPoint][n] * basisFuncCellIntegration[iGaussPoint][m] ;
            } 
            cout << "n:  " << n << "  m:   " << m << "    tmpMatrix[n][m]   " << tmpMatrix[n][m] << endl;
        }
    }
}

void HighOrderCellSol::ComputeGradOfBasisFunctionIntegrationPoint(HighOrderMatrix & orthogonalizationMatrix, int iGaussPoint, 
    RDouble lenScaleOfCell)
{
    //! dxBasisFuncCellIntegration;
    RDouble oneDLenScal = 1.0/lenScaleOfCell;
    RDouble dxdBasis[35] = {0};    //! number=10:P2,number=20:P3
    RDouble dydBasis[35] = {0};
    RDouble dzdBasis[35] = {0};

    if (nPolySol == 0)
    {
        dxdBasis[0] = 0.0;
        dydBasis[0] = 0.0;
        dzdBasis[0] = 0.0;
    }
    else if (nPolySol == 1)
    {
        dxdBasis[0] = 0.0;
        dxdBasis[1] = oneDLenScal;
        dxdBasis[2] = 0.0;
        dxdBasis[3] = 0.0;

        dydBasis[0] = 0.0;
        dydBasis[1] = 0.0;
        dydBasis[2] = oneDLenScal;
        dydBasis[3] = 0.0;

        dzdBasis[0] = 0.0;
        dzdBasis[1] = 0.0;
        dzdBasis[2] = 0.0;
        dzdBasis[3] = oneDLenScal;
    }
    else if (nPolySol == 2)
    {
        dxdBasis[0] = 0.0;
        dxdBasis[1] = oneDLenScal;
        dxdBasis[2] = 0.0;
        dxdBasis[3] = 0.0;
        dxdBasis[4] = 2.0 * basisFuncCellIntegration[iGaussPoint][1] * oneDLenScal;
        dxdBasis[5] = 0.0;
        dxdBasis[6] = 0.0;
        dxdBasis[7] =       basisFuncCellIntegration[iGaussPoint][2] * oneDLenScal;
        dxdBasis[8] =       basisFuncCellIntegration[iGaussPoint][3] * oneDLenScal;
        dxdBasis[9] = 0.0;

        dydBasis[0] = 0.0;
        dydBasis[1] = 0.0;
        dydBasis[2] = oneDLenScal;
        dydBasis[3] = 0.0;
        dydBasis[4] = 0.0;
        dydBasis[5] = 2.0 * basisFuncCellIntegration[iGaussPoint][2] * oneDLenScal;
        dydBasis[6] = 0.0;
        dydBasis[7] =       basisFuncCellIntegration[iGaussPoint][1] * oneDLenScal;
        dydBasis[8] = 0.0;
        dydBasis[9] =       basisFuncCellIntegration[iGaussPoint][3] * oneDLenScal;

        dzdBasis[0] = 0.0;
        dzdBasis[1] = 0.0;
        dzdBasis[2] = 0.0;
        dzdBasis[3] = oneDLenScal;
        dzdBasis[4] = 0.0;
        dzdBasis[5] = 0.0;
        dzdBasis[6] = 2.0 * basisFuncCellIntegration[iGaussPoint][3] * oneDLenScal;
        dzdBasis[7] = 0.0;
        dzdBasis[8] =       basisFuncCellIntegration[iGaussPoint][1] * oneDLenScal;
        dzdBasis[9] =       basisFuncCellIntegration[iGaussPoint][2] * oneDLenScal;
    }

    for (int idof = 0; idof < nDOFsSol; ++idof)
    {
        dxBasisFuncCellIntegration[iGaussPoint][idof] = 0.0;
        dyBasisFuncCellIntegration[iGaussPoint][idof] = 0.0;
        dzBasisFuncCellIntegration[iGaussPoint][idof] = 0.0;

        for (int jdof = 0; jdof < nDOFsSol;++jdof)
        {
            dxBasisFuncCellIntegration[iGaussPoint][idof] += orthogonalizationMatrix(idof, jdof) * dxdBasis[jdof];
            dyBasisFuncCellIntegration[iGaussPoint][idof] += orthogonalizationMatrix(idof, jdof) * dydBasis[jdof];
            dzBasisFuncCellIntegration[iGaussPoint][idof] += orthogonalizationMatrix(idof, jdof) * dzdBasis[jdof];
        }
    }
}

void HighOrderFaceSol::ComputeCoorOfIntegrationPoint(int nodenumber, const unsigned long *nodeID, const RDouble *pointX, const RDouble *pointY, const RDouble *pointZ,
    const vector < vector < RDouble > > & shapeFunctionFacePoint)
{
    //! The dof of shape function is equal to the number of nodes;
    const int dofOfShapeFunction = nodenumber;

    for (int iGaussPoint = 0; iGaussPoint < integNumberOfGaussPoint; ++iGaussPoint)
    {
        coorXIntegrationPoints[iGaussPoint] = 0.0;
        coorYIntegrationPoints[iGaussPoint] = 0.0;
        coorZIntegrationPoints[iGaussPoint] = 0.0;

        for (int iDofShapeFun = 0; iDofShapeFun < dofOfShapeFunction; ++iDofShapeFun)
        {
            coorXIntegrationPoints[iGaussPoint] += pointX[nodeID[iDofShapeFun]]*shapeFunctionFacePoint[iGaussPoint][iDofShapeFun];
            coorYIntegrationPoints[iGaussPoint] += pointY[nodeID[iDofShapeFun]]*shapeFunctionFacePoint[iGaussPoint][iDofShapeFun];
            coorZIntegrationPoints[iGaussPoint] += pointZ[nodeID[iDofShapeFun]]*shapeFunctionFacePoint[iGaussPoint][iDofShapeFun];
        }    
    }
}

void HighOrderFaceSol::ComputeBasisFunction(int iGaussPoint, int leftCell, int rightCell, int nTotalCell, int nPloyOfLeftCell, int nPloyOfRightCell, 
    RDouble lenScaleOfLeftCell, RDouble lenScaleOfRightCell, const RDouble *Xcc, const RDouble *Ycc, const RDouble *Zcc, 
    const RDouble *leftOrthMatrix, const RDouble *rightOrthMatrix)
{
    //! leftBasisFuncFaceIntegration;
    RDouble oneDLeftLenScal  = 1.0/lenScaleOfLeftCell;
    RDouble oneDRightLenScal = 1.0/lenScaleOfRightCell;
    //! 10 if only for P1/P2, if accuracy is P3, 10 should be 20 
    RDouble refCoord[3],leftBasis[10], rightBasis[10], orthMatrix[10][10];
    int leftDof = 0, rightDof = 0, counter = 0;
    if (nPloyOfLeftCell == 0)
    {
        leftDof = 1;

        leftBasis[0] = 1.0;
    }
    else if (nPloyOfLeftCell == 1)
    {
        leftDof = 4;

        leftBasis[0] = 1.0;
        leftBasis[1] = (coorXIntegrationPoints[iGaussPoint] - Xcc[leftCell]) * oneDLeftLenScal;
        leftBasis[2] = (coorYIntegrationPoints[iGaussPoint] - Ycc[leftCell]) * oneDLeftLenScal;
        leftBasis[3] = (coorZIntegrationPoints[iGaussPoint] - Zcc[leftCell]) * oneDLeftLenScal;
    }
    else if (nPloyOfLeftCell == 2)
    {
        leftDof = 10;
        refCoord[0] = (coorXIntegrationPoints[iGaussPoint] - Xcc[leftCell]) * oneDLeftLenScal;
        refCoord[1] = (coorYIntegrationPoints[iGaussPoint] - Ycc[leftCell]) * oneDLeftLenScal;
        refCoord[2] = (coorZIntegrationPoints[iGaussPoint] - Zcc[leftCell]) * oneDLeftLenScal;

        leftBasis[0] = 1.0;
        leftBasis[1] = refCoord[0];
        leftBasis[2] = refCoord[1];
        leftBasis[3] = refCoord[2];

        leftBasis[4] = refCoord[0] * refCoord[0];
        leftBasis[5] = refCoord[1] * refCoord[1];
        leftBasis[6] = refCoord[2] * refCoord[2];
        leftBasis[7] = refCoord[0] * refCoord[1];
        leftBasis[8] = refCoord[0] * refCoord[2];
        leftBasis[9] = refCoord[1] * refCoord[2];
    }

    //! left side of face
    counter = 0;
    for (int iDof = 0; iDof < leftDof; ++iDof)
    {
        for (int jDof = 0; jDof < leftDof; ++jDof)
        {
            orthMatrix[iDof][jDof] = leftOrthMatrix[counter ++];
        }
    }

    for (int iDof = 0; iDof < leftDof; ++iDof)
    {
        leftBasisFuncFaceIntegration[iGaussPoint][iDof] = 0.0;
        for (int jDof = 0; jDof < leftDof; ++jDof)
        {
            leftBasisFuncFaceIntegration[iGaussPoint][iDof] += orthMatrix[iDof][jDof] * leftBasis[jDof];
        }
    }

    if (rightCell < nTotalCell)
    {
        if (nPloyOfRightCell == 0) 
        {
            rightDof = 1;

            rightBasis[0] = 1.0;
        }
        else if (nPloyOfRightCell == 1) 
        {
            rightDof = 4;

            rightBasis[0] = 1.0;
            rightBasis[1] = (coorXIntegrationPoints[iGaussPoint] - Xcc[rightCell]) * oneDRightLenScal;
            rightBasis[2] = (coorYIntegrationPoints[iGaussPoint] - Ycc[rightCell]) * oneDRightLenScal;
            rightBasis[3] = (coorZIntegrationPoints[iGaussPoint] - Zcc[rightCell]) * oneDRightLenScal;
        }
        else if (nPloyOfRightCell == 2)
        {
            rightDof = 10;

            refCoord[0] = (coorXIntegrationPoints[iGaussPoint] - Xcc[rightCell]) * oneDRightLenScal;
            refCoord[1] = (coorYIntegrationPoints[iGaussPoint] - Ycc[rightCell]) * oneDRightLenScal;
            refCoord[2] = (coorZIntegrationPoints[iGaussPoint] - Zcc[rightCell]) * oneDRightLenScal;

            rightBasis[0] = 1.0;
            rightBasis[1] = refCoord[0];
            rightBasis[2] = refCoord[1];
            rightBasis[3] = refCoord[2];

            rightBasis[4] = refCoord[0] * refCoord[0];
            rightBasis[5] = refCoord[1] * refCoord[1];
            rightBasis[6] = refCoord[2] * refCoord[2];
            rightBasis[7] = refCoord[0] * refCoord[1];
            rightBasis[8] = refCoord[0] * refCoord[2];
            rightBasis[9] = refCoord[1] * refCoord[2];            
        }

        //! right side of face
        counter = 0;
        for (int iDof = 0; iDof < rightDof; ++iDof)
        {
            for (int jDof = 0; jDof < rightDof; ++jDof)
            {
                orthMatrix[iDof][jDof] = rightOrthMatrix[counter ++];
            }
        }
        for (int iDof = 0; iDof < rightDof; ++iDof)
        {
            rightBasisFuncFaceIntegration[iGaussPoint][iDof] = 0.0;
            for (int jDof = 0; jDof < rightDof; ++jDof)
            {
                rightBasisFuncFaceIntegration[iGaussPoint][iDof] += orthMatrix[iDof][jDof] * rightBasis[jDof];
            }
        }
    }
    else
    {
        for (int iDof = 0; iDof < rightDof; ++iDof)
        {
            rightBasisFuncFaceIntegration[iGaussPoint][iDof] = leftBasisFuncFaceIntegration[iGaussPoint][iDof];
        }
    }
}

void HighOrderFaceSol::ComputeMetricNormalsFace(int iGaussPoint, int iFaceVTKType, int nodenumber, const unsigned long * nodeID, 
    const RDouble * pointX, const RDouble * pointY, const RDouble * pointZ, RDouble * gradXiShapeFunction, RDouble * gradEtShapeFunction)
{
    RDouble dxDxi = 0.0;
    RDouble dxDet = 0.0;
    RDouble dyDxi = 0.0;
    RDouble dyDet = 0.0;
    RDouble dzDxi = 0.0;
    RDouble dzDet = 0.0;
    RDouble totalJacobin = 0.0;

    //! The dof of shape function is equal to the number of nodes;
    metricNormalsFace[iGaussPoint][0] = 0.0;
    metricNormalsFace[iGaussPoint][1] = 0.0;
    metricNormalsFace[iGaussPoint][2] = 0.0;

    for (int iLocalNode = 0; iLocalNode < nodenumber; ++ iLocalNode)
    {
        dxDxi += pointX[nodeID[iLocalNode]] * gradXiShapeFunction[iLocalNode];
        dxDet += pointX[nodeID[iLocalNode]] * gradEtShapeFunction[iLocalNode];

        dyDxi += pointY[nodeID[iLocalNode]] * gradXiShapeFunction[iLocalNode];
        dyDet += pointY[nodeID[iLocalNode]] * gradEtShapeFunction[iLocalNode];

        dzDxi += pointZ[nodeID[iLocalNode]] * gradXiShapeFunction[iLocalNode];
        dzDet += pointZ[nodeID[iLocalNode]] * gradEtShapeFunction[iLocalNode];
    }

    metricNormalsFace[iGaussPoint][0] = dyDxi * dzDet - dyDet * dzDxi;
    metricNormalsFace[iGaussPoint][1] = dzDxi * dxDet - dzDet * dxDxi;
    metricNormalsFace[iGaussPoint][2] = dxDxi * dyDet - dxDet * dyDxi;
    totalJacobin = sqrt(pow(metricNormalsFace[iGaussPoint][0], 2) + pow(metricNormalsFace[iGaussPoint][1], 2) + pow(metricNormalsFace[iGaussPoint][2], 2));

    metricNormalsFace[iGaussPoint][0] /= totalJacobin;
    metricNormalsFace[iGaussPoint][1] /= totalJacobin;
    metricNormalsFace[iGaussPoint][2] /= totalJacobin;

    switch(iFaceVTKType)
    {
        case HOUnstruct::TRIANGLE:
            JacDetFaceIntegration[iGaussPoint] = totalJacobin/2.0;
            break;
        case HOUnstruct::QUADRILATERAL:
            JacDetFaceIntegration[iGaussPoint] = totalJacobin/1.0;
            break;
        default:
            break;
    }
}

//! @Application memory for high order cell and  copy the geometry of grid to high order cell,containing 
//! @param[in]   grid                      the input unstructgrid
//! @param[out]                            initialize HighOrderCellSol and build HighOrderStandardCell   
void HighOrderGrid::CopyGeometryInfoToHighOrderCell(UnstructGrid * grid)
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    highOrderCells.resize(nTotal);

    const int *pFaceNumberOfCells = grid->GetFaceNumberOfEachCell();
    int **ppFacesOfCells = grid->GetCell2Face();
    int *pNodeNumberOfCell = grid->GetNodeNumberOfEachCell();
    //int ** ppNodesOfCells = grid->GetCell2NodeArray();
    int *pNodesOfCells = grid->GetCell2Node();
    int cellToNodeIndex = 0;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        highOrderCells[iCell].ID = iCell;

        highOrderCells[iCell].nFaces = static_cast<unsigned short>(pFaceNumberOfCells[iCell]);
        for (int ilocalface = 0; ilocalface < highOrderCells[iCell].nFaces; ++ ilocalface)
        {
            highOrderCells[iCell].faceIDs.push_back(ppFacesOfCells[iCell][ilocalface]);
        }

        highOrderCells[iCell].nNodes = static_cast<unsigned short>(pNodeNumberOfCell[iCell]);
        for (int ilocalnode = 0; ilocalnode < highOrderCells[iCell].nNodes; ++ ilocalnode)
        {
            //highOrderCells[iCell].nodeIDsGrid.push_back(ppNodesOfCells[iCell][ilocalnode]);
            highOrderCells[iCell].nodeIDs.push_back(pNodesOfCells[cellToNodeIndex++]);
        }
        if (highOrderCells[iCell].nNodes == 4)
        {
            highOrderCells[iCell].nPolyGrid = 1;
            highOrderCells[iCell].VTK_Type = HOUnstruct::TETRAHEDRON;
        }
        else if (highOrderCells[iCell].nNodes == 5)
        {
            highOrderCells[iCell].nPolyGrid = 1;
            highOrderCells[iCell].VTK_Type = HOUnstruct::PYRAMID;
        }
        else if (highOrderCells[iCell].nNodes == 6)
        {
            highOrderCells[iCell].nPolyGrid = 1;
            highOrderCells[iCell].VTK_Type = HOUnstruct::PRISM;
        }
        else if (highOrderCells[iCell].nNodes == 8)
        {
            highOrderCells[iCell].nPolyGrid = 1;
            highOrderCells[iCell].VTK_Type = HOUnstruct::HEXAHEDRON;
        }
        else if (highOrderCells[iCell].nNodes == 10)
        {
            highOrderCells[iCell].nPolyGrid = 2;
            highOrderCells[iCell].VTK_Type = HOUnstruct::TETRAHEDRON;
        }
        else if (highOrderCells[iCell].nNodes == 14)
        {
            highOrderCells[iCell].nPolyGrid = 1;
            highOrderCells[iCell].VTK_Type = HOUnstruct::PYRAMID;
        }
        else if (highOrderCells[iCell].nNodes == 18)
        {
            highOrderCells[iCell].nPolyGrid = 2;
            highOrderCells[iCell].VTK_Type = HOUnstruct::PRISM;
        }
        else if (highOrderCells[iCell].nNodes == 27)
        {
            highOrderCells[iCell].nPolyGrid = 2;
            highOrderCells[iCell].VTK_Type = HOUnstruct::HEXAHEDRON;
        }

        highOrderCells[iCell].nDOFsGrid = highOrderCells[iCell].nNodes;

        //highOrderCells[iCell].volume = grid->GetVolume()[iCell];
        //highOrderCells[iCell].centerCoord[0] = grid->GetCellCenterX()[iCell];
        //highOrderCells[iCell].centerCoord[1] = grid->GetCellCenterY()[iCell];
        //highOrderCells[iCell].centerCoord[2] = grid->GetCellCenterZ()[iCell];

        //highOrderCells[iCell].lenScale = pow(abs(highOrderCells[iCell].volume), 1.0/3.0);
    }

    //ghost cell value , we should attention parallel BC
    int *left_cell_of_face  = grid->GetLeftCellOfFace();
    int *right_cell_of_face = grid->GetRightCellOfFace();
    
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    int le = 0, re = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            le = left_cell_of_face [ iFace ];
            re = right_cell_of_face[ iFace ];
            // if is parallel BC, we should get the physical cell 
            if (IsInterface(bcType)) //interface
            {
                continue;
            }
            else
            {
                highOrderCells[re].ID = re;

                highOrderCells[re].nFaces = highOrderCells[le].nFaces;
                highOrderCells[re].nNodes = highOrderCells[le].nNodes;

                highOrderCells[re].nPolyGrid = highOrderCells[le].nPolyGrid;
                highOrderCells[re].VTK_Type = highOrderCells[le].VTK_Type;
                highOrderCells[re].nDOFsGrid =  highOrderCells[le].nDOFsGrid;

                //highOrderCells[re].volume = grid->GetVolume()[re];
                //highOrderCells[re].centerCoord[0] = grid->GetCellCenterX()[re];
                //highOrderCells[re].centerCoord[1] = grid->GetCellCenterY()[re];
                //highOrderCells[re].centerCoord[2] = grid->GetCellCenterZ()[re];

                //highOrderCells[re].lenScale = pow(abs(highOrderCells[re].volume), 1.0/3.0);
            }
        }
    }
}

//! @Application memory for high order face and  copy the geometry of grid to high order face,containing 
void HighOrderGrid::CopyGeometryInfoToHighOrderFace(UnstructGrid *grid)
{
    int nTotalFace = grid->GetNTotalFace();

    highOrderFaces.resize(nTotalFace);
    int * pNodeNumberOfFaces = grid->GetNodeNumberOfEachFace();
    int ** ppNodesOfFaces = grid->GetFace2NodeArray();

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        highOrderFaces[iFace].ID = iFace;

        highOrderFaces[iFace].nNodes = pNodeNumberOfFaces[iFace];
        for (unsigned int ilocalnode = 0; ilocalnode < highOrderFaces[iFace].nNodes; ++ ilocalnode)
        {
            highOrderFaces[iFace].nodeIDs.push_back(ppNodesOfFaces[iFace][ilocalnode]);
        }

        if (highOrderFaces[iFace].nNodes == 3)
        {
            highOrderFaces[iFace].nPolyGrid = 1;
            highOrderFaces[iFace].VTK_Type = TRIANGLE;
        }
        else if (highOrderFaces[iFace].nNodes == 4)
        {
            highOrderFaces[iFace].nPolyGrid = 1;
            highOrderFaces[iFace].VTK_Type = QUADRILATERAL;
        }
        else if (highOrderFaces[iFace].nNodes == 6)
        {
            highOrderFaces[iFace].nPolyGrid = 2;
            highOrderFaces[iFace].VTK_Type = TRIANGLE;
        }
        else if (highOrderFaces[iFace].nNodes == 9)
        {
            highOrderFaces[iFace].nPolyGrid = 2;
            highOrderFaces[iFace].VTK_Type = QUADRILATERAL;
        }

        highOrderFaces[iFace].nDOFsGrid = static_cast<unsigned short>(highOrderFaces[iFace].nNodes);
        highOrderFaces[iFace].cellID0 = grid->GetLeftCellOfFace()[iFace];
        highOrderFaces[iFace].cellID1 = grid->GetRightCellOfFace()[iFace];

        //highOrderFaces[iFace].area = grid->GetFaceArea()[iFace];
        //highOrderFaces[iFace].centerCoord[0] = grid->GetFaceCenterX()[iFace];
        //highOrderFaces[iFace].centerCoord[1] = grid->GetFaceCenterY()[iFace];
        //highOrderFaces[iFace].centerCoord[2] = grid->GetFaceCenterZ()[iFace];
    }
}

//! @Application memory for high order node and  copy the geometry of grid to high order node,containing 
void HighOrderGrid::CopyGeometryInfoToHighOrderNode(UnstructGrid *grid)
{
    int nTotalNode = grid->GetNTotalNode();

    highOrderNodes.resize(nTotalNode);

    int **ppCellsOfNodes = grid->GetNode2CellArray();

    /*  Xu-Gang
    int * ppCellsOfNodes = grid->GetCell2Node();
    int * NodeNumberOfEachCell = grid->GetNodeNumberOfEachCell();
    */
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        highOrderNodes[iNode].ID = iNode;

        highOrderNodes[iNode].coor[0] = grid->GetX()[iNode];
        highOrderNodes[iNode].coor[1] = grid->GetY()[iNode];
        highOrderNodes[iNode].coor[2] = grid->GetZ()[iNode];

        highOrderNodes[iNode].cellIDsOfPoint.resize(grid->GetCellNumberOfEachNode()[iNode]);

        for (int ilocalcell = 0; ilocalcell < highOrderNodes[iNode].cellIDsOfPoint.size(); ++ ilocalcell)
        {
            highOrderNodes[iNode].cellIDsOfPoint[ilocalcell] = ppCellsOfNodes[iNode][ilocalcell];
        }
       // highOrderNodes[iNode].cellIDsOfPoint.resize(grid->GetCellNumberOfEachNode()[iNode]);
       //
       // for (int ilocalcell = 0; ilocalcell < highOrderNodes[iNode].cellIDsOfPoint.size(); ++ ilocalcell)
       // {
       //     highOrderNodes[iNode].cellIDsOfPoint[ilocalcell] = ppCellsOfNodes[iNode][ilocalcell];
       // }
    }
}

//! @Initialize HighOrderCellSol for different p_multigrid 
//! @param[in]   grid_in                   the input unstructgrid
//! @param[in]   dgSolOrder                the accuracy of DG solution
//! @param[in]   pMultiGrid                pmultigrid number
//! @param[out]                            initialize HighOrderCellSol and build HighOrderStandardCell
void HighOrderGrid::InitHighOrderCellSol(UnstructGrid *grid_in, int dgSolOrder, int pMultiGridLevel, int isUnsteady, int viscousType)
{     
    UnstructGrid *grid = UnstructGridCast(grid_in);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    const RDouble * pointX = grid->GetX();
    const RDouble * pointY = grid->GetY();
    const RDouble * pointZ = grid->GetZ();
    RDouble * Xcc = grid->GetCellCenterX();
    RDouble * Ycc = grid->GetCellCenterY();
    RDouble * Zcc = grid->GetCellCenterZ();
    RDouble * vol = grid->GetCellVolume();

    int pMultiGridOrder = 0;                            //! @p_multigrid order
    unsigned long positionOfDof = 0;                    //! @the dof position of curret cell,using for
    unsigned short int elementNumberOfOrth = 0;         //! @the number of orth matirx
    unsigned short int integNumberOfGauss = 0; 
    unsigned short int accuracyOfVolumeIntegration = 0; //! @the volume integation accuracy
    unsigned short int iCellVTKType = 0;
    unsigned short int nCellNodes = 0;
    unsigned short int degreeOfFreedom = 0;
    RDouble lenScaleOfCell = 0.0;

    positionOfDof = 0;

    highOrderCellsSol[pMultiGridLevel].resize(nTotal);

    int nDim = grid->GetDim();
    int neqn = 5;

    if (nDim == 3)
    {
        neqn = 5;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {

        highOrderCellsSol[pMultiGridLevel][iCell].ID = iCell;
        highOrderCellsSol[pMultiGridLevel][iCell].nPolySol = static_cast<unsigned short>(dgSolOrder - pMultiGridLevel); //polynomial order
        // the volume integation accuracy
        accuracyOfVolumeIntegration = 2 * highOrderCellsSol[pMultiGridLevel][iCell].nPolySol + 2;
        // unsigned short iCellVTKType = 
        iCellVTKType = highOrderCells[iCell].VTK_Type;
        //highOrderCells[iCell].VTK_Type = highOrderCells[iCell].VTK_Type;
        // get nodes of cell
        nCellNodes = highOrderCells[iCell].nNodes;
        unsigned short int iStandElemIndex = 0;                                   //! @the HigOrderStandardCell type

        for (iStandElemIndex = 0; iStandElemIndex < standardCellElements.size(); ++iStandElemIndex)
        {
            if (standardCellElements[iStandElemIndex].SameStandardElement(iCellVTKType, highOrderCells[iCell].nPolyGrid, accuracyOfVolumeIntegration))
            {
                break;
            }
        }

        if (iStandElemIndex == standardCellElements.size())
        {
            standardCellElements.push_back(HighOrderStandardCell(iCellVTKType, nCellNodes, accuracyOfVolumeIntegration));
        }
        //! highOrderCellsSol[pMultiGrid][iCell].nPolySol = dgSolOrder - pMultiGrid;
        pMultiGridOrder = highOrderCellsSol[pMultiGridLevel][iCell].nPolySol;
        highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol = static_cast<unsigned short>(((pMultiGridOrder + 1) * (pMultiGridOrder + 2) * (pMultiGridOrder + 3)) / 6);
        highOrderCellsSol[pMultiGridLevel][iCell].offsetDOFsSolLocal = positionOfDof;

        positionOfDof += highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol;

        highOrderCellsSol[pMultiGridLevel][iCell].integNumberOfGaussPoint = static_cast<unsigned short>(standardCellElements[iStandElemIndex].GetNumberOfIntegrationPoint());
        integNumberOfGauss = highOrderCellsSol[pMultiGridLevel][iCell].integNumberOfGaussPoint;
        const unsigned long * nodeID = highOrderCells[iCell].GetnodeIDsGrid();

        //apply memory to save coord of integration
        highOrderCellsSol[pMultiGridLevel][iCell].coorXIntegrationPoints.resize(integNumberOfGauss);
        highOrderCellsSol[pMultiGridLevel][iCell].coorYIntegrationPoints.resize(integNumberOfGauss);
        highOrderCellsSol[pMultiGridLevel][iCell].coorZIntegrationPoints.resize(integNumberOfGauss);

        //compute Coor of volume integration
        highOrderCellsSol[pMultiGridLevel][iCell].ComputeCoorOfIntegrationPoint(standardCellElements[iStandElemIndex].GetNnode(), nodeID, pointX, pointY, pointZ, 
            standardCellElements[iStandElemIndex].GetShapeFunctionsIntegration());

        //compute basis function of volume integration
        //1:compute Jacobian determinant of Guass point
        highOrderCellsSol[pMultiGridLevel][iCell].JacDetCellIntegration.resize(integNumberOfGauss);

        RDouble * gradXiShapeFunction;
        RDouble * gradEtShapeFunction;
        RDouble * gradZtShapeFunction;
        for (int iGaussPoint = 0; iGaussPoint < integNumberOfGauss; ++iGaussPoint)
        {
            gradXiShapeFunction = standardCellElements[iStandElemIndex].GetdXiShapeFunctionAtIntegration(iGaussPoint);
            gradEtShapeFunction = standardCellElements[iStandElemIndex].GetdEtShapeFunctionAtIntegration(iGaussPoint);
            gradZtShapeFunction = standardCellElements[iStandElemIndex].GetdZtShapeFunctionAtIntegration(iGaussPoint);
            highOrderCellsSol[pMultiGridLevel][iCell].ComputeJacobianDeterminant(iGaussPoint, iCellVTKType, standardCellElements[iStandElemIndex].GetNnode(), nodeID, pointX, pointY, pointZ,
                gradXiShapeFunction, gradEtShapeFunction, gradZtShapeFunction);
        }
        const RDouble *weightCoeff = 0;
        
        if (pMultiGridLevel == 0)
        {
            RDouble volumJacobian = 0.0;
            weightCoeff = standardCellElements[iStandElemIndex].GetWtIntegration();

            double xcc = 0, ycc = 0, zcc = 0;
            for (int iGaussPoint = 0; iGaussPoint < integNumberOfGauss; ++ iGaussPoint)
            {
                volumJacobian += weightCoeff[iGaussPoint] * highOrderCellsSol[pMultiGridLevel][iCell].JacDetCellIntegration[iGaussPoint];
                xcc += weightCoeff[iGaussPoint] * highOrderCellsSol[pMultiGridLevel][iCell].JacDetCellIntegration[iGaussPoint] * highOrderCellsSol[pMultiGridLevel][iCell].coorXIntegrationPoints[iGaussPoint];
                ycc += weightCoeff[iGaussPoint] * highOrderCellsSol[pMultiGridLevel][iCell].JacDetCellIntegration[iGaussPoint] * highOrderCellsSol[pMultiGridLevel][iCell].coorYIntegrationPoints[iGaussPoint];
                zcc += weightCoeff[iGaussPoint] * highOrderCellsSol[pMultiGridLevel][iCell].JacDetCellIntegration[iGaussPoint] * highOrderCellsSol[pMultiGridLevel][iCell].coorZIntegrationPoints[iGaussPoint];
            }

            vol[iCell] = volumJacobian;
            Xcc[iCell] = xcc / volumJacobian;
            Ycc[iCell] = ycc / volumJacobian;
            Zcc[iCell] = zcc / volumJacobian;
            
            highOrderCells[iCell].volume = volumJacobian;
            highOrderCells[iCell].centerCoord[0] = Xcc[iCell];
            highOrderCells[iCell].centerCoord[1] = Ycc[iCell];
            highOrderCells[iCell].centerCoord[2] = Zcc[iCell];

            highOrderCells[iCell].lenScale = pow(abs(highOrderCells[iCell].volume), 1.0/3.0);
        }

        //2:compute original basis function
        degreeOfFreedom = highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol;
        lenScaleOfCell = highOrderCells[iCell].lenScale;

        highOrderCellsSol[pMultiGridLevel][iCell].basisFuncCellIntegration.resize(integNumberOfGauss);
        for (int iGaussPoint = 0; iGaussPoint < integNumberOfGauss; ++iGaussPoint)
        {
            highOrderCellsSol[pMultiGridLevel][iCell].basisFuncCellIntegration[iGaussPoint].resize(degreeOfFreedom);
        }

        highOrderCellsSol[pMultiGridLevel][iCell].ComputeOriginalBasisFunction(pMultiGridOrder, lenScaleOfCell, Xcc[iCell], Ycc[iCell], Zcc[iCell]);

        //3:compute original massmatrix
        HighOrderMatrix taylorBasisMassMatrix(degreeOfFreedom, degreeOfFreedom);
        taylorBasisMassMatrix.setConstant(0.0);          
        highOrderCellsSol[pMultiGridLevel][iCell].ComputeOriginalMassmatrix(taylorBasisMassMatrix, weightCoeff);

        //4:compute orthogonalizationMatrix for basis function
        HighOrderMatrix orthogonalizationMatrix(degreeOfFreedom, degreeOfFreedom);
        orthogonalizationMatrix.setConstant(0.0);

        Orthogonalize(degreeOfFreedom, taylorBasisMassMatrix, orthogonalizationMatrix);

        //5:apply memory to save orthogonal matrix, get the orthogonalizedMatrix
        elementNumberOfOrth = highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol * highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol;
        highOrderCellsSol[pMultiGridLevel][iCell].orthogonalizedMatrix.resize(elementNumberOfOrth);
        unsigned short int numberOfMatrix = 0;
        for (int idof = 0; idof < degreeOfFreedom; ++idof)
        {
            for (int jdof = 0; jdof < degreeOfFreedom; ++jdof)
            {
                highOrderCellsSol[pMultiGridLevel][iCell].orthogonalizedMatrix[ numberOfMatrix++ ] = orthogonalizationMatrix(idof,jdof);
            }
        }

        //6:compute the grad of basis function  for volume gauss point.
        highOrderCellsSol[pMultiGridLevel][iCell].dxBasisFuncCellIntegration.resize(integNumberOfGauss);
        highOrderCellsSol[pMultiGridLevel][iCell].dyBasisFuncCellIntegration.resize(integNumberOfGauss);
        highOrderCellsSol[pMultiGridLevel][iCell].dzBasisFuncCellIntegration.resize(integNumberOfGauss);

        for (int iGaussPoint = 0; iGaussPoint < integNumberOfGauss; ++iGaussPoint)
        {                
            highOrderCellsSol[pMultiGridLevel][iCell].dxBasisFuncCellIntegration[iGaussPoint].resize(degreeOfFreedom);
            highOrderCellsSol[pMultiGridLevel][iCell].dyBasisFuncCellIntegration[iGaussPoint].resize(degreeOfFreedom);
            highOrderCellsSol[pMultiGridLevel][iCell].dzBasisFuncCellIntegration[iGaussPoint].resize(degreeOfFreedom);

            highOrderCellsSol[pMultiGridLevel][iCell].ComputeGradOfBasisFunctionIntegrationPoint(orthogonalizationMatrix, iGaussPoint, lenScaleOfCell);
        }
        
        //7:compute new basis function for volume gauss point,the new basis function is orthogonal
        highOrderCellsSol[pMultiGridLevel][iCell].ComputeOrthogonalBasisFunction(degreeOfFreedom, orthogonalizationMatrix);
        //highOrderCellsSol[pMultiGrid][iCell].TestOrthogonalBasisFunction(weightCoeff);

        //8:finish init HighOrderCellSol

        //allocate memory for q and res
        highOrderCellsSol[pMultiGridLevel][iCell].dofQ.resize(neqn);
        highOrderCellsSol[pMultiGridLevel][iCell].res.resize(neqn);
        highOrderCellsSol[pMultiGridLevel][iCell].q.resize(neqn);
        highOrderCellsSol[pMultiGridLevel][iCell].qAverage.resize(neqn);

        for (int ieqn = 0; ieqn < neqn; ++ieqn)
        {
            highOrderCellsSol[pMultiGridLevel][iCell].dofQ[ieqn].resize(highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol);
            highOrderCellsSol[pMultiGridLevel][iCell].res[ieqn].resize(highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol);
            highOrderCellsSol[pMultiGridLevel][iCell].q[ieqn].resize(highOrderCellsSol[pMultiGridLevel][iCell].integNumberOfGaussPoint);            
        }
    }

    for (int iCell = nTotalCell; iCell < nTotal; ++ iCell)
    {
        highOrderCellsSol[pMultiGridLevel][iCell].ID = iCell;
        highOrderCellsSol[pMultiGridLevel][iCell].nPolySol = static_cast<unsigned short>(dgSolOrder - pMultiGridLevel); //polynomial order

        //! highOrderCellsSol[pMultiGrid][iCell].nPolySol = dgSolOrder - pMultiGrid;
        pMultiGridOrder = highOrderCellsSol[pMultiGridLevel][iCell].nPolySol;
        highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol = static_cast<unsigned short>(((pMultiGridOrder + 1) * (pMultiGridOrder + 2) * (pMultiGridOrder + 3)) / 6);
        highOrderCellsSol[pMultiGridLevel][iCell].offsetDOFsSolLocal = positionOfDof;
        positionOfDof += highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol;

        //allocate memory for q and res
        highOrderCellsSol[pMultiGridLevel][iCell].dofQ.resize(neqn);
        highOrderCellsSol[pMultiGridLevel][iCell].qAverage.resize(neqn);

        for (int ieqn = 0; ieqn < neqn; ++ieqn)
        {
            highOrderCellsSol[pMultiGridLevel][iCell].dofQ[ieqn].resize(highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol);
        }
    }

    //HighOrderCellµÄtype
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        //! Determine the standardElementType, pMultiGrid
        iCellVTKType = highOrderCells[iCell].VTK_Type;
        int nPolyGrid = highOrderCells[iCell].nPolyGrid;

        for (int iStandElemIndex = 0; iStandElemIndex < standardCellElements.size(); ++iStandElemIndex)
        {
            if (iCellVTKType == standardCellElements[iStandElemIndex].GetVTK_Type() && 
                nPolyGrid    == standardCellElements[iStandElemIndex].GetNpoly())
            {
                highOrderCells[iCell].standardElementIndex = static_cast<unsigned short>(iStandElemIndex);
            }
        }
    }

    if (isUnsteady)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            highOrderCellsSol[pMultiGridLevel][iCell].dofQ1.resize(neqn);
            highOrderCellsSol[pMultiGridLevel][iCell].dofQ2.resize(neqn);
            for (int ieqn = 0; ieqn < neqn; ++ieqn)
            {
                highOrderCellsSol[pMultiGridLevel][iCell].dofQ1[ieqn].resize(highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol);
                highOrderCellsSol[pMultiGridLevel][iCell].dofQ2[ieqn].resize(highOrderCellsSol[pMultiGridLevel][iCell].nDOFsSol);
            }
        }
    }

    if (viscousType == LAMINAR)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            highOrderCellsSol[pMultiGridLevel][iCell].visl.resize(highOrderCellsSol[pMultiGridLevel][iCell].integNumberOfGaussPoint);
        }
    }
    else if (viscousType > LAMINAR)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            highOrderCellsSol[pMultiGridLevel][iCell].visl.resize(highOrderCellsSol[pMultiGridLevel][iCell].integNumberOfGaussPoint);
            highOrderCellsSol[pMultiGridLevel][iCell].vist.resize(highOrderCellsSol[pMultiGridLevel][iCell].integNumberOfGaussPoint);
        }
    }
}

//! @Initialize HighOrderFaceSol for different p_multigrid
//! @param[in]   grid_in                   the input unstructgrid
//! @param[in]   dgSolOrder                the accuracy of DG solution
//! @param[in]   pMultiGrid                pmultigrid number
//! @param[out]                            initialize HighOrderfaceSol and build HighOrderStandardface
void HighOrderGrid::InitHighOrderFaceSol(UnstructGrid *grid, int pMultiGrid, int viscousType)
{     
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    const RDouble *pointX = grid->GetX();
    const RDouble *pointY = grid->GetY();
    const RDouble *pointZ = grid->GetZ();
    const RDouble *Xcc = grid->GetCellCenterX();
    const RDouble *Ycc = grid->GetCellCenterY();
    const RDouble *Zcc = grid->GetCellCenterZ();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();


    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    unsigned long positionOfDof = 0;                    //! @the dof position of curret cell,using for 
    unsigned short int integNumberOfGauss = 0; 
    unsigned short int accuracyOfIntigration = 0;       //! @the face integation accuracy
    unsigned short int iFaceVTKType = 0;
    unsigned short int iStandElemIndex = 0;             //! @the type of face
    unsigned short int nFaceNodes = 0;
    unsigned short int leftCell = 0, rightCell = 0;
    unsigned short int dofOfLeftCell = 0, dofOfRightCell = 0;
    unsigned short int nPloyOfLeftCell = 0, nPloyOfRightCell = 0;
    RDouble lenScaleOfLeftCell = 0.0, lenScaleOfRightCell = 0.0;
    RDouble *leftOrthMatrix, *rightOrthMatrix; 
    positionOfDof = 0;
    int nDim = grid->GetDim();
    int neqn = 5;
    unsigned long pTotalGaussPoint = 0;

    if (nDim == 3)
    {
        neqn = 5;
    }

    highOrderFacesSol[pMultiGrid].resize(nTotalFace);

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == INTERFACE)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                highOrderFacesSol[pMultiGrid][iFace].globalID = iFace;

                leftCell  = static_cast<unsigned short>(leftCellOfFace[iFace]);
                rightCell = static_cast<unsigned short>(rightCellOfFace[iFace]);
                highOrderFacesSol[pMultiGrid][iFace].offsetGaussPoint = pTotalGaussPoint;

                nPloyOfLeftCell  = highOrderCellsSol[pMultiGrid][leftCell].nPolySol;
                nPloyOfRightCell = highOrderCellsSol[pMultiGrid][rightCell].nPolySol;

                highOrderFacesSol[pMultiGrid][iFace].nPolySol = max(nPloyOfLeftCell, nPloyOfRightCell);    //! polynomial order
                //! the face integation accuracy
                highOrderFacesSol[pMultiGrid][iFace].accuracyOfFaceIntigration = 2 * highOrderFacesSol[pMultiGrid][iFace].nPolySol + 2;
                accuracyOfIntigration = highOrderFacesSol[pMultiGrid][iFace].accuracyOfFaceIntigration;

                //! unsigned short iFaceVTKType 
                iFaceVTKType = highOrderFaces[iFace].VTK_Type;

                //! get nodes of face
                nFaceNodes = static_cast<unsigned short>(highOrderFaces[iFace].nNodes);

                //! Determine the VTK type
                for (iStandElemIndex = 0; iStandElemIndex < standardFaceElements.size(); ++iStandElemIndex) 
                {
                    if (standardFaceElements[iStandElemIndex].SameStandardElement(iFaceVTKType, highOrderFaces[iFace].nPolyGrid, accuracyOfIntigration))
                    {
                        break;
                    }
                }

                if (iStandElemIndex == standardFaceElements.size())
                    standardFaceElements.push_back(HighOrderStandardFace (iFaceVTKType, nFaceNodes, accuracyOfIntigration));

                highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint = static_cast<unsigned short>(standardFaceElements[iStandElemIndex].GetNumberOfIntegrationPoint());
                integNumberOfGauss = highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint;
                pTotalGaussPoint += integNumberOfGauss;
                //! apply memory to save coord of integration
                highOrderFacesSol[pMultiGrid][iFace].coorXIntegrationPoints.resize(integNumberOfGauss);
                highOrderFacesSol[pMultiGrid][iFace].coorYIntegrationPoints.resize(integNumberOfGauss);
                highOrderFacesSol[pMultiGrid][iFace].coorZIntegrationPoints.resize(integNumberOfGauss);
                const unsigned long * nodeID = highOrderFaces[iFace].GetnodeIDsGrid();

                //! 1:compute Coor of face integration
                highOrderFacesSol[pMultiGrid][iFace].ComputeCoorOfIntegrationPoint(standardFaceElements[iStandElemIndex].GetNnode(), nodeID, pointX, pointY, pointZ, 
                    standardFaceElements[iStandElemIndex].GetShapeFuncFaceIntegration());

                //! compute basis function of face integration
                //! 2:compute orthgonal basis function    
                dofOfLeftCell = highOrderCellsSol[pMultiGrid][leftCell].nDOFsSol;
                dofOfRightCell = highOrderCellsSol[pMultiGrid][rightCell].nDOFsSol;
                lenScaleOfLeftCell  = highOrderCells[leftCell].lenScale;
                lenScaleOfRightCell = highOrderCells[rightCell].lenScale;
                leftOrthMatrix = highOrderCellsSol[pMultiGrid][leftCell].GetOrthogonalizedMatrix();
                rightOrthMatrix = highOrderCellsSol[pMultiGrid][rightCell].GetOrthogonalizedMatrix();

                highOrderFacesSol[pMultiGrid][iFace].leftBasisFuncFaceIntegration.resize(integNumberOfGauss);
                highOrderFacesSol[pMultiGrid][iFace].rightBasisFuncFaceIntegration.resize(integNumberOfGauss);
                for (int iGaussPoint = 0; iGaussPoint < integNumberOfGauss; ++iGaussPoint)
                {                
                    highOrderFacesSol[pMultiGrid][iFace].leftBasisFuncFaceIntegration[iGaussPoint].resize(dofOfLeftCell); 
                    highOrderFacesSol[pMultiGrid][iFace].rightBasisFuncFaceIntegration[iGaussPoint].resize(dofOfRightCell);

                    highOrderFacesSol[pMultiGrid][iFace].ComputeBasisFunction(iGaussPoint, leftCell, rightCell, nTotalCell, nPloyOfLeftCell, nPloyOfRightCell, 
                        lenScaleOfLeftCell, lenScaleOfRightCell, Xcc, Ycc, Zcc, leftOrthMatrix, rightOrthMatrix);
                }

                //! 3:compute the normals in the integration points of the face
                highOrderFacesSol[pMultiGrid][iFace].metricNormalsFace.resize(integNumberOfGauss);
                highOrderFacesSol[pMultiGrid][iFace].JacDetFaceIntegration.resize(integNumberOfGauss);
                //! 3-1: compute dx/dxi,dx/det,   dy/dxi,dy/det,   dz/dxi,dz/det
                RDouble * gradXiShapeFunction;
                RDouble * gradEtShapeFunction;

                for (int iGaussPoint = 0; iGaussPoint < integNumberOfGauss; ++ iGaussPoint)
                {
                    gradXiShapeFunction = standardFaceElements[iStandElemIndex].GetdXiShapeFunctionAtIntegration(iGaussPoint);
                    gradEtShapeFunction = standardFaceElements[iStandElemIndex].GetdEtShapeFunctionAtIntegration(iGaussPoint);

                    //! 0:nx,1:ny,2:nz
                    highOrderFacesSol[pMultiGrid][iFace].metricNormalsFace[iGaussPoint].resize(3);

                    highOrderFacesSol[pMultiGrid][iFace].ComputeMetricNormalsFace(iGaussPoint, iFaceVTKType, standardFaceElements[iStandElemIndex].GetNnode(), nodeID, pointX, pointY, pointZ,
                        gradXiShapeFunction, gradEtShapeFunction);
                }

            }
        }
        else if (bcType != INTERFACE)    //! need modify to fit parallel
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;

                highOrderFacesSol[pMultiGrid][iFace].globalID = iFace;

                leftCell  = static_cast<unsigned short>(leftCellOfFace[iFace]);
                rightCell = static_cast<unsigned short>(rightCellOfFace[iFace]);
                highOrderFacesSol[pMultiGrid][iFace].offsetGaussPoint = pTotalGaussPoint;


                nPloyOfLeftCell  = highOrderCellsSol[pMultiGrid][leftCell].nPolySol;
                nPloyOfRightCell = nPloyOfLeftCell;
                highOrderFacesSol[pMultiGrid][iFace].nPolySol = nPloyOfLeftCell;    //!polynomial order
                //! the face integation accuracy
                highOrderFacesSol[pMultiGrid][iFace].accuracyOfFaceIntigration = 2 * highOrderFacesSol[pMultiGrid][iFace].nPolySol + 2;
                accuracyOfIntigration = highOrderFacesSol[pMultiGrid][iFace].accuracyOfFaceIntigration;

                //! unsigned short iFaceVTKType
                iFaceVTKType = highOrderFaces[iFace].VTK_Type;

                //! get nodes of face
                nFaceNodes = static_cast<unsigned short>(highOrderFaces[iFace].nNodes);

                //! Determine the VTK type
                for (iStandElemIndex = 0; iStandElemIndex < standardFaceElements.size(); ++ iStandElemIndex)
                {
                    if (standardFaceElements[iStandElemIndex].SameStandardElement(iFaceVTKType, highOrderFaces[iFace].nPolyGrid, accuracyOfIntigration))
                    {
                        break;
                    }
                }

                if (iStandElemIndex == standardFaceElements.size())
                    standardFaceElements.push_back(HighOrderStandardFace (iFaceVTKType, nFaceNodes, accuracyOfIntigration));

                highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint = static_cast<unsigned short>(standardFaceElements[iStandElemIndex].GetNumberOfIntegrationPoint());
                integNumberOfGauss = highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint;
                pTotalGaussPoint += integNumberOfGauss;
                //! apply memory to save coord of integration
                highOrderFacesSol[pMultiGrid][iFace].coorXIntegrationPoints.resize(integNumberOfGauss);
                highOrderFacesSol[pMultiGrid][iFace].coorYIntegrationPoints.resize(integNumberOfGauss);
                highOrderFacesSol[pMultiGrid][iFace].coorZIntegrationPoints.resize(integNumberOfGauss);
                const unsigned long *nodeID = highOrderFaces[iFace].GetnodeIDsGrid();

                //! 1:compute Coor of face integration
                highOrderFacesSol[pMultiGrid][iFace].ComputeCoorOfIntegrationPoint(standardFaceElements[iStandElemIndex].GetNnode(), nodeID, pointX, pointY, pointZ, 
                    standardFaceElements[iStandElemIndex].GetShapeFuncFaceIntegration());

                //! compute basis function of face integration
                //! 2:compute orthgonal basis function    
                dofOfLeftCell = highOrderCellsSol[pMultiGrid][leftCell].nDOFsSol;
                dofOfRightCell = dofOfLeftCell;
                lenScaleOfLeftCell  = highOrderCells[leftCell].lenScale;
                lenScaleOfRightCell = lenScaleOfLeftCell;
                leftOrthMatrix = highOrderCellsSol[pMultiGrid][leftCell].GetOrthogonalizedMatrix();
                rightOrthMatrix = leftOrthMatrix;

                highOrderFacesSol[pMultiGrid][iFace].leftBasisFuncFaceIntegration.resize(integNumberOfGauss);
                highOrderFacesSol[pMultiGrid][iFace].rightBasisFuncFaceIntegration.resize(integNumberOfGauss);

                for (int iGaussPoint = 0; iGaussPoint < integNumberOfGauss; ++ iGaussPoint)
                {
                    highOrderFacesSol[pMultiGrid][iFace].leftBasisFuncFaceIntegration[iGaussPoint].resize(dofOfLeftCell);
                    highOrderFacesSol[pMultiGrid][iFace].rightBasisFuncFaceIntegration[iGaussPoint].resize(dofOfRightCell);
                    highOrderFacesSol[pMultiGrid][iFace].ComputeBasisFunction(iGaussPoint, leftCell, rightCell, nTotalCell, nPloyOfLeftCell, nPloyOfRightCell, 
                        lenScaleOfLeftCell, lenScaleOfRightCell, Xcc, Ycc, Zcc, leftOrthMatrix, rightOrthMatrix);
                }

                //! 3:compute the normals in the integration points of the face
                highOrderFacesSol[pMultiGrid][iFace].metricNormalsFace.resize(integNumberOfGauss);
                highOrderFacesSol[pMultiGrid][iFace].JacDetFaceIntegration.resize(integNumberOfGauss);

                //! 3-1: compute dx/dxi,dx/det, dy/dxi,dy/det, dz/dxi,dz/det
                RDouble * gradXiShapeFunction;
                RDouble * gradEtShapeFunction;

                for (int iGaussPoint = 0; iGaussPoint < integNumberOfGauss; ++ iGaussPoint)
                {
                    gradXiShapeFunction = standardFaceElements[iStandElemIndex].GetdXiShapeFunctionAtIntegration(iGaussPoint);
                    gradEtShapeFunction = standardFaceElements[iStandElemIndex].GetdEtShapeFunctionAtIntegration(iGaussPoint);

                    //! 0:nx,1:ny,2:nz
                    highOrderFacesSol[pMultiGrid][iFace].metricNormalsFace[iGaussPoint].resize(3);

                    highOrderFacesSol[pMultiGrid][iFace].ComputeMetricNormalsFace(iGaussPoint, iFaceVTKType, standardFaceElements[iStandElemIndex].GetNnode(), 
                        nodeID, pointX, pointY, pointZ, gradXiShapeFunction, gradEtShapeFunction);
                }

            }

        }
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        highOrderFacesSol[pMultiGrid][iFace].globalID = iFace;

        leftCell  = static_cast<unsigned short>(leftCellOfFace[iFace]);
        rightCell = static_cast<unsigned short>(rightCellOfFace[iFace]);

        nPloyOfLeftCell  = highOrderCellsSol[pMultiGrid][leftCell].nPolySol;
        nPloyOfRightCell = highOrderCellsSol[pMultiGrid][rightCell].nPolySol;

        highOrderFacesSol[pMultiGrid][iFace].nPolySol = max(nPloyOfLeftCell, nPloyOfRightCell);   //polynomial order
        //! the face integation accuracy
        highOrderFacesSol[pMultiGrid][iFace].accuracyOfFaceIntigration = 2 * highOrderFacesSol[pMultiGrid][iFace].nPolySol + 2;
        accuracyOfIntigration = highOrderFacesSol[pMultiGrid][iFace].accuracyOfFaceIntigration;

        //! unsigned short iFaceVTKType
        iFaceVTKType = highOrderFaces[iFace].VTK_Type;

        //! get nodes of face
        nFaceNodes = static_cast<unsigned short>(highOrderFaces[iFace].nNodes);

        //! Determine the VTK type
        for (iStandElemIndex = 0; iStandElemIndex < standardFaceElements.size(); ++iStandElemIndex)
        {
            if (standardFaceElements[iStandElemIndex].SameStandardElement(iFaceVTKType, highOrderFaces[iFace].nPolyGrid, accuracyOfIntigration))
            {
                break;
            }
        }

        if (iStandElemIndex == standardFaceElements.size())
            standardFaceElements.push_back(HighOrderStandardFace (iFaceVTKType, nFaceNodes, accuracyOfIntigration));

        highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint = static_cast<unsigned short>(standardFaceElements[iStandElemIndex].GetNumberOfIntegrationPoint());
        integNumberOfGauss = highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint;
        highOrderFacesSol[pMultiGrid][iFace].offsetGaussPoint = pTotalGaussPoint;
        pTotalGaussPoint += integNumberOfGauss;
        //! apply memory to save coord of integration
        highOrderFacesSol[pMultiGrid][iFace].coorXIntegrationPoints.resize(integNumberOfGauss);
        highOrderFacesSol[pMultiGrid][iFace].coorYIntegrationPoints.resize(integNumberOfGauss);
        highOrderFacesSol[pMultiGrid][iFace].coorZIntegrationPoints.resize(integNumberOfGauss);
        const unsigned long *nodeID = highOrderFaces[iFace].GetnodeIDsGrid();

        //! 1:compute Coor of face integration
        highOrderFacesSol[pMultiGrid][iFace].ComputeCoorOfIntegrationPoint(standardFaceElements[iStandElemIndex].GetNnode(), nodeID, pointX, pointY, pointZ, 
            standardFaceElements[iStandElemIndex].GetShapeFuncFaceIntegration());

        //! compute basis function of face integration
        //! 2:compute orthgonal basis function    
        dofOfLeftCell = highOrderCellsSol[pMultiGrid][leftCell].nDOFsSol;
        dofOfRightCell = highOrderCellsSol[pMultiGrid][rightCell].nDOFsSol;
        lenScaleOfLeftCell  = highOrderCells[leftCell].lenScale;
        lenScaleOfRightCell = highOrderCells[rightCell].lenScale;
        leftOrthMatrix = highOrderCellsSol[pMultiGrid][leftCell].GetOrthogonalizedMatrix();
        rightOrthMatrix = highOrderCellsSol[pMultiGrid][rightCell].GetOrthogonalizedMatrix();

        highOrderFacesSol[pMultiGrid][iFace].leftBasisFuncFaceIntegration.resize(integNumberOfGauss);
        highOrderFacesSol[pMultiGrid][iFace].rightBasisFuncFaceIntegration.resize(integNumberOfGauss);
        for (int iGaussPoint = 0; iGaussPoint < integNumberOfGauss; ++ iGaussPoint)
        {
            highOrderFacesSol[pMultiGrid][iFace].leftBasisFuncFaceIntegration[iGaussPoint].resize(dofOfLeftCell);
            highOrderFacesSol[pMultiGrid][iFace].rightBasisFuncFaceIntegration[iGaussPoint].resize(dofOfRightCell);
            highOrderFacesSol[pMultiGrid][iFace].ComputeBasisFunction(iGaussPoint, leftCell, rightCell, nTotalCell, nPloyOfLeftCell, nPloyOfRightCell, 
                lenScaleOfLeftCell, lenScaleOfRightCell, Xcc, Ycc, Zcc, leftOrthMatrix, rightOrthMatrix);
        }

        //! 3:compute the normals in the integration points of the face
        highOrderFacesSol[pMultiGrid][iFace].metricNormalsFace.resize(integNumberOfGauss);
        highOrderFacesSol[pMultiGrid][iFace].JacDetFaceIntegration.resize(integNumberOfGauss);
        //! 3-1: compute dx/dxi,dx/det,   dy/dxi,dy/det,   dz/dxi,dz/det
        RDouble *gradXiShapeFunction;
        RDouble *gradEtShapeFunction;

        for (int iGaussPoint = 0; iGaussPoint < integNumberOfGauss; ++ iGaussPoint)
        { 
            gradXiShapeFunction = standardFaceElements[iStandElemIndex].GetdXiShapeFunctionAtIntegration(iGaussPoint);
            gradEtShapeFunction = standardFaceElements[iStandElemIndex].GetdEtShapeFunctionAtIntegration(iGaussPoint);

            //! 0:nx,1:ny,2:nz
            highOrderFacesSol[pMultiGrid][iFace].metricNormalsFace[iGaussPoint].resize(3);

            highOrderFacesSol[pMultiGrid][iFace].ComputeMetricNormalsFace(iGaussPoint, iFaceVTKType, standardFaceElements[iStandElemIndex].GetNnode(), nodeID, pointX, pointY, pointZ,
                gradXiShapeFunction, gradEtShapeFunction);
        }
    }

    totalFaceGaussPoint[pMultiGrid] = pTotalGaussPoint;

    //! HighOrderFaceµÄstandardFaceElements type
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        //! Determine the standardElementType, pMultiGrid
        iFaceVTKType = highOrderFaces[iFace].VTK_Type;
        int nPolyGrid = highOrderFaces[iFace].nPolyGrid;
        for (iStandElemIndex = 0; iStandElemIndex < standardFaceElements.size(); ++iStandElemIndex)
        {
            if (iFaceVTKType == standardFaceElements[iStandElemIndex].GetVTK_Type() && 
                nPolyGrid   == standardFaceElements[iStandElemIndex].GetNpoly())
            {
                highOrderFaces[iFace].standardElementIndex = static_cast<unsigned short>(iStandElemIndex);
            }
        } 
        if (pMultiGrid == 0)
        {
            RDouble areaTmp = 0.0, xfc = 0.0, yfc = 0.0, zfc = 0.0;
            const RDouble *weightCoeff;
            int iCurStandElemIndex = highOrderFaces[iFace].standardElementIndex;
            weightCoeff = standardFaceElements[iCurStandElemIndex].GetWtIntegration();
            for (int iGaussPoint = 0; iGaussPoint < highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint; ++ iGaussPoint)
            {
                areaTmp += weightCoeff[iGaussPoint] * highOrderFacesSol[pMultiGrid][iFace].JacDetFaceIntegration[iGaussPoint];
                xfc += weightCoeff[iGaussPoint] * highOrderFacesSol[pMultiGrid][iFace].coorXIntegrationPoints[iGaussPoint] * highOrderFacesSol[pMultiGrid][iFace].JacDetFaceIntegration[iGaussPoint]; 
                yfc += weightCoeff[iGaussPoint] * highOrderFacesSol[pMultiGrid][iFace].coorYIntegrationPoints[iGaussPoint] * highOrderFacesSol[pMultiGrid][iFace].JacDetFaceIntegration[iGaussPoint];
                zfc += weightCoeff[iGaussPoint] * highOrderFacesSol[pMultiGrid][iFace].coorZIntegrationPoints[iGaussPoint] * highOrderFacesSol[pMultiGrid][iFace].JacDetFaceIntegration[iGaussPoint];
            }         
            highOrderFaces[iFace].area = areaTmp;
            highOrderFaces[iFace].centerCoord[0] = xfc / areaTmp;
            highOrderFaces[iFace].centerCoord[1] = yfc / areaTmp;
            highOrderFaces[iFace].centerCoord[2] = zfc / areaTmp;
        }
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        highOrderFacesSol[pMultiGrid][iFace].q.resize(neqn);
        highOrderFacesSol[pMultiGrid][iFace].visl.resize(highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint);

        for (int ieqn = 0; ieqn < neqn; ++ieqn)
        {
            //[neqn][2*gausspoint], [neqn][0]:L,[neqn][0+nTotalFace]:R
            highOrderFacesSol[pMultiGrid][iFace].q[ieqn].resize(2 * highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint);
        }
    }

    if (viscousType == LAMINAR)
    {
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            highOrderFacesSol[pMultiGrid][iFace].visl.resize(2 * highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint);
        }
    }
    else if (viscousType > LAMINAR)
    {
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            highOrderFacesSol[pMultiGrid][iFace].visl.resize(2 * highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint);
            highOrderFacesSol[pMultiGrid][iFace].vist.resize(2 * highOrderFacesSol[pMultiGrid][iFace].integNumberOfGaussPoint);
        }
    }
}

void HighOrderGrid::InitDofOfCell(UnstructGrid * grid, vector < RDouble > & conserValue, int neqn)
{
    int nTotalCell = grid->GetNTotalCell();
    const RDouble *weightCoeff;
    RDouble dofQTmp[35], weigthMultJac = 0.0; //35------>P4
    RDouble volumTmp = 0.0;

    for (int ieqn = 0; ieqn < neqn; ++ ieqn)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {

            for (int iDof = 0; iDof < highOrderCellsSol[0][iCell].nDOFsSol; ++ iDof)
            {
                dofQTmp[iDof] = 0.0;
            }
            volumTmp = 0.0;
            weightCoeff = standardCellElements[highOrderCells[iCell].standardElementIndex].GetWtIntegration();

            for (int iGaussPoint = 0; iGaussPoint < highOrderCellsSol[0][iCell].integNumberOfGaussPoint; ++ iGaussPoint)
            {
                weigthMultJac = weightCoeff[iGaussPoint] * highOrderCellsSol[0][iCell].JacDetCellIntegration[iGaussPoint];

                for (int iDof = 0; iDof < highOrderCellsSol[0][iCell].nDOFsSol; ++ iDof)
                {
                    dofQTmp[iDof] += weigthMultJac * highOrderCellsSol[0][iCell].basisFuncCellIntegration[iGaussPoint][iDof];
                }
                volumTmp += weigthMultJac;
            }
            for (int iDof = 0; iDof < highOrderCellsSol[0][iCell].nDOFsSol; ++ iDof)
            {
                highOrderCellsSol[0][iCell].dofQ[ieqn][iDof] = dofQTmp[iDof] * conserValue[ieqn] / volumTmp;
            }
        }
    }

    //init res of cell
    for (int ieqn = 0; ieqn < neqn; ++ ieqn)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            for (int iDof = 0;iDof < highOrderCellsSol[0][iCell].nDOFsSol; ++iDof)
            {
                highOrderCellsSol[0][iCell].res[ieqn][iDof] = 0.0;
            }
        }
    }

    //init gauss point value of cell
    for (int ieqn = 0; ieqn < neqn; ++ ieqn)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            for (int iGaussPoint = 0;iGaussPoint < highOrderCellsSol[0][iCell].integNumberOfGaussPoint; ++iGaussPoint)
            {
                highOrderCellsSol[0][iCell].q[ieqn][iGaussPoint] = conserValue[ieqn];
            }
        }
    }
}

void HighOrderGrid::InitDofOfCellUnsteady(UnstructGrid * grid, vector < RDouble > & conserValue, int neqn)
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    for (int ieqn = 0; ieqn < neqn; ++ ieqn)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            highOrderCellsSol[0][iCell].dofQ1[ieqn][0] = conserValue[ieqn];
            highOrderCellsSol[0][iCell].dofQ2[ieqn][0] = conserValue[ieqn];

            for (int iDof = 1;iDof < highOrderCellsSol[0][iCell].nDOFsSol; ++ iDof)
            {
                highOrderCellsSol[0][iCell].dofQ1[ieqn][iDof] = 0.0;
                highOrderCellsSol[0][iCell].dofQ2[ieqn][iDof] = 0.0;
            }
        }
    }
}

void HighOrderGrid::InitValueOfFace(const UnstructGrid * grid, const vector < RDouble > & conserValue, int neqn)
{
    //int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    //int nBoundFace = grid->GetNBoundFace();
    //int nTotal = nTotalCell + nBoundFace;

    //init gauss point value of cell
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        HighOrderFaceSol & faceSol = highOrderFacesSol[0][iFace];

        for (int iEqn = 0; iEqn < neqn; ++ iEqn)
        {
            for (int iGaussPoint = 0; iGaussPoint < (2 * faceSol.integNumberOfGaussPoint); ++ iGaussPoint)
            {
                faceSol.q[iEqn][iGaussPoint] = conserValue[iEqn];
            }
        }
    }
}

void HighOrderGrid::InitValueOfFace(const UnstructGrid * grid, const RDouble * conserValue, int neqn)
{
    //int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    //int nBoundFace = grid->GetNBoundFace();
    //int nTotal = nTotalCell + nBoundFace;

    //init gauss point value of cell
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        HighOrderFaceSol & faceSol = highOrderFacesSol[0][iFace];

        for (int iEqn = 0; iEqn < neqn; ++ iEqn)
        {
            for (int iGaussPoint = 0; iGaussPoint < (2 * faceSol.integNumberOfGaussPoint); ++ iGaussPoint)
            {
                faceSol.q[iEqn][iGaussPoint] = conserValue[iEqn];
            }
        }
    }
}

void HighOrderGrid::GetNumberOfDOFs(int nCells, int * dofsOfCells)
{
    for (int iCell = 0; iCell < nCells; ++ iCell)
    {
        dofsOfCells[ iCell ] = highOrderCellsSol[0][iCell].nDOFsSol;
    }
}

void HighOrderGrid::GetQDOFs(const UnstructGrid * grid, int neqn, RDouble *** q)
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        for (int iEqn = 0; iEqn < neqn; ++ iEqn)
        {
            for (int iDof = 0;iDof < highOrderCellsSol[0][iCell].nDOFsSol; ++iDof)
            {
                q[iCell][iEqn][iDof] = highOrderCellsSol[0][iCell].dofQ[iEqn][iDof];
            }
        }
    }
}

}