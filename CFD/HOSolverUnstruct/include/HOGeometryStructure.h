//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      HOGeometryStructure.h
//! @brief     High Order Geometry_structure
//! @author    Li Ming, Gong Xiaoquan, Zhang Jian, Wan yunbo, Xu gang.

#pragma once

#include "HODefine.h"
#include "Geo_UnstructGrid.h"
#include "HOStandardElement.h"

#include <vector>
#pragma warning(disable:4100)
namespace HOUnstruct
{

//class HighOrderCell;
//class HighOrderFace;
//class HO_Point;

/*!
 * \class HO_Point
 * \brief Class to a point for the unstr high-order solver.
 * \author 
 * \version 
 */
class HighOrderNode
{
public:
    unsigned long ID;    /*!< \brief The global ID of this point in the grid. */

    RDouble coor[3];         /*!< \brief Array with the coordinates of the node. */

    vector<unsigned long> cellIDsOfPoint; /*!< \brief Vector with the cell IDs of this point. */

    /*!
    * \brief Default constructor of the class. Initialize the coordinates to zero
    to avoid a valgrind warning in two space dimensions.
    */
    HighOrderNode(void) {}

    /*!
    * \brief Destructor of the class. Nothing to be done.
    */
    ~HighOrderNode(void) {}

    /*!
    * \brief Copy constructor of the class.
    */
    HighOrderNode(const HighOrderNode &other) {}

    /*!
    * \brief Assignment operator of the class.
    */
    HighOrderNode& operator=(const HighOrderNode &other) { return * this; }

    /*!
    * \brief Less than operator of the class. Needed for the sorting.
    */
    bool operator<(const HighOrderNode &other) const;

    /*!
    * \brief Equal operator of the class. Needed for the removal of double entities.
    */
    bool operator==(const HighOrderNode &other) const;

private:
    /*!
    * \brief Copy function. Needed for the copy constructor and assignment operator.
    */
    void Copy(const HighOrderNode &other);
};

/*!
 * \class HighOrderFace
 * \brief Class to store a cell for the unstr high-order solver.
 * \only geometry information,using for h-multigrid
 * \author 
 * \version
 */
class HighOrderFace
{
public:
    unsigned long ID;    /*!< \brief The global ID of this face in the grid. */
    unsigned int nNodes;
    vector<unsigned long> nodeIDs; /*!< \brief Vector with the node IDs of the grid for this element.
                                                 In this vector the original sequence of the grid file is stored. */
    unsigned short nPolyGrid;    /*!< \brief Polynomial degree for the geometry of the element. */
    unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
                                            /*!< different VTK_Type,  VTK_Type only mark the type, no linear or quadratic type */
    unsigned short nDOFsGrid;    /*!< \brief Number of DOFs for the geometry of the element. */

    //unsigned short indStandardElement; /*!< \brief Index in the vector of standard face elements. */

    unsigned long cellID0;              /*!< \brief cell ID adjacent to side 0 of the face. */
    unsigned long cellID1;              /*!< \brief cell ID adjacent to side 1 of the face. */

    unsigned short standardElementIndex;     /*!< \brief mark standard Element type , standardFaceElements[standardElementType]

    //unsigned short localIDInCell0;   /*!< \brief face ID in adjacent cell 0. */
    //unsigned short localIDInCell1;   /*!< \brief face ID in adjacent cell 1. */

    RDouble area;
    RDouble centerCoord[3];

    /*!
    * \brief Constructor of the class. Initialize some pointers to NULL.
    */
    HighOrderFace(void) {}

    /*!
    * \brief Destructor of the class. Nothing to be done.
    */
    ~HighOrderFace(void) {}

    /*!
    * \brief Copy constructor of the class.
    */
    HighOrderFace(const HighOrderFace &other) {}

    /*!
    * \brief Assignment operator of the class.
    */
    HighOrderFace& operator=(const HighOrderFace &other) { return * this; }

    /*!
    * \brief Less than operator of the class. Needed for the sorting.
    The criterion for comparison are the standard element and
    adjacent volume ID's.
    */
    bool operator<(const HighOrderFace &other) const;
    const unsigned long * GetnodeIDsGrid(void) const { return nodeIDs.data();}

private:
    /*!
    * \brief Copy function. Needed for the copy constructor and assignment operator.
    */
    void Copy(const HighOrderFace &other);
};

/*!
 * \class HighOrderFaceSol
 * \brief Class to store a Face for the unstr high-order solver.
 * \only solution information, using for Hp-multigrid
 * \author 
 * \version
 */
class HighOrderFaceSol
{
public:
    unsigned long globalID;                 /*!< \brief The global ID of this face in the grid. */
    unsigned short nPolySol;                /*!< \brief Polynomial degree for this face. */
    //unsigned short nDOFsSol;                /*!< \brief Number of DOFs for the solution of the element of this face. */
    unsigned short integNumberOfGaussPoint;
    unsigned long offsetGaussPoint;         //! the location of gausspoint in all face gausspoint
    unsigned short accuracyOfFaceIntigration;  //! @the max accuracy of face for left and right gauss numerical integration
    vector< vector< RDouble > > metricNormalsFace;     /*!< \brief The normals in the integration points of the face.
                                               The normals point from side 0 to side 1. 0-2:nx,ny,nz,[iGausspoint][normal]*/
                   
    vector< RDouble > JacDetFaceIntegration;   /*!< \brief The Jacobian determinant for the integration points of this face. */
    vector< RDouble > coorXIntegrationPoints;  /*!< \brief Coordinates for the integration points of this face. */
    vector< RDouble > coorYIntegrationPoints;  /*!< \brief Coordinates for the integration points of this face. */
    vector< RDouble > coorZIntegrationPoints;  /*!< \brief Coordinates for the integration points of this face. */
    vector< vector< RDouble > > leftBasisFuncFaceIntegration; /*!< \brief the cell basis functions at the integration points. */
    vector< vector< RDouble > > rightBasisFuncFaceIntegration; /*!< \brief the cell basis functions at the integration points. */

   //vector<RDouble> gridVelocities;         /*!< \brief Grid velocities in the integration points of this face. */
    vector< RDouble > wallDistance;           /*!< \brief The wall distance to the viscous walls for
                                              the integration points of this face. */
    //[nequ][GaussPointIndex]
    vector< vector< RDouble > > q;           /*!< \conservation variable of gauss point,cuurect step n. */
    //[GaussPointIndex]
    vector< RDouble > visl;             /*!< \the visl in face gauss point*/
    vector< RDouble > vist;             /*!< \the vist in face gauss point*/
    /*!
    * \brief Constructor of the class. Initialize some pointers to NULL.
    */
    HighOrderFaceSol(void) {}

    /*!
    * \brief Destructor of the class. Nothing to be done.
    */
    ~HighOrderFaceSol(void) {}

    /*!
    * \brief Copy constructor of the class.
    */
    HighOrderFaceSol(const HighOrderFaceSol &other) {}

    /*!
    * \brief Assignment operator of the class.
    */
    HighOrderFaceSol& operator=(const HighOrderFaceSol &other) { return * this; }

    /*!
    * \brief Less than operator of the class. Needed for the sorting.
    The criterion for comparison are the standard element and
    adjacent volume ID's.
    */
    bool operator<(const HighOrderFaceSol &other) const;

    void ComputeCoorOfIntegrationPoint(int nodeNumber, const unsigned long *nodeID, const RDouble *pointX, const RDouble *pointY, const RDouble *pointZ,
         const vector < vector < RDouble > > & shapeFunctionFacePoint);
    
    void ComputeBasisFunction(int iGaussPoint, int leftCell, int rightCell, int nTotalCell, int nPloyOfLeftCell, int nPloyOfRightCell, 
         RDouble lenScaleOfLeftCell, RDouble lenScaleOfRightCell, const RDouble *Xcc, const RDouble *Ycc, const RDouble *Zcc, const RDouble *leftOrthMatrix, const RDouble *rightOrthMatrix);

    void ComputeMetricNormalsFace(int iGaussPoint, int iFaceVTKType, int nodenumber, const unsigned long *nodeID, const RDouble *pointX, const RDouble *pointY, const RDouble *pointZ,
                                                        RDouble *gradXiShapeFunction, RDouble *gradEtShapeFunction);

private:
    /*!
    * \brief Copy function. Needed for the copy constructor and assignment operator.
    */
    void Copy(const HighOrderFaceSol &other);
};

/*!
 * \class HighOrderCell
 * \brief Class to store a cell for the unstr high-order solver.
 * \only geometry information,using for h-multigrid
 * \author 
 * \version
 */
class HighOrderCell 
{
public:
    unsigned long ID;             /*!< \brief Global element ID of this element. */
    
    unsigned short nFaces;                  /*!< \brief Number of faces of the element. */
    vector<unsigned long> faceIDs;      /*!< \brief Vector with the node IDs of the grid for this element. */
    unsigned short nNodes;
    vector<unsigned long> nodeIDs;      /*!< \brief Vector with the node IDs of the grid for this element. */

    //unsigned short indStandardElement;      /*!< \brief Index in the vector of standard elements. */
    //bool IsOwnedCell;                       /*!< \brief Whether or not this is an owned element. */
    //bool JacIsConsideredConstant;           /*!< \brief Whether or not the Jacobian of the transformation to the standard element is considered constant. */
    
    unsigned short nPolyGrid;               /*!< \brief Polynomial degree for the geometry of the element. */
    unsigned short VTK_Type;                /*!< \brief Element type using the VTK convention. */
    unsigned short nDOFsGrid;               /*!< \brief Number of DOFs for the geometry of the element. */
   
    unsigned short standardElementIndex;     /*!< \brief mark standard Element type , standardCellElements[standardElementType]
                                            /*!< different VTK_Type,  VTK_Type only mark the type, no linear or quadratic type */
    RDouble volume;
    RDouble centerCoord[3];
    RDouble lenScale;                /*!< \brief Length scale of the element. */

    //RDouble shockSensorValue;               /*!< \brief Value for sensing a shock */
    //RDouble shockArtificialViscosity;       /*!< \brief Artificial viscosity for a shock */
       
    /*!
    * \brief Constructor of the class. Initialize the pointers to NULL.
    */
    HighOrderCell(void) {}

    /*!
    * \brief Destructor of the class. Nothing to be done.
    */
    ~HighOrderCell(void) {}

    /*!
    * \brief Get all the corner points of all the faces of this element. It must be made sure
    that the numbering of the faces is identical to the numbering used for the
    standard elements.
    * \param[out] nFaces         - Number of faces of this element.
    * \param[out] nPointsPerFace - Number of corner points for each of the faces.
    * \param[out] faceConn       - Global IDs of the corner points of the faces.
    */
    //void GetCornerPointsAllFaces(unsigned short &numFaces, unsigned short nPointsPerFace[], unsigned long  faceConn[6][4]);
    const unsigned long * GetnodeIDsGrid(void) const { return nodeIDs.data();}
};

/*!
 * \class HighOrderCellSol
 * \brief Class to store a cell for the unstr high-order solver.
 * \only solution information, using for Hp-multigrid
 * \author 
 * \version
 */

class HighOrderCellSol
{
public:
    //unsigned short VTK_Type;             /*!< \brief Element type using the VTK convention. */
    unsigned long ID;                      /*!< \brief Global element ID of this element. */
    unsigned short nPolySol;               /*!< \brief Polynomial degree for the solution of the element. */
    unsigned short nFVPolySol;             /*!< \brief Polynomial degree for the solution of the element. */
    unsigned short nDOFsSol;               /*!< \brief Number of DOFs for the solution of the element. */
    unsigned short integNumberOfGaussPoint;      
    
    unsigned long offsetDOFsSolGlobal;     /*!< \brief Global offset of the solution DOFs of this element. */
    unsigned long offsetDOFsSolLocal;      /*!< \brief Local offset of the solution DOFs of this element. */

    vector<RDouble> massMatrix;             /*!< \brief Mass matrix for this element. */
    vector<RDouble> invMassMatrix;          /*!< \brief Inverse mass matrix for this element. */
    vector<RDouble> orthogonalizedMatrix;   /*!< \brief orthogonalized matrix for basis function. */

    vector<RDouble> coorXIntegrationPoints; /*!< \brief The coordinates X of the integration points of this element. */
    vector<RDouble> coorYIntegrationPoints; /*!< \brief The coordinates Y of the integration points of this element. */
    vector<RDouble> coorZIntegrationPoints; /*!< \brief The coordinates Z of the integration points of this element. */
    
    //[GaussPointIndex][basisFuncIndex]
    vector< vector<RDouble> > basisFuncCellIntegration;   /*!< \brief the cell basis functions at the integration points. */
    vector< vector<RDouble> > dxBasisFuncCellIntegration; /*!< \brief x-derivatives of the cell basis functions at the integration points. */
    vector< vector<RDouble> > dyBasisFuncCellIntegration; /*!< \brief y-derivatives of the cell basis functions at the integration points. */
    vector< vector<RDouble> > dzBasisFuncCellIntegration; /*!< \brief z-derivatives of the cell basis functions at the integration points. */
    vector<RDouble> JacDetCellIntegration;                /*!< \brief The Jacobian determinant for the integration points of this cell. */

    vector<RDouble> wallDistance;            /*!< \brief The wall distance to the viscous walls for
                                              the integration points of this element. */

    vector<RDouble> reconstructMatrix;       /*!< \brief reconstruction matrix for this element. */
    //[nequ][NDOF]
    vector< vector< RDouble > > dofQ;        /*!< \dof of conservation variable,cuurect step n. */
    vector< vector< RDouble > > dofQ1;       /*!< \dof of conservation variable, step n-1. */
    vector< vector< RDouble > > dofQ2;       /*!< \dof of conservation variable, step n-2. */
    vector< vector< RDouble > > res;         /*!< \res of cuurect step n. */
    //[nequ][GaussPointIndex]
    vector< vector< RDouble > > q;           /*!< \conservation variable of gauss point,cuurect step n. */
    //[nequ]
    vector< RDouble > qAverage;              /*!< \conservation variable of CellSol */
    //[GaussPointIndex]
    vector< RDouble > visl;                  /*!< \the visl  in volume gauss point*/
    vector< RDouble > vist;                  /*!< \the vist  in volume gauss point*/
    RDouble vislAverage;                     /*!< \average of visl */
    RDouble vistAverage;                     /*!< \average of vist */

    typedef Eigen::Matrix<RDouble, Eigen::Dynamic, Eigen::Dynamic> HighOrderMatrix;

    //! @Compute Coordinate of integration point  
    //! @param[in]   nodenumber                               node number of iCell
    //! @param[in]   nodeID                                   node ID of iCell
    //! @param[in]   pointX                                   the coord_X of grid
    //! @param[in]   pointY                                   the coord_Y of grid
    //! @param[in]   pointZ                                   the coord_Z of grid
    //! @param[in]   shapeFunctionVolumPoint                  the shape Function of volum gauss point
    //! @param[out]  return                                   Coordinate of integration point    
    void ComputeCoorOfIntegrationPoint(int nodeNumber, const unsigned long *nodeID, const RDouble *pointX, const RDouble *pointY, const RDouble *pointZ,
         const vector < vector < RDouble > > & shapeFunctionVolumPoint);
    
    void ComputeOriginalBasisFunction(int currentOrder, RDouble lenScaleOfCell, const RDouble Xcc, const RDouble Ycc, const RDouble Zcc);

    void ComputeOriginalMassmatrix(HighOrderMatrix & taylorBasisMassMatrix, const RDouble *weightCoeff);

    void ComputeJacobianDeterminant(int iGaussPoint, int iCellVTKType, int nodenumber, const unsigned long *nodeID, const RDouble *pointX, const RDouble *pointY, const RDouble *pointZ,
        RDouble *gradXiShapeFunction, RDouble *gradEtShapeFunction, RDouble *gradZtShapeFunction);

    void ComputeOrthogonalBasisFunction(int degreeOfFreedom, HighOrderMatrix & orthogonalizationMatrix);
    void TestOrthogonalBasisFunction(const RDouble *weight);

    void ComputeGradOfBasisFunctionIntegrationPoint(HighOrderMatrix & orthogonalizationMatrix, int iGaussPoint, 
        RDouble lenScaleOfCell);

    RDouble * GetOrthogonalizedMatrix()
    {
        return orthogonalizedMatrix.data();
    }
    /*!
    * \brief Constructor of the class. Initialize the pointers to NULL.
    */
    //HighOrderCellSol(void) {}
    HighOrderCellSol(void);

    /*!
    * \brief Destructor of the class. Nothing to be done.
    */
    ~HighOrderCellSol(void) {}
};

//class HighOrderStandardCell;

class HighOrderGrid
{
public:
    //[iCell]
    std::vector<HighOrderCell> highOrderCells;
    std::vector<HighOrderFace> highOrderFaces;
    std::vector<HighOrderNode> highOrderNodes;

    //[pLevel][iCell]
    typedef std::vector< std::vector<HighOrderCellSol> > CellSolType;
    typedef std::vector< std::vector<HighOrderFaceSol> > FaceSolType;
    typedef Eigen::Matrix<RDouble, Eigen::Dynamic, Eigen::Dynamic> HighOrderMatrix;

    //[pLevel][iCell]
    CellSolType highOrderCellsSol;
    FaceSolType highOrderFacesSol;
    //[pLevel]
    std::vector< int > totalFaceGaussPoint;
    
    vector< HighOrderStandardCell > standardCellElements;
    vector< HighOrderStandardFace > standardFaceElements;

public:
    HighOrderGrid(void) {}
    ~HighOrderGrid(void) {}

    void CopyGeometryInfoToHighOrderCell(UnstructGrid * grid);
    void CopyGeometryInfoToHighOrderFace(UnstructGrid * grid);
    void CopyGeometryInfoToHighOrderNode(UnstructGrid * grid);

    //! @Initialize HighOrderCellSol for different p_multigrid 
    //! @param[in]   grid                      the input unstructgrid
    //! @param[in]   dgSolOrder                the accuracy of DG solution
    //! @param[in]   pMultiGrid                pmultigrid number
    //! @param[out]                            initialize HighOrderCellSol and build HighOrderStandardCell 
    void InitHighOrderCellSol(UnstructGrid *grid, int dgSolOrder, int pMultiGrid, int isUnsteady, int viscousType);
    
    //! @Initialize HighOrderFaceSol for different p_multigrid 
    //! @param[in]   grid                      the input unstructgrid
    //! @param[in]   dgSolOrder                the accuracy of DG solution
    //! @param[in]   pMultiGrid                pmultigrid number
    //! @param[out]                            initialize HighOrderFaceSol and build HighOrderStandardFace
    void InitHighOrderFaceSol(UnstructGrid *grid, int pMultiGrid, int viscousType);
    void InitDofOfCell(UnstructGrid *grid, vector < RDouble > & conserValue, int nequ);
    void InitDofOfCellUnsteady(UnstructGrid * grid, vector < RDouble > & conserValue, int nequ);
    void InitValueOfFace(const UnstructGrid * grid, const vector < RDouble > & conserValue, int nequ);
    void InitValueOfFace(const UnstructGrid * grid, const RDouble * conserValue, int nequ);

    void GetNumberOfDOFs(int nCells, int * dofsOfCells);

    //q[iCell][iEqn][iDof]
    void GetQDOFs(const UnstructGrid *grid, int nequ, RDouble *** q);
};

}