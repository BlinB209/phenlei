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
//! @file      HOStandardElement.h
//! @brief     Cell class definition for High Order Solver.
//! @author    Li Ming, Gong Xiaoquan, Zhang Jian, Wan yunbo, Xu gang.

#pragma once

#include "HODefine.h"
#include "Geo_UnstructGrid.h"

#include <vector>

using namespace PHSPACE;
using namespace std;

namespace HOUnstruct
{

class HighOrderStandardElementBase;
class HighOrderStandardCell;
class HighOrderStandardFace;

//! HighOrderStandardElementBase get nIntegrationPoint /location of Integration_Point/ the weights of the integration points
class HighOrderStandardElementBase 
{
protected:
    unsigned short VTK_Type;                 //!  @StandardElement Type
    unsigned short integrationOrder;         //!  @order of integration
    unsigned short nIntegrationPoint;        //!  @integration point's number for this standard element
    bool constJacobian;             //!  @brief Whether or not the element has a constant Jacobian 

    vector< RDouble > xiIntegration; //! @< \brief xi-location of the integration points for this standard element. */
    vector< RDouble > etIntegration; //! @< \brief et-location of the integration points for this standard element, if needed. */
    vector< RDouble > ztIntegration; //! @< \brief zt-location of the integration points for this standard element, if needed. */
    vector< RDouble > wtIntegration; //! @< \brief The weights of the integration points for this standard element. */

public:
    HighOrderStandardElementBase() {}
    virtual ~HighOrderStandardElementBase() {}
    HighOrderStandardElementBase(unsigned short iVTK_Type,unsigned short iIntegrationOrder);  //! import:VTK_Type and integrationOrder，define:xiIntegration,wtIntegration

    int GetVTK_Type(void) const {return VTK_Type; }; //Get StandardElementBase type
    int GetNumberOfIntegrationPoint(void) const {return nIntegrationPoint;}           //Get integration point's number of StandardElementBase
    int GetOrderExact(void){return integrationOrder;}                                 //Get integrationOrder of  StandardElementBase
    const RDouble * GetXiIntegrationCoor(void) const { return xiIntegration.data();}  //Get coord[Xi] Integration_point of StandardElementBase
    const RDouble * GetEtIntegrationCoor(void) const { return etIntegration.data();}  //Get coord[Et] Integration_point of StandardElementBase
    const RDouble * GetZtIntegrationCoor(void) const { return ztIntegration.data();}  //Get coord[Zt] Integration_point of StandardElementBase
    const RDouble * GetWtIntegration(void) const { return wtIntegration.data();}      //Get weights of the integration points for StandardElementBase

    /*!
    * \brief Function, which determines the integration points for a triangle/quadrilateral/tetrahedron/pyramid/prism/hexahedron
    such that polynomials of orderExact are integrated exactly.
    */
    void IntegrationPointsTriangle(void);
    void IntegrationPointsQuadrilateral(void);
    void IntegrationPointsTetrahedron(void);
    void IntegrationPointsPyramid(void);
    void IntegrationPointsPrism(void);
    void IntegrationPointsHexahedron(void);

    /*!
    * \brief Function, which determines the 1D Gauss Legendre integration points and weights.
    * \param[in,out] GLPoints  - The location of the Gauss-Legendre integration points.
    * \param[in,out] GLWeights - The weights of the Gauss-Legendre integration points.
    */
    void GaussLegendrePoints1D(vector<RDouble> & GLPoints, vector<RDouble> & GLWeights);

    /*!
    * \brief Function, which copies the data of the given object into the current object.
    * \param[in] other - Object, whose data is copied.
    */
    void Copy(const HighOrderStandardElementBase & other);

    void ComputeOrthogoMatrix();
};

void TestOrthogoMatrix();

//StandardElement 
class HighOrderStandardCell : public HighOrderStandardElementBase
{
protected:
    //bool constJacobian;           //brief Whether or not the element has a constant Jacobian 

    unsigned short nPoly;        /*!< \brief Polynomial degree of the element. */
    unsigned short nNodes;        /*!< \brief Number of DOFs of the element, the total vertex number. */

    //vector<RDouble> BasisFunction, nNodes defining element;
    vector< RDouble > vertexXi;     /*!< \brief r-location of vertexs for this standard element. */
    vector< RDouble > vertexEt;     /*!< \brief s-location of vertexs for this standard element, if needed. */
    vector< RDouble > vertexZt;     /*!< \brief t-location of vertexs for this standard element, if needed. */

    //[GaussPointsIdx][ShapeDOFsIdx]
    vector< vector< RDouble > > shapeFunctionAtIntegration;       /*!< \brief Lagrangian basis functions in the integration points. */
    vector< vector< RDouble > > dXiShapeFunctionAtIntegration; /*!< \brief xi-derivatives of the Lagrangian basis functions in the integration points. */
    vector< vector< RDouble > > dEtShapeFunctionAtIntegration; /*!< \brief et-derivatives of the Lagrangian basis functions in the integration points. */
    vector< vector< RDouble > > dZtShapeFunctionAtIntegration; /*!< \brief zt-derivatives of the Lagrangian basis functions in the integration points. */

    vector< unsigned short > connFace0; /*!< \brief Local connectivity of face 0 of the element. The numbering of the DOFs is
                                        such that the element is to the left of the face. */
    vector< unsigned short > connFace1; /*!< \brief Local connectivity of face 1 of the element. The numbering of the DOFs is
                                        such that the element is to the left of the face. */
    vector< unsigned short > connFace2; /*!< \brief Local connectivity of face 2 of the element, if present. The numbering
                                        of the DOFs is such that the element is to the left of the face. */
    vector< unsigned short > connFace3; /*!< \brief Local connectivity of face 3 of the element, if present. The numbering
                                        of the DOFs is such that the element is to the left of the face. */
    vector< unsigned short > connFace4; /*!< \brief Local connectivity of face 4 of the element, if present. The numbering
                                        of the DOFs is such that the element is to the left of the face. */
    vector< unsigned short > connFace5; /*!< \brief Local connectivity of face 5 of the element, if present. The numbering
                                        of the DOFs is such that the element is to the left of the face. */    
public:
    HighOrderStandardCell() {}
    virtual ~HighOrderStandardCell() {}
    HighOrderStandardCell(unsigned short VTK_Type, unsigned short nCellNodes, unsigned short integrationOrder); //! import:VTK_Type and integrationOrder，export:integrationOrder、nIntegrationPoint\
    //! position of xi,et,zt

    HighOrderStandardCell(const HighOrderStandardCell &other) : HighOrderStandardElementBase(other) {Copy(other);}

    HighOrderStandardCell& operator= (const HighOrderStandardCell &other) {Copy(other); return (*this);}

    void ShapeFunctionsAtGaussPoints();  //compute shape functions of all points
    void ShapeFunctionsAndDerivativesAtPoint();

    const vector< vector< RDouble > > & GetShapeFunctionsIntegration() { return shapeFunctionAtIntegration; }
    //RDouble * GetShapeFunctionsIntegration(int iGaussPoint);
    RDouble * GetdXiShapeFunctionAtIntegration(int iGaussPoint)
    {
        return dXiShapeFunctionAtIntegration[iGaussPoint].data();
    }

    RDouble * GetdEtShapeFunctionAtIntegration(int iGaussPoint)
    {
        return dEtShapeFunctionAtIntegration[iGaussPoint].data();
    }

    RDouble * GetdZtShapeFunctionAtIntegration(int iGaussPoint)
    {
        return dZtShapeFunctionAtIntegration[iGaussPoint].data();
    }

    /*!
    * \brief Function, which makes available the connectivity of face faceIdx.
    * \return  The pointer to data, which stores the connectivity of face faceIdx.
    */
    unsigned short * GetConnFace(short int faceIdx);

    /*!
    * \brief Function, which makes available the number of DOFs for this standard element.
    * \return  The number of DOFs of this standard element.
    */
    unsigned short GetNnode(void) const { return nNodes; };
    /*!
    * \brief Function, which makes available the polynomial degree for this standard element.
    * \return  The polynomial degree of this standard element.
    */
    unsigned short GetNpoly(void) const { return nPoly; };

    //! @brief Function, which checks if the function arguments correspond to this standard element.
    //! @param[in] val_VTK_Type           - Type of the face using the VTK convention.
    //! @param[in] nCellNodes             - the nodes of cell (linear, quadratic).
    //! @param[in] valIntegrationOrder    - integration order of cell .
    //! @return Whether or not the function arguments correspond to this standard element.
    bool SameStandardElement(unsigned short val_VTK_Type, unsigned short val_nPoly, unsigned short valIntegrationOrder);

    const vector<RDouble> *GetVertexXi(void) const;
    const vector<RDouble> *GetVertexEt(void) const;
    const vector<RDouble> *GetVertexZt(void) const;

    void Copy(const HighOrderStandardCell &other);

    //! @Determine shape function and its derivativew
    void DataStandardTetrahedron(void);
    void DataStandardPyramid(void);
    void DataStandardPrism(void);
    void DataStandardHexahedron(void);
};

//StandardFace
class HighOrderStandardFace : public HighOrderStandardElementBase
{
protected:
    unsigned short nPoly;             //! @the face grid is linear or quadratic
    unsigned short nNodes;            //! @the number of nodes on face

    vector< RDouble > vertexXiFace;   /*!< \brief r-location of the DOFs on side 0 of the face. the vertex of face using for shapefunction  */
    vector< RDouble > vertexEtFace;   /*!< \brief s-location of the DOFs on side 0 of the face */
    //[GaussPointsIdx][ShapeDOFsIdx]
    vector< vector< RDouble > > shapeFuncAtFaceIntegration;
                             
    vector< vector< RDouble > > dXiShapeFuncFaceIntegration;
    vector< vector< RDouble > > dEtShapeFuncFaceIntegration;

public:
    HighOrderStandardFace(){}
    ~HighOrderStandardFace(){}
    HighOrderStandardFace(unsigned short VTK_Type, unsigned short nFaceNodes,unsigned short integrationOrder);  //输入VTK_Type以及integrationOrder，输出constJacobian以及单元积分点个数、位置及权重系数
    HighOrderStandardFace(const HighOrderStandardFace &other) : HighOrderStandardElementBase(other) {Copy(other);}
    HighOrderStandardFace& operator= (const HighOrderStandardFace &other) {Copy(other); return (*this);}
    void Copy(const HighOrderStandardFace &other);

    /*!
    * \brief Function, which makes available the r-derivatives of the elements
    basis functions of side 0 in the integration points.
    * \return  The pointer to data, which stores this information.
    */
    RDouble * GetdXiShapeFunctionAtIntegration(int iGaussPoint)
    {
        return dXiShapeFuncFaceIntegration[iGaussPoint].data();
    }

    RDouble * GetdEtShapeFunctionAtIntegration(int iGaussPoint)
    {
        return dEtShapeFuncFaceIntegration[iGaussPoint].data();
    }
    /*!
    * \brief Function, which makes available the face basis functions of side 0
    in the integration points.
    * \return  The pointer to data, which stores this information.
    */
    const vector< vector< RDouble > > & GetShapeFuncFaceIntegration(void) const { return shapeFuncAtFaceIntegration; }
    /*!
    * \brief Function, which makes available the polynomial degree for this standard element.
    * \return  The polynomial degree of this standard element.
    */
    unsigned short GetNpoly(void) const { return nPoly; };
    unsigned short GetNnode(void) const { return nNodes; };
    
    //! @brief Function, which checks if the function arguments correspond to this standard element.
    //! @param[in] val_VTK_Type        - Type of the face using the VTK convention.
    //! @param[in] nFaceNodes          - the nodes of face (linear, quadratic).
    //! @param[in] integrationOrder    - integration order of face .
    //! @return Whether or not the function arguments correspond to this standard element.
    
    bool SameStandardElement(unsigned short val_VTK_Type, unsigned short val_nPoly, unsigned short integrationOrder);
    //! @Determine shape function and its derivativew
    void DataStandardTriangle(void);
    void DataStandardQuadrilateral(void);
};
typedef Eigen::Matrix<RDouble, Eigen::Dynamic, Eigen::Dynamic> HighOrderMatrix;
void Orthogonalize(int degreeOfFreedom, const HighOrderMatrix & taylorBasisMassMatrix, 
                   HighOrderMatrix & orthogonalizationMatrix);

}