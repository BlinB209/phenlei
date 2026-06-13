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
//! @file      FaceProxy.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "TypeDefine.h"
#include "GeneralFieldProxy.h"
#include "Geo_Grid.h"
#include "GeomProxy.h"
#include "Constants.h"

namespace PHSPACE
{
class NSFaceValue
{
private:
    GeneralFieldProxy *rho_ds, *hint_s;
    GeneralFieldProxy *prim, *t;
    RDouble *kcp, *mul, *mut;
public:
    NSFaceValue(int neqn, int ns, int nlen);
    ~NSFaceValue();
    GeneralFieldProxy * GetRhoDS() { return rho_ds; }
    GeneralFieldProxy * GetHintS() { return hint_s; }
    GeneralFieldProxy * GetPrim()  { return prim; }
    GeneralFieldProxy * GetT()     { return t; }

    //! Heat conductivity coefficient, usually computed by Cp.
    RDouble * GetKCP() { return kcp; }
    RDouble * GetMUL() { return mul; }
    RDouble * GetMUT() { return mut; }
};

class TurbFaceValue
{
public:
    RDouble **prim, **mlt;
    RDouble *mul, *mut;
public:
    TurbFaceValue(int neqn, int nlen);
    ~TurbFaceValue();
};

class TransitionFaceValue
{
public:
    RDouble **prim, **mlt;
    RDouble *mul, *mut;
public:
    TransitionFaceValue(int neqn, int nlen);
    ~TransitionFaceValue();
};

class FaceProxy
{
public:
    //! Type definitions.
    typedef RDouble **value_type;
private:
    value_type ql, qr, flux;
    value_type qlc,qrc; //GMRESPassQC , pass the value at cell center
    RDouble **tl, **tr;
    RDouble** dRdq; // GMRES
    RDouble** dRdq1st; // GMRESJac1st
    RDouble** dDdP;  // GMRESBoundary
    RDouble *deltl, *deltr;
    RDouble *gamal, *gamar;
    RDouble *lrtl, *lrtr;    //! Pressure factor.
    int* leftcellindexofFace, *rightcellindexofFace; //GMRES
    int* qlsign, *qrsign; //GMRESnolim GMRESSign , according the sign for limiter
    RDouble *tcl, *tcr;
    RDouble *pcl, *pcr;
    NSFaceValue   *ns_facevar;
    TurbFaceValue *turb_facevar;
    TransitionFaceValue *transition_facevar;
    GeomProxy     *geom_proxy;
    int nlen, nsize, neqn;
    int localStart, localEnd, nBoundFace, nTotalCell;  // GMRESBoundary
    vector<int>* wallFaceIndex; // GMRESBCorrection , this is for the correction of wall bc
    vector<int> JacobianAI; // GMRESCSR
    vector<int> JacobianAJ; // GMRESCSR
    vector<int> JacobianAI1st; // GMRESJac1st
    vector<int> JacobianAJ1st; // GMRESJac1st
    int JacOrder; // GMRESJac1st
    int isViscous;  // GMRESBCorrection
    FaceProxy *next;
public:
    FaceProxy();
    ~FaceProxy();
public:
    FaceProxy * GetNext() { return next; }
    void SetNext(FaceProxy *next) { this->next = next; }
public:
    void Create(int nsize, int neqn);
    value_type GetQL();
    value_type GetQR();
    value_type GetQLC(); // GMRESPassQC
    value_type GetQRC(); // GMRESPassQC
    value_type GetFlux();

    RDouble ** GetLeftTemperature();
    RDouble ** GetRightTemperature();
    RDouble ** GetJacobianMatrix(); //GMRES
    void SetJacobianMatrix(RDouble** dRdqfromGrid); //GMRES

    //! GMRESJac1st
    RDouble ** GetJacobianMatrix1st();
    void SetJacobianMatrix1st(RDouble** dRdq1stfromGrid);

    RDouble ** GetdDdPMatrix(); // GMRESBoundary
    void SetdDdPMatrix(RDouble** dDdPfromGrid); // GMRESBoundary

    // GMRESBoundary
    int GetlocalStart() {return localStart;};  
    int GetlocalEnd() {return localEnd;};      
    int GetnBoundFace() {return nBoundFace;};  
    int GetnTotalCell() {return nTotalCell;};  
    int GetWallType() {return isViscous;}; // GMRESBCorrection
    vector<int>* GetWallFaceIndex() const { return wallFaceIndex;}; // GMRESBCorrection
    void SetlocalStart(int localStart) {this->localStart = localStart;}; 
    void SetlocalEnd(int localEnd) {this->localEnd = localEnd;};         
    void SetnBoundFace(int nBoundFace) {this->nBoundFace = nBoundFace;};
    void SetnTotalCell(int nTotalCell) {this->nTotalCell = nTotalCell;};
    void SetWallType(int isViscous) {this->isViscous = isViscous;}; // GMRESBCorrection
    void SetWallFaceIndex(vector<int>* wallFaceIndex){this->wallFaceIndex = wallFaceIndex;};// GMRESBCorrection


    // GMRESCSR
    vector<int> GetJacobianAI4GMRES() const { return JacobianAI; };
    vector<int> GetJacobianAJ4GMRES() const { return JacobianAJ; };
    void SetJacobianAI4GMRES(vector<int> JacobianAI) { this->JacobianAI = JacobianAI; };
    void SetJacobianAJ4GMRES(vector<int> JacobianAJ) { this->JacobianAJ = JacobianAJ; };
    //! GMRESJac1st
    vector<int> GetJacobianAI1st4GMRES() const { return JacobianAI1st; };
    vector<int> GetJacobianAJ1st4GMRES() const { return JacobianAJ1st; };
    void SetJacobianAI1st4GMRES(vector<int> JacobianAI1st) { this->JacobianAI1st = JacobianAI1st; };
    void SetJacobianAJ1st4GMRES(vector<int> JacobianAJ1st) { this->JacobianAJ1st = JacobianAJ1st; };

    int GetJacobianOrder() const { return JacOrder; };
    void SetJacobianOrder(int JacOrder) { this->JacOrder = JacOrder; };

    int size() { return nlen; }
    int capacity() { return nsize; }
    void setsize(int nlen) { this->nlen = nlen; }

    // GMRES
    int* GetLeftCellIndexOfFace();
    int* GetRightCellIndexOfFace();

    // GMRESnolim GMRESSign
    int* GetQLSign();
    int* GetQRSign();

    //! Get left cell gama of the face.
    RDouble * GetGamaL();

    //! Get right cell gama of the face.
    RDouble * GetGamaR();

    //! Get left cell pressure coefficient of the face.
    RDouble * GetPressureCoefficientL();

    //! Get right cell pressure coefficient of the face.
    RDouble * GetPressureCoefficientR();

    //! Get left cell time coefficient for precondition of the face.
    RDouble * GetTimeCoefficientL();

    //! Get right cell time coefficient for precondition of the face.
    RDouble * GetTimeCoefficientR();

    //! Get left cell precondition coefficient of the face.
    RDouble * GetPreconCoefficientL();

    //! Get right cell precondition coefficient of the face.
    RDouble * GetPreconCoefficientR();

    //! Get left cell weight of the face.
    RDouble * GetWeightL();

    //! Get right cell weight of the face.
    RDouble * GetWeightR();

    GeomProxy * GetGeomProxy();
    void SetGeomProxy(GeomProxy *geom_proxy);

    NSFaceValue * GetNSFaceValue();
    void SetNSFaceValue(NSFaceValue *facevar);

    TurbFaceValue * GetTurbFaceValue();
    void SetTurbFaceValue(TurbFaceValue *facevar);

    TransitionFaceValue * GetTransitionFaceValue();
    void SetTransitionFaceValue(TransitionFaceValue *facevar);
};

//! Check whether the density and pressure are positive.
inline bool PositiveCheckForDensityAndPressure(RDouble *qTry)
{
    if (qTry[IDX::IR] <= 0.0 || qTry[IDX::IP] <= 0.0) return false;
    return true;
}

//! ReConstruct the left and right cell variables (Ql and Qr) of the face for structured grid.
void ReConstructQlQr_STR(Grid *grid_in, FaceProxy *face_proxy, int neqn, int nst, int ned);
}
