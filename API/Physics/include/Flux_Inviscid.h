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
//! @file      Flux_Inviscid.h
//! @brief     Inviscid (convection) flux computing, including number of inviscid scheme,
//!            such as Roe, Steger-Warmming, Van Leer, HLLC, AUSM, etc.
//! @author    Bell, He Xin, Wang Boqian.

#pragma once
#include "PHVector3.h"
#include "FaceProxy.h"
#include "Geo_UnstructBC.h"    //! GMRES3D
#ifdef USE_GMRESSOLVER
#include "Sacado.hpp"
#endif
namespace PHSPACE
{
class InviscidSchemeParameter;
#ifdef USE_GMRESSOLVER
//! GMRESAD
typedef Sacado::Fad::DFad<RDouble> ADReal;
#endif
typedef void (*INVScheme) (FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

int GetSchemeID(const string &scheme_name);

INVScheme GetINVScheme(int scheme_id);

//! This Roe scheme is coding by the conservative variable way, in which the DF is expressed by dH.
void Roe_Scheme_ConservativeForm(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

#ifdef USE_GMRESSOLVER
//! GMRES, the flux is calculated in Cal_GMRES_Roe_Scheme_Flux 
//! assembling the Jacobian matrix using primitive variables
//! GMRESAD GMRESCSR
void GMRES_Roe_Scheme_ConservativeForm(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

//! GMRESJac1st GMRES3D
void GMRES_INVScheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

//! GMRESnolim
void GMRES_Roe_Scheme_ConservativeForm_nolim_2nd(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

//! GMRESnolim
void GMRES_Roe_Scheme_ConservativeForm_nolim_Matrix(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

//! GMRESAD_WITH_CSR GMRESCSR
void GMRES_Roe_Scheme_ConservativeForm_CSR(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

//! GMRESPV
void GMRES_Roe_Scheme_ConservativeForm_FD_PV(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

//! _CV : assembling the Jacobian matrix using conservative variables(consider the contribution from bc)
void GMRES_Roe_Scheme_ConservativeForm_CV(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

//! _FD: assembling the Jacobian matrix using conservative variables(not considering the contribution from bc)
void GMRES_Roe_Scheme_ConservativeForm_FD(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

//! _wholeJacobian: for outputing the wholeJacobian matrix, so remove the finite difference part
void GMRES_Roe_Scheme_ConservativeForm_wholeJacobian(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

//! _rhow: modification for the components of rhow when assembling the Jacobian matrix
void GMRES_Roe_Scheme_ConservativeForm_rhow(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);
void GMRES_Roe_Scheme_ConservativeForm_test(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

//! GMRESPV
template <typename T>
void Cal_GMRES_Roe_Scheme_Flux_PV(T* primL, T* primR,int iFace, T* flux,
    FaceProxy* face_proxy, InviscidSchemeParameter* invSchemePara);

//! GMRES3D 
template <typename T>
void Cal_GMRES_Steger_Scheme_Flux_PV(T* primL, T* primR,int iLen, T* flux,
    FaceProxy* face_proxy, InviscidSchemeParameter* invSchemePara);

//! GMRESnolim  primcc: primitive variables at the cell center, primf: primitive variables at the face
void cal_ReConstructFaceValue(RDouble* primcc, RDouble* dqdx, RDouble* dqdy, RDouble* dqdz,
                              RDouble dx, RDouble dy, RDouble dz, RDouble * primf,
                              InviscidSchemeParameter* invSchemePara);
//! using conservative variables
void Cal_GMRES_Roe_Scheme_Flux(RDouble* primL, RDouble* primR, RDouble* qCvL, RDouble* qCvR, int iFace, RDouble* flux,
    FaceProxy* face_proxy, InviscidSchemeParameter* invSchemePara);

//! GMRES CSR
int GetIndexOfBlockJacobianMatrix(vector<int> &AI, vector<int> &AJ, int nEquation, int rCell, int qCell);
#endif

//! This Roe scheme is coding by the primitive variable way, in which the DF is split into three parts.
void Roe_Scheme_PrimitiveForm(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Steger_Scheme (FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Ausmdv_Scheme (FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Ausmpw_Scheme (FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void AusmpwPlus_Scheme (FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Hlle_Scheme   (FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Vanleer_Scheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Rotate_Vanleer_Scheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Gks_Scheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);
//gks所需函数
void local(RDouble* loc, const RDouble* var, const Vector3D e1, const Vector3D e2, const Vector3D e3);
void localGrad(RDouble r[5][3], const RDouble grad[5][3], const Vector3D e1, const Vector3D e2, const Vector3D e3);
void chainRule(const RDouble gama, RDouble* swout, const RDouble* prim, const RDouble swin[5][3], const int& dir);
void calcMomentUvwxi(const RDouble& ck, const int& dims, const RDouble* prim, RDouble* mu, RDouble* mv, RDouble* mw, RDouble* mxi, RDouble* mu_l, RDouble* mu_r);
void GmUvwxipsi(RDouble result[5], const RDouble mu[7], const RDouble mv[7], const RDouble mw[7], const RDouble mxi[3],const std::vector<int> param);
void set_mt(const RDouble& dt, const RDouble gama, const RDouble prim[5], const RDouble primleft[5], const RDouble primright[5],int iFace, FaceProxy* faceProxy, RDouble *mt);
void microSlope(RDouble ms[5], const RDouble& ck, const RDouble& dims, const RDouble prim[5], const RDouble sw[5]);
void GmMsuvwxipsi(RDouble result[5], const RDouble ms[5], const RDouble mu[7], const RDouble mv[7], const RDouble mw[7], const RDouble mxi[3], const std::vector<int> param);
RDouble GmUvwxi1(const RDouble mu[7], const RDouble mv[7], const RDouble mw[7], const RDouble mxi[3], const std::vector<int> param);
RDouble GmMsuvwxi1(const RDouble ms[5], const RDouble mu[7], const RDouble mv[7], const RDouble mw[7], const RDouble mxi[3], const std::vector<int> param);
RDouble GmUvwxipsiq(const RDouble prim[5], const RDouble mu[7], const RDouble mv[7], const RDouble mw[7], const RDouble mxi[3], const std::vector<int> param);
RDouble GmMsuvwxipsiq(const RDouble prim[5], const RDouble ms[5], const RDouble mu[7], const RDouble mv[7], const RDouble mw[7], const RDouble mxi[3], const std::vector<int> param);

RDouble FMSplit1(const RDouble &mach, const RDouble &ipn);
RDouble FMSplit2(const RDouble &mach, const RDouble &ipn);
RDouble FMSplit4(const RDouble &mach, const RDouble &beta, const RDouble &ipn);
RDouble FPSplit5(const RDouble &mach, const RDouble &alpha, const RDouble &ipn);
RDouble FMmpn(const RDouble &mach, const RDouble &ipn);
RDouble FMppn(const RDouble &mach, const RDouble &ipn);

void RoeAveragedVariables(RDouble rhoLeft, RDouble uLeft, RDouble vLeft, RDouble wLeft, RDouble pLeft, RDouble gamaLeft, 
                          RDouble rhoRight, RDouble uRight, RDouble vRight, RDouble wRight, RDouble pRight, RDouble gamaRight,
                          RDouble &rhoRoe, RDouble &uRoe, RDouble &vRoe, RDouble &wRoe, RDouble &pRoe);

//! @brief The parameters that inviscid flux computing need.
//!        These data are compacted into this class to fill the need of OpenMP parallel.
//!        If not, would lead to data completion problems.
class InviscidSchemeParameter
{
public:
    InviscidSchemeParameter();
    ~InviscidSchemeParameter();

private:
    int nLen;
    int nm;
    int nLaminar;
    int nEquation;
    int nChemical;
    int nPrecondition;
    int isUnsteady;
    int nTemperatureModel;
    int nEntrFix;
    int nNitrogenIndex;
    int nElectronIndex;
    int nIdealState;
    RDouble kPrec;
    int limiterType;    //! GMRESnolim
	RDouble dt;    //!< time step, used in the inviscid flux computing.

    //! GMRESVenkat
    int limitModel;
    int limitVector;

    RDouble AusmpwPlusLimiter;

    RDouble refMachNumber;

    RDouble entrFixCoefLeft;
    RDouble entrFixCoefRight;

    RDouble *gamal; RDouble *gamar;
    RDouble *lrtl; RDouble *lrtr;
    RDouble *tcl; RDouble *tcr;
    RDouble *pcl; RDouble *pcr;
    RDouble **ql; RDouble **qr;
    RDouble **qlc; RDouble **qrc;    //! GMRESPassQC
    RDouble **tl; RDouble **tr;    //! left and right temperature.
    RDouble **flux;
    RDouble *xfn; RDouble *yfn; RDouble *zfn;
    RDouble* xfv; RDouble* yfv; RDouble* zfv;    //! GKS
    RDouble *xfc; RDouble *yfc; RDouble *zfc;    //! GMRESnolim
    RDouble *xcc; RDouble *ycc; RDouble *zcc;    //! GMRESnolim
    RDouble *nsx; RDouble *nsy; RDouble *nsz;    //! GMRESnolim absolute index
    RDouble *ns; RDouble *vol;   //! GMRESnolim absolute index
    RDouble ** prim; RDouble *gama;    //! GMRESnolim
    vector<int>* neighborCells;
    vector<int>* neighborFaces;
    vector<int>* neighborLR;
    vector<int> BCLeftCells;    //! GMRES2ndCorrection
    vector<int> BCRightCells;   //! GMRES2ndCorrection
    vector<int> BCFaces;        //! GMRES2ndCorrection
    RDouble **dqdx; RDouble **dqdy; RDouble **dqdz;    //! GMRESnolim
    int * qlsign; int * qrsign;    //! GMRESnolim
    int* nodeNumOfEachFace;    //!GKS
    UnstructBCSet* unstructBCSet;    //! GMRES3D
    RDouble *area;
    RDouble *vgn;
    RDouble averageVolume;    //! GMRESVenkat

public:
    void SetLength(const int nLen) {this->nLen = nLen;}
    void SetNumberOfLaminar(const int nLaminar) {this->nLaminar = nLaminar;}
    void Setnm(const int nm) {this->nm = nm;}
    void SetIfPrecondition(const int nPrecondition) {this->nPrecondition = nPrecondition;}
    void SetIfIsUnsteady(const int isUnsteady) {this->isUnsteady = isUnsteady;}
    void SetIfIdealState(const int nIdealState) { this->nIdealState = nIdealState; }
    void SetPreconCoefficient(const RDouble kPrec) { this->kPrec = kPrec; }
    void SetNumberOfChemical(const int nChemical) {this->nChemical = nChemical;}
    void SetNumberOfTemperatureModel(const int nTemperatureModel) {this->nTemperatureModel = nTemperatureModel;}
    void SetNumberOfTotalEquation(const int nEquation) {this->nEquation = nEquation;}
    void SetIndexOfNitrogen(const int nIndex) {this->nNitrogenIndex = nIndex;}
    void SetIndexOfElectron(const int nIndex) {this->nElectronIndex = nIndex;}
    void SetAusmpwPlusLimiter(const RDouble Limiter) {this->AusmpwPlusLimiter = Limiter;}
    void SetLimiterType(const int limiterType){ this->limiterType = limiterType;}    //!GMRESnolim
    void SetLimitModel(const int limitModel){ this->limitModel = limitModel;}    //! GMRESVenkat
    void SetLimitVector(const int limitVector){ this->limitVector = limitVector;}    //! GMRESVenkat
    void SetAverageVolume(const RDouble averageVolume){ this->averageVolume = averageVolume; }    //! GMRESVenkat
	void SetTimeStep(const RDouble dt) { this->dt = dt; }    //! time step, used in the inviscid flux computing.

    void SetEntropyFixMethod(const int nEntrFix) {this->nEntrFix = nEntrFix;}
    void SetEntropyFixCoefficients(const RDouble entrFixCoefLeft, const RDouble entrFixCoefRight)
    {this->entrFixCoefLeft = entrFixCoefLeft; this->entrFixCoefRight = entrFixCoefRight;}
    void SetMachNumber(const RDouble refMachNumber) {this->refMachNumber = refMachNumber;}
    void SetLeftAndRightGama(RDouble *gamal, RDouble *gamar) {this->gamal = gamal; this->gamar = gamar;}
    void SetLeftAndRightPressureCoefficient(RDouble *lrtl, RDouble *lrtr) {this->lrtl = lrtl; this->lrtr = lrtr;}
    void SetLeftAndRightTimeCoefficient(RDouble *tcl, RDouble *tcr) {this->tcl = tcl; this->tcr = tcr;}
    void SetLeftAndRightPreconCoefficient(RDouble *pcl, RDouble *pcr) {this->pcl = pcl; this->pcr = pcr;}
    void SetLeftAndRightQ(RDouble **ql, RDouble **qr) {this->ql = ql; this->qr = qr;}
    void SetLeftAndRightQC(RDouble **qlc, RDouble **qrc) {this->qlc = qlc; this->qrc = qrc;}    //! GMRESPassQC
    void SetLeftAndRightQSign(int* qlsign, int* qrsign) {this->qlsign = qlsign; this->qrsign = qrsign;}    //! GMRESnolim GMRESSign
	void SetNodeNumOfEachFace(int* nodeNumOfEachFace) {this->nodeNumOfEachFace = nodeNumOfEachFace;}    //! GKS
    void SetLeftAndRightTemperature(RDouble **tl, RDouble **tr) {this->tl = tl; this->tr = tr;}
    void SetFlux(RDouble **flux) {this->flux = flux;}
    void SetFaceNormal(RDouble *xfn, RDouble *yfn, RDouble *zfn) {this->xfn = xfn; this->yfn = yfn; this->zfn = zfn;}
    void SetFaceArea(RDouble *area) {this->area = area;}
    void SetFaceVelocity(RDouble *vgn) {this->vgn = vgn;}
    //! GKS
    //void SetFaceVector(RDouble *xfv, RDouble *yfv, RDouble *zfv) {this->xfv = xfv; this->yfv = yfv; this->zfv = zfv; }
    /*void SetFaceVector(RDouble* xfv, RDouble* yfv, RDouble* zfv, int nTotalFace) {
        this->xfv = new RDouble[nTotalFace];
        this->yfv = new RDouble[nTotalFace];
        this->zfv = new RDouble[nTotalFace];
        //xfv.resize(nTotalFace);
        //yfv.resize(nTotalFace);
        //zfv.resize(nTotalFace);
        for (int i = 0; i < nTotalFace; ++i) {
            this->xfv[i] = xfv[i];
            this->yfv[i] = yfv[i];
            this->zfv[i] = zfv[i];
        }
    }
    void DelFaceVector()
    {
        delete[] xfv;    xfv = NULL;
        delete[] yfv;    yfv = NULL;
        delete[] zfv;    zfv = NULL;
    }*/
    //! GMRESnolim
    void SetFaceCenter(RDouble *xfc, RDouble *yfc, RDouble *zfc) {this->xfc = xfc; this->yfc = yfc; this->zfc = zfc;}
    void SetCellCenter(RDouble *xcc, RDouble *ycc, RDouble *zcc) {this->xcc = xcc; this->ycc = ycc; this->zcc = zcc;}
    void SetGradient(RDouble ** dqdx, RDouble ** dqdy, RDouble ** dqdz) {this->dqdx = dqdx; this->dqdy = dqdy; this->dqdz = dqdz;}
    void SetFaceNormalAbs(RDouble *nsx, RDouble *nsy, RDouble *nsz) {this->nsx = nsx; this->nsy = nsy; this->nsz = nsz;}
    void SetFaceAreaAbs(RDouble *ns) {this->ns = ns;}
    void SetVolume(RDouble *vol) {this->vol = vol;}
    void SetPrimitive(RDouble **prim){this->prim = prim;}
    void SetGama(RDouble* gama){this->gama = gama;}
    void SetNeighborCells(vector<int>* neighborCells){this->neighborCells=neighborCells;}
    void SetNeighborFaces(vector<int>* neighborFaces){this->neighborFaces=neighborFaces;}
    void SetNeighborLR(vector<int>* neighborLR){this->neighborLR=neighborLR;}

    //! GMRES2ndCorrection
    void SetBCLeftCells(vector<int> BCLeftCells){this->BCLeftCells = BCLeftCells;}
    void SetBCRightCells(vector<int> BCRightCells){this->BCRightCells = BCRightCells;}
    void SetBCFaces(vector<int> BCFaces){this->BCFaces = BCFaces;}

    //! GMRES3D
    void SetUnstructBCSet(UnstructBCSet* unstructBCSet){this->unstructBCSet = unstructBCSet;}


    int GetLength() {return this->nLen;}
    int Getnm() {return this->nm;}
    int GetNumberOfLaminar() {return this->nLaminar;}
    int GetNumberOfTotalEquation() {return this->nEquation;}
    int GetNumberOfChemical() {return this->nChemical;}
    int GetNumberOfTemperatureModel() {return this->nTemperatureModel;}
    int GetIndexOfNitrogen() {return this->nNitrogenIndex;}
    int GetIndexOfElectron() {return this->nElectronIndex;}

    int GetEntropyFixMethod() {return this->nEntrFix;}
    int GetIfPrecondition() {return this->nPrecondition;}
    int GetIfIsUnsteady() {return this->isUnsteady;}
    int GetIfIdealState() { return this->nIdealState; }
    RDouble GetPreconCoefficient() { return this->kPrec; }
    int GetLimiterType(){ return this->limiterType;}    //! GMRESnolim
    int GetLimitModel(){ return this->limitModel;}    //! GMRESVenkat
    int GetLimitVector(){ return this->limitVector;}    //! GMRESVenkat
    RDouble GetAverageVolume(){ return this->averageVolume; }    //! GMRESVenkat
	RDouble GetTimeStep() { return this->dt; }    //! time step, used in the inviscid flux computing.

    RDouble GetAusmpwPlusLimiter() {return this->AusmpwPlusLimiter;}

    RDouble GetLeftEntropyFixCoefficient() {return this->entrFixCoefLeft;}
    RDouble GetRightEntropyFixCoefficient() {return this->entrFixCoefRight;}
    RDouble GetMachNumber() {return this->refMachNumber;}
    RDouble *GetLeftGama() {return this->gamal;}
    RDouble *GetRightGama() {return this->gamar;}

    RDouble *GetLeftPressureCoefficient() {return this->lrtl;}
    RDouble *GetRightPressureCoefficient() {return this->lrtr;}

    RDouble *GetLeftTimeCoefficient() {return this->tcl;}
    RDouble *GetRightTimeCoefficient() {return this->tcr;}

    RDouble *GetLeftPreconCoefficient() {return this->pcl;}
    RDouble *GetRightPreconCoefficient() {return this->pcr;}

    RDouble **GetLeftQ() {return this->ql;}
    RDouble **GetRightQ() {return this->qr;}

    RDouble **GetLeftQC() {return this->qlc;}    //! GMRESPassQC
    RDouble **GetRightQC() {return this->qrc;}   //! GMRESPassQC

    int* GetLeftQSign() {return this->qlsign;}    //! GMRESSign GMRESnolim
    int* GetRightQSign() {return this->qrsign;}   //! GMRESSign GMRESnolim

	int* GetNodeNumOfEachFace() {return this->nodeNumOfEachFace; }    //! GKS

    RDouble ** GetLeftTemperature() {return this->tl;}
    RDouble ** GetRightTemperature() {return this->tr;}

    RDouble **GetFlux() {return this->flux;}
    RDouble *GetFaceNormalX() {return this->xfn;}
    RDouble *GetFaceNormalY() {return this->yfn;}
    RDouble *GetFaceNormalZ() {return this->zfn;}
    RDouble *GetFaceArea() {return this->area;}
    RDouble *GetFaceVelocity() {return this->vgn;}
    //! GKS
    RDouble *GetFaceVectorX() {return this->xfv;}  // 法向量X分量（非结构网格的面向量）
    RDouble *GetFaceVectorY() {return this->yfv;}  // 法向量Y分量
    RDouble *GetFaceVectorZ() {return this->zfv;}  // 法向量Z分量
    //! GMRESnolim
    RDouble *GetFaceCenterX() {return this->xfc;}
    RDouble *GetFaceCenterY() {return this->yfc;}
    RDouble *GetFaceCenterZ() {return this->zfc;}
    RDouble *GetCellCenterX() {return this->xcc;}
    RDouble *GetCellCenterY() {return this->ycc;}
    RDouble *GetCellCenterZ() {return this->zcc;}
    RDouble **GetGradientX() {return this->dqdx;}
    RDouble **GetGradientY() {return this->dqdy;}
    RDouble **GetGradientZ() {return this->dqdz;}
    RDouble *GetFaceNormalXAbs() {return this->nsx;}
    RDouble *GetFaceNormalYAbs() {return this->nsy;}
    RDouble *GetFaceNormalZAbs() {return this->nsz;}
    RDouble *GetFaceAreaAbs(){return this->ns;}
    RDouble *GetVolume(){return this->vol;}
    RDouble **GetPrimitive(){return this->prim;}
    RDouble *GetGama(){return this->gama;}
    vector<int>* GetNeighborCells() {return this->neighborCells;}
    vector<int>* GetNeighborFaces() {return this->neighborFaces;}
    vector<int>* GetNeighborLR() {return this->neighborLR;}

    //! GMRES2ndCorrection
    vector<int> GetBCLeftCells() {return this-> BCLeftCells;}
    vector<int> GetBCRightCells() {return this-> BCRightCells;}
    vector<int> GetBCFaces() {return this-> BCFaces;}

    //! GMRES3D
    UnstructBCSet* GetUnstructBCSet() {return this->unstructBCSet;}
};

}
