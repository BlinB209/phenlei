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
//! @file      Constants.h
//! @brief     Define some global const value.
//! @author    He Xin.

#pragma once
#include "Precision.h"

namespace PHSPACE
{
const int SOLVER_BASED        = 0;
const int GRID_BASED          = 1;

const int DENSITY_BASED_SOLVER = 0;
const int PRESSURE_BASED_SOLVER = 1;

const int FV_METHOD           = 0;
const int FD_METHOD           = 1;

const int SOLVE_FIELD         = 0;
const int CREATE_GRID         = 1;
const int CAL_WALL_DIST       = 2;
const int PARTITION_GRID      = 3;
const int REPOSITORY          = 4;  //! Useless.
const int OVERSETGRID_VIEW    = 5;
const int OVERSET_CONFIG      = 6;
const int CONVERT_PARTITION   = 7;
const int DSMC_SOLVER         = 8;  //! Useless.
const int HELICOPTER          = 9;  //! Useless.
const int UGKS_SOLVER         = 10; //! Useless.
const int HO_SOLVER           = 12;
const int LBM_SOLVER_MPI      = 13;
const int LBM_SOLVER_OMP      = 17;
const int INTEGRATIVE_SOLVER  = 14;
const int GRIDCHECK_SOLVER    = 15;
const int GRID_SOLVER         = 16;

const int SPECDIFFHYB_SOLVER  = 51;
const int SPEC_SOLVER         = 52;

const int DERIVATIVE_IDENITFY = 98;
const int POST_PROCESSING     = 99;

const int NOLES_SOLVER        = 0;
const int LES_SOLVER          = 1;

const int GRID_CONVERSION     = 1;
const int GRID_REFINE         = 2;
const int GRID_MERGING        = 3;
const int GRID_DEFORMATION    = 4;
const int GRID_REPAIR         = 5;
const int GRID_MIRROR         = 6;
const int GRID_STRTOUNS       = 7;

//! Type of  rotating reference frame.
const int STATIONARY_FRAME = 0;
const int TRANSLATIONAL_FRAME = 1;
const int ROTATIONAL_FRAME = 2;

//! Type of InitFlow.
const int INITINFLOW      = 0;
const int RESTARTFLIE     = 1;
const int INTERPOLATEFLIE = 2;

//! periodicType
const int NO_PERIODICITY = 0;
const int TRANSLATIONAL_PERIODICITY = 1;
const int ROTATIONAL_PERIODICITY = 2;

const RDouble half            = 0.5;
const RDouble third           = 1.0 / 3.0;
const RDouble two3rd          = 2.0 / 3.0;
const RDouble four3rd         = 4.0 / 3.0;
const RDouble fourth          = 0.25;
const RDouble sixth           = 1.0 / 6.0;
const RDouble eighth          = 1.0 / 8.0;
const RDouble tenth           = 0.1;
const RDouble zero            = 0.0;
const RDouble one             = 1.0;
const RDouble two             = 2.0;
const RDouble three           = 3.0;
const RDouble four            = 4.0;
const RDouble five            = 5.0;
const RDouble six             = 6.0;
const RDouble seven           = 7.0;
const RDouble eight           = 8.0;
const RDouble nine            = 9.0;
const RDouble ten             = 10.0;
const RDouble eleven          = 11.0;
const RDouble twelve          = 12.0;
const RDouble thirteen        = 13.0;
const RDouble fourteen        = 14.0;
const RDouble fifteen         = 15.0;
const RDouble sixteen         = 16.0;
const RDouble seventeen       = 17.0;
const RDouble eighteen        = 18.0;
const RDouble nineteen        = 19.0;
const RDouble twenty          = 20.0;

const RDouble PI              = 3.141592653589793;
const RDouble straightAngle   = 180.0;

const int SEGCTION_LENGTH     = 10240;
const int MAXCELLLAYER        = 1024;

const int NS_EQUATION         = 1;
const int TURBULENCE          = 2;
const int TRANSITION          = 3;

const int NONDIMENSIONCONDITION = 0;
const int FLIGHTCONDITION       = 1;
const int EXPERIMENTCONDITION   = 2;
const int SUBSONICBCCONDITION   = 3;    //! Useless.
const int TEMPERATURE_DENSITY   = 4;
const int TEMPERATURE_PRESSURE  = 5;
const int MACH_TEMP_PRE         = 6;
const int WEATHERCONDITION      = 7;
const int WINDSPEEDPROFILE      = 8;
const int PRIMITIVE_VARIABLES   = 9;
const int TRAJECTORY            = 10;

const int POWERLAW              = 0;
const int LOGARITHMICLAW        = 1;

const int VIS_STD               = 0;
const int VIS_AVER              = 1;
const int VIS_TEST              = 2;
const int VIS_NEW1              = 3;
const int VIS_NEW2              = 4;
const int VIS_ERIC              = 5;

const int ISCHEME_ROE            = 1;
const int ISCHEME_VANLEER        = 2;
const int ISCHEME_STEGER         = 3;
const int ISCHEME_HLLE           = 4;
const int ISCHEME_LAX_FRIEDRICHS = 5;
const int ISCHEME_AUSM           = 6;
const int ISCHEME_AUSMP          = 7;
const int ISCHEME_AUSMPUP        = 8;
const int ISCHEME_AUSMDV         = 9;
const int ISCHEME_AUSM_W         = 10;
const int ISCHEME_AUSMPW         = 11;
const int ISCHEME_Rotate         = 12;
const int ISCHEME_ROE_MODIFIED   = 13;
const int ISCHEME_AUSMPW_PLUS    = 14;
const int ISCHEME_GMRES_ROE      = 15;    //! GMRES modified Roe scheme for assembling Jacobian matrix using FD
const int ISCHEME_GMRES_Steger   = 16;    //! GMRES modified Steger scheme  GMRES3D
const int ISCHEME_GKS            = 17;

const int FIRST_ORDER            = 1;
const int SECOND_ORDER           = 2;

const int ILMT_NOLIM             = 1;
const int ILMT_FIRST             = 2;
const int ILMT_BARTH             = 3;
const int ILMT_VENCAT            = 4;
const int ILMT_STRUCT            = 10;
const int ILMT_MINMOD            = 11;
const int ILMT_VAN_ALBADA        = 12;
const int ILMT_VAN_ALBADA_CLZ    = 13;
const int ILMT_VANLEER           = 14;
const int ILMT_SMOOTH            = 15;
const int ILMT_MINMODTEST        = 16;
const int ILMT_MIN_VAN           = 17;
const int ILMT_3rd_SMOOTH        = 18;
const int ILMT_3rd_Minmod_SMOOTH = 19;
const int ILMT_WENO3_JS          = 20;
const int ILMT_WENN3_PRM211      = 21;
const int ILMT_WENN3_ZM          = 22;
const int ILMT_WENN3_ZES2        = 23;
const int ILMT_WENN3_ZES3        = 24;

const int INVISCID               = 0;
const int LAMINAR                = 1;
const int ALGEBRAIC              = 2;
const int ONE_EQU                = 3;
const int TWO_EQU                = 4;

const int IREGAMA                = 2;
const int BALDWIN_LOMAX          = 10;
const int SPALART_ALLMARAS       = 11;
const int KW_SST                 = 12;
const int KEPSILON               = 13;

const int INCOMPRESSIBLE         = 0;
const int COMPRESSIBLE           = 1;

const int MULTI_STAGE            = 1;
const int LU_SGS                 = 4;
const int Block_LU_SGS           = 5;
const int JACOBIAN_ITERATION     = 6;
const int Line_LU_SGS            = 7;
const int Matrix_LU_SGS          = 8;
const int GMRES                  = 9;  

const int OneRefineOneDeform     = 1;
const int DeformAfterAllRefine   = 2;

const int X_DIR    = 1;
const int Y_DIR    = 2;
const int Z_DIR    = 3;

const int PHINT    = 1;
const int PHFLOAT  = 2;
const int PHDOUBLE = 3;
const int PHSTRING = 4;
const int PHBOOL   = 5;

const RDouble SAME_POINT_TOL = static_cast<RDouble> (1.0e-12);

const int PROBESMONITOR   = 0;
const int LINESMONITOR    = 1;
const int SURFACESMONITOR = 2;

const int NEARESTCELLDATA = 0;
const int REALCELLDATA    = 1;

const int CELLVALUE       = 0;
const int CELLSTOPROBE    = 1;
const int NODESTOPROBE    = 2;

const int NON_GHOST       = 0;    //! These are interior cells.
const int ONE_GHOST_LAYER = 1;    //! Their have one-layer ghost cells.
const int TWO_GHOST_LAYER = 2;    //! Their have two-layer ghost cells.
const int THREE_GHOST_LAYER = 3;    //! Their have three-layer ghost cells.

const int GlobalCoordinate = 0;
const int LocalCoordinate  = 1;

//! These variables are for porous media algorithm.
const int NON_POROUS = 0;
const int POROUS = 1;

//! these variables are for overset algorithm.
const int INACTIVE        = 0;
const int ACTIVE          = 1;
const int INTERPOLATION   = -1;
const int VARIABLEACTIVE  = 2;
const int DONORCELL       = 3;
const int EXPLICIT_CONFIG = 0;
const int IMPLICIT_CONFIG = 1;
namespace IDX
{
const int IR   = 0;
const int IU   = 1;
const int IV   = 2;
const int IW   = 3;
const int IP   = 4;

const int IRU  = 1;
const int IRV  = 2;
const int IRW  = 3;
const int IRE  = 4;

const int IA   = 0;
const int IX   = 1;
const int IY   = 2;
const int IZ   = 3;

const int IKES = 1;
const int IETA = 2;
const int IZTA = 3;

const int ISA  = 0;

const int IKE  = 0;
const int IKW  = 1;

const int IGAMA = 0;
const int IRECT = 1;

const int ITT  = 0;
const int ITV  = 1;
const int ITE  = 2;

// these variables are for SIMPLE algorithm.
const int S_IU = 0;
const int S_IV = 1;
const int S_IW = 2;
const int S_IP = 3;
const int S_IR = 4;
const int S_ITEMP  = 5;
const int S_ITURBK = 6;
const int S_IEPSILON  = 7;

const int S_IOMEGA = 8;
const int S_ISA = 9;
const int S_IEDDYVISC = 10;

const int S_IAIR            = 15;
const int S_IH2             = 16;
const int S_IH2_l           = 17;
const int S_INH3_l          = 18;
const int S_ICH4            = 19;
const int S_IC2H4           = 20;
const int S_IC3H8           = 21;
const int S_ING             = 22;
const int S_ILNG            = 23;
const int S_IH2S            = 24;
const int S_ICl2            = 25;
const int S_INH3            = 26;
const int S_IC2H6           = 27;

const int MOMX_EQ    = 0;
const int MOMY_EQ    = 1;
const int MOMZ_EQ    = 2;
const int CONT_EQ    = 3;
const int ENER_EQ    = 4;
const int ZNUT_EQ    = 5;
const int TURBK_EQ   = 6;
const int EPSILON_EQ = 7;
const int OMEGA_EQ   = 8;
}
const int MAX_SPECIES_NUM  = 20;     //! The maximum number of species.
const int MAX_REACTION_NUM = 50;    //! The maximum number of reactions.

namespace WALLFUNCTION
{
    const int NONE     = 0;
    const int STANDARD = 1;
    const int PAB3D    = 2;
}

namespace ABLATION
{
const RDouble FormationHeat = 13.083E+07;
}

#define BOLTZMANN_CONSTANT  1.380649E-23      //! The Boltzmann constant, Unit: J/K.
#define AVOGADRO_CONSTANT   6.0225E+23        //! The Avogadro constant, Unit: mol-1.
#define VACUUM_PERMITTIVITY 8.854187817E-12   //! The permittivity of vacuum,Unit: C/(V*m).
#define STANDARD_PRESSURE   101325.0          //! The standard atmosphere pressure,Unit: Pa.
#define ELECTRIC_QUANTITY   1.602176634E-19   //! The electric quantity, Unit: C.
#define GeneralGasConstant  8.31434           //! The general gas constant, J/(kg * K).
}

//---------------------------
// Pressure-based enum
//----------------------------
typedef enum  
{
    Projection = 1,
    SIMPLE     = 2,
    SIMPLEC    = 3,
    PISO       = 4
}InCom_method;

typedef enum
{
    HYPRELib = 0,
    UNAPLib  = 1,
    YHAMGLib = 2
}mathLibName;

typedef enum  
{
    GaussSeidelMethod = 0,
    GMRESMethod       = 1,
    BicgstabMethod    = 2,
    AMGMethod         = 3,
    PCGMethod         = 4,
    FlexGMRESMethod   = 5,
    LGMRESMethod      = 6,
    COGMRESMethod     = 7,
    CGNRMethod        = 8,
    MGRMethod         = 9,
    ILUMethod         = 10,
    AMGDDMethod       = 11,
    HybridMethod       = 12,
    PipeBicgstabMethod = 13,
    PipeCGMethod = 14,
    PipeGMRESMethod = 15
}LinearSolverID;

typedef enum  
{
    L1Res = 1,
    L2Res    = 2,
    LinfRes  = 3,
}ResidualID;

typedef enum
{
    iNone = 0,
    iILU = 1,
    iAMG = 2,
    iEuclid = 3,
    DiagScale = 4,
    PILUT = 5,
    GSMG = 6,
    ParaSails = 7,
    LduDiag = 8,
    LduDIC = 9,
    LduDILU = 10,
    LduBJacobi = 11,
    iBJILU =   12,
    iBJSOR =   13,
    iJacobi =  14,
    iSOR    =  15,
    iAMGDD = 16
}PrecondTypeID;

typedef enum
{
    None                  = 0,
    MinimumCorrection     = 1,
    OrthogonalCorrection  = 2,
    OverrelaxedCorrection = 3
} InCom_nonOrthMethod;

typedef enum
{
    firstUpWind   = 0,
    secondUpWind  = 1,
    QUICKMTHOD    = 2,
    SecondCentral = 3,
} invFluxScheme;

typedef enum
{
    iDefScheme            = 0,
    iCentralScalar        = 1,
    iCentralScalarReconst = 5,
    iCentralMatrix        = 2,
    iUpwindAUSMPlus       = 3,
    iUpwindRoe            = 4
}SpaceScheme;

typedef enum
{
    iEuler      = 0,
    iNSLam      = 1,
    iSAWallF    = 2,
    iKEStandard = 3,
    iSSTWallF   = 4
}Model;

typedef enum
{
    WriteBoundary,
    WriteBlock,
    Writecomplete
} plotFieldType;

typedef enum
{
    Tecplot_Binary,
    Tecplot_ASCII,
    Ensight_Binary,
    Ensight_ASCII,
    Paraview
} VisualfileType;

namespace TEC_SPACE
{
const int ORDERED         = 0;
const int FELINESEG       = 1;
const int FETRIANGLE      = 2;
const int FEQUADRILATERAL = 3;
const int FETETRAHEDRON   = 4;
const int FEBRICK         = 5;
const int FEPOLYGON       = 6;
const int FEPOLYHEDRON    = 7;

const int BoundaryVisual  = 0;
const int FieldVisual     = 1;
const int BlockVisual     = 2;

const int WriteBoundary = 0;
const int WriteBlock = 1;
const int Writecomplete = 2;
}

