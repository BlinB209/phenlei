//! @file      main.cpp.
//! @brief     PHengLEI's main entrance.
//! @author    Bell.
//!                     1. Introduction of PHengLEI
//!    PHengLEI is a Computational Fluid Dynamics (CFD) software platform developed
//! by China Aerodynamics Research and Development Center (CARDC).
//! PHengLEI is a structured/unstructured
//! general brand CFD software, meaning the structured and unstructured solvers
//! can run on structured and unstructured grid respectively. PHengLEI is developed on
//! the basis of existing in-house codes and software in CARDC, and making full use of
//! years of technical reserves on the developing of professional software.
//! The range of applied flow velocity covers low-subsonic-transonic-hypersonic speed.
//! Modern software engineering method was employed as the guidance in the process of
//! the software development, and the software architecture for next generation was designed
//! in consideration of the CFD industrial features.

//!                     2. Development history
//!    Since 2003, a structured/unstructured hybrid prototype codes named 'Fantasy' was developed
//! by Dr. Eric (He Xin) and Dr. Zhang Laiping, supported by National Basic Research Program of China,
//! (No. 2009CB723802).
//!    Since 2011, a professional CFD software development team was assembled to develop general
//! CFD software which is named 'HyperFLOW' firstly, supported by National Natural Science Foundation
//! of China (No. 11272339).
//!    2013.12.28, the general CFD software 'PHengLEI' was released to CARDC users.
//!    2016.04.28, PHengLEI (version 1.0) was officially released to CFD users in the field of aeronautic
//! and aerospace in CHINA. Several versions were released in the next few yeas till 2018.12.
//!    2018.12.07, Computational Fluid Dynamics Software Engineering Center was set up in CARDC,
//! to develop National Numerical Wind-tunnel (NNW) project. PHengLEI, as the most important part of NNW,
//! is entrusted with a new mission.
//!    2020, Open source.
//!    2023, the PHengLEI development team relocated to Chengdu.

//!                     3. Acknowledgements
//!    During the development of prototype before 2003, we got a lot of help from Dr. Wang and other
//! researchers. PHengLEI is developed from HyperFLOW v1.0, which is mostly developed by Eric.
//!    Thanks for all!

//! Declaration: If the user uses the "PHengLEI" software for academic research or engineering application,
//! the software based on the National Numerical Windtunnel should be marked in the prominent position of 
//! the paper results, and references related to the "PHengLEI" software (e.g. [1] and [2]) should be cited.
//! 
//! [1] Zhao Z, et al. Design of general CFD software PHengLEI[J]. Computer Engineering & Science, 2020, 42(2): 
//!     210-219.
//! [2] Zhao Z, et al. PHengLEI: A Large Scale Parallel CFD Framework for Arbitrary Grids[J]. Chinese Journal 
//!     of Computers, 2018, 42(11): 2368-2383.
#include "CFDSimulation.h"
#include "TK_Time.h"
using namespace std;
using namespace PHSPACE;

int main()
{
    //! Do some prepare for time statistic.
    setvbuf(stdout, NULL, _IONBF, 0);
    PHSPACE::TIME_SPACE::GetWallTimeForResidual();
    PHSPACE::TIME_SPACE::GetWallTime();
    PHSPACE::TIME_SPACE::GetCPUTime(0);

    //! Entrance of PHengLEI.
    RunPHengLEI();

    //! Print job over information.
    string infor = "Job Done.\n";
    WriteLogFile(infor);
    PrintToWindow(infor);

    return 0;
}
