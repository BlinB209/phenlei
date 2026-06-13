#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
#          P   P  H   H  E      NN   N  G      L      E       I          +
#          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
#          P      H   H  E      N  N N  G   G  L      E       I          +
#          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
#------------------------------------------------------------------------+
#          Platform for Hybrid Engineering Simulation of Flows           +
#           China Aerodynamics Research and Development Center           +
#                     (C) Copyright, Since 2010                          +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! @file      CMakeLists.txt
#! @brief     It defines the project.
#! @author    Bell, Xu Qingxin, He Kun, He Xin.
cmake_minimum_required(VERSION 3.5)

#set working directory:
set(EXEC_WORKING_DIR PHengLEIv3d0)

#set cmake environment:
set(CMAKE_USE_RELATIVE_PATHS ON)

#output message:
message(STATUS "Welcome to Build Project ${CMAKE_PROJECT_NAME}.")
message(STATUS "The OS is ${CMAKE_SYSTEM_NAME}, version ${CMAKE_SYSTEM_VERSION}.")
message(STATUS "The working directory is ${CMAKE_BINARY_DIR}/${EXEC_WORKING_DIR}/${DLL_WORKING_DIR}.")
message(STATUS "--- Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

# Define some pre-definition.
include(${PROJECT_SOURCE_DIR}/../PHengLEI/PHMacro.cmake)
if(USE_OCC)
  message(STATUS "USE_OCC is turned on!")
  add_definitions(-DUSE_OCC)
endif()
if(USE_DOUBLE_PRECISION)
  message(STATUS "USE_DOUBLE_PRECISION is turned on!")
  add_definitions(-DPH_DOUBLE_PREC)
endif()
if(USE_SMARTARRAY)
  message(STATUS "USE_SMARTARRAY is turned on!")
  add_definitions(-DUSE_SMARTARRAY)
endif()
if(WIN32)
if(USE_Windows_X64)
  message(STATUS "USE_Windows_X64 is turned on!")
  add_definitions(-DUSE_WINDOWS_X64)
endif(USE_Windows_X64)
endif(WIN32)
if(BUILD_SpecDiffHybSolver)
  message(STATUS "BUILD_SpecDiffHybSolver is turned on!")
  add_definitions(-DUSE_SpecDiffHybSolver)
endif()

if(BUILD_IncomSolver)
  message(STATUS "BUILD_IncomSolver is turned on!")
  add_definitions(-DUSE_INCOMSOLVER)
  if(USE_HypreLib)
    add_definitions(-DUSE_HypreLib)
  endif()
  
  if(USE_UnapLib)
    add_definitions(-DUSE_UnapLib)
  endif()
  
  if(USE_YhamgLib)
    add_definitions(-DUSE_YhamgLib)
  endif()
endif()


#1. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory <CFD>, Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

if(WIN32)
if(USE_Windows_X64)
  set(LIB_API_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(CMAKE_EXE_LINKER_FLAGS "/machine:x64")
else(USE_Windows_X64)
  set(LIB_API_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(USE_Windows_X64)
else(WIN32)
  set(LIB_API_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)
link_directories("${LIB_API_EXPORT_PATH}")

set(CFD_TARGET_NAME "PHCFD")
set(LIBRARY_OUTPUT_PATH)
if(WIN32)
if(USE_Windows_X64)
 set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(CMAKE_EXE_LINKER_FLAGS "/machine:x64")
else(USE_Windows_X64)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(USE_Windows_X64)
else(WIN32)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)
message(STATUS "CFD library path is: ${LIBRARY_OUTPUT_PATH}.")

#set options
if (NOT BUILD_LIBRARY_TYPE)
  set (BUILD_LIBRARY_TYPE "Static" CACHE STRING "The type of PHengLEI libraries" FORCE)
  SET_PROPERTY(CACHE BUILD_LIBRARY_TYPE PROPERTY STRINGS Shared Static)
endif()

#3. set environment for MPICH and 3rdParty lib directory:
if(WIN32)
  include_directories("${MSMPI_HOME_DIRECTORY}/Include")
  link_directories("${MSMPI_HOME_DIRECTORY}/lib")
  if(USE_Windows_X64)
    link_directories("${MSMPI_HOME_DIRECTORY}/Lib/x64")
    set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x64")
  else(USE_Windows_X64) 
    link_directories("${MSMPI_HOME_DIRECTORY}/Lib/x86")
    set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x86")
  endif(USE_Windows_X64)
else(WIN32)
  set(CMAKE_CXX_COMPILER  mpicxx)
  set(CMAKE_C_COMPILER    mpicc)
  if(${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR} MATCHES "Linux-x86_64")
    if(USE_INTEL_COMPILER)
      set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/3rdparty/liblinux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/intel")
      message(STATUS "USE_INTEL_COMPILER is turned on!")
    else(USE_INTEL_COMPILER)
      set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/3rdparty/liblinux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/gcc")
      message(STATUS "USE_INTEL_COMPILER is turned off!")
    endif()
  else(${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR} MATCHES "Linux-aarch64")
    set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/3rdparty/liblinux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
  endif()
endif(WIN32)
# specify directory in which the linker will look for libraries for myself libraries:
message(STATUS "The 3rdparty path is: ${THIRDPARTY_LIB_DIR}.") 
link_directories("${THIRDPARTY_LIB_DIR}")

includeHeadFiles()

# select HDF5 H file path with corresponding environment
if(WIN32)
  if(USE_VS2019)
    set(HDF5_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/windows/vs19/win-x64/include)
    add_definitions(-DUSE_VS2019)
  else(USE_VS2019)
    set(HDF5_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/windows/vs12/win-x64/include)
  endif()
else(WIN32)
  if(USE_INTEL_COMPILER)
    set(HDF5_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/linux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/intel/include)
  else(USE_INTEL_COMPILER)
    set(HDF5_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/linux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/gcc/include)
  endif()
endif(WIN32)
include_directories("${HDF5_H_PATH}")
if(BUILD_GPU) 
  add_definitions(-DUSE_CUDA)
  #add_definitions(-DCUDAUNITTEST)
  #add_definitions(-DUNITTESTOUTPUT)

  set(GPU_H_PATH   
  ${PROJECT_SOURCE_DIR}/../GPU/Manage/include
  ${PROJECT_SOURCE_DIR}/../GPU/Kernels/include
  )
  include_directories("${GPU_H_PATH}")
endif(BUILD_GPU)

if(USE_TecplotLib)
set(TECIO_H_PATH   ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/TECIO/include)
include_directories("${TECIO_H_PATH}")
endif(USE_TecplotLib)

set(MEITIS_H_PATH  ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/metis/include)
include_directories("${MEITIS_H_PATH}")
# select CGNS H file path with corresponding environment
if(WIN32)
  if(USE_VS2019)
    set(CGNS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/cgns/windows/vs19/win-x64/include)
  else(USE_VS2019)
    set(CGNS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/cgns/windows/vs12/win-x64/include)
  endif()
else(WIN32)
  if(USE_INTEL_COMPILER)
    set(CGNS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/cgns/linux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/intel/include)
  else(USE_INTEL_COMPILER)
    set(CGNS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/cgns/linux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/gcc/include)
  endif()
endif(WIN32)
include_directories("${CGNS_H_PATH}")
set(MGRID_H_PATH    ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/mgrid/include)
include_directories("${MGRID_H_PATH}")


if(USE_HypreLib)
  set(HYPRE_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/hypre/include)
  include_directories("${HYPRE_H_PATH}")
endif(USE_HypreLib)

if(USE_UnapLib)
  set(UNAP_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/unapLib/include)
  include_directories("${UNAP_H_PATH}")
endif(USE_UnapLib)

if(USE_YhamgLib)
  set(YHAMG_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/yhamg/include)
  include_directories("${YHAMG_H_PATH}")
endif(USE_YhamgLib)


#set(API_H_DIR ${PROJECT_SOURCE_DIR}/3rdparty)

#5. set current directory:
# set directory of Solver
set(SOLVER_H_PATH    Solver/include)
set(SOLVER_H_LIST
${SOLVER_H_PATH}/CFDSolver.h
${SOLVER_H_PATH}/Controller.h
${SOLVER_H_PATH}/Param_CFDSolver.h
)
source_group(Solver\\include FILES ${SOLVER_H_LIST}) 
include_directories(${SOLVER_H_PATH})

set(SOLVER_CPP_PATH      Solver/src)
set(SOLVER_CPP_LIST
${SOLVER_CPP_PATH}/CFDSolver.cpp
${SOLVER_CPP_PATH}/Controller.cpp
${SOLVER_CPP_PATH}/Param_CFDSolver.cpp
)
source_group(Solver\\src FILES ${SOLVER_CPP_LIST})

# set directory of NSSolver:
set(NSSOLVER_H_PATH    NSSolver/include)
set(NSSOLVER_H_LIST
${NSSOLVER_H_PATH}/NSSolver.h
${NSSOLVER_H_PATH}/NSSolverStruct.h
${NSSOLVER_H_PATH}/NSSolverStructFD.h
${NSSOLVER_H_PATH}/NSSolverStructFV.h
${NSSOLVER_H_PATH}/NSSolverUnstruct.h
${NSSOLVER_H_PATH}/MultiGridOperation.h
${NSSOLVER_H_PATH}/Param_NSSolver.h
${NSSOLVER_H_PATH}/Param_NSSolverUnstruct.h
${NSSOLVER_H_PATH}/Param_NSSolverStruct.h
)
source_group(NSSolver\\include FILES ${NSSOLVER_H_LIST}) 
include_directories(${NSSOLVER_H_PATH})

set(NSSOLVER_CPP_PATH      NSSolver/src)
set(NSSOLVER_CPP_LIST
${NSSOLVER_CPP_PATH}/NSSolver.cpp
${NSSOLVER_CPP_PATH}/NSSolverStruct.cpp
${NSSOLVER_CPP_PATH}/NSSolverStructFD.cpp
${NSSOLVER_CPP_PATH}/NSSolverStructFV.cpp
${NSSOLVER_CPP_PATH}/NSSolverUnstruct.cpp
${NSSOLVER_CPP_PATH}/Param_NSSolver.cpp
${NSSOLVER_CPP_PATH}/Param_NSSolverUnstruct.cpp
${NSSOLVER_CPP_PATH}/Param_NSSolverStruct.cpp
)
source_group(NSSolver\\src FILES ${NSSOLVER_CPP_LIST})

# file(GLOB_RECURSE NSSOLVER_H_LIST "NSSolver/*.h")
# source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${NSSOLVER_H_LIST}) 
# FIND_INCLUDE_DIR(NSSOLVER_H_PATH NSSolver)
# include_directories(${NSSOLVER_H_PATH})

# file(GLOB_RECURSE NSSOLVER_CPP_LIST "NSSolver/*.cpp")
# source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${NSSOLVER_CPP_LIST}) 

# set directory of TurbSolver:
set(TurbSOLVER_H_PATH        TurbSolver/include)
set(TurbSOLVER_H_LIST
${TurbSOLVER_H_PATH}/TurbSolver.h
${TurbSOLVER_H_PATH}/TurbSolverStruct.h
${TurbSOLVER_H_PATH}/TransitionSolver.h
${TurbSOLVER_H_PATH}/TransitionSolverStruct.h
${TurbSOLVER_H_PATH}/TurbSolverUnstr.h
${TurbSOLVER_H_PATH}/TransitionSolverUnstr.h
${TurbSOLVER_H_PATH}/Transition.h
${TurbSOLVER_H_PATH}/LESSolver.h
${TurbSOLVER_H_PATH}/LESSolverStruct.h
${TurbSOLVER_H_PATH}/LESSolverUnstr.h
${TurbSOLVER_H_PATH}/TurbSolverStructFD.h
${TurbSOLVER_H_PATH}/Param_TurbSolver.h
${TurbSOLVER_H_PATH}/Param_TransitionSolver.h
${TurbSOLVER_H_PATH}/Param_TurbSolverStruct.h
${TurbSOLVER_H_PATH}/Param_TurbSolverUnstruct.h
${TurbSOLVER_H_PATH}/Param_LESSolver.h
${TurbSOLVER_H_PATH}/Param_LESSolverStruct.h
${TurbSOLVER_H_PATH}/Param_LESSolverUnstruct.h
)
source_group(TurbSolver\\include FILES ${TurbSOLVER_H_LIST}) 
include_directories(${TurbSOLVER_H_PATH})

set(TurbSOLVER_CPP_PATH      TurbSolver/src)
set(TurbSOLVER_CPP_LIST
${TurbSOLVER_CPP_PATH}/TurbSolver.cpp
${TurbSOLVER_CPP_PATH}/TransitionSolver.cpp
${TurbSOLVER_CPP_PATH}/TurbSolverStruct.cpp
${TurbSOLVER_CPP_PATH}/TransitionSolverStruct.cpp
${TurbSOLVER_CPP_PATH}/TurbSolverStructFD.cpp
${TurbSOLVER_CPP_PATH}/TurbSolverUnstr.cpp
${TurbSOLVER_CPP_PATH}/TransitionSolverUnstr.cpp
${TurbSOLVER_CPP_PATH}/Transition.cpp
${TurbSOLVER_CPP_PATH}/LESSolver.cpp
${TurbSOLVER_CPP_PATH}/LESSolverStruct.cpp
${TurbSOLVER_CPP_PATH}/LESSolverUnstr.cpp
${TurbSOLVER_CPP_PATH}/Param_TurbSolver.cpp
${TurbSOLVER_CPP_PATH}/Param_TransitionSolver.cpp
${TurbSOLVER_CPP_PATH}/Param_TurbSolverStruct.cpp
${TurbSOLVER_CPP_PATH}/Param_TurbSolverUnstruct.cpp
${TurbSOLVER_CPP_PATH}/Param_LESSolver.cpp
${TurbSOLVER_CPP_PATH}/Param_LESSolverStruct.cpp
${TurbSOLVER_CPP_PATH}/Param_LESSolverUnstruct.cpp
)
source_group(TurbSolver\\src FILES ${TurbSOLVER_CPP_LIST})

# set directory of HOSolverUnstruct:
set(HOSOLVERUNSTRUCT_H_PATH        HOSolverUnstruct/include)
set(HOSOLVERUNSTRUCT_H_LIST
${HOSOLVERUNSTRUCT_H_PATH}/HOBasisFunction.h
${HOSOLVERUNSTRUCT_H_PATH}/HOBoundary.h
${HOSOLVERUNSTRUCT_H_PATH}/HODefine.h
${HOSOLVERUNSTRUCT_H_PATH}/HOGaussJacobiQuadrature.h
${HOSOLVERUNSTRUCT_H_PATH}/HOGeometryStructure.h
${HOSOLVERUNSTRUCT_H_PATH}/HOInviscidFlux.h
${HOSOLVERUNSTRUCT_H_PATH}/HOShapeFunction.h
${HOSOLVERUNSTRUCT_H_PATH}/HOSimulation.h
${HOSOLVERUNSTRUCT_H_PATH}/HOSolverUnstruct.h
${HOSOLVERUNSTRUCT_H_PATH}/HOSolverUnstructParam.h
${HOSOLVERUNSTRUCT_H_PATH}/HOStandardElement.h
)
source_group(HOSolverUnstruct\\include FILES ${HOSOLVERUNSTRUCT_H_LIST}) 
include_directories(${HOSOLVERUNSTRUCT_H_PATH})

set(HOSOLVERUNSTRUCT_CPP_PATH      HOSolverUnstruct/src)
set(HOSOLVERUNSTRUCT_CPP_LIST
${HOSOLVERUNSTRUCT_CPP_PATH}/HOBasisFunction.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOBoundary.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOGaussJacobiQuadrature.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOGeometryDtructure.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOInviscidFlux.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOShapeFunction.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOSimulation.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOSolverUnstruct.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOSolverUnstructParam.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOStandardElement.cpp
)
source_group(HOSolverUnstruct\\src FILES ${HOSOLVERUNSTRUCT_CPP_LIST})

# Build In-house version:
# set directory of IncomSolver:
if(BUILD_IncomSolver)
set(INCOMSOLVER_H_PATH IncomSolver/include)
set(INCOMSOLVER_H_LIST
${INCOMSOLVER_H_PATH}/MomEqCalculator.h
${INCOMSOLVER_H_PATH}/PPEqCalculator.h
${INCOMSOLVER_H_PATH}/IncomCalculator.h
${INCOMSOLVER_H_PATH}/IncomScalarEqCalculator.h
${INCOMSOLVER_H_PATH}/IncomEnergyEqCalculator.h
${INCOMSOLVER_H_PATH}/IncomKETurbKEqCalculator.h
${INCOMSOLVER_H_PATH}/IncomSATurbKEqCalculator.h
${INCOMSOLVER_H_PATH}/IncomKETurbEpsilonEqCalculator.h
${INCOMSOLVER_H_PATH}/IncomSpeciesEqCalculator.h
${INCOMSOLVER_H_PATH}/IncomUnstructSolver.h
${INCOMSOLVER_H_PATH}/hypreCalculator.h
${INCOMSOLVER_H_PATH}/unapCalculator.h
${INCOMSOLVER_H_PATH}/yhamgCalculator.h
${INCOMSOLVER_H_PATH}/LinearEquationSystemCalculator.h
${INCOMSOLVER_H_PATH}/Param_INCompSolverUnstruct.h
)
source_group(IncomSolver\\include FILES ${INCOMSOLVER_H_LIST})
include_directories(${INCOMSOLVER_H_PATH})

set(INCOMSOLVER_CPP_PATH IncomSolver/src)
set(INCOMSOLVER_CPP_LIST
${INCOMSOLVER_CPP_PATH}/MomEqCalculator.cpp
${INCOMSOLVER_CPP_PATH}/PPEqCalculator.cpp
${INCOMSOLVER_CPP_PATH}/IncomCalculator.cpp
${INCOMSOLVER_CPP_PATH}/IncomScalarEqCalculator.cpp
${INCOMSOLVER_CPP_PATH}/IncomEnergyEqCalculator.cpp
${INCOMSOLVER_CPP_PATH}/IncomKETurbEpsilonEqCalculator.cpp
${INCOMSOLVER_CPP_PATH}/IncomKETurbKEqCalculator.cpp
${INCOMSOLVER_CPP_PATH}/IncomSATurbKEqCalculator.cpp
${INCOMSOLVER_CPP_PATH}/IncomSpeciesEqCalculator.cpp
${INCOMSOLVER_CPP_PATH}/IncomUnstructSolver.cpp
${INCOMSOLVER_CPP_PATH}/hypreCalculator.cpp
${INCOMSOLVER_CPP_PATH}/unapCalculator.cpp
${INCOMSOLVER_CPP_PATH}/yhamgCalculator.cpp
${INCOMSOLVER_CPP_PATH}/LinearEquationSystemCalculator.cpp
${INCOMSOLVER_CPP_PATH}/Param_INCompSolverUnstruct.cpp
)
source_group(IncomSolver\\src FILES ${INCOMSOLVER_CPP_LIST})
endif(BUILD_IncomSolver)

# set directory of SpecDiffHybSolver:
if(BUILD_SpecDiffHybSolver)
set(SpecDiffHybSolver_H_PATH SpecDiffHybSolver/include)
set(SpecDiffHybSolver_H_LIST
${SpecDiffHybSolver_H_PATH}/SpecDiffHybSolver.h
${SpecDiffHybSolver_H_PATH}/SpecGrid.h
${SpecDiffHybSolver_H_PATH}/SpecDiffHybGrid.h
${SpecDiffHybSolver_H_PATH}/ExplicitDifferenceBoundary.h
${SpecDiffHybSolver_H_PATH}/CompactDifferenceFirstDerivative.h
${SpecDiffHybSolver_H_PATH}/CompactDifferenceSecondDerivative.h
${SpecDiffHybSolver_H_PATH}/PoissonSolver.h
${SpecDiffHybSolver_H_PATH}/CorrectWallDivergence.h
${SpecDiffHybSolver_H_PATH}/Complex.h
${SpecDiffHybSolver_H_PATH}/RescaleField.h
${SpecDiffHybSolver_H_PATH}/Statistics.h
${SpecDiffHybSolver_H_PATH}/p3dfft.h
${SpecDiffHybSolver_H_PATH}/fftw3.h
${SpecDiffHybSolver_H_PATH}/SpecDiffHybFieldView.h
${SpecDiffHybSolver_H_PATH}/LES.h
${SpecDiffHybSolver_H_PATH}/SpecSolver.h
${SpecDiffHybSolver_H_PATH}/Spectrum.h
${SpecDiffHybSolver_H_PATH}/Param_SpecSolver.h
${SpecDiffHybSolver_H_PATH}/Param_SpecDiffHybSolver.h
)
source_group(SpecDiffHybSolver\\include FILES ${SpecDiffHybSolver_H_LIST})
include_directories(${SpecDiffHybSolver_H_PATH})

set(SpecDiffHybSolver_CPP_PATH SpecDiffHybSolver/src)
set(SpecDiffHybSolver_CPP_LIST
${SpecDiffHybSolver_CPP_PATH}/SpecDiffHybSolver.cpp
${SpecDiffHybSolver_CPP_PATH}/SpecGrid.cpp
${SpecDiffHybSolver_CPP_PATH}/SpecDiffHybGrid.cpp
${SpecDiffHybSolver_CPP_PATH}/ExplicitDifferenceBoundary.cpp
${SpecDiffHybSolver_CPP_PATH}/CompactDifferenceFirstDerivative.cpp
${SpecDiffHybSolver_CPP_PATH}/CompactDifferenceSecondDerivative.cpp
${SpecDiffHybSolver_CPP_PATH}/PoissonSolver.cpp
${SpecDiffHybSolver_CPP_PATH}/CorrectWallDivergence.cpp
${SpecDiffHybSolver_CPP_PATH}/RescaleField.cpp
${SpecDiffHybSolver_CPP_PATH}/Statistics.cpp
${SpecDiffHybSolver_CPP_PATH}/SpecDiffHybFieldView.cpp
${SpecDiffHybSolver_CPP_PATH}/LES.cpp
${SpecDiffHybSolver_CPP_PATH}/SpecSolver.cpp
${SpecDiffHybSolver_CPP_PATH}/Spectrum.cpp
${SpecDiffHybSolver_CPP_PATH}/Param_SpecSolver.cpp
${SpecDiffHybSolver_CPP_PATH}/Param_SpecDiffHybSolver.cpp
)
source_group(SpecDiffHybSolver\\src FILES ${SpecDiffHybSolver_CPP_LIST})
endif(BUILD_SpecDiffHybSolver)

# set directory of DemoSolver:
if(BUILD_DEMOSOLVER)
set(DEMOSOLVER_H_PATH DemoSolver/include)
set(DEMOSOLVER_H_LIST
${DEMOSOLVER_H_PATH}/DemoSolver.h
${DEMOSOLVER_H_PATH}/Param_DemoSolver.h
)
source_group(DemoSolver\\include FILES ${DEMOSOLVER_H_LIST})
include_directories(${DEMOSOLVER_H_PATH})

set(DEMOSOLVER_CPP_PATH DemoSolver/src)
set(DEMOSOLVER_CPP_LIST
${DEMOSOLVER_CPP_PATH}/DemoSolver.cpp
${DEMOSOLVER_CPP_PATH}/Param_DemoSolver.cpp
)
source_group(DemoSolver\\src FILES ${DEMOSOLVER_CPP_LIST})
endif(BUILD_DEMOSOLVER)

#LBMSolverMPI
if(BUILD_LBMSolverMPI)
set(LBMSolverMPI_H_PATH LBMSolverMPI/include)
set(LBMSolverMPI_H_LIST
${LBMSolverMPI_H_PATH}/LBMSolverMPI.hpp
${LBMSolverMPI_H_PATH}/vector3.hpp
)
source_group(LBMSolverMPI\\include FILES ${LBMSolverMPI_H_LIST})
include_directories(${LBMSolverMPI_H_PATH})

set(LBMSolverMPI_CPP_PATH LBMSolverMPI/src)
set(LBMSolverMPI_CPP_LIST
${LBMSolverMPI_CPP_PATH}/bcs.cpp
${LBMSolverMPI_CPP_PATH}/inout.cpp
${LBMSolverMPI_CPP_PATH}/lbmMPI.cpp
${LBMSolverMPI_CPP_PATH}/lbmMRT.cpp
${LBMSolverMPI_CPP_PATH}/lbmSolver.cpp
${LBMSolverMPI_CPP_PATH}/LBMSolverMPISimulation.cpp
${LBMSolverMPI_CPP_PATH}/vector3.cpp
)
source_group(LBMSolverMPI\\src FILES ${LBMSolverMPI_CPP_LIST})
endif(BUILD_LBMSolverMPI)

#LBMSolverOMP
if(BUILD_LBMSolverOMP)
set(LBMSolverOMP_H_PATH LBMSolverOMP/include)
set(LBMSolverOMP_H_LIST
${LBMSolverOMP_H_PATH}/LBMSolverOMP.hpp
${LBMSolverOMP_H_PATH}/vector3.hpp
)
source_group(LBMSolverOMP\\include FILES ${LBMSolverOMP_H_LIST})
include_directories(${LBMSolverOMP_H_PATH})

set(LBMSolverOMP_CPP_PATH LBMSolverOMP/src)
set(LBMSolverOMP_CPP_LIST
${LBMSolverOMP_CPP_PATH}/bcs.cpp
${LBMSolverOMP_CPP_PATH}/bouzidi.cpp
${LBMSolverOMP_CPP_PATH}/collisionstream.cpp
${LBMSolverOMP_CPP_PATH}/exchange.cpp
${LBMSolverOMP_CPP_PATH}/LBMSolverOMPSimulation.cpp
${LBMSolverOMP_CPP_PATH}/MRT.cpp
${LBMSolverOMP_CPP_PATH}/read.cpp
${LBMSolverOMP_CPP_PATH}/readGeo.cpp
${LBMSolverOMP_CPP_PATH}/vector3.cpp
${LBMSolverOMP_CPP_PATH}/write.cpp
)
source_group(LBMSolverOMP\\src FILES ${LBMSolverOMP_CPP_LIST})

set(IBFSI_H_PATH IBFSI/include)
set(IBFSI_H_LIST
${IBFSI_H_PATH}/IB.hpp
${IBFSI_H_PATH}/vector.hpp
)
source_group(IBFSI\\include FILES ${IBFSI_H_LIST})
include_directories(${IBFSI_H_PATH})

set(IBFSI_CPP_PATH IBFSI/src)
set(IBFSI_CPP_LIST
${IBFSI_CPP_PATH}/IBFSI.cpp
${IBFSI_CPP_PATH}/vector.cpp
${IBFSI_CPP_PATH}/deployPart.cpp
${IBFSI_CPP_PATH}/structrueSolver.cpp
)
source_group(IBFSI\\src FILES ${IBFSI_CPP_LIST})
endif(BUILD_LBMSolverOMP)


#20.1 set header file list
set(CFD_HEAD_LIST 
${SOLVER_H_LIST}
${NSSOLVER_H_LIST}
${TurbSOLVER_H_LIST}
${SpecDiffHybSolver_H_LIST}
${HOSOLVERUNSTRUCT_H_LIST}
)

# add INCOMSOLVER_H_LIST
if(BUILD_IncomSolver)
list(APPEND CFD_HEAD_LIST ${INCOMSOLVER_H_LIST})
endif()

# add SpecDiffHybSolver_H_LIST
if(BUILD_SpecDiffHybSolver)
list(APPEND CFD_HEAD_LIST ${SpecDiffHybSolver_H_LIST})
endif()

# add DEMOSOLVER_H_LIST
if(BUILD_DEMOSOLVER)
list(APPEND CFD_HEAD_LIST ${DEMOSOLVER_H_LIST})
endif()

# add LBMSolverMPI_H_LIST
if(BUILD_LBMSolverMPI)
list(APPEND CFD_HEAD_LIST ${LBMSolverMPI_H_LIST})
endif()

# add LBMSolverOMP_H_LIST
if(BUILD_LBMSolverOMP)
list(APPEND CFD_HEAD_LIST ${LBMSolverOMP_H_LIST})
list(APPEND CFD_HEAD_LIST ${IBFSI_H_LIST})
endif()
#20.2 set source file list:
set(CFD_CPP_LIST 
${SOLVER_CPP_LIST}
${NSSOLVER_CPP_LIST}
${TurbSOLVER_CPP_LIST}
${SpecDiffHybSolver_CPP_LIST}
${HOSOLVERUNSTRUCT_CPP_LIST}
${DEMOSOLVER_CPP_LIST}
)
# add SpecDiffHybSolver_CPP_LIST
if(BUILD_SpecDiffHybSolver)
list(APPEND CFD_CPP_LIST ${SpecDiffHybSolver_CPP_LIST})
endif()

# add DEMOSOLVER_CPP_LIST
if(BUILD_DEMOSOLVER)
list(APPEND CFD_CPP_LIST ${DEMOSOLVER_CPP_LIST})
endif()

# add LBMSolverMPI_CPP_LIST
if(BUILD_LBMSolverMPI)
list(APPEND CFD_CPP_LIST ${LBMSolverMPI_CPP_LIST})
endif()

# add LBMSolverOMP_CPP_LIST
if(BUILD_LBMSolverOMP)
list(APPEND CFD_CPP_LIST ${LBMSolverOMP_CPP_LIST})
list(APPEND CFD_CPP_LIST ${IBFSI_CPP_LIST})
endif()

# add INCOMSOLVER_CPP_LIST
if(BUILD_IncomSolver)
list(APPEND CFD_CPP_LIST ${INCOMSOLVER_CPP_LIST})
endif()

#21. adds flags to the compiler command line for sources in the current directory and below:
add_definitions(-DPH_PARALLEL -DPH_CGNS)

# For OCC.
set(OCC_LIBS)
if(USE_OCC)
  message(STATUS "Using OpenCasCade ...... ")
  set(OCC_LIBS_REQUIRED
      # subset of DataExchange
      TKSTEP TKSTEP209 TKSTEPAttr TKSTEPBase TKIGES TKSTL TKXSBase
      # ModelingAlgorithms
      TKOffset TKFeat TKFillet TKBool TKShHealing TKMesh TKHLR TKBO TKPrim
      TKTopAlgo TKGeomAlgo
      # ModelingData
      TKBRep TKGeomBase TKG3d TKG2d
      # FoundationClasses
      # TKAdvTools -- not necessary? (and removed from OCC 6.8)
      TKMath TKernel)
  
  # Set the path of the OCC library.
  set(3RDPARTY_LIB_PATH ${PROJECT_SOURCE_DIR}/3rdparty/occ730-vc11-bit64-release)
  set(OCC_LIB_PATH ${3RDPARTY_LIB_PATH} CACHE PATH "path to OpenCasCade library.")
  message(STATUS "Path to OpenCasCade library: " ${OCC_LIB_PATH})
  
  # Check if all the library exist -- it's need to be checked. 
  foreach(OCC ${OCC_LIBS_REQUIRED})     
      list(APPEND OCC_LIBS ${OCC})
      unset(OCC CACHE)
  endforeach(OCC)
  
  if(WIN32)           
    MESSAGE(STATUS "Using Win32 OpenCasCade.")
    MESSAGE(STATUS "${OCC_LIB_PATH}/lib")
    MESSAGE(STATUS "${OCC_LIB_PATH}/bin")
    link_directories("${OCC_LIB_PATH}/lib")
    link_directories("${OCC_LIB_PATH}/bin") 
    include_directories("${OCC_LIB_PATH}/inc")
  else(WIN32)
    MESSAGE(STATUS "Using NON-Win32 OpenCasCade.")
    MESSAGE(STATUS "${OCC_LIB_PATH}/lin64/gcc/lib")
    MESSAGE(STATUS "${OCC_LIB_PATH}/inc")
    link_directories("${OCC_LIB_PATH}/lin64/gcc/lib")
    include_directories("${OCC_LIB_PATH}/inc")
  endif(WIN32)
endif()

# Some pre-definitions that occ need.
if(USE_OCC)
  #target_link_libraries(${EXECUTABLE_TARGET_NAME} ${OCC_LIBS})  
  add_definitions(-DWNT)
  add_definitions(-DNOTCL)
endif()

# build executable program or library to the project using the specified source files and header files:
if ("${BUILD_LIBRARY_TYPE}" STREQUAL "Shared")
  # Building Shared Library.   
  message(STATUS "Building Shared Library ${CFD_TARGET_NAME} ...")
  add_library(${CFD_TARGET_NAME} SHARED ${CFD_CPP_LIST} ${CFD_HEAD_LIST})
elseif ("${BUILD_LIBRARY_TYPE}" STREQUAL "Static")
  # Building Static Library.
  message(STATUS "Building Static Library ${CFD_TARGET_NAME} ...")
  add_library(${CFD_TARGET_NAME} STATIC ${CFD_CPP_LIST} ${CFD_HEAD_LIST}) 
endif() 

# link API library:
set(API_TARGET_NAME "PHAPI")

if(WIN32)
  if(USE_VS2019)
    list(APPEND LINK_LIBRARY_LIST ${API_TARGET_NAME} metis parmetis cgns libhdf5 libhdf5_hl)
  else(USE_VS2019)
list(APPEND LINK_LIBRARY_LIST ${API_TARGET_NAME} metis parmetis cgns hdf5)
  endif()
else(WIN32)
  list(APPEND LINK_LIBRARY_LIST ${API_TARGET_NAME} metis parmetis cgns hdf5 hdf5_hl)
endif(WIN32)

if(USE_HypreLib)
  list(APPEND LINK_LIBRARY_LIST hypre)
endif(USE_HypreLib)

if(USE_UnapLib)
  list(APPEND LINK_LIBRARY_LIST unap utilities)
endif(USE_UnapLib)

if(USE_YhamgLib)
  list(APPEND LINK_LIBRARY_LIST yhamg)
endif(USE_YhamgLib)


if(USE_TecplotLib)
  list(APPEND LINK_LIBRARY_LIST tecio)
endif(USE_TecplotLib)

if(WIN32)
  list(APPEND LINK_LIBRARY_LIST ${mpiname})
endif(WIN32)

if(BUILD_SpecDiffHybSolver)
  list(APPEND LINK_LIBRARY_LIST p3dfft fftw3 mpifort gfortran)
endif(BUILD_SpecDiffHybSolver)

if(BUILD_GMRESSolver)
  # sacado library
  list(APPEND LINK_LIBRARY_LIST ${Sacado_LIBRARIES})
  # PETSc library
  list(APPEND LINK_LIBRARY_LIST ${PETSC_Libraries})
endif()

target_link_libraries(${CFD_TARGET_NAME} ${LINK_LIBRARY_LIST})

message(STATUS "Link library list: ${LINK_LIBRARY_LIST}.")


# For OpenMP
if(USE_OMP)
find_package (OpenMP)
if (OpenMP_FOUND)

if(WIN32)
target_include_directories (${CFD_TARGET_NAME}
    PRIVATE
        ${OpenMP_CXX_INCLUDE_DIRS}
)

target_compile_options (${CFD_TARGET_NAME}
    PRIVATE
        ${OpenMP_CXX_FLAGS}
)

target_compile_features (${CFD_TARGET_NAME}
    PUBLIC 
        cxx_std_14
)

target_link_libraries (${CFD_TARGET_NAME}
    PRIVATE 
        ${OpenMP_CXX_LIBRARIES}
)
else(WIN32) 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(WIN32)
endif (OpenMP_FOUND)
endif(USE_OMP)

#24. exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")

