cmake_minimum_required(VERSION 3.5)

#build the project:
project(PHengLEIv3d0 CXX)

if(WIN32)
  add_compile_options("/utf-8")
endif()

#set working directory:
set(EXEC_WORKING_DIR PHengLEIv3d0)

#set COMPILER:
if(WIN32)
  if(CMAKE_GENERATOR MATCHES "Visual Studio 11 2012")
    SET(THIRD_PARTY_COMPILER_FLAG "libwindows/libvs12")
    SET(IS_VS2019 OFF)
    add_definitions(-DVS12)
  elseif(CMAKE_GENERATOR MATCHES "Visual Studio 12 2013")
    SET(THIRD_PARTY_COMPILER_FLAG "libwindows/libvs12")
    SET(IS_VS2019 OFF)
    add_definitions(-DVS13)
  elseif(CMAKE_GENERATOR MATCHES "Visual Studio 14 2015")
    SET(THIRD_PARTY_COMPILER_FLAG "libwindows/libvs19")
    SET(IS_VS2019 ON)
  elseif(CMAKE_GENERATOR MATCHES "Visual Studio 15 2017")
    SET(THIRD_PARTY_COMPILER_FLAG "libwindows/libvs19")
    SET(IS_VS2019 ON)
  elseif(CMAKE_GENERATOR MATCHES "Visual Studio 16 2019")
    SET(THIRD_PARTY_COMPILER_FLAG "libwindows/libvs19")
    SET(IS_VS2019 ON)
  elseif(CMAKE_GENERATOR MATCHES "Visual Studio 17 2022")
    SET(THIRD_PARTY_COMPILER_FLAG "libwindows/libvs19")
    SET(IS_VS2019 ON)
  elseif(CMAKE_GENERATOR MATCHES "Visual Studio 18 2026")
    SET(THIRD_PARTY_COMPILER_FLAG "libwindows/libvs19")
    SET(IS_VS2019 ON)
  endif()
else(WIN32)
  SET(THIRD_PARTY_COMPILER_FLAG "Unix")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ldl")
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ldl")
  endif()
endif()

include(${PROJECT_SOURCE_DIR}/../PHengLEI/PHMacro.cmake)

#set cmake environment:
set(CMAKE_USE_RELATIVE_PATHS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "CHOOSE" FORCE)
endif()
message(STATUS "The current CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}.")

#output message:
message(STATUS "Welcome to Build Project ${CMAKE_PROJECT_NAME}.")
message(STATUS "The OS is ${CMAKE_SYSTEM_ID}, version ${CMAKE_SYSTEM_VERSION}.")
message(STATUS "The COMPILER is ${THIRD_PARTY_COMPILER_FLAG}.")
message(STATUS "The working directory is ${CMAKE_BINARY_DIR}/${EXEC_WORKING_DIR}/${DLL_WORKING_DIR}.")
message(STATUS "--- Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "The current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "The current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

# Define some pre-definition.
opt_build(API "Build PHengLEI framework or not" ON) #Build API or not.
opt_build(CFD "Build CFD module or not" ON) #Build CFD or not.
opt_build(MESH "Build Mesh module or not" ON) #Build Mesh or not.
opt_build(DemoSolver "Build Solver Template or not" OFF)# Build Solver Template option.
opt_build(IncomSolver "Build Pressure Based Solver or not" ON)# Build IncomSolver option.
if(BUILD_IncomSolver)
  opt(HypreLib "Build Pressure Based Solver or not" ON)# Build IncomSolver option.
  opt(UnapLib "Build Pressure Based Solver or not" OFF)# Build IncomSolver option.
  opt(YhamgLib "Build Pressure Based Solver or not" OFF)# Build IncomSolver option.
endif()
opt_build(GMRESSolver "Build GMRES Solver or not" OFF)# Build GMRESSolver option.
opt_build(SpecDiffHybSolver "Build Spectra-Difference Hybrid Solver or not" OFF) #Build Spectra-Difference Hybrid Solver or not.
opt(DOUBLE_PRECISION "Use double precision or not" ON)
opt(ERROR_DATA_DEBUG "Build error data debug or not" OFF)
opt(OCC "Use OpenCasCade Library" OFF) #Use occ or not.
opt(OMP "Use OpenMP or not" OFF) #Use OPenMP or not.
opt(SMARTARRAY "Use smart array more efficenet than origial array" ON) #Use smart array or not.
opt(TecplotLib "Use tecplot lib to dump binary flow file" ON) #Use tecplotlib or not.
opt_build(LBMSolverMPI "Build LBM-MPI module or not" OFF) #Build LBM-MPI or not.
opt_build(LBMSolverOMP "Build LBM-OMP module or not" OFF) #Build LBM-OMP or not.
opt_build(GPU "Build GPU module or not" OFF) #Build GPU or not.
opt(Windows_X64 "Compile PHengLEI X64 on Windows platform" ON) ###
opt(AI_RandomForestRegressor "Compile data generation interface for AI predictions of RandomForestRegressor " OFF) ###
if(WIN32)
  opt(VS2019 "Compile PHengLEI X64 with Microsoft Visual Studio 16 2019 generator" ${IS_VS2019}) # Choose VS2019 or VS2012 to build this project.
  opt(INTELMPI "Use INTELMPI module or not" OFF) #Use INTELMPI or not.
  if(USE_INTELMPI)
    set(mpiname impi)
    set(MSMPI_HOME_DIRECTORY "C:/Program Files (x86)/Intel/IMPI/intel64")
  else(USE_INTELMPI)
    set(mpiname msmpi)
    set(MSMPI_HOME_DIRECTORY "C:/Program Files (x86)/Microsoft SDKs/MPI")
  endif()
  message(STATUS "\n\nThe current MSMPI_HOME_DIRECTORY is ${MSMPI_HOME_DIRECTORY}.")
endif(WIN32)

set(GPROF OFF)
if(GPROF)
  set(CMAKE_CXX_FLAGS "$ENV{CMAKE_CXX_FLAGS} -pg -ldl")
  set(CMAKE_C_FLAGS "$ENV{CMAKE_C_FLAGS} -pg")
  set(CMAKE_EXE_LINKER_FLAGS "$ENV{CMAKE_EXE_LINKER_FLAGS} -pg")
  set(CMAKE_SHARED_LINKER_FLAGS "$ENV{CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

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

if(BUILD_SpecDiffHybSolver)
  message(STATUS "BUILD_SpecDiffHybSolver is turned on!")
  add_definitions(-DUSE_SpecDiffHybSolver)
endif()

if(BUILD_DemoSolver)
  message(STATUS "BUILD_DemoSolver is turned on!")
  add_definitions(-DUSE_DEMOSOLVER)
endif()

if(BUILD_GMRESSolver)
  message(STATUS "BUILD_GMRESSolver is turned on!")
  add_definitions(-DUSE_GMRESSOLVER)

  #For Sacado
  SET(SACADO_DIR ${PROJECT_SOURCE_DIR}/3rdparty/Sacado)
  FIND_PATH(Sacado_INCLUDE_DIR Sacado.hpp
  HINTS ${SACADO_DIR}
  PATH_SUFFIXES include
  )
  FIND_LIBRARY(Sacado_LIBRARY NAMES sacado
  HINTS ${SACADO_DIR}
  PATH_SUFFIXES lib
  )
  set(Sacado_LIBRARIES ${Sacado_LIBRARY})
  include_directories("${Sacado_INCLUDE_DIR}")

  #For petsc
  set(PETSC ${PROJECT_SOURCE_DIR}/3rdparty/PETSc/)
  set(PETSC_Include_Dir ${PROJECT_SOURCE_DIR}/3rdparty/PETSc/include)
  include_directories("${PETSC_Include_Dir}")

  FIND_LIBRARY(PETSC_Lib NAMES petsc
    HINTS ${PETSC}
    PATH_SUFFIXES lib
 )

 FIND_LIBRARY(F2cblas_Lib NAMES f2cblas
    HINTS ${PETSC}
    PATH_SUFFIXES lib
 )
 FIND_LIBRARY(F2clapack_Lib NAMES f2clapack
    HINTS ${PETSC}
    PATH_SUFFIXES lib
 )
 set(PETSC_Libraries ${PETSC_Lib} ${F2cblas_Lib} ${F2clapack_Lib})

endif()

if(USE_OMP)
  message(STATUS "USE_OMP is turned on!")
  add_definitions(-DUSE_OMP)
  find_package(OpenMP)
  if(OpenMP_FOUND)
    message(STATUS "OpenMP_VERSION = ${OpenMP_VERSION}")
    message(STATUS "OpenMP_FOUND   = ${OpenMP_FOUND}")
    message(STATUS "OpenMP_CXX_FOUND = ${OpenMP_CXX_FOUND}")
    message(STATUS "OpenMP_CXX_FLAGS = ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP_CXX_INCLUDE_DIRS = ${OpenMP_CXX_INCLUDE_DIRS}")
    message(STATUS "OpenMP_CXX_LIB_NAMES = ${OpenMP_CXX_LIB_NAMES}")
    message(STATUS "OpenMP_CXX_LIBRARIES = ${OpenMP_CXX_LIBRARIES}")
    message(STATUS "OpenMP_LIBRARY = ${OpenMP_LIBRARY}")
    message(STATUS "OpenMP_CXX_SPEC_DATE = ${OpenMP_CXX_SPEC_DATE}")
    message(STATUS "OpenMP_CXX_VERSION_MAJOR = ${OpenMP_CXX_VERSION_MAJOR}")
    message(STATUS "OpenMP_CXX_VERSION_MINOR = ${OpenMP_CXX_VERSION_MINOR}")
    message(STATUS "OpenMP_CXX_VERSION = ${OpenMP_CXX_VERSION}") 
  endif(OpenMP_FOUND)
endif(USE_OMP)

find_package(MPI REQUIRED)
include_directories ("${MPI_CXX_INCLUDE_DIRS}")
message (STATUS "MPI_CXX_INCLUDE_DIRS = ${MPI_CXX_INCLUDE_DIRS}") 
message (STATUS "MPI_CXX_LIBRARIES = ${MPI_CXX_LIBRARIES}") 

if(USE_ERROR_DATA_DEBUG)
  message(STATUS "USE_ERROR_DATA_DEBUG is turned on!")
  add_definitions(-DUSE_ERROR_DATA_DEBUG)
endif(USE_ERROR_DATA_DEBUG)

if(BUILD_IncomSolver)
  message(STATUS "BUILD_IncomSolver is turned on!")
  add_definitions(-DUSE_INCOMSOLVER)
  if(USE_HypreLib)
    add_definitions(-DUSE_HypreLib)
  endif(USE_HypreLib)
  
  if(USE_UnapLib)
    add_definitions(-DUSE_UnapLib)
  endif(USE_UnapLib)
  
  if(USE_YhamgLib)
    add_definitions(-DUSE_YhamgLib)
  endif(USE_YhamgLib)
endif(BUILD_IncomSolver)

if(BUILD_LBM)
  message(STATUS "BUILD_LBM is turned on!")
  add_definitions(-DUSE_LBM_CODE)
  add_compile_options(-std=c++14)
endif(BUILD_LBM)

if(USE_TecplotLib)
  message(STATUS "USE_TecplotLib is turned on!")
  add_definitions(-DUSE_TecplotLib)
endif(USE_TecplotLib)

if(WIN32)
  if(USE_Windows_X64)
    message(STATUS "USE_Windows_X64 is turned on!")
    add_definitions(-DUSE_WINDOWS_X64)
  endif(USE_Windows_X64)
endif(WIN32)

if(WIN32)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
endif()

if(BUILD_LBMSolverMPI)
  message(STATUS "BUILD_LBM-MPI is turned on!")
  add_definitions(-DUSE_LBMSolverMPI)
endif(BUILD_LBMSolverMPI)

if(BUILD_LBMSolverOMP)
  message(STATUS "BUILD_LBM-OMP is turned on!")
  add_definitions(-DUSE_LBMSolverOMP)
endif(BUILD_LBMSolverOMP)
# Add subdirectory or link library PHAPI.
set(API_TARGET_NAME "PHAPI")

if(WIN32)
  if(USE_AI_RandomForestRegressor)
    add_definitions(-DAI_RandomForestRegressor)
  endif(USE_AI_RandomForestRegressor)
endif()

includeHeadFiles()

# select HDF5 H file path with corresponding environment
if(WIN32)
  if(IS_VS2019)
    set(HDF5_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/windows/vs19/win-x64/include)
    add_definitions(-DUSE_VS2019)
  else(IS_VS2019)
    set(HDF5_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/windows/vs12/win-x64/include)
  endif()
else(WIN32)
  if(USE_INTEL_COMPILER)
    set(HDF5_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/linux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/intel/include)
  else(USE_INTEL_COMPILER)
    set(HDF5_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/linux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/gcc/include)
  endif()
  add_definitions(-DUSE_VS2019)
endif(WIN32)
include_directories("${HDF5_H_PATH}")
set(GMSH_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/Gmsh/include)
include_directories("${GMSH_H_PATH}")
set(PARMETIS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/parmetis/include)
include_directories("${PARMETIS_H_PATH}")
if(USE_TecplotLib)
  set(TECIO_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/TECIO/include)
  include_directories("${TECIO_H_PATH}")
endif(USE_TecplotLib)

set(MEITIS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/metis/include)
include_directories("${MEITIS_H_PATH}")
# select CGNS H file path with corresponding environment
if(WIN32)
  if(IS_VS2019)
    set(CGNS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/cgns/windows/vs19/win-x64/include)
  else(IS_VS2019)
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
set(MGRID_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/mgrid/include)
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

# Add subdirectory or link library PHMesh.
set(MESH_TARGET_NAME "PHMesh")
if(WIN32)
  if(USE_Windows_X64)
    set(LIB_MESH_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../Mesh/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  else(USE_Windows_X64)
    set(LIB_MESH_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../Mesh/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  endif(USE_Windows_X64)
else(WIN32)
  set(LIB_MESH_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../Mesh/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)

set(MESH_H_PATH ${PROJECT_SOURCE_DIR}/../Mesh/include)
include_directories(${MESH_H_PATH})
if(BUILD_MESH)
  message(STATUS "Build MESH!")
  add_subdirectory(../Mesh ${EXEC_WORKING_DIR}/Mesh)
  if(WIN32)
    set_target_properties(${Mesh_TARGET_NAME} PROPERTIES COMPILE_FLAGS "/MP")
  endif(WIN32)
else(BUILD_MESH)
  message(STATUS "Link MESH!") 
  link_directories("${LIB_MESH_EXPORT_PATH}")  
endif(BUILD_MESH)
# Add subdirectory or link library PHCFD.
set(CFD_TARGET_NAME "PHCFD")
if(WIN32)
    if(USE_Windows_X64)
      set(LIB_Grid_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../Grider/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
    else(USE_Windows_X64)
      set(LIB_Grid_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../Grider/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
    endif(USE_Windows_X64)
else(WIN32)
  set(LIB_CFD_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)

set(CFD_H_PATH ${PROJECT_SOURCE_DIR}/../CFD)
# FIND_INCLUDE_DIR(NSSOLVER_H_PATH ${CFD_H_PATH}/NSSolver)
# include_directories(${NSSOLVER_H_PATH})

include_directories(${CFD_H_PATH}/NSSolver/include)
include_directories(${CFD_H_PATH}/Solver/include)
include_directories(${CFD_H_PATH}/TurbSolver/include)
include_directories(${CFD_H_PATH}/PostProcess/include)
include_directories(${CFD_H_PATH}/HOSolverUnstruct/include)
if(BUILD_LBMSolverMPI)
include_directories(${CFD_H_PATH}/LBMSolverMPI/include)
endif(BUILD_LBMSolverMPI)
if(BUILD_LBMSolverOMP)
include_directories(${CFD_H_PATH}/LBMSolverOMP/include)
include_directories(${CFD_H_PATH}/IBFSI/include)
endif(BUILD_LBMSolverOMP)

if(BUILD_SpecDiffHybSolver)
  include_directories(${CFD_H_PATH}/SpecDiffHybSolver/include)
endif(BUILD_SpecDiffHybSolver)
if(BUILD_IncomSolver)
include_directories(${CFD_H_PATH}/IncomSolver/include)
endif(BUILD_IncomSolver)
if(BUILD_DemoSolver)
  include_directories(${CFD_H_PATH}/DemoSolver/include)
endif(BUILD_DemoSolver)

# Add subdirectory or link library PHGPU.
if(BUILD_GPU)
  SET(CUDA_SEPARABLE_COMPILATION ON)
  set(GPU_TARGET_NAME "PHGPU")
  if(WIN32)
    if(USE_Windows_X64)
      set(LIB_GPU_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../GPU/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
    else(USE_Windows_X64)
      set(LIB_GPU_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../GPU/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
    endif(USE_Windows_X64)
  else(WIN32)
    set(LIB_GPU_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../GPU/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  endif(WIN32)

  set(GPU_H_PATH   ${PROJECT_SOURCE_DIR}/../GPU)
  include_directories(${GPU_H_PATH}/Manage/include)
  include_directories(${GPU_H_PATH}/Kernels/include)
endif(BUILD_GPU)

if(BUILD_API)
  message(STATUS "Build API!")
  add_subdirectory(../API ${EXEC_WORKING_DIR}/API)
  if(WIN32)
    set_target_properties(${API_TARGET_NAME} PROPERTIES COMPILE_FLAGS "/MP /W4")
    ADD_COMPILE_OPTIONS(/wd26451)
  endif(WIN32)
else(BUILD_API)
  message(STATUS "Link API!") 
  link_directories("${LIB_API_EXPORT_PATH}")  
endif(BUILD_API)

if(BUILD_CFD)
  message(STATUS "Build CFD!")
  add_subdirectory(../CFD ${EXEC_WORKING_DIR}/CFD)
  if(WIN32)
    set_target_properties(${CFD_TARGET_NAME} PROPERTIES COMPILE_FLAGS "/MP /W4")
    ADD_COMPILE_OPTIONS(/wd26451)
  endif(WIN32)
else(BUILD_CFD)
  message(STATUS "Link CFD!") 
  link_directories("${LIB_CFD_EXPORT_PATH}")  
endif(BUILD_CFD)

if(BUILD_GPU)
  message(STATUS "Build GPU!")
  add_definitions(-DUSE_CUDA)
  add_subdirectory(../GPU ${EXEC_WORKING_DIR}/GPU)
  set(GPUMODE ON)
  if(WIN32)
    set_target_properties(${GPU_TARGET_NAME} PROPERTIES COMPILE_FLAGS "/MP")
  endif(WIN32)
endif(BUILD_GPU)

# define compile macro for nvcc, open only one mode
set(UNITTESTMODE OFF)

if (GPUMODE)
    message(STATUS "GPUMODE is turned on!")
    #define architecture
    set(V100MODE ON) #for Nvidia Tesla V100, else for Nvidia Tesla K80
    set(ORGMPIMODE ON)
    set(NEWMPIMODE OFF)
    set(TESTMPIMODE OFF)
    add_definitions(-DUSE_CUDA)
    add_definitions(-DLOOPMODEOPT) #loop mode optimization for NSSolverUnstruct::ComputeNodeValue and Limiter::MinMaxDiff
    add_definitions(-DNESTLOOPOPT) # split nested loop into separate loops
  add_definitions(-DCUDADEBUGATOMIC)
    add_definitions(-DCUDACOPYBACK)

    #define compile macro for multi-color LU-SGS
    set(MCLUSGS OFF)
    if (MCLUSGS)
        add_definitions(-DMCLUSGS) #Just for multi-color LU-SGS, only works at tscheme=4
        add_definitions(-DMCLUSEPSGS) #L-SGS and U-SGS are separated
    else(MCLUSGS) #MCLUSGS is OFF
        add_definitions(-DGPUCPURUN) #only works in LU-SGS scheme, right hand side is calculated by GPU. LU-SGS part is calculated by CPU.
    endif(MCLUSGS)

    #define compile macro for gprof
    set(GPROF OFF)
    if (GPROF)
        set(CMAKE_CXX_FLAGS "$ENV{CMAKE_CXX_FLAGS} -pg")
        set(CMAKE_C_FLAGS "$ENV{CMAKE_C_FLAGS} -pg")
        set(CMAKE_EXE_LINKER_FLAGS "$ENV{CMAKE_EXE_LINKER_FLAGS} -pg")
        set(CMAKE_SHARED_LINKER_FLAGS "$ENV{CMAKE_SHARED_LINKER_FLAGS} -pg")
    endif()

    if (ORGMPIMODE)
        add_definitions(-DORGMPI) #for large MPI buffer
        add_definitions(-DMPIBUFFERLARGE) #for large MPI buffer
        add_definitions(-DREPLACEBUFFER) #Use only with MPIBUFFERLARGE
  endif()

    if (NEWMPIMODE)
        add_definitions(-DNEWMPI)
        add_definitions(-DMPIBUFFERLARGE) #for large Host MPI buffer
        add_definitions(-DISENDIRECVMPI) #MPI_Isend and MPI_Irecv are used. Otherwise, MPI_Sendrecv is used.
    endif()

    if (TESTMPIMODE)
        add_definitions(-DORGMPI) #for large MPI buffer
        add_definitions(-DNEWMPI)
        add_definitions(-DMPIBUFFERLARGE) #for large MPI buffer
        add_definitions(-DGPUBUFFERSMALL) #for small GPU buffer, every variable owns ngbZone's buffer. Only used with ISENDIRECVMPI and MPIOVERLAP. In the future, CUDA_AWARE_MPI may be considered.
        add_definitions(-DISENDIRECVMPI) #MPI_Isend and MPI_Irecv are used. Otherwise, MPI_Sendrecv is used.
        add_definitions(-DMPIOVERLAP) #MPIOVERLAP mode cannot be test correctly, because test function consumes a lot time that MPI message may be passed.
        add_definitions(-DHDCOVERLAP) # HtoD, DtoH and Computing can be overlapped together.
        add_definitions(-DTESTNEWORGMPI) #Test CUDA_AWARE_MPI.
    endif()
endif()

if (UNITTESTMODE)
  message(STATUS "UNITTESTMODE is turned on!")
    add_definitions(-DCPURUN)
    add_definitions(-DUSE_CUDA)
    set(V100MODE ON) #for Nvidia Tesla V100, else for Nvidia Tesla K80
    add_definitions(-DCUDAUNITTEST)
    add_definitions(-DDEBUGPOW) #FOR DEBUG, cancel all of pow operations, even cpu code is modified.
    add_definitions(-DUNITTESTOUTPUT) #Most are the successfull information in unit test.
    add_definitions(-DDEBUGSolution) #SolutionFix makes GPUUpdateFlowFieldM1_S1 execute series.
  set(ORGMPIMODE OFF)
    if (ORGMPIMODE)
        add_definitions(-DORGMPI) #for large MPI buffer
        add_definitions(-DMPIBUFFERLARGE) #for large MPI buffer
        add_definitions(-DREPLACEBUFFER) #Use only with MPIBUFFERLARGE
  endif()

    #define compile macro for multi-color LU-SGS
    set(MCLUSGS OFF)
    if (MCLUSGS)
        add_definitions(-DMCLUSGS) #Just for multi-color LU-SGS, only works at tscheme=4
        add_definitions(-DMCLUSEPSGS) #L-SGS and U-SGS are separated
    else(MCLUSGS) #MCLUSGS is OFF
        add_definitions(-DGPUCPURUN) #only works in LU-SGS scheme, right hand side is calculated by GPU. LU-SGS part is calculated by CPU.
    endif(MCLUSGS)
endif()

if(WIN32)
  if(USE_Windows_X64)
    link_directories("${PROJECT_SOURCE_DIR}/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x64")
  else(USE_Windows_X64)
    link_directories("${PROJECT_SOURCE_DIR}/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x86")
  endif(USE_Windows_X64)
  message(STATUS "CMAKE_CXX_COMPILER! ${CMAKE_CXX_COMPILER}")
else(WIN32)
  if(${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR} MATCHES "Linux-x86_64")
    opt(INTEL_COMPILER "Use intel compiler or not" OFF) #Use INTEL COMPILER or not.
    if(USE_INTEL_COMPILER)
      link_directories("${PROJECT_SOURCE_DIR}/3rdparty/liblinux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/intel")
      message(STATUS "USE_INTEL_COMPILER is turned on!")
    else(USE_INTEL_COMPILER)
      link_directories("${PROJECT_SOURCE_DIR}/3rdparty/liblinux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/gcc")
      message(STATUS "USE_INTEL_COMPILER is turned off!")
    endif()
  else(${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR} MATCHES "Linux-aarch64")
    link_directories("${PROJECT_SOURCE_DIR}/3rdparty/liblinux/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
  endif()
endif(WIN32)

# look for CMakeLists.txt in subdirectory <codes> and build an executable target:
add_subdirectory(../Main ${EXEC_WORKING_DIR})


#exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
message(STATUS "Finished build project ${PROJECT_NAME}! ")

