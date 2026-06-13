#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
#          P   P  H   H  E      NN   N  G      L      E       I          +
#          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
#          P      H   H  E      N  N N  G   G  L      E       I          +
#          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
#------------------------------------------------------------------------+
#          Platform for Hybrid Engineering Simulation of Flows           +
#          China Aerodynamics Research and Development Center            +
#                     (C) Copyright, Since 2010                          +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! @file      CMakeLists.txt
#! @brief     It defines the Mesh module, such as parallel mesh generation.
#! @author    Bell.

#. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory <MESH>, Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

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


set(MESH_TARGET_NAME "PHMesh")
set(LIBRARY_OUTPUT_PATH)

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

if(WIN32)
if(USE_Windows_X64)
 set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../Mesh/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(CMAKE_EXE_LINKER_FLAGS "/machine:x64")
else(USE_Windows_X64)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../Mesh/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(USE_Windows_X64)
else(WIN32)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../Mesh/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)
message(STATUS "Mesh library path is: ${LIBRARY_OUTPUT_PATH}.")

#set options
if (NOT BUILD_LIBRARY_TYPE)
  set (BUILD_LIBRARY_TYPE "Static" CACHE STRING "The type of PHengLEI libraries" FORCE)
  SET_PROPERTY(CACHE BUILD_LIBRARY_TYPE PROPERTY STRINGS Shared Static)
endif()

#3. set environment for MPICH:
if(WIN32)
  set(MSMPI_HOME_DIRECTORY "C:/Program Files (x86)/Microsoft SDKs/MPI" CACHE PATH "path to installed MPICH")
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

if(USE_TecplotLib)
set(TECIO_H_PATH  ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/TECIO/include)
include_directories("${TECIO_H_PATH}")
endif(USE_TecplotLib)

set(MEITIS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/metis/include)
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
# set(EIGEN_PATH   ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/eigen)
# include_directories("${EIGEN_PATH}")
set(MGRID_H_PATH   ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/mgrid/include)
include_directories("${MGRID_H_PATH}")

# set current directory:
# set(EIGEN_H_PATH   ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty)
# include_directories(${EIGEN_H_PATH})

# set subdirectory <include>:
set(OVERSET_H_PATH   include)
set(OVERSET_H_LIST
${OVERSET_H_PATH}/UnstructuredOversetConfig.h
${OVERSET_H_PATH}/OversetKDTree.h
${OVERSET_H_PATH}/DynamicGrid.h
${OVERSET_H_PATH}/Mesh_Deformation.h
${OVERSET_H_PATH}/Mesh_DeformationRBF.h
${OVERSET_H_PATH}/Mesh_DeformationSPRING.h
${OVERSET_H_PATH}/ComputationalGrid.h
${OVERSET_H_PATH}/InterfaceLinkStructure.h
${OVERSET_H_PATH}/Mesh_CellTopo.h
${OVERSET_H_PATH}/Mesh_FaceTopo.h
${OVERSET_H_PATH}/Mesh_Refine.h
${OVERSET_H_PATH}/Mesh_Refine_Struct.h
${OVERSET_H_PATH}/Mesh_Refine_Unstruct.h
${OVERSET_H_PATH}/Mesh_RefineController.h
${OVERSET_H_PATH}/Mesh_ConstructInterfaceLink.h
)
source_group(include FILES ${OVERSET_H_LIST}) 
include_directories(${OVERSET_H_PATH})


# set directory <src>:
set(OVERSET_CPP_PATH  src)
set(OVERSET_CPP_LIST
${OVERSET_CPP_PATH}/UnstructuredOversetConfig.cpp
${OVERSET_CPP_PATH}/OversetKDTree.cpp
${OVERSET_CPP_PATH}/DynamicGrid.cpp
${OVERSET_CPP_PATH}/Mesh_Deformation.cpp
${OVERSET_CPP_PATH}/Mesh_DeformationRBF.cpp
${OVERSET_CPP_PATH}/Mesh_DeformationSPRING.cpp
${OVERSET_CPP_PATH}/ComputationalGrid.cpp
${OVERSET_CPP_PATH}/InterfaceLinkStructure.cpp
${OVERSET_CPP_PATH}/Mesh_CellTopo.cpp
${OVERSET_CPP_PATH}/Mesh_FaceTopo.cpp
${OVERSET_CPP_PATH}/Mesh_Refine.cpp
${OVERSET_CPP_PATH}/Mesh_Refine_Struct.cpp
${OVERSET_CPP_PATH}/Mesh_Refine_Unstruct.cpp
${OVERSET_CPP_PATH}/Mesh_RefineController.cpp
${OVERSET_CPP_PATH}/Mesh_ConstructInterfaceLink.cpp
)
source_group(src FILES ${OVERSET_CPP_LIST}) 
include_directories(${OVERSET_CPP_PATH})


# set header file list:
set(HEAD_LIST
${OVERSET_H_LIST}
)
# set source file list:
set(CPP_LIST 
${OVERSET_CPP_LIST}
)
# adds flags to the compiler command line for sources in the current directory and below:
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
  set(3RDPARTY_LIB_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/occ730-vc11-bit64-release)
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

# build executable program or library to the project using the specified source files and header files:
if ("${BUILD_LIBRARY_TYPE}" STREQUAL "Shared")
  # Building Shared Library.   
  message(STATUS "Building Shared Library ${MESH_TARGET_NAME} ...")
  add_library(${MESH_TARGET_NAME} SHARED ${CPP_LIST} ${HEAD_LIST})
elseif ("${BUILD_LIBRARY_TYPE}" STREQUAL "Static")
  # Building Static Library.
  message(STATUS "Building Static Library ${MESH_TARGET_NAME} ...")
  add_library(${MESH_TARGET_NAME} STATIC ${CPP_LIST} ${HEAD_LIST}) 
endif() 

# link API library:
set(API_TARGET_NAME "PHAPI")
list(APPEND LINK_LIBRARY_LIST ${API_TARGET_NAME})
target_link_libraries(${MESH_TARGET_NAME} ${LINK_LIBRARY_LIST})

# For OpenMP
if(USE_OMP)
find_package (OpenMP)
if (OpenMP_FOUND)

if(WIN32)
target_include_directories (${MESH_TARGET_NAME}
  PRIVATE
    ${OpenMP_CXX_INCLUDE_DIRS}
)

target_compile_options (${MESH_TARGET_NAME}
  PRIVATE
    ${OpenMP_CXX_FLAGS}
)

target_compile_features (${MESH_TARGET_NAME}
  PUBLIC
    cxx_std_14
)

target_link_libraries (${MESH_TARGET_NAME}
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

# exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")

