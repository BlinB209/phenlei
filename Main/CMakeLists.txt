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
#! @brief     It defines the project.
#! @author    Bell, Xu Qingxin, He Kun, He Xin.
cmake_minimum_required(VERSION 3.5)

#1. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory <Main>, Start to call ${CMAKE_CURRENT_LIST_FILE} ---")
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

#2. set executable name:
set(EXECUTABLE_TARGET_NAME "${CMAKE_PROJECT_NAME}")
message(STATUS "The executable name is ${EXECUTABLE_TARGET_NAME}.")

# Define some pre-definition.
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

#5. set current directory:

#6. set subdirectory <include>:
set(MAIN_H_PATH   include)
set(MAIN_H_LIST
${MAIN_H_PATH}/CFDSimulation.h
${MAIN_H_PATH}/Simulation.h
${MAIN_H_PATH}/GridFactory.h
${MAIN_H_PATH}/FieldViewIn.h
${MAIN_H_PATH}/GridPartition.h
${MAIN_H_PATH}/GridView.h
${MAIN_H_PATH}/FluentIn.h
${MAIN_H_PATH}/Gmesh.h
${MAIN_H_PATH}/GridgenIn.h
${MAIN_H_PATH}/SolversList.h
)
source_group(Main\\include FILES ${MAIN_H_LIST}) 
include_directories(${MAIN_H_PATH})

#7. set directory <src>:
set(MAIN_CPP_PATH  src)
set(MAIN_CPP_LIST
${MAIN_CPP_PATH}/main.cpp
${MAIN_CPP_PATH}/CFDSimulation.cpp
${MAIN_CPP_PATH}/Simulation.cpp
${MAIN_CPP_PATH}/GridFactory.cpp
${MAIN_CPP_PATH}/FieldViewIn.cpp
${MAIN_CPP_PATH}/GridPartition.cpp
${MAIN_CPP_PATH}/GridView.cpp
${MAIN_CPP_PATH}/FluentIn.cpp
${MAIN_CPP_PATH}/Gmesh.cpp
${MAIN_CPP_PATH}/GridgenIn.cpp
)

source_group(Main\\src FILES ${MAIN_CPP_LIST}) 
include_directories(${MAIN_CPP_PATH})

#20.1 set header file list
set(HEAD_LIST 
${MAIN_H_LIST}
)

#20.2 set source file list:
set(CPP_LIST 
${MAIN_CPP_LIST}
)

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
endif(USE_OCC)

#22. build executable program:

# Building Executable Program.
message(STATUS "Building Executable Program ${EXECUTABLE_TARGET_NAME} ...")
add_executable(${EXECUTABLE_TARGET_NAME} ${CPP_LIST} ${HEAD_LIST})  
 if(WIN32) 
    set_target_properties(${EXECUTABLE_TARGET_NAME} PROPERTIES COMPILE_FLAGS "/MP")
 endif(WIN32)
# Collect the linraries to be linked.
# link the API/MESH etc. libraries, no matter build them or not.
# Important: the top library should be list firstly!!!
# For example, the Mesh is previous API.

set(API_TARGET_NAME "PHAPI")
set(MESH_TARGET_NAME "PHMesh")
set(CFD_TARGET_NAME "PHCFD")
set(Grider_TARGET_NAME "PHGrider")

message(STATUS "Linking ${MESH_TARGET_NAME} ...")
list(APPEND LINK_LIBRARY_LIST ${MESH_TARGET_NAME})

message(STATUS "Linking ${CFD_TARGET_NAME} ...")
list(APPEND LINK_LIBRARY_LIST ${CFD_TARGET_NAME})

if(BUILD_GPU) 
  add_definitions(-DUSE_CUDA)
  #add_definitions(-DCUDAUNITTEST)
  #add_definitions(-DUNITTESTOUTPUT) 
  message(STATUS "Linking ${GPU_TARGET_NAME} ...")
  list(APPEND LINK_LIBRARY_LIST ${GPU_TARGET_NAME})
endif(BUILD_GPU)

list(APPEND LINK_LIBRARY_LIST ${API_TARGET_NAME})

message(STATUS "Link library list: ${LINK_LIBRARY_LIST}.")

if(BUILD_GMRESSolver)
  #for Sacado
  list(APPEND LINK_LIBRARY_LIST ${Sacado_LIBRARIES})
  #for PETSc
  list(APPEND LINK_LIBRARY_LIST ${PETSC_Libraries})

  list(APPEND LINK_LIBRARY_LIST ${CMAKE_DL_LIBS})
endif()
#23. link a target to given libraries:
target_link_libraries(${EXECUTABLE_TARGET_NAME} ${LINK_LIBRARY_LIST})

if(USE_TecplotLib)
  if(WIN32)
    if(USE_Windows_X64)
      add_custom_command(TARGET ${EXECUTABLE_TARGET_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PROJECT_SOURCE_DIR}/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x64/tecio.dll" $<TARGET_FILE_DIR:${EXECUTABLE_TARGET_NAME}>)
    else(USE_Windows_X64) 
      add_custom_command(TARGET ${EXECUTABLE_TARGET_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PROJECT_SOURCE_DIR}/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x86/tecio.dll" $<TARGET_FILE_DIR:${EXECUTABLE_TARGET_NAME}>)
    endif(USE_Windows_X64)
  endif(WIN32)
endif(USE_TecplotLib)

# Some pre-definitions that occ need.
if(USE_OCC)
  target_link_libraries(${EXECUTABLE_TARGET_NAME} ${OCC_LIBS})  
  add_definitions(-DWNT)
  add_definitions(-DNOTCL)
endif(USE_OCC)

# For OpenMP
if(USE_OMP)
find_package (OpenMP)
if (OpenMP_FOUND)

if(WIN32)
target_include_directories (${API_TARGET_NAME}
    PRIVATE
        ${OpenMP_CXX_INCLUDE_DIRS}
)

target_compile_options (${API_TARGET_NAME}
    PRIVATE
        ${OpenMP_CXX_FLAGS}
)

target_compile_features (${API_TARGET_NAME}
    PUBLIC 
        cxx_std_14
)

target_link_libraries (${API_TARGET_NAME}
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

# For UGKS

#24. exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")

