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
#! @file     CMakeLists.txt
#! @brief    It defines the GPU module, using Cuda to accelerate computing.
#! @author   Yong Zhang, Xi Zhang.
cmake_minimum_required(VERSION 3.5)

#set working directory:
set(EXEC_WORKING_DIR PHengLEIv3d0)

#set cmake environment:
set(CMAKE_USE_RELATIVE_PATHS ON)

#find nvcc compiler
find_package(CUDA REQUIRED)

#1. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory <GPU>, Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

set(GPU_TARGET_NAME "PHGPU")
set(LIBRARY_OUTPUT_PATH)
set(LIB_API_EXPORT_PATH)
if(WIN32)
if(USE_Windows_X64)
 set(LIB_API_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../GPU/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(CMAKE_EXE_LINKER_FLAGS "/machine:x64")
else(USE_Windows_X64)
  set(LIB_API_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../GPU/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(USE_Windows_X64)
else(WIN32)
  set(LIB_API_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../GPU/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)
message(STATUS "API library path is: ${LIB_API_EXPORT_PATH}.")
link_directories("${LIB_API_EXPORT_PATH}")
message(STATUS "GPU library path is: ${LIBRARY_OUTPUT_PATH}.")

#set options
if (NOT BUILD_LIBRARY_TYPE)
  set (BUILD_LIBRARY_TYPE "Static" CACHE STRING "The type of PHengLEI libraries" FORCE)
  SET_PROPERTY(CACHE BUILD_LIBRARY_TYPE PROPERTY STRINGS Shared Static)
endif()

# define compile macro for nvcc
#add_definitions(-DCUDAUNITTEST)
#add_definitions(-DUNITTESTOUTPUT)
# add warings ignore
add_definitions(-Wattributes)

#3. set environment for MPICH:
set(THIRDPARTY_LIB_DIR)
if(WIN32)
  set(MSMPI_HOME_DIRECTORY "C:/Program Files (x86)/Microsoft SDKs/MPI" CACHE PATH "path to installed MPICH")
  include_directories("${MSMPI_HOME_DIRECTORY}/Include")
  if(USE_Windows_X64)
    link_directories("${MSMPI_HOME_DIRECTORY}/Lib/x64")
    set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/3rdparty/win-x64")
  else(USE_Windows_X64) 
    link_directories("${MSMPI_HOME_DIRECTORY}/Lib/x86")
    set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/3rdparty/win-x86")
  endif(USE_Windows_X64)
else(WIN32)

set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/3rdparty/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif(WIN32)
# specify directory in which the linker will look for libraries for myself libraries:
message(STATUS "The 3rdparty path is: ${THIRDPARTY_LIB_DIR}.") 
link_directories("${THIRDPARTY_LIB_DIR}")

set(CFD_H_PATH   
${PROJECT_SOURCE_DIR}/../CFD/Solver/include
${PROJECT_SOURCE_DIR}/../CFD/TurbSolver/include
${PROJECT_SOURCE_DIR}/../CFD/NSSolver/include
)
include_directories("${CFD_H_PATH}")


#5. set current directory:
# set directory of Manage
set(Manage_H_PATH    Manage/include)
set(Manage_H_LIST
${Manage_H_PATH}/cudaErrorHandle.h
${Manage_H_PATH}/GPUDeviceControl.h
${Manage_H_PATH}/TimeMeasure.h
${Manage_H_PATH}/OutputDebug.h
)
source_group(Manage\\include FILES ${Manage_H_LIST}) 
include_directories(${Manage_H_PATH})

set(Manage_CPP_PATH       Manage/src)
set(Manage_CPP_LIST
${Manage_CPP_PATH}/cudaErrorHandle.cu
${Manage_CPP_PATH}/GPUDeviceControl.cu
${Manage_CPP_PATH}/TimeMeasure.cu
${Manage_CPP_PATH}/OutputDebug.cu
)
source_group(Manage\\src FILES ${Manage_CPP_LIST})


# set directory of Kernels:
set(Kernels_H_PATH        Kernels/include)
set(Kernels_H_LIST
${Kernels_H_PATH}/GPUCompNodeVar.h
${Kernels_H_PATH}/GPUFixBCNodeVar.h
${Kernels_H_PATH}/GPUCompGradientGGNode.h
${Kernels_H_PATH}/GPUStoreBoundGrad.h
${Kernels_H_PATH}/GPUCompViscousCoef.h
${Kernels_H_PATH}/GPUFaceWeight.h
${Kernels_H_PATH}/GPUGetVisFaceValue.h
${Kernels_H_PATH}/GPUCompVisfluxTEST.h
${Kernels_H_PATH}/GPULoadFlux.h
${Kernels_H_PATH}/GPUSpectrumRadius1eq.h
${Kernels_H_PATH}/GPUSourceFlux_1eq_Original.h
${Kernels_H_PATH}/GPUTurbAhead.h
${Kernels_H_PATH}/GPUFaceColor.h
${Kernels_H_PATH}/GPUTestFunctions.h
${Kernels_H_PATH}/GPUKernelTest.h
${Kernels_H_PATH}/GPUKernelTestAPI.h
${Kernels_H_PATH}/GPUKernelTestPart2.h
${Kernels_H_PATH}/BasicDeviceVariables.h
${Kernels_H_PATH}/TemporaryOperations.h
${Kernels_H_PATH}/GPUInviscidFlux.h
${Kernels_H_PATH}/TemporaryOperationsPart2.h
${Kernels_H_PATH}/GPUBasicFunctions.h
${Kernels_H_PATH}/GPUBasicFunctionsPart2.h
${Kernels_H_PATH}/GPUNSSolver.h
${Kernels_H_PATH}/GPUInvVisSpectrum.h
${Kernels_H_PATH}/GPUCompGradientGGCell.h
${Kernels_H_PATH}/GPUComputeNodeValue.h
${Kernels_H_PATH}/GPUCompGradientGGNodeNew.h
${Kernels_H_PATH}/GPUMPICommunication.h
${Kernels_H_PATH}/GPUPortability.h
${Kernels_H_PATH}/GPUSpectrumRadius.h
${Kernels_H_PATH}/GPULUSGS.h
${Kernels_H_PATH}/GPUNSUnsteady.h
)
source_group(Kernels\\include FILES ${Kernels_H_LIST}) 
include_directories(${Kernels_H_PATH})

set(Kernels_CPP_PATH      Kernels/src)
set(Kernels_CPP_LIST
${Kernels_CPP_PATH}/GPUCompNodeVar.cu
${Kernels_CPP_PATH}/GPUFixBCNodeVar.cu
${Kernels_CPP_PATH}/GPUCompGradientGGNode.cu
${Kernels_CPP_PATH}/GPUStoreBoundGrad.cu
${Kernels_CPP_PATH}/GPUCompViscousCoef.cu
${Kernels_CPP_PATH}/GPUFaceWeight.cu
${Kernels_CPP_PATH}/GPUGetVisFaceValue.cu
${Kernels_CPP_PATH}/GPUCompVisfluxTEST.cu
${Kernels_CPP_PATH}/GPULoadFlux.cu
${Kernels_CPP_PATH}/GPUSpectrumRadius1eq.cu
${Kernels_CPP_PATH}/GPUSourceFlux_1eq_Original.cu
${Kernels_CPP_PATH}/GPUTurbAhead.cu
${Kernels_CPP_PATH}/GPUFaceColor.cu
${Kernels_CPP_PATH}/GPUTestFunctions.cu
${Kernels_CPP_PATH}/GPUKernelTest.cu
${Kernels_CPP_PATH}/GPUKernelTestAPI.cu
${Kernels_CPP_PATH}/GPUKernelTestPart2.cu
${Kernels_CPP_PATH}/BasicDeviceVariables.cu
${Kernels_CPP_PATH}/TemporaryOperations.cu
${Kernels_CPP_PATH}/GPUInviscidFlux.cu
${Kernels_CPP_PATH}/TemporaryOperationsPart2.cu
${Kernels_CPP_PATH}/GPUBasicFunctions.cu
${Kernels_CPP_PATH}/GPUBasicFunctionsPart2.cu
${Kernels_CPP_PATH}/GPUNSSolver.cu
${Kernels_CPP_PATH}/GPUInvVisSpectrum.cu
${Kernels_CPP_PATH}/GPUCompGradientGGCell.cu
${Kernels_CPP_PATH}/GPUComputeNodeValue.cu
${Kernels_CPP_PATH}/GPUCompGradientGGNodeNew.cu
${Kernels_CPP_PATH}/GPUMPICommunication.cu
${Kernels_CPP_PATH}/GPUPortability.cu
${Kernels_CPP_PATH}/GPUSpectrumRadius.cu
${Kernels_CPP_PATH}/GPULUSGS.cu
${Kernels_CPP_PATH}/GPUNSUnsteady.cu
)
source_group(Kernels\\src FILES ${Kernels_CPP_LIST})


#20.1 set header file list
set(GPU_HEAD_LIST
${Manage_H_LIST}
${Kernels_H_LIST}
)

#20.2 set source file list:
set(GPU_CPP_LIST 
${Manage_CPP_LIST}
${Kernels_CPP_LIST}
)

#21. adds flags to the compiler command line for sources in the current directory and below:
add_definitions(-DPH_PARALLEL -DPH_CGNS)

# cuda settings
set(V100MODE ON)
if (V100MODE)
    MESSAGE(STATUS "Using Nvida Tesla V100")    
#for Nvidia Tesla V100
    if (UNITTESTMODE)
        MESSAGE(STATUS "UNITTESTMODE")    
    #for UNITTESTMODE
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -arch=sm_75 --fmad=false;) #for V100
    else (UNITTESTMODE)
        MESSAGE(STATUS "GPUMODE")    
    #for GPUMODE
        #set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -arch=sm_70;) #for V100
        #set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -arch=sm_70; -pg;) #for V100
        if (GPROF)
            MESSAGE(STATUS "GPROF")    
        #-pg is open
            set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -arch=sm_70; -pg;) #for V100
        else (GPROF)
            set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -arch=sm_70;) #for V100
        endif()
    endif()
else (V100MODE)
        MESSAGE(STATUS "Using Nvida Tesla K80")    
#for Nvidia Tesla K80 or the other low version
    if (UNITTESTMODE)
        MESSAGE(STATUS "UNITTESTMODE")    
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -arch=sm_35 --fmad=false;) #for K80
    else (UNITTESTMODE)
        MESSAGE(STATUS "GPUMODE")    
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -arch=sm_35;) #for K80
    endif()
endif()
# turn on the seperate compile for seperate device function file
SET(CUDA_SEPARABLE_COMPILATION ON)

# Add subdirectory or link library PHAPI.
# Some pre-definitions that occ need.
if(USE_OCC)
  #target_link_libraries(${EXECUTABLE_TARGET_NAME} ${OCC_LIBS})  
  add_definitions(-DWNT)
  add_definitions(-DNOTCL)
endif()

# build executable program or library to the project using the specified source files and header files:
if ("${BUILD_LIBRARY_TYPE}" STREQUAL "Shared")
  # Building Shared Library.   
  message(STATUS "Building Shared Library ${GPU_TARGET_NAME} ...")
  cuda_add_library(${GPU_TARGET_NAME} SHARED ${GPU_CPP_LIST} ${GPU_HEAD_LIST})
elseif ("${BUILD_LIBRARY_TYPE}" STREQUAL "Static")
  # Building Static Library.
  message(STATUS "Building Static Library ${GPU_TARGET_NAME} ...")
  cuda_add_library(${GPU_TARGET_NAME} STATIC ${GPU_CPP_LIST} ${GPU_HEAD_LIST}) 
endif() 


# link API library:
set(API_TARGET_NAME "PHAPI")
message(STATUS "API library name is: ${API_TARGET_NAME}.")
list(APPEND LINK_LIBRARY_LIST ${API_TARGET_NAME})
target_link_libraries(${GPU_TARGET_NAME} ${LINK_LIBRARY_LIST})

list(APPEND LINK_LIBRARY_LIST metis parmetis cgns hdf5 tecio)
if(WIN32)
  list(APPEND LINK_LIBRARY_LIST msmpi)
endif(WIN32)
message(STATUS "Link library list: ${LINK_LIBRARY_LIST}.")

#24. exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
