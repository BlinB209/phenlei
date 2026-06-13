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

#1. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory <API>, Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

set(API_TARGET_NAME "PHAPI")
set(LIBRARY_OUTPUT_PATH)
if(WIN32)
if(CMAKE_CL_64)
 set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(CMAKE_EXE_LINKER_FLAGS "/machine:x64")
 message(STATUS "USE_Windows_X64 is turned on!")
 add_definitions(-DUSE_WINDOWS_X64)
else(CMAKE_CL_64)
  message(STATUS "Windows_X86 is turned on!")
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(CMAKE_CL_64)
else(WIN32)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)

message(STATUS "API library path is: ${LIBRARY_OUTPUT_PATH}.")

#set options
if(NOT BUILD_LIBRARY_TYPE)
  set (BUILD_LIBRARY_TYPE "Static" CACHE STRING "The type of PHengLEI libraries" FORCE)
  SET_PROPERTY(CACHE BUILD_LIBRARY_TYPE PROPERTY STRINGS Shared Static)
endif()

#3. set environment for MPICH and 3rdParty lib directory:
if(WIN32)
  include_directories("${MSMPI_HOME_DIRECTORY}/Include")
  link_directories("${MSMPI_HOME_DIRECTORY}/lib")
  if(CMAKE_CL_64)
    link_directories("${MSMPI_HOME_DIRECTORY}/Lib/x64")
    set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x64")
  else(CMAKE_CL_64) 
    link_directories("${MSMPI_HOME_DIRECTORY}/Lib/x86")
    set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x86")
  endif(CMAKE_CL_64)
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

set(MESH_H_PATH   ${PROJECT_SOURCE_DIR}/../Mesh/include)
include_directories(${MESH_H_PATH})

# specify directory in which the linker will look for libraries for myself libraries:
message(STATUS "The 3rdparty path is: ${THIRDPARTY_LIB_DIR}.") 
link_directories("${THIRDPARTY_LIB_DIR}")

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
set(TECIO_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/TECIO/include)
include_directories("${TECIO_H_PATH}")
endif(USE_TecplotLib)

set(PARMETIS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/parmetis/include)
include_directories("${PARMETIS_H_PATH}")

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

set(MGRID_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/mgrid/include)
include_directories("${MGRID_H_PATH}")

add_definitions(-DPH_PARALLEL -DPH_CGNS)
#set(API_H_DIR ${PROJECT_SOURCE_DIR}/3rdparty)

# set directory of DataStruct
set(DataStruct_H_PATH DataStruct/include)
set(DataStruct_H_LIST
${DataStruct_H_PATH}/DataStruct_Array.h
${DataStruct_H_PATH}/DataStruct_AdtNode.h
${DataStruct_H_PATH}/DataStruct_AdtTree.h
${DataStruct_H_PATH}/DataStruct_KDNode.h
${DataStruct_H_PATH}/DataStruct_KDTree.h
${DataStruct_H_PATH}/DataStruct_BinaryTree.h
${DataStruct_H_PATH}/DataStruct_Sort.h
${DataStruct_H_PATH}/DataStruct_Range.h
${DataStruct_H_PATH}/PHET.h
${DataStruct_H_PATH}/PHETBase.h
${DataStruct_H_PATH}/PHExpr.h
${DataStruct_H_PATH}/PHExprWrap.h
${DataStruct_H_PATH}/PHFastArrayIterator.h
${DataStruct_H_PATH}/PHMemblock.h
${DataStruct_H_PATH}/PHOps.h
${DataStruct_H_PATH}/PHProduct.h
${DataStruct_H_PATH}/PHSimpleArray.h
${DataStruct_H_PATH}/PHStorage.h
${DataStruct_H_PATH}/PHSlice.h
${DataStruct_H_PATH}/Pointer.h
${DataStruct_H_PATH}/PHVector3.h
${DataStruct_H_PATH}/PHMatrix.h
${DataStruct_H_PATH}/TypeDefine.h
${DataStruct_H_PATH}/HyList.h
${DataStruct_H_PATH}/DataContainer.h
${DataStruct_H_PATH}/PHHeapSort.h
)
source_group(DataStruct\\include FILES ${DataStruct_H_LIST}) 
include_directories(${DataStruct_H_PATH})

set(DataStruct_CPP_PATH DataStruct/src)
set(DataStruct_CPP_LIST
${DataStruct_CPP_PATH}/PHStorage.cpp
${DataStruct_CPP_PATH}/DataStruct_BinaryTree.cpp
${DataStruct_CPP_PATH}/DataStruct_KDTree.cpp
${DataStruct_CPP_PATH}/DataContainer.cpp
)
source_group(DataStruct\\src FILES ${DataStruct_CPP_LIST})

# set directory of Data:
set(Data_H_PATH Data/include)
set(Data_H_LIST
${Data_H_PATH}/Data_SafeData.h
${Data_H_PATH}/Data_Util.h
${Data_H_PATH}/Data_IndexData.h
${Data_H_PATH}/Data_Param.h
${Data_H_PATH}/Data_Field.h
${Data_H_PATH}/Data_ParamFieldSuite.h
${Data_H_PATH}/GlobalDataBase.h
)
source_group(Data\\include FILES ${Data_H_LIST}) 
include_directories(${Data_H_PATH})

set(Data_CPP_PATH Data/src)
set(Data_CPP_LIST
${Data_CPP_PATH}/Data_SafeData.cpp
${Data_CPP_PATH}/Data_Util.cpp
${Data_CPP_PATH}/Data_IndexData.cpp
${Data_CPP_PATH}/Data_Param.cpp
${Data_CPP_PATH}/Data_Field.cpp
${Data_CPP_PATH}/Data_ParamFieldSuite.cpp
${Data_CPP_PATH}/GlobalDataBase.cpp
)
source_group(Data\\src FILES ${Data_CPP_LIST})

# set directory of Math:
set(Math_H_PATH Math/include)
set(Math_H_LIST
${Math_H_PATH}/Math_BasisFunction.h
${Math_H_PATH}/Math_Limiter.h
)
source_group(Math\\include FILES ${Math_H_LIST}) 
include_directories(${Math_H_PATH})

set(Math_CPP_PATH Math/src)
set(Math_CPP_LIST
${Math_CPP_PATH}/Math_BasisFunction.cpp
)
source_group(Math\\src FILES ${Math_CPP_LIST})


# set directory of MPI:
set(MPI_H_PATH FYMPI/include)
set(MPI_H_LIST
${MPI_H_PATH}/PHMpi.h
${MPI_H_PATH}/Parallel.h
${MPI_H_PATH}/ActionKey.h
)
source_group(MPI\\include FILES ${MPI_H_LIST}) 
include_directories(${MPI_H_PATH})

set(MPI_CPP_PATH FYMPI/src)
set(MPI_CPP_LIST
${MPI_CPP_PATH}/PHMpi.cpp
${MPI_CPP_PATH}/Parallel.cpp
${MPI_CPP_PATH}/ActionKey.cpp
)
source_group(MPI\\src FILES ${MPI_CPP_LIST})

# set directory of Toolkit:
set(Toolkit_H_PATH Toolkit/include)
set(Toolkit_H_LIST
${Toolkit_H_PATH}/TK_Log.h
${Toolkit_H_PATH}/TK_Exit.h
${Toolkit_H_PATH}/TK_Parse.h
${Toolkit_H_PATH}/TK_Time.h
${Toolkit_H_PATH}/TK_Warning.h
)
source_group(Toolkit\\include FILES ${Toolkit_H_LIST}) 
include_directories(${Toolkit_H_PATH})

set(Toolkit_CPP_PATH Toolkit/src)
set(Toolkit_CPP_LIST
${Toolkit_CPP_PATH}/TK_Log.cpp
${Toolkit_CPP_PATH}/TK_Exit.cpp
${Toolkit_CPP_PATH}/TK_Parse.cpp
${Toolkit_CPP_PATH}/TK_Time.cpp
${Toolkit_CPP_PATH}/TK_Warning.cpp
)
source_group(Toolkit\\src FILES ${Toolkit_CPP_LIST})

# set directory of IO:
set(IO_H_PATH IO/include)
set(IO_H_LIST
${IO_H_PATH}/PHIO.h
${IO_H_PATH}/IO_VirtualFile.h
${IO_H_PATH}/IO_FileReader.h
${IO_H_PATH}/IO_HDF5File.h
${IO_H_PATH}/IO_FileName.h
)
source_group(IO\\include FILES ${IO_H_LIST}) 
include_directories(${IO_H_PATH})

set(IO_CPP_PATH IO/src)
set(IO_CPP_LIST
${IO_CPP_PATH}/PHIO.cpp
${IO_CPP_PATH}/IO_VirtualFile.cpp
${IO_CPP_PATH}/IO_FileReader.cpp
${IO_CPP_PATH}/IO_HDF5File.cpp
${IO_CPP_PATH}/IO_FileName.cpp
)
source_group(IO\\src FILES ${IO_CPP_LIST})

# set directory of Geometry.
set(Geometry_H_PATH Geometry/include)
set(Geometry_H_LIST
${Geometry_H_PATH}/Geo_Grid.h
${Geometry_H_PATH}/Geo_GridManager.h
${Geometry_H_PATH}/GridType.h
${Geometry_H_PATH}/Geo_Interface.h
${Geometry_H_PATH}/Geo_Interpoint.h
${Geometry_H_PATH}/Geo_SimpleGrid.h
${Geometry_H_PATH}/Geo_StructGrid.h
${Geometry_H_PATH}/Geo_UnstructGrid.h
${Geometry_H_PATH}/MultigridFactory.h
${Geometry_H_PATH}/Geo_CellMetrics_Struct.h
${Geometry_H_PATH}/Geo_CellMetrics_Unstruct.h
${Geometry_H_PATH}/Geo_CellTopo_Unstruct.h
${Geometry_H_PATH}/Geo_DynamicGridMetrics_Struct.h
${Geometry_H_PATH}/Geo_DynamicGridMetrics_Unstruct.h
${Geometry_H_PATH}/Geo_FaceMetrics_Struct.h
${Geometry_H_PATH}/Geo_FaceMetrics_Unstruct.h
${Geometry_H_PATH}/Geo_FaceTopo_Unstruct.h
${Geometry_H_PATH}/Geo_LineTopo_Unstruct.h
${Geometry_H_PATH}/Geo_GridIndex.h
${Geometry_H_PATH}/Geo_LSQWeight_Unstruct.h
${Geometry_H_PATH}/Geo_MultiGridInfo_Struct.h
${Geometry_H_PATH}/Geo_NodeTopo_Struct.h
${Geometry_H_PATH}/Geo_NodeTopo_Unstruct.h
${Geometry_H_PATH}/Geo_OversetGridTopo_Struct.h
${Geometry_H_PATH}/Geo_SimpleBC.h
${Geometry_H_PATH}/Geo_StructBC.h
${Geometry_H_PATH}/Geo_UnstructBC.h
${Geometry_H_PATH}/Geo_OversetGrid.h
${Geometry_H_PATH}/OversetGridFactory.h
${Geometry_H_PATH}/OversetInformation.h
${Geometry_H_PATH}/MeshAgglomeration.h
${Geometry_H_PATH}/Geometry.h
${Geometry_H_PATH}/LinkStruct.h
${Geometry_H_PATH}/Geo_Point.h
${Geometry_H_PATH}/Geo_Face.h
${Geometry_H_PATH}/GRHash.h
${Geometry_H_PATH}/OverLappingGrid.h
${Geometry_H_PATH}/GeometryUnit.h
${Geometry_H_PATH}/Geo_PointFactory.h
${Geometry_H_PATH}/Geo_Element.h
${Geometry_H_PATH}/Geo_Sample.h
${Geometry_H_PATH}/GeomProxy.h
${Geometry_H_PATH}/FieldProxy.h
${Geometry_H_PATH}/GeneralFieldProxy.h
${Geometry_H_PATH}/FaceProxy.h
${Geometry_H_PATH}/Geo_SimpleVC.h
${Geometry_H_PATH}/Geo_StructMerge.h
${Geometry_H_PATH}/Geo_UnstructMerge.h
)
source_group(Geometry\\include FILES ${Geometry_H_LIST})
include_directories(${Geometry_H_PATH})

set(Geometry_CPP_PATH Geometry/src)
set(Geometry_CPP_LIST
${Geometry_CPP_PATH}/Geo_Grid.cpp
${Geometry_CPP_PATH}/Geo_GridManager.cpp
${Geometry_CPP_PATH}/GridType.cpp
${Geometry_CPP_PATH}/Geo_SimpleGrid.cpp
${Geometry_CPP_PATH}/Geo_StructGrid.cpp
${Geometry_CPP_PATH}/Geo_UnstructGrid.cpp
${Geometry_CPP_PATH}/StructGridFD.cpp
${Geometry_CPP_PATH}/Geo_CellMetrics_Struct.cpp
${Geometry_CPP_PATH}/Geo_CellMetrics_Unstruct.cpp
${Geometry_CPP_PATH}/Geo_CellTopo_Unstruct.cpp
${Geometry_CPP_PATH}/Geo_DynamicGridMetrics_Struct.cpp
${Geometry_CPP_PATH}/Geo_DynamicGridMetrics_Unstruct.cpp
${Geometry_CPP_PATH}/Geo_FaceMetrics_Struct.cpp
${Geometry_CPP_PATH}/Geo_FaceMetrics_Unstruct.cpp
${Geometry_CPP_PATH}/Geo_FaceTopo_Unstruct.cpp
${Geometry_CPP_PATH}/Geo_LineTopo_Unstruct.cpp
${Geometry_CPP_PATH}/Geo_GridIndex.cpp
${Geometry_CPP_PATH}/Geo_LSQWeight_Unstruct.cpp
${Geometry_CPP_PATH}/Geo_MultiGridInfo_Struct.cpp
${Geometry_CPP_PATH}/Geo_NodeTopo_Struct.cpp
${Geometry_CPP_PATH}/Geo_NodeTopo_Unstruct.cpp
${Geometry_CPP_PATH}/Geo_OversetGridTopo_Struct.cpp
${Geometry_CPP_PATH}/Geo_SimpleBC.cpp
${Geometry_CPP_PATH}/Geo_StructBC.cpp
${Geometry_CPP_PATH}/Geo_UnstructBC.cpp
${Geometry_CPP_PATH}/Geo_OversetGrid.cpp
${Geometry_CPP_PATH}/OversetGridFactory.cpp
${Geometry_CPP_PATH}/OversetInformation.cpp
${Geometry_CPP_PATH}/MeshAgglomeration.cpp
${Geometry_CPP_PATH}/Geometry.cpp
${Geometry_CPP_PATH}/Geo_Interface.cpp
${Geometry_CPP_PATH}/Geo_Interpoint.cpp
${Geometry_CPP_PATH}/MultigridFactory.cpp
${Geometry_CPP_PATH}/GRHash.cpp
${Geometry_CPP_PATH}/OverLappingGrid.cpp
${Geometry_CPP_PATH}/GeometryUnit.cpp
${Geometry_CPP_PATH}/Geo_PointFactory.cpp
${Geometry_CPP_PATH}/Geo_Element.cpp
${Geometry_CPP_PATH}/Geo_Sample.cpp
${Geometry_CPP_PATH}/GeomProxy.cpp
${Geometry_CPP_PATH}/FieldProxy.cpp
${Geometry_CPP_PATH}/GeneralFieldProxy.cpp
${Geometry_CPP_PATH}/FaceProxy.cpp
${Geometry_CPP_PATH}/Geo_SimpleVC.cpp
${Geometry_CPP_PATH}/Geo_StructMerge.cpp
${Geometry_CPP_PATH}/Geo_UnstructMerge.cpp
)
source_group(Geometry\\src FILES ${Geometry_CPP_LIST})

# set directory of Physics.
set(Physics_H_PATH Physics/include)
set(Physics_H_LIST
${Physics_H_PATH}/AleModel.h
${Physics_H_PATH}/AleManager.h
${Physics_H_PATH}/AleForceManager.h
${Physics_H_PATH}/MovingMesh.h
${Physics_H_PATH}/DeformingSolverManager.h
${Physics_H_PATH}/SixDofManager.h
${Physics_H_PATH}/PHSuperMatrix.h
${Physics_H_PATH}/SixDOF.h
${Physics_H_PATH}/Force.h
${Physics_H_PATH}/Zone.h
${Physics_H_PATH}/Residual.h
${Physics_H_PATH}/GradientOperation.h
${Physics_H_PATH}/AerodynamicsBasicRelations.h
${Physics_H_PATH}/Limiter.h
${Physics_H_PATH}/Task.h
${Physics_H_PATH}/Task_Simple.h
${Physics_H_PATH}/Task_ServerUpdateInterface.h
${Physics_H_PATH}/Task_ServerUpdateInterpoint.h
${Physics_H_PATH}/Task_WriteBinaryTecplot.h
${Physics_H_PATH}/Region.h
${Physics_H_PATH}/IO_GridReader.h
${Physics_H_PATH}/GridGroup.h
${Physics_H_PATH}/PerfectGas.h
${Physics_H_PATH}/Gas.h
${Physics_H_PATH}/Solver.h
${Physics_H_PATH}/Flux_Inviscid.h
${Physics_H_PATH}/Flux_RoeEntropyCorrection.h
${Physics_H_PATH}/Jacobian.h
${Physics_H_PATH}/Gradient.h
${Physics_H_PATH}/MatrixLUSGS.h
)

if(BUILD_IncomSolver)
list(APPEND Physics_H_LIST ${Physics_H_PATH}/IncomGas.h)
endif(BUILD_IncomSolver)

source_group(Physics\\include FILES ${Physics_H_LIST}) 
include_directories(${Physics_H_PATH})

set(Physics_CPP_PATH Physics/src)
set(Physics_CPP_LIST
${Physics_CPP_PATH}/Zone.cpp
${Physics_CPP_PATH}/Force.cpp
${Physics_CPP_PATH}/AleForceManager.cpp
${Physics_CPP_PATH}/AleModel.cpp
${Physics_CPP_PATH}/AleManager.cpp
${Physics_CPP_PATH}/MovingMesh.cpp
${Physics_CPP_PATH}/DeformingSolverManager.cpp
${Physics_CPP_PATH}/SixDofManager.cpp
${Physics_CPP_PATH}/SixDOF.cpp
${Physics_CPP_PATH}/Residual.cpp
${Physics_CPP_PATH}/GradientOperation.cpp
${Physics_CPP_PATH}/AerodynamicsBasicRelations.cpp
${Physics_CPP_PATH}/Limiter.cpp
${Physics_CPP_PATH}/Task.cpp
${Physics_CPP_PATH}/Task_Simple.cpp
${Physics_CPP_PATH}/Task_ServerUpdateInterface.cpp
${Physics_CPP_PATH}/Task_ServerUpdateInterpoint.cpp
${Physics_CPP_PATH}/Task_WriteBinaryTecplot.cpp
${Physics_CPP_PATH}/Region.cpp
${Physics_CPP_PATH}/IO_GridReader.cpp
${Physics_CPP_PATH}/GridGroup.cpp
${Physics_CPP_PATH}/PerfectGas.cpp
${Physics_CPP_PATH}/Solver.cpp
${Physics_CPP_PATH}/Flux_Inviscid.cpp
${Physics_CPP_PATH}/Flux_RoeEntropyCorrection.cpp
${Physics_CPP_PATH}/Jacobian.cpp
${Physics_CPP_PATH}/Gradient.cpp
${Physics_CPP_PATH}/Gas.cpp
${Physics_CPP_PATH}/MatrixLUSGS.cpp
)

if(BUILD_IncomSolver)
list(APPEND Physics_CPP_LIST ${Physics_CPP_PATH}/IncomGas.cpp)
endif(BUILD_IncomSolver)

source_group(Physics\\src FILES ${Physics_CPP_LIST})

# set directory of PreProcess:
set(PreProcess_H_PATH PreProcess/include)
set(PreProcess_H_LIST
${PreProcess_H_PATH}/Pre_Block_Struct.h
${PreProcess_H_PATH}/Pre_UnstructGridPartition.h
${PreProcess_H_PATH}/Pre_WalldistCompute.h
${PreProcess_H_PATH}/Pre_GridConversion.h
${PreProcess_H_PATH}/Pre_FieldViewToPHengLEI_Unstruct.h
${PreProcess_H_PATH}/Pre_CGNSToPHengLEI_Unstruct.h
${PreProcess_H_PATH}/Pre_CGNSConversion_Struct.h
${PreProcess_H_PATH}/Pre_CGNSToPlot3D_Struct.h
${PreProcess_H_PATH}/Pre_CGNSToPHengLEI_Struct.h
${PreProcess_H_PATH}/Pre_Plot3DToPHengLEI_Struct.h
${PreProcess_H_PATH}/Pre_FluentToPHengLEI_Unstruct.h
${PreProcess_H_PATH}/Pre_HDF5File.h
${PreProcess_H_PATH}/Pre_GridMirror.h
${PreProcess_H_PATH}/Pre_StructGridPartition.h
${PreProcess_H_PATH}/Pre_StructToUnstruct.h
${PreProcess_H_PATH}/Pre_GridTranslator.h
${PreProcess_H_PATH}/Pre_GridCombine.h
${PreProcess_H_PATH}/Pre_GridBase.h
${PreProcess_H_PATH}/MixGrid.h
${PreProcess_H_PATH}/Pre_Reorder.h
)
source_group(PreProcess\\include FILES ${PreProcess_H_LIST}) 
include_directories(${PreProcess_H_PATH})

set(PreProcess_CPP_PATH PreProcess/src)
set(PreProcess_CPP_LIST
${PreProcess_CPP_PATH}/Pre_Block_Struct.cpp
${PreProcess_CPP_PATH}/Pre_UnstructGridPartition.cpp
${PreProcess_CPP_PATH}/Pre_WalldistCompute.cpp
${PreProcess_CPP_PATH}/Pre_GridConversion.cpp
${PreProcess_CPP_PATH}/Pre_FieldViewToPHengLEI_Unstruct.cpp
${PreProcess_CPP_PATH}/Pre_CGNSToPHengLEI_Unstruct.cpp
${PreProcess_CPP_PATH}/Pre_CGNSConversion_Struct.cpp
${PreProcess_CPP_PATH}/Pre_CGNSToPlot3D_Struct.cpp
${PreProcess_CPP_PATH}/Pre_CGNSToPHengLEI_Struct.cpp
${PreProcess_CPP_PATH}/Pre_Plot3DToPHengLEI_Struct.cpp
${PreProcess_CPP_PATH}/Pre_FluentToPHengLEI_Unstruct.cpp
${PreProcess_CPP_PATH}/Pre_HDF5File.cpp
${PreProcess_CPP_PATH}/Pre_GridMirror.cpp
${PreProcess_CPP_PATH}/Pre_StructGridPartition.cpp
${PreProcess_CPP_PATH}/Pre_StructToUnstruct.cpp
${PreProcess_CPP_PATH}/Pre_GridTranslator.cpp
${PreProcess_CPP_PATH}/Pre_GridCombine.cpp
${PreProcess_CPP_PATH}/Pre_GridBase.cpp
${PreProcess_CPP_PATH}/MixGrid.cpp
${PreProcess_CPP_PATH}/Pre_Reorder.cpp
)
source_group(PreProcess\\src FILES ${PreProcess_CPP_LIST})

# set directory of PostProcess:
set(PostProcess_H_PATH PostProcess/include)
set(PostProcess_H_LIST
${PostProcess_H_PATH}/Post_Visual.h
${PostProcess_H_PATH}/Post_WriteVisualFile.h
${PostProcess_H_PATH}/Post_WriteTecplot.h
${PostProcess_H_PATH}/Post_WriteEnsight.h
${PostProcess_H_PATH}/Post_ForceMoment.h
${PostProcess_H_PATH}/PostProcess.h
${PostProcess_H_PATH}/IO_FlowReader.h
${PostProcess_H_PATH}/ForceProcessStruct.h
${PostProcess_H_PATH}/ForceProcessUnstruct.h
${PostProcess_H_PATH}/PostProcessData.h
${PostProcess_H_PATH}/Post_Probes.h
${PostProcess_H_PATH}/Post_ProbesWrite.h
${PostProcess_H_PATH}/Post_WriteParaview.h
)
source_group(PostProcess\\include FILES ${PostProcess_H_LIST}) 
include_directories(${PostProcess_H_PATH})

set(PostProcess_CPP_PATH PostProcess/src)
set(PostProcess_CPP_LIST
${PostProcess_CPP_PATH}/Post_Visual.cpp
${PostProcess_CPP_PATH}/Post_WriteVisualFile.cpp
${PostProcess_CPP_PATH}/Post_WriteTecplot.cpp
${PostProcess_CPP_PATH}/Post_WriteEnsight.cpp
${PostProcess_CPP_PATH}/Post_ForceMoment.cpp
${PostProcess_CPP_PATH}/PostProcess.cpp
${PostProcess_CPP_PATH}/IO_FlowReader.cpp
${PostProcess_CPP_PATH}/ForceProcessStruct.cpp
${PostProcess_CPP_PATH}/ForceProcessUnstruct.cpp
${PostProcess_CPP_PATH}/PostProcessData.cpp
${PostProcess_CPP_PATH}/Post_Probes.cpp
${PostProcess_CPP_PATH}/Post_ProbesWrite.cpp
${PostProcess_CPP_PATH}/Post_WriteParaview.cpp
)
source_group(PostProcess\\src FILES ${PostProcess_CPP_LIST})

# set directory of Common:
set(Common_H_PATH Common/include)
set(Common_H_LIST
${Common_H_PATH}/AMRDef.h
${Common_H_PATH}/Constants.h
${Common_H_PATH}/PHHeader.h
${Common_H_PATH}/Precision.h
${Common_H_PATH}/Glb_Dimension.h
${Common_H_PATH}/LIB_Macro.h
)
source_group(Common\\include FILES ${Common_H_LIST}) 
include_directories(${Common_H_PATH})

set(Common_CPP_PATH Common/src)
set(Common_CPP_LIST
${Common_CPP_PATH}/Glb_Dimension.cpp
)
source_group(Common\\src FILES ${Common_CPP_LIST})

#20.1 set header file list
set(API_HEAD_LIST 
${DataStruct_H_LIST}
${Data_H_LIST}
${Math_H_LIST}
${MPI_H_LIST}
${Toolkit_H_LIST}
${IO_H_LIST}
${Geometry_H_LIST}
${Physics_H_LIST}
${PreProcess_H_LIST}
${PostProcess_H_LIST}
${Common_H_LIST}
)

#20.2 set source file list:
set(API_CPP_LIST 
${DataStruct_CPP_LIST}
${Data_CPP_LIST}
${Math_CPP_LIST}
${MPI_CPP_LIST}
${Toolkit_CPP_LIST}
${IO_CPP_LIST}
${Geometry_CPP_LIST}
${Physics_CPP_LIST}
${PreProcess_CPP_LIST}
${PostProcess_CPP_LIST}
${Common_CPP_LIST}
)

# set directory of MultiSpeciesGas.
set(MultiSpecies_H_PATH MultiSpecies/include)
set(MultiSpecies_H_LIST
${MultiSpecies_H_PATH}/ChemModel.h
${MultiSpecies_H_PATH}/MixedGas.h
)
source_group(MultiSpecies\\include FILES ${MultiSpecies_H_LIST}) 
include_directories(${MultiSpecies_H_PATH})

set(MultiSpecies_CPP_PATH MultiSpecies/src)
set(MultiSpecies_CPP_LIST
${MultiSpecies_CPP_PATH}/ChemModel.cpp
${MultiSpecies_CPP_PATH}/MixedGas.cpp
)
source_group(MultiSpecies\\src FILES ${MultiSpecies_CPP_LIST})
list(APPEND API_HEAD_LIST ${MultiSpecies_H_LIST})
list(APPEND API_CPP_LIST ${MultiSpecies_CPP_LIST})


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

if(BUILD_GPU) 
  add_definitions(-DUSE_CUDA)
  set(GPU_H_PATH   
  ${PROJECT_SOURCE_DIR}/../GPU/Manage/include
  ${PROJECT_SOURCE_DIR}/../GPU/Kernels/include
  )
  include_directories("${GPU_H_PATH}")
endif(BUILD_GPU)

# Some pre-definitions that occ need.
if(USE_OCC)
  #target_link_libraries(${EXECUTABLE_TARGET_NAME} ${OCC_LIBS})  
  add_definitions(-DWNT)
  add_definitions(-DNOTCL)
endif()

# build executable program or library to the project using the specified source files and header files:
if ("${BUILD_LIBRARY_TYPE}" STREQUAL "Shared")
  # Building Shared Library.   
  message(STATUS "Building Shared Library ${API_TARGET_NAME} ...")
  add_library(${API_TARGET_NAME} SHARED ${API_CPP_LIST} ${API_HEAD_LIST})
elseif ("${BUILD_LIBRARY_TYPE}" STREQUAL "Static")
  # Building Static Library.
  message(STATUS "Building Static Library ${API_TARGET_NAME} ...")
  add_library(${API_TARGET_NAME} STATIC ${API_CPP_LIST} ${API_HEAD_LIST}) 
endif()

# link API library:
if(BUILD_GPU)
  message(STATUS "Linking ${GPU_TARGET_NAME} ...")
  list(APPEND LINK_LIBRARY_LIST ${GPU_TARGET_NAME})
endif(BUILD_GPU)
if(WIN32)
  if(USE_VS2019)
    list(APPEND LINK_LIBRARY_LIST metis parmetis cgns libhdf5 libhdf5_hl mgrid IMlib)
  else(USE_VS2019)
    list(APPEND LINK_LIBRARY_LIST metis parmetis cgns hdf5 mgrid IMlib)
  endif()
else(WIN32)
  list(APPEND LINK_LIBRARY_LIST metis parmetis cgns hdf5 hdf5_hl mgrid)
endif(WIN32)

if(USE_TecplotLib)
  list(APPEND LINK_LIBRARY_LIST tecio)
endif(USE_TecplotLib)

if(WIN32)
  list(APPEND LINK_LIBRARY_LIST ${mpiname})
endif(WIN32)

if(BUILD_GMRESSolver)
  # sacado library
  list(APPEND LINK_LIBRARY_LIST ${Sacado_LIBRARIES})
  # PETSc library
  list(APPEND LINK_LIBRARY_LIST ${PETSC_Libraries})
endif()

target_link_libraries(${API_TARGET_NAME} ${LINK_LIBRARY_LIST})

message(STATUS "Link library list: ${LINK_LIBRARY_LIST}.")

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
endif (USE_OMP)

#24. exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")

