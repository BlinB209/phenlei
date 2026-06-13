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
#! @file      PHMacro.cmake
#! @brief     It defines the macro of USE_XX precondition.
#! @author    Bell.

#Mac function to define pre-definition.
macro(opt OPTION HELP VALUE)
  option(USE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item USE_${OPTION}\n${HELP} (default: ${VALUE})")
  if(${VALUE} STREQUAL "ON")
    message(STATUS "${OPTION} is turned on by default")
  else()
    message(STATUS "${OPTION} is turned off by default")
  endif()
endmacro(opt)

macro(opt_build OPTION HELP VALUE)
  option(BUILD_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item BUILD_${OPTION}\n${HELP} (default: ${VALUE})")
  if(${VALUE} STREQUAL "ON")
    message(STATUS "${OPTION} is turned on by default")
  else()
    message(STATUS "${OPTION} is turned off by default")
  endif()
endmacro(opt_build)

macro(includeAPIHeadFiles path)
message(STATUS "includeAPIHeadFiles is ${path}")
include_directories(${path}/DataStruct/include)
include_directories(${path}/Data/include)
include_directories(${path}/Math/include)
include_directories(${path}/FYMPI/include)
include_directories(${path}/Toolkit/include)
include_directories(${path}/IO/include)
include_directories(${path}/Geometry/include)
include_directories(${path}/Physics/include)
include_directories(${path}/PreProcess/include)
include_directories(${path}/PostProcess/include)
include_directories(${path}/Common/include)
endmacro(includeAPIHeadFiles)

macro(includeHeadFiles)
message(STATUS "includeHeadFiles is ${PROJECT_SOURCE_DIR}")
set(DataStruct_H_PATH   ${PROJECT_SOURCE_DIR}/../API/DataStruct/include)
include_directories(${DataStruct_H_PATH})
set(Data_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Data/include)
include_directories(${Data_H_PATH})
set(Math_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Math/include)
include_directories(${Math_H_PATH})
set(PHMPI_H_PATH   ${PROJECT_SOURCE_DIR}/../API/FYMPI/include)
include_directories(${PHMPI_H_PATH})
set(Toolkit_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Toolkit/include)
include_directories(${Toolkit_H_PATH})
set(IO_H_PATH   ${PROJECT_SOURCE_DIR}/../API/IO/include)
include_directories(${IO_H_PATH})
set(Geometry_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Geometry/include)
include_directories(${Geometry_H_PATH})
set(Physics_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Physics/include)
include_directories(${Physics_H_PATH})
set(PreProcess_H_PATH   ${PROJECT_SOURCE_DIR}/../API/PreProcess/include)
include_directories(${PreProcess_H_PATH})
set(PostProcess_H_PATH   ${PROJECT_SOURCE_DIR}/../API/PostProcess/include)
include_directories(${PostProcess_H_PATH})
set(Common_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Common/include)
include_directories(${Common_H_PATH})
set(HyFLOW_H_PATH   ${PROJECT_SOURCE_DIR}/../HyFLOW/include)
include_directories(${HyFLOW_H_PATH})
set(MultiSpecies_H_PATH   ${PROJECT_SOURCE_DIR}/../API/MultiSpecies/include)
include_directories(${MultiSpecies_H_PATH})
endmacro(includeHeadFiles)

macro(include3rpartyHeadFiles)
message(STATUS "include3rpartyHeadFiles is ${PROJECT_SOURCE_DIR}")
set(HDF5_H_PATH    ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/include)
include_directories("${HDF5_H_PATH}")

if(USE_TecplotLib)
set(TECIO_H_PATH   ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/TECIO/include)
include_directories("${TECIO_H_PATH}")
endif(USE_TecplotLib)

set(PARMETIS_H_PATH  ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/parmetis/include)
include_directories("${PARMETIS_H_PATH}")

set(MEITIS_H_PATH    ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/metis/include)
include_directories("${MEITIS_H_PATH}")

set(CGNS_H_PATH      ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/cgns/include)
include_directories("${CGNS_H_PATH}")

set(MGRID_H_PATH    ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/mgrid/include)
include_directories("${MGRID_H_PATH}")
endmacro(include3rpartyHeadFiles)

macro(FIND_INCLUDE_DIR result curdir)                                    #定义函数,2个参数:存放结果result；指定路径curdir；
    file(GLOB_RECURSE children "${curdir}/*.h")                          #遍历获取{curdir}中*.hpp和*.h文件列表
    set(dirlist "")                                                      #定义dirlist中间变量，并初始化
    foreach(child ${children})                                           #for循环
        string(REGEX REPLACE "(.*)/.*" "\\1" LIB_NAME ${child})          #字符串替换,用/前的字符替换/*h
        if(IS_DIRECTORY ${LIB_NAME})                                     #判断是否为路径
            list (FIND dirlist ${LIB_NAME} list_index)                   #去重，查找dirlist中是否有${LIB_NAME}指定的值，可以区分字符串相同数子后缀不同的路径：例如/app/test_1和/app/test_2
            if(${list_index} LESS 0)                                     #若没找到则代表列表中没有该路径
                LIST(APPEND dirlist ${LIB_NAME})                         #将合法的路径加入dirlist变量中  
            endif()                                                       #结束判断
        endif()                                                            
    endforeach()                                                          #结束for循环
    set(${result} ${dirlist})                                            #dirlist结果放入result变量中
endmacro()

