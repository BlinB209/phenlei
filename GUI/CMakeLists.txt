#! This framework is developed based on QT 5.12 and VTK 7.1
#! If you cannot compile the project, configure the corresponding library files for the QT and VTK versions
PROJECT(PHengLEI)
set(AppName PHengLEI)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.0)

# Add qt lib
if(UNIX)
    set(path_to_qt5_root /vol8/home/tool/qt5)
else()
    set(path_to_qt5_root C:/Qt/Qt5.12.12/5.12.12/msvc2015_64)
endif()

SET(CMAKE_PREFIX_PATH
    ${path_to_qt5_root}
	${path_to_qt5_root}/lib/cmake/Qt5Core
	${path_to_qt5_root}/lib/cmake/Qt5Gui
	${path_to_qt5_root}/lib/cmake/Qt5Widgets)

SET(CMAKE_INCLUDE_CURRENT_DIR ON)
SET(CMAKE_AUTOMOC ON)
SET(CMAKE_AUTOUIC ON)
SET(CMAKE_AUTORCC ON)
set(AUTOGEN_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/temp)

# Set the compiler parameters
if(MSVC)
   set(CMAKE_CXX_FLAGS "/utf-8 ${CMAKE_CXX_FLAGS}")
endif()

# Add an output path
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Add the path of the vtk in release mode
set(VTKPath ${PROJECT_SOURCE_DIR}/3rdparty/VTK/x64)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set(VTKPath ${PROJECT_SOURCE_DIR}/3rdparty/VTK/x64)
endif()
message("vtk path ${VTKPath}")
IF(NOT VTK_BINARY_DIR)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        MESSAGE( " ^^^^^^^^^^ x64")
        FIND_PACKAGE(VTK  REQUIRED
                   NO_MODULE
                   PATHS ${VTKPath}
                   NO_DEFAULT_PATH
        )
    else()
        MESSAGE( " ********** x86")
        FIND_PACKAGE(VTK  REQUIRED
                    NO_MODULE
                    PATHS ${VTKPath}
                    NO_DEFAULT_PATH
        )
    endif()
IF(NOT VTK_DIR)
    MESSAGE(FATAL_ERROR "Please set VTK_DIR.")
ENDIF(NOT VTK_DIR)
    INCLUDE(${VTK_USE_FILE})
ENDIF(NOT VTK_BINARY_DIR)

# use what QVTK built with
SET(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")

find_package(Qt5 COMPONENTS
                    Widgets
                    Core
                    Gui
                 REQUIRED)
if(WIN32)
    find_package(Qt5 COMPONENTS WinExtras REQUIRED)
endif()

LINK_DIRECTORIES(3rdparty)
add_definitions(-D VS)

# Add all the source files
SET(MainWindow_CPP_PATH ./MainWindow/src)
SET(MainWindow_SOURCES 
	${MainWindow_CPP_PATH}/GuiPHengLEI.cpp	
	${MainWindow_CPP_PATH}/main.cpp
	${MainWindow_CPP_PATH}/DefaultWidget.cpp
	${MainWindow_CPP_PATH}/TextView.cpp
)
source_group(src FILES ${MainWindow_SOURCES})
include_directories(${MainWindow_CPP_PATH})

# Add all the header files
SET(MainWindow_H_PATH ./MainWindow/include)
SET(MainWindow_HEADERS
	${MainWindow_H_PATH}/GuiPHengLEI.h
	${MainWindow_H_PATH}/DefaultWidget.h
	${MainWindow_H_PATH}/TextView.h
)
source_group(include FILES ${MainWindow_HEADERS})
include_directories(${MainWindow_H_PATH})
# Add the .ui files
SET(PHengLEI_FORMS 
	MainWindow/src/GuiPHengLEI.ui
)
# Add a resource file 
SET(PHengLEI_RESOURCES 
  ./Resources/images/images.qrc
  ./Resources/images/png.qrc
)
SET(EXE_ICON
  ./Resources/logo.rc
)

# The rest should just work 
QT5_WRAP_UI(UISrcs ${PHengLEI_FORMS})
qt5_add_resources(ResourceSrcs ${PHengLEI_RESOURCES})

SOURCE_GROUP("Resources" FILES
  ${PHengLEI_RESOURCES}
  ${EXE_ICON}
)

ADD_DEFINITIONS(-DQT_GUI_LIBS -DQT_CORE_LIB -DQT3_SUPPORT)

ADD_EXECUTABLE( ${AppName}
    ${PHengLEI_FORMS}
    ${MainWindow_SOURCES}
    ${MainWindow_HEADERS}
    ${ResourceSrcs}
    ${BigResourceSrcs}
    ${EXE_ICON}
    ${UI_FILES}
)

set(LinkList
    Qt5::Core
    Qt5::Gui
    Qt5::Widgets
    ${VTK_LIBRARIES}
    )
TARGET_LINK_LIBRARIES(${AppName}
    ${LinkList}
 )
