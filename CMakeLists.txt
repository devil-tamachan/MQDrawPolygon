cmake_minimum_required(VERSION 2.8)
project (MQDrawPolygon)

add_definitions(-D_UNICODE -DUNICODE)

find_package(CGAL REQUIRED COMPONENTS Core)
include( ${CGAL_USE_FILE} )

# find_package(GMP REQUIRED)
# include_directories( ${GMP_INCLUDE_DIR} )

# find_package(MPFR REQUIRED)
# include_directories( ${MPFR_INCLUDE_DIR} )

find_package( OpenCV REQUIRED )

add_definitions( -DBOOST_ALL_NO_LIB )
find_package(Boost 1.64.0 COMPONENTS program_options filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIR})

set(MQSDK_DIR "" CACHE PATH "MQInit.cpp")

if(EXISTS ${MQSDK_DIR}/mqsdk/MQInit.cpp)
  set(MQSDK_INCLUDE_DIR ${MQSDK_DIR}/mqsdk/)
elseif(EXISTS ${MQSDK_DIR}/MQInit.cpp)
  set(MQSDK_INCLUDE_DIR ${MQSDK_DIR})
endif()

set(EIGEN3_DIR "" CACHE PATH "")
if(EXISTS ${EIGEN3_DIR}/Eigen/Dense)
  set(EIGEN3_INCLUDE_DIR ${EIGEN3_DIR})
elseif(EXISTS ${EIGEN3_DIR}/Dense)
  set(EIGEN3_INCLUDE_DIR ${EIGEN3_DIR}/../)
elseif(EXISTS ${EIGEN3_DIR}/include/eigen3/Eigen/Dense)
  set(EIGEN3_INCLUDE_DIR ${EIGEN3_DIR}/include/eigen3/)
endif()

if(NOT DEFINED MQSDK_INCLUDE_DIR)
  message(FATAL_ERROR "MQSDK_DIR not found!")
endif()

INCLUDE_DIRECTORIES( "${MQSDK_INCLUDE_DIR}" )

if(NOT DEFINED EIGEN3_INCLUDE_DIR)
  message(FATAL_ERROR "EIGEN3_DIR not found!")
endif()

INCLUDE_DIRECTORIES( "${EIGEN3_INCLUDE_DIR}" )

add_library(MQDrawPolygon SHARED
  MQDrawPolygon.cpp
  ${MQSDK_INCLUDE_DIR}/MQInit.cpp
  ${MQSDK_INCLUDE_DIR}/MQPlugin.cpp
  ${MQSDK_INCLUDE_DIR}/MQWidget.cpp
)
target_link_libraries( MQDrawPolygon ${GMP_LIBRARIES} )
target_link_libraries( MQDrawPolygon ${MPFR_LIBRARIES} )
target_link_libraries( MQDrawPolygon ${OpenCV_LIBS} )
target_link_libraries( MQDrawPolygon ${Boost_LIBRARIES} )