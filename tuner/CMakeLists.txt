CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)
PROJECT(tuner)
INCLUDE(ExternalProject)

# CMake module(s) path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

set(CMAKE_C_FLAGS "-O3")

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX .)

# Add Ceres solver
ExternalProject_Add(
ceres
URL http://ceres-solver.googlecode.com/files/ceres-solver-1.8.0.tar.gz
TIMEOUT 300
# Disable install step
INSTALL_COMMAND ""
CMAKE_ARGS
	-DLAPACK=OFF
	-DSUITESPARSE=OFF
	-DCXSPARSE=OFF
	-DGFLAGS=OFF
	-DMINIGLOG=ON
# Wrap download, configure and build steps in a script to log output
LOG_DOWNLOAD ON
LOG_CONFIGURE ON
LOG_BUILD ON)

ExternalProject_Get_Property(ceres source_dir)
SET(ceres_dir ${source_dir})

INCLUDE_DIRECTORIES(../include ${ceres_dir}/include
	${ceres_dir}/internal/ceres/miniglog)

ADD_LIBRARY(tuner SHARED tuner.cpp)

ADD_DEPENDENCIES(tuner ukf ceres)

ExternalProject_Get_Property(ceres binary_dir)
SET(cereslib_dir ${binary_dir})

TARGET_LINK_LIBRARIES(tuner
	${cereslib_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}ceres${CMAKE_STATIC_LIBRARY_SUFFIX}
	${cereslib_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}miniglog${CMAKE_STATIC_LIBRARY_SUFFIX}
	ukf)
