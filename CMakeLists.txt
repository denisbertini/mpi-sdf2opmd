cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(SDFROOT C CXX Fortran)

##### PATH OF CMAKE PACKAGES
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Compilers
enable_language(Fortran CXX)

set(CMAKE_Fortran_COMPILER mpif90)

# MPI is mandatory
find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX -fpermissive)

# OpenMP
find_package(OpenMP REQUIRED)


# EPOCH
set(EPOCH_PATH $ENV{EPOCH_ROOT})
set(SDF_PATH $ENV{EPOCH_ROOT}/SDF/FORTRAN)
set(SDF_LIB sdf)

# Compilation options
set(CMAKE_CXX_FLAGS "-fPIC -Wshadow -Wall -W -Woverloaded-virtual -fsigned-char -Wextra -std=c++17 -g")

message(STATUS "MPI_INCLUDE_PATH : ${MPI_INCLUDE_PATH}")
message(STATUS "C++ compiler is : ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ compiler Flags are : ${CMAKE_CXX_FLAGS}")
message(STATUS "C++ compiler identification is : ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "EPOCH     PATH is : ${EPOCH_PATH}")
message(STATUS "EPOCH SDF PATH is : ${SDF_PATH}")


add_subdirectory( sdf2opmd_2d )


