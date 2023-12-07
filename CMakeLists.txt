#specefiy required version - this is the newest one
cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#using intel compiler
set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_C_COMPILER icc)

enable_language(CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

enable_testing()

find_package(HDF5 REQUIRED COMPONENTS CXX)
#find_package(MKL CONFIG REQUIRED)
#find_package(MPI REQUIRED)

set(project_sources
        src/eliashbergChain.cpp
        src/fileHandling.cpp
        src/processing/localOneGap.cpp
        )

set(header_files
        include/fileHandling.h
        parameters/globals.h
        include/processing/localOneGap.h
        )

include_directories("./include")
include_directories("./include/processing")
include_directories("./parameters")
#include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(${HDF5_INCLUDE_DIRS})
message(STATUS "---------------------------------------------------------------------")

#add compiler flags
set(C_flags "-O3 -fvar-tracking -Wall -Wextra -pedantic -qmkl -qopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${C_flags}")

#set project name
project(EliashbergChain)
#add an executable - using variables
add_executable(${PROJECT_NAME}.out ${project_sources} ${header_files})

target_link_libraries(${PROJECT_NAME}.out ${HDF5_LIBRARIES})
#target_link_libraries(${PROJECT_NAME}.out ${MPI_CXX_LIBRARIES})


