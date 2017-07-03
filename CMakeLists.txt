cmake_minimum_required(VERSION 2.8)
project(simplicial CXX)

# set(CMAKE_CXX_COMPILER_ID "Clang")
# message("${CMAKE_CXX_COMPILER_ID}")

# requirements: CHANGE LOCATIONS HERE!!!
#   * Gudhi under /usr/include/gudhi
#   * Eigen3 under /usr/include/eigen3
#   * boost under /usr/include/boost
#   * CGAL under /usr/include/CGAL
#
# in case Eigen, and/or CGAL, and/or GUDHI are not found by find_package
# try setting the paths to them manually

#
# set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3")
# include_directories(${EIGEN3_INCLUDE_DIR})
#
# set(CGAL_INCLUDE_DIR "/usr/local/include/eigen3")
# include_directories(${CGAL_INCLUDE_DIR})
#
# GUDHI is being included manually further down
# version should be 1.3 or higher, provided they don't break anything being used
#
# you will also need to have libyaml-cpp installed

find_package(Boost REQUIRED)
message(STATUS "Boost directory: ${Boost_INCLUDE_DIRS}")
include_directories(${Boost_INCLUDE_DIRS})

find_package(Eigen3 3.1.2 REQUIRED)
message(STATUS "Eigen3 version: ${EIGEN3_VERSION}")
message(STATUS "Eigen3 directory: ${EIGEN3_INCLUDE_DIR}")
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(CGAL REQUIRED)
message(STATUS "CGAL version: ${CGAL_VERSION}")
message(STATUS "CGAL directory: ${CGAL_INCLUDE_DIRS}")
include_directories(${CGAL_INCLUDE_DIRS})

set(GUDHI_INCLUDE_DIR "/usr/local/include/gudhi")
if(EXISTS ${GUDHI_INCLUDE_DIR})
    include_directories(${GUDHI_INCLUDE_DIR})
else(EXISTS ${GUDHI_INCLUDE_DIR})
    message(FATAL_ERROR "Gudhi not found")
endif(EXISTS ${GUDHI_INCLUDE_DIR})

if(NOT EXISTS "/usr/include/yaml-cpp")
    message(FATAL_ERROR "libyaml (yaml-cpp) not found")
endif(NOT EXISTS "/usr/include/yaml-cpp")

# find_library(HAVE_YAML NAME YAML_CPP)

set(CMAKE_C_COMPILER             "/usr/bin/clang-3.8")
set(CMAKE_C_FLAGS                "-Wall -std=c99")
set(CMAKE_C_FLAGS_DEBUG          "-g")
set(CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

set(CMAKE_CXX_COMPILER             "/usr/bin/clang++-3.8")
set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG          "-g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

set(CMAKE_AR      "/usr/bin/llvm-ar")
set(CMAKE_LINKER  "/usr/bin/llvm-ld")
set(CMAKE_NM      "/usr/bin/llvm-nm")
set(CMAKE_OBJDUMP "/usr/bin/llvm-objdump")
set(CMAKE_RANLIB  "/usr/bin/llvm-ranlib")

set(CMAKE_BUILD_TYPE Release)

include_directories("./lib") # where the header files are

add_library(scomplex SHARED "./lib/scomplex/simplicial_complex.cpp")
add_library(pathsnap SHARED "./lib/scomplex/path_snapper.cpp")

add_executable(yamltest "src/testing_facility.cpp")
target_link_libraries(yamltest yaml-cpp scomplex pathsnap)
