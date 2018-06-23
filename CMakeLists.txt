cmake_minimum_required(VERSION 3.0)
project(coefficient-flow CXX)

# requirements:
#   * Gudhi under /usr/include/gudhi
#   * Eigen3 under /usr/include/eigen3
#   * boost under /usr/include/boost
#   * libyaml under /usr/local/include/libyaml
#

# in case Eigen, and/or GUDHI are not found by find_package
# try setting the paths to them manually
#
# set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3")
# include_directories(${EIGEN3_INCLUDE_DIR})
#

# GUDHI version should be 1.3 or higher, provided they don't break the
# SimplexTree interface
#
# libyaml-cpp this should actually be found in "/usr/include" so there is no
# real need for the find_package, but it's always nice to make sure
#

# TODO remove
add_subdirectory(KDTree)


add_subdirectory(tinyply)
include_directories("./tinyply/source/")

add_subdirectory(tinyobjloader)
include_directories("./tinyobjloader")

find_package(Boost REQUIRED)
message(STATUS "Boost version: ${Boost_VERSION}")
message(STATUS "Boost directory: ${Boost_INCLUDE_DIRS}")
message(STATUS "boost libraries: ${Boost_LIBRARIES}")

find_package(Eigen3 3.1.2 REQUIRED QUIET)
message(STATUS "Eigen3 version: ${EIGEN3_VERSION}")
message(STATUS "Eigen3 directory: ${EIGEN3_INCLUDE_DIR}")
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(yaml-cpp REQUIRED QUIET)
message(STATUS "yaml-cpp version: ${YAML_CPP_VERSION}")
message(STATUS "yaml-cpp directory: ${YAML_CPP_INCLUDE_DIR}")
include_directories(${YAML_CPP_INCLUDE_DIR})

find_package(GUDHI 1.3 REQUIRED QUIET)
message(STATUS "Gudhi version: ${GUDHI_VERSION}")

set(CMAKE_CXX_COMPILER             "/usr/bin/clang++")
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
# set(CMAKE_BUILD_TYPE Debug)

include_directories(
    "./lib" # where all the library files are
    "KDTree" # for the KDTree library
    )


add_library(scomplex SHARED "./lib/scomplex/simplicial_complex.cpp")
add_library(pathsnap SHARED "./lib/scomplex/path_snapper.cpp")
target_link_libraries(pathsnap KDTree)

add_executable(qhull2ply "src/make_mesh.cpp")

add_executable(yamltest "src/testing_facility.cpp")
target_link_libraries(yamltest yaml-cpp scomplex pathsnap tinyply tinyobjloader)

message(INFO ${CMAKE_CURRENT_SOURCE_DIR})
message(INFO ${CMAKE_CURRENT_BINARY_DIR})
# TEST target
add_custom_target(timing_test
    COMMAND "cp" "${CMAKE_CURRENT_SOURCE_DIR}/test/run_dumb_tests.sh" "${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND "cp" "${CMAKE_CURRENT_SOURCE_DIR}/test/dumbexample.yaml" "${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND "${CMAKE_CURRENT_BINARY_DIR}/run_dumb_tests.sh" "1"
    DEPENDS yamltest qhull2ply "test/run_dumb_tests.sh" "test/dumbexample.yaml")

