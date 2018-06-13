cmake_minimum_required(VERSION 3.0)
project(coefficient-flow CXX)

# requirements:
#   * Gudhi under /usr/include/gudhi
#   * Eigen3 under /usr/include/eigen3
#   * boost under /usr/include/boost

# in case Eigen, and/or GUDHI are not found by find_package
# try setting the paths to them manually
#
# set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3")
# include_directories(${EIGEN3_INCLUDE_DIR})
#

# GUDHI version should be 1.3 or higher, provided they don't break the
# SimplexTree interface
#

add_subdirectory(yaml-cpp)
include_directories("./yaml-cpp/include/")

add_subdirectory(KDTree)
include_directories("./KDTree")

add_subdirectory(tinyply)
include_directories("./tinyply/source/")

add_subdirectory(tinyobjloader)
include_directories("./tinyobjloader")

find_package(Boost REQUIRED)

find_package(Eigen3 3.3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(GUDHI 1.3 REQUIRED)
include_directories(${GUDHI_INCLUDE_DIRS})

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

include_directories( "./lib" )

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

# this compiles it against python 2.7
#find_package(PythonInterp REQUIRED)
#find_package(PythonLibs REQUIRED)
#find_package(Boost REQUIRED COMPONENTS python3)

#set(PYTHON_INCLUDE_DIRS /usr/include/python3.6m )
#set(PYTHON_LIBRARIES /usr/lib/libpython3.6m.so /usr/lib/libpython3.6m.so.1.0 /usr/lib/libpython3.so)

#include_directories(${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
#link_directories(${Boost_LIBRARY_DIR})

#message(STATUS "python version: ${PYTHON_VERSION} ")
#message(STATUS "python dirs: ${PYTHON_INCLUDE_DIRS}")
#message(STATUS "python libs: ${PYTHON_LIBRARIES}")
#message(STATUS "boost dirs: ${Boost_INCLUDE_DIRS}")

find_package(pybind11 REQUIRED)
add_library(coeffflow MODULE python/bindings.cpp)
target_link_libraries(coeffflow PRIVATE scomplex pybind11::module)

set_target_properties(coeffflow PROPERTIES PREFIX "" SUFFIX .so)

#set_target_properties(coeffflow PROPERTIES SUFFIX .so)
#set_target_properties(coeffflow PROPERTIES PREFIX "")

#target_link_libraries(coeffflow "python2.7" "boost_python" ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} scomplex)


# This should really work, but it doesnt
#set(Python_ADDITIONAL_VERSIONS 3.6)
#find_package(PythonLibs REQUIRED)
#
#set(PYTHON_INCLUDE_DIR /usr/include/python3.6m)
#set(PYTHON_LIBRARY libpython3.6m.so libpython3.6m.so.1.0 libpython3.so)
#
#include_directories(${PYTHON_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
#python_add_module(coeffflow python/bindings.cpp)
#target_link_libraries(coeffflow "boost_python3" "python3.6m" ${Boost_LIBRARIES} ${PYTHON_LIBRARY} scomplex)
#
