cmake_minimum_required(VERSION 3.14)
project(AutoSPH)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
option(GIT_SUBMODULES_SSH "Use SSH for git submodules instead of HTTPS" OFF)

include(autopas)      # loads cmake/modules/autopas.cmake

#file(GLOB_RECURSE MY_SRC "${PROJECT_SOURCE_DIR}/*.cpp" "${PROJECT_SOURCE_DIR}/*.h")

# Disable everything we don't need
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Disable yaml-cpp tests")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Disable yaml-cpp contrib")
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Disable yaml-cpp tools")
set(YAML_CPP_INSTALL OFF CACHE BOOL "Disable yaml-cpp install rules")

add_subdirectory("${CMAKE_BINARY_DIR}/_deps/autopasfetch-src/libs/yaml-cpp")

find_package(MPI REQUIRED)

add_executable(sph-exec main.cpp)
target_link_libraries(sph-exec
  PUBLIC
    ${AUTOPAS_LIB}
    MPI::MPI_CXX
)

target_link_libraries(sph-exec PRIVATE yaml-cpp)

target_compile_features(sph-exec PRIVATE cxx_std_17)
