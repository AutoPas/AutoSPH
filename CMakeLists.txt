cmake_minimum_required(VERSION 3.14)
project(AutoSPH)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
option(GIT_SUBMODULES_SSH "Use SSH for git submodules instead of HTTPS" OFF)

include(autopas)      # loads cmake/modules/autopas.cmake

#file(GLOB_RECURSE MY_SRC "${PROJECT_SOURCE_DIR}/*.cpp" "${PROJECT_SOURCE_DIR}/*.h")

find_package(MPI REQUIRED)

add_executable(sph-exec main.cpp)
target_link_libraries(sph-exec
  PUBLIC
    ${AUTOPAS_LIB}
    MPI::MPI_CXX
)
target_compile_features(sph-exec PRIVATE cxx_std_17)
