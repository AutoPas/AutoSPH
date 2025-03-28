cmake_minimum_required(VERSION 3.21)
project(AutoSPH)

set(CMAKE_CXX_STANDARD 17)

# Include FetchContent logic from another file
include(cmake/modules/autopas.cmake)

# Create an executable
add_executable(AutoSPH src/sph.cpp)

# Link the fetched library
target_link_libraries(AutoSPH PRIVATE autopasfetch)

# Include the fetched library headers
target_include_directories(AutoSPH PRIVATE autopasfetch)
