# autopas library
message(STATUS "Using AutoPas.")

# Select https (default) or ssh path.
set(autopasRepoPath https://github.com/AutoPas/AutoPas.git)
if (GIT_SUBMODULES_SSH)
set(autopasRepoPath git@github.com:AutoPas/AutoPas.git)
endif ()

set(AUTOPAS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(AUTOPAS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(AUTOPAS_ENABLE_DYNAMIC_CONTAINERS OFF CACHE BOOL "")
set(AUTOPAS_ENABLE_ENERGY_MEASUREMENTS OFF CACHE BOOL "")
set(spdlog_ForceBundled ON CACHE BOOL "" FORCE)
set(Eigen3_ForceBundled ON CACHE BOOL "" FORCE)

# Master as on March 27th, 2026. Includes PMT energy measurements, dynamic rebuilding and 3-body interactions.
set(AUTOPAS_TAG 0338ff17d55cc770a403efee69740a69753d64ad CACHE STRING "AutoPas Git tag or commit id to use")

# Enable ExternalProject CMake module
include(FetchContent)

FetchContent_Declare(
    autopasfetch
    GIT_REPOSITORY ${autopasRepoPath}
    GIT_TAG master
)

# Get autopas source and binary directories from CMake project
FetchContent_GetProperties(autopasfetch)

if (NOT autopasfetch_POPULATED)
	FetchContent_Populate(autopasfetch)
	add_subdirectory(${autopasfetch_SOURCE_DIR} ${autopasfetch_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()

# Linking AutoPas
set(AUTOPAS_LIB "autopas")
