# Enable ExternalProject CMake module
include(FetchContent)

# Select https (default) or ssh path.
set(autopasRepoPath https://github.com/AutoPas/AutoPas.git)

set(AUTOPAS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(AUTOPAS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(AUTOPAS_ENABLE_ADDRESS_SANITIZER ${ENABLE_ADDRESS_SANITIZER} CACHE BOOL "" FORCE)
set(AUTOPAS_OPENMP ${OPENMP} CACHE BOOL "" FORCE)
set(spdlog_ForceBundled ON CACHE BOOL "" FORCE)
set(Eigen3_ForceBundled ON CACHE BOOL "" FORCE)

# Version from March 19th, 2025. Includes PMT energy measurements, dynamic rebuilding and 3-body interactions.
set(AUTOPAS_TAG 0338ff17d55cc770a403efee69740a69753d64ad CACHE STRING "AutoPas Git tag or commit id to use")

# Fetch external library from GitHub
FetchContent_Declare(
  autopasfetch
  GIT_REPOSITORY ${autopasRepoPath}
  GIT_TAG ${AUTOPAS_TAG}
)

# Download and configure the fetched content
FetchContent_MakeAvailable(autopasfetch)


