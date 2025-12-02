# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(project VERSION 1.0)

# Set policy for OpenGL to suppress warning
set(OpenGL_GL_PREFERENCE GLVND)

# ============================================================================
# MPI Support Option (OFF by default)
# ============================================================================
option(USE_MPI "Enable MPI parallelization for path planning" OFF)

if(USE_MPI)
    find_package(MPI REQUIRED)
    message(STATUS "MPI enabled - path planner will use parallel algorithms")
    add_definitions(-DUSE_MPI)
else()
    message(STATUS "MPI disabled - using single-threaded path planner")
endif()

# Find required packages
find_package(OpenGL REQUIRED)
find_package(MAVSDK REQUIRED)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define a common output directory
set(COMMON_OUTPUT_DIR "${CMAKE_BINARY_DIR}/output")
if(WIN32)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${COMMON_OUTPUT_DIR}/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG "${COMMON_OUTPUT_DIR}/lib")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${COMMON_OUTPUT_DIR}/bin")

    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${COMMON_OUTPUT_DIR}/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE "${COMMON_OUTPUT_DIR}/lib")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${COMMON_OUTPUT_DIR}/bin")
elseif(UNIX)
    # Set output directories for different types of binaries
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${COMMON_OUTPUT_DIR}/bin")  # For executables
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${COMMON_OUTPUT_DIR}/lib")  # For shared libraries
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${COMMON_OUTPUT_DIR}/lib")  # For static libraries
endif()

# Don't build GLFW examples, tests, or documentation
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)

add_subdirectory(TerrainPathPlanner)
add_subdirectory(TerrainServer)
add_subdirectory(TerrainVisualizer)

# Source Files
set(MAIN_SOURCES
    main.cpp
    PX4Interface.cpp
)

add_executable(main ${MAIN_SOURCES})

target_link_libraries(main
    terrain_server_lib
    terrain_path_planner_lib
    terrain_visualizer_lib
    MAVSDK::mavsdk
)

# Link MPI if enabled
if(USE_MPI)
    target_link_libraries(main MPI::MPI_CXX)
    target_compile_definitions(main PRIVATE USE_MPI)
endif()

target_include_directories(main PRIVATE
    ${PROJECT_SOURCE_DIR}/TerrainVisualizer/include
    ${PROJECT_SOURCE_DIR}/TerrainPathPlanner/include
    ${PROJECT_SOURCE_DIR}/TerrainServer/include
)

# ============================================================================
# Installation
# ============================================================================
install(TARGETS main
    RUNTIME DESTINATION bin
)

# Print configuration summary
message(STATUS "")
message(STATUS "=== Configuration Summary ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "MPI support: ${USE_MPI}")
if(USE_MPI)
    message(STATUS "MPI C++ compiler: ${MPI_CXX_COMPILER}")
endif()
message(STATUS "Output directory: ${COMMON_OUTPUT_DIR}")
message(STATUS "=============================")
message(STATUS "")