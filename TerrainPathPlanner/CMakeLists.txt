cmake_minimum_required(VERSION 3.16)
project(TerrainPathPlanner LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ============================================================================
# MPI Support (Optional)
# ============================================================================
option(USE_MPI "Enable MPI parallelization" OFF)

if(USE_MPI)
    find_package(MPI REQUIRED)
    message(STATUS "MPI enabled - building with parallel support")
    add_definitions(-DUSE_MPI)
else()
    message(STATUS "MPI disabled - building single-threaded version")
endif()

# ============================================================================
# Source Files
# ============================================================================
set(PLANNER_SOURCES
    src/TerrainPathPlanner.cpp
)

# Add MPI implementation only if MPI is enabled
if(USE_MPI)
    list(APPEND PLANNER_SOURCES 
        src/TerrainPathPlannerMPI.cpp
        src/TerrainPathPlanner_MPIStubs.cpp
    )
endif()

# ============================================================================
# Main Library Target
# ============================================================================
add_library(terrain_path_planner_lib ${PLANNER_SOURCES})

# Include directories
target_include_directories(terrain_path_planner_lib PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/../TerrainServer/include
)

# Link libraries
target_link_libraries(terrain_path_planner_lib PUBLIC
    pthread
)

# Link MPI if enabled
if(USE_MPI)
    target_link_libraries(terrain_path_planner_lib PUBLIC
        MPI::MPI_CXX
    )
    target_compile_definitions(terrain_path_planner_lib PUBLIC USE_MPI)
endif()

# Optional compiler warnings
if(MSVC)
    target_compile_options(terrain_path_planner_lib PRIVATE /W4)
else()
    target_compile_options(terrain_path_planner_lib PRIVATE -Wall -Wextra -Wpedantic)
endif()

# ============================================================================
# Installation (optional)
# ============================================================================
install(TARGETS terrain_path_planner_lib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)