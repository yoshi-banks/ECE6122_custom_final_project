cmake_minimum_required(VERSION 3.16)
project(TerrainServer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find packages needed by TerrainServer
find_package(CURL REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(GDAL REQUIRED)

# ============================================================================
# Main Library Target
# ============================================================================
add_library(terrain_server_lib 
    src/ITerrainServer.cpp
    src/Logger.cpp
    src/TerrainServerGeoTIFF.cpp
    src/TerrainServerOpenTopo.cpp
)

# Include directories local to TerrainServer
target_include_directories(terrain_server_lib PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
)

# Link libraries
target_link_libraries(terrain_server_lib PUBLIC
    nlohmann_json::nlohmann_json
    CURL::libcurl
    ${GDAL_LIBRARIES}
    pthread
)

# Optional compiler warnings
if (MSVC)
    target_compile_options(terrain_server_lib PRIVATE /W4)
else()
    target_compile_options(terrain_server_lib PRIVATE -Wall -Wextra -Wpedantic)
endif()

# ============================================================================
# Tests (Optional)
# ============================================================================
option(BUILD_TERRAIN_TESTS "Build terrain server tests" ON)

if(BUILD_TERRAIN_TESTS)
    add_subdirectory(test)
endif()