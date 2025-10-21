cmake_minimum_required(VERSION 3.16)
project(TerrainServer LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Find required packages
find_package(CURL REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)

# Source files
set(SOURCES
    src/TerrainServerOpenTopo.cpp
    src/main.cpp
)

# Create executable
add_executable(terrain_server ${SOURCES})

# Link libraries
target_link_libraries(terrain_server PRIVATE CURL::libcurl nlohmann_json::nlohmann_json pthread)

# Optional compiler warnings
if (MSVC)
    target_compile_options(terrain_server PRIVATE /W4)
else()
    target_compile_options(terrain_server PRIVATE -Wall -Wextra -Wpedantic)
endif()
