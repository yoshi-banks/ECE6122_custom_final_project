cmake_minimum_required(VERSION 3.16)
project(TerrainVisualizer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)
find_package(X11 REQUIRED)

# Source files
set(VISUALIZER_SOURCES
    src/Camera.cpp
    src/Cube.cpp
    src/Shader.cpp
    src/TerrainMesh.cpp
)

message(STATUS "VISUALIZER_SOURCES: ${VISUALIZER_SOURCES}")
message(STATUS "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Create library
add_library(terrain_visualizer_lib ${VISUALIZER_SOURCES})

# Include directories
target_include_directories(terrain_visualizer_lib PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

# Link libraries
target_link_libraries(terrain_visualizer_lib PRIVATE
    terrain_server_lib
    terrain_path_planner_lib
    OpenGL::GL
    GLEW::GLEW
    glfw
    glm::glm
    ${CMAKE_DL_LIBS}
    ${X11_LIBRARIES}
    Xrandr
    Xcursor
    Xi
    Xinerama
    Xxf86vm
)

# Optional compiler warnings
if (MSVC)
    target_compile_options(terrain_visualizer_lib PRIVATE /W4)
else()
    target_compile_options(terrain_visualizer_lib PRIVATE -Wall -Wextra -Wpedantic)
endif()

# If building static GLFW (libglfw3.a), explicitly include required system libs:
if (X11_FOUND)
    target_include_directories(terrain_visualizer_lib PRIVATE ${X11_INCLUDE_DIR})
endif()
