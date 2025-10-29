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
    src/main.cpp
    src/Camera.cpp
    src/Shader.cpp
    src/TerrainMesh.cpp
)

# Create executable
add_executable(terrain_visualizer ${VISUALIZER_SOURCES})

# Include directories
target_include_directories(terrain_visualizer PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/../TerrainServer/include
)

# Link libraries
target_link_libraries(terrain_visualizer PRIVATE
    terrain_server_lib
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
    target_compile_options(terrain_visualizer PRIVATE /W4)
else()
    target_compile_options(terrain_visualizer PRIVATE -Wall -Wextra -Wpedantic)
endif()

# If building static GLFW (libglfw3.a), explicitly include required system libs:
if (X11_FOUND)
    target_include_directories(terrain_visualizer PRIVATE ${X11_INCLUDE_DIR})
endif()
