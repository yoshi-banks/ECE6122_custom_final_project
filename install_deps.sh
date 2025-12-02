#!/bin/bash

set -e

# Update package lists
sudo apt update

# Install build tools and CMake >= 3.16 (cmake from apt is 3.22+ on 22.04)
sudo apt install -y build-essential cmake

# Install OpenGL dependencies
sudo apt install -y libgl1-mesa-dev libglu1-mesa-dev

# Install GLEW
sudo apt install -y libglew-dev

# Install GLFW3
sudo apt install -y libglfw3-dev

# Install GLM
sudo apt install -y libglm-dev

# Install X11 dependencies
sudo apt install -y libx11-dev libxrandr-dev libxcursor-dev libxi-dev libxinerama-dev libxxf86vm-dev

# Install CURL
sudo apt install -y libcurl4-openssl-dev

# Install pthreads (provided by libc6-dev)
sudo apt install -y libc6-dev

# Install nlohmann_json (>= 3.2.0)
sudo apt install -y nlohmann-json3-dev

echo "All required development packages have been installed."
