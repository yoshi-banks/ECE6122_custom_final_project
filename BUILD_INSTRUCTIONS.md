# PX4/OpenGL/A* Pathfinding Project

## Build

### Prerequisites
Before building the project, ensure you have the following prerequisites installed:
- A compatible operating system (Only tested on Ubuntu 22.04 and WSL2 with Ubuntu 22.04)
- Git
- CMake
- A C++ compiler (GCC 7 or later)
- OpenGL development libraries
- PX4 Development Environment (refer to the PX4 Dev Guide)

#### WSL 

This was tested on Windows Subsystem for Linux (WSL) version 2 with Ubuntu 22.04. Make sure to install WSL2 and set up Ubuntu from the Microsoft Store. Follow the official Microsoft documentation for detailed instructions on setting up WSL2.

### Setup PX4
To set up the PX4 development environment, follow these steps:
1. **Clone the PX4 Firmware Repository**:
    - Tested with the v1.16.0 release.
2. **Install Required Dependencies**:
    - Follow the instructions on the [PX4 Dev Guide](https://dev.px4.io/master/en/setup/dev_env.html) to install all necessary dependencies for your operating system.
3. **Set Up the Build Environment**:
    - Source the PX4 environment setup script:
      ```bash
      source ~/px4_firmware/Tools/setup/ubuntu.sh
      ```

### Build cpp software
To build the C++ software for the PX4/OpenGL/A* Pathfinding project, follow these steps:
1. **Navigate to the Project Directory**:
2. **Install Project Dependencies**:
    - Use the provided `install_dependencies.sh` script to install all required libraries and dependencies:
      ```bash
      ./install_dependencies.sh
      ```
3. **Build the Project**:
    - Run the build.sh script provided in the project:

The executables will be generated in the `build/output/bin` directory.

### Install QGroundControl
To install QGroundControl, follow these steps:
1. **Download QGroundControl**:
    - Visit the [QGroundControl Download Page](https://docs.qgroundcontrol.com/en/getting_started/download_and_install.html) and download the appropriate version for your operating system.
2. **Install QGroundControl**:
    - Follow the installation instructions specific to your operating system provided on the download page.

## Run

To run the PX4/OpenGL/A* Pathfinding project, follow these steps:
1. **Start the PX4 Simulation**:
    - Navigate to the project directory and run the simulation script:
      ```bash
      ./launch_px4.sh
      ```

      This assumes that you cloned PX4-Autopilot in the project directory. I.e. the directory structure should look like this:
      ```
      /path/to/project/
      ├── px4-Autopilot
      └── other_project_files
      ```
2. **Launch QGroundControl**:
    - Open QGroundControl from your applications menu or by running the executable you installed earlier.
    - Connect to the PX4 simulation by selecting the appropriate connection method (usually UDP, see QGroundControl documentation for details).

3. **Run the OpenGL Visualization**:
    - In a separate terminal, navigate to the and run the main executable for OpenGL visualization:
      ```bash
      ./build/output/bin/main
      ```

      Note: the project will cache topography files in terrain_cache in the directory you call the binary from. Make sure you have write permissions to this directory. And call from the same directory if repeating runs to take advantage of the cache. 