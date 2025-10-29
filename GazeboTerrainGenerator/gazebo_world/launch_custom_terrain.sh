#!/bin/bash
# Launch script for PX4 SITL with custom terrain

# Set home position (origin of your GPS data)
export PX4_HOME_LAT=37.7
export PX4_HOME_LON=-119.7
export PX4_HOME_ALT=1903.0

# Launch PX4 SITL with Gazebo
# Make sure you're in your PX4-Autopilot directory
# and the world file is in Tools/simulation/gazebo-classic/sitl_gazebo-classic/worlds/

echo "Starting PX4 SITL with custom terrain..."
echo "Home Position: LAT=$PX4_HOME_LAT, LON=$PX4_HOME_LON, ALT=$PX4_HOME_ALT"

# Option 1: If world is in the worlds directory
make px4_sitl gazebo-classic_iris__custom_terrain

# Option 2: Specify world file directly
# GAZEBO_WORLD_PATH=/path/to/your/custom_terrain.world make px4_sitl gazebo-classic_iris
