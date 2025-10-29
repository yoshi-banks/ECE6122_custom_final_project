#!/usr/bin/env python3
"""
GPS to Gazebo World Converter
Converts JSON GPS data (lat, lon, alt) to Gazebo world with terrain heightmap
"""

import json
import math
import numpy as np
from PIL import Image
import argparse
import os


def gps_to_local(lat, lon, alt, origin_lat, origin_lon, origin_alt):
    """
    Convert GPS coordinates to local ENU (East-North-Up) coordinates
    
    Args:
        lat, lon, alt: Target GPS coordinates
        origin_lat, origin_lon, origin_alt: Origin GPS coordinates
    
    Returns:
        (x, y, z) in meters relative to origin
    """
    # Earth radius in meters
    R = 6378137.0
    
    # Convert to radians
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)
    origin_lat_rad = math.radians(origin_lat)
    origin_lon_rad = math.radians(origin_lon)
    
    # Calculate local coordinates (ENU)
    x = R * (lon_rad - origin_lon_rad) * math.cos(origin_lat_rad)  # East
    y = R * (lat_rad - origin_lat_rad)  # North
    z = alt - origin_alt  # Up
    
    return x, y, z


def create_heightmap_from_points(points, resolution=512, interpolation='linear'):
    """
    Create a heightmap image from GPS points
    
    Args:
        points: List of local coordinate dicts with 'x', 'y', 'z'
        resolution: Output image resolution (pixels)
        interpolation: 'linear' or 'nearest'
    
    Returns:
        numpy array representing heightmap, bounds info
    """
    if not points:
        raise ValueError("No points provided")
    
    # Extract coordinates
    xs = [p['x'] for p in points]
    ys = [p['y'] for p in points]
    zs = [p['z'] for p in points]
    
    # Get bounds
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    z_min, z_max = min(zs), max(zs)
    
    # Add padding
    x_range = x_max - x_min
    y_range = y_max - y_min
    padding = 0.1  # 10% padding
    
    x_min -= x_range * padding
    x_max += x_range * padding
    y_min -= y_range * padding
    y_max += y_range * padding
    
    # Create grid
    heightmap = np.zeros((resolution, resolution))
    
    # Simple grid interpolation
    for i in range(resolution):
        for j in range(resolution):
            # Map grid position to world coordinates
            x = x_min + (x_max - x_min) * j / (resolution - 1)
            y = y_min + (y_max - y_min) * i / (resolution - 1)
            
            # Find nearest points or interpolate
            if interpolation == 'nearest':
                # Find nearest point
                distances = [(p['x'] - x)**2 + (p['y'] - y)**2 for p in points]
                nearest_idx = distances.index(min(distances))
                heightmap[i, j] = points[nearest_idx]['z']
            else:
                # Inverse distance weighting
                total_weight = 0
                weighted_z = 0
                for p in points:
                    dist = math.sqrt((p['x'] - x)**2 + (p['y'] - y)**2)
                    if dist < 0.01:  # Very close
                        weighted_z = p['z']
                        total_weight = 1
                        break
                    weight = 1 / (dist + 1)  # Inverse distance
                    weighted_z += weight * p['z']
                    total_weight += weight
                
                if total_weight > 0:
                    heightmap[i, j] = weighted_z / total_weight
    
    # Normalize to 0-1 range
    if z_max - z_min > 0:
        heightmap = (heightmap - z_min) / (z_max - z_min)
    
    bounds = {
        'x_min': x_min,
        'x_max': x_max,
        'y_min': y_min,
        'y_max': y_max,
        'z_min': z_min,
        'z_max': z_max,
        'x_size': x_max - x_min,
        'y_size': y_max - y_min,
        'z_size': z_max - z_min
    }
    
    return heightmap, bounds


def save_heightmap_image(heightmap, output_path):
    """Save heightmap as 16-bit grayscale PNG"""
    # Convert to 16-bit for better precision
    heightmap_16bit = (heightmap * 65535).astype(np.uint16)
    img = Image.fromarray(heightmap_16bit, mode='I;16')
    img.save(output_path)
    print(f"Heightmap saved to: {output_path}")


def create_gazebo_world(world_name, heightmap_filename, bounds, origin_gps, output_path):
    """Create Gazebo SDF world file"""
    
    x_size = bounds['x_size']
    y_size = bounds['y_size']
    z_size = bounds['z_size']
    
    # Use larger of x or y for square terrain
    terrain_size = max(x_size, y_size) * 1.2  # Add 20% buffer
    
    world_content = f"""<?xml version="1.0" ?>
<sdf version="1.6">
  <world name="{world_name}">
    
    <!-- Physics -->
    <physics type="ode">
      <max_step_size>0.001</max_step_size>
      <real_time_factor>1.0</real_time_factor>
      <real_time_update_rate>1000</real_time_update_rate>
    </physics>
    
    <!-- Lighting -->
    <include>
      <uri>model://sun</uri>
    </include>
    
    <light name="sun_light" type="directional">
      <cast_shadows>true</cast_shadows>
      <pose>0 0 10 0 0 0</pose>
      <diffuse>0.8 0.8 0.8 1</diffuse>
      <specular>0.2 0.2 0.2 1</specular>
      <attenuation>
        <range>1000</range>
        <constant>0.9</constant>
        <linear>0.01</linear>
        <quadratic>0.001</quadratic>
      </attenuation>
      <direction>-0.5 0.1 -0.9</direction>
    </light>
    
    <!-- Terrain Model -->
    <model name="terrain">
      <static>true</static>
      <link name="link">
        <collision name="collision">
          <geometry>
            <heightmap>
              <uri>file://{heightmap_filename}</uri>
              <size>{terrain_size:.2f} {terrain_size:.2f} {z_size:.2f}</size>
              <pos>0 0 {bounds['z_min']:.2f}</pos>
            </heightmap>
          </geometry>
        </collision>
        <visual name="visual">
          <geometry>
            <heightmap>
              <uri>file://{heightmap_filename}</uri>
              <size>{terrain_size:.2f} {terrain_size:.2f} {z_size:.2f}</size>
              <pos>0 0 {bounds['z_min']:.2f}</pos>
              <texture>
                <diffuse>file://media/materials/textures/dirt_diffusespecular.png</diffuse>
                <normal>file://media/materials/textures/flat_normal.png</normal>
                <size>10</size>
              </texture>
              <blend>
                <min_height>{bounds['z_min']:.2f}</min_height>
                <fade_dist>5</fade_dist>
              </blend>
            </heightmap>
          </geometry>
        </visual>
      </link>
    </model>
    
    <!-- Ground Plane (below terrain) -->
    <model name="ground_plane">
      <static>true</static>
      <link name="link">
        <collision name="collision">
          <geometry>
            <plane>
              <normal>0 0 1</normal>
              <size>10000 10000</size>
            </plane>
          </geometry>
        </collision>
        <visual name="visual">
          <geometry>
            <plane>
              <normal>0 0 1</normal>
              <size>10000 10000</size>
            </plane>
          </geometry>
          <material>
            <ambient>0.5 0.5 0.5 1</ambient>
          </material>
        </visual>
      </link>
    </model>
    
    <!-- Scene settings -->
    <scene>
      <ambient>0.4 0.4 0.4 1</ambient>
      <background>0.7 0.7 0.7 1</background>
      <shadows>true</shadows>
    </scene>
    
  </world>
</sdf>
"""
    
    with open(output_path, 'w') as f:
        f.write(world_content)
    
    print(f"Gazebo world file saved to: {output_path}")


def create_px4_launch_script(world_name, origin_gps, output_path):
    """Create a bash script to launch PX4 SITL with the world"""
    
    script_content = f"""#!/bin/bash
# Launch script for PX4 SITL with custom terrain

# Set home position (origin of your GPS data)
export PX4_HOME_LAT={origin_gps['lat']}
export PX4_HOME_LON={origin_gps['lon']}
export PX4_HOME_ALT={origin_gps['alt']}

# Launch PX4 SITL with Gazebo
# Make sure you're in your PX4-Autopilot directory
# and the world file is in Tools/simulation/gazebo-classic/sitl_gazebo-classic/worlds/

echo "Starting PX4 SITL with custom terrain..."
echo "Home Position: LAT=$PX4_HOME_LAT, LON=$PX4_HOME_LON, ALT=$PX4_HOME_ALT"

# Option 1: If world is in the worlds directory
make px4_sitl gazebo-classic_iris__{world_name}

# Option 2: Specify world file directly
# GAZEBO_WORLD_PATH=/path/to/your/{world_name}.world make px4_sitl gazebo-classic_iris
"""
    
    with open(output_path, 'w') as f:
        f.write(script_content)
    
    # Make executable
    os.chmod(output_path, 0o755)
    print(f"PX4 launch script saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Convert GPS JSON data to Gazebo world for PX4 drone simulation'
    )
    parser.add_argument(
        'json_file',
        help='Path to JSON file containing GPS points'
    )
    parser.add_argument(
        '-o', '--output-dir',
        default='gazebo_world',
        help='Output directory for generated files (default: gazebo_world)'
    )
    parser.add_argument(
        '-n', '--name',
        default='custom_terrain',
        help='Name for the world (default: custom_terrain)'
    )
    parser.add_argument(
        '-r', '--resolution',
        type=int,
        default=512,
        help='Heightmap resolution in pixels (default: 512)'
    )
    parser.add_argument(
        '--interpolation',
        choices=['linear', 'nearest'],
        default='linear',
        help='Interpolation method (default: linear)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load JSON data
    print(f"Loading GPS data from: {args.json_file}")
    with open(args.json_file, 'r') as f:
        data = json.load(f)
    
    if 'points' not in data:
        raise ValueError("JSON must contain 'points' array")
    
    gps_points = data['points']
    print(f"Loaded {len(gps_points)} GPS points")
    
    # Use first point as origin
    origin = gps_points[0]
    origin_gps = {
        'lat': origin['lat'],
        'lon': origin['lon'],
        'alt': origin['alt']
    }
    
    print(f"\nOrigin GPS: LAT={origin_gps['lat']}, LON={origin_gps['lon']}, ALT={origin_gps['alt']}m")
    
    # Convert all points to local coordinates
    print("\nConverting GPS to local coordinates...")
    local_points = []
    for i, p in enumerate(gps_points):
        x, y, z = gps_to_local(
            p['lat'], p['lon'], p['alt'],
            origin_gps['lat'], origin_gps['lon'], origin_gps['alt']
        )
        local_points.append({'x': x, 'y': y, 'z': z})
        if i < 5 or i == len(gps_points) - 1:  # Show first 5 and last
            print(f"  Point {i+1}: GPS({p['lat']}, {p['lon']}, {p['alt']}m) -> Local({x:.2f}, {y:.2f}, {z:.2f}m)")
        elif i == 5:
            print(f"  ... ({len(gps_points) - 6} more points) ...")
    
    # Create heightmap
    print(f"\nGenerating {args.resolution}x{args.resolution} heightmap...")
    heightmap, bounds = create_heightmap_from_points(
        local_points,
        resolution=args.resolution,
        interpolation=args.interpolation
    )
    
    print(f"Terrain bounds:")
    print(f"  X: {bounds['x_min']:.2f}m to {bounds['x_max']:.2f}m ({bounds['x_size']:.2f}m)")
    print(f"  Y: {bounds['y_min']:.2f}m to {bounds['y_max']:.2f}m ({bounds['y_size']:.2f}m)")
    print(f"  Z: {bounds['z_min']:.2f}m to {bounds['z_max']:.2f}m ({bounds['z_size']:.2f}m)")
    
    # Save heightmap
    heightmap_filename = f"{args.name}_heightmap.png"
    heightmap_path = os.path.join(args.output_dir, heightmap_filename)
    save_heightmap_image(heightmap, heightmap_path)
    
    # Create Gazebo world file
    world_filename = f"{args.name}.world"
    world_path = os.path.join(args.output_dir, world_filename)
    create_gazebo_world(args.name, heightmap_filename, bounds, origin_gps, world_path)
    
    # Create launch script
    launch_script_path = os.path.join(args.output_dir, f"launch_{args.name}.sh")
    create_px4_launch_script(args.name, origin_gps, launch_script_path)
    
    # Save local coordinates for reference
    local_coords_path = os.path.join(args.output_dir, f"{args.name}_local_coords.json")
    with open(local_coords_path, 'w') as f:
        json.dump({
            'origin_gps': origin_gps,
            'points': local_points,
            'bounds': bounds
        }, f, indent=2)
    print(f"Local coordinates saved to: {local_coords_path}")
    
    print("\n" + "="*60)
    print("SUCCESS! Files generated:")
    print("="*60)
    print(f"1. Heightmap: {heightmap_path}")
    print(f"2. World file: {world_path}")
    print(f"3. Launch script: {launch_script_path}")
    print(f"4. Local coords: {local_coords_path}")
    
    print("\n" + "="*60)
    print("Next Steps:")
    print("="*60)
    print("1. Copy the world file to your PX4-Autopilot worlds directory:")
    print(f"   cp {world_path} ~/PX4-Autopilot/Tools/simulation/gazebo-classic/sitl_gazebo-classic/worlds/")
    print()
    print("2. Copy the heightmap to the same directory:")
    print(f"   cp {heightmap_path} ~/PX4-Autopilot/Tools/simulation/gazebo-classic/sitl_gazebo-classic/worlds/")
    print()
    print("3. Run the launch script from your PX4-Autopilot directory:")
    print(f"   cd ~/PX4-Autopilot")
    print(f"   {launch_script_path}")
    print()
    print("Or manually launch with:")
    print(f"   export PX4_HOME_LAT={origin_gps['lat']}")
    print(f"   export PX4_HOME_LON={origin_gps['lon']}")
    print(f"   export PX4_HOME_ALT={origin_gps['alt']}")
    print(f"   make px4_sitl gazebo-classic_iris__{args.name}")
    print("="*60)


if __name__ == '__main__':
    main()