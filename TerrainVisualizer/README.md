# Terrain Visualizer

A modern OpenGL application for visualizing real-world terrain data fetched from OpenTopoData API.

## Features

### 🎨 Rendering

- **Modern OpenGL 3.3+** with shaders
- **Phong lighting** for realistic terrain shading
- **Height-based coloring**: Valleys (green) → Hills (brown) → Mountains (gray) → Peaks (white/snow)
- **Smooth normals** for natural terrain appearance
- **4x MSAA** for anti-aliasing

### 🎮 Camera Controls

- **Left Mouse + Drag**: Rotate camera around terrain
- **Right Mouse + Drag**: Pan camera
- **Mouse Wheel**: Zoom in/out
- **ESC**: Exit application

### ⚡ Performance

- **Smart caching** via TerrainServer
- **Efficient mesh generation** with indexed triangles
- **Vertex Array Objects (VAO)** for fast rendering
- **60+ FPS** on modern hardware

## Building

### Dependencies

```bash
# Ubuntu/Debian
sudo apt-get install libglew-dev libglfw3-dev libglm-dev

# macOS
brew install glew glfw glm
```

### Compile

```bash
cd build
cmake ..
make
```

## Usage

### Basic Usage

```bash
./output/bin/terrain_visualizer
```

Loads Grand Canyon area by default (36.0°N to 36.3°N, -112.3°W to -112.0°W)

### Custom Region

```bash
./output/bin/terrain_visualizer <latStart> <lonStart> <latEnd> <lonEnd> [rows] [cols]
```

### Examples

**Mount Everest Region:**

```bash
./output/bin/terrain_visualizer 27.8 86.7 28.1 87.0 60 60
```

**Yosemite National Park:**

```bash
./output/bin/terrain_visualizer 37.7 -119.7 38.0 -119.4 50 50
```

**Grand Canyon (High Resolution):**

```bash
./output/bin/terrain_visualizer 36.0 -112.3 36.3 -112.0 100 100
```

**Death Valley:**

```bash
./output/bin/terrain_visualizer 36.2 -117.1 36.6 -116.7 50 50
```

## Architecture

### Class Structure

```
Camera
├── Orbit camera around terrain
├── Mouse-based rotation and panning
└── Smooth zoom control

TerrainMesh
├── Generates mesh from TerrainServer data
├── Calculates smooth normals
├── Height-based coloring
└── OpenGL buffer management (VAO/VBO/EBO)

Shader
├── Vertex and fragment shaders
├── Phong lighting model
└── Uniform management

main.cpp
├── Window and input management
├── Render loop
└── Integrates all components
```

### Rendering Pipeline

```
TerrainServer → TerrainPoint[] → Mesh Generation
                                        ↓
                                  Vertices + Indices
                                        ↓
                                   GPU Buffers
                                        ↓
                         Vertex Shader (Transform + Lighting)
                                        ↓
                          Fragment Shader (Color + Shading)
                                        ↓
                                   Screen Display
```

## Color Scheme

The terrain uses realistic elevation-based coloring:

| Height  | Color                    | Represents        |
| ------- | ------------------------ | ----------------- |
| 0-20%   | Dark Green → Light Green | Valleys, lowlands |
| 20-50%  | Green → Brown            | Hills, foothills  |
| 50-70%  | Brown → Gray             | Mountains         |
| 70-100% | Gray → White             | High peaks, snow  |

## Performance Tips

### For Large Areas

1. **Start with coarse grids** (20x20 or 30x30) to build cache
2. **Increase resolution** gradually (50x50, then 100x100)
3. **Use smart caching** - overlapping regions use interpolation

### For Best Quality

- Use **higher grid resolution** (100x100 or more)
- Smaller geographic areas give better detail
- Mountainous regions benefit from denser grids

### Grid Size Guidelines

- **20x20**: Quick preview (400 points, ~1 second)
- **50x50**: Good balance (2,500 points, ~3-5 seconds first time)
- **100x100**: High detail (10,000 points, ~10-15 seconds first time)
- **200x200**: Maximum quality (40,000 points, ~30-60 seconds first time)

_Note: Subsequent requests use cache and are nearly instant!_

## Keyboard Shortcuts (Future Enhancement)

Future versions may include:

- **W**: Wireframe toggle
- **L**: Lighting toggle
- **C**: Color scheme cycle
- **R**: Reset camera
- **F**: Fullscreen toggle
- **Space**: Pause/resume animation

## Technical Details

### Vertex Format

```cpp
struct Vertex {
    glm::vec3 position;  // World coordinates
    glm::vec3 normal;    // Surface normal for lighting
    glm::vec3 color;     // Height-based color
};
```

### Lighting Model

Uses **Phong reflection model**:

- **Ambient**: 30% base illumination
- **Diffuse**: Angle-based lighting
- **Specular**: 20% highlights on peaks

### Mesh Generation

1. Fetch elevation data from TerrainServer
2. Normalize elevations to [0, 10] range
3. Generate grid of vertices
4. Create triangle indices (2 triangles per quad)
5. Calculate per-vertex normals
6. Upload to GPU

## Troubleshooting

### Black Screen

- Check OpenGL version: `glxinfo | grep "OpenGL version"`
- Requires OpenGL 3.3 or higher

### Slow Performance

- Reduce grid resolution (use smaller rows/cols)
- Check terrain cache is being used (should be fast on 2nd run)
- Update graphics drivers

### Missing Terrain

- Check internet connection (first fetch requires API access)
- Verify coordinates are valid
- Check terrain_cache/ directory for cached data

### Compilation Errors

```bash
# Missing GLEW
sudo apt-get install libglew-dev

# Missing GLFW
sudo apt-get install libglfw3-dev

# Missing GLM
sudo apt-get install libglm-dev
```

## Examples Gallery

### Grand Canyon

```bash
./terrain_visualizer 36.0 -112.3 36.3 -112.0 80 80
```

Features: Deep canyon walls, plateaus, dramatic elevation changes

### Himalayas (Mount Everest)

```bash
./terrain_visualizer 27.8 86.7 28.1 87.0 60 60
```

Features: Extreme elevations, snow-capped peaks, valleys

### Yosemite Valley

```bash
./terrain_visualizer 37.7 -119.7 38.0 -119.4 70 70
```

Features: Granite cliffs, valleys, Half Dome

### Hawaiian Volcanoes

```bash
./terrain_visualizer 19.3 -155.6 19.6 -155.2 60 60
```

Features: Volcanic cones, lava flows, coastal slopes

## Future Enhancements

- [ ] Texture mapping support
- [ ] Water plane rendering
- [ ] Fog effects for distance
- [ ] Shadow mapping
- [ ] Export to OBJ/STL formats
- [ ] Multiple terrain overlays
- [ ] Vegetation/land cover data
- [ ] Real-time terrain editing
- [ ] VR support

## Credits

- **OpenTopoData API**: Free elevation data source
- **SRTM Dataset**: NASA Shuttle Radar Topography Mission
- **Libraries**: GLFW, GLEW, GLM
