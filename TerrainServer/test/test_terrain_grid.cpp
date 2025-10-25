#include <iostream>
#include <chrono>
#include "TerrainServerOpenTopo.hpp"

void testGridSize(TerrainServerOpenTopo& server, const BoundingBox& bbox, 
                  int grid_lat, int grid_lon)
{
    std::cout << "\n--- Testing " << grid_lat << "x" << grid_lon << " grid ---\n";
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<TerrainPoint> terrain = server.getTerrainData(bbox, grid_lat, grid_lon);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Retrieved " << terrain.size() << " points in " 
              << duration.count() << " ms\n";
    
    if (!terrain.empty())
    {
        double min_elev = terrain[0].alt;
        double max_elev = terrain[0].alt;
        
        for (const auto& pt : terrain)
        {
            if (pt.alt < min_elev) min_elev = pt.alt;
            if (pt.alt > max_elev) max_elev = pt.alt;
        }
        
        std::cout << "Elevation range: " << min_elev << " to " << max_elev << " m\n";
    }
}

int main()
{
    std::cout << "=== Grid Performance Test ===\n";
    std::cout << "Testing different grid sizes for caching performance\n";

    TerrainServerOpenTopo server;

    // Yosemite National Park area
    BoundingBox bbox;
    bbox.min_lat = 37.7;
    bbox.max_lat = 38.0;
    bbox.min_lon = -119.7;
    bbox.max_lon = -119.4;

    // Test different grid sizes
    testGridSize(server, bbox, 10, 10);   // Small - 100 points
    testGridSize(server, bbox, 25, 25);   // Medium - 625 points
    testGridSize(server, bbox, 50, 50);   // Large - 2500 points
    
    std::cout << "\n--- Testing cache performance ---\n";
    std::cout << "Re-running 50x50 grid to test cache...\n";
    testGridSize(server, bbox, 50, 50);   // Should be instant from cache

    std::cout << "\n--- Testing getElevationGrid (batched) ---\n";
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<TerrainPoint> grid = server.getElevationGrid(
        37.7, -119.7,  // start lat, lon
        38.0, -119.4,  // end lat, lon
        30, 30         // 30x30 grid
    );
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Retrieved " << grid.size() << " points in " 
              << duration.count() << " ms using getElevationGrid\n";

    std::cout << "\n=== Test Complete ===\n";

    return 0;
}