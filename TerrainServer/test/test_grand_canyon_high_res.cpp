#include <iostream>
#include <chrono>
#include "TerrainServerOpenTopo.hpp"

void printStats(const std::vector<TerrainPoint>& terrain)
{
    if (terrain.empty()) return;
    
    double min_elev = terrain[0].alt;
    double max_elev = terrain[0].alt;
    double sum_elev = 0.0;
    
    for (const auto& pt : terrain)
    {
        if (pt.alt < min_elev) min_elev = pt.alt;
        if (pt.alt > max_elev) max_elev = pt.alt;
        sum_elev += pt.alt;
    }
    
    std::cout << "  Points: " << terrain.size() << "\n";
    std::cout << "  Elevation: " << min_elev << " to " << max_elev << " m\n";
    std::cout << "  Average: " << (sum_elev / terrain.size()) << " m\n";
}

int main()
{
    std::cout << "=== Smart Cache & Interpolation Test ===\n\n";

    TerrainServerOpenTopo server;
    
    // Set interpolation tolerance (0.01 degrees ≈ 1.1km)
    server.setInterpolationTolerance(0.01);

    std::cout << "Step 1: Fetch initial dataset (Grand Canyon area)\n";
    std::cout << "------------------------------------------------\n";
    
    BoundingBox bbox1;
    bbox1.min_lat = 36.0;
    bbox1.max_lat = 36.3;
    bbox1.min_lon = -112.3;
    bbox1.max_lon = -112.0;
    
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<TerrainPoint> terrain1 = server.getTerrainData(bbox1, 30, 30);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Time: " << duration.count() << " ms\n";
    printStats(terrain1);

    std::cout << "\n\nStep 2: Fetch overlapping dataset (should use interpolation)\n";
    std::cout << "-------------------------------------------------------------\n";
    
    BoundingBox bbox2;
    bbox2.min_lat = 36.05;  // Overlaps with bbox1
    bbox2.max_lat = 36.25;
    bbox2.min_lon = -112.25;
    bbox2.max_lon = -112.05;
    
    start = std::chrono::high_resolution_clock::now();
    std::vector<TerrainPoint> terrain2 = server.getTerrainData(bbox2, 500, 500);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Time: " << duration.count() << " ms (should be much faster!)\n";
    printStats(terrain2);

    // Export some results
    server.exportToCSV(terrain2, "high_res_grandcanyon_interpolated.csv");

    std::cout << "\n\n=== Summary ===\n";
    std::cout << "The smart caching system:\n";
    std::cout << "1. Loads all existing cache files on startup\n";
    std::cout << "2. Checks for exact cache matches first\n";
    std::cout << "3. Interpolates from nearby cached data when possible\n";
    std::cout << "4. Only fetches missing data from the API\n";
    std::cout << "5. Saves new data to cache for future use\n";
    std::cout << "\nThis dramatically reduces API calls and improves performance!\n";

    std::cout << "\n=== Test Complete ===\n";

    return 0;
}