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
    std::vector<TerrainPoint> terrain2 = server.getTerrainData(bbox2, 25, 25);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Time: " << duration.count() << " ms (should be much faster!)\n";
    printStats(terrain2);

    std::cout << "\n\nStep 3: Fetch adjacent dataset (partial overlap)\n";
    std::cout << "-------------------------------------------------\n";
    
    BoundingBox bbox3;
    bbox3.min_lat = 36.2;   // Partially overlaps
    bbox3.max_lat = 36.5;
    bbox3.min_lon = -112.3;
    bbox3.max_lon = -112.0;
    
    start = std::chrono::high_resolution_clock::now();
    std::vector<TerrainPoint> terrain3 = server.getTerrainData(bbox3, 30, 30);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Time: " << duration.count() << " ms (partial interpolation + API fetch)\n";
    printStats(terrain3);

    std::cout << "\n\nStep 4: Fetch non-overlapping dataset (different region)\n";
    std::cout << "---------------------------------------------------------\n";
    
    BoundingBox bbox4;
    bbox4.min_lat = 37.7;   // Yosemite area
    bbox4.max_lat = 38.0;
    bbox4.min_lon = -119.7;
    bbox4.max_lon = -119.4;
    
    start = std::chrono::high_resolution_clock::now();
    std::vector<TerrainPoint> terrain4 = server.getTerrainData(bbox4, 20, 20);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Time: " << duration.count() << " ms (new region, full API fetch)\n";
    printStats(terrain4);

    std::cout << "\n\nStep 5: Re-fetch Step 2 dataset (should be instant from exact cache)\n";
    std::cout << "---------------------------------------------------------------------\n";
    
    start = std::chrono::high_resolution_clock::now();
    std::vector<TerrainPoint> terrain5 = server.getTerrainData(bbox2, 25, 25);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Time: " << duration.count() << " ms (exact cache match)\n";
    printStats(terrain5);

    std::cout << "\n\nStep 6: High-resolution query within cached area\n";
    std::cout << "-------------------------------------------------\n";
    
    BoundingBox bbox6;
    bbox6.min_lat = 36.10;  // Small area within bbox1
    bbox6.max_lat = 36.15;
    bbox6.min_lon = -112.20;
    bbox6.max_lon = -112.15;
    
    start = std::chrono::high_resolution_clock::now();
    std::vector<TerrainPoint> terrain6 = server.getTerrainData(bbox6, 50, 50);  // High resolution
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Time: " << duration.count() << " ms (2500 points via interpolation!)\n";
    printStats(terrain6);

    // Export some results
    server.exportToCSV(terrain6, "high_res_interpolated.csv");

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