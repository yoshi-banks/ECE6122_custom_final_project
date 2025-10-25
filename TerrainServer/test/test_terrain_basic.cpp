#include <iostream>
#include "TerrainServerOpenTopo.hpp"

int main()
{
    std::cout << "=== Basic Terrain Data Retrieval Test ===\n\n";

    TerrainServerOpenTopo server;

    // Grand Canyon area
    BoundingBox bbox;
    bbox.min_lat = 36.056361;
    bbox.max_lat = 36.137089;
    bbox.min_lon = -112.175607;
    bbox.max_lon = -112.105038;

    // Get terrain data (50x50 grid)
    std::cout << "Requesting 50x50 grid...\n";
    std::vector<TerrainPoint> terrain = server.getTerrainData(bbox, 50, 50);

    std::cout << "\nRetrieved " << terrain.size() << " terrain points\n";

    // Print first few points
    std::cout << "\nFirst 5 points:\n";
    for (int i = 0; i < std::min(5, (int)terrain.size()); ++i)
    {
        std::cout << "  Lat: " << terrain[i].lat
                  << ", Lon: " << terrain[i].lon
                  << ", Elevation: " << terrain[i].alt << " m\n";
    }

    // Calculate statistics
    if (!terrain.empty())
    {
        double min_elev = terrain[0].alt;
        double max_elev = terrain[0].alt;
        double sum_elev = 0.0;

        for (const auto& pt : terrain)
        {
            if (pt.alt < min_elev) min_elev = pt.alt;
            if (pt.alt > max_elev) max_elev = pt.alt;
            sum_elev += pt.alt;
        }

        std::cout << "\nElevation Statistics:\n";
        std::cout << "  Min: " << min_elev << " m\n";
        std::cout << "  Max: " << max_elev << " m\n";
        std::cout << "  Avg: " << (sum_elev / terrain.size()) << " m\n";
    }

    // Export to CSV
    server.exportToCSV(terrain, "test_terrain_data.csv");

    std::cout << "\n=== Test Complete ===\n";

    return 0;
}