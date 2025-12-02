/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file main.cpp
 * @brief Main function for testing the TerrainServerOpenTopo and TerrainServerGeoTIFF classes
 */

#include <iostream>

#include "TerrainServerOpenTopo.hpp"

/**
 * @brief Main function
 */
int main()
{
    // TerrainServerOpenTopo server;
    TerrainServerGeoTIFF server;

    BoundingBox bbox;
    bbox.min_lat = 35.991364;
    bbox.max_lat = 36.332758;
    bbox.min_lon = -112.289476;
    bbox.max_lon = -111.894833;

    // Get terrain data (50x50 grid)
    std::vector<TerrainPoint> terrain = server.getTerrainData(bbox, 50, 50);

    std::cout << "\nRetrieved " << terrain.size() << " terrain points" << std::endl;

    // Print first few points
    std::cout << "\nFirst 5 points:" << std::endl;
    for (int i = 0; i < std::min(5, (int)terrain.size()); ++i)
    {
        std::cout << "Lat: " << terrain[i].lat
                    << ", Lon: " << terrain[i].lon
                    << ", Elevation: " << terrain[i].alt << " m" << std::endl;
    }

    // Export to CSV
    server.exportToCSV(terrain, "terrain_data.csv");

    return 0;
}