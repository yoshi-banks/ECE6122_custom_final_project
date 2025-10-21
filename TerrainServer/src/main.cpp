#include <iostream>

#include "TerrainServerOpenTopo.hpp"

int main()
{
    TerrainServerOpenTopo topoServer("srtm90m");

    // Example 1: Single Point
    double elev = topoServer.getElevation(34.012, -84.321);
    std::cout << "Elevation at (34.012, -84.321): " << elev << " meters\n";

    // Example 2: Start a TCP socket server
    topoServer.startServer(8080);

    return 0;
}