#pragma once

#include <string>
#include <vector>

struct TerrainPoint 
{
    double lat;
    double lon;
    double elevation;
};

class ITerrainServer
{
    public:
        virtual ~ITerrainServer() = default;

        // Query elevation for a single point
        virtual double getElevation(double lat, double lon) = 0;

        // Query multiple elevations (optional override)
        virtual std::vector<TerrainPoint> getElevationGrid(double latStart, 
            double lonStart, double latEnd, double lonEnd, int numLatSamples, int numLonSamples) = 0;

        // Optional: Start a TCP server to respond to requests
        virtual void startServer(int port) = 0;
};