#pragma once

#include <string>
#include <mutex>

#include "ITerrainServer.hpp"

class TerrainServerOpenTopo : public ITerrainServer
{
    public:
        TerrainServerOpenTopo(const std::string& = "srtm90m");
        ~TerrainServerOpenTopo() override = default;

        double getElevation(double lat, double lon) override;
        std::vector<TerrainPoint> getElevationGrid(double latStart,
            double lonStart, double latEnd, double lonEnd, int numLatSamples, int numLonSamples) override;

        void startServer(int port) override;

    private:
        std::string queryAPI(const std::string& url);
        std::string dataset_;
        std::mutex ioMutex_;
};