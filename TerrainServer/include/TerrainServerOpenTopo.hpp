#pragma once

#include <string>
#include <mutex>

#include "ITerrainServer.hpp"

struct CachedDataset
{
    BoundingBox bbox;
    int grid_lat;
    int grid_lon;
    std::vector<TerrainPoint> points;
    std::string filename;
};

class TerrainServerOpenTopo : public ITerrainServer
{
    public:
        TerrainServerOpenTopo(const std::string& cache_directory = "./terrain_cache",
                              const std::string& api_base_url = "https://api.opentopodata.org/v1/",
                              const std::string& dataset_name = "srtm90m");

        ~TerrainServerOpenTopo() override;
        
        std::vector<TerrainPoint> getTerrainData(const BoundingBox& bbox, 
                                                 int grid_points_lat = 50, 
                                                 int grid_points_lon = 50);
        void exportToCSV(const std::vector<TerrainPoint>& points, const std::string& filename);

        // Set interpolation tolerance (in degrees)
        void setInterpolationTolerance(double tolerance) { interpolationTolerance_ = tolerance; }

        // ITerrainServer interface implementation
        double getElevation(double lat, double lon) override;
        std::vector<TerrainPoint> getElevationGrid(double latStart,
            double lonStart, double latEnd, double lonEnd, 
            int numLatSamples, int numLonSamples) override;
        void startServer(int port) override;

    private:
        std::string cache_dir;
        std::string api_url;
        std::string dataset;
        std::mutex ioMutex_;
        double interpolationTolerance_; // Max distance for interpolation (degrees)

        // Cache management
        std::vector<CachedDataset> loadedCaches_;

        // Callback for CURL
        static size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp);

        // Cache file operations
        std::string getCacheFilename(const BoundingBox& bbox, int grid_points_lat, int grid_points_lon);
        bool loadFromCache(const std::string& filename, std::vector<TerrainPoint>& points);
        void saveToCache(const std::string& filename, const std::vector<TerrainPoint>& points);

        // Multi-cache operations
        std::vector<std::string> findAllCacheFiles();
        void loadAllCaches();
        BoundingBox parseBBoxFromFilename(const std::string& filename);

        // Interpolation functions
        bool canInterpolatePoint(double lat, double lon, const CachedDataset& cache);
        double interpolateElevation(double lat, double lon, const CachedDataset& cache);
        double bilinearInteropolation(double lat, double lon, const CachedDataset& cache);

        // Smart data fetching
        std::vector<TerrainPoint> getTerrainDataSmart(const BoundingBox& bbox,
                                                      int grid_points_lat,
                                                      int grid_points_lon);

        // API query methods
        std::string queryAPISingle(const std::string& url);
        std::vector<TerrainPoint> queryAPI(const std::vector<std::pair<double, double>>& locations);
        std::vector<TerrainPoint> queryAPIBatched(const std::vector<std::pair<double, double>>& locations, 
                                                  int batch_size = 100);
};