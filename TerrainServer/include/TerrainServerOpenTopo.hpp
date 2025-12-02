/** 
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file TerrainServerOpenTopo.hpp
 * @brief TerrainServerOpenTopo class declaration
*/

#pragma once

#include <string>
#include <vector>

#include "ITerrainServer.hpp"

/**
 * @brief CachedDataset struct
 * @details This struct represents a cached dataset of terrain data.
 */
struct CachedDataset
{
    BoundingBox bbox;
    int grid_lat;
    int grid_lon;
    std::vector<TerrainPoint> points;
    std::string filename;
};

/**
 * @brief TerrainServerOpenTopo class
 * @details This class represents a terrain data server that fetches data from the OpenTopoData API.
 */
class TerrainServerOpenTopo : public ITerrainServer
{
    public:
        TerrainServerOpenTopo(const std::string& cache_directory = "./terrain_cache",
                              const std::string& api_base_url = "https://api.opentopodata.org/v1",
                              const std::string& dataset_name = "srtm90m");

        ~TerrainServerOpenTopo() override;

        // ITerrainServer interface - these are the ONLY public methods for getting elevation data
        double getElevation(double lat, double lon) override;
        std::vector<TerrainPoint> getElevationGrid(double latStart, double lonStart, 
                                                    double latEnd, double lonEnd, 
                                                    int numLatSamples, int numLonSamples) override;
        std::vector<TerrainPoint> getTerrainData(BoundingBox bbox, int numLatSamples, int numLonSamples);

        // Configuration
        void setBatchSize(int batchSize) { batch_size_ = batchSize; }

    private:
        std::string cache_dir_;
        std::string api_url_;
        std::string dataset_;
        int batch_size_;

        // Cache management
        std::vector<CachedDataset> loaded_caches_;

        // ===================================================================
        // Core fetching method - everything goes through this
        // ===================================================================
        std::vector<TerrainPoint> fetchElevations(const std::vector<std::pair<double, double>>& locations,
                                                   bool save_to_cache = false,
                                                   const BoundingBox* bbox = nullptr,
                                                   int grid_lat = 0, int grid_lon = 0);

        // ===================================================================
        // Cache operations
        // ===================================================================
        std::string getCacheFilename(const BoundingBox& bbox, int grid_points_lat, int grid_points_lon);
        bool loadFromCache(const std::string& filename, std::vector<TerrainPoint>& points);
        void saveToCache(const std::string& filename, const std::vector<TerrainPoint>& points);
        
        void loadAllCaches();
        std::vector<std::string> findAllCacheFiles();
        BoundingBox parseBBoxFromFilename(const std::string& filename);

        // ===================================================================
        // Cache interpolation
        // ===================================================================
        bool canInterpolatePoint(double lat, double lon, const CachedDataset& cache);
        double interpolateElevation(double lat, double lon, const CachedDataset& cache);
        double bilinearInterpolation(double lat, double lon, const CachedDataset& cache);

        // ===================================================================
        // API communication
        // ===================================================================
        static size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp);
        std::vector<TerrainPoint> queryAPI(const std::vector<std::pair<double, double>>& locations);
};