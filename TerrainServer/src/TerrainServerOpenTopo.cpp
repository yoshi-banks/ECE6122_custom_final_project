/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file TerrainServerOpenTopo.cpp
 * 
 * @brief TerrainServerOpenTopo class implementation
 */

#include "TerrainServerOpenTopo.hpp"
#include "Logger.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <iomanip>
#include <thread>
#include <chrono>

using json = nlohmann::json;
namespace fs = std::filesystem;

extern Logger logger;

// ============================================================================
// Constructor / Destructor
// ============================================================================

/**
 * @brief Constructor for TerrainServerOpenTopo class
 */
TerrainServerOpenTopo::TerrainServerOpenTopo(const std::string& cache_directory,
                                             const std::string& api_base_url,
                                             const std::string& dataset_name)
    : cache_dir_(cache_directory)
    , api_url_(api_base_url)
    , dataset_(dataset_name)
    , batch_size_(100)
{
    logger.log("[INFO] TerrainServerOpenTopo initialized with dataset: " + dataset_);
    
    if (!fs::exists(cache_dir_))
    {
        fs::create_directories(cache_dir_);
        logger.log("[INFO] Created cache directory: " + cache_dir_);
    }

    curl_global_init(CURL_GLOBAL_DEFAULT);
    
    loadAllCaches();
    logger.log("[INFO] Loaded " + std::to_string(loaded_caches_.size()) + " cached datasets");
}

/**
 * @brief Destructor for TerrainServerOpenTopo class
 */
TerrainServerOpenTopo::~TerrainServerOpenTopo()
{
    curl_global_cleanup();
}

// ============================================================================
// PUBLIC INTERFACE - Only two methods for getting elevation data
// ============================================================================

/**
 * @brief Get elevation at given latitude and longitude
 * @param lat Latitude in degrees
 * @param lon Longitude in degrees
 * @return Elevation in meters
 */
double TerrainServerOpenTopo::getElevation(double lat, double lon)
{
    std::vector<std::pair<double, double>> locations = {{lat, lon}};
    std::vector<TerrainPoint> results = fetchElevations(locations);
    
    if (results.empty())
    {
        logger.log("[ERROR] Failed to get elevation for (" + std::to_string(lat) + 
                   ", " + std::to_string(lon) + ")");
        return -9999.0;
    }
    
    return results[0].alt;
}

/**
 * @brief Get elevation grid for a given bounding box
 * @param latStart Starting latitude in degrees
 * @param lonStart Starting longitude in degrees
 * @param latEnd Ending latitude in degrees
 * @param lonEnd Ending longitude in degrees
 * @param numLatSamples Number of latitude samples
 * @param numLonSamples Number of longitude samples
 * @return Vector of TerrainPoints
 */
std::vector<TerrainPoint> TerrainServerOpenTopo::getElevationGrid(double latStart, double lonStart,
                                                                   double latEnd, double lonEnd,
                                                                   int numLatSamples, int numLonSamples)
{
    logger.log("[INFO] [TerrainServerOpenTopo::getElevationGrid] Getting elevation grid: [" + std::to_string(latStart) + ", " +
               std::to_string(lonStart) + "] to [" + std::to_string(latEnd) + ", " +
               std::to_string(lonEnd) + "] with " + std::to_string(numLatSamples) + "x" +
               std::to_string(numLonSamples) + " samples");

    // Check for exact cache match first
    BoundingBox bbox{latStart, latEnd, lonStart, lonEnd};
    std::string cache_file = getCacheFilename(bbox, numLatSamples, numLonSamples);
    
    std::vector<TerrainPoint> points;
    if (loadFromCache(cache_file, points))
    {
        logger.log("[INFO] [TerrainServerOpenTopo::getElevationGrid] Loaded exact grid match from cache");
        return points;
    }

    // Generate grid locations
    std::vector<std::pair<double, double>> locations;
    double lat_step = (latEnd - latStart) / (numLatSamples - 1);
    double lon_step = (lonEnd - lonStart) / (numLonSamples - 1);

    for (int i = 0; i < numLatSamples; ++i)
    {
        for (int j = 0; j < numLonSamples; ++j)
        {
            double lat = latStart + i * lat_step;
            double lon = lonStart + j * lon_step;
            locations.push_back({lat, lon});
        }
    }

    // Fetch elevations and save to cache
    points = fetchElevations(locations, true, &bbox, numLatSamples, numLonSamples);

    return points;
}

// ============================================================================
// CORE FETCHING METHOD - Everything goes through this
// ============================================================================

/**
 * @brief Fetch elevation data for a given set of locations
 * @param locations Vector of latitude-longitude pairs
 * @param save_to_cache Whether to save results to cache
 * @param bbox Optional bounding box for caching
 * @param grid_lat Number of latitude samples in grid (for caching)
 * @param grid_lon Number of longitude samples in grid (for caching)
 * @return Vector of TerrainPoints  
 */
std::vector<TerrainPoint> TerrainServerOpenTopo::fetchElevations(
    const std::vector<std::pair<double, double>>& locations,
    bool save_to_cache,
    const BoundingBox* bbox,
    int grid_lat,
    int grid_lon)
{
    std::vector<TerrainPoint> results;
    results.reserve(locations.size());
    
    std::vector<std::pair<size_t, std::pair<double, double>>> missing_locations;
    
    logger.log("[INFO] Fetching " + std::to_string(locations.size()) + " elevation points");
    
    // Step 1: Try to interpolate from existing caches
    for (size_t i = 0; i < locations.size(); ++i)
    {
        const auto& loc = locations[i];
        bool found = false;
        
        for (const auto& cache : loaded_caches_)
        {
            if (canInterpolatePoint(loc.first, loc.second, cache))
            {
                TerrainPoint pt;
                pt.lat = loc.first;
                pt.lon = loc.second;
                pt.alt = interpolateElevation(loc.first, loc.second, cache);
                results.push_back(pt);
                found = true;
                break;
            }
        }
        
        if (!found)
        {
            missing_locations.push_back({i, loc});
            // Add placeholder
            results.push_back(TerrainPoint{loc.first, loc.second, 0.0});
        }
    }
    
    int cached_count = locations.size() - missing_locations.size();
    if (cached_count > 0)
    {
        logger.log("[INFO] Found " + std::to_string(cached_count) + " points in cache");
    }
    
    // Step 2: Fetch missing points from API in batches
    if (!missing_locations.empty())
    {
        logger.log("[INFO] Fetching " + std::to_string(missing_locations.size()) + 
                   " points from API");
        
        size_t points_so_far = 0;
        for (size_t i = 0; i < missing_locations.size(); i += batch_size_)
        {
            size_t end = std::min(i + batch_size_, missing_locations.size());
            
            std::vector<std::pair<double, double>> batch_locations;
            std::vector<size_t> batch_indices;
            
            for (size_t j = i; j < end; ++j)
            {
                batch_locations.push_back(missing_locations[j].second);
                batch_indices.push_back(missing_locations[j].first);
            }
            
            points_so_far += batch_locations.size();
            std::cout << "API batch " << (i / batch_size_ + 1)
                      << " (" << batch_locations.size() << " points) "
                      << points_so_far << "/" << missing_locations.size() << " fetched" 
                      << std::endl;
            
            std::vector<TerrainPoint> batch_results = queryAPI(batch_locations);
            
            // Fill in results at correct indices
            for (size_t j = 0; j < batch_results.size() && j < batch_indices.size(); ++j)
            {
                results[batch_indices[j]] = batch_results[j];
            }
            
            // Rate limiting
            if (end < missing_locations.size())
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(500));
            }
        }
    }
    
    // Step 3: Save to cache if requested
    if (save_to_cache && !results.empty() && bbox != nullptr && grid_lat > 0 && grid_lon > 0)
    {
        std::string cache_file = getCacheFilename(*bbox, grid_lat, grid_lon);
        saveToCache(cache_file, results);
        
        // Add to loaded caches for future interpolation
        CachedDataset new_cache;
        new_cache.bbox = *bbox;
        new_cache.grid_lat = grid_lat;
        new_cache.grid_lon = grid_lon;
        new_cache.points = results;
        new_cache.filename = cache_file;
        loaded_caches_.push_back(new_cache);
        
        logger.log("[INFO] Saved " + std::to_string(results.size()) + " points to cache");
    }
    
    return results;
}

// ============================================================================
// CACHE OPERATIONS
// ============================================================================

/**
 * @brief Generate cache filename based on bounding box and grid size
 * @param bbox Bounding box of the dataset
 * @param grid_points_lat Number of latitude samples in grid
 * @param grid_points_lon Number of longitude samples in grid
 * @return Cache filename as string
 */
std::string TerrainServerOpenTopo::getCacheFilename(const BoundingBox& bbox,
                                                    int grid_points_lat,
                                                    int grid_points_lon)
{
    std::ostringstream oss;
    oss << cache_dir_ << "/terrain_"
        << std::fixed << std::setprecision(6)
        << bbox.min_lat << "_" << bbox.max_lat << "_"
        << bbox.min_lon << "_" << bbox.max_lon << "_"
        << grid_points_lat << "x" << grid_points_lon << ".json";
    return oss.str();
}

/**
 * @brief Load cached dataset from file
 * @param filename Cache filename
 * @param points Vector to populate with loaded TerrainPoints
 * @return true if loaded successfully, false otherwise
 */
bool TerrainServerOpenTopo::loadFromCache(const std::string& filename, 
                                          std::vector<TerrainPoint>& points)
{
    if (!fs::exists(filename))
    {
        return false;
    }

    try
    {
        std::ifstream file(filename);
        json j;
        file >> j;

        points.clear();
        for (const auto& item : j["points"])
        {
            TerrainPoint pt;
            pt.lat = item["lat"];
            pt.lon = item["lon"];
            pt.alt = item["alt"];
            points.push_back(pt);
        }

        logger.log("[INFO] [TerrainServerOpenTopo::loadFromCache] Loaded " + std::to_string(points.size()) + 
                   " points from cache: " + filename);
        return true;
    }
    catch (const std::exception& e)
    {
        logger.log("[ERROR] Error loading cache: " + std::string(e.what()));
        return false;
    }
}

/**
 * @brief Save dataset to cache file
 * @param filename Cache filename
 * @param points Vector of TerrainPoints to save
 */
void TerrainServerOpenTopo::saveToCache(const std::string& filename, 
                                        const std::vector<TerrainPoint>& points)
{
    try 
    {
        json j;
        j["points"] = json::array();

        for (const auto& pt : points)
        {
            j["points"].push_back({
                {"lat", pt.lat},
                {"lon", pt.lon},
                {"alt", pt.alt}
            });
        }

        std::ofstream file(filename);
        file << j.dump(2);

        logger.log("[INFO] [TerrainServerOpenTopo::saveToCache]Saved " + std::to_string(points.size()) + 
                   " points to cache: " + filename);
    }
    catch (const std::exception& e)
    {
        logger.log("[ERROR] Error saving cache: " + std::string(e.what()));
    }
}

/**
 * @brief Load all cached datasets from cache directory
 */
void TerrainServerOpenTopo::loadAllCaches()
{
    loaded_caches_.clear();
    std::vector<std::string> cache_files = findAllCacheFiles();

    for (const auto& filename : cache_files)
    {
        CachedDataset cache;
        if (loadFromCache(filename, cache.points))
        {
            cache.bbox = parseBBoxFromFilename(filename);
            cache.filename = filename;

            // Extract grid size from filename
            size_t grid_pos = filename.find_last_of("_");
            size_t ext_pos = filename.find_last_of(".");
            if (grid_pos != std::string::npos && ext_pos != std::string::npos)
            {
                std::string grid_str = filename.substr(grid_pos + 1, ext_pos - grid_pos - 1);
                size_t x_pos = grid_str.find('x');
                if (x_pos != std::string::npos)
                {
                    cache.grid_lat = std::stoi(grid_str.substr(0, x_pos));
                    cache.grid_lon = std::stoi(grid_str.substr(x_pos + 1));
                }
            }

            loaded_caches_.push_back(cache);
        }
    }
}

/**
 * @brief Find all cache files in the cache directory
 * @return Vector of cache filenames
 */
std::vector<std::string> TerrainServerOpenTopo::findAllCacheFiles()
{
    std::vector<std::string> cache_files;

    if (!fs::exists(cache_dir_))
    {
        return cache_files;
    }

    for (const auto& entry : fs::directory_iterator(cache_dir_))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".json")
        {
            cache_files.push_back(entry.path().string());
        }
    }

    return cache_files;
}

/**
 * @brief Parse bounding box from filename
 * @param filename Cache filename
 * @return Bounding box
 */
BoundingBox TerrainServerOpenTopo::parseBBoxFromFilename(const std::string& filename)
{
    BoundingBox bbox = {0, 0, 0, 0};
    
    size_t start = filename.find("terrain_");
    if (start == std::string::npos) return bbox;
    
    start += 8; // Skip "terrain_"
    
    std::vector<double> values;
    std::string num_str;
    bool negative = false;
    
    for (size_t i = start; i < filename.size(); ++i)
    {
        char c = filename[i];
        
        if (c == '-')
        {
            if (!num_str.empty())
            {
                values.push_back(negative ? -std::stod(num_str) : std::stod(num_str));
                num_str.clear();
            }
            negative = true;
        }
        else if (c == '_' || c == 'x' || c == '.')
        {
            if (!num_str.empty())
            {
                values.push_back(negative ? -std::stod(num_str) : std::stod(num_str));
                num_str.clear();
                negative = false;
            }
            if (c == 'x' || c == '.') break;
        }
        else if (std::isdigit(c) || c == '.')
        {
            num_str += c;
        }
    }
    
    if (values.size() >= 4)
    {
        bbox.min_lat = values[0];
        bbox.max_lat = values[1];
        bbox.min_lon = values[2];
        bbox.max_lon = values[3];
    }
    
    return bbox;
}

// ============================================================================
// CACHE INTERPOLATION
// ============================================================================

/**
 * @brief Check if point can be interpolated from cache
 * @param lat Latitude of the point
 * @param lon Longitude of the point
 * @param cache Cached dataset
 * @return True if point can be interpolated
 */
bool TerrainServerOpenTopo::canInterpolatePoint(double lat, double lon, 
                                                 const CachedDataset& cache)
{
    return (lat >= cache.bbox.min_lat - interpolationTolerance_ &&
            lat <= cache.bbox.max_lat + interpolationTolerance_ &&
            lon >= cache.bbox.min_lon - interpolationTolerance_ &&
            lon <= cache.bbox.max_lon + interpolationTolerance_);
}

/**
 * @brief Interpolate elevation from cache
 * @param lat Latitude of the point
 * @param lon Longitude of the point
 * @param cache Cached dataset
 * @return Elevation in meters
 */
double TerrainServerOpenTopo::interpolateElevation(double lat, double lon, 
                                                    const CachedDataset& cache)
{
    return bilinearInterpolation(lat, lon, cache);
}

/**
 * @brief Perform bilinear interpolation to get elevation at given lat/lon
 * @param lat Latitude in degrees
 * @param lon Longitude in degrees
 * @param cache Cached dataset
 * @return Elevation in meters
 */
double TerrainServerOpenTopo::bilinearInterpolation(double lat, double lon, 
                                                     const CachedDataset& cache)
{
    double lat_step = (cache.bbox.max_lat - cache.bbox.min_lat) / (cache.grid_lat - 1);
    double lon_step = (cache.bbox.max_lon - cache.bbox.min_lon) / (cache.grid_lon - 1);
    
    double lat_idx = (lat - cache.bbox.min_lat) / lat_step;
    double lon_idx = (lon - cache.bbox.min_lon) / lon_step;
    
    int lat_i0 = static_cast<int>(std::floor(lat_idx));
    int lat_i1 = static_cast<int>(std::ceil(lat_idx));
    int lon_j0 = static_cast<int>(std::floor(lon_idx));
    int lon_j1 = static_cast<int>(std::ceil(lon_idx));
    
    // Clamp to bounds
    lat_i0 = std::max(0, std::min(lat_i0, cache.grid_lat - 1));
    lat_i1 = std::max(0, std::min(lat_i1, cache.grid_lat - 1));
    lon_j0 = std::max(0, std::min(lon_j0, cache.grid_lon - 1));
    lon_j1 = std::max(0, std::min(lon_j1, cache.grid_lon - 1));
    
    auto getElevationAt = [&](int i, int j) -> double {
        int idx = i * cache.grid_lon + j;
        if (idx >= 0 && idx < static_cast<int>(cache.points.size()))
        {
            return cache.points[idx].alt;
        }
        return 0.0;
    };
    
    double z00 = getElevationAt(lat_i0, lon_j0);
    double z01 = getElevationAt(lat_i0, lon_j1);
    double z10 = getElevationAt(lat_i1, lon_j0);
    double z11 = getElevationAt(lat_i1, lon_j1);
    
    if (lat_i0 == lat_i1 && lon_j0 == lon_j1)
    {
        return z00;
    }
    
    double lat_weight = lat_idx - lat_i0;
    double lon_weight = lon_idx - lon_j0;
    
    double z0 = z00 * (1 - lon_weight) + z01 * lon_weight;
    double z1 = z10 * (1 - lon_weight) + z11 * lon_weight;
    double z = z0 * (1 - lat_weight) + z1 * lat_weight;
    
    return z;
}

// ============================================================================
// API COMMUNICATION
// ============================================================================

/**
 * @brief Callback function for libcurl to write data to a string
 * @param contents Pointer to the string to write to
 * @param size Size of each element
 * @param nmemb Number of elements
 * @param userp Pointer to the string
 * @return Number of elements written
 */
size_t TerrainServerOpenTopo::WriteCallback(void* contents, size_t size, 
                                             size_t nmemb, void* userp)
{
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

/**
 * @brief Query OpenTopoData API for elevation data
 * @param locations Vector of latitude-longitude pairs
 * @return Vector of TerrainPoints
 */
std::vector<TerrainPoint> TerrainServerOpenTopo::queryAPI(
    const std::vector<std::pair<double, double>>& locations)
{
    std::vector<TerrainPoint> points;

    // Build locations string
    std::ostringstream loc_stream;
    for (size_t i = 0; i < locations.size(); ++i)
    {
        if (i > 0) loc_stream << "|";
        loc_stream << std::fixed << std::setprecision(6)
                   << locations[i].first << "," << locations[i].second;
    }

    // Build URL
    std::string url = api_url_ + "/" + dataset_ + "?locations=" + loc_stream.str();

    // Make HTTP request
    CURL* curl = curl_easy_init();
    std::string response;

    if (curl)
    {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);

        CURLcode res = curl_easy_perform(curl);

        if (res != CURLE_OK)
        {
            logger.log("[ERROR] CURL error: " + std::string(curl_easy_strerror(res)));
        }
        else
        {
            try 
            {
                json j = json::parse(response);

                if (j["status"] == "OK")
                {
                    for (const auto& result : j["results"])
                    {
                        TerrainPoint pt;
                        pt.lat = result["location"]["lat"];
                        pt.lon = result["location"]["lng"];
                        pt.alt = result["elevation"].is_null() ? 0.0 : result["elevation"].get<double>();
                        points.push_back(pt);
                    }
                }
                else
                {
                    logger.log("[ERROR] API error: " + j.value("error", "Unknown error"));
                }
            }
            catch (const std::exception& e)
            {
                logger.log("[ERROR] JSON parsing error: " + std::string(e.what()));
            }
        }

        curl_easy_cleanup(curl);
    }

    return points;
}

/**
 * @brief Get terrain data for a given bounding box
 * @param bbox Bounding box
 * @param numLatSamples Number of latitude samples
 * @param numLonSamples Number of longitude samples
 * @return Vector of TerrainPoints
 */
std::vector<TerrainPoint> TerrainServerOpenTopo::getTerrainData(BoundingBox bbox, int numLatSamples, int numLonSamples) 
{ 
    return getElevationGrid(bbox.min_lat, bbox.min_lon, bbox.max_lat, bbox.max_lon, numLatSamples, numLonSamples); 
}
