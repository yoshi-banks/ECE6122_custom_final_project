#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <iomanip>
#include <netinet/in.h>
#include <unistd.h>
#include <thread>
#include <chrono>

#include "Logger.hpp"
#include "TerrainServerOpenTopo.hpp"

using json = nlohmann::json;
namespace fs = std::filesystem;

Logger logger("TerrainServerOpenTopo.log");

TerrainServerOpenTopo::TerrainServerOpenTopo(const std::string& cache_directory,
                                             const std::string& api_base_url,
                                             const std::string& dataset_name)
    : cache_dir(cache_directory)
    , api_url(api_base_url)
    , dataset(dataset_name)
    , interpolationTolerance_(0.01) // Default tolerance ~1km
    
{
    logger.log("[INFO] TerrainServerOpenTopo initialized with dataset: " + dataset);
    if (!fs::exists(cache_dir))
    {
        fs::create_directories(cache_dir);
        logger.log("[INFO] Created cache directory: " + cache_dir);
    }

    // Initialize CURL
    curl_global_init(CURL_GLOBAL_DEFAULT);

    // Load all existing caches
    loadAllCaches();
    logger.log("[INFO] Loaded " + std::to_string(loadedCaches_.size()) + " cached datasets.");
}

TerrainServerOpenTopo::~TerrainServerOpenTopo()
{
    // Cleanup CURL
    curl_global_cleanup();
}

size_t TerrainServerOpenTopo::WriteCallback(void* contents, size_t size, size_t nmemb, void* userp)
{
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

std::string TerrainServerOpenTopo::getCacheFilename(const BoundingBox& bbox,
                                                    int grid_points_lat,
                                                    int grid_points_lon)
{
    std::ostringstream oss;
    oss << cache_dir << "/terrain_"
        << std::fixed << std::setprecision(6)
        << bbox.min_lat << "_" << bbox.max_lat << "_"
        << bbox.min_lon << "_" << bbox.max_lon << "_"
        << grid_points_lat << "x" << grid_points_lon << ".json";
    return oss.str();
}

bool TerrainServerOpenTopo::loadFromCache(const std::string& filename, std::vector<TerrainPoint>& points)
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

        std::cout << "Loaded " << points.size() << " points from cache: " << filename << std::endl;
        logger.log("[INFO] Loaded " + std::to_string(points.size()) + " points from cache: " + filename);
        return true;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading cache: " << e.what() << std::endl;
        logger.log("[ERROR] Error loading cache: " + std::string(e.what()));
        return false;
    }
}

void TerrainServerOpenTopo::saveToCache(const std::string& filename, const std::vector<TerrainPoint>& points)
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

        std::cout << "Saved " << points.size() << " points to cache: " << filename << std::endl;
        logger.log("[INFO] Saved " + std::to_string(points.size()) + " points to cache: " + filename);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error saving cache: " << e.what() << std::endl;
        logger.log("[ERROR] Error saving cache: " + std::string(e.what()));
    }
}

std::vector<TerrainPoint> TerrainServerOpenTopo::queryAPI(const std::vector<std::pair<double, double>>& locations)
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
    std::string url = api_url + "/" + dataset + "?locations=" + loc_stream.str();

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
            std::cerr << "CURL error: " << curl_easy_strerror(res) << std::endl;
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

                        if (result["elevation"].is_null())
                        {
                            pt.alt = 0.0; // or use nan
                        }
                        else
                        {
                            pt.alt = result["elevation"];
                        }

                        points.push_back(pt);
                    }
                }
                else
                {
                    std::cerr << "API error: " << j.value("error", "Unknown error") << std::endl;
                    logger.log("[ERROR] API error: " + j.value("error", "Unknown error"));
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "JSON parsing error: " << e.what() << std::endl;
                logger.log("[ERROR] JSON parsing error: " + std::string(e.what()));
            }
        }

        curl_easy_cleanup(curl);
    }

    return points;
}

std::vector<TerrainPoint> TerrainServerOpenTopo::queryAPIBatched(
    const std::vector<std::pair<double, double>>& locations,
    int batch_size)
{
    std::vector<TerrainPoint> all_points;

    for (size_t i = 0; i < locations.size(); i += batch_size)
    {
        size_t end = std::min(i + batch_size, locations.size());
        std::vector<std::pair<double, double>> batch(locations.begin() + i, locations.begin() + end);

        std::cout << "Querying batch " << (i / batch_size + 1)
                  << " (" << batch.size() << " points)..." << std::endl;

        std::vector<TerrainPoint> batch_points = queryAPI(batch);
        all_points.insert(all_points.end(), batch_points.begin(), batch_points.end());

        // Be nice to the API - add a small delay between batches
        if (end < locations.size())
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
    }

    return all_points;
}

std::vector<TerrainPoint> TerrainServerOpenTopo::getTerrainData(const BoundingBox& bbox, 
                                                                int grid_points_lat, 
                                                                int grid_points_lon)
{
    return getTerrainDataSmart(bbox, grid_points_lat, grid_points_lon);
}

std::vector<TerrainPoint> TerrainServerOpenTopo::pygetTerrainDataSmart(const BoundingBox& bbox,
                                                                     int grid_points_lat,
                                                                     int grid_points_lon)
{
    std::string cache_file = getCacheFilename(bbox, grid_points_lat, grid_points_lon);

    // Try to load exact match from cache
    std::vector<TerrainPoint> points;
    if (loadFromCache(cache_file, points))
    {
        std::cout << "Loaded exact match from cache" << std::endl;
        return points;
    }

    // Generate requested grid points
    std::vector<std::pair<double, double>> requested_locations;
    double lat_step = (bbox.max_lat - bbox.min_lat) / (grid_points_lat - 1);
    double lon_step = (bbox.max_lon - bbox.min_lon) / (grid_points_lon - 1);

    for (int i = 0; i < grid_points_lat; ++i)
    {
        for (int j = 0; j < grid_points_lon; ++j)
        {
            double lat = bbox.min_lat + i * lat_step;
            double lon = bbox.min_lon + j * lon_step;
            requested_locations.push_back({lat, lon});
        }
    }

    // Try to interpolate from existing caches
    std::vector<TerrainPoint> interpolated_points;
    std::vector<std::pair<double, double>> missing_locations;
    int interpolated_count = 0;

    for (const auto& loc : requested_locations)
    {
        bool found = false;

        // Check all loaded caches
        for (const auto& cache : loadedCaches_)
        {
            if (canInterpolatePoint(loc.first, loc.second, cache))
            {
                TerrainPoint pt;
                pt.lat = loc.first;
                pt.lon = loc.second;
                pt.alt = interpolateElevation(loc.first, loc.second, cache);
                interpolated_points.push_back(pt);
                interpolated_count++;
                found = true;
                break;
            }
        }

        if (!found)
        {
            missing_locations.push_back(loc);
        }
    }

    std::cout << "Interpolated " << interpolated_count << " points from "
              << loadedCaches_.size() << " cached datasets\n";
    logger.log("[INFO] Interpolated " + std::to_string(interpolated_count) + " points from cache");

    // Fetch missing data from API
    if (!missing_locations.empty())
    {
        std::cout << "Fetching " << missing_locations.size() << " missing points from API...\n";
        logger.log("[INFO] Fetching " + std::to_string(missing_locations.size()) + " missing points from API");

        std::vector<TerrainPoint> fetched_points = queryAPIBatched(missing_locations);
        interpolated_points.insert(interpolated_points.end(),
                                   fetched_points.begin(),
                                   fetched_points.end());
    }

    // Sort points to match original grid order
    std::sort(interpolated_points.begin(), interpolated_points.end(),
        [](const TerrainPoint& a, const TerrainPoint& b) {
            if (std::abs(a.lat - b.lat) < 1e-6)
                return a.lon < b.lon;
            return a.lat < b.lat;
        });

    // Save to cache for future use
    if (!interpolated_points.empty())
    {
        saveToCache(cache_file, interpolated_points);

        // Add to loaded caches
        CachedDataset new_cache;
        new_cache.bbox = bbox;
        new_cache.grid_lat = grid_points_lat;
        new_cache.grid_lon = grid_points_lon;
        new_cache.points = interpolated_points;
        new_cache.filename = cache_file;
        loadedCaches_.push_back(new_cache);
    }

    return interpolated_points;
}

void TerrainServerOpenTopo::exportToCSV(const std::vector<TerrainPoint>& points, const std::string& filename)
{
    std::ofstream file(filename);
    file << "latitude,longitude,elevation\n";

    for (const auto& pt: points)
    {
        file << std::fixed << std::setprecision(6)
             << pt.lat << "," << pt.lon << "," << pt.alt << "\n";
    }

    std::cout << "Exported to CSV: " << filename << std::endl;
    logger.log("[INFO] Exported to CSV: " + filename);
}

std::string TerrainServerOpenTopo::queryAPISingle(const std::string& url)
{
    logger.log("[INFO] Querying OpenTopo API: " + url);

    CURL* curl = curl_easy_init();
    std::string readBuffer;
    if (curl)
    {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
        CURLcode res = curl_easy_perform(curl);
        if (res != CURLE_OK)
        {
            logger.log("[ERROR] CURL request failed: " + std::string(curl_easy_strerror(res)));
            std::cerr << "CURL request failed: " << curl_easy_strerror(res) << "\n";
        }
        else
        {
            logger.log("[INFO] CURL request successful.");
        }
        curl_easy_cleanup(curl);
    }
    else
    {
        logger.log("[ERROR] Failed to initialize CURL.");
        std::cerr << "Failed to initialize CURL.\n";
    }
    return readBuffer;
}

double TerrainServerOpenTopo::getElevation(double lat, double lon)
{
    std::lock_guard<std::mutex> lock(ioMutex_);
    std::string url = api_url + "/" + dataset +
                      "?locations=" + std::to_string(lat) + "," + std::to_string(lon);

    logger.log("[INFO] Getting elevation for (" + std::to_string(lat) + ", " + std::to_string(lon) + ")");
    std::string response = queryAPISingle(url);

    try 
    {
        auto j = json::parse(response);
        auto& elevationField = j["results"][0]["elevation"];
        if (elevationField.is_null())
        {
            logger.log("[WARN] Elevation data not available for the given coordinates.");
            return -9999.0; // Indicate no data
        }
        double elevation = elevationField.get<double>();
        logger.log("[INFO] Retrieved elevation: " + std::to_string(elevation) + " m");
        return elevation;
    }
    catch (const std::exception& e)
    {
        logger.log("[ERROR] JSON parsing failed: " + std::string(e.what()));
        std::cerr << "JSON parsing failed: " << e.what() << "\n";
        return -9999.0; // Indicate error
    }
    catch (...)
    {
        logger.log("[ERROR] Unknown error during JSON parsing.");
        std::cerr << "Unknown error during JSON parsing.\n";
        return -9999.0; // Indicate error
    }
}

std::vector<TerrainPoint> TerrainServerOpenTopo::getElevationGrid(double latStart,
    double lonStart, double latEnd, double lonEnd, int numLatSamples, int numLonSamples)
{
    logger.log("[INFO] Getting elevation grid from (" + std::to_string(latStart) + ", " +
                std::to_string(lonStart) + ") to (" + std::to_string(latEnd) + ", " +
                std::to_string(lonEnd) + ") with " + std::to_string(numLatSamples) +
                "x" + std::to_string(numLonSamples) + " samples.");

    // Use the new batched querying method
    BoundingBox bbox;
    bbox.min_lat = latStart;
    bbox.max_lat = latEnd;
    bbox.min_lon = lonStart;
    bbox.max_lon = lonEnd;

    return getTerrainData(bbox, numLatSamples, numLonSamples);
}

void TerrainServerOpenTopo::startServer(int port)
{
    int server_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (server_fd == 0)
    {
        logger.log("[ERROR] Socket creation failed.");
        std::cerr << "Socket creation failed.\n";
        return;
    }

    sockaddr_in address{};
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(port);

    if (bind(server_fd, (struct sockaddr*)&address, sizeof(address)) < 0)
    {
        logger.log("[ERROR] Socket bind failed.");
        std::cerr << "Socket bind failed.\n";
        return;
    }

    if (listen(server_fd, 3) < 0)
    {
        logger.log("[ERROR] Socket listen failed.");
        std::cerr << "Socket listen failed.\n";
        return;
    }

    logger.log("[INFO] TerrainServerOpenTopo listening on port " + std::to_string(port));
    std::cout << "OpenTopo Terrain Server listening on port " << port << "...\n";

    while (true)
    {
        int clientSocket = accept(server_fd, nullptr, nullptr);
        if (clientSocket < 0)
        {
            logger.log("[ERROR] Socket accept failed.");
            std::cerr << "Socket accept failed.\n";
            continue;
        }

        std::thread([this, clientSocket]() {
            char buffer[1024] = {0};
            read(clientSocket, buffer, sizeof(buffer));
            double lat, lon;
            if (sscanf(buffer, "%lf,%lf", &lat, &lon) != 2)
            {
                std::string errorMsg = "Invalid request format. Use: <lat>,<lon>\n";
                send(clientSocket, errorMsg.c_str(), errorMsg.size(), 0);
                close(clientSocket);
                return;
            }
            
            logger.log("[INFO] Received request for elevation at (" + std::to_string(lat) + ", " + std::to_string(lon) + ")");

            double elev = getElevation(lat, lon);
            std::string response = "Elevation: " + std::to_string(elev) + "\n";
            send(clientSocket, response.c_str(), response.size(), 0);

            logger.log("[INFO] Sent elevation response: " + std::to_string(elev));
            close(clientSocket);
        }).detach();
    }
}

// ============================================================================
// Cache Management Functions
// ============================================================================

std::vector<std::string> TerrainServerOpenTopo::findAllCacheFiles()
{
    std::vector<std::string> cache_files;

    if (!fs::exists(cache_dir))
    {
        return cache_files;
    }

    for (const auto& entry : fs::directory_iterator(cache_dir))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".json")
        {
            cache_files.push_back(entry.path().string());
        }
    }

    return cache_files;
}

void TerrainServerOpenTopo::loadAllCaches()
{
    loadedCaches_.clear();
    std::vector<std::string> cache_files = findAllCacheFiles();

    for (const auto& filename : cache_files)
    {
        CachedDataset cache;
        if (loadFromCache(filename, cache.points))
        {
            cache.bbox = parseBBoxFromFilename(filename);
            cache.filename = filename;

            // Extract grid size from filename (format: terrain_lat1_lat2_lon1_lon2_NxM.json)
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

            loadedCaches_.push_back(cache);
        }
    }
}

BoundingBox TerrainServerOpenTopo::parseBBoxFromFilename(const std::string& filename)
{
    BoundingBox bbox = {0, 0, 0, 0};
    
    // Extract bbox from filename format: terrain_lat1_lat2_lon1_lon2_NxM.json
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
            if (c == 'x' || c == '.') break; // Stop at grid size
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
// Interpolation Functions
// ============================================================================

bool TerrainServerOpenTopo::canInterpolatePoint(double lat, double lon, const CachedDataset& cache)
{
    // Check if point is within bbox with tolerance
    if (lat < cache.bbox.min_lat - interpolationTolerance_ ||
        lat > cache.bbox.max_lat + interpolationTolerance_ ||
        lon < cache.bbox.min_lon - interpolationTolerance_ ||
        lon > cache.bbox.max_lon + interpolationTolerance_)
    {
        return false;
    }

    return true;
}

double TerrainServerOpenTopo::interpolateElevation(double lat, double lon, const CachedDataset& cache)
{
    return bilinearInteropolation(lat, lon, cache);
}

double TerrainServerOpenTopo::bilinearInteropolation(double lat, double lon, const CachedDataset& cache)
{
    // Find the four nearest grid points
    double lat_step = (cache.bbox.max_lat - cache.bbox.min_lat) / (cache.grid_lat - 1);
    double lon_step = (cache.bbox.max_lon - cache.bbox.min_lon) / (cache.grid_lon - 1);
    
    // Find grid indices
    double lat_idx = (lat - cache.bbox.min_lat) / lat_step;
    double lon_idx = (lon - cache.bbox.min_lon) / lon_step;
    
    int lat_i0 = static_cast<int>(std::floor(lat_idx));
    int lat_i1 = static_cast<int>(std::ceil(lat_idx));
    int lon_j0 = static_cast<int>(std::floor(lon_idx));
    int lon_j1 = static_cast<int>(std::ceil(lon_idx));
    
    // Clamp to grid bounds
    lat_i0 = std::max(0, std::min(lat_i0, cache.grid_lat - 1));
    lat_i1 = std::max(0, std::min(lat_i1, cache.grid_lat - 1));
    lon_j0 = std::max(0, std::min(lon_j0, cache.grid_lon - 1));
    lon_j1 = std::max(0, std::min(lon_j1, cache.grid_lon - 1));
    
    // Get elevations at four corners
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
    
    // If all corners are the same (at edge), return that value
    if (lat_i0 == lat_i1 && lon_j0 == lon_j1)
    {
        return z00;
    }
    
    // Compute interpolation weights
    double lat_weight = lat_idx - lat_i0;
    double lon_weight = lon_idx - lon_j0;
    
    // Bilinear interpolation
    double z0 = z00 * (1 - lon_weight) + z01 * lon_weight;
    double z1 = z10 * (1 - lon_weight) + z11 * lon_weight;
    double z = z0 * (1 - lat_weight) + z1 * lat_weight;
    
    return z;
}