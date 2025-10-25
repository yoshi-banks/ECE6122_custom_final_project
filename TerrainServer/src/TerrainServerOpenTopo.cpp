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

#include "Logger.hpp"
#include "TerrainServerOpenTopo.hpp"

using json = nlohmann::json;
namespace fs = std::filesystem;

Logger logger("TerrainServerOpenTopo.log");

TerrainServerOpenTopo::TerrainServerOpenTopo(const std::string& cache_directory,
                                             const std::string& api_base_url,
                                             const std::string& dataset_name)
    : cache_dir(cache_directory), api_url(api_base_url), dataset(dataset_name)
{
    logger.log("[INFO] TerrainServerOpenTopo initialized with dataset: " + dataset);
    if (!fs::exists(cache_dir))
    {
        fs::create_directories(cache_dir);
        logger.log("[INFO] Created cache directory: " + cache_dir);
    }

    // Initialize CURL
    curl_global_init(CURL_GLOBAL_DEFAULT);
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
    std::string cache_file = getCacheFilename(bbox, grid_points_lat, grid_points_lon);

    // Try to load from cache
    std::vector<TerrainPoint> points;
    if (loadFromCache(cache_file, points))
    {
        return points;
    }

    // Generate grid of lat/lon points
    std::vector<std::pair<double, double>> locations;
    double lat_step = (bbox.max_lat - bbox.min_lat) / (grid_points_lat - 1);
    double lon_step = (bbox.max_lon - bbox.min_lon) / (grid_points_lon - 1);

    for (int i = 0; i < grid_points_lat; ++i)
    {
        for (int j = 0; j < grid_points_lon; ++j)
        {
            double lat = bbox.min_lat + i * lat_step;
            double lon = bbox.min_lon + j * lon_step;
            locations.push_back({lat, lon});
        }
    }

    std::cout << "Querying " << locations.size() << " points from API..." << std::endl;
    logger.log("[INFO] Querying " + std::to_string(locations.size()) + " points from API...");

    // Query API in batches
    points = queryAPIBatched(locations);

    // Save to cache
    if (!points.empty())
    {
        saveToCache(cache_file, points);
    }

    return points;
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