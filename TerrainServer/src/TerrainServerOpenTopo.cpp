#include <iostream>
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <netinet/in.h>
#include <unistd.h>
#include <thread>

#include "TerrainServerOpenTopo.hpp"

using json = nlohmann::json;

static size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* output)
{
    size_t totalSize = size * nmemb;
    output->append((char*)contents, totalSize);
    return totalSize;
}

TerrainServerOpenTopo::TerrainServerOpenTopo(const std::string& dataset)
    : dataset_(dataset)
{

}

std::string TerrainServerOpenTopo::queryAPI(const std::string& url)
{
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
            std::cerr << "CURL request failed: " << curl_easy_strerror(res) << "\n";
        }
        curl_easy_cleanup(curl);
    }
    return readBuffer;
}

double TerrainServerOpenTopo::getElevation(double lat, double lon)
{
    std::lock_guard<std::mutex> lock(ioMutex_);
    std::string url = "https://api.opentopodata.org/v1/" + dataset_ +
                      "?locations=" + std::to_string(lat) + "," + std::to_string(lon);

    std::string response = queryAPI(url);

    try 
    {
        auto j = json::parse(response);
        return j["results"][0]["elevation"].get<double>();
    }
    catch (...)
    {
        std::cerr << "Failed to parse JSON response or extract elevation.\n";
        return -9999.0; // Indicate error
    }
}

std::vector<TerrainPoint> TerrainServerOpenTopo::getElevationGrid(double latStart,
    double lonStart, double latEnd, double lonEnd, int numLatSamples, int numLonSamples)
{
    std::vector<TerrainPoint> grid;
    double dLat = (latEnd - latStart) / (numLatSamples - 1);
    double dLon = (lonEnd - lonStart) / (numLonSamples - 1);

    for (int i = 0; i < numLatSamples; ++i)
    {
        for (int j = 0; j < numLonSamples; ++j)
        {
            double lat = latStart + i * dLat;
            double lon = lonStart + j * dLon;
            double elev = getElevation(lat, lon);
            grid.push_back({lat, lon, elev});
        }
    }
    return grid;
}

void TerrainServerOpenTopo::startServer(int port)
{
    int server_fd = socket(AF_INET, SOCK_STREAM, 0);
    sockaddr_in address{};
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(port);

    bind(server_fd, (struct sockaddr*)&address, sizeof(address));
    listen(server_fd, 3);

    std::cout << "OpenTopo Terrain Server listening on port " << port << "...\n";

    while (true)
    {
        int clientSocket = accept(server_fd, nullptr, nullptr);
        std::thread([this, clientSocket]() {
            char buffer[1024] = {0};
            read(clientSocket, buffer, sizeof(buffer));
            double lat, lon;
            sscanf(buffer, "%lf,%lf", &lat, &lon);

            double elev = getElevation(lat, lon);
            std::string response = "Elevation: " + std::to_string(elev) + "\n";
            send(clientSocket, response.c_str(), response.size(), 0);

            close(clientSocket);
        }).detach();
    }
}