/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file ITerrainServer.cpp
 * @brief ITerrainServer class implementation
 */

#include "ITerrainServer.hpp"
#include "Logger.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <netinet/in.h>
#include <unistd.h>
#include <thread>
#include <cstring>

extern Logger logger;

/**
 * @brief Starts the terrain server to listen for client requests
 * @param port Port number to listen on
 */
void ITerrainServer::startServer(int port)
{
    int server_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (server_fd == 0)
    {
        logger.log("[ERROR] Socket creation failed.");
        std::cerr << "Socket creation failed.\n";
        return;
    }

    // Set socket options to reuse address
    int opt = 1;
    if (setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt)))
    {
        logger.log("[WARN] Failed to set SO_REUSEADDR");
    }

    sockaddr_in address{};
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(port);

    if (bind(server_fd, (struct sockaddr*)&address, sizeof(address)) < 0)
    {
        logger.log("[ERROR] Socket bind failed.");
        std::cerr << "Socket bind failed.\n";
        close(server_fd);
        return;
    }

    if (listen(server_fd, 3) < 0)
    {
        logger.log("[ERROR] Socket listen failed.");
        std::cerr << "Socket listen failed.\n";
        close(server_fd);
        return;
    }

    logger.log("[INFO] TerrainServer listening on port " + std::to_string(port));
    std::cout << "Terrain Server listening on port " << port << "...\n";

    while (true)
    {
        int clientSocket = accept(server_fd, nullptr, nullptr);
        if (clientSocket < 0)
        {
            logger.log("[ERROR] Socket accept failed.");
            std::cerr << "Socket accept failed.\n";
            continue;
        }

        // Handle client in a separate thread
        std::thread([this, clientSocket]() {
            handleClient(clientSocket);
        }).detach();
    }

    close(server_fd);
}

/**
 * @brief Handles a client connection
 * @param clientSocket Socket descriptor for the client
 */
void ITerrainServer::handleClient(int clientSocket)
{
    char buffer[1024] = {0};
    ssize_t bytes_read = read(clientSocket, buffer, sizeof(buffer) - 1);
    
    if (bytes_read <= 0)
    {
        close(clientSocket);
        return;
    }

    double lat, lon;
    if (sscanf(buffer, "%lf,%lf", &lat, &lon) != 2)
    {
        std::string errorMsg = "Invalid request format. Use: <lat>,<lon>\n";
        send(clientSocket, errorMsg.c_str(), errorMsg.size(), 0);
        close(clientSocket);
        return;
    }
    
    logger.log("[INFO] Received request for elevation at (" + 
              std::to_string(lat) + ", " + std::to_string(lon) + ")");

    double elev = getElevation(lat, lon);
    std::string response = "Elevation: " + std::to_string(elev) + "\n";
    send(clientSocket, response.c_str(), response.size(), 0);

    logger.log("[INFO] Sent elevation response: " + std::to_string(elev));
    close(clientSocket);
}

/**
 * @brief Exports terrain points to a CSV file
 * @param points Vector of TerrainPoint to export
 * @param filename Output CSV filename
 */
void ITerrainServer::exportToCSV(const std::vector<TerrainPoint>& points, const std::string& filename)
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

/**
 * @brief Checks if a point is within the bounding box
 * @param lat Latitude of the point
 * @param lon Longitude of the point
 */
bool ITerrainServer::isPointInBounds(double lat, double lon) const
{
    return (lat >= bbox_.min_lat && lat <= bbox_.max_lat &&
            lon >= bbox_.min_lon && lon <= bbox_.max_lon);
}

/**
 * @brief Checks if a point is interpolatable based on tolerance
 * @param lat Latitude of the point
 * @param lon Longitude of the point
 */
bool ITerrainServer::isPointInterpolatable(double lat, double lon) const
{
    // Check if point is within bounds plus tolerance
    return (lat >= bbox_.min_lat - interpolationTolerance_ && 
            lat <= bbox_.max_lat + interpolationTolerance_ &&
            lon >= bbox_.min_lon - interpolationTolerance_ && 
            lon <= bbox_.max_lon + interpolationTolerance_);
}