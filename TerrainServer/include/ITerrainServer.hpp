/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file ITerrainServer.hpp
 * @brief ITerrainServer class declaration
 */

#pragma once

#include <string>
#include <vector>
#include <mutex>

/**
 * @brief TerrainPoint struct
 * @details This struct represents a point on the terrain with latitude, longitude, and altitude information.
 */
struct TerrainPoint 
{
    double lat;
    double lon;
    double alt;
};

/**
 * @brief BoundingBox struct
 * @details This struct represents a bounding box defined by minimum and maximum latitude and longitude.
 */
struct BoundingBox 
{
    double min_lat;
    double max_lat;
    double min_lon;
    double max_lon;
};

/**
 * @brief ITerrainServer interface
 * 
 * @details This interface defines the methods for a terrain data server.
 */
class ITerrainServer
{
    public:
        virtual ~ITerrainServer() = default;

        // Pure virtual methods - must be implemented by derived classes
        virtual double getElevation(double lat, double lon) = 0;
        virtual std::vector<TerrainPoint> getElevationGrid(double latStart, 
            double lonStart, double latEnd, double lonEnd, 
            int numLatSamples, int numLonSamples) = 0;

        // Common virtual methods with default implementations
        virtual void startServer(int port);
        virtual void setInterpolationTolerance(double tolerance) { interpolationTolerance_ = tolerance; }
        virtual double getInterpolationTolerance() const { return interpolationTolerance_; }
        virtual BoundingBox getBoundingBox() const { return bbox_; }
        virtual void exportToCSV(const std::vector<TerrainPoint>& points, const std::string& filename);

    protected:
        // Common member variables available to all derived classes
        std::mutex ioMutex_;
        double interpolationTolerance_ = 0.01; // Default tolerance ~1km
        BoundingBox bbox_ = {0, 0, 0, 0};

        // Common helper methods
        virtual void handleClient(int clientSocket);
        bool isPointInBounds(double lat, double lon) const;
        bool isPointInterpolatable(double lat, double lon) const;
};