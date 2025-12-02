/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file TerrainMesh.hpp
 * 
 * @brief TerrainMesh class definition
 * @details This file declares the TerrainMesh class for generating and rendering
 * terrain meshes from terrain data servers in the Terrain Visualizer application.
 */

#pragma once

#include <vector>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <iostream>

#include "TerrainServerOpenTopo.hpp"
#include "NEDPoint.hpp"

/**
 * @brief TerrainPoint struct
 * 
 * @details This struct represents a point on the terrain with position, normal, and color information.
 */
struct Vertex
{
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec3 color;
};

/**
 * @brief TerrainMesh class
 * 
 * @details This class represents a terrain mesh generated from terrain data.
 */
class TerrainMesh
{
public:
    TerrainMesh();
    
    ~TerrainMesh();

    void generateFromServer(ITerrainServer& server,
                           double latStart, double lonStart,
                           double latEnd, double lonEnd,
                           int rows, int cols);

    void generateDataOnly(ITerrainServer& server, 
                          double latStart, double lonStart,
                          double latEnd, double lonEnd,
                          int rows, int cols);

    void render(bool wireframe = false);

    glm::vec3 getCenterPoint() const { return centerPoint_; }

    glm::vec3 latLonToMeshCoords(double lat, double lon, double alt) const;

    double getElevationAtLatLon(double lat, double lon) const;

    double getRefLat() const { return ref_lat_; }
    double getRefLon() const { return ref_lon_; }
    double getRefAlt() const { return ref_alt_; }
    double getLatStart() const { return latStart_; }
    double getLatEnd() const { return latEnd_; }
    double getLonStart() const { return lonStart_; }
    double getLonEnd() const { return lonEnd_; }
    int getRows() const { return rows_; }
    int getCols() const { return cols_; }
    double getMinElev() const { return minElev_; }
    double getMaxElev() const { return maxElev_; }

    // Get terrain height grid for path planning
    std::vector<std::vector<double>> getHeightGrid() const { return height_grid_; }
    double getNorthMin() const { return north_min_; }
    double getNorthMax() const { return north_max_; }
    double getEastMin() const { return east_min_; }
    double getEastMax() const { return east_max_; }

    // Convert lat/lon to NEDPoint
    NEDPoint latLonToNED(double lat, double lon, double alt) const;
    
    // Convert NEDPoint to lat/lon
    void nedToLatLon(const NEDPoint& ned, double& lat, double& lon, double& alt) const;

    // Helper to create visible marker at lat/lon with altitude offset
    NEDPoint latLonToNEDWithOffset(double lat, double lon, double offset_agl) const;

    void verifyCoordinateTransform(double lat, double lon, double alt) const;

private:
    void calculateNormals(std::vector<Vertex>& vertices,
                         const std::vector<unsigned int>& indices,
                         int cols);

    glm::vec3 getTerrainColor(float height);

    void createBuffers(const std::vector<Vertex>& vertices,
                      const std::vector<unsigned int>& indices);

    void cleanup();

    GLuint vao_;
    GLuint vbo_;
    GLuint ebo_;
    unsigned int indexCount_;
    glm::vec3 centerPoint_;

    double latStart_, lonStart_, latEnd_, lonEnd_;
    int rows_, cols_;
    double minElev_, maxElev_;
    std::vector<TerrainPoint> points_; // store terrain grid points

    std::vector<std::vector<double>> height_grid_;  // Height at each grid point
    double north_min_, north_max_, east_min_, east_max_;
    double ref_lat_, ref_lon_, ref_alt_;  // Reference point for NED conversion
};