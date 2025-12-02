/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file TerrainMesh.cpp
 * @brief TerrainMesh class implementation
 * @details This file implements the TerrainMesh class for generating and rendering
 * terrain meshes from terrain data servers in the Terrain Visualizer application.
 */

#include <iomanip>
#include <iostream>

#include "TerrainMesh.hpp"
#include "constants.hpp"
#include "NED.hpp"

/**
 * @brief Constructor for TerrainMesh class
 */
TerrainMesh::TerrainMesh() 
    : vao_(0)
    , vbo_(0)
    , ebo_(0)
    , indexCount_(0) 
{

}

/**
 * @brief Destructor for TerrainMesh class
 */
TerrainMesh::~TerrainMesh()
{
    cleanup();
}

/**
 * @brief Generates a terrain mesh from terrain data server
 * @param server Terrain data server
 * @param latStart Starting latitude
 * @param lonStart Starting longitude
 * @param latEnd Ending latitude
 * @param lonEnd Ending longitude
 * @param rows Number of rows in the grid
 * @param cols Number of columns in the grid
 */
void TerrainMesh::generateFromServer(ITerrainServer& server,
                        double latStart, double lonStart,
                        double latEnd, double lonEnd,
                        int rows, int cols)
{
    std::cout << "Fetching terrain data...\n";
    std::vector<TerrainPoint> points = server.getElevationGrid(latStart, lonStart, latEnd, lonEnd, rows, cols);

    // Calculate refernece point (center of region)
    double refLat = (latStart + latEnd) / 2.0;
    double refLon = (lonStart + lonEnd) / 2.0;
    
    std::cout << "Finding elevation range" << std::endl;
    // Find elevation range for coloring
    double minElev = 1e9, maxElev = -1e9;
    for (const auto& p : points)
    {
        minElev = std::min(minElev, p.alt);
        maxElev = std::max(maxElev, p.alt);
    }

    // Use reference altitude as the mean elevation
    double refAlt = (minElev + maxElev) / 2.0;

    double elevRange = maxElev - minElev;
    std::cout << "Elevation range: " << minElev << " to " << maxElev << " m\n";
    std::cout << "Reference point: lat=" << refLat << " lon=" << refLon << " alt=" << refAlt << "\n";

    minElev_ = minElev;
    maxElev_ = maxElev;
    latStart_ = latStart;
    lonStart_ = lonStart;
    latEnd_ = latEnd;
    lonEnd_ = lonEnd;
    rows_ = rows;
    cols_ = cols;
    points_ = points;

    // Store reference for later conversions
    ref_lat_ = refLat;
    ref_lon_ = refLon;
    ref_alt_ = refAlt;
    
    // Build height grid for path planning
    height_grid_.resize(rows, std::vector<double>(cols));
    
    // Convert all points to NED and build height grid
    std::vector<NED> nedPoints(rows * cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int idx = i * cols + j;
            const auto& pt = points[idx];
            nedPoints[idx] = latLonAltToNED(pt.lat, pt.lon, pt.alt, refLat, refLon, refAlt);
            
            // Store the down component (altitude) in grid
            height_grid_[i][j] = nedPoints[idx].down;
        }
    }
    
    // Track NED bounds
    north_min_ = nedPoints[0].north;
    north_max_ = nedPoints[0].north;
    east_min_ = nedPoints[0].east;
    east_max_ = nedPoints[0].east;
    
    for (const auto& ned : nedPoints) {
        north_min_ = std::min(north_min_, ned.north);
        north_max_ = std::max(north_max_, ned.north);
        east_min_ = std::min(east_min_, ned.east);
        east_max_ = std::max(east_max_, ned.east);
    }

    // Smooth the down (altitude component to remove spikes)
    std::cout << "Smoothing altitude data...\n";
    std::vector<double> smoothedDown = medianFilterAltitude(nedPoints, rows, cols, 3); // 3x3 kernel

    // Generate vertices using smoothed NED coordinates
    std::vector<Vertex> vertices;
    vertices.reserve(rows * cols);

    // Track NED bounds for centering
    double minNorth = 1e9, maxNorth = -1e9;
    double minEast = 1e9, maxEast = -1e9;
    double minDown = 1e9, maxDown = -1e9;

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            int idx = i * cols + j;
            const auto& pt = points[idx];
            const NED& ned = nedPoints[idx];

            // Use smoothed down vlaue
            double smoothDown = smoothedDown[idx];

            // Track bounds
            minNorth = std::min(minNorth, ned.north);
            maxNorth = std::max(maxNorth, ned.north);
            minEast = std::min(minEast, ned.east);
            maxEast = std::max(maxEast, ned.east);
            minDown = std::min(minDown, ned.down);
            maxDown = std::max(maxDown, ned.down);

            // Use NED directly - note: in graphics typcially Y is up
            // So we map: North->X, -Down->Y (up), East->Z
            float x = static_cast<float>(ned.north);
            float y = static_cast<float>(-smoothDown); // Negative because down is positive downward
            float z = static_cast<float>(ned.east);

            // Height-based coloring
            float normHeight = (elevRange > 0) ? (pt.alt - minElev) / elevRange : 0.5f;
            glm::vec3 color = getTerrainColor(normHeight);

            vertices.push_back({
                glm::vec3(x, y, z),
                glm::vec3(0.0f, 1.0f, 0.0f), // Will calculate normals later
                color
            });
        }
    }

    std::cout << "NED bounds:\n"
              << "  North: [" << minNorth << ", " << maxNorth << "] m\n"
              << "  East:  [" << minEast << ", " << maxEast << "] m\n"
              << "  Down:  [" << minDown << ", " << maxDown << "] m\n";

    // Generate indices for triangles
    std::vector<unsigned int> indices;
    indices.reserve((rows - 1) * (cols - 1) * 6);

    for (int i = 0; i < rows - 1; ++i)
    {
        for (int j = 0; j < cols - 1; ++j)
        {
            unsigned int topLeft = i * cols + j;
            unsigned int topRight = topLeft + 1;
            unsigned int bottomLeft = (i + 1) * cols + j;
            unsigned int bottomRight = bottomLeft + 1;

            // First triangle
            indices.push_back(topLeft);
            indices.push_back(bottomLeft);
            indices.push_back(topRight);

            // Second triangle
            indices.push_back(topRight);
            indices.push_back(bottomLeft);
            indices.push_back(bottomRight);
        }
    }

    // Calculate normals
    calculateNormals(vertices, indices, cols);

    // Create OpenGL buffers
    createBuffers(vertices, indices);

    // Center point in NED Coordaintes
    centerPoint_ = glm::vec3(
        (minNorth + maxNorth) / 2.0f,
        -(minDown + maxDown) / 2.0f,
        (minEast + maxEast) / 2.0f
    );

    std::cout << "Generated terrain mesh: " << vertices.size() << " vertices, "
                << indices.size() / 3 << " triangles\n";
    std::cout << "Center point (NED): " << centerPoint_.x << ", "
              << centerPoint_.y << ", " << centerPoint_.z << "m\n";
}

void TerrainMesh::generateDataOnly(ITerrainServer& server,
                                  double latStart, double lonStart,
                                  double latEnd, double lonEnd,
                                  int rows, int cols) {
    std::cout << "Fetching terrain data (data-only mode, no OpenGL)...\n";
    std::vector<TerrainPoint> points = server.getElevationGrid(latStart, lonStart, latEnd, lonEnd, rows, cols);

    // Calculate reference point (center of region)
    double refLat = (latStart + latEnd) / 2.0;
    double refLon = (lonStart + lonEnd) / 2.0;
    
    std::cout << "Finding elevation range" << std::endl;
    // Find elevation range
    double minElev = 1e9, maxElev = -1e9;
    for (const auto& p : points)
    {
        minElev = std::min(minElev, p.alt);
        maxElev = std::max(maxElev, p.alt);
    }

    // Use reference altitude as the mean elevation
    double refAlt = (minElev + maxElev) / 2.0;

    std::cout << "Elevation range: " << minElev << " to " << maxElev << " m\n";
    std::cout << "Reference point: lat=" << refLat << " lon=" << refLon << " alt=" << refAlt << "\n";

    minElev_ = minElev;
    maxElev_ = maxElev;
    latStart_ = latStart;
    lonStart_ = lonStart;
    latEnd_ = latEnd;
    lonEnd_ = lonEnd;
    rows_ = rows;
    cols_ = cols;
    points_ = points;

    // Store reference for later conversions
    ref_lat_ = refLat;
    ref_lon_ = refLon;
    ref_alt_ = refAlt;
    
    // Build height grid for path planning
    height_grid_.resize(rows, std::vector<double>(cols));
    
    // Convert all points to NED and build height grid
    std::vector<NED> nedPoints(rows * cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int idx = i * cols + j;
            const auto& pt = points[idx];
            nedPoints[idx] = latLonAltToNED(pt.lat, pt.lon, pt.alt, refLat, refLon, refAlt);
            
            // Store the down component (altitude) in grid
            height_grid_[i][j] = nedPoints[idx].down;
        }
    }
    
    // Track NED bounds
    north_min_ = nedPoints[0].north;
    north_max_ = nedPoints[0].north;
    east_min_ = nedPoints[0].east;
    east_max_ = nedPoints[0].east;
    
    for (const auto& ned : nedPoints) {
        north_min_ = std::min(north_min_, ned.north);
        north_max_ = std::max(north_max_, ned.north);
        east_min_ = std::min(east_min_, ned.east);
        east_max_ = std::max(east_max_, ned.east);
    }

    std::cout << "NED bounds:\n"
              << "  North: [" << north_min_ << ", " << north_max_ << "] m\n"
              << "  East:  [" << east_min_ << ", " << east_max_ << "] m\n";

    // Calculate center point (for consistency, though not used for rendering)
    double minDown = 1e9, maxDown = -1e9;
    for (const auto& ned : nedPoints) {
        minDown = std::min(minDown, ned.down);
        maxDown = std::max(maxDown, ned.down);
    }
    
    centerPoint_ = glm::vec3(
        (north_min_ + north_max_) / 2.0f,
        -(minDown + maxDown) / 2.0f,
        (east_min_ + east_max_) / 2.0f
    );

    std::cout << "Generated terrain data (no mesh): " << points.size() << " points\n";
    std::cout << "Center point (NED): " << centerPoint_.x << ", "
              << centerPoint_.y << ", " << centerPoint_.z << "m\n";
}

/**
 * @brief Convert lat/lon/alt to mesh coordinates
 * @param lat Latitude in degrees
 * @param lon Longitude in degrees
 * @param alt Altitude in meters
 * @return glm::vec3 Mesh coordinates
 */
NEDPoint TerrainMesh::latLonToNED(double lat, double lon, double alt) const 
{
    NED ned = latLonAltToNED(lat, lon, alt, ref_lat_, ref_lon_, ref_alt_);
    return NEDPoint(ned.north, ned.east, ned.down);
}

/**
 * @brief Convert lat/lon with altitude offset to NEDPoint
 * @param lat Latitude in degrees
 * @param lon Longitude in degrees
 * @param offset_agl Altitude offset in meters above ground level
 * @return NEDPoint Corresponding NED coordinates
 */
NEDPoint TerrainMesh::latLonToNEDWithOffset(double lat, double lon, double offset_agl) const {
    // Get terrain elevation at this lat/lon
    double terrain_alt = getElevationAtLatLon(lat, lon);
    
    // Calculate absolute altitude MSL (Mean Sea Level)
    double marker_alt_msl = terrain_alt + offset_agl;
    
    // Convert to NED
    NED ned = latLonAltToNED(lat, lon, marker_alt_msl, ref_lat_, ref_lon_, ref_alt_);
    
    // IMPORTANT: Verify the sign is correct
    // If offset_agl is positive (above ground), ned.down should be NEGATIVE
    // because we're going UP from the reference altitude
    
    std::cout << "DEBUG latLonToNEDWithOffset:\n"
              << "  Input: lat=" << lat << ", lon=" << lon << ", offset_agl=" << offset_agl << "\n"
              << "  Terrain alt: " << terrain_alt << " m MSL\n"
              << "  Marker alt: " << marker_alt_msl << " m MSL\n"
              << "  Reference alt: " << ref_alt_ << " m MSL\n"
              << "  NED down: " << ned.down << " (should be negative if above reference)\n";
    
    return NEDPoint(ned.north, ned.east, ned.down);
}

/**
 * @brief Convert NED displacements to lat/lon/alt changes
 * @param ned NED displacements
 * @param lat Reference latitude in degrees
 * @param lon Reference longitude in degrees
 * @param alt Reference altitude in meters
 * @return lat, lon, alt Corresponding latitude, longitude, and altitude
 */
void TerrainMesh::nedToLatLon(const NEDPoint& ned, double& lat, double& lon, double& alt) const 
{    
    double lat_rad = ref_lat_ * M_PI / 180.0;
    double lon_rad = ref_lon_ * M_PI / 180.0;

    // Radius of curvature in meridian
    double sin_lat = std::sin(lat_rad);
    double N = TerrainConstants::EARTH_RADIUS / std::sqrt(1.0 - TerrainConstants::EARTH_ECCENTRICITY_SQ * sin_lat * sin_lat);
    double M = TerrainConstants::EARTH_RADIUS * (1.0 - TerrainConstants::EARTH_ECCENTRICITY_SQ) 
        / std::pow(1.0 - TerrainConstants::EARTH_ECCENTRICITY_SQ * sin_lat * sin_lat, 1.5);
    
    // Convert NED displacements to lat/lon/alt changes
    double d_lat_rad = ned.north / M;
    double d_lon_rad = ned.east / (N * std::cos(lat_rad));
    double d_alt = -ned.down;  // Down is positive downward, altitude is positive upward
    
    lat = ref_lat_ + d_lat_rad * 180.0 / M_PI;
    lon = ref_lon_ + d_lon_rad * 180.0 / M_PI;
    alt = ref_alt_ + d_alt;
}

/**
 * @brief Render the mesh
 * @param wireframe Render wireframe if true, solid if false
 */
void TerrainMesh::render(bool wireframe)
{
    if (vao_ == 0) return;

    glBindVertexArray(vao_);
    
    if (wireframe)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
    else
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    glDrawElements(GL_TRIANGLES, indexCount_, GL_UNSIGNED_INT, 0);
    
    glBindVertexArray(0);
}

/**
 * @brief Calculate normals for the terrain mesh
 * @param vertices Vertex list
 * @param indices Index list
 * @param cols Number of columns in the grid
 */
void TerrainMesh::calculateNormals(std::vector<Vertex>& vertices,
                        const std::vector<unsigned int>& indices,
                        int cols)
{
    // Reset normals
    for (auto& v : vertices)
    {
        v.normal = glm::vec3(0.0f);
    }

    // Accumulate face normals
    for (size_t i = 0; i < indices.size(); i += 3)
    {
        unsigned int i0 = indices[i];
        unsigned int i1 = indices[i + 1];
        unsigned int i2 = indices[i + 2];

        glm::vec3 v0 = vertices[i0].position;
        glm::vec3 v1 = vertices[i1].position;
        glm::vec3 v2 = vertices[i2].position;

        glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));

        vertices[i0].normal += normal;
        vertices[i1].normal += normal;
        vertices[i2].normal += normal;
    }

    // Normalize
    for (auto& v : vertices)
    {
        v.normal = glm::normalize(v.normal);
    }
}

/**
 * @brief Get the color of the terrain based on height
 * @param height Height of the terrain (0.0 to 1.0 normalized)
 */
glm::vec3 TerrainMesh::getTerrainColor(float height)
{
    // Realistic terrain coloring
    if (height < 0.2f)
        return glm::mix(glm::vec3(0.2f, 0.5f, 0.3f), glm::vec3(0.3f, 0.6f, 0.3f), height / 0.2f); // Dark to light green
    else if (height < 0.5f)
        return glm::mix(glm::vec3(0.3f, 0.6f, 0.3f), glm::vec3(0.6f, 0.5f, 0.3f), (height - 0.2f) / 0.3f); // Green to brown
    else if (height < 0.7f)
        return glm::mix(glm::vec3(0.6f, 0.5f, 0.3f), glm::vec3(0.5f, 0.5f, 0.5f), (height - 0.5f) / 0.2f); // Brown to gray
    else
        return glm::mix(glm::vec3(0.5f, 0.5f, 0.5f), glm::vec3(1.0f, 1.0f, 1.0f), (height - 0.7f) / 0.3f); // Gray to white (snow)
}

/**
 * @brief Create buffers for the mesh
 * @param vertices Vertex list
 * @param indices Index list
 */
void TerrainMesh::createBuffers(const std::vector<Vertex>& vertices,
                    const std::vector<unsigned int>& indices)
{
    cleanup(); // Clean up old buffers if any

    indexCount_ = indices.size();

    // Generate and bind VAO
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);

    // Generate and bind VBO
    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex),
                vertices.data(), GL_STATIC_DRAW);

    // Generate and bind EBO
    glGenBuffers(1, &ebo_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int),
                indices.data(), GL_STATIC_DRAW);

    // Position attribute
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                        (void*)offsetof(Vertex, position));

    // Normal attribute
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                        (void*)offsetof(Vertex, normal));

    // Color attribute
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                        (void*)offsetof(Vertex, color));

    glBindVertexArray(0);
}

/**
 * @brief Clean up OpenGL buffers
 */
void TerrainMesh::cleanup()
{
    if (vbo_) glDeleteBuffers(1, &vbo_);
    if (ebo_) glDeleteBuffers(1, &ebo_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
    vbo_ = ebo_ = vao_ = 0;
}

/**
 * @brief Get the elevation at a given latitude and longitude
 * @param lat Latitude in degrees
 * @param lon Longitude in degrees
 * @return Elevation in meters
 */
double TerrainMesh::getElevationAtLatLon(double lat, double lon) const
{
    // Simple nearest-neighbor for now
    double minDist = 1e9;
    double alt = minElev_;
    for (const auto& p : points_)
    {
        double d = (p.lat - lat)*(p.lat - lat) + (p.lon - lon)*(p.lon - lon);
        if (d < minDist)
        {
            minDist = d;
            alt = p.alt;
        }
    }
    return alt;
}

/**
 * @brief Verify coordinate transformations by round-tripping a lat/lon/alt point
 * @param lat Latitude in degrees
 * @param lon Longitude in degrees
 * @param alt Altitude in meters
 */
void TerrainMesh::verifyCoordinateTransform(double lat, double lon, double alt) const
{
    std::cout << "\n=== COORDINATE TRANSFORM VERIFICATION ===\n";
    std::cout << std::fixed << std::setprecision(8);
    
    // Forward: LLA → NED
    NED ned = latLonAltToNED(lat, lon, alt, ref_lat_, ref_lon_, ref_alt_);
    std::cout << "Input LLA:\n"
              << "  Lat: " << lat << "°\n"
              << "  Lon: " << lon << "°\n"
              << "  Alt: " << alt << " m MSL\n\n";
    
    std::cout << "NED Result:\n"
              << "  North: " << ned.north << " m\n"
              << "  East:  " << ned.east << " m\n"
              << "  Down:  " << ned.down << " m\n\n";
    
    // Backward: NED → LLA
    double lat_back, lon_back, alt_back;
    NEDPoint ned_pt(ned.north, ned.east, ned.down);
    nedToLatLon(ned_pt, lat_back, lon_back, alt_back);
    
    std::cout << "Round-trip LLA:\n"
              << "  Lat: " << lat_back << "°\n"
              << "  Lon: " << lon_back << "°\n"
              << "  Alt: " << alt_back << " m MSL\n\n";
    
    // Calculate errors
    double lat_error_m = (lat_back - lat) * 111000.0;
    double lon_error_m = (lon_back - lon) * 111000.0 * std::cos(lat * M_PI / 180.0);
    double alt_error_m = alt_back - alt;
    
    std::cout << "Errors:\n"
              << "  Latitude:  " << lat_error_m << " m\n"
              << "  Longitude: " << lon_error_m << " m\n"
              << "  Altitude:  " << alt_error_m << " m\n";
    
    double total_error = std::sqrt(lat_error_m*lat_error_m + 
                                   lon_error_m*lon_error_m + 
                                   alt_error_m*alt_error_m);
    std::cout << "  Total 3D:  " << total_error << " m\n";
    
    if (total_error > 1.0) {
        std::cout << "\n⚠️  WARNING: Round-trip error > 1m! Transform may be incorrect.\n";
    } else {
        std::cout << "\n✓ Transform appears correct (error < 1m)\n";
    }
    std::cout << "=========================================\n\n";
}