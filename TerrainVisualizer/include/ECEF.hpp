/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file ECEF.hpp
 * @brief ECEF class definition
 * @details This file declares the ECEF class for representing a point in Earth Centered Earth Fixed (ECEF) coordinates.
 */

#pragma once

#include <cmath>

#include "functions.hpp"
#include "constants.hpp"

/**
 * @brief ECEF struct
 * 
 * @details This struct represents a point in Earth Centered Earth Fixed (ECEF) coordinates.
 */
struct ECEF
{
    double x;
    double y;
    double z;
};

/**
 * @brief Convert lat/lon/alt to ECEF
 * 
 * @details Converts latitude, longitude, and altitude to ECEF coordinates.
 * @param latDeg Latitude in degrees.
 * @param lonDeg Longitude in degrees.
 * @param alt Altitude in meters.
 */
inline ECEF latLonAltToECEF(double latDeg, double lonDeg, double alt)
{
    double latRad = degToRad(latDeg);
    double lonRad = degToRad(lonDeg);

    double sinLat = std::sin(latRad);
    double cosLat = std::cos(latRad);
    double sinLon = std::sin(lonRad);
    double cosLon = std::cos(lonRad);

    // Radius of curvature in the prime vertical
    double N = TerrainConstants::EARTH_RADIUS / std::sqrt(1.0 - TerrainConstants::EARTH_ECCENTRICITY_SQ * sinLat * sinLat);

    ECEF ecef;
    ecef.x = (N + alt) * cosLat * cosLon;
    ecef.y = (N + alt) * cosLat * sinLon;
    ecef.z = (N * (1.0 - TerrainConstants::EARTH_ECCENTRICITY_SQ) + alt) * sinLat;

    return ecef;
}
