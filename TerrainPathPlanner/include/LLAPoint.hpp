/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file LLAPoint.hpp
 * @brief LLAPoint struct definition
 */

#pragma once

#include <cmath>

/**
 * @brief LLAPoint struct
 * @details This struct represents a point on the terrain with latitude, longitude, and altitude information.
 */
struct LLAPoint
{
    double lat;
    double lon;
    double alt;
    
    LLAPoint() : lat(0), lon(0), alt(0) {}
    LLAPoint(double lat_, double lon_, double alt_) : lat(lat_), lon(lon_), alt(alt_) {}
    
    bool operator==(const LLAPoint& other) const
    {
        return std::abs(lat - other.lat) < 1e-6 &&
               std::abs(lon - other.lon) < 1e-6 &&
               std::abs(alt - other.alt) < 1e-6;
    }
};