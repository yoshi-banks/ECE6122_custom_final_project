/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file NEDPoint.hpp
 * @brief NEDPoint struct definition
 */

#pragma once

#include <cmath>

/**
 * @brief NEDPoint struct
 * @details This struct represents a point in North-East-Down (NED) coordinates.
 */
struct NEDPoint {
    double north;
    double east;
    double down;
    
    NEDPoint(double n = 0, double e = 0, double d = 0) 
        : north(n), east(e), down(d) {}
    
    bool operator==(const NEDPoint& other) const {
        return north == other.north && east == other.east && down == other.down;
    }
    
    double distanceTo(const NEDPoint& other) const {
        double dn = north - other.north;
        double de = east - other.east;
        double dd = down - other.down;
        return std::sqrt(dn*dn + de*de + dd*dd);
    }
};