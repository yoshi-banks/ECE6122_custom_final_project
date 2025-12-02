/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file constants.hpp
 * @brief Terrain-related constants
 */

#pragma once

namespace TerrainConstants
{
    constexpr double EARTH_RADIUS = 6378137.0;
    constexpr double EARTH_FLATTENING = 1.0 / 298.257223563;
    constexpr double EARTH_ECCENTRICITY_SQ = 2.0 * EARTH_FLATTENING - EARTH_FLATTENING * EARTH_FLATTENING;
}