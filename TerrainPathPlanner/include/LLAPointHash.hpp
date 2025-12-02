/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file LLAPointHash.hpp
 * @brief LLAPointHash struct definition
 */

#pragma once

#include <vector>
#include <cmath>

#include "LLAPoint.hpp"

// Hash function for LLAPoint (for use in unordered_map)
struct LLAPointHash
{
    std::size_t operator()(const LLAPoint& p) const
    {
        // Discretize to grid cells for hashing
        long long lat_grid = static_cast<long long>(p.lat * 1000000);
        long long lon_grid = static_cast<long long>(p.lon * 1000000);
        long long alt_grid = static_cast<long long>(p.alt * 10);
        
        std::size_t h1 = std::hash<long long>{}(lat_grid);
        std::size_t h2 = std::hash<long long>{}(lon_grid);
        std::size_t h3 = std::hash<long long>{}(alt_grid);
        
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};