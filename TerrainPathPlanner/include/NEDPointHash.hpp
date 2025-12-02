/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file NEDPointHash.hpp
 * @brief NEDPointHash struct definition
 */

#pragma once

#include <cmath>
#include <vector>

#include "NEDPoint.hpp"

/**
 * @brief NEDPointHash struct
 * @details This struct represents a hash function for NEDPoint (for use in unordered_map).
 */
struct NEDPointHash {
    size_t operator()(const NEDPoint& p) const {
        // Quantize to grid resolution for hashing
        int64_t n = static_cast<int64_t>(p.north * 10.0);
        int64_t e = static_cast<int64_t>(p.east * 10.0);
        int64_t d = static_cast<int64_t>(p.down * 10.0);
        return std::hash<int64_t>()(n) ^ 
               (std::hash<int64_t>()(e) << 1) ^ 
               (std::hash<int64_t>()(d) << 2);
    }
};