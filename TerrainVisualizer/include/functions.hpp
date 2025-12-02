/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file functions.hpp
 * @brief Utility functions for converting degrees to radians
 */

#pragma once

#include <cmath>

/**
 * @brief Convert degrees to radians
 * @details Converts degrees to radians.
 * @param degrees Degrees to convert.
 * @return Radians.
 */
inline double degToRad(double degrees) { return degrees * M_PI / 180.0; }