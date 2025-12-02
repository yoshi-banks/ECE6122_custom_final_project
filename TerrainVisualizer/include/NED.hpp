/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file NED.hpp
 * @brief NED class definition
 */

#pragma once

#include <cmath>
#include <vector>
#include <algorithm>

#include "constants.hpp"
#include "functions.hpp"
#include "ECEF.hpp"

/**
 * @brief NED struct
 * @details This struct represents a point in North-East-Down (NED) coordinates.
 */
struct NED
{
    double north;
    double east;
    double down;

    NED() : north(0), east(0), down(0) {}
    NED(double n, double e, double d) : north(n), east(e), down(d) {}
};

/**
 * @brief Convert ECEF to NED
 * @details Converts a point from ECEF coordinates to NED coordinates relative to a reference point.
 * @param point The point in ECEF coordinates.
 * @param reference The reference point in ECEF coordinates.
 * @param refLat The latitude of the reference point in degrees.
 * @param refLon The longitude of the reference point in degrees.
 * @return The point in NED coordinates.
 */
inline NED ecefToNed(const ECEF& point, const ECEF& reference, double refLat, double refLon)
{
    double latRad = degToRad(refLat);
    double lonRad = degToRad(refLon);

    double sinLat = std::sin(latRad);
    double cosLat = std::cos(latRad);
    double sinLon = std::sin(lonRad);
    double cosLon = std::cos(lonRad);

    // Difference3 vector in ECEF
    double dx = point.x - reference.x;
    double dy = point.y - reference.y;
    double dz = point.z - reference.z;

    // Rotation matrix from ECEF to NED
    NED ned;
    ned.north = -sinLat * cosLon * dx - sinLat * sinLon * dy + cosLat * dz;
    ned.east  = -sinLon * dx + cosLon * dy; 
    ned.down  = -cosLat * cosLon * dx - cosLat * sinLon * dy - sinLat * dz;

    return ned;
}

/**
 * @brief Convert lat/lon/alt to NED
 * @details Converts latitude, longitude, and altitude to NED coordinates relative to a reference location.
 * @param latDeg Latitude in degrees.
 * @param lonDeg Longitude in degrees.
 * @param alt Altitude in meters.
 * @param refLatDeg Reference latitude in degrees.
 * @param refLonDeg Reference longitude in degrees.
 * @param refAlt Reference altitude in meters.
 */
inline NED latLonAltToNED(double latDeg, double lonDeg, double alt,
                   double refLatDeg, double refLonDeg, double refAlt)
{
    ECEF point = latLonAltToECEF(latDeg, lonDeg, alt);
    ECEF reference = latLonAltToECEF(refLatDeg, refLonDeg, refAlt);
    return ecefToNed(point, reference, refLatDeg, refLonDeg);
}

/**
 * @brief Median filter for altitude
 * 
 * @details Applies a median filter to the down (altitude) component of NED points to smooth out noise.
 * @param nedPoints Vector of NED points.
 * @param rows Number of rows in the grid.
 * @param cols Number of columns in the grid.
 * @param kernelSize Size of the median filter kernel.
 */
inline std::vector<double> medianFilterAltitude(const std::vector<NED>& nedPoints,
                                         int rows, int cols, int kernelSize)
{
    std::vector<double> smoothed(rows * cols);
    int halfKernel = kernelSize / 2;

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            int idx = i * cols + j;

            // Collect neighborhood values
            std::vector<double> neighborhood;
            neighborhood.reserve(kernelSize * kernelSize);

            for (int di = -halfKernel; di <= halfKernel; ++di)
            {
                for (int dj = -halfKernel; dj <= halfKernel; ++dj)
                {
                    int ni = i + di;
                    int nj = j + dj;

                    // Clamp to boundaries
                    ni = std::max(0, std::min(rows - 1, ni));
                    nj = std::max(0, std::min(cols - 1, nj));

                    int nidx = ni * cols + nj;
                    neighborhood.push_back(nedPoints[nidx].down);
                }
            }

            // Find median
            std::nth_element(neighborhood.begin(), neighborhood.begin() + neighborhood.size() / 2, neighborhood.end());
            smoothed[idx] = neighborhood[neighborhood.size() / 2];
        }
    }

    return smoothed;
}

/**
 * @brief Gaussian filter for altitude
 * 
 * @details Applies a Gaussian filter to the down (altitude) component of NED points to smooth out noise.
 * @param nedPoints Vector of NED points.
 * @param rows Number of rows in the grid.
 * @param cols Number of columns in the grid.
 * @param sigma Standard deviation of the Gaussian kernel.
 */
inline std::vector<double> gaussianFilterAltitude(const std::vector<NED>& nedPoints,
                                          int rows, int cols, double sigma)
{
    std::vector<double> smoothed(rows * cols);
    
    // Create Gaussian kernel
    int kernelSize = static_cast<int>(std::ceil(3.0 * sigma)) * 2 + 1;
    int halfKernel = kernelSize / 2;
    std::vector<double> kernel(kernelSize);
    double sum = 0.0;
    
    for (int i = 0; i < kernelSize; ++i) {
        int x = i - halfKernel;
        kernel[i] = std::exp(-(x * x) / (2.0 * sigma * sigma));
        sum += kernel[i];
    }
    
    // Normalize kernel
    for (auto& k : kernel) {
        k /= sum;
    }
    
    // Apply separable filter (horizontal then vertical for efficiency)
    std::vector<double> temp(rows * cols);
    
    // Horizontal pass
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double value = 0.0;
            double weightSum = 0.0;
            
            for (int k = -halfKernel; k <= halfKernel; ++k) {
                int nj = j + k;
                if (nj >= 0 && nj < cols) {
                    int nidx = i * cols + nj;
                    value += nedPoints[nidx].down * kernel[k + halfKernel];
                    weightSum += kernel[k + halfKernel];
                }
            }
            
            temp[i * cols + j] = value / weightSum;
        }
    }
    
    // Vertical pass
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double value = 0.0;
            double weightSum = 0.0;
            
            for (int k = -halfKernel; k <= halfKernel; ++k) {
                int ni = i + k;
                if (ni >= 0 && ni < rows) {
                    int nidx = ni * cols + j;
                    value += temp[nidx] * kernel[k + halfKernel];
                    weightSum += kernel[k + halfKernel];
                }
            }
            
            smoothed[i * cols + j] = value / weightSum;
        }
    }
    
    return smoothed;
}