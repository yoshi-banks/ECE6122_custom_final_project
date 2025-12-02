/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file TerrainServerGeoTIFF.cpp
 * @brief TerrainServerGeoTIFF class implementation
 */

#include "TerrainServerGeoTIFF.hpp"
#include "Logger.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>

extern Logger logger;

/**
 * @brief Constructor for TerrainServerGeoTIFF class
 */
TerrainServerGeoTIFF::TerrainServerGeoTIFF(const std::string& geotiff_path,
                                           double interpolation_tolerance)
    : geotiff_path_(geotiff_path)
    , dataset_(nullptr)
    , band_(nullptr)
    , raster_width_(0)
    , raster_height_(0)
    , nodata_value_(-9999.0)
    , has_nodata_(false)
{
    interpolationTolerance_ = interpolation_tolerance;
    
    logger.log("[INFO] Initializing TerrainServerGeoTIFF with file: " + geotiff_path);
    
    // Initialize GDAL
    GDALAllRegister();
    
    // Load the GeoTIFF
    if (!loadGeoTIFF(geotiff_path))
    {
        logger.log("[ERROR] Failed to load GeoTIFF file: " + geotiff_path);
        std::cerr << "Failed to load GeoTIFF file: " << geotiff_path << std::endl;
    }
    else
    {
        logger.log("[INFO] GeoTIFF loaded successfully. Bounds: (" + 
                  std::to_string(bbox_.min_lat) + ", " + std::to_string(bbox_.min_lon) + ") to (" +
                  std::to_string(bbox_.max_lat) + ", " + std::to_string(bbox_.max_lon) + ")");
                  
    }
}

/**
 * @brief Destructor for TerrainServerGeoTIFF class
 */
TerrainServerGeoTIFF::~TerrainServerGeoTIFF()
{
    if (dataset_ != nullptr)
    {
        GDALClose(dataset_);
        logger.log("[INFO] GeoTIFF dataset closed.");
    }
}

/**
 * @brief Load GeoTIFF file and initialize members
 * @param path Path to GeoTIFF file
 * @return true if loaded successfully, false otherwise
 */
bool TerrainServerGeoTIFF::loadGeoTIFF(const std::string& path)
{
    dataset_ = (GDALDataset*)GDALOpen(path.c_str(), GA_ReadOnly);
    if (dataset_ == nullptr)
    {
        return false;
    }

    // Get raster dimensions
    raster_width_ = dataset_->GetRasterXSize();
    raster_height_ = dataset_->GetRasterYSize();

    logger.log("[INFO] Raster dimensions: " + std::to_string(raster_width_) + " x " + 
              std::to_string(raster_height_));

    // Get geotransform
    if (dataset_->GetGeoTransform(geoTransform_) != CE_None)
    {
        logger.log("[ERROR] Failed to get geotransform from GeoTIFF");
        return false;
    }

    // Get first raster band (assuming single-band elevation data)
    band_ = dataset_->GetRasterBand(1);
    if (band_ == nullptr)
    {
        logger.log("[ERROR] Failed to get raster band from GeoTIFF");
        return false;
    }

    // Check for nodata value
    int success;
    nodata_value_ = band_->GetNoDataValue(&success);
    has_nodata_ = (success != 0);
    
    if (has_nodata_)
    {
        logger.log("[INFO] NoData value: " + std::to_string(nodata_value_));
    }

    // Calculate bounding box
    calculateBoundingBox();

    calculateMaxElevation();

    return true;
}

/**
 * @brief Calculate the bounding box of the GeoTIFF
 */
void TerrainServerGeoTIFF::calculateBoundingBox()
{
    // Calculate lat/lon bounds from geotransform
    // GeoTransform: [0] = top left x, [1] = w-e pixel resolution, [2] = rotation (0 if north up)
    //               [3] = top left y, [4] = rotation (0 if north up), [5] = n-s pixel resolution (negative)
    
    double min_x = geoTransform_[0];
    double max_x = geoTransform_[0] + raster_width_ * geoTransform_[1];
    double max_y = geoTransform_[3];
    double min_y = geoTransform_[3] + raster_height_ * geoTransform_[5];

    // Assuming the GeoTIFF is in lat/lon coordinates (WGS84)
    // If it's in a different projection, we'd need to transform
    bbox_.min_lon = min_x;
    bbox_.max_lon = max_x;
    bbox_.min_lat = min_y;
    bbox_.max_lat = max_y;
}

/**
 * @brief Calculate the maximum elevation in the GeoTIFF
 */
void TerrainServerGeoTIFF::calculateMaxElevation()
{
    if (band_ == NULL)
    {
        logger.log("[WARNING] band_ not set");
        max_elevation_ = 0.0;
    }

    double minVal, maxVal, meanVal, stdDev;
    if (band_->GetStatistics(FALSE, TRUE, &minVal, &maxVal, &meanVal, &stdDev) == CE_None)
    {
        max_elevation_ = maxVal;
        logger.log("[INFO] Max elevation found: " + std::to_string(max_elevation_) + " m");
    }
    else
    {
        logger.log("[WARN] Could not compute statistics; using default max elevation of 0");
        max_elevation_ = 0.0;
    }
}

/**
 * @brief Convert lat/lon to pixel coordinates
 * @param lat Latitude in degrees
 * @param lon Longitude in degrees
 * @param pixel_x Pixel x-coordinate
 * @param pixel_y Pixel y-coordinate
 * @return true if successful, false otherwise
 */
bool TerrainServerGeoTIFF::latLonToPixel(double lat, double lon, double& pixel_x, double& pixel_y) const
{
    // Convert lat/lon to pixel coordinates using inverse geotransform
    // Assuming no rotation (geoTransform_[2] and [4] are 0)
    
    pixel_x = (lon - geoTransform_[0]) / geoTransform_[1];
    pixel_y = (lat - geoTransform_[3]) / geoTransform_[5];
    
    return true;
}

/**
 * @brief Convert pixel coordinates to lat/lon
 * @param pixel_x Pixel x-coordinate
 * @param pixel_y Pixel y-coordinate
 * @param lat Latitude in degrees
 * @param lon Longitude in degrees
 * @return true if successful, false otherwise
 */
bool TerrainServerGeoTIFF::pixelToLatLon(double pixel_x, double pixel_y, double& lat, double& lon) const
{
    lon = geoTransform_[0] + pixel_x * geoTransform_[1];
    lat = geoTransform_[3] + pixel_y * geoTransform_[5];
    
    return true;
}

/**
 * @brief Perform bilinear interpolation to get elevation at fractional pixel coordinates
 * @param pixel_x Pixel x-coordinate
 * @param pixel_y Pixel y-coordinate
 * @return Elevation in meters
 */
double TerrainServerGeoTIFF::bilinearInterpolation(double pixel_x, double pixel_y)
{
    // Get integer pixel coordinates
    int x0 = static_cast<int>(std::floor(pixel_x));
    int y0 = static_cast<int>(std::floor(pixel_y));
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    // Clamp to raster bounds
    x0 = std::max(0, std::min(x0, raster_width_ - 1));
    x1 = std::max(0, std::min(x1, raster_width_ - 1));
    y0 = std::max(0, std::min(y0, raster_height_ - 1));
    y1 = std::max(0, std::min(y1, raster_height_ - 1));

    // Read elevation values at four corners
    float z00, z01, z10, z11;
    band_->RasterIO(GF_Read, x0, y0, 1, 1, &z00, 1, 1, GDT_Float32, 0, 0);
    band_->RasterIO(GF_Read, x1, y0, 1, 1, &z01, 1, 1, GDT_Float32, 0, 0);
    band_->RasterIO(GF_Read, x0, y1, 1, 1, &z10, 1, 1, GDT_Float32, 0, 0);
    band_->RasterIO(GF_Read, x1, y1, 1, 1, &z11, 1, 1, GDT_Float32, 0, 0);

    // Check for nodata values
    if (has_nodata_)
    {
        if (std::abs(z00 - nodata_value_) < 1e-6 || 
            std::abs(z01 - nodata_value_) < 1e-6 ||
            std::abs(z10 - nodata_value_) < 1e-6 || 
            std::abs(z11 - nodata_value_) < 1e-6)
        {
            return 0.0; // Return 0 if any corner has nodata
        }
    }

    // If at exact pixel boundary, return that value
    if (x0 == x1 && y0 == y1)
    {
        return static_cast<double>(z00);
    }

    // Calculate interpolation weights
    double dx = pixel_x - x0;
    double dy = pixel_y - y0;

    // Bilinear interpolation
    double z0 = z00 * (1.0 - dx) + z01 * dx;
    double z1 = z10 * (1.0 - dx) + z11 * dx;
    double z = z0 * (1.0 - dy) + z1 * dy;

    return z;
}

/**
 * @brief Get elevation at given pixel coordinates
 * @param pixel_x Pixel x-coordinate
 * @param pixel_y Pixel y-coordinate
 * @return Elevation in meters
 */
double TerrainServerGeoTIFF::getElevationAtPixel(double pixel_x, double pixel_y)
{
    // Check if pixel is within raster bounds (with some tolerance for interpolation)
    if (pixel_x < -0.5 || pixel_x > raster_width_ - 0.5 ||
        pixel_y < -0.5 || pixel_y > raster_height_ - 0.5)
    {
        logger.log("[WARN] Pixel out of bounds - returning max elevation");
        return max_elevation_;
    }

    return bilinearInterpolation(pixel_x, pixel_y);
}

/**
 * @brief Get elevation from Topo server
 * @param lat Latitude in degrees
 * @param lon Longitude in degrees
 * @return Elevation in meters
 */
double TerrainServerGeoTIFF::getElevationFromTopo(double lat, double lon)
{
    if (pServerTopo_ != nullptr)
    {
        logger.log("[INFO] Getting elevation from Topo server");
        return pServerTopo_->getElevation(lat, lon);
    }
    else 
    {
        logger.log("[WARN] Topo server isn't assigned in geotiff server");
        return 0.0;
    }
}

/**
 * @brief Get elevation at given latitude and longitude
 * @param lat Latitude in degrees
 * @param lon Longitude in degrees
 * @return Elevation in meters
 */
double TerrainServerGeoTIFF::getElevation(double lat, double lon)
{
    std::lock_guard<std::mutex> lock(ioMutex_);
    
    if (dataset_ == nullptr)
    {
        logger.log("[ERROR] GeoTIFF dataset not loaded");
        return getElevationFromTopo(lat, lon);
    }

    // Check if point is within interpolatable range (using base class method)
    if (!isPointInterpolatable(lat, lon))
    {
        logger.log("[WARN] Point (" + std::to_string(lat) + ", " + std::to_string(lon) + 
                  ") is outside interpolation range");
        return getElevationFromTopo(lat, lon);
    }

    // Convert lat/lon to pixel coordinates
    double pixel_x, pixel_y;
    if (!latLonToPixel(lat, lon, pixel_x, pixel_y))
    {
        logger.log("[ERROR] Failed to convert lat/lon to pixel coordinates");
        return getElevationFromTopo(lat, lon);
    }

    double elevation = getElevationAtPixel(pixel_x, pixel_y);

    // the geotiff has elevation at 0 for the not canyon surroundings
    // use the topo server to get that data
    if (elevation == 0)
    {
        elevation = getElevationFromTopo(lat, lon);
    }
    
    logger.log("[INFO] Elevation at (" + std::to_string(lat) + ", " + std::to_string(lon) + 
              "): " + std::to_string(elevation) + " m");
    
    return elevation;
}

/**
 * @brief Get a grid of elevation values for a given latitude and longitude range
 * @param latStart Starting latitude in degrees
 * @param lonStart Starting longitude in degrees
 * @param latEnd Ending latitude in degrees
 * @param lonEnd Ending longitude in degrees
 * @param numLatSamples Number of latitude samples
 * @param numLonSamples Number of longitude samples
 */
std::vector<TerrainPoint> TerrainServerGeoTIFF::getElevationGrid(double latStart,
    double lonStart, double latEnd, double lonEnd, int numLatSamples, int numLonSamples)
{
    logger.log("[INFO] [TerrainServerGeoTiff::getElevationGrid] Getting elevation grid from (" + std::to_string(latStart) + ", " +
              std::to_string(lonStart) + ") to (" + std::to_string(latEnd) + ", " +
              std::to_string(lonEnd) + ") with " + std::to_string(numLatSamples) +
              "x" + std::to_string(numLonSamples) + " samples.");

    std::vector<TerrainPoint> points;
    points.reserve(numLatSamples * numLonSamples);

    if (dataset_ == nullptr)
    {
        logger.log("[ERROR] GeoTIFF dataset not loaded, using OpenTopo");
        
        // Use OpenTopo's grid method which batches requests AND handles caching
        if (pServerTopo_ != nullptr)
        {
            return pServerTopo_->getElevationGrid(latStart, lonStart, latEnd, lonEnd, 
                                                  numLatSamples, numLonSamples);
        }
        return points;
    }

    double lat_step = (latEnd - latStart) / (numLatSamples - 1);
    double lon_step = (lonEnd - lonStart) / (numLonSamples - 1);

    // First pass: get all GeoTIFF data and identify missing points
    std::vector<std::pair<int, std::pair<double, double>>> missing_points;
    
    for (int i = 0; i < numLatSamples; ++i)
    {
        for (int j = 0; j < numLonSamples; ++j)
        {
            double lat = latStart + i * lat_step;
            double lon = lonStart + j * lon_step;

            TerrainPoint pt;
            pt.lat = lat;
            pt.lon = lon;
            
            // Try to get from GeoTIFF first
            bool needsTopoData = false;
            
            if (isPointInterpolatable(lat, lon))
            {
                double pixel_x, pixel_y;
                if (latLonToPixel(lat, lon, pixel_x, pixel_y))
                {
                    double elevation = getElevationAtPixel(pixel_x, pixel_y);
                    
                    // Check if we got valid data
                    if (has_nodata_ && std::abs(elevation - nodata_value_) < 1e-6)
                    {
                        needsTopoData = true;
                    }
                    else if (elevation == 0.0)  // Your specific case
                    {
                        needsTopoData = true;
                    }
                    else
                    {
                        pt.alt = elevation;
                    }
                }
                else
                {
                    needsTopoData = true;
                }
            }
            else
            {
                needsTopoData = true;
            }
            
            if (needsTopoData)
            {
                int index = i * numLonSamples + j;
                missing_points.push_back({index, {lat, lon}});
                pt.alt = 0.0;  // Placeholder
            }
            
            points.push_back(pt);
        }
    }

    // Second pass: batch fetch missing data from OpenTopo
    if (!missing_points.empty() && pServerTopo_ != nullptr)
    {
        logger.log("[INFO] [TerrainServerGeoTiff::getElevationGrid] Fetching " + std::to_string(missing_points.size()) + 
                  " missing points from OpenTopo");
        
        // Calculate bounding box for missing points
        double min_lat = missing_points[0].second.first;
        double max_lat = missing_points[0].second.first;
        double min_lon = missing_points[0].second.second;
        double max_lon = missing_points[0].second.second;
        
        for (const auto& mp : missing_points)
        {
            min_lat = std::min(min_lat, mp.second.first);
            max_lat = std::max(max_lat, mp.second.first);
            min_lon = std::min(min_lon, mp.second.second);
            max_lon = std::max(max_lon, mp.second.second);
        }
        
        // Use OpenTopo's getElevationGrid which handles caching and batching
        // This will check cache first, then fetch only what's needed
        std::vector<TerrainPoint> topo_grid = pServerTopo_->getElevationGrid(
            min_lat, min_lon, max_lat, max_lon,
            numLatSamples, numLonSamples);
        
        logger.log("[INFO] [TerrainServerGeoTiff::getElevationGrid] Fetched missing points from OpenTopo");

        // Map the topo_grid data to our missing points
        // Since we're using the same grid resolution, we can match by coordinates
        for (const auto& mp : missing_points)
        {
            int index = mp.first;
            double target_lat = mp.second.first;
            double target_lon = mp.second.second;
            
            // Find matching point in topo_grid
            for (const auto& topo_pt : topo_grid)
            {
                if (std::abs(topo_pt.lat - target_lat) < 1e-6 &&
                    std::abs(topo_pt.lon - target_lon) < 1e-6)
                {
                    points[index].alt = topo_pt.alt;
                    break;
                }
            }
        }
        
        logger.log("[INFO] [TerrainServerGeoTiff::getElevationGrid] Successfully filled missing points from OpenTopo");
    }

    logger.log("[INFO] [TerrainServerGeoTiff::getElevationGrid] Generated " + std::to_string(points.size()) + " elevation points");
    return points;
}

/**
 * @brief Set the Topo server for fallback elevation data
 * @param server Unique pointer to TerrainServerOpenTopo instance
 * @return true if server is set successfully, false otherwise
 */
bool TerrainServerGeoTIFF::setTopoServer(std::unique_ptr<TerrainServerOpenTopo> server)
{
    if (!server) return false; 

    pServerTopo_ = std::move(server);  // transfer ownership
    return true;
}