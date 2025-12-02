/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file TerrainServerGeoTIFF.hpp
 * @brief TerrainServerGeoTIFF class declaration
 */

#pragma once

#include <string>
#include <vector>
#include <gdal/gdal_priv.h>
#include <gdal/ogr_spatialref.h>

#include "ITerrainServer.hpp"
#include "TerrainServerOpenTopo.hpp"

/**
 * @brief TerrainServerGeoTIFF class
 * @details This class represents a terrain data server that reads from a GeoTIFF (Geographic Tiff) file.
 */
class TerrainServerGeoTIFF : public ITerrainServer
{
    public:
        // Constructor accepts path to GeoTIFF DEM file
        explicit TerrainServerGeoTIFF(const std::string& geotiff_path,
                                      double interpolation_tolerance = 0.01);
        
        ~TerrainServerGeoTIFF() override;

        // ITerrainServer interface implementation
        double getElevation(double lat, double lon) override;
        std::vector<TerrainPoint> getElevationGrid(double latStart,
            double lonStart, double latEnd, double lonEnd,
            int numLatSamples, int numLonSamples) override;

        // Additional utility functions
        bool isLoaded() const { return dataset_ != nullptr; }

        bool setTopoServer(std::unique_ptr<TerrainServerOpenTopo> server);

    private:
        std::string geotiff_path_;
        GDALDataset* dataset_;
        GDALRasterBand* band_;
        double geoTransform_[6];
        int raster_width_;
        int raster_height_;
        double nodata_value_;
        bool has_nodata_;
        double max_elevation_;
        std::unique_ptr<TerrainServerOpenTopo> pServerTopo_;

        // Internal helper functions
        bool loadGeoTIFF(const std::string& path);
        void calculateBoundingBox();
        void calculateMaxElevation();
        bool latLonToPixel(double lat, double lon, double& pixel_x, double& pixel_y) const;
        bool pixelToLatLon(double pixel_x, double pixel_y, double& lat, double& lon) const;
        double getElevationAtPixel(double pixel_x, double pixel_y);
        double bilinearInterpolation(double pixel_x, double pixel_y);

        double getElevationFromTopo(double lat, double lon);
};