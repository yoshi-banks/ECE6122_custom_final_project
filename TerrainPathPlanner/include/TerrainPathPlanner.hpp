/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file TerrainPathPlanner.hpp
 * @brief TerrainPathPlanner class declaration with optional MPI support
 */

#ifndef TERRAIN_PATH_PLANNER_HPP
#define TERRAIN_PATH_PLANNER_HPP

#include <vector>
#include <memory>
#include <functional>
#include <cmath>
#include <chrono>

#include "NEDPoint.hpp"

#ifdef USE_MPI
#include <mpi.h>
#endif

/**
 * @brief A* path planning node
 */
struct PathNode {
    NEDPoint point;
    std::shared_ptr<PathNode> parent;
    double g_cost;  // Cost from start
    double h_cost;  // Heuristic cost to goal
    double f_cost;  // Total cost (g + h)
    
    PathNode(const NEDPoint& p) : point(p), parent(nullptr), g_cost(0), h_cost(0), f_cost(0) {}
};

/**
 * @brief Terrain-aware path planner
 */
class TerrainPathPlanner {
public:
    /**
     * @brief Configuration parameters
     */
    struct Config {
        double grid_resolution = 10.0;          // Grid spacing (meters)
        double altitude_step = 10.0;            // Altitude discretization (meters)
        double min_altitude_agl = 20.0;         // Minimum altitude AGL (meters)
        double max_altitude_agl = 100.0;        // Maximum altitude AGL (meters)
        double vertical_safety_margin = 5.0;    // Vertical clearance (meters)
        double horizontal_safety_radius = 10.0; // Horizontal clearance radius (meters)
        double max_terrain_gradient = 0.7;      // Max slope (rise/run)
        
        bool use_diagonal_moves = true;         // Allow 8-directional movement
        int max_iterations = 50000;             // A* iteration limit
        
        // Cost function weights
        double distance_weight = 1.0;
        double altitude_change_penalty = 0.5;
        double terrain_clearance_weight = 0.2;
    };
    
    /**
     * @brief Path statistics
     */
    struct PathStats {
        double total_distance = 0.0;
        double total_altitude_change = 0.0;
        double avg_agl = 0.0;
        double max_agl = 0.0;
        int num_waypoints = 0;
        double computation_time = 0.0;
    };
    
    /**
     * @brief Constructor
     */
    TerrainPathPlanner(const std::vector<std::vector<double>>& terrain_heights,
                      double north_min, double north_max,
                      double east_min, double east_max,
                      const Config& config);
    
    /**
     * @brief Compute path using A* (single-threaded)
     */
    std::vector<NEDPoint> computePath(const NEDPoint& start, const NEDPoint& goal);
    
#ifdef USE_MPI
    /**
     * @brief Compute path using MPI-parallelized A*
     * @note Only available when compiled with USE_MPI
     */
    std::vector<NEDPoint> computePathMPI(const NEDPoint& start, const NEDPoint& goal);
    
    /**
     * @brief Compute path using parallel RRT
     * @note Only available when compiled with USE_MPI
     */
    std::vector<NEDPoint> computePathRRTParallel(const NEDPoint& start, 
                                                 const NEDPoint& goal, 
                                                 int max_iterations);
    
    // Allow MPI implementation to access private members
    friend class TerrainPathPlannerMPI;
#endif
    
    /**
     * @brief Compute path using RRT (single-threaded)
     */
    std::vector<NEDPoint> computePathRRT(const NEDPoint& start, 
                                        const NEDPoint& goal, 
                                        int max_iterations);
    
    /**
     * @brief Get statistics of last computed path
     */
    const PathStats& getLastPathStats() const;
    
    /**
     * @brief Check if a point is valid
     */
    bool isValidPoint(const NEDPoint& point) const;
    
    /**
     * @brief Debug function for validating points
     */
    void debugValidPoint(const NEDPoint& point, const std::string& label) const;
    
private:
    // Terrain data
    std::vector<std::vector<double>> terrain_heights_;
    double north_min_, north_max_, north_range_;
    double east_min_, east_max_, east_range_;
    int grid_rows_, grid_cols_;
    
    // Configuration
    Config config_;
    
    // Statistics
    PathStats last_stats_;
    
    // Helper methods
    double getTerrainHeight(double north, double east) const;
    double heuristic(const NEDPoint& a, const NEDPoint& b) const;
    double computeCost(const NEDPoint& from, const NEDPoint& to) const;
    bool checkHorizontalClearance(const NEDPoint& point) const;
    std::vector<NEDPoint> getNeighbors(const NEDPoint& point, const NEDPoint& goal) const;
    std::vector<NEDPoint> reconstructPath(std::shared_ptr<PathNode> node);
    std::vector<NEDPoint> smoothPath(const std::vector<NEDPoint>& raw_path);
    void calculatePathStats(const std::vector<NEDPoint>& path);
    
public:
    // Public accessors for MPI implementation
    double getGridResolution() const { return config_.grid_resolution; }
    int getMaxIterations() const { return config_.max_iterations; }
};

#endif // TERRAIN_PATH_PLANNER_HPP