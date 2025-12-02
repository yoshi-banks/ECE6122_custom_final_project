/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file TerrainPathPlanner.cpp
 * @brief TerrainPathPlanner class implementation
 */

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <cmath>

#include "TerrainPathPlanner.hpp"
#include "NEDPointHash.hpp"

/**
 * @brief TerrainPathPlanner constructor
 * @details Initializes the TerrainPathPlanner with terrain data and configuration.
 * @param terrain_heights 2D grid of terrain heights (meters)
 * @param north_min Minimum north coordinate (meters)
 * @param north_max Maximum north coordinate (meters)
 * @param east_min Minimum east coordinate (meters)
 * @param east_max Maximum east coordinate (meters)
 * @param config Configuration parameters
 */
TerrainPathPlanner::TerrainPathPlanner(const std::vector<std::vector<double>>& terrain_heights,
                               double north_min, double north_max,
                               double east_min, double east_max,
                               const Config& config)
    : terrain_heights_(terrain_heights)
    , north_min_(north_min)
    , north_max_(north_max)
    , east_min_(east_min)
    , east_max_(east_max)
    , config_(config)
{
    north_range_ = north_max_ - north_min_;
    east_range_ = east_max_ - east_min_;
    grid_rows_ = terrain_heights_.size();
    grid_cols_ = terrain_heights_.empty() ? 0 : terrain_heights_[0].size();
    
    std::cout << "[PathPlanner] Initialized with " << grid_rows_ << "x" << grid_cols_ 
              << " terrain grid\n";
}

/**
 * @brief Get the terrain height at a given NED point
 * @param north North coordinate (meters)
 * @param east East coordinate (meters)
 * @return Terrain height (meters)
 */
double TerrainPathPlanner::getTerrainHeight(double north, double east) const {
    // Convert NED to grid indices
    double norm_n = (north - north_min_) / north_range_;
    double norm_e = (east - east_min_) / east_range_;
    
    if (norm_n < 0 || norm_n > 1 || norm_e < 0 || norm_e > 1) {
        return 0; // Out of bounds
    }
    
    // Bilinear interpolation
    double row_f = norm_n * (grid_rows_ - 1);
    double col_f = norm_e * (grid_cols_ - 1);
    
    int row = static_cast<int>(row_f);
    int col = static_cast<int>(col_f);
    
    row = std::max(0, std::min(grid_rows_ - 2, row));
    col = std::max(0, std::min(grid_cols_ - 2, col));
    
    double dr = row_f - row;
    double dc = col_f - col;
    
    double h00 = terrain_heights_[row][col];
    double h10 = terrain_heights_[row + 1][col];
    double h01 = terrain_heights_[row][col + 1];
    double h11 = terrain_heights_[row + 1][col + 1];
    
    double h0 = h00 * (1 - dr) + h10 * dr;
    double h1 = h01 * (1 - dr) + h11 * dr;
    
    return h0 * (1 - dc) + h1 * dc;
}

/**
 * @brief Check if a NED point is valid
 * @param point NED point to check
 * @return true if valid, false otherwise
 */
bool TerrainPathPlanner::isValidPoint(const NEDPoint& point) const {
    // Check bounds
    if (point.north < north_min_ || point.north > north_max_ ||
        point.east < east_min_ || point.east > east_max_) {
        return false;
    }
    
    // Get terrain height (down is positive downward)
    double terrain_down = getTerrainHeight(point.north, point.east);
    
    // Calculate altitude above ground (negative because down is positive downward)
    double agl = terrain_down - point.down;

    // VERTICAL SAFETY: Must be above terrain + vertical margin
    if (agl < config_.vertical_safety_margin) {
        return false;  // Too close to terrain
    }
    
    // Check altitude constraints
    if (agl < config_.min_altitude_agl || agl > config_.max_altitude_agl) {
        return false;
    }
    
    // HORIZONTAL SAFETY: Check terrain gradient around the point
    if (!checkHorizontalClearance(point)) {
        return false;  // Too close to steep terrain
    }
    
    // Check altitude constraints
    return (agl >= config_.min_altitude_agl && agl <= config_.max_altitude_agl);
}

/**
 * @brief Debug function to print information about a valid NED point
 * @param point NED point to debug
 * @param label Label for the point (e.g., "START", "GOAL")
 * @details This function is used for debugging and provides information about a valid NED point.
 * It prints the point coordinates, terrain height, and altitude above ground.
 * If the point is invalid, it provides reasons why.
 */
void TerrainPathPlanner::debugValidPoint(const NEDPoint& point, const std::string& label) const {
    std::cout << "\n[DEBUG] Validating " << label << ":\n";
    std::cout << "  Point: (" << point.north << ", " << point.east << ", " << point.down << ")\n";
    
    // Check bounds
    bool in_bounds = (point.north >= north_min_ && point.north <= north_max_ &&
                      point.east >= east_min_ && point.east <= east_max_);
    std::cout << "  In bounds: " << (in_bounds ? "YES" : "NO") << "\n";
    
    if (!in_bounds) return;
    
    // Check terrain clearance
    double terrain_down = getTerrainHeight(point.north, point.east);
    double agl = terrain_down - point.down;
    std::cout << "  Terrain down: " << terrain_down << "\n";
    std::cout << "  Point down: " << point.down << "\n";
    std::cout << "  AGL: " << agl << " (min: " << config_.min_altitude_agl 
              << ", max: " << config_.max_altitude_agl << ")\n";
    std::cout << "  Vertical margin: " << config_.vertical_safety_margin 
              << " (need AGL >= " << config_.vertical_safety_margin << ")\n";
    
    if (agl < config_.vertical_safety_margin) {
        std::cout << "  FAILED: AGL too low for vertical safety margin\n";
        return;
    }
    
    if (agl < config_.min_altitude_agl || agl > config_.max_altitude_agl) {
        std::cout << "  FAILED: AGL outside bounds\n";
        return;
    }
    
    // Check horizontal clearance
    bool horiz_clear = checkHorizontalClearance(point);
    std::cout << "  Horizontal clearance: " << (horiz_clear ? "PASS" : "FAIL") << "\n";
    
    if (isValidPoint(point)) {
        std::cout << "  OVERALL: VALID\n";
    } else {
        std::cout << "  OVERALL: INVALID\n";
    }
}

/**
 * @brief Heuristic function for A* (Euclidean distance)
 * @param a First NED point
 * @param b Second NED point
 * @return Heuristic cost
 */
double TerrainPathPlanner::heuristic(const NEDPoint& a, const NEDPoint& b) const {
    // Euclidean distance (admissible)
    return a.distanceTo(b);
}

/**
 * @brief Compute cost from one point to another
 * @param from Starting NED point
 * @param to Target NED point
 * @return Cost value
 */
double TerrainPathPlanner::computeCost(const NEDPoint& from, const NEDPoint& to) const {
    double distance = from.distanceTo(to);
    double distance_cost = distance * config_.distance_weight;
    
    // Altitude change penalty
    double alt_change = std::abs(to.down - from.down);
    double altitude_cost = alt_change * config_.altitude_change_penalty;
    
    // Terrain clearance (prefer lower altitudes)
    double terrain_down = getTerrainHeight(to.north, to.east);
    double agl = terrain_down - to.down;
    double clearance_cost = (agl / config_.max_altitude_agl) * 
                           config_.terrain_clearance_weight * distance;
    
    return distance_cost + altitude_cost + clearance_cost;
}

/**
 * @brief Check horizontal clearance around a NED point
 * @param point NED point to check
 * @return true if clearance is sufficient, false otherwise
 */
bool TerrainPathPlanner::checkHorizontalClearance(const NEDPoint& point) const {
    // Sample terrain heights in a circle around the point
    double center_terrain = getTerrainHeight(point.north, point.east);
    double center_agl = center_terrain - point.down;
    
    // Check 8 directions around the point
    const int num_samples = 8;
    for (int i = 0; i < num_samples; ++i) {
        double angle = 2.0 * M_PI * i / num_samples;
        double sample_north = point.north + config_.horizontal_safety_radius * cos(angle);
        double sample_east = point.east + config_.horizontal_safety_radius * sin(angle);
        
        // Check if sample point is in bounds
        if (sample_north < north_min_ || sample_north > north_max_ ||
            sample_east < east_min_ || sample_east > east_max_) {
            continue;  // Skip out-of-bounds samples
        }
        
        double sample_terrain = getTerrainHeight(sample_north, sample_east);
        
        // Check if we're too close to higher terrain (cliff/steep slope)
        double terrain_diff = center_terrain - sample_terrain;  // Positive if sample is higher
        
        // If nearby terrain is significantly higher, we might be too close to a cliff
        if (terrain_diff > config_.vertical_safety_margin) {
            return false;  // Nearby terrain rises above our safety margin
        }
        
        // Check gradient (slope steepness)
        double gradient = std::abs(terrain_diff) / config_.horizontal_safety_radius;
        if (gradient > config_.max_terrain_gradient) {
            return false;  // Terrain is too steep nearby
        }
        
        // Also check if the flight path at this radius would be below terrain
        double sample_agl = sample_terrain - point.down;
        if (sample_agl < config_.vertical_safety_margin) {
            return false;  // Would be too close to terrain at this radius
        }
    }
    
    return true;
}

/**
 * @brief Get neighbors of a NED point
 * @param point NED point to get neighbors for
 * @param goal Target NED point
 * @return Vector of neighbors
 */
std::vector<NEDPoint> TerrainPathPlanner::getNeighbors(const NEDPoint& point, 
                                                    const NEDPoint& goal) const {
    std::vector<NEDPoint> neighbors;
    
    // Horizontal directions
    std::vector<std::pair<int, int>> directions = {
        {1, 0}, {-1, 0}, {0, 1}, {0, -1}
    };
    
    if (config_.use_diagonal_moves) {
        directions.push_back({1, 1});
        directions.push_back({1, -1});
        directions.push_back({-1, 1});
        directions.push_back({-1, -1});
    }
    
    for (const auto& [dn, de] : directions) {
        double new_north = point.north + dn * config_.grid_resolution;
        double new_east = point.east + de * config_.grid_resolution;
        
        // Get terrain height at new position
        double terrain_down = getTerrainHeight(new_north, new_east);
        
        // SIMPLIFIED: Try only 3 altitude candidates instead of 9+
        std::vector<double> altitude_candidates;
        
        // 1. Keep current altitude (terrain following)
        altitude_candidates.push_back(point.down);
        
        // 2. Optimal altitude (middle of AGL range)
        double optimal_agl = (config_.min_altitude_agl + config_.max_altitude_agl) * 0.5;
        altitude_candidates.push_back(terrain_down - optimal_agl);
        
        // 3. Goal altitude if nearby (helps reach goal)
        double dist_to_goal = std::sqrt((new_north - goal.north) * (new_north - goal.north) +
                                       (new_east - goal.east) * (new_east - goal.east));
        if (dist_to_goal < config_.grid_resolution * 5) {
            altitude_candidates.push_back(goal.down);
        }
        
        for (double down : altitude_candidates) {
            NEDPoint neighbor(new_north, new_east, down);
            if (isValidPoint(neighbor)) {
                neighbors.push_back(neighbor);
            }
        }
    }
    
    return neighbors;
}

/**
 * @brief Reconstruct path from goal to start
 * @param node Goal PathNode
 * @return Vector of NEDPoints representing the path
 */
std::vector<NEDPoint> TerrainPathPlanner::reconstructPath(std::shared_ptr<PathNode> node) {
    std::vector<NEDPoint> path;
    while (node != nullptr) {
        path.push_back(node->point);
        node = node->parent;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

/**
 * @brief Compute path from start to goal using A* algorithm
 * @param start Starting NED point
 * @param goal Goal NED point
 * @return Vector of NEDPoints representing the computed path
 */
std::vector<NEDPoint> TerrainPathPlanner::computePath(const NEDPoint& start, const NEDPoint& goal) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "[PathPlanner] Computing A* path...\n";
    
    if (!isValidPoint(start)) {
        std::cerr << "[PathPlanner] ERROR: Start point invalid\n";
        return {};
    }
    
    if (!isValidPoint(goal)) {
        std::cerr << "[PathPlanner] ERROR: Goal point invalid\n";
        return {};
    }
    
    std::priority_queue<std::shared_ptr<PathNode>,
                       std::vector<std::shared_ptr<PathNode>>,
                       std::function<bool(std::shared_ptr<PathNode>, std::shared_ptr<PathNode>)>>
        open_set([](const auto& a, const auto& b) { return a->f_cost > b->f_cost; });
    
    std::unordered_map<NEDPoint, double, NEDPointHash> g_scores;
    std::unordered_set<NEDPoint, NEDPointHash> closed_set;
    
    auto start_node = std::make_shared<PathNode>(start);
    start_node->g_cost = 0;
    start_node->h_cost = heuristic(start, goal);
    start_node->f_cost = start_node->h_cost;
    
    open_set.push(start_node);
    g_scores[start] = 0;
    
    int iterations = 0;
    
    while (!open_set.empty() && iterations < config_.max_iterations) {
        iterations++;
        
        auto current = open_set.top();
        open_set.pop();
        
        // Check if reached goal
        double dist_to_goal = std::sqrt(
            (current->point.north - goal.north) * (current->point.north - goal.north) +
            (current->point.east - goal.east) * (current->point.east - goal.east));
        
        if (dist_to_goal < config_.grid_resolution * 2) {
            std::cout << "[PathPlanner] Path found in " << iterations << " iterations\n";
            
            auto raw_path = reconstructPath(current);
            // auto smoothed_path = smoothPath(raw_path);
            auto smoothed_path = raw_path;
            
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            
            calculatePathStats(smoothed_path);
            last_stats_.computation_time = elapsed.count();
            
            std::cout << "[PathPlanner] Computed in " << elapsed.count() << "s, "
                     << smoothed_path.size() << " waypoints, "
                     << last_stats_.total_distance << "m distance\n";
            
            return smoothed_path;
        }
        
        if (closed_set.count(current->point) > 0) {
            continue;
        }
        
        closed_set.insert(current->point);
        
        auto neighbors = getNeighbors(current->point, goal);
        
        for (const auto& neighbor_point : neighbors) {
            if (closed_set.count(neighbor_point) > 0) {
                continue;
            }
            
            double tentative_g = current->g_cost + computeCost(current->point, neighbor_point);
            
            if (g_scores.find(neighbor_point) == g_scores.end() ||
                tentative_g < g_scores[neighbor_point]) {
                auto neighbor_node = std::make_shared<PathNode>(neighbor_point);
                neighbor_node->parent = current;
                neighbor_node->g_cost = tentative_g;
                neighbor_node->h_cost = heuristic(neighbor_point, goal);
                neighbor_node->f_cost = neighbor_node->g_cost + neighbor_node->h_cost;
                
                g_scores[neighbor_point] = tentative_g;
                open_set.push(neighbor_node);
            }
        }
        
        if (iterations % 1000 == 0) {
            std::cout << "Iteration " << iterations << ", dist to goal: " 
                     << dist_to_goal << "m       \r" << std::flush;
        }
    }
    
    std::cerr << "[PathPlanner] Path not found after " << iterations << " iterations\n";
    return {};
}

/**
 * @brief Compute path using RRT algorithm
 * @param start Starting NED point
 * @param goal Goal NED point
 * @return Vector of NEDPoints representing the computed path
 */
std::vector<NEDPoint> TerrainPathPlanner::computePathRRT(const NEDPoint& start, const NEDPoint& goal, int max_iterations)
{
    // TODO
    return {};
}

/**
 * @brief Get statistics of last computed path
 * @return PathStats structure
 */
const TerrainPathPlanner::PathStats& TerrainPathPlanner::getLastPathStats() const
{
    return last_stats_;
}

/**
 * @brief Smooth the computed path
 * @param raw_path Raw path as vector of NEDPoints
 * @return Smoothed path as vector of NEDPoints
 */
std::vector<NEDPoint> TerrainPathPlanner::smoothPath(const std::vector<NEDPoint>& raw_path) {
    if (raw_path.size() <= 2) {
        return raw_path;
    }
    
    std::vector<NEDPoint> smoothed;
    smoothed.push_back(raw_path[0]);
    
    size_t i = 0;
    while (i < raw_path.size() - 1) {
        size_t farthest = i + 1;
        
        for (size_t j = i + 2; j < raw_path.size(); ++j) {
            bool can_skip = true;
            
            int num_samples = 20;
            for (int k = 1; k < num_samples; ++k) {
                double t = static_cast<double>(k) / num_samples;
                NEDPoint sample(
                    raw_path[i].north + t * (raw_path[j].north - raw_path[i].north),
                    raw_path[i].east + t * (raw_path[j].east - raw_path[i].east),
                    raw_path[i].down + t * (raw_path[j].down - raw_path[i].down)
                );
                
                if (!isValidPoint(sample)) {
                    can_skip = false;
                    break;
                }
            }
            
            if (can_skip) {
                farthest = j;
            } else {
                break;
            }
        }
        
        smoothed.push_back(raw_path[farthest]);
        i = farthest;
    }
    
    return smoothed;
}

/**
 * @brief Calculate statistics of a path
 * @param path Path as vector of NEDPoints
 */
void TerrainPathPlanner::calculatePathStats(const std::vector<NEDPoint>& path) {
    last_stats_ = PathStats();
    
    if (path.empty()) return;
    
    last_stats_.num_waypoints = path.size();
    
    double total_dist = 0;
    double total_alt_change = 0;
    double total_agl = 0;
    double max_agl = 0;
    
    for (size_t i = 0; i < path.size(); ++i) {
        double terrain_down = getTerrainHeight(path[i].north, path[i].east);
        double agl = terrain_down - path[i].down;
        
        total_agl += agl;
        max_agl = std::max(max_agl, agl);
        
        if (i > 0) {
            total_dist += path[i].distanceTo(path[i-1]);
            total_alt_change += std::abs(path[i].down - path[i-1].down);
        }
    }
    
    last_stats_.total_distance = total_dist;
    last_stats_.total_altitude_change = total_alt_change;
    last_stats_.avg_agl = total_agl / path.size();
    last_stats_.max_agl = max_agl;
}