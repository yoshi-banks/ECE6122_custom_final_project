/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file PathNode.hpp
 * @brief PathNode struct definition
 */

#pragma once

#include <memory>

#include "NEDPoint.hpp"

/**
 * @brief PathNode struct
 * @details This struct represents a node in the pathfinding algorithm,
 * containing NED coordinates, cost values, and a pointer to the parent node.
 */
struct PathNode {
    NEDPoint point;
    double g_cost;
    double h_cost;
    double f_cost;
    std::shared_ptr<PathNode> parent;
    
    PathNode(const NEDPoint& p) 
        : point(p), g_cost(0), h_cost(0), f_cost(0), parent(nullptr) {}
    
    bool operator>(const PathNode& other) const {
        return f_cost > other.f_cost;
    }
};