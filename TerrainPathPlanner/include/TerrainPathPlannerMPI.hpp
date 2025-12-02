/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file TerrainPathPlannerMPI.hpp
 * @brief MPI-specific extensions for TerrainPathPlanner
 * 
 * @note This file is only included when USE_MPI is defined
 * 
 * This header provides MPI-parallelized versions of path planning algorithms.
 * It extends the base TerrainPathPlanner class with distributed computing capabilities.
 */

#ifndef TERRAIN_PATH_PLANNER_MPI_HPP
#define TERRAIN_PATH_PLANNER_MPI_HPP

#ifdef USE_MPI

#include <mpi.h>
#include <vector>
#include <memory>
#include <chrono>
#include "TerrainPathPlanner.hpp"
#include "NEDPoint.hpp"

// MPI message tags
#define TAG_NEIGHBOR_REQUEST 100
#define TAG_NEIGHBOR_RESPONSE 101
#define TAG_TERMINATE 103
#define TAG_VALIDATION_REQUEST 104
#define TAG_VALIDATION_RESPONSE 105

/**
 * @brief MPI-specific extensions to TerrainPathPlanner
 * 
 * @details This class provides MPI-parallelized implementations of path planning
 * algorithms. It uses a master-worker pattern where:
 * - Rank 0 (master) manages the search algorithm and open set
 * - Ranks 1-N (workers) perform neighbor expansion and validation
 * 
 * All methods in this class are only available when compiled with USE_MPI.
 */
class TerrainPathPlannerMPI {
public:
    /**
     * @brief Compute path using MPI-parallelized A*
     * 
     * @param planner Reference to the base planner (contains terrain data and config)
     * @param start Starting NED point
     * @param goal Goal NED point
     * @return Vector of NEDPoints representing the computed path
     * 
     * @details Master process manages the priority queue and closed set.
     * Worker processes expand neighbors in parallel using round-robin assignment.
     * Automatically falls back to sequential A* if only 1 process is available.
     */
    static std::vector<NEDPoint> computePathMPI(
        TerrainPathPlanner& planner,
        const NEDPoint& start, 
        const NEDPoint& goal);
    
    /**
     * @brief Compute path using parallel RRT
     * 
     * @param planner Reference to the base planner
     * @param start Starting NED point
     * @param goal Goal NED point
     * @param max_iterations Maximum iterations per process
     * @return Vector of NEDPoints representing the best computed path
     * 
     * @details Each process builds an independent RRT tree with different random seed.
     * After all processes complete, the master selects and broadcasts the best path.
     */
    static std::vector<NEDPoint> computePathRRTParallel(
        TerrainPathPlanner& planner,
        const NEDPoint& start, 
        const NEDPoint& goal, 
        int max_iterations);

private:
    /**
     * @brief Master process for MPI A* search
     * 
     * @param planner Reference to the base planner
     * @param start Starting NED point
     * @param goal Goal NED point
     * @param num_processes Total number of MPI processes
     * @param start_time Timestamp when computation started
     * @return Vector of NEDPoints representing the computed path
     * 
     * @details Manages the A* algorithm:
     * - Maintains priority queue (open set)
     * - Tracks visited nodes (closed set)
     * - Distributes neighbor expansion to workers
     * - Reconstructs and smooths final path
     */
    static std::vector<NEDPoint> computePathMaster(
        TerrainPathPlanner& planner,
        const NEDPoint& start,
        const NEDPoint& goal,
        int num_processes,
        std::chrono::high_resolution_clock::time_point start_time);
    
    /**
     * @brief Worker process for MPI A* search
     * 
     * @param planner Reference to the base planner
     * @param goal Goal NED point (used for heuristic)
     * 
     * @details Worker main loop:
     * 1. Wait for neighbor request from master
     * 2. Compute valid neighbors for given point
     * 3. Send neighbors back to master
     * 4. Repeat until termination signal
     */
    static void computePathWorker(
        TerrainPathPlanner& planner,
        const NEDPoint& goal);
    
    /**
     * @brief MPI-parallelized path smoothing
     * 
     * @param planner Reference to the base planner
     * @param raw_path Raw path from A* or RRT
     * @param num_processes Total number of MPI processes
     * @return Smoothed path with fewer waypoints
     * 
     * @details Distributes path segment validation across workers:
     * - Master tries to skip waypoints by validating direct paths
     * - Workers validate line-of-sight between waypoints
     * - Falls back to sequential smoothing if only 1 process
     */
    static std::vector<NEDPoint> smoothPathMPI(
        TerrainPathPlanner& planner,
        const std::vector<NEDPoint>& raw_path, 
        int num_processes);
    
    /**
     * @brief Validate a path segment (used by workers)
     * 
     * @param planner Reference to the base planner
     * @param from Starting point of segment
     * @param to Ending point of segment
     * @param num_samples Number of points to check along segment
     * @return true if entire segment is valid, false otherwise
     * 
     * @details Checks multiple points along the line segment to ensure
     * the entire path from 'from' to 'to' is collision-free.
     */
    static bool validateSegment(
        TerrainPathPlanner& planner,
        const NEDPoint& from,
        const NEDPoint& to,
        int num_samples = 20);
    
    /**
     * @brief Send neighbor data to master
     * 
     * @param neighbors Vector of neighbor points
     * @param dest_rank Destination rank (typically 0 for master)
     */
    static void sendNeighbors(
        const std::vector<NEDPoint>& neighbors,
        int dest_rank);
    
    /**
     * @brief Receive neighbor data from worker
     * 
     * @param source_rank Source rank to receive from
     * @return Vector of neighbor points
     */
    static std::vector<NEDPoint> receiveNeighbors(int source_rank);
    
    /**
     * @brief Send path data (used for broadcasting best RRT result)
     * 
     * @param path Path to send
     * @param dest_rank Destination rank (-1 for broadcast)
     */
    static void sendPath(const std::vector<NEDPoint>& path, int dest_rank = -1);
    
    /**
     * @brief Receive path data
     * 
     * @param source_rank Source rank to receive from
     * @return Received path
     */
    static std::vector<NEDPoint> receivePath(int source_rank);
};

#endif // USE_MPI

#endif // TERRAIN_PATH_PLANNER_MPI_HPP