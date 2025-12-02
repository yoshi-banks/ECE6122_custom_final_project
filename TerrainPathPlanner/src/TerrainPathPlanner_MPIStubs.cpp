/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file TerrainPathPlanner_MPIStubs.cpp
 * @brief MPI method stubs for TerrainPathPlanner
 * 
 * @details This file provides the implementation of MPI methods in the
 * TerrainPathPlanner class. These methods delegate to TerrainPathPlannerMPI.
 * 
 * Add this to your TerrainPathPlanner.cpp or compile separately.
 */

#ifdef USE_MPI

#include "TerrainPathPlanner.hpp"
#include "TerrainPathPlannerMPI.hpp"

/**
 * @brief MPI-parallelized A* implementation (delegates to TerrainPathPlannerMPI)
 */
std::vector<NEDPoint> TerrainPathPlanner::computePathMPI(
    const NEDPoint& start, 
    const NEDPoint& goal)
{
    return TerrainPathPlannerMPI::computePathMPI(*this, start, goal);
}

/**
 * @brief Parallel RRT implementation (delegates to TerrainPathPlannerMPI)
 */
std::vector<NEDPoint> TerrainPathPlanner::computePathRRTParallel(
    const NEDPoint& start,
    const NEDPoint& goal,
    int max_iterations)
{
    return TerrainPathPlannerMPI::computePathRRTParallel(*this, start, goal, max_iterations);
}

#endif // USE_MPI