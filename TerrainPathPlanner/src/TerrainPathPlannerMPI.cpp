/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file TerrainPathPlannerMPI.cpp
 * @brief MPI-parallelized terrain path planner implementation
 * 
 * @note This file is only compiled when USE_MPI is defined
 */

#ifdef USE_MPI

#include "TerrainPathPlannerMPI.hpp"
#include "NEDPointHash.hpp"

#include <algorithm>
#include <iostream>
#include <random>
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <limits>
#include <ctime>

/**
 * @brief Compute path using MPI-parallelized A*
 */
std::vector<NEDPoint> TerrainPathPlannerMPI::computePathMPI(
    TerrainPathPlanner& planner,
    const NEDPoint& start, 
    const NEDPoint& goal) 
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    if (rank == 0) {
        std::cout << "[PathPlanner] Computing MPI A* path with " << size << " process(es)...\n";
    }
    
    // Validate start and goal on all processes
    if (!planner.isValidPoint(start)) {
        if (rank == 0) std::cerr << "[PathPlanner] ERROR: Start point invalid\n";
        return {};
    }
    
    if (!planner.isValidPoint(goal)) {
        if (rank == 0) std::cerr << "[PathPlanner] ERROR: Goal point invalid\n";
        return {};
    }
    
    // If only one process, fall back to sequential version
    if (size == 1) {
        if (rank == 0) {
            std::cout << "[PathPlanner] Single process detected, using sequential A*\n";
        }
        return planner.computePath(start, goal);
    }
    
    if (rank == 0) {
        // MASTER PROCESS: Manages A* search
        return computePathMaster(planner, start, goal, size, start_time);
    } else {
        // WORKER PROCESSES: Expand neighbors
        computePathWorker(planner, goal);
        return {};
    }
}

/**
 * @brief Master process for MPI A* search
 */
std::vector<NEDPoint> TerrainPathPlannerMPI::computePathMaster(
    TerrainPathPlanner& planner,
    const NEDPoint& start,
    const NEDPoint& goal,
    int num_processes,
    std::chrono::high_resolution_clock::time_point start_time)
{
    // Priority queue for open set (min-heap based on f_cost)
    std::priority_queue<std::shared_ptr<PathNode>,
                       std::vector<std::shared_ptr<PathNode>>,
                       std::function<bool(std::shared_ptr<PathNode>, std::shared_ptr<PathNode>)>>
        open_set([](const auto& a, const auto& b) { return a->f_cost > b->f_cost; });
    
    std::unordered_map<NEDPoint, double, NEDPointHash> g_scores;
    std::unordered_set<NEDPoint, NEDPointHash> closed_set;
    
    // Initialize start node
    auto start_node = std::make_shared<PathNode>(start);
    start_node->g_cost = 0;
    start_node->h_cost = start.distanceTo(goal);
    start_node->f_cost = start_node->h_cost;
    
    open_set.push(start_node);
    g_scores[start] = 0;
    
    int iterations = 0;
    int next_worker = 1;  // Round-robin worker assignment
    
    // Get config from planner (need to access private members through public interface)
    const int max_iterations = 300000; // Use planner's max_iterations
    const double grid_resolution = 30.0; // Use planner's grid_resolution
    
    while (!open_set.empty() && iterations < max_iterations) {
        iterations++;
        
        auto current = open_set.top();
        open_set.pop();
        
        // Check if reached goal
        double dist_to_goal = current->point.distanceTo(goal);
        
        if (dist_to_goal < grid_resolution * 2) {
            std::cout << "[PathPlanner] Path found in " << iterations << " iterations\n";
            
            // Signal all workers to terminate USING TAG_NEIGHBOR_REQUEST with special marker
            for (int i = 1; i < num_processes; ++i) {
                double terminate_marker[3] = {-999999.0, -999999.0, -999999.0};
                MPI_Send(terminate_marker, 3, MPI_DOUBLE, i, TAG_NEIGHBOR_REQUEST, MPI_COMM_WORLD);
            }
            
            // Reconstruct path
            std::vector<NEDPoint> raw_path;
            auto node = current;
            while (node != nullptr) {
                raw_path.push_back(node->point);
                node = node->parent;
            }
            std::reverse(raw_path.begin(), raw_path.end());
            
            // Smooth path
            auto smoothed_path = smoothPathMPI(planner, raw_path, num_processes);
            
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            
            std::cout << "[PathPlanner] Computed in " << elapsed.count() << "s, "
                     << smoothed_path.size() << " waypoints\n";
            
            return smoothed_path;
        }
        
        // Skip if already visited
        if (closed_set.count(current->point) > 0) {
            continue;
        }
        
        closed_set.insert(current->point);
        
        // Send current point to worker for neighbor expansion
        double point_data[3] = {current->point.north, current->point.east, current->point.down};
        MPI_Send(point_data, 3, MPI_DOUBLE, next_worker, TAG_NEIGHBOR_REQUEST, MPI_COMM_WORLD);
        
        // Receive neighbors from worker
        int num_neighbors;
        MPI_Status status;
        MPI_Recv(&num_neighbors, 1, MPI_INT, next_worker, TAG_NEIGHBOR_RESPONSE, 
                MPI_COMM_WORLD, &status);
        
        if (num_neighbors > 0) {
            std::vector<double> neighbor_data(num_neighbors * 3);
            MPI_Recv(neighbor_data.data(), num_neighbors * 3, MPI_DOUBLE, 
                    next_worker, TAG_NEIGHBOR_RESPONSE, MPI_COMM_WORLD, &status);
            
            // Process received neighbors
            for (int i = 0; i < num_neighbors; ++i) {
                NEDPoint neighbor_point(neighbor_data[i*3], neighbor_data[i*3+1], neighbor_data[i*3+2]);
                
                if (closed_set.count(neighbor_point) > 0) {
                    continue;
                }
                
                // Compute cost (simple Euclidean distance for now)
                double tentative_g = current->g_cost + current->point.distanceTo(neighbor_point);
                
                if (g_scores.find(neighbor_point) == g_scores.end() ||
                    tentative_g < g_scores[neighbor_point]) {
                    auto neighbor_node = std::make_shared<PathNode>(neighbor_point);
                    neighbor_node->parent = current;
                    neighbor_node->g_cost = tentative_g;
                    neighbor_node->h_cost = neighbor_point.distanceTo(goal);
                    neighbor_node->f_cost = neighbor_node->g_cost + neighbor_node->h_cost;
                    
                    g_scores[neighbor_point] = tentative_g;
                    open_set.push(neighbor_node);
                }
            }
        }
        
        // Round-robin worker assignment
        next_worker = (next_worker % (num_processes - 1)) + 1;
        
        if (iterations % 1000 == 0) {
            std::cout << "Iteration " << iterations << ", dist to goal: " 
                     << dist_to_goal << "m       \r" << std::flush;
        }
    }
    
    // Terminate workers - ALWAYS DO THIS
    for (int i = 1; i < num_processes; ++i) {
        double terminate_marker[3] = {-999999.0, -999999.0, -999999.0};
        MPI_Send(terminate_marker, 3, MPI_DOUBLE, i, TAG_NEIGHBOR_REQUEST, MPI_COMM_WORLD);
    }
    
    std::cerr << "[PathPlanner] Path not found after " << iterations << " iterations\n";
    return {};
}

/**
 * @brief Worker process for MPI A* search
 */
void TerrainPathPlannerMPI::computePathWorker(
    TerrainPathPlanner& planner,
    const NEDPoint& goal) 
{
    while (true) {
        // Wait for neighbor request
        double point_data[3];
        MPI_Status status;
        MPI_Recv(point_data, 3, MPI_DOUBLE, 0, TAG_NEIGHBOR_REQUEST, 
                MPI_COMM_WORLD, &status);
        
        // Check for termination marker
        if (point_data[0] == -999999.0 && point_data[1] == -999999.0 && point_data[2] == -999999.0) {
            break;  // Exit worker loop
        }
        
        NEDPoint current(point_data[0], point_data[1], point_data[2]);
        
        // Compute neighbors
        auto neighbors = planner.getNeighbors(current, goal);
        
        // Send results...
        int num_neighbors = neighbors.size();
        MPI_Send(&num_neighbors, 1, MPI_INT, 0, TAG_NEIGHBOR_RESPONSE, MPI_COMM_WORLD);
        
        if (num_neighbors > 0) {
            std::vector<double> neighbor_data;
            neighbor_data.reserve(num_neighbors * 3);
            for (const auto& n : neighbors) {
                neighbor_data.push_back(n.north);
                neighbor_data.push_back(n.east);
            }
            MPI_Send(neighbor_data.data(), num_neighbors * 3, MPI_DOUBLE, 
                    0, TAG_NEIGHBOR_RESPONSE, MPI_COMM_WORLD);
        }
    }
}

/**
 * @brief MPI-parallelized path smoothing
 */
std::vector<NEDPoint> TerrainPathPlannerMPI::smoothPathMPI(
    TerrainPathPlanner& planner,
    const std::vector<NEDPoint>& raw_path,
    int num_processes)
{
    if (raw_path.size() <= 2 || num_processes <= 1) {
        return planner.smoothPath(raw_path);  // Fall back to sequential
    }
    
    std::vector<NEDPoint> smoothed;
    smoothed.push_back(raw_path[0]);
    
    size_t i = 0;
    int worker_idx = 0;
    
    while (i < raw_path.size() - 1) {
        size_t farthest = i + 1;
        
        // Try to skip ahead as far as possible
        for (size_t j = i + 2; j < raw_path.size(); ++j) {
            // Use workers to validate segments in parallel
            int worker = (worker_idx % (num_processes - 1)) + 1;
            worker_idx++;
            
            // Send segment to worker for validation
            double segment[6] = {
                raw_path[i].north, raw_path[i].east, raw_path[i].down,
                raw_path[j].north, raw_path[j].east, raw_path[j].down
            };
            MPI_Send(segment, 6, MPI_DOUBLE, worker, TAG_VALIDATION_REQUEST, MPI_COMM_WORLD);
            
            // Receive validation result
            int valid;
            MPI_Status status;
            MPI_Recv(&valid, 1, MPI_INT, worker, TAG_VALIDATION_RESPONSE, MPI_COMM_WORLD, &status);
            
            if (valid) {
                farthest = j;
            } else {
                break;
            }
        }
        
        smoothed.push_back(raw_path[farthest]);
        i = farthest;
    }
    
    // Signal workers to stop validation
    for (int w = 1; w < num_processes; ++w) {
        double dummy[6] = {0};
        MPI_Send(dummy, 6, MPI_DOUBLE, w, TAG_TERMINATE, MPI_COMM_WORLD);
    }
    
    return smoothed;
}

/**
 * @brief Validate a path segment
 */
bool TerrainPathPlannerMPI::validateSegment(
    TerrainPathPlanner& planner,
    const NEDPoint& from,
    const NEDPoint& to,
    int num_samples)
{
    for (int k = 1; k < num_samples; ++k) {
        double t = static_cast<double>(k) / num_samples;
        NEDPoint sample(
            from.north + t * (to.north - from.north),
            from.east + t * (to.east - from.east),
            from.down + t * (to.down - from.down)
        );
        
        if (!planner.isValidPoint(sample)) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Parallel RRT implementation using MPI
 */
std::vector<NEDPoint> TerrainPathPlannerMPI::computePathRRTParallel(
    TerrainPathPlanner& planner,
    const NEDPoint& start,
    const NEDPoint& goal,
    int max_iterations)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == 0) {
        std::cout << "[PathPlanner] Computing parallel RRT with " << size << " process(es)...\n";
    }
    
    // Each process runs RRT with different random seed
    std::mt19937 rng(rank * 12345 + time(nullptr));
    
    // Run local RRT (sequential implementation)
    std::vector<NEDPoint> local_path = planner.computePathRRT(start, goal, max_iterations);
    
    // Calculate local path cost
    double local_cost = std::numeric_limits<double>::max();
    if (!local_path.empty()) {
        local_cost = 0.0;
        for (size_t i = 1; i < local_path.size(); ++i) {
            local_cost += local_path[i].distanceTo(local_path[i-1]);
        }
    }
    
    // Gather all path costs to master
    std::vector<double> all_costs(size);
    MPI_Gather(&local_cost, 1, MPI_DOUBLE, all_costs.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    int best_rank = 0;
    if (rank == 0) {
        // Master finds best path
        best_rank = std::distance(all_costs.begin(), 
                                 std::min_element(all_costs.begin(), all_costs.end()));
        std::cout << "[PathPlanner] Best path found by rank " << best_rank 
                  << " with cost " << all_costs[best_rank] << "\n";
    }
    
    // Broadcast best rank
    MPI_Bcast(&best_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Best process broadcasts its path size
    int path_size = local_path.size();
    MPI_Bcast(&path_size, 1, MPI_INT, best_rank, MPI_COMM_WORLD);
    
    // Other processes allocate memory
    if (rank != best_rank) {
        local_path.resize(path_size);
    }
    
    // Broadcast path data
    std::vector<double> path_data(path_size * 3);
    if (rank == best_rank) {
        for (size_t i = 0; i < local_path.size(); ++i) {
            path_data[i*3] = local_path[i].north;
            path_data[i*3+1] = local_path[i].east;
            path_data[i*3+2] = local_path[i].down;
        }
    }
    
    MPI_Bcast(path_data.data(), path_size * 3, MPI_DOUBLE, best_rank, MPI_COMM_WORLD);
    
    // Reconstruct path on non-best processes
    if (rank != best_rank) {
        local_path.clear();
        for (int i = 0; i < path_size; ++i) {
            local_path.emplace_back(path_data[i*3], path_data[i*3+1], path_data[i*3+2]);
        }
    }
    
    return local_path;
}

/**
 * @brief Send neighbor data to master
 */
void TerrainPathPlannerMPI::sendNeighbors(
    const std::vector<NEDPoint>& neighbors,
    int dest_rank)
{
    int num_neighbors = neighbors.size();
    MPI_Send(&num_neighbors, 1, MPI_INT, dest_rank, TAG_NEIGHBOR_RESPONSE, MPI_COMM_WORLD);
    
    if (num_neighbors > 0) {
        std::vector<double> data;
        data.reserve(num_neighbors * 3);
        for (const auto& n : neighbors) {
            data.push_back(n.north);
            data.push_back(n.east);
            data.push_back(n.down);
        }
        MPI_Send(data.data(), num_neighbors * 3, MPI_DOUBLE, 
                dest_rank, TAG_NEIGHBOR_RESPONSE, MPI_COMM_WORLD);
    }
}

/**
 * @brief Receive neighbor data from worker
 */
std::vector<NEDPoint> TerrainPathPlannerMPI::receiveNeighbors(int source_rank)
{
    int num_neighbors;
    MPI_Status status;
    MPI_Recv(&num_neighbors, 1, MPI_INT, source_rank, TAG_NEIGHBOR_RESPONSE, 
            MPI_COMM_WORLD, &status);
    
    std::vector<NEDPoint> neighbors;
    if (num_neighbors > 0) {
        std::vector<double> data(num_neighbors * 3);
        MPI_Recv(data.data(), num_neighbors * 3, MPI_DOUBLE, 
                source_rank, TAG_NEIGHBOR_RESPONSE, MPI_COMM_WORLD, &status);
        
        neighbors.reserve(num_neighbors);
        for (int i = 0; i < num_neighbors; ++i) {
            neighbors.emplace_back(data[i*3], data[i*3+1], data[i*3+2]);
        }
    }
    
    return neighbors;
}

/**
 * @brief Send path data
 */
void TerrainPathPlannerMPI::sendPath(
    const std::vector<NEDPoint>& path,
    int dest_rank)
{
    int path_size = path.size();
    
    if (dest_rank == -1) {
        // Broadcast
        MPI_Bcast(&path_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        MPI_Send(&path_size, 1, MPI_INT, dest_rank, 0, MPI_COMM_WORLD);
    }
    
    if (path_size > 0) {
        std::vector<double> data(path_size * 3);
        for (size_t i = 0; i < path.size(); ++i) {
            data[i*3] = path[i].north;
            data[i*3+1] = path[i].east;
            data[i*3+2] = path[i].down;
        }
        
        if (dest_rank == -1) {
            MPI_Bcast(data.data(), path_size * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else {
            MPI_Send(data.data(), path_size * 3, MPI_DOUBLE, dest_rank, 0, MPI_COMM_WORLD);
        }
    }
}

/**
 * @brief Receive path data
 */
std::vector<NEDPoint> TerrainPathPlannerMPI::receivePath(int source_rank)
{
    int path_size;
    MPI_Status status;
    
    if (source_rank == -1) {
        // Broadcast receive
        MPI_Bcast(&path_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        MPI_Recv(&path_size, 1, MPI_INT, source_rank, 0, MPI_COMM_WORLD, &status);
    }
    
    std::vector<NEDPoint> path;
    if (path_size > 0) {
        std::vector<double> data(path_size * 3);
        
        if (source_rank == -1) {
            MPI_Bcast(data.data(), path_size * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else {
            MPI_Recv(data.data(), path_size * 3, MPI_DOUBLE, source_rank, 0, 
                    MPI_COMM_WORLD, &status);
        }
        
        path.reserve(path_size);
        for (int i = 0; i < path_size; ++i) {
            path.emplace_back(data[i*3], data[i*3+1], data[i*3+2]);
        }
    }
    
    return path;
}

#endif // USE_MPI