/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file main.cpp
 * 
 * @brief Main entry point for the Terrain Visualizer application.
 * 
 * @details This application visualizes terrain data and simulates UAV flight paths
 * using OpenGL. It connects to a PX4-based UAV via MAVSDK to retrieve real-time telemetry data.
 * It also supports uploading and commanding missions for the UAV.
 */

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <iomanip>
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "TerrainServerOpenTopo.hpp"
#include "TerrainServerGeoTIFF.hpp"
#include "Camera.hpp"
#include "TerrainMesh.hpp"
#include "Shader.hpp"
#include "Cube.hpp"
#include "TerrainPathPlanner.hpp"
#include "MouseState.hpp"
#include "PX4Interface.hpp"
#include "NED.hpp"

#ifdef USE_MPI
bool is_mpi_initialized = false;
int mpi_rank = 0;
int mpi_size = 1;
#endif

// Global state for mouse input
MouseState mouseState;

Camera camera;

float deltaTime = 0.0f;
float lastFrame = 0.0f;

const int WINDOW_WIDTH = 1200;
const int WINDOW_HEIGHT = 800;

// Add global PX4 state (after other globals)
std::unique_ptr<PX4Interface> px4_interface;
std::vector<Cube> vehicleCubes;
std::mutex vehicle_mutex; // Protect vehicle cubes from concurrent access

NED px4_to_terrain_offset(0, 0, 0);  // Offset from PX4 home to terrain center

// Callbacks
/**
 * @brief Framebuffer size callback
 * 
 * @details Called when the window is resized. Updates the viewport to match the new size.
 * @param window The GLFW window.
 * @param width The new width of the window.
 * @param height The new height of the window.
 */
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

/**
 * @brief Mouse movement callback
 * 
 * @details Called when the mouse moves. Updates the camera's position based on mouse input.
 * @param window The GLFW window.
 * @param xpos The new x-position of the mouse.
 * @param ypos The new y-position of the mouse.
 */
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (mouseState.firstMouse)
    {
        mouseState.lastX = xpos;
        mouseState.lastY = ypos;
        mouseState.firstMouse = false;
    }

    double xoffset = xpos - mouseState.lastX;
    double yoffset = mouseState.lastY - ypos; // Reversed since y-coordinates go from bottom to top

    mouseState.lastX = xpos;
    mouseState.lastY = ypos;

    if (mouseState.leftButtonPressed)
    {
        camera.processMouseMovement(xoffset, yoffset);
    }
    else if (mouseState.rightButtonPressed)
    {
        camera.pan(-xoffset, yoffset);
    }
}

/**
 * @brief Mouse button callback
 * 
 * @details Called when a mouse button is pressed or released. Updates the mouse state accordingly.
 * @param window The GLFW window.   
 * @param button The button that was pressed or released.
 * @param action GLFW_PRESS or GLFW_RELEASE.
 * @param mods Any modifier keys that were held down when the button was pressed.
 */
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        mouseState.leftButtonPressed = (action == GLFW_PRESS);
    }
    else if (button == GLFW_MOUSE_BUTTON_RIGHT)
    {
        mouseState.rightButtonPressed = (action == GLFW_PRESS);
    }
    else if (button == GLFW_MOUSE_BUTTON_MIDDLE)
    {
        mouseState.middleButtonPressed = (action == GLFW_PRESS);
    }
}

/**
 * @brief Mouse scroll callback
 * 
 * @details Called when the mouse wheel is scrolled. Zooms the camera in or out.
 * @param window The GLFW window.
 * @param xoffset The scroll offset in the x-direction.
 * @param yoffset The scroll offset in the y-direction.
 */
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.processMouseScroll(yoffset);
}

/**
 * @brief Keyboard key callback
 * 
 * @details Called when a key is pressed or released. Handles escape key and speed adjustment.
 * @param window The GLFW window.
 * @param key The key that was pressed or released.
 * @param scancode The scancode of the key.
 * @param action GLFW_PRESS or GLFW_RELEASE.
 * @param mods Any modifier keys that were held down when the key was pressed.
 */
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, true);
    }

    // Speed adjustment with Shift
    if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_PRESS)
    {
        camera.setMovementSpeed(400.0f);  // Fast mode
    }
    if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_RELEASE)
    {
        camera.setMovementSpeed(100.0f);   // Normal mode
    }
    
    // Reset camera view
    if (key == GLFW_KEY_R && action == GLFW_PRESS)
    {
        camera.setDistance(50.0f);
    }
}

/**
 * @brief Process input from the keyboard
 * 
 * @details Checks for continuous key presses and updates the camera accordingly.
 * @param window The GLFW window.
 */
void processInput(GLFWwindow* window)
{
    // Movement keys
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_W, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_S, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_A, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_D, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_Q, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_E, deltaTime);
    
    // Rotation keys
    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_RIGHT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_UP, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_DOWN, deltaTime);
    
    // Zoom keys
    if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_Z, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
        camera.processKeyboard(GLFW_KEY_X, deltaTime);
}

/**
 * @brief Vehicle position update callback
 * 
 * @details Called when the vehicle position changes. Updates the vehicle position in the scene.
 * @param state The vehicle state.
 */
void onVehiclePositionUpdate(const PX4Interface::VehicleState& state)
{
    std::lock_guard<std::mutex> lock(vehicle_mutex);
    
    if (!vehicleCubes.empty()) {
        // Vehicle state is in PX4's NED frame (relative to start location)
        // Need to convert to terrain center frame for visualization
        NEDPoint vehicle_pos(
            state.north - px4_to_terrain_offset.north,  // Add offset to convert frames
            state.east - px4_to_terrain_offset.east,
            state.down - px4_to_terrain_offset.down
        );
        vehicleCubes[0].setPosition(vehicle_pos);
    }
}

/**
 * @brief Process PX4 input commands
 * 
 * @details Checks for specific key presses to control the PX4 vehicle.
 * @param window The GLFW window.
 */
void processPX4Input(GLFWwindow* window)
{
    if (!px4_interface || !px4_interface->isConnected()) return;
    
    // '1' key - Arm and takeoff
    static bool key1_pressed = false;
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS && !key1_pressed) {
        key1_pressed = true;
        std::cout << "\n[Control] Arming and taking off...\n";
        px4_interface->armAndTakeoff(10.0f);
    }
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_RELEASE) {
        key1_pressed = false;
    }
    
    // '2' key - Start mission
    static bool key2_pressed = false;
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS && !key2_pressed) {
        key2_pressed = true;
        std::cout << "\n[Control] Starting mission...\n";
        px4_interface->startMission();
    }
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_RELEASE) {
        key2_pressed = false;
    }
    
    // '3' key - Pause mission
    static bool key3_pressed = false;
    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS && !key3_pressed) {
        key3_pressed = true;
        std::cout << "\n[Control] Pausing mission...\n";
        px4_interface->pauseMission();
    }
    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_RELEASE) {
        key3_pressed = false;
    }
    
    // '4' key - Return to launch
    static bool key4_pressed = false;
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS && !key4_pressed) {
        key4_pressed = true;
        std::cout << "\n[Control] Returning to launch...\n";
        px4_interface->returnToLaunch();
    }
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_RELEASE) {
        key4_pressed = false;
    }
    
    // '5' key - Land
    static bool key5_pressed = false;
    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS && !key5_pressed) {
        key5_pressed = true;
        std::cout << "\n[Control] Landing...\n";
        px4_interface->land();
    }
    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_RELEASE) {
        key5_pressed = false;
    }
}

/**
 * @brief Print help instructions
 */
void printHelp()
{
    std::cout << "\n=== Terrain Visualizer Controls ===\n";
    std::cout << "MOUSE:\n";
    std::cout << "  Left Click + Drag:    Rotate camera\n";
    std::cout << "  Right/Middle + Drag:  Pan camera\n";
    std::cout << "  Mouse Wheel:          Zoom in/out\n";
    std::cout << "\nKEYBOARD:\n";
    std::cout << "  W/A/S/D:     Move forward/left/back/right\n";
    std::cout << "  Q/E:         Move down/up\n";
    std::cout << "  Arrow Keys:  Rotate camera\n";
    std::cout << "  Z/X:         Zoom in/out\n";
    std::cout << "  Left Shift:  Hold for faster movement\n";
    std::cout << "  R:           Reset camera distance\n";
    std::cout << "\nPX4 CONTROLS (Fixed-Wing):\n";
    std::cout << "  1:           Arm and takeoff (catapult launch)\n";
    std::cout << "  2:           Start mission\n";
    std::cout << "  3:           Pause mission (enter loiter)\n";
    std::cout << "  4:           Return to launch\n";
    std::cout << "  5:           Land (approach pattern)\n";
    std::cout << "  ESC:         Exit\n";
    std::cout << "\nNOTE: Fixed-wing aircraft cannot hover!\n";
    std::cout << "      Vehicle will circle when paused.\n";
    std::cout << "====================================\n\n";
}

/**
 * @brief Diagnose path for PX4
 * 
 * @details Diagnose the path to be used with PX4. Prints out key information about waypoints,
 *          home position, and home position in NED frame.
 * 
 * @param terrain Terrain mesh
 * @param path Path to diagnose
 * @param home_lat Home latitude
 * @param home_lon Home longitude
 * @param home_alt Home altitude
 */
void diagnosePathForPX4(const TerrainMesh& terrain, 
                        const std::vector<NEDPoint>& path,
                        double home_lat, double home_lon, double home_alt)
{
    std::cout << "\n=== PX4 MISSION DIAGNOSTICS ===\n";
    std::cout << std::fixed << std::setprecision(8);
    
    std::cout << "Home Position:\n"
              << "  Lat: " << home_lat << "°\n"
              << "  Lon: " << home_lon << "°\n"
              << "  Alt: " << home_alt << " m MSL\n\n";
    
    // Verify home position should be NED origin
    NED home_ned = latLonAltToNED(home_lat, home_lon, home_alt,
                                   terrain.getRefLat(), terrain.getRefLon(), terrain.getRefAlt());
    std::cout << "Home in terrain NED frame:\n"
              << "  North: " << home_ned.north << " m\n"
              << "  East:  " << home_ned.east << " m\n"
              << "  Down:  " << home_ned.down << " m\n\n";
    
    if (std::abs(home_ned.north) > 10 || std::abs(home_ned.east) > 10) {
        std::cout << "⚠️  WARNING: Home position is far from terrain reference!\n"
                  << "   This may cause alignment issues.\n\n";
    }
    
    std::cout << "Waypoint Analysis:\n";
    std::cout << "  Total waypoints: " << path.size() << "\n\n";
    
    for (size_t i = 0; i < std::min(path.size(), size_t(5)); ++i) {
        const auto& wp = path[i];
        
        // Convert back to LLA for verification
        double wp_lat, wp_lon, wp_alt;
        terrain.nedToLatLon(wp, wp_lat, wp_lon, wp_alt);
        
        // Get terrain elevation at this point
        double terrain_elev = terrain.getElevationAtLatLon(wp_lat, wp_lon);
        double agl = wp_alt - terrain_elev;
        
        std::cout << "Waypoint " << i << ":\n"
                  << "  NED: (" << wp.north << ", " << wp.east << ", " << wp.down << ")\n"
                  << "  LLA: (" << wp_lat << "°, " << wp_lon << "°, " << wp_alt << " m MSL)\n"
                  << "  Terrain: " << terrain_elev << " m MSL\n"
                  << "  AGL: " << agl << " m\n";
        
        if (agl < 10.0) {
            std::cout << "  ⚠️  WARNING: Less than 10m AGL!\n";
        }
        std::cout << "\n";
    }
    
    if (path.size() > 5) {
        std::cout << "... (" << (path.size() - 5) << " more waypoints)\n\n";
    }
    
    std::cout << "================================\n\n";
}

/**
 * @brief Debug altitude visualization
 * 
 * @details Prints out key information about the terrain mesh and path for altitude visualization.
 * @param terrain Terrain mesh
 * @param path Path to analyze
 */
void debugAltitudeVisualization(const TerrainMesh& terrain, const std::vector<NEDPoint>& path) {
    std::cout << "\n=== ALTITUDE VISUALIZATION DEBUG ===\n";
    std::cout << std::fixed << std::setprecision(3);
    
    // Terrain reference (NED origin)
    std::cout << "Terrain Reference (NED origin):\n";
    std::cout << "  LLA: (" << terrain.getRefLat() << ", " 
              << terrain.getRefLon() << ", " << terrain.getRefAlt() << " m MSL)\n";
    std::cout << "  NED: (0, 0, 0) by definition\n";
    std::cout << "  OpenGL Y: 0 (after centering)\n\n";
    
    // Terrain center point (visualization origin)
    glm::vec3 center = terrain.getCenterPoint();
    std::cout << "Terrain Center Point (OpenGL origin after translate):\n";
    std::cout << "  NED: (" << center.x << ", " << center.z << ", " 
              << -center.y << ") [note: Y is flipped]\n";
    std::cout << "  This gets subtracted from all positions!\n\n";
    
    // First waypoint
    if (!path.empty()) {
        const auto& wp0 = path[0];
        
        std::cout << "First Waypoint Analysis:\n";
        std::cout << "  Stored NED: (" << wp0.north << ", " << wp0.east 
                  << ", " << wp0.down << ")\n";
        
        // Convert to LLA
        double lat, lon, alt;
        terrain.nedToLatLon(wp0, lat, lon, alt);
        std::cout << "  LLA: (" << lat << ", " << lon << ", " << alt << " m MSL)\n";
        
        // What OpenGL sees BEFORE centering
        float opengl_x = wp0.north;
        float opengl_y = -wp0.down;  // Flip for Y-up
        float opengl_z = wp0.east;
        std::cout << "  OpenGL (before centering): (" << opengl_x << ", " 
                  << opengl_y << ", " << opengl_z << ")\n";
        
        // What OpenGL sees AFTER centering
        float final_x = opengl_x - center.x;
        float final_y = opengl_y - center.y;
        float final_z = opengl_z - center.z;
        std::cout << "  OpenGL (after centering): (" << final_x << ", " 
                  << final_y << ", " << final_z << ")\n";
        
        // Altitude analysis
        double terrain_alt = terrain.getElevationAtLatLon(lat, lon);
        double agl = alt - terrain_alt;
        std::cout << "  Terrain elevation: " << terrain_alt << " m MSL\n";
        std::cout << "  AGL: " << agl << " m\n";
        std::cout << "  Expected 'down' value: " << (terrain.getRefAlt() - alt) << "\n";
        std::cout << "  Actual 'down' value: " << wp0.down << "\n";
        double down_error = std::abs((terrain.getRefAlt() - alt) - wp0.down);
        if (down_error > 1.0) {
            std::cout << "  ⚠️  ERROR: Down value mismatch by " << down_error << "m!\n";
        }
    }
    
    std::cout << "\n=== Key Insight ===\n";
    std::cout << "In OpenGL visualization:\n";
    std::cout << "  - Terrain reference (home) appears at Y = " << -center.y << "\n";
    std::cout << "  - If center.y is NOT zero, home is offset!\n";
    std::cout << "  - Vehicle at 'down=0' should appear at Y = " << -center.y << "\n";
    std::cout << "  - Vehicle at 'down=-252' should appear at Y = " 
              << (-center.y + 252) << "\n";
    std::cout << "===================\n\n";
}

/**
 * @brief Main entry point
 * 
 * @param argc Number of command-line arguments
 * @param argv Array of command-line arguments
 */
int main(int argc, char** argv)
{
    #ifdef USE_MPI
        // Initialize MPI
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        is_mpi_initialized = true;
        
        if (mpi_rank == 0) {
            std::cout << "[MPI] Running with " << mpi_size << " process(es)\n";
            printHelp();
        }
        
        // Only rank 0 does visualization, all ranks do path planning
        bool do_visualization = (mpi_rank == 0);
    #else
        std::cout << "[Single-threaded] Running without MPI\n";
        printHelp();
        bool do_visualization = true;
    #endif
    
    // ============================================================================================
    // Initialize OpenGL (only on rank 0)
    // ============================================================================================
    GLFWwindow* window = nullptr;
        
    if (do_visualization) {
        // [Keep all your existing GLFW/OpenGL initialization code]
        if (!glfwInit()) {
            std::cerr << "Failed to initialize GLFW\n";
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#endif
            return -1;
        }
        
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_SAMPLES, 4);
        
        window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT,
                                  "Terrain Visualizer", nullptr, nullptr);
        if (!window) {
            std::cerr << "Failed to create GLFW window\n";
            glfwTerminate();
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#endif
            return -1;
        }
        
        glfwMakeContextCurrent(window);
        glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
        glfwSetCursorPosCallback(window, mouse_callback);
        glfwSetMouseButtonCallback(window, mouse_button_callback);
        glfwSetScrollCallback(window, scroll_callback);
        glfwSetKeyCallback(window, key_callback);
        
        glewExperimental = GL_TRUE;
        if (glewInit() != GLEW_OK) {
            std::cerr << "Failed to initialize GLEW\n";
#ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, -1);
#endif
            return -1;
        }
        
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_MULTISAMPLE);
        glClearColor(0.53f, 0.81f, 0.92f, 1.0f);
    }

    // ============================================================================================
    // Terrain Setup (ALL ranks need this for path planning)
    // ============================================================================================

#ifdef USE_MPI
    if (mpi_rank == 0) {
        std::cout << "Initializing terrain server...\n";
    }
#else
    std::cout << "Initializing terrain server...\n";
#endif

    std::unique_ptr<TerrainServerOpenTopo> pTopoServer = 
        std::make_unique<TerrainServerOpenTopo>("./terrain_cache", 
                                               "https://api.opentopodata.org/v1", 
                                               "srtm90m");
    std::string geotiff_path = "/home/yoshi/ECE6122_custom_final_project/TerrainServer/data/DEME_Zone3_2021_lla.tif";
    TerrainServerGeoTIFF server(geotiff_path);
    server.setInterpolationTolerance(0.1);
    pTopoServer->setBatchSize(100);
    server.setTopoServer(std::move(pTopoServer));
    
    // Create terrain mesh (all ranks)
    TerrainMesh terrain;

    // Parse command line arguments or use defaults
    double latStart =  36.381030;
    double lonStart = -111.857478;
    double latEnd =  36.400147;
    double lonEnd = -111.842598;
    int rows = 50;
    int cols = 50;

    if (argc >= 5)
    {
        latStart = std::stod(argv[1]);
        lonStart = std::stod(argv[2]);
        latEnd = std::stod(argv[3]);
        lonEnd = std::stod(argv[4]);
    }
    if (argc >= 7)
    {
        rows = std::stoi(argv[5]);
        cols = std::stoi(argv[6]);
    }

#ifdef USE_MPI
    if (mpi_rank == 0) {
#endif
        std::cout << "Loading terrain: [" << latStart << ", " << lonStart << "] to ["
                << latEnd << ", " << lonEnd << "] with " << rows << "x" << cols << " grid\n";
#ifdef USE_MPI
    }
#endif

#ifdef USE_MPI
    if (mpi_rank == 0) {
        // Full terrain mesh with OpenGL for visualization
        terrain.generateFromServer(server, latStart, lonStart, latEnd, lonEnd, rows, cols);
    } else {
        // Workers only need height data, not OpenGL mesh
        terrain.generateDataOnly(server, latStart, lonStart, latEnd, lonEnd, rows, cols);
    }
#else
    terrain.generateFromServer(server, latStart, lonStart, latEnd, lonEnd, rows, cols);
#endif

    // ============================================================================================
    // Path Planning with MPI
    // ============================================================================================

    double start_lat = 36.384218;
    double start_lon = -111.853969;
    double terrain_alt_at_start = terrain.getElevationAtLatLon(start_lat, start_lon);
    
    double goal_lat = 36.397386;
    double goal_lon = -111.854242;
    double goal_terrain_alt = terrain.getElevationAtLatLon(goal_lat, goal_lon);
    
    double flight_agl = 100.0;
    NEDPoint start_ned = terrain.latLonToNEDWithOffset(start_lat, start_lon, flight_agl);
    NEDPoint goal_ned = terrain.latLonToNEDWithOffset(goal_lat, goal_lon, flight_agl);
    
    // Create path planner (all ranks)
    TerrainPathPlanner::Config planner_config;
    planner_config.max_altitude_agl = 200.0;
    planner_config.min_altitude_agl = -20.0;
    planner_config.grid_resolution = 30.0;
    planner_config.altitude_step = 30.0;
    planner_config.use_diagonal_moves = true;
    planner_config.max_iterations = 300000;
    planner_config.vertical_safety_margin = 0.0;
    planner_config.horizontal_safety_radius = 0.0;
    planner_config.max_terrain_gradient = 100.0;
    
    TerrainPathPlanner planner(
        terrain.getHeightGrid(),
        terrain.getNorthMin(),
        terrain.getNorthMax(),
        terrain.getEastMin(),
        terrain.getEastMax(),
        planner_config
    );
    
#ifdef USE_MPI
    if (mpi_rank == 0) {
#endif
        planner.debugValidPoint(start_ned, "START");
        planner.debugValidPoint(goal_ned, "GOAL");
#ifdef USE_MPI
    }
#endif
    
    // Compute path (all ranks participate)
    std::vector<NEDPoint> path;
    
#ifdef USE_MPI
    if (mpi_size > 1) {
        // Use MPI-parallelized version
        path = planner.computePathMPI(start_ned, goal_ned);
    } else {
        // Fall back to sequential even with MPI compiled
        if (mpi_rank == 0) {
            path = planner.computePath(start_ned, goal_ned);
        }
    }
#else
    // Non-MPI build uses sequential version
    path = planner.computePath(start_ned, goal_ned);
#endif

    // ============================================================================================
    // Visualization Setup (only rank 0)
    // ============================================================================================
    std::vector<Cube> terrainCubes;
    std::vector<Cube> pathCubes;
    
    if (do_visualization && !path.empty()) {
        // [Keep all your existing visualization code]
        auto stats = planner.getLastPathStats();
        std::cout << "\n=== Path Statistics ===\n";
        std::cout << "Total distance: " << stats.total_distance << " m\n";
        std::cout << "Altitude change: " << stats.total_altitude_change << " m\n";
        std::cout << "Average AGL: " << stats.avg_agl << " m\n";
        std::cout << "Max AGL: " << stats.max_agl << " m\n";
        std::cout << "Waypoints: " << stats.num_waypoints << "\n";
        std::cout << "Computation time: " << stats.computation_time << " s\n";
        
        // Create marker cubes
        NEDPoint start_marker = terrain.latLonToNEDWithOffset(start_lat, start_lon, 100.0);
        NEDPoint goal_marker = terrain.latLonToNEDWithOffset(goal_lat, goal_lon, 100.0);
        
        terrainCubes.emplace_back(start_marker, glm::vec3(1.0f, 0.0f, 0.0f), 50.0f);
        terrainCubes.emplace_back(goal_marker, glm::vec3(0.0f, 1.0f, 0.0f), 50.0f);
        
        // Create path cubes
        pathCubes.emplace_back(path[0], glm::vec3(1.0f, 1.0f, 0.0f), 30.0f);
        
        for (size_t i = 1; i < path.size() - 1; ++i) {
            float t = static_cast<float>(i) / (path.size() - 1);
            glm::vec3 color(1.0f, 1.0f - t * 0.5f, 0.5f * t);
            pathCubes.emplace_back(path[i], color, 5.0f);
        }
        
        if (path.size() > 1) {
            pathCubes.emplace_back(path[path.size() - 1], glm::vec3(0.0f, 1.0f, 1.0f), 30.0f);
        }
        
        for (auto& cube : pathCubes) {
            terrainCubes.push_back(std::move(cube));
        }

        std::cout << "\nVisualization:\n";
        std::cout << "  - RED cube (50m): Start location marker (100m above ground)\n";
        std::cout << "  - GREEN cube (50m): Goal location marker (100m above ground)\n";
        std::cout << "  - YELLOW cube (30m): Path start waypoint\n";
        std::cout << "  - " << (path.size() - 2) << " small cubes: Path waypoints\n";
        std::cout << "  - CYAN cube (30m): Path end waypoint\n";
    } 
    else if (path.empty())
    {
        std::cerr << "\n!!! FAILED TO FIND PATH !!!\n";
        std::cerr << "Possible reasons:\n";
        std::cerr << "  1. Start/goal altitudes outside min/max AGL bounds\n";
        std::cerr << "  2. Terrain too steep between points\n";
        std::cerr << "  3. Grid resolution too coarse\n";
        std::cerr << "  4. Max iterations too low\n";
        std::cerr << "\nYou should still see RED (start) and GREEN (goal) marker cubes!\n";
    }

    // ============================================================================================
    // Render Loop (only rank 0)
    // ============================================================================================

    if (do_visualization && window) {
        // [Keep all your existing render loop code]
        camera.setTarget(terrain.getCenterPoint());
        Shader shader = Shader::createDefaultTerrainShader();
        Shader cubeShader = Shader::createCubeShader();
        
        // ... [Keep all PX4 integration and render loop code] ...

        
        // ============================================================================================
        // Camera
        // ============================================================================================
        // Setup camera to look at terrain center
        camera.setTarget(terrain.getCenterPoint());

        // ============================================================================================
        // PX4 Integration
        // ============================================================================================
        
        std::cout << "\n=== PX4 Integration ===\n";
        
        px4_interface = std::make_unique<PX4Interface>();
        
        // Connect to PX4 on port 14541 (QGC uses 14550, SITL default is 14540)
        // We'll use 14541 to avoid conflicts
        if (px4_interface->connect("udp://:14541")) {
            
            // Upload the computed path as a mission
            if (!path.empty()) {
                std::cout << "\n=== Converting path to start-relative coordinates ===\n";
                
                // Calculate offset between terrain center and start location
                NED offset = latLonAltToNED(
                    start_lat, start_lon, terrain_alt_at_start,                     // Start location (NEW home)
                    terrain.getRefLat(), terrain.getRefLon(), terrain.getRefAlt()  // Terrain center (OLD reference)
                );

                px4_to_terrain_offset = offset;  // ← ADD THIS LINE
                px4_to_terrain_offset.north = -offset.north;  // ← NEGATE to go opposite direction
                px4_to_terrain_offset.east = -offset.east;
                px4_to_terrain_offset.down = -offset.down;
                
                std::cout << "Offset from start to terrain center:\n";
                std::cout << "  North: " << offset.north << " m\n";
                std::cout << "  East: " << offset.east << " m\n";
                std::cout << "  Down: " << offset.down << " m\n\n";
                
                // Convert all path waypoints to be relative to start location
                std::vector<NEDPoint> path_from_start;
                path_from_start.reserve(path.size());
                
                for (const auto& wp : path) {
                    // Subtract the offset to convert from terrain-center-relative to start-relative
                    NEDPoint wp_from_start;
                    wp_from_start.north = wp.north - offset.north;
                    wp_from_start.east = wp.east - offset.east;
                    wp_from_start.down = wp.down - offset.down;
                    path_from_start.push_back(wp_from_start);
                }
                
                std::cout << "Converted " << path_from_start.size() << " waypoints to start-relative frame\n";
                std::cout << "First waypoint (terrain-relative): (" << path[0].north << ", " 
                        << path[0].east << ", " << path[0].down << ")\n";
                std::cout << "First waypoint (start-relative): (" << path_from_start[0].north << ", " 
                        << path_from_start[0].east << ", " << path_from_start[0].down << ")\n\n";
                
                // Run diagnostics with start-relative path
                diagnosePathForPX4(terrain, path_from_start, start_lat, start_lon, terrain_alt_at_start);

                // Use terrain center or start position as home
                bool mission_uploaded = px4_interface->uploadMission(
                    path_from_start, 
                    start_lat,   // Home latitude
                    start_lon,   // Home longitude
                    terrain_alt_at_start    // Home altitude
                );
                
                if (mission_uploaded) {
                    std::cout << "[PX4] Mission uploaded successfully!\n";
                }
            } else {
                std::cout << "[PX4] No path to upload (path planning failed)\n";
            }
            
            // Create vehicle visualization cube (BLUE, 20m size)
            NEDPoint initial_pos(0, 0, 0); // Will be updated by telemetry
            vehicleCubes.emplace_back(initial_pos, glm::vec3(0.0f, 0.0f, 1.0f), 20.0f);
            
            // Set up real-time position callback
            px4_interface->setPositionCallback(onVehiclePositionUpdate);
            
            std::cout << "[PX4] Ready! Use keys 1-5 to control vehicle\n";
            std::cout << "  BLUE cube: Real-time vehicle position\n";
            
        } else {
            std::cout << "[PX4] Failed to connect - continuing without PX4\n";
            px4_interface.reset();
        }

        // ============================================================================================
        // Render Loop (UPDATE)
        // ============================================================================================
        
        // Render loop
        std::cout << "\nStarting render loop...\n";
        
        // Timing variables (move outside loop)
        float lastFrame = 0.0f;
        double lastFPSTime = glfwGetTime();
        int frameCount = 0;
        
        while (!glfwWindowShouldClose(window))
        {
            // Calculate delta time
            float currentFrame = static_cast<float>(glfwGetTime());
            deltaTime = currentFrame - lastFrame;
            lastFrame = currentFrame;
            
            // Process continuous input (held keys)
            processInput(window);

            // Process PX4 controls
            processPX4Input(window);
            
            // Calculate FPS
            frameCount++;
            if (currentFrame - lastFPSTime >= 1.0)
            {
                std::cout << "FPS: " << frameCount << " | Camera Pos: (" 
                        << std::fixed << std::setprecision(1)
                        << camera.getPosition().x << ", "
                        << camera.getPosition().y << ", "
                        << camera.getPosition().z << ")    \r" << std::flush;
                frameCount = 0;
                lastFPSTime = currentFrame;
            }

            // Clear buffers
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // ========================================================================
            // Render TERRAIN
            // ========================================================================
            shader.use();

            // Set up transforms
            glm::mat4 model = glm::mat4(1.0f);
            model = glm::translate(model, -terrain.getCenterPoint());
            
            glm::mat4 view = camera.getViewMatrix();
            
            glm::mat4 projection = glm::perspective(
                glm::radians(45.0f),
                static_cast<float>(WINDOW_WIDTH) / static_cast<float>(WINDOW_HEIGHT),
                0.1f, 10000.0f  // Increased far plane for larger terrain
            );

            shader.setMat4("model", model);
            shader.setMat4("view", view);
            shader.setMat4("projection", projection);

            // Set lighting
            glm::vec3 lightPos = camera.getPosition() + glm::vec3(10.0f, 20.0f, 10.0f);
            shader.setVec3("lightPos", lightPos);
            shader.setVec3("viewPos", camera.getPosition());
            shader.setVec3("lightColor", glm::vec3(1.0f, 1.0f, 0.95f));

            // Render terrain
            terrain.render(false);

            // ========================================================================
            // Render CUBES
            // ========================================================================
            cubeShader.use();
            
            // Reuse same model transform (terrain centering)
            glm::mat4 cubeModel = glm::mat4(1.0f);
            cubeModel = glm::translate(cubeModel, -terrain.getCenterPoint());
            
            cubeShader.setMat4("view", view);
            cubeShader.setMat4("projection", projection);
            cubeShader.setMat4("model", cubeModel);  // Set base model matrix
            
            // Render all cubes
            for (auto& cube : terrainCubes) 
            {
                cube.render(cubeModel);
            }

            // Render vehicle cube (thread-safe)
            {
                std::lock_guard<std::mutex> lock(vehicle_mutex);
                for (auto& cube : vehicleCubes) 
                {
                    cube.render(cubeModel);
                }
            }

            // Swap buffers and poll events
            glfwSwapBuffers(window);
            glfwPollEvents();
        }

        // Cleanup
        if (px4_interface) 
        {
            px4_interface->disconnect();
        }

        // Cleanup
        glfwTerminate();
    }

#ifdef USE_MPI
    // Cleanup MPI
    if (is_mpi_initialized) {
        MPI_Finalize();
    }
#endif
    
    std::cout << "\nExiting...\n";
    return 0;
}