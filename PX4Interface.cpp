/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @brief PX4 Interface Implementation. 
 * 
 * @details This file implements the PX4Interface class for connecting to a PX4-based UAV,
 * uploading missions, controlling flight, and retrieving telemetry data using MAVSDK.
 */

#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <future>
#include <cmath>

#include "PX4Interface.hpp"
#include "NEDPoint.hpp"

using namespace mavsdk;

/**
 * @brief Constructor
 * 
 * @details Creates a PX4Interface instance.
 */
PX4Interface::PX4Interface()
    : connected_(false)
    , home_lat_(0)
    , home_lon_(0)
    , home_alt_(0)
    , home_set_(false)
    , stop_telemetry_(false)
{
    // CRITICAL: Use GCS system ID, NOT autopilot ID
    // System ID 255 = Ground Control Station
    // Component ID 190 = MAV_COMP_ID_MISSIONPLANNER
    Mavsdk::Configuration config(
        255,    // system_id (was 1 - WRONG!)
        190,    // component_id (was 1 - WRONG!)
        true    // always_send_heartbeats
    );
    mavsdk_ = std::make_shared<Mavsdk>(config);
}

/**
 * @brief Destructor
 * 
 * @details Destroys the PX4Interface instance.
 */
PX4Interface::~PX4Interface()
{
    disconnect();
}

/**
 * @brief Get Vehicle Type
 * 
 * @details Returns the type of the connected vehicle.
 * 
 * @return VehicleType Enum indicating vehicle type.
 */
PX4Interface::VehicleType PX4Interface::getVehicleType() const
{
    if (!connected_) return VehicleType::UNKNOWN;
    
    // Get vehicle type from telemetry
    // This requires MAVSDK to expose vehicle type
    // For now, you can manually set this based on what you launched
    return VehicleType::FIXED_WING; // Or MULTIROTOR
}

/**
 * @brief Connect to a PX4-based UAV
 * 
 * @details Establishes a connection to a PX4-based UAV using MAVSDK.
 * 
 * @param connection_url The connection URL (default: "udp://:14540").
 * @return true if connection is successful, false otherwise.
 */
bool PX4Interface::connect(const std::string& connection_url)
{
    std::cout << "[PX4] Connecting to " << connection_url << "...\n";
    
    ConnectionResult connection_result = mavsdk_->add_any_connection(connection_url);
    
    if (connection_result != ConnectionResult::Success) {
        std::cerr << "[PX4] Connection failed: " << connection_result << "\n";
        return false;
    }

    std::cout << "[PX4] Waiting for system to connect...\n";
    
    // Wait for system with autopilot (increased timeout to 15 seconds)
    auto start_time = std::chrono::steady_clock::now();
    std::shared_ptr<System> found_system = nullptr;
    
    while (std::chrono::steady_clock::now() - start_time < std::chrono::seconds(15)) {
        // Check all discovered systems
        for (auto system : mavsdk_->systems()) {
            if (system->has_autopilot()) {
                found_system = system;
                std::cout << "[PX4] Discovered autopilot (System ID: " 
                          << (int)system->get_system_id() << ")\n";
                break;
            }
        }
        
        if (found_system) break;
        
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    
    if (!found_system) {
        std::cerr << "[PX4] No autopilot found (timeout)\n";
        std::cerr << "[PX4] Discovered " << mavsdk_->systems().size() << " systems:\n";
        for (auto sys : mavsdk_->systems()) {
            std::cerr << "  - System ID " << (int)sys->get_system_id() 
                      << " (has autopilot: " << sys->has_autopilot() << ")\n";
        }
        return false;
    }

    system_ = found_system;

    // Instantiate plugins
    action_ = std::make_shared<Action>(system_);
    mission_ = std::make_shared<Mission>(system_);
    telemetry_ = std::make_shared<Telemetry>(system_);

    setupTelemetrySubscriptions();
    
    connected_ = true;
    std::cout << "[PX4] Connected successfully!\n";
    
    return true;
}

/**
 * @brief Disconnect from the PX4-based UAV
 * 
 * @details Cleans up and disconnects from the UAV.
 */
void PX4Interface::disconnect()
{
    if (!connected_) return;
    
    stop_telemetry_ = true;
    if (telemetry_thread_.joinable()) {
        telemetry_thread_.join();
    }
    
    connected_ = false;
    std::cout << "[PX4] Disconnected\n";
}

/**
 * @brief Set the home position
 * 
 * @details Sets the home position for NED coordinate calculations.
 * 
 * @param lat Latitude in degrees.
 * @param lon Longitude in degrees.
 * @param alt Altitude in meters (AMSL).
 */
void PX4Interface::setHomePosition(double lat, double lon, double alt)
{
    home_lat_ = lat;
    home_lon_ = lon;
    home_alt_ = alt;
    home_set_ = true;
    std::cout << "[PX4] Home position set: " << lat << ", " << lon << ", " << alt << "m\n";
}

/**
 * @brief Set up subscriptions for telemetry data
 * 
 * @details Sets up subscriptions for position, attitude, and velocity data from the UAV.
 */
void PX4Interface::setupTelemetrySubscriptions()
{
    // Set telemetry rate
    telemetry_->set_rate_position(10.0); // 10 Hz
    telemetry_->set_rate_attitude_euler(10.0);
    telemetry_->set_rate_velocity_ned(10.0);

    // Start telemetry thread
    stop_telemetry_ = false;
    telemetry_thread_ = std::thread(&PX4Interface::telemetryLoop, this);
}

/**
 * @brief Telemetry loop
 * 
 * @details Continuously retrieves telemetry data and updates the vehicle state.
 */
void PX4Interface::telemetryLoop()
{
    while (!stop_telemetry_) {
        VehicleState state;
        state.timestamp = std::chrono::system_clock::now();

        // Get position (no .value() needed in older API)
        auto pos = telemetry_->position();
        state.latitude = pos.latitude_deg;
        state.longitude = pos.longitude_deg;
        state.altitude_amsl = pos.absolute_altitude_m;
        state.altitude_agl = pos.relative_altitude_m;
        
        // Convert to NED if home is set
        if (home_set_) {
            updateNEDFromGPS(pos.latitude_deg, pos.longitude_deg, 
                           pos.absolute_altitude_m, state);
        }

        // Get attitude (no .value() needed)
        auto att = telemetry_->attitude_euler();
        state.roll = att.roll_deg;
        state.pitch = att.pitch_deg;
        state.yaw = att.yaw_deg;

        // Get velocity (no .value() needed)
        auto vel = telemetry_->velocity_ned();
        state.velocity_north = vel.north_m_s;
        state.velocity_east = vel.east_m_s;
        state.velocity_down = vel.down_m_s;

        // Get armed/in-air status
        state.is_armed = telemetry_->armed();
        state.in_air = telemetry_->in_air();

        // Update current state (thread-safe)
        {
            std::lock_guard<std::mutex> lock(state_mutex_);
            current_state_ = state;
        }

        // Call callback if set
        if (position_callback_) {
            position_callback_(state);
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // 10 Hz
    }
}

/**
 * @brief Update NED coordinates from GPS
 * 
 * @details Converts GPS coordinates to NED coordinates relative to home position.
 */
void PX4Interface::updateNEDFromGPS(double lat, double lon, double alt, VehicleState& state)
{
    // Convert lat/lon to NED relative to home position
    double dlat = (lat - home_lat_) * DEG_TO_RAD;
    double dlon = (lon - home_lon_) * DEG_TO_RAD;
    
    // Approximate conversion (flat Earth, good for small distances)
    state.north = dlat * EARTH_RADIUS;
    state.east = dlon * EARTH_RADIUS * std::cos(home_lat_ * DEG_TO_RAD);
    state.down = -(alt - home_alt_); // NED down is negative up
}

/**
 * @brief Get the current vehicle state
 */
PX4Interface::VehicleState PX4Interface::getVehicleState() const
{
    std::lock_guard<std::mutex> lock(state_mutex_);
    return current_state_;
}

/**
 * @brief Set the position callback
 * 
 * @details Sets the callback function to be called when the vehicle position changes.
 * 
 * @param callback The callback function.
 */
void PX4Interface::setPositionCallback(std::function<void(const VehicleState&)> callback)
{
    position_callback_ = callback;
}

/**
 * @brief Upload a mission to the vehicle
 * 
 * @details Uploads a mission to the vehicle using the PX4 command interface.
 * 
 * @param waypoints The waypoints to upload.
 * @param home_lat The latitude of the home position.
 * @param home_lon The longitude of the home position.
 * @param home_alt The altitude of the home position.
 */
bool PX4Interface::uploadMission(const std::vector<NEDPoint>& waypoints, 
                                 double home_lat, double home_lon, double home_alt)
{
    if (!connected_ || waypoints.empty()) {
        std::cerr << "[PX4] Cannot upload mission: not connected or empty waypoints\n";
        return false;
    }

    float speed = 10.0f;
    float acceptance_radius = 10.0f;
    
    auto vehicle_type = getVehicleType();
    if (vehicle_type == VehicleType::FIXED_WING) {
        speed = 20.0f;
        acceptance_radius = 50.0f;
        std::cout << "[PX4] Using fixed-wing mission parameters\n";
    } else if (vehicle_type == VehicleType::MULTIROTOR) {
        speed = 10.0f;
        acceptance_radius = 5.0f;
        std::cout << "[PX4] Using multirotor mission parameters\n";
    }

    std::cout << "[PX4] Uploading mission with " << waypoints.size() << " waypoints...\n";

    // Set home position for NED conversion
    setHomePosition(home_lat, home_lon, home_alt);

    std::vector<Mission::MissionItem> mission_items;

    // 1. Add the Path Waypoints
    for (size_t i = 0; i < waypoints.size(); ++i) {
        const auto& wp = waypoints[i];
        
        // Convert NED back to GPS coordinates
        double lat = home_lat + (wp.north / EARTH_RADIUS) / DEG_TO_RAD;
        double lon = home_lon + (wp.east / (EARTH_RADIUS * std::cos(home_lat * DEG_TO_RAD))) / DEG_TO_RAD;
        double alt_amsl = home_alt - wp.down; // NED down is negative up

        Mission::MissionItem item;
        item.latitude_deg = lat;
        item.longitude_deg = lon;
        item.relative_altitude_m = alt_amsl - home_alt; // Relative to home
        item.speed_m_s = speed; 
        item.is_fly_through = true;
        item.loiter_time_s = 0.0f; 
        item.gimbal_pitch_deg = -90.0f;
        item.gimbal_yaw_deg = 0.0f;
        item.camera_action = Mission::MissionItem::CameraAction::None;
        item.camera_photo_interval_s = 0.0;
        item.acceptance_radius_m = acceptance_radius; 
        
        mission_items.push_back(item);

        std::cout << "  WP" << i << ": lat=" << std::fixed << std::setprecision(6) << lat 
                  << ", lon=" << lon 
                  << ", rel_alt=" << std::setprecision(1) << (alt_amsl - home_alt) << "m"
                  << ", speed=" << item.speed_m_s << "m/s\n";
    }

    // ==========================================================
    // 2. Add Landing Sequence (Return to Home + Land)
    // ==========================================================

    // Optional: Add approach waypoint at home with some altitude
    Mission::MissionItem approach_item;
    approach_item.latitude_deg = home_lat;
    approach_item.longitude_deg = home_lon;
    approach_item.relative_altitude_m = 50.0f; // Approach altitude
    approach_item.speed_m_s = (vehicle_type == VehicleType::FIXED_WING) ? 15.0f : 5.0f;
    approach_item.is_fly_through = true;
    approach_item.acceptance_radius_m = acceptance_radius;
    approach_item.vehicle_action = Mission::MissionItem::VehicleAction::None;

    mission_items.push_back(approach_item);

    // Then add the land item as shown above
    if (!waypoints.empty()) {
        // Land at HOME position (QGC typically requires this)
        Mission::MissionItem land_item;
        land_item.latitude_deg = home_lat;      // HOME, not last waypoint
        land_item.longitude_deg = home_lon;     // HOME, not last waypoint
        land_item.relative_altitude_m = 0.0f;   // Ground level
        land_item.speed_m_s = (vehicle_type == VehicleType::FIXED_WING) ? 12.0f : 2.0f;
        land_item.is_fly_through = false;       // Must stop here
        land_item.acceptance_radius_m = acceptance_radius;
        
        // THIS IS THE KEY LINE YOU WERE MISSING:
        land_item.vehicle_action = Mission::MissionItem::VehicleAction::Land;
        
        land_item.loiter_time_s = 0.0f;
        land_item.gimbal_pitch_deg = 0.0f;
        land_item.gimbal_yaw_deg = 0.0f;
        land_item.camera_action = Mission::MissionItem::CameraAction::None;

        mission_items.push_back(land_item);
        std::cout << "[PX4] Added LANDING waypoint at HOME (Lat: " << home_lat 
                << ", Lon: " << home_lon << ") with VehicleAction::Land\n";
    }

    Mission::MissionPlan mission_plan{};
    mission_plan.mission_items = mission_items;

    const Mission::Result upload_result = mission_->upload_mission(mission_plan);

    if (upload_result != Mission::Result::Success) {
        std::cerr << "[PX4] Mission upload failed: " << upload_result << "\n";
        return false;
    }

    std::cout << "[PX4] Mission uploaded successfully!\n";
    return true;
}

/**
 * @brief Arms the vehicle and takes it off
 * 
 * @details Arms the vehicle and commands it to take off to a specified altitude.
 */
bool PX4Interface::armAndTakeoff(float altitude_m)
{
    if (!connected_) {
        std::cerr << "[PX4] Not connected\n";
        return false;
    }

    std::cout << "[PX4] Arming...\n";
    const Action::Result arm_result = action_->arm();
    if (arm_result != Action::Result::Success) {
        std::cerr << "[PX4] Arming failed: " << arm_result << "\n";
        return false;
    }

    // For fixed-wing: Set takeoff altitude first
    std::cout << "[PX4] Setting takeoff altitude to " << altitude_m << "m...\n";
    action_->set_takeoff_altitude(altitude_m);

    std::cout << "[PX4] Taking off...\n";
    const Action::Result takeoff_result = action_->takeoff();
    if (takeoff_result != Action::Result::Success) {
        std::cerr << "[PX4] Takeoff failed: " << takeoff_result << "\n";
        return false;
    }

    std::cout << "[PX4] Takeoff command sent (fixed-wing will use catapult/runway)\n";
    return true;
}

/**
 * @brief Starts the mission
 */
bool PX4Interface::startMission()
{
    if (!connected_) {
        std::cerr << "[PX4] Not connected\n";
        return false;
    }

    std::cout << "[PX4] Starting mission...\n";
    const Mission::Result result = mission_->start_mission();
    
    if (result != Mission::Result::Success) {
        std::cerr << "[PX4] Mission start failed: " << result << "\n";
        return false;
    }

    std::cout << "[PX4] Mission started!\n";
    return true;
}

/**
 * @brief Pauses the mission
 */
bool PX4Interface::pauseMission()
{
    if (!connected_) return false;
    
    std::cout << "[PX4] Pausing mission...\n";
    const Mission::Result result = mission_->pause_mission();
    return result == Mission::Result::Success;
}

/**
 * @brief Command returns to launch
 */
bool PX4Interface::returnToLaunch()
{
    if (!connected_) return false;
    
    std::cout << "[PX4] Returning to launch...\n";
    const Action::Result result = action_->return_to_launch();
    return result == Action::Result::Success;
}

/**
 * @brief Command landing
 */
bool PX4Interface::land()
{
    if (!connected_) return false;
    
    std::cout << "[PX4] Landing...\n";
    const Action::Result result = action_->land();
    return result == Action::Result::Success;
}