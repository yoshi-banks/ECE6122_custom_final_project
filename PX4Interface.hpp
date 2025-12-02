/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * @file PX4Interface.hpp
 * @brief PX4Interface class definition
 * 
 * @details This file declares the PX4Interface class for connecting to a PX4-based UAV,
 * uploading missions, controlling flight, and retrieving telemetry data using MAVSDK.
 */

#pragma once

#include <mavsdk/mavsdk.h>
#include <mavsdk/plugins/action/action.h>
#include <mavsdk/plugins/mission/mission.h>
#include <mavsdk/plugins/telemetry/telemetry.h>
#include <mavsdk/plugins/offboard/offboard.h>
#include <mavsdk/system.h>
#include <mavsdk/plugin_base.h>
#include <vector>
#include <memory>
#include <atomic>
#include <thread>
#include <functional>
#include <mutex>

struct NEDPoint;

/**
 * @brief PX4Interface class
 */
class PX4Interface
{
public:
    struct VehicleState {
        double latitude;
        double longitude;
        double altitude_amsl;  // Altitude above mean sea level
        double altitude_agl;   // Altitude above ground level (if available)
        double north;          // NED coordinates relative to home
        double east;
        double down;
        double roll;           // degrees
        double pitch;          // degrees
        double yaw;            // degrees
        double velocity_north;
        double velocity_east;
        double velocity_down;
        bool is_armed;
        bool in_air;
        std::chrono::system_clock::time_point timestamp;
        
        VehicleState() : latitude(0), longitude(0), altitude_amsl(0), altitude_agl(0),
                        north(0), east(0), down(0), roll(0), pitch(0), yaw(0),
                        velocity_north(0), velocity_east(0), velocity_down(0),
                        is_armed(false), in_air(false) {}
    };

    PX4Interface();
    ~PX4Interface();

    enum class VehicleType {
        UNKNOWN,
        MULTIROTOR,
        FIXED_WING,
        VTOL
    };

    VehicleType getVehicleType() const;

    // Connection
    bool connect(const std::string& connection_url = "udp://:14540");
    void disconnect();
    bool isConnected() const { return connected_; }

    // Mission upload
    bool uploadMission(const std::vector<NEDPoint>& waypoints, 
                      double home_lat, double home_lon, double home_alt);
    
    // Mission control
    bool armAndTakeoff(float altitude_m = 10.0f);
    bool startMission();
    bool pauseMission();
    bool returnToLaunch();
    bool land();

    // Telemetry - thread-safe access
    VehicleState getVehicleState() const;
    
    // Set home position for NED calculations
    void setHomePosition(double lat, double lon, double alt);

    // Callback for real-time updates (called from telemetry thread)
    void setPositionCallback(std::function<void(const VehicleState&)> callback);

private:
    std::shared_ptr<mavsdk::Mavsdk> mavsdk_;
    std::shared_ptr<mavsdk::System> system_;
    std::shared_ptr<mavsdk::Action> action_;
    std::shared_ptr<mavsdk::Mission> mission_;
    std::shared_ptr<mavsdk::Telemetry> telemetry_;
    
    std::atomic<bool> connected_;
    mutable std::mutex state_mutex_;
    VehicleState current_state_;
    
    // Home position for NED conversion
    double home_lat_;
    double home_lon_;
    double home_alt_;
    bool home_set_;
    
    std::function<void(const VehicleState&)> position_callback_;
    std::thread telemetry_thread_;
    std::atomic<bool> stop_telemetry_;

    // Helper functions
    void setupTelemetrySubscriptions();
    void telemetryLoop();
    void updateNEDFromGPS(double lat, double lon, double alt, VehicleState& state);
    
    // Coordinate conversion helpers
    // TODO move these to a constants file
    constexpr static double EARTH_RADIUS = 6378137.0; // meters
    constexpr static double DEG_TO_RAD = M_PI / 180.0;
};