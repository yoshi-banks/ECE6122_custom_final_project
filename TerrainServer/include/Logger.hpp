/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file Logger.hpp
 */

#pragma once

#include <string>
#include <fstream>
#include <mutex>

/**
 * @brief Logger class
 * 
 * @details This class provides simple logging functionality to log messages to a file with timestamps.
 */
class Logger
{
    public:
        Logger(const std::string& filename);

        void log(const std::string& message);

    private:
        std::string timestamp();

        std::ofstream logFile_;
        std::mutex mutex_;
};