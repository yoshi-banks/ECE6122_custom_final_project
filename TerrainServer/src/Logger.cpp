/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file Logger.cpp
 * @brief Logger class implementation
 */

#include "Logger.hpp"

#include <sstream>
#include <fstream>
#include <string>
#include <iostream>

/**
 * @brief Check if file exists
 * @param filename - Name of file to check
 */
static bool fileExists(const std::string& filename)
{
    std::ifstream infile(filename);
    return infile.good();
}

/**
 * @brief Constructor for Logger class
 * @param baseFilename - Base name of log file
 */
Logger::Logger(const std::string& baseFilename) 
{
    std::string filename = baseFilename;
    size_t count = 1;

    // Check if file exists and increment filename
    while (fileExists(filename))
    {
        // Insert (cout) before file extension
        // Assume extension like ".log"
        auto dotPos = baseFilename.rfind('.');
        std::ostringstream oss;
        if (dotPos == std::string::npos)
        {
            oss << baseFilename << count; // No extension
        }
        else
        {
            oss << baseFilename.substr(0, dotPos) << count << baseFilename.substr(dotPos);
        }
        filename = oss.str();
        count++;
    }

    logFile_.open(filename, std::ios::app); // Now open new file
}

/**
 * @brief Log message to console and file
 * @param message - Message to log
 */
void Logger::log(const std::string& message)
{
    std::cout << timestamp() << " " << message << "\n";

    std::lock_guard<std::mutex> lock(mutex_);
    if (logFile_.is_open())
    {
        logFile_ << timestamp() << " " << message << "\n";
        logFile_.flush();
    }
}

/**
 * @brief Get current timestamp as string
 * @return Timestamp string in format [YYYY-MM-DD HH:MM:SS]
 */
std::string Logger::timestamp()
{
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf;
    localtime_r(&now_c, &tm_buf);
    char buf[20];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tm_buf);
    return std::string("[") + buf + "]";
}

Logger logger("TerrainLog.log");