#pragma once

#include <fstream>
#include <mutex>
#include <chrono>
#include <ctime>
#include <sstream>

class Logger
{
    public:
        Logger(const std::string& filename) : logFile_(filename, std::ios::app) {}

        void log(const std::string& message);

    private:
        std::string timestamp();

        std::ofstream logFile_;
        std::mutex mutex_;
};