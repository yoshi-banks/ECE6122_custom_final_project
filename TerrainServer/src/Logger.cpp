#include "Logger.hpp"

void Logger::log(const std::string& message)
{
    std::lock_guard<std::mutex> lock(mutex_);
    if (logFile_.is_open())
    {
        logFile_ << timestamp() << " " << message << "\n";
        logFile_.flush();
    }
}

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