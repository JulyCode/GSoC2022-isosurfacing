#pragma once

#include <iostream>
#include <chrono>

class ScopeTimer {
    using TimePoint = std::chrono::_V2::system_clock::time_point;
public:
    ScopeTimer() : start(std::chrono::high_resolution_clock::now()), msg("Duration") {}

    explicit ScopeTimer(const std::string& msg) : start(std::chrono::high_resolution_clock::now()), msg(msg) {
        std::cout << msg << "..." << std::endl;
    }

    ~ScopeTimer() {
        TimePoint end = std::chrono::high_resolution_clock::now();
        int64_t duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << msg << ": " << duration << " ms" << std::endl;
    }

private:
    const TimePoint start;
    const std::string msg;
};
