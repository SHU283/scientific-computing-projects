#pragma once
#include <chrono>

// timer
struct Timer {
    // Monotonically increasing clock
    using clock = std::chrono::steady_clock;   
    // The most recent start() time point                     
    clock::time_point t0;                            
    // cumulative time               
    double acc = 0.0;                                               

    // Start timing
    void start() { t0 = clock::now(); }   
    // End timing                          
    void stop()  {                                                  
        auto t1 = clock::now();
        acc += std::chrono::duration<double>(t1 - t0).count();
    }
    // Reset cumulative time
    void reset() { acc = 0.0; }                                     
};
