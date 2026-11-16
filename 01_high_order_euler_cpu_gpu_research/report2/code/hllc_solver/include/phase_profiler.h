#pragma once

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

namespace fp {

struct PhaseProfile {
    using Clock = std::chrono::steady_clock;
    using TimePoint = Clock::time_point;

    bool enabled = false;
    double total_s = 0.0;
    double setup_s = 0.0;
    double cfl_s = 0.0;
    double workspace_s = 0.0;
    double flux_x_s = 0.0;
    double flux_y_s = 0.0;
    double spatial_fused_s = 0.0;
    double rhs_assembly_s = 0.0;
    double state_update_s = 0.0;
    double admissibility_s = 0.0;
    double transfer_s = 0.0;
    int retries = 0;

    static bool requested() {
        const char* value = std::getenv("FP_PHASE_PROFILE");
        return value != nullptr && std::string(value) != "0";
    }

    static double elapsed(const TimePoint& start) {
        return std::chrono::duration<double>(Clock::now() - start).count();
    }

    double measured_s() const {
        return setup_s + cfl_s + workspace_s + flux_x_s + flux_y_s +
               spatial_fused_s + rhs_assembly_s + state_update_s +
               admissibility_s + transfer_s;
    }

    void print(const std::string& method,
               const std::string& backend,
               int nx,
               int ny,
               int steps) const {
        if (!enabled) {
            return;
        }
        const double unattributed = std::max(0.0, total_s - measured_s());
        std::cerr << std::setprecision(12)
                  << "PHASE_PROFILE"
                  << " method=" << method
                  << " backend=" << backend
                  << " nx=" << nx
                  << " ny=" << ny
                  << " steps=" << steps
                  << " retries=" << retries
                  << " total_s=" << total_s
                  << " setup_s=" << setup_s
                  << " cfl_s=" << cfl_s
                  << " workspace_s=" << workspace_s
                  << " flux_x_s=" << flux_x_s
                  << " flux_y_s=" << flux_y_s
                  << " spatial_fused_s=" << spatial_fused_s
                  << " rhs_assembly_s=" << rhs_assembly_s
                  << " state_update_s=" << state_update_s
                  << " admissibility_s=" << admissibility_s
                  << " transfer_s=" << transfer_s
                  << " unattributed_s=" << unattributed
                  << "\n";
    }
};

}  // namespace fp
