#pragma once
#include <cstdlib>

struct Params {
    // grid
    int nx = 500;
    int ny = 197;
    int ng = 2;

    // gas
    double gamma = 1.4;         // specific heat ratio

    // physical constants
    double rho_air = 1.29;      // air density
    double rho_he  = 0.214;     // density of helium gas
    double p0      = 1.01325e5; // initial atmospheric pressure

    // domain
    double Lx = 0.225;          
    double Ly = 0.089;

    // The center position and radius of the bubble
    double xb = 0.035;          
    double yb = 0.0445;   
    double rb = 0.025;

    // shock
    double M_shock = 1.22;      // mach number
    double x_shock = 0.005;     // position

    // time
    double CFL = 0.4;
    //double T_final = 2e-4;    // Minor verification
    //double T_final = 5e-6;    // Fixed dt experiment, very short    
    double T_final = 1.1741e-3; // Final time

    // output step
    int output_every = 200;

    // weak scaling: duplicate bubble in y
    int weakN = 1;   // default: no duplication

    // Make it convenient to change the resolution
    static int getenv_int(const char* key, int fallback) {
        if (const char* v = ::getenv(key)) return ::atoi(v);
        return fallback;
    }
    static double getenv_double(const char* key, double fallback) {
        if (const char* v = ::getenv(key)) return ::atof(v);
        return fallback;
    }

    // Apply environment variable overrides
    void apply_env_overrides() {
        nx = getenv_int("HRSC_NX", nx);
        ny = getenv_int("HRSC_NY", ny);
        T_final = getenv_double("HRSC_TFINAL", T_final);
    }

    Params() {
        apply_env_overrides();
    }
};
