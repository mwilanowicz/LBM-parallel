#pragma once
#include <array>

/*
File: Constants.hpp
Description: D2Q9 model constants and simulation parameters for Lattice Boltzmann Method.
Author: Marcel Wilanowicz
Date: 2026-06-10
*/

namespace LBM {
    // --------------- D2Q9 model constants (Compile-time) --------------- 
    // Discrete velocity vectors (e_i) for D2Q9 model [He & Luo, 1997, Eq. 17].
    inline constexpr int Q = 9; // Number of discrete velocity and distribution function directions
    inline constexpr std::array<int, Q> ex = {0, 1, 0, -1, 0, 1, -1, -1, 1}; // x-components of velocity vectors (ex_i)
    inline constexpr std::array<int, Q> ey = {0, 0, 1, 0, -1, 1, 1, -1, -1}; // y-components of velocity vectors (ey_i)
    /* Vortices are ordered as follows:

        6  2  5
        3  0  1
        7  4  8
    
        where:

        {0, 0}, // Rest vertex
        {1, 0}, // Right vertex
        {0, 1}, // Up vertex
        {-1, 0}, // Left vertex
        {0, -1}, // Down vertex
        {1, 1}, // Up-Right vertex
        {-1, 1}, // Up-Left vertex
        {-1, -1}, // Down-Left vertex
        {1, -1} // Down-Right vertex
    */

    // Opposite directions for LBM bounce-back [Mohamad, 8.3.1 Bounce-Back]: 
    // 0->0, 1->3, 2->4, 3->1, 4->2, 5->7, 6->8, 7->5, 8->6
    inline constexpr std::array<int, Q> opposite = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    // Weight scalars (w_i) as function distribution coefficients [He & Luo, 1997, Eq. 15]
    inline constexpr std::array<double, Q> w = { 
        4.0 / 9.0, // Rest vertex
        1.0 / 9.0, // Right vertex
        1.0 / 9.0, // Up vertex
        1.0 / 9.0, // Left vertex
        1.0 / 9.0, // Down vertex
        1.0 / 36.0, // Up-Right vertex
        1.0 / 36.0, // Up-Left vertex
        1.0 / 36.0, // Down-Left vertex
        1.0 / 36.0 // Down-Right vertex
    };

    // Speed of sound (squared) in the lattice (cs^2 = c^2 / 3, where c=1) [He & Luo, 1997, R6334]
    inline constexpr double cs2 = 1.0 / 3.0; 

    // ---------------  Simulation parameters (Compile-time) ---------------
    struct Config {
        /*
        Reynolds number: Re = (u * L)/nu [White, p.27 Eq. 1.24],
        
        where:
        - u is the flow velocity (in our case u_lid)
        - L is characterstic length (in our case height)
        - nu is kinematic viscosity
        
        To get target Re = 100:
        nu = (2 * tau − 1)/6 = (2 * 0.8 - 1)/6 = 0.1
        u_lid = (Re * nu) / H = (100 * 0.1) / 100 = 0.1
        */

        // Kinematic viscosity: nu = (2*tau − 1)/6, therefore tau must be > 0.5 for stability [Chen & Doolen, 1998, p.335].
        static constexpr double tau = 0.8; 

        // Lid velocity: we use 0.1 for low mach number
        static constexpr double u_lid = 0.1;

        // Tested grid sizes: 100x100, 256x256, 512x512
        static constexpr int width = 100; // Width of the simulation grid
        static constexpr int height = 100; // Height of the simulation grid
        static constexpr int max_time_steps = 10000; // Number of simulation time steps
    };

}
