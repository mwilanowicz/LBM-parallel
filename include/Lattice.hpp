# pragma once
#include "Constants.hpp"
#include <array>
#include <vector>
#include <utility>

#ifdef USE_MPI
#include <mpi.h>
#endif

/*
File: Lattice.hpp
Description: Lattice class definition for Lattice Boltzmann Method simulation.
It utilizes SoA (Structure of Arrays) layout for distribution functions to optimize memory 
access.

Author: Marcel Wilanowicz
Date: 2026-05-21
*/

class Lattice {
private:
    // Distribution functions for each velocity direction (Structure of Arrays)
    std::array<std::vector<double>, LBM::Q> f_old; // Input Buffer (current state of the lattice)
    std::array<std::vector<double>, LBM::Q> f_new; // Output Buffer (result of computation for the next time step)

    int local_width; // LBM::Config::width
    int local_height; // Band height + 2 Ghost layers (rows)
    int y_start;
    int band_height; // Computational band height (without Ghost layers)

public:
    // Empty constructor: allocation is deferred to initialize()
    Lattice() : local_width(0), local_height(0), y_start(0), band_height(0) {}

    // Map (linearize) the index based on x and y coordinates of the lattice (row-major order) for quick memory access
    inline size_t map_idx(int x, int y) const {
        return static_cast<size_t>(y) * local_width + x; 
    }

    // Swap the distribution function arrays for the next iteration
    void swap () {
        std::swap(f_old, f_new); 
    }

    // Test function to check memory allocation
    size_t get_component_size (int q) const {
        return f_old[q].size();
    }

    // Setter for writing access to f_new
    inline void set_f_new(int x, int y, int q, double val) {
        size_t idx = map_idx(x, y);
        f_new[q][idx] = val;
    }

    // Getter for reading access to f_old
    inline double get_f_old(int x, int y, int q) {
        int idx = map_idx(x, y);
        return f_old[q][idx];
    }

    // Initialization of the grid
    void initialize(int width, int height, int start_y, bool is_parallel);

    // Getter for reading summed rho-value of current cell (total mass conservation) 
    // [Chen & Doolen, 1998, Eq. 3].
    inline double get_rho(int x, int y) {
        double rho = 0.0;
        size_t idx = map_idx(x, y); // Map 2D coordinates to 1D linear memory index
    
        for (int q = 0; q < LBM::Q; ++q) {
            rho += f_old[q][idx];
        }
        return rho;
    }

    // Getter for reading summed u-value of current cell (total momentum conservation)
    // [Chen & Doolen, 1998, Eq. 3].
    inline std::pair<double, double> get_u(int x, int y) {
        std::pair<double, double> u = {0.0, 0.0};
        double rho = 0.0;
        size_t idx = map_idx(x, y);
        for (int q = 0; q < LBM::Q; ++q) {
            rho += f_old[q][idx];
            u.first += f_old[q][idx] * LBM::ex[q];
            u.second += f_old[q][idx] * LBM::ey[q];
        }
        u.first = u.first / rho;
        u.second = u.second / rho;

        return u;
    }

    // The general form of the equilibrium distribution function [Chen & Doolen, 1998, Eq. 22].
    inline double get_equilibrium(int q, double rho, std::pair<double, double> u) const {
        return LBM::w[q] * rho * (
            1.0 
            + (LBM::ex[q] * u.first + LBM::ey[q] * u.second) / LBM::cs2 
            + ((LBM::ex[q] * u.first + LBM::ey[q] * u.second) 
            * (LBM::ex[q] * u.first + LBM::ey[q] * u.second)) / (2.0 * LBM::cs2 * LBM::cs2) 
            - (u.first * u.first + u.second * u.second) / (2.0 * LBM::cs2)
        );
    }

    // Performs a single discrete time step (delta t) of the LBM evolution
    // (includes local boundaries). 
    void time_step() {
        int start_y;
        int end_y;

        if (local_height > band_height) { // Parallel (MPI) mode
            start_y = 1; // Avoiding lower Ghost layer
            end_y = local_height - 1; // Avoiding upper Ghost layer (last index)
        } else { // Sequential mode
            start_y = 0;
            end_y = local_height;
        }

        // (x, y): actual position
        for (int y = start_y; y < end_y; ++y) {
            for (int x = 0; x < local_width; ++x) {
                size_t idx = map_idx(x, y); // Map actual position

                double rho = 0.0;
                double ux = 0.0;
                double uy = 0.0;

                // Macroscopic (summed up) p,u values of the fluid [Chen & Doolen, 1998, Eq. 3].
                for (int q = 0; q < LBM::Q; ++q) {
                    // Fluid total mass conservation (left equation)
                    rho += f_old[q][idx]; 

                    // Fluid total momentum conservation (right equation)
                    ux += f_old[q][idx] * LBM::ex[q];
                    uy += f_old[q][idx] * LBM::ey[q];
                }

                // We devide by rho to get the actual u value
                std::pair<double, double> u = {ux / rho, uy / rho};

                // LBGK Collision & Streaming step [Chen & Doolen, 1998, Eq. 14].
                for (int q = 0; q < LBM::Q; ++q) {

                    // Right side of the equation
                    double f_eq = get_equilibrium(q, rho, u);
                    double f_coll = f_old[q][idx] - (1.0 / LBM::Config::tau) * (f_old[q][idx] - f_eq);
                    
                    // Left side of the equation (target position): applying bounce-back boundary condition
                    int tx = (x + LBM::ex[q]);
                    int ty = (y + LBM::ey[q]);

                    int global_y = y_start + y - start_y;
                    int global_ty = global_y + LBM::ey[q];

                    // Check for hitting the boundary: Did hit == wall
                    if (tx < 0 || tx >= LBM::Config::width || global_ty < 0 || global_ty >= LBM::Config::height) {
                        
                        // If the wall is the moving lid (top)
                        if (global_ty >= LBM::Config::height) {
                            // Bounce-Back at the boundary [Mohamad, eq. 8.17]:
                            double momentum_correction =  2 * ((rho * LBM::w[q]) / LBM::cs2) * (LBM::Config::u_lid * LBM::ex[q]);
                            f_new[LBM::opposite[q]][idx] = f_coll - momentum_correction;
                        }
                        
                        // Stationary walls (bottom, left, right)
                        else {
                            f_new[LBM::opposite[q]][idx] = f_coll;
                        }
                        
                    }
                    // Didn't hit == fluid
                    else {
                        // Time Step (mapping target position)
                        f_new[q][map_idx(tx, ty)] = f_coll;
                    }

                }                
            }
        }
        swap(); // Swap the distribution function arrays for the next iteration
    }

    // Exports full velocity field (point data) for ParaView visualization
    void save_vtk(int step);

    // Exports velocity field for Ghia et al. benchmark and data analysis
    void save_csv(int step);
};
