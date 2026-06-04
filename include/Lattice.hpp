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
Date: 2026-05-27
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
    bool use_halo_exchange;

public:
    // Empty constructor: allocation is deferred to initialize()
    Lattice() : local_width(0), local_height(0), y_start(0), band_height(0), use_halo_exchange(false) {}

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
    void initialize(int width, int height, int start_y, bool is_parallel, bool halo_on);

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

    void exchange_halo();
    
    // Performs a single discrete time step (delta t) of the LBM evolution
    void time_step();

    // Exports full velocity field (point data) for ParaView visualization
    void save_vtk(int step);

    // Exports velocity field for Ghia et al. benchmark and data analysis
    void save_csv(int step);
};
