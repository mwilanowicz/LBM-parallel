#include "Constants.hpp"
#include "Lattice.hpp"
#include <iostream>
#include <iomanip>
#include <ios>

/*
File: main_seq.cpp
Description: Sequential implementation of the D2Q9 LBM simulation
Author: Marcel Wilanowicz
Date: 2026-05-07
*/

int main() {
    LBM::Config cfg;
    std::cout << "Simulation width: " << cfg.width << std::endl;
    std::cout << "Simulation height: " << cfg.height << std::endl;
    std::cout << "Relaxation time (tau): " << cfg.tau << std::endl;
    std::cout << "Lid Velocity: " << cfg.u_lid << std::endl;
    std::cout << "Number of time steps: " << cfg.max_time_steps << std::endl;
    std::cout << "Reynolds Number: " << (cfg.u_lid * cfg.height) 
    / ((2 * cfg.tau - 1) / 6)  << std::endl;

    Lattice simulation;

    simulation.initialize(cfg.width, cfg.height, 0, false, false); // No mpi, no halo

    // Initial mass
    double initial_mass = 0.0; 
    for (int y = 0; y < cfg.height; ++y) {
        for (int x = 0; x < cfg.width; ++x) {
            initial_mass += simulation.get_rho(x, y); // 100x100 = 10000
        }
    }
    std::cout << "\nInitial Mass: " << std::fixed << std::setprecision(10) << initial_mass << std::endl;

    for (int t = 1; t <= cfg.max_time_steps; ++t) {
        simulation.time_step();

        // Mass monitoring (after 1000 steps each)
        if (t % 1000 == 0) {
            double current_mass = 0.0;
            for (int y = 0; y < cfg.height; ++y) {
                for (int x = 0; x < cfg.width; ++x) {
                    current_mass += simulation.get_rho(x, y);
                }
            }

            // Computing relative error
            double rel_error = std::abs(current_mass - initial_mass) / initial_mass;

            std::cout << "Step: " << std::setw(6) << t 
            << " | Mass: " << std::fixed << std::setprecision(10) << current_mass 
            << " | Rel. Error: " << std::scientific << std::setprecision(3) << rel_error
            << "\n";
        }
    }

    // Export simulation data
    simulation.save_vtk();
    simulation.save_csv();

    return 0;
}