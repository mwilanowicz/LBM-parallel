#include "Constants.hpp"
#include "Lattice.hpp"
#include <iostream>
#include <iomanip>
#include <ios>
#include <mpi.h>

/*
File: main_parallel.cpp
Description: Parallel implementation of the D2Q9 LBM simulation.

Use: mpirun -np <num_procs> ./main_par

Grid decomposistion examples:

64x64 grid:
- mpirun -np 1: 64x64 grid will be divided into 1 band of 64 rows (sequential).
- mpirun -np 2: 64x64 grid will be divided into 2 bands of 32 rows each.
- mpirun -np 4: 64x64 grid will be divided into 4 bands of 16 rows each.

128x128 grid:
- mpirun -np 1: 128x128 grid will be divided into 1 band of 128 rows (sequential).
- mpirun -np 2: 128x128 grid will be divided into 2 bands of 64 rows each.
- mpirun -np 4: 128x128 grid will be divided into 4 bands of 32 rows each.

256x256 grid:
- mpirun -np 1: 256x256 grid will be divided into 1 band of 256 rows (sequential).
- mpirun -np 2: 256x256 grid will be divided into 2 bands of 128 rows each.
- mpirun -np 4: 256x256 grid will be divided into 4 bands of 64 rows each.

Author: Marcel Wilanowicz
Date: 2026-11-11
*/

int main(int argc, char** argv) {
    LBM::cfg = LBM::loadConfig("."); // Load config.csv parameters

    int rank, nprocs;
    // initialise MPI
    MPI_Init(&argc,&argv);

    // save rank of this process
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // save number of processes
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    
    // --- Domain Decomposition ---
    int rank_ny; // Local band height
    int rank_ystart; // Global offset
    int NY = LBM::cfg.height; // Global height

    // if ranks are less than remainder (non-ideal grid division)
    if(rank < NY % nprocs) { 
        rank_ny = NY / nprocs + 1; // Executes only for processes with extra rows
        rank_ystart = rank * rank_ny; // Starting point for initial processes
    }
    else { // Ideal grid division (for higher rank processes)
        rank_ny = NY / nprocs; // Standard, even number of rows (pure quotient)
        rank_ystart = NY - (nprocs - rank) * rank_ny; // Starting point from the "end" of grid
    }
    // ----------------------------

    if (rank == 0) {
        std::cout << "Processes: " << nprocs << std::endl;
        std::cout << "HALO exchange: " << (LBM::cfg.is_halo ? "ON" : "OFF") << std::endl; 
        std::cout << "Simulation width: " << LBM::cfg.width << std::endl;
        std::cout << "Simulation height: " << LBM::cfg.height << std::endl;
        std::cout << "Relaxation time (tau): " << LBM::cfg.tau << std::endl;
        std::cout << "Lid Velocity: " << LBM::cfg.u_lid << std::endl;
        std::cout << "Number of time steps: " << LBM::cfg.max_time_steps << std::endl;
        std::cout << "Reynolds Number: " << (LBM::cfg.u_lid * LBM::cfg.height) 
        / ((2 * LBM::cfg.tau - 1) / 6)  << std::endl;
    }
    
    Lattice simulation;

    // Initialization
    simulation.initialize(LBM::cfg.width, rank_ny, rank_ystart, true, LBM::cfg.is_halo); 

    // Local Initial mass
    double local_initial_mass = 0.0;
    for (int y = 1; y <= rank_ny; ++y) { // Avoids Ghost layer (y = 0)
        for (int x = 0; x < LBM::cfg.width; ++x) {
            local_initial_mass += simulation.get_rho(x, y);
        }
    }

    double global_initial_mass = 0.0;

    // Collective communication (everyone with everyone)
    MPI_Allreduce(&local_initial_mass, &global_initial_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

    if (rank == 0) {
        std::cout << "\nInitial Mass: " << std::fixed << std::setprecision(10) << global_initial_mass << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD); // Synchronization between MPI processes in a group (start)
    double start_time = MPI_Wtime(); // Simulation time measurement (start)

    for (int t = 1; t <= LBM::cfg.max_time_steps; ++t) {
        simulation.time_step();

        // Mass monitoring (after 1000 steps each)
        if (t % 1000 == 0) {
            double local_current_mass = 0.0;
            // Iterate ONLY over the rows of this concrete process
            for (int y = 1; y <= rank_ny; ++y) { // Avoids Ghost layer (y = 0)
                for (int x = 0; x < LBM::cfg.width; ++x) {
                    local_current_mass += simulation.get_rho(x, y);
                }
            }

            double global_current_mass = 0.0;

            // Collective communication (all to one): Reduction for all processes, so each knows initial mass for computing error later
            MPI_Reduce(&local_current_mass, &global_current_mass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            if (rank == 0) {
                // Computing relative error
                double rel_error = std::abs(global_current_mass - global_initial_mass) / global_initial_mass;

                std::cout << "Step: " << std::setw(6) << t 
                << " | Mass: " << std::fixed << std::setprecision(10) << global_current_mass 
                << " | Rel. Error: " << std::scientific << std::setprecision(3) << rel_error
                << "\n";
            }
            
        }
    }

    MPI_Barrier(MPI_COMM_WORLD); // Synchronization between MPI processes in a group (end)
    double end_time = MPI_Wtime(); // Simulation time measurement (end)

    double elapsed_time = end_time - start_time;

    if (rank == 0) {
        std::cout << "\nSimulation finished in: " << std::fixed << std::setprecision(4) 
        << elapsed_time << " seconds.\n";
    }

    // Export simulation data (methods use MPI_Gather internally)
    simulation.save_vtk();
    simulation.save_csv();

    // Export simulation metadata
    simulation.save_metadata(LBM::cfg.width, LBM::cfg.height, nprocs, elapsed_time);

    MPI_Finalize(); // End all processes safely and release resources

    return 0;
}
