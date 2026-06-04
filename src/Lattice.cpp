#include "Constants.hpp"
#include "Lattice.hpp"
#include <fstream>
#include <string>
#include <filesystem>
#include <iostream>

/*
File: Lattice.cpp
Description: Implementation of the Lattice class methods for LBM simulation and data export.
Author: Marcel Wilanowicz
Date: 2026-05-28
*/

void Lattice::initialize(int width, int height, int start_y, bool is_parallel, bool halo_on) {
    local_width = width;
    band_height = height;
    y_start = start_y;
    use_halo_exchange = halo_on;

    if (is_parallel) {
        local_height = height + 2; // Band + 1 Ghost layer on top + 1 Ghost layer at the bottom
    } else {
        local_height = height;
    }

    // Total number of cells in the local grid
    size_t total_cells = static_cast<size_t>(local_width) * local_height; 

    // Resizing f_old and f_new
    for (int q = 0; q < LBM::Q; ++q) {
        f_old[q].resize(total_cells, 0.0);
        f_new[q].resize(total_cells, 0.0);
    }

    int start_loop_y;
    int end_loop_y;

    if (is_parallel) {
        start_loop_y = 1; // We start from 1, cause 0 is lower Ghost layer in parallel version
        end_loop_y = local_height - 1;
    } else {
        start_loop_y = 0; // We start from 0, casue no Ghost layers in sequential version
        end_loop_y = local_height;
    }

    // Row-Major Order traversal for optimization
    for (int y = start_loop_y; y < end_loop_y; ++y) { // Iterate over rows first for sped up
        for (int x = 0; x < local_width; ++x) { // Iterate over columns (elements of current row)
            size_t idx = map_idx(x, y); // Map 2D coordinates to 1D linear memory index
            for (int q = 0; q < LBM::Q; ++q) {
                // Equilibrium state at u=0, p=1 [Chen & Doolen, 1998, Eq. 22].
                f_old[q][idx] = LBM::w[q]; 
                f_new[q][idx] = LBM::w[q];
            }
        }
    }
}


void Lattice::exchange_halo() {
    #ifdef USE_MPI
    if (use_halo_exchange == false) {
        return; 
    }

    int rank = 0;
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (local_height <= band_height) return;

    int bottom_row = 1; 
    int top_row = local_height - 2; 

    int ghost_bottom = 0; 
    int ghost_top = local_height - 1; 

    int target_bottom = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int target_top = (rank == nprocs - 1) ? MPI_PROC_NULL : rank + 1;

    const int up_dirs[3] = {2, 5, 6};
    const int down_dirs[3] = {4, 7, 8};

    int up_tag = 1;
    int down_tag = 2;

    // Temporary buffer to accumulate data (prevents destroying boundary Bounce-Backs)
    std::vector<double> recv_temp(local_width, 0.0);

    // ---------- UP-Exchange ----------
    for (int i = 0; i < 3; ++i) {
        int q = up_dirs[i];
        double* send_buffer = &f_new[q][map_idx(0, ghost_top)];

        MPI_Sendrecv(
            send_buffer, local_width, MPI_DOUBLE, target_top, up_tag,
            recv_temp.data(), local_width, MPI_DOUBLE, target_bottom, up_tag,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE 
        );

        if (target_bottom != MPI_PROC_NULL) {
            for (int x = 0; x < local_width; ++x) {
                // prevents mass leakage via vertical walls
                f_new[q][map_idx(x, bottom_row)] += recv_temp[x];
            }
        }
    }

    // ---------- DOWN-Exchange ----------
    for (int i = 0; i < 3; ++i) {
        int q = down_dirs[i];
        double* send_buffer = &f_new[q][map_idx(0, ghost_bottom)];

        MPI_Sendrecv(
            send_buffer, local_width, MPI_DOUBLE, target_bottom, down_tag,
            recv_temp.data(), local_width, MPI_DOUBLE, target_top, down_tag,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE 
        );

        if (target_top != MPI_PROC_NULL) {
            for (int x = 0; x < local_width; ++x) {
                // prevents mass leakage via vertical walls
                f_new[q][map_idx(x, top_row)] += recv_temp[x]; 
            }
        }
    }
    #endif
}

// Performs a single discrete time step (delta t) of the LBM evolution
// (includes local boundaries). 
void Lattice::time_step() {

    for (int q = 0; q < LBM::Q; ++q) {
        std::fill(f_new[q].begin(), f_new[q].end(), 0.0);
    }

    int start_y;
    int end_y;
    int ghost_offset;

    if (local_height > band_height) { // Parallel (MPI) mode
        start_y = 1; // Avoiding lower Ghost layer (first index)
        end_y = local_height - 1; // Avoiding upper Ghost layer (last index)
        ghost_offset = 1;
        
    } else { // Sequential mode
        start_y = 0;
        end_y = local_height;
        ghost_offset = 0;
    }

    // (x, y): actual position
    for (int y = start_y; y < end_y; ++y) {

        int global_y = y_start + y - ghost_offset;

        if (global_y < 0 || global_y >= LBM::Config::height) {
            continue; 
        }

        for (int x = 0; x < local_width; ++x) {
            size_t idx = map_idx(x, y); 

            double rho = 0.0;
            double ux = 0.0;
            double uy = 0.0;
                        
            // Macroscopic (summed up) p, u values of the fluid [Chen & Doolen, 1998, Eq. 3].
            for (int q = 0; q < LBM::Q; ++q) {
                // Fluid total mass conservation (left equation)
                rho += f_old[q][idx]; 

                // Fluid total momentum conservation (right equation)
                ux += f_old[q][idx] * LBM::ex[q];
                uy += f_old[q][idx] * LBM::ey[q];
            }

            if (rho <= 0.0) rho = 1.0;
            std::pair<double, double> u = {ux / rho, uy / rho};

            for (int q = 0; q < LBM::Q; ++q) {
                double f_eq = get_equilibrium(q, rho, u);
                double f_coll = f_old[q][idx] - (1.0 / LBM::Config::tau) * (f_old[q][idx] - f_eq);
                
                int tx = x + LBM::ex[q];
                int ty = y + LBM::ey[q];
                int global_ty = global_y + LBM::ey[q];

                bool hit_wall = (tx < 0 || tx >= LBM::Config::width || global_ty < 0 || global_ty >= LBM::Config::height);

                // If HALO is disabled, we treat artificial MPI boundaries as physical walls (Bounce-Back)
                if (use_halo_exchange == false && local_height > band_height) {
                    if (ty < ghost_offset || ty >= local_height - ghost_offset) {
                        hit_wall = true;
                    }
                }

                // Check for hitting the boundary: Did hit == wall
                if (hit_wall) {
                    // If the wall is the moving lid (top)
                    if (global_ty >= LBM::Config::height) {
                        double momentum_correction = 2 * ((rho * LBM::w[q]) / LBM::cs2) * (LBM::Config::u_lid * LBM::ex[q]);
                        f_new[LBM::opposite[q]][idx] = f_coll - momentum_correction;
                    } 
                    
                    // Stationary walls (bottom, left, right)
                    else {
                        f_new[LBM::opposite[q]][idx] = f_coll;
                    }  
                }
                // Didn't hit == fluid
                else {
                    if (ty >= 0 && ty < local_height) {
                        f_new[q][map_idx(tx, ty)] = f_coll;
                    }
                }
            }                
        }
    }

    if (local_height > band_height && use_halo_exchange) {
        exchange_halo();
    }

    swap(); 
}


void Lattice::save_vtk(int step) {
    int rank = 0;
    
    #ifdef USE_MPI
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    #endif

    int start_y = (local_height > band_height) ? 1 : 0;

    // Allocate local buffers for computational data (excluding Ghost layers)
    size_t local_points = static_cast<size_t>(local_width) * band_height;
    std::vector<double> local_ux(local_points);
    std::vector<double> local_uy(local_points);
    std::vector<double> local_rho(local_points);

    // Extract data from local process memory
    for (int y = 0; y < band_height; ++y) {
        for (int x = 0; x < local_width; ++x) {
            int local_y = start_y + y;
            auto [ux, uy] = get_u(x, local_y);
            double rho = get_rho(x, local_y);

            size_t idx = static_cast<size_t>(y) * local_width + x;
            local_ux[idx] = ux;
            local_uy[idx] = uy;
            local_rho[idx] = rho;
        }
    }

    // Prepare global buffers on rank 0
    size_t global_points = static_cast<size_t>(LBM::Config::width) * LBM::Config::height;
    std::vector<double> global_ux;
    std::vector<double> global_uy;
    std::vector<double> global_rho;

    if (rank == 0) {
        global_ux.resize(global_points);
        global_uy.resize(global_points);
        global_rho.resize(global_points);
    }

    // MPI Communication - gather all bands into the global buffer on rank 0
    #ifdef USE_MPI
    MPI_Gather(local_ux.data(), local_points, MPI_DOUBLE, global_ux.data(), local_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(local_uy.data(), local_points, MPI_DOUBLE, global_uy.data(), local_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(local_rho.data(), local_points, MPI_DOUBLE, global_rho.data(), local_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #else
    global_ux = local_ux;
    global_uy = local_uy;
    global_rho = local_rho;
    #endif

    // File export executed exclusively by rank 0
    if (rank == 0) {
        std::filesystem::create_directory("outputs");
        std::string mpi_suffix = "";
        std::string halo_suffix = "";
        #ifdef USE_MPI
        mpi_suffix = "_mpi_" + std::to_string(nprocs);  
        halo_suffix = use_halo_exchange ? "_halo_on" : "_halo_off";      
        #endif
        std::string filename = "outputs/output_" + std::to_string(step) + mpi_suffix + halo_suffix + ".vtk";
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Unable to open a file!" << std::endl;
            return;
        }

        // VTK Legacy format header
        file << "# vtk DataFile Version 3.0\n";
        file << "LBM_LidDrivenCavity_Output\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << LBM::Config::width << " " << LBM::Config::height << " 1\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING 1 1 1\n\n";

        // Write velocity vectors (X, Y, Z components - Z is set to 0.0)
        file << "POINT_DATA " << global_points << "\n";
        file << "VECTORS velocity double\n";
        for (size_t i = 0; i < global_points; ++i) {
            file << global_ux[i] << " " << global_uy[i] << " 0.0\n";
        }

        // Write density scalars
        file << "\nSCALARS density double 1\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < global_points; ++i) {
            file << global_rho[i] << "\n";
        }

        file.close();
    }
}

void Lattice::save_csv(int step) {
    int rank = 0;
    #ifdef USE_MPI
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    #endif

    int start_y = (local_height > band_height) ? 1 : 0;

    // Every process prepares its local data (without Ghosts)
    std::vector<double> local_ux(local_width * band_height);
    std::vector<double> local_uy(local_width * band_height);

    for (int y = 0; y < band_height; ++y) {
        for (int x = 0; x < local_width; ++x) {
            int local_y = start_y + y;
            auto [ux, uy] = get_u(x, local_y);
            size_t idx = static_cast<size_t>(y) * local_width + x;
            local_ux[idx] = ux;
            local_uy[idx] = uy;
        }
    }

    // Global buffers on process 0
    std::vector<double> global_ux;
    std::vector<double> global_uy;

    if (rank == 0) {
        global_ux.resize(LBM::Config::width * LBM::Config::height);
        global_uy.resize(LBM::Config::width * LBM::Config::height);
    }

    #ifdef USE_MPI
    MPI_Gather(local_ux.data(), local_width * band_height, MPI_DOUBLE,
               global_ux.data(), local_width * band_height, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    MPI_Gather(local_uy.data(), local_width * band_height, MPI_DOUBLE,
               global_uy.data(), local_width * band_height, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
    #else
    global_ux = local_ux;
    global_uy = local_uy;
    #endif

    // Saving to a file only by process 0
    if (rank == 0) {
        std::filesystem::create_directory("outputs");
        std::string mpi_suffix = "";
        std::string halo_suffix = "";
        #ifdef USE_MPI
        mpi_suffix = "_mpi_" + std::to_string(nprocs);
        halo_suffix = use_halo_exchange ? "_halo_on" : "_halo_off";
        #endif
        std::string filename = "outputs/output_" + std::to_string(step) + mpi_suffix + halo_suffix + ".csv";
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Unable to open a file!" << std::endl;
            return;
        }

        file << "x,y,ux,uy\n";

        // Saving in Row-Major
        for (int y = 0; y < LBM::Config::height; ++y) {
            for (int x = 0; x < LBM::Config::width; ++x) {
                size_t idx = static_cast<size_t>(y) * LBM::Config::width + x;
                file << x << "," << y << "," << global_ux[idx] << "," << global_uy[idx] << "\n";
            }
        }
        file.close();
    }
}
