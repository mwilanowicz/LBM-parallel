"""
File: benchmark2.py
Description: Parallel Verification Benchmark
Author: Marcel Wilanowicz
Date: 19.11.26
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import re

def validate_par(file_seq, file_par):
    """
    Loads sequential and parallel simulation data, computes node-by-node absolute errors for 
    velocity components and their magnitude, and visualizes all three.
    """
    seq_file = file_seq
    par_file = file_par

    if not seq_file or not par_file:
        print("Error: Missing required CSV files for error mapping.")
        return
    
    # print(f"Comparing: {seq_file} VS {par_file}")
    df_seq = pd.read_csv(seq_file)
    df_par = pd.read_csv(par_file)

    # Reconstruct flat CSV rows into 2D grid shapes (H x W)
    ux_seq = df_seq.pivot(index='y', columns='x', values='ux').to_numpy()
    uy_seq = df_seq.pivot(index='y', columns='x', values='uy').to_numpy()
    ux_par = df_par.pivot(index='y', columns='x', values='ux').to_numpy()
    uy_par = df_par.pivot(index='y', columns='x', values='uy').to_numpy()

    # Calculate velocity magnitude: |U| = sqrt(ux^2 + uy^2)
    u_mag_seq = np.sqrt(ux_seq**2 + uy_seq**2)
    u_mag_par = np.sqrt(ux_par**2 + uy_par**2)

    # Calculate absolute error maps (|Sequential - Parallel|)
    diff_ux = np.abs(ux_seq - ux_par)
    diff_uy = np.abs(uy_seq - uy_par)
    diff_umag = np.abs(u_mag_seq - u_mag_par)

    # Pack structures for the 3-channel iterative plotting loop
    channels = ['Streamwise Velocity ($u_x$)', 'Transverse Velocity ($u_y$)', 'Velocity Magnitude ($|U|$)']
    seq_data = [ux_seq, uy_seq, u_mag_seq]
    par_data = [ux_par, uy_par, u_mag_par]
    diff_data = [diff_ux, diff_uy, diff_umag]
    
    return channels, seq_data, par_data, diff_data
    
def plot_par(channels, seq_data, par_data, diff_data, par_file):
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))

    for i in range(3):
        im1 = axes[i, 0].imshow(seq_data[i], cmap='jet', origin='lower')
        axes[i, 0].set_title(f"Sequential: {channels[i]}")
        axes[i, 0].set_xlabel("x / W")
        axes[i, 0].set_ylabel("y / H")
        fig.colorbar(im1, ax=axes[i, 0], label="$u / u_{lid}$", shrink=0.7)

        im2 = axes[i, 1].imshow(par_data[i], cmap='jet', origin='lower')
        axes[i, 1].set_title(f"MPI Parallel: {channels[i]}")
        axes[i, 1].set_xlabel("x / W")
        axes[i, 1].set_ylabel("y / H")
        fig.colorbar(im2, ax=axes[i, 1], label="$u / u_{lid}$", shrink=0.7)

        im3 = axes[i, 2].imshow(diff_data[i], cmap='inferno', origin='lower')
        axes[i, 2].set_title(f"Absolute Error: {channels[i]}")
        axes[i, 2].set_xlabel("x / W")
        axes[i, 2].set_ylabel("y / H")
        fig.colorbar(im3, ax=axes[i, 2], label="Absolute Diff", shrink=0.7)

    match = re.search(r"mpi_(\d+)", par_file) # MPI process count extraction
    if match and int(match.group(1)) > 1:
        num_procs = match.group(1)
    else:
        num_procs = "1 (Sequential)"
        
    halo = "ON" if "halo_on" in par_file else "OFF"

    plt.suptitle(f"Parallel Verification & Error Analysis (Step: 10000, MPI Procs: {num_procs}, HALO: {halo})", fontsize=16, y=1.02)
    plt.tight_layout()
    plt.show()