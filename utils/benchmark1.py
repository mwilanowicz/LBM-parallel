"""
File: benchmark1.py
Description: Physical Validation with Ghia et. al. Benchmark
Author: Marcel Wilanowicz
Date: 19.11.26
"""

from pathlib import Path
import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import re

# Anchors to parallel-lbm/ regardless of where the notebook is executed
PROJECT_ROOT = Path(__file__).resolve().parent.parent

def load_data(W, H, nprocs=1, halo_on=False):
    """
    Finds CSV output file (compressed or uncompressed) for given grid size, 
    process count, and halo mode.
    """
    halo_dir = "halo_on" if halo_on else "halo_off"
    mpi_dir = f"mpi_{nprocs}"
    grid = f"{W}x{H}"
    
    target_dir = PROJECT_ROOT / "outputs" / halo_dir / mpi_dir / grid
    
    if not target_dir.exists():
        print(f"Error: Directory '{target_dir}' not found")
        return None
    
    matches = glob.glob(str(target_dir / "output*.csv*"))
    if not matches:
        print(f"Error: No output CSV found in '{target_dir}'")
        return None
        
    return matches[0]

def validate_ghia(file, W, H, U_LID):
    
    target_file = file

    if not target_file:
        print("Error: Not found any files in 'outputs' directory.")
        return
    
    # print(f"Loading data from: {target_file}")
    data = pd.read_csv(target_file)

    # Centerline grid coordinates
    mid_x = (W - 1) // 2
    mid_y = (H - 1) // 2

    # Extracting vertical profile for ux comparison
    v_line = data[data['x'] == mid_x].sort_values('y')

    # Extracting horizontal profile for uy comparison
    h_line = data[data['y'] == mid_y].sort_values('x')

    # Normalizing the data (velocity profiles and dimension of the grid): 
    # we want dimensionless values for correct benchmark [Chen & Doolen, 1998, p.342].
    ux_norm = v_line['ux'] / U_LID
    y_norm = v_line['y'] / (H - 1)
    uy_norm = h_line['uy'] / U_LID
    x_norm = h_line['x'] / (W - 1)
    
    return ux_norm, y_norm, uy_norm, x_norm

def plot_ghia(ux_norm, y_norm, uy_norm, x_norm, par_file):  
      
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    # ------------ 1. Plotting vertical profile: ux(y) ------------

    # Simulation data
    ax1.plot(ux_norm, y_norm, 'b-', label='LBM Simulation Data', linewidth=2)

    # Ghia data
    # Results for u-velocity along Vertical Line through Geometric Center of Cavity 
    # [Ghia et al. p.398, Table I]
    ghia_y =  [0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 
            0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1.0000]

    ghia_ux = [0.00000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, 
            -0.20581, -0.14098, -0.01447, 0.17527, 0.55892, 0.61756, 0.68439, 0.75837, 1.00000]

    # Results for u-velocity along Vertical Line through Geometric Center of Cavity 
    #[Ghia et al. p.399, Table II]
    ghia_x = [0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5000, 
            0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1.0000]

    ghia_uy = [0.00000, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507, 0.17527, 
            0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0.00000]    
    
    ax1.plot(ghia_ux, ghia_y, 'ro', label='Ghia et al. (1982) Reference', markersize=6)

    ax1.set_title(f"Vertical Profile $u_x$ (Step: 10000)")
    ax1.set_xlabel("$u_x / u_{lid}$")
    ax1.set_ylabel("$y / H$")
    ax1.set_xlim([-1.0, 1.0])
    ax1.set_ylim([-0.1, 1.1])
    ax1.grid(True, which='both', linestyle='--', alpha=0.5)
    ax1.legend()

    # ------------ 2. Plotting horizontal profile: uy(x) ------------

    # LBM data
    ax2.plot(uy_norm, x_norm, 'b-', label='LBM Simulation Data', linewidth=2)

    # Ghia data
    ax2.plot(ghia_uy, ghia_x, 'ro', label='Ghia et al. (1982) Reference', markersize=6)

    ax2.set_title(f"Horizontal Profile $u_y$ (Step: 10000)")
    ax2.set_xlabel("$u_y / u_{lid}$")
    ax2.set_ylabel("$x / W$")
    ax2.set_xlim([-1.0, 1.0])
    ax2.set_ylim([-0.1, 1.1])
    ax2.grid(True, which='both', linestyle='--', alpha=0.5)
    ax2.legend()

    # Extract MPI process count directly from the par_file path string
    match = re.search(r"mpi_(\d+)", par_file)
    if match and int(match.group(1)) > 1:
        num_procs = match.group(1)
    else:
        num_procs = "1 (Sequential)"

    plt.suptitle(f"LBM Velocity Profiles Validation vs Ghia et al. Reference (Step: 10000, MPI Procs: {num_procs})", fontsize=16, y=1.02)
    plt.tight_layout()
    plt.show()