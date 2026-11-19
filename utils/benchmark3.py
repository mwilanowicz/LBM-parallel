"""
File: benchmark3.py
Description: Scaling Benchmark
Author: Marcel Wilanowicz
Date: 17.11.26
"""

import matplotlib.pyplot as plt
import pandas as pd
import glob
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def load_benchmark(GRID):
    """
    Finds and merges all metadata.csv files across halo modes and MPI process counts for a given grid.
    """
    pattern = str(PROJECT_ROOT / "outputs" / "halo_*" / "mpi_*" / GRID / "metadata.csv")
    files = glob.glob(pattern)

    if not files:
        print(f"Error: No metadata.csv files found for grid={GRID}")
        return None

    dfs = [pd.read_csv(f) for f in files]
    return pd.concat(dfs, ignore_index=True)


def validate_benchmark(loader):
    """
    Aggregates execution times and computes parallel speedup and efficiency
    against a shared single-process baseline.
    """
    if loader is None:
        print("Error: No benchmark data provided.")
        return None

    df = loader if isinstance(loader, pd.DataFrame) else pd.read_csv(loader)

    # Compute average execution time for each MPI configuration
    df_avg = (
        df.groupby(["halo_mode", "nprocs"])["elapsed_time"]
        .mean()
        .reset_index()
    )

    # Shared sequential baseline (Halo off, 1 process; fallback to any 1-process run)
    baseline_rows = df_avg[(df_avg["halo_mode"] == "off") & (df_avg["nprocs"] == 1)]
    if baseline_rows.empty:
        baseline_rows = df_avg[df_avg["nprocs"] == 1]

    if baseline_rows.empty:
        print("Error: Missing nprocs=1 baseline run required to compute speedup.")
        return None

    t1 = baseline_rows["elapsed_time"].values[0]

    # Speedup = T1 / Tp
    df_avg["speedup"] = t1 / df_avg["elapsed_time"]

    # Parallel efficiency = (Speedup / P) * 100
    df_avg["efficiency"] = (df_avg["speedup"] / df_avg["nprocs"]) * 100

    return df_avg


def plot_benchmark(GRID, df_avg):
    """
    Plots execution time, speedup, and parallel efficiency from validated benchmark metrics.
    """
    if df_avg is None or df_avg.empty:
        return

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4))

    for halo in sorted(df_avg["halo_mode"].unique()):
        data = df_avg[df_avg["halo_mode"] == halo].sort_values("nprocs")

        ax1.plot(
            data["nprocs"],
            data["elapsed_time"],
            marker="o",
            linewidth=2,
            markersize=6,
            label=f"Halo {halo}"
        )

        ax2.plot(
            data["nprocs"],
            data["speedup"],
            marker="o",
            linewidth=2,
            markersize=6,
            label=f"Halo {halo}"
        )

        ax3.plot(
            data["nprocs"],
            data["efficiency"],
            marker="o",
            linewidth=2,
            markersize=6,
            label=f"Halo {halo}"
        )

    nprocs_values = sorted(df_avg["nprocs"].unique())

    # Execution time plot
    ax1.set_title("Elapsed time")
    ax1.set_xlabel("MPI processes")
    ax1.set_ylabel("Time [s]")
    ax1.set_xticks(nprocs_values)
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend()

    # Speedup plot
    ax2.plot(nprocs_values, nprocs_values, "k--", linewidth=1, label="Ideal")
    ax2.set_title("Speedup")
    ax2.set_xlabel("MPI processes")
    ax2.set_ylabel("Speedup")
    ax2.set_xticks(nprocs_values)
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend()

    # Parallel efficiency plot
    ax3.axhline(y=100, color="k", linestyle="--", linewidth=1, label="Ideal")
    ax3.set_title("Parallel efficiency")
    ax3.set_xlabel("MPI processes")
    ax3.set_ylabel("Efficiency [%]")
    ax3.set_xticks(nprocs_values)
    ax3.set_ylim([0, 140])
    ax3.grid(True, linestyle="--", alpha=0.5)
    ax3.legend()

    plt.suptitle(f"Benchmark — grid: {GRID}", fontsize=14)
    plt.tight_layout()
    plt.show()
