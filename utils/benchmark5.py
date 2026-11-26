"""
File: benchmark4.py
Description: Strong Scaling Benchmark
Author: Marcel Wilanowicz
Date: 26.11.26
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import glob
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]

DATA_DIR = PROJECT_ROOT / "data" / "rosi_v3"


def load_benchmark(RANK, GRID):
    """
    Finds and merges all metadata.csv files across halo modes and MPI process counts for a given grid.
    """
    pattern = str(DATA_DIR / "halo_*" / f"mpi_{RANK}" / GRID / "metadata.csv")
    files = glob.glob(pattern)

    if not files:
        print(f"Error: No metadata.csv files found for grid={GRID}")
        return None

    dfs = [pd.read_csv(f) for f in files]
    df = pd.concat(dfs, ignore_index=True)
    df["grid"] = GRID
    return df


def validate_benchmark(loader):
    """
    Aggregates execution times and computes Strong scaling efficiency
    against the shared single-process reference run.
    """
    if loader is None:
        print("Error: No benchmark data provided.")
        return None

    df = loader if isinstance(loader, pd.DataFrame) else pd.read_csv(loader)

    # Compute average execution time for each grid and MPI configuration
    df_avg = (
        df.groupby(["grid", "halo_mode", "nprocs"])["elapsed_time"]
        .mean()
        .reset_index()
    )

    # Shared reference run for Strong scaling (same grid, 1 process)
    baseline_rows = df_avg[df_avg["nprocs"] == 1]
    if baseline_rows.empty:
        print("Error: Missing nprocs=1 baseline run required to compute efficiency.")
        return None

    # T1 taken per grid and halo mode, not averaged across all grids
    t1 = baseline_rows.set_index(["grid", "halo_mode"])["elapsed_time"]

    df_avg["speedup"] = [
        t1.loc[(g, h)] / t
        for g, h, t in zip(df_avg["grid"], df_avg["halo_mode"], df_avg["elapsed_time"])
    ]

    # Strong scaling efficiency = (T1 / (Tp * P)) * 100
    df_avg["efficiency"] = df_avg["speedup"] / df_avg["nprocs"] * 100

    return df_avg


def plot_benchmark(df_avg):
    """
    Plots execution time and Strong scaling efficiency across process counts.
    """
    if df_avg is None or df_avg.empty:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

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
            data["efficiency"],
            marker="o",
            linewidth=2,
            markersize=6,
            label=f"Halo {halo}"
        )

    # Execution time plot
    ax1.set_title("Elapsed time")
    ax1.set_xlabel("Ranks")
    ax1.set_ylabel("Time [s]")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend()

    # Strong scaling efficiency plot
    ax2.axhline(y=100, color="k", linestyle="--", linewidth=1, label="Ideal")
    ax2.set_title("Strong scaling efficiency")
    ax2.set_xlabel("Ranks")
    ax2.set_ylabel("Efficiency [%]")
    ax2.set_ylim([0, 140])
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend()

    plt.suptitle("Strong scaling — constant grid per ranks", fontsize=14)
    plt.tight_layout()

    # Saving (option)
    fig.savefig("strong_scaling.png", dpi=150)
    plt.close(fig)

def main():
    RANKS = [1, 4, 16, 64]
    GRID = "256x256"
    frames = []
    for rank in RANKS:
        df = load_benchmark(rank, GRID)
        if df is not None:
            frames.append(df)

    df_avg = validate_benchmark(pd.concat(frames, ignore_index=True))
    plot_benchmark(df_avg)


if __name__ == "__main__":
    main()
