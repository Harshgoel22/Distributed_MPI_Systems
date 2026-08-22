#!/usr/bin/env python3
"""
plot_timing_boxplot.py

Reads a whitespace-separated timing log (e.g. timing.txt) with columns:
    Run  P  M  D1  D2  T  Seed  MaxD1  MaxD2  Time

The file may contain repeated header lines, blank lines, and comment
lines starting with '#' (e.g. "# ---- New Script Run ... ----").
Those are all handled automatically.

It produces a box-plot of `Time` showing variation across repeated
`Run`s, grouped by number of processes (P) and data size (M).

Usage:
    python3 plot_timing_boxplot.py timing.txt
    python3 plot_timing_boxplot.py timing.txt -o my_plot.png
"""

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def load_timing_data(path: str) -> pd.DataFrame:
    """Robustly load a timing.txt file that has repeated headers,
    comment lines, and blank lines scattered through it."""

    rows = []
    header = None

    with open(path, "r") as f:
        for raw_line in f:
            line = raw_line.strip()

            # Skip blank lines and comment/banner lines
            if not line or line.startswith("#"):
                continue

            fields = line.split()

            # Detect header line (first field is non-numeric, e.g. "Run")
            if fields[0].lower() == "run":
                header = fields
                continue

            # Any other line is treated as a data row
            rows.append(fields)

    if header is None:
        sys.exit("Could not find a header row (line starting with 'Run ...').")

    df = pd.DataFrame(rows, columns=header)

    # Convert everything except we keep P and M as categorical-friendly ints
    numeric_cols = ["Run", "P", "M", "D1", "D2", "T", "Seed",
                     "MaxD1", "MaxD2", "Time"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["P", "M", "Time"])
    df["P"] = df["P"].astype(int)
    df["M"] = df["M"].astype(int)

    return df


def make_boxplots(df: pd.DataFrame, out_path: str):
    sns.set_theme(style="whitegrid")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # --- Plot 1: Time vs P, one box per M (data size) ---
    sns.boxplot(
        data=df, x="P", y="Time", hue="M",
        ax=axes[0], palette="Set2"
    )
    axes[0].set_title("Timing variation vs. number of processes (P)")
    axes[0].set_xlabel("Number of Processes (P)")
    axes[0].set_ylabel("Time (s)")
    axes[0].legend(title="Data Size (M)")

    # --- Plot 2: Time vs M, one box per P (process count) ---
    sns.boxplot(
        data=df, x="M", y="Time", hue="P",
        ax=axes[1], palette="Set1"
    )
    axes[1].set_title("Timing variation vs. data size (M)")
    axes[1].set_xlabel("Data Size (M)")
    axes[1].set_ylabel("Time (s)")
    axes[1].legend(title="Processes (P)")

    fig.suptitle("Timing Variation Across Runs, Processes, and Data Sizes",
                  fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    fig.savefig(out_path, dpi=150)
    print(f"Saved plot to {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Box-plot timing.txt data")
    parser.add_argument("input", help="Path to timing.txt")
    parser.add_argument("-o", "--output", default="timing_boxplot.png",
                         help="Output image path (default: timing_boxplot.png)")
    args = parser.parse_args()

    df = load_timing_data(args.input)
    print(df.groupby(["P", "M"])["Time"].describe())

    make_boxplots(df, args.output)


if __name__ == "__main__":
    main()
