## Course Instructor (CS633 Parallel Computing)

**Prof. Preeti Malalkar**

---

## Team Members (Group 22)

- **Harsh Goel** — 251110034
- **Kishan Kumar Mishra** — 251110044
- **Aditya Pushkar** — 251110004
- **Govind Gupta** — 251110033

---

# Parallel Data Exchange using MPI

## Overview

This project implements a **parallel data exchange pattern using MPI**, where each process communicates with other processes located at fixed distances (`D1` and `D2`) in a ring-like topology. Every process exchanges data with its `D1` and `D2` neighbours over `T` iterations, and the goal is to analyze how **communication overhead** affects execution time as the system is scaled across:

- Number of processes: `P = 8, 16, 32`
- Data sizes: `M = 262144, 1048576` (doubles per message)
- 5 repeated runs per `(P, M)` configuration, to capture run-to-run timing variability

Results are visualized as **boxplots**, which let us see not just the average execution time per configuration but also its spread — a much better picture of real-world scalability than a single number per data point.

---

## Communication Pattern

Each process:
- Sends data **forward** to processes at distance `D1` and `D2`
- Receives updated data **back** from those same processes
- Handles boundary conditions when destination ranks exceed the total process count

This creates a **bidirectional communication workflow** with synchronization and message passing. To maximize parallelism, processes are grouped into **independent communication blocks** (see `Group22.pdf`, Section 1.2) so that non-interfering sends/receives can proceed concurrently, and only dependent blocks are serialized.

At the end of the run, each process computes the maximum value it received from its `D1` and `D2` neighbours (`MaxD1`, `MaxD2`), and these are reduced to rank 0 for reporting alongside the elapsed time.

---

## Project Structure

```
├── src.c                          # MPI implementation
├── plot_timing_boxplot.py         # Python/Matplotlib boxplot generator
├── timing.txt                     # Raw execution timing results
├── timing_boxplot.png             # Final boxplot output
├── job_output_files               # SLURM job output files directory
├── job_scripts                    # SLURM job scripts directory
├── Group22.pdf                    # Detailed project report
└── README.md                      # Project documentation
```

---

## ⚙️ Compilation & Execution

### 1. Compile the MPI program

```bash
mpicc -o src src.c -lm
```

### 2. Run with different process counts (locally)

```bash
mpirun -np 8  ./src
mpirun -np 16 ./src
mpirun -np 32 ./src
```

Execution times are appended to `timing.txt` in the format:

```
Run  P  M  D1  D2  T  Seed  MaxD1  MaxD2  Time
```

### 3. Running on the supercomputer (SLURM)

```bash
sbatch job_scripts/run_8.sh
sbatch job_scripts/run_16.sh
sbatch job_scripts/run_32.sh
```

Each job appends a new `<Run, P, M, D1, D2, T, Seed, MaxD1, MaxD2, Time>` row to `timing.txt` and writes `mpi_<job_id>.out` / `.err` logs to `job_output_files/`.

---

## Generating the Boxplot

**Python/Matplotlib (recommended, used for the results below):**

```bash
pip install pandas matplotlib seaborn
python3 plot_timing_boxplot.py timing.txt -o boxplot_time_vs_process.png
```

## 📈 Results

**Table 1 — Execution time (s) by run, process count, and data size**

| P | M | Run 1 | Run 2 | Run 3 | Run 4 | Run 5 |
|---|---|-------|-------|-------|-------|-------|
| 8  | 262144  | 0.1170 | 0.2666 | 0.2687 | 0.2818 | 0.2860 |
| 8  | 1048576 | 0.5966 | 0.6002 | 0.5921 | 0.5939 | 0.6138 |
| 16 | 262144  | 0.3304 | 0.4928 | 0.4835 | 0.5063 | 0.4798 |
| 16 | 1048576 | 0.9415 | 0.9213 | 0.9514 | 0.9323 | 0.9538 |
| 32 | 262144  | 0.7258 | 0.8560 | 0.8542 | 0.8532 | 0.7729 |
| 32 | 1048576 | 1.2766 | 1.3103 | 1.2914 | 1.2884 | 1.3314 |

**Summary statistics (mean ± std, seconds):**

| P | M = 262144 | M = 1048576 |
|---|---|---|
| 8  | 0.244 ± 0.071 | 0.599 ± 0.009 |
| 16 | 0.459 ± 0.072 | 0.940 ± 0.014 |
| 32 | 0.812 ± 0.060 | 1.300 ± 0.022 |

### Boxplot

![Timing variation across processes and data sizes](https://github.com/Harshgoel22/Parallel_Computing_Assignments/blob/e3ff5f677a69bbac7f580939e675fd614d4d1e47/timing_boxplot.png)

### What the boxplot shows

- **Execution time rises monotonically with P.** For both data sizes, median time roughly doubles going from `P = 8 → 16` and increases further from `16 → 32`. This is the signature of **communication overhead dominating** as more processes are added — each process participates in more `D1`/`D2` message exchanges as the system scales.
- **Larger data size (M) consistently costs more time**, since larger messages take longer to serialize and transfer over MPI, independent of process count. The `M = 1048576` boxes sit clearly above the `M = 262144` boxes at every P.
- **Spread (box height) grows with P.** At `P = 8` the five runs are tightly clustered (std ≈ 0.07 s or less), but at `P = 32` the spread is visibly wider. This reflects greater sensitivity to scheduling jitter, network contention, and synchronization delays as more processes coordinate.
- **Outliers appear at low P.** The `P = 8, M = 262144` box shows an outlier point near 0.12 s (Run 1) — a run that happened to have less contention or a faster cold start, which is exactly the kind of anomaly a boxplot is designed to surface at a glance, unlike a bar chart of averages.
- **Median vs. mean**: because the distributions are slightly right-skewed at higher P (a few slow runs pull the tail up), the boxplot's median line is a more robust "typical performance" indicator than a plain average.

---

## Key Concepts Used

1. MPI point-to-point communication (`MPI_Send`, `MPI_Recv`)
2. Independent-block scheduling to maximize concurrent communication
3. Distance-based (`D1`, `D2`) parallel data exchange pattern
4. Performance measurement across repeated runs
5. Boxplot visualization for variability analysis (Python/Matplotlib)
6. SLURM job scheduling for cluster execution

---

## Conclusions

- Execution time increases with the number of processes due to communication overhead, even though more processes parallelize the underlying computation.
- For smaller process counts (`P = 8`), execution times are low and stable across runs.
- As process count increases (`P = 16, 32`), both execution time and run-to-run variability increase, driven by MPI communication and synchronization delays.
- The algorithm scales reasonably but is sensitive to the distance-based communication pattern (`D1`, `D2`), which introduces extra synchronization cost.
- Larger data sizes add a roughly constant time overhead across all P, confirming that transfer cost and computation cost scale somewhat independently in this implementation.

---

## Possible Improvements

1. Use non-blocking MPI (`MPI_Isend`, `MPI_Irecv`) to reduce waiting and overlap communication with computation
2. Improve scalability for larger `P` (e.g., collective operations instead of manual point-to-point loops)

---

## Report

Detailed explanation of the code, communication pattern construction, and results are provided in **`Group22.pdf`**.
