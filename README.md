
---

## 👨‍🏫 Course Instructor ( CS633 Parallel Computing)

**Prof. Preeti Malalkar**

---

## 👥 Team Members (Group 22)

- **Harsh Goel** — 251110034  
- **Kishan Kumar Mishra** — 251110044  
- **Aditya Pushkar** — 251110004  
- **Govind Gupta** — 251110033  

---

# Parallel Data Exchange using MPI

## Overview
This project implements a **parallel data exchange pattern using MPI**, where each process communicates with other processes located at fixed distances (`D1` and `D2`) in a ring-like topology. The goal is to analyze how communication overhead affects execution time when scaling the number of processes.

The program measures execution time for different:
- Number of processes (`P = 8, 16, 32`)
- Data sizes (`M = 262144, 1048576`)
- Multiple runs per configuration

The results are visualized using **boxplots** to study performance variability and communication overhead.

---

## Communication Pattern

Each process:
- Sends data **forward** to processes at distance `D1` and `D2`
- Receives updated data **back** from those same processes
- Handles boundary conditions when destination ranks exceed total processes

This creates a **bidirectional communication workflow** with synchronization and message passing.

---

## Project Structure

├── src.c / mpi_program.c # MPI implementation

├── plot_script.c # Generates gnuplot script from timing data

├── boxplot.sh # Script to compile & generate boxplot

├── timing.txt # Raw execution timing results

├── plot_data.dat # Processed data for plotting

├── boxplot.gnu # Auto-generated gnuplot script

├── boxplot_time_vs_process.png # Final boxplot output

├── result.txt # SLURM job output (if run on cluster)

└── README.md # Project documentation

---

## ⚙️ Compilation & Execution

### 1. Compile MPI Program

```
mpicc src.c -o mpi_program
```
### 2. Run with different process counts
```
mpirun -np 8  ./mpi_program
mpirun -np 16 ./mpi_program
mpirun -np 32 ./mpi_program
```
Execution times will be stored in:
```
timing.txt
```

---

## 📊 Generate Boxplot
### Compile plot generator
```
gcc plot_script.c -o plot_script
```

Run
```
./plot_script
```

Output:
```
boxplot_time_vs_process.png
```

---

## 🖥 Running on Supercomputer (SLURM)

### Submit job:
```
sbatch job_scripts/boxplot.sh
```

Outputs will be stored inside:
```
job_output_files/

```
---

## 📈 Results Summary

Execution time increases as number of processes increases

Communication overhead dominates for larger P

Larger data size (M = 1048576) shows higher execution time than smaller data size

---

## 🧠 Key Concepts Used

1. MPI Point-to-Point Communication (MPI_Send, MPI_Recv)

2. Process communication scheduling

3. Parallel data exchange

4. Performance measurement

5. Boxplot visualization using gnuplot

6. SLURM job scheduling

---

## 📊 Performance Visualization

The generated plot:
```
boxplot_time_vs_process.png

```

Shows:

1. X-axis → Number of Processes (P = 8, 16, 32)

2. Y-axis → Execution Time (seconds)

3. Median, spread, and variability of execution times

---

## 📜 Report

Detailed explanation, communication pattern, results, and conclusions are provided in:

```
Group22.pdf
```

---

## 🚀 Possible Improvements

1. Use non-blocking MPI (MPI_Isend, MPI_Irecv) to reduce waiting

2. Optimize communication ordering

3. Reduce synchronization overhead

4. Improve scalability for larger P

---

🏁 Conclusion

The experiment demonstrates that while parallelization enables distributed computation, communication cost becomes the dominant factor as the number of processes increases. Efficient communication scheduling is critical for achieving scalability in MPI-based systems.
Boxplot visualizes variability and median execution time across runs
