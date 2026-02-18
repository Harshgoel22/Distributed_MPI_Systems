#!/bin/bash
#SBATCH --job-name=g22_main_run_P16
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=cpu
#SBATCH --time=00:10:00
#SBATCH --output=job_output_files/mpi_%j.out
#SBATCH --error=job_output_files/mpi_%j.err

# create folder if not exists
OUTPUT_DIR=job_output_files
mkdir -p $OUTPUT_DIR

module load compiler/oneapi-2024/mpi

# generate hostfile
scontrol show hostnames $SLURM_NODELIST > hostfile
HOSTFILE=hostfile

# ---------------- CONFIG ----------------
P=16
EXEC=$HOME/Group22/src

# Fixed parameters
D1=2
D2=4
T=10
SEED=1000

# M values
M_VALUES=(262144 1048576)

TIMING_FILE="timing.txt"

# ----------- HEADER (WRITE ONLY ONCE) ----------
if [ ! -f "$TIMING_FILE" ]; then
  printf "%-4s %-3s %-8s %-3s %-3s %-3s %-6s %-8s %-8s %-10s\n" \
  "Run" "P" "M" "D1" "D2" "T" "Seed" "MaxD1" "MaxD2" "Time" > "$TIMING_FILE"
fi

# ----------- MARK NEW SCRIPT RUN ----------
echo "" >> "$TIMING_FILE"
echo "# ---- New Script Run $(date) ----" >> "$TIMING_FILE"

# ---------------- RUN ----------------
for run in {1..5}
do
  for M in "${M_VALUES[@]}"
  do
    echo "Run=$run | P=$P | M=$M | D1=$D1 | D2=$D2 | T=$T | seed=$SEED"

    OUTPUT=$(mpirun -np $P -f $HOSTFILE $EXEC $M $D1 $D2 $T $SEED)

    LAST_LINE=$(echo "$OUTPUT" | tail -1)

    MaxD1=$(echo "$LAST_LINE" | awk '{print $1}')
    MaxD2=$(echo "$LAST_LINE" | awk '{print $2}')
    Time=$(echo  "$LAST_LINE" | awk '{print $3}')

    printf "%-4s %-3s %-8s %-3s %-3s %-3s %-6s %-8s %-8s %-10s\n" \
    "$run" "$P" "$M" "$D1" "$D2" "$T" "$SEED" "$MaxD1" "$MaxD2" "$Time" >> "$TIMING_FILE"

  done
done
