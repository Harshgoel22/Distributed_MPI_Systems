#!/bin/bash
#SBATCH --job-name=iso_run
#SBATCH -N 2
#SBATCH --ntasks=96
#SBATCH --partition=cpu
#SBATCH --time=00:10:00
#SBATCH --output=job_output_files/mpi_%j.out
#SBATCH --error=job_output_files/mpi_%j.err

module load compiler/oneapi-2024/mpi

EXEC=$HOME/Group22/src

# Fixed params
d=7; T=5; seed=1000; F=2; iso=500

# Ensure output folder exists
mkdir -p job_output_files

# Initialize timing file with header
TIMING_FILE="timing.txt"
echo "Run,P,N,px,py,pz,Time(s)" > $TIMING_FILE

echo "===== MPI Runs Started ====="
echo ""

for N in 120 240
do
  for cfg in \
    "32 4 4 2" \
    "48 6 4 2" \
    "64 4 4 4" \
    "96 6 4 4"
  do
    set -- $cfg
    P=$1; px=$2; py=$3; pz=$4

    for run in {1..5}
    do
      echo "-------------------------------------------"
      echo "Running: Run=$run P=$P N=$N px=$px py=$py pz=$pz"

      # Capture output
      OUTPUT=$(mpirun -np $P $EXEC $d $P $px $py $pz $N $N $N $T $seed $F $iso)
      echo "$OUTPUT"

      # Extract last float line (the timing value)
      TIME=$(echo "$OUTPUT" | grep -E '^[0-9]+\.[0-9]+$' | tail -1)

      # Append row to timing file
      echo "$run,$P,$N,$px,$py,$pz,$TIME" >> $TIMING_FILE

      echo ""
    done
  done
done

echo "===== All MPI Runs Completed ====="
echo ""
echo "Timing data written to: $TIMING_FILE"
