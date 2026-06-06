#!/bin/bash
# SLURM job script for bertini2 parallel benchmark.
# Customize the #SBATCH directives and paths below for your cluster.

#SBATCH --job-name=bertini2_benchmark
#SBATCH --nodes=4                  # total number of nodes
#SBATCH --ntasks-per-node=1        # one MPI rank per node
#SBATCH --cpus-per-task=8          # CPU cores per rank (= OMP_NUM_THREADS)
#SBATCH --time=02:00:00
#SBATCH --output=benchmark_%j.log
#SBATCH --error=benchmark_%j.err

# --- Cluster-specific setup ---
# Uncomment / adjust as needed for your environment:
# module load mpi/openmpi-x86_64
# module load python/3.10
# source /path/to/your/conda/etc/profile.d/conda.sh && conda activate b2-ubuntu

# --- Paths (edit these) ---
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BERTINI2="${REPO_ROOT}/build/core/bertini2"
BENCHMARK_SCRIPT="${REPO_ROOT}/benchmark/run_benchmark.py"
INPUT="${REPO_ROOT}/benchmark/inputs/large.b2"

# Output CSV tagged with SLURM job ID
OUTPUT="${REPO_ROOT}/benchmark/results_${SLURM_JOB_ID:-local}.csv"

# --- Sweep configuration ---
# Ranks to test. Should not exceed (--nodes * --ntasks-per-node).
RANKS="1 2 4"

# Threads per rank. Should not exceed --cpus-per-task.
THREADS="1 2 4 8"

# --- Run ---
python "${BENCHMARK_SCRIPT}" \
    --bertini2 "${BERTINI2}" \
    --input    "${INPUT}" \
    --ranks    ${RANKS} \
    --threads  ${THREADS} \
    --output   "${OUTPUT}" \
    --repeats  3

echo "Benchmark complete. Results: ${OUTPUT}"
