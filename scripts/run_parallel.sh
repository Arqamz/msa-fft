#!/usr/bin/env bash
# Parallel run script for the project using OpenMPI

# Check if the executable exists
if [ ! -f bin/msafft_parallel ]; then
  echo "Error: Parallel executable not found. Please build it first."
  exit 1
fi

# Check if NP (number of processes) is provided
if [ -z "$NP" ]; then
  echo "Number of processes (NP) not set. Using default: 4."
  NP=4  # Default to 4 processes if not provided
fi

# Set default MPI flags if not provided
if [ -z "$MPI_FLAGS" ]; then
    echo "MPI flags not set. Using default: --oversubscribe."
    MPI_FLAGS="--oversubscribe --use-hwthread-cpus"
fi
# Check if HOSTFILE is provided
if [ -z "$HOSTFILE" ]; then
    echo "Hostfile not set. Using default: hostfile."
    HOSTFILE="hostfile"  # Default to hosts.txt if not provided
fi

# Run the parallel program using mpirun with dynamic NP, MPI_FLAGS, and HOSTFILE
mpirun -np "$NP" $MPI_FLAGS --hostfile "$HOSTFILE" bin/msafft_parallel

echo "Parallel program run complete."
