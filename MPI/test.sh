#!/bin/bash

# Set the number of OpenMP threads per MPI process
export OMP_NUM_THREADS=2

# Compile the program
mpicc -fopenmp heat_mpi.c -o heat_mpi -lm

# Define the problem size and number of steps
# Use smaller values for testing locally
SIZE=100
STEPS=100
OUTPUT="output_hybrid.bmp"

# Run the program with 2 MPI processes
mpirun -np 2 ./heat_mpi $SIZE $STEPS $OUTPUT

echo "Simulation completed. Output saved to $OUTPUT"