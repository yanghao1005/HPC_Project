#!/bin/bash
# filepath: d:\UDL\HPC\HPC_Project\CUDA\heat_cuda_test.sh

# Single test parameters
SIZE=100
STEPS=100

# Make sure the executable is compiled
make clean
make

# Directory for results
RESULTS_DIR="results_cuda"
mkdir -p $RESULTS_DIR

# Create job name and output filename
JOB_NAME="heat_cuda_${SIZE}_${STEPS}"
OUTPUT_FILE="${RESULTS_DIR}/${JOB_NAME}.bmp"
LOG_FILE="${RESULTS_DIR}/${JOB_NAME}.log"

echo "Running: ${JOB_NAME} (Size: ${SIZE}, Steps: ${STEPS})"

# Run the simulation directly and redirect output to log file
./heat_cuda ${SIZE} ${STEPS} ${OUTPUT_FILE} > ${LOG_FILE} 2>&1

echo "Completed: ${JOB_NAME}"

echo "Test simulation completed. Results are in ${RESULTS_DIR}/"
echo "Output file: ${OUTPUT_FILE}"
echo "Log file: ${LOG_FILE}"