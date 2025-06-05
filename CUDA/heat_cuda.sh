#!/bin/bash
# filepath: d:\UDL\HPC\HPC_Project\CUDA\heat_cuda.sh

# Arrays of parameters to run
SIZES=(100 1000 2000)
STEPS=(100 1000 10000 100000)

# Make sure the executable is compiled
make clean
make

# Directory for results
RESULTS_DIR="results_cuda"
mkdir -p $RESULTS_DIR

# Loop through all combinations
for SIZE in "${SIZES[@]}"; do
    for STEP in "${STEPS[@]}"; do
        # Create a unique job name and output filename
        JOB_NAME="heat_cuda_${SIZE}_${STEP}"
        OUTPUT_FILE="${RESULTS_DIR}/${JOB_NAME}.bmp"
        LOG_FILE="${RESULTS_DIR}/${JOB_NAME}.log"
        
        echo "Running: ${JOB_NAME} (Size: ${SIZE}, Steps: ${STEP})"
        
        # Run the simulation directly and redirect output to log file
        ./heat_cuda ${SIZE} ${STEP} ${OUTPUT_FILE} > ${LOG_FILE} 2>&1
        
        echo "Completed: ${JOB_NAME}"
        
        # Optional: wait a bit between runs to let the GPU cool down
        sleep 2
    done
done

echo "All simulations completed. Results are in ${RESULTS_DIR}/"