#!/bin/bash
# filepath: d:\UDL\HPC\HPC_Project\CUDA\heat_cuda_8x8.sh

# Test configuration: 8x8 blocks (64 threads per block)
BLOCK_X=8
BLOCK_Y=8
CONFIG_NAME="block_8x8"

# Arrays of parameters to run
SIZES=(100 1000 2000)
STEPS=(100 1000 10000 100000)



# Directory for this configuration
RESULTS_DIR="cuda_benchmark_results/${CONFIG_NAME}"
mkdir -p $RESULTS_DIR

echo "========== Testing Block Size: ${BLOCK_X}x${BLOCK_Y} (${CONFIG_NAME}) =========="
echo "Threads per block: $((BLOCK_X * BLOCK_Y))"
echo "Results directory: ${RESULTS_DIR}"

# Create performance log header
PERF_LOG="${RESULTS_DIR}/performance_summary.txt"
echo "Block Configuration: ${BLOCK_X}x${BLOCK_Y} (64 threads per block)" > $PERF_LOG
echo "========================================" >> $PERF_LOG
echo "Size\tSteps\tTime(s)\tGFLOPs\tBandwidth(GB/s)" >> $PERF_LOG
echo "----------------------------------------" >> $PERF_LOG

# Loop through all combinations
for SIZE in "${SIZES[@]}"; do
    for STEP in "${STEPS[@]}"; do
        JOB_NAME="heat_${CONFIG_NAME}_${SIZE}_${STEP}"
        OUTPUT_FILE="${RESULTS_DIR}/${JOB_NAME}.bmp"
        LOG_FILE="${RESULTS_DIR}/${JOB_NAME}.log"
        
        echo "Running: ${JOB_NAME} (Size: ${SIZE}x${SIZE}, Steps: ${STEP})"
        
        # Run the simulation
        ./heat_cuda ${SIZE} ${STEP} ${BLOCK_X} ${BLOCK_Y} ${OUTPUT_FILE} > ${LOG_FILE} 2>&1
        
        # Extract performance data and append to summary
        if [ -f ${LOG_FILE} ]; then
            TIME=$(grep "Execution Time:" ${LOG_FILE} | awk '{print $3}')
            GFLOPS=$(grep "Performance:" ${LOG_FILE} | awk '{print $2}')
            BANDWIDTH=$(grep "Memory Bandwidth:" ${LOG_FILE} | awk '{print $3}')
            echo -e "${SIZE}\t${STEP}\t${TIME}\t${GFLOPS}\t${BANDWIDTH}" >> $PERF_LOG
        fi
        
        echo "Completed: ${JOB_NAME}"
        sleep 1
    done
done

echo "Block ${CONFIG_NAME} testing completed!"