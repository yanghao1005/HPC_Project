#!/bin/bash
# filepath: d:\UDL\HPC\HPC_Project\CUDA\heat_cuda_16x16.sh

# Test configuration: 16x16 blocks (256 threads per block)
BLOCK_X=16
BLOCK_Y=16
CONFIG_NAME="block_16x16"

SIZES=(100 1000 2000)
STEPS=(100 1000 10000 100000)


RESULTS_DIR="cuda_benchmark_results/${CONFIG_NAME}"
mkdir -p $RESULTS_DIR

echo "========== Testing Block Size: ${BLOCK_X}x${BLOCK_Y} (${CONFIG_NAME}) =========="
echo "Threads per block: $((BLOCK_X * BLOCK_Y))"

PERF_LOG="${RESULTS_DIR}/performance_summary.txt"
echo "Block Configuration: ${BLOCK_X}x${BLOCK_Y} (256 threads per block)" > $PERF_LOG
echo "========================================" >> $PERF_LOG
echo "Size\tSteps\tTime(s)\tGFLOPs\tBandwidth(GB/s)" >> $PERF_LOG
echo "----------------------------------------" >> $PERF_LOG

for SIZE in "${SIZES[@]}"; do
    for STEP in "${STEPS[@]}"; do
        JOB_NAME="heat_${CONFIG_NAME}_${SIZE}_${STEP}"
        OUTPUT_FILE="${RESULTS_DIR}/${JOB_NAME}.bmp"
        LOG_FILE="${RESULTS_DIR}/${JOB_NAME}.log"
        
        echo "Running: ${JOB_NAME} (Size: ${SIZE}x${SIZE}, Steps: ${STEP})"
        
        ./heat_cuda ${SIZE} ${STEP} ${BLOCK_X} ${BLOCK_Y} ${OUTPUT_FILE} > ${LOG_FILE} 2>&1
        
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