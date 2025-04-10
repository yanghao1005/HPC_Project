#!/bin/bash

# Arrays of parameters to run
SIZES=(100 1000 2000)
STEPS=(100 1000 10000 100000)
THREADS=(2 4 8)

# Make sure the executable is compiled
make clean
make

# Directory for results
RESULTS_DIR="results_omp"
mkdir -p $RESULTS_DIR

# Loop through all combinations
for SIZE in "${SIZES[@]}"; do
    for STEP in "${STEPS[@]}"; do
        for THREAD in "${THREADS[@]}"; do
            # Create a unique job name and output filename
            JOB_NAME="heat_omp_${SIZE}_${STEP}_${THREAD}t"
            OUTPUT_FILE="${RESULTS_DIR}/${JOB_NAME}.bmp"
            LOG_FILE="${RESULTS_DIR}/${JOB_NAME}.log"
            
            # Create a temporary job script for this specific run
            cat > temp_job_${JOB_NAME}.sh << EOF
#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -pe smp 4
#$ -cwd
#$ -N ${JOB_NAME}
#$ -m e
#$ -M hy2@alumnes.udl.cat
#$ -o ${LOG_FILE}
#$ -j y
#$ -v OMP_NUM_THREADS=${THREAD}

# Run the heat diffusion simulation with specific parameters
./heat_omp ${SIZE} ${STEP} ${OUTPUT_FILE}
EOF

            # Submit the job to the SGE queue
            echo "Submitting job: ${JOB_NAME} (Size: ${SIZE}, Steps: ${STEP}, Threads: ${THREAD})"
            qsub temp_job_${JOB_NAME}.sh
            
            # Optional: wait a bit between submissions to avoid overwhelming the scheduler
            sleep 1
        done
    done
done

# Clean up temporary job scripts (uncomment if desired)
rm -f temp_job_heat_omp_*.sh

echo "All jobs submitted. Results will be in ${RESULTS_DIR}/"