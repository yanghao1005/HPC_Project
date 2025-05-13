#!/bin/bash

# Arrays of parameters to run
SIZES=(100 1000 2000)
STEPS=(100 1000 10000 100000)
PROCESSES=(2 4 8 16 32)

# Make sure the executable is compiled
make clean
make

# Directory for results
RESULTS_DIR="results_mpi"
mkdir -p $RESULTS_DIR

# Loop through all combinations
for SIZE in "${SIZES[@]}"; do
    for STEP in "${STEPS[@]}"; do
        for PROC in "${PROCESSES[@]}"; do
            # Create a unique job name and output filename
            JOB_NAME="heat_mpi_${SIZE}_${STEP}_${PROC}p"
            OUTPUT_FILE="${RESULTS_DIR}/${JOB_NAME}.bmp"
            LOG_FILE="${RESULTS_DIR}/${JOB_NAME}.log"
            
            # Create a temporary job script for this specific run
            cat > temp_job_${JOB_NAME}.sh << EOF
#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -pe mpich ${PROC}
#$ -cwd
#$ -N ${JOB_NAME}
#$ -m
#$ -M hy2@alumnes.udl.cat
#$ -o ${LOG_FILE}
#$ -j y

# Set up MPI machine file
MPICH_MACHINES=\$TMPDIR/mpich_machines
cat \$PE_HOSTFILE | awk '{print \$1":"\$2}' > \$MPICH_MACHINES

# Run the heat diffusion simulation with specific parameters
mpiexec -f \$MPICH_MACHINES -n ${PROC} ./heat_mpi ${SIZE} ${STEP} ${OUTPUT_FILE}

# Clean up
rm -rf \$MPICH_MACHINES
EOF

            # Submit the job to the SGE queue
            echo "Submitting job: ${JOB_NAME} (Size: ${SIZE}, Steps: ${STEP}, Processes: ${PROC})"
            qsub temp_job_${JOB_NAME}.sh
            
            # Optional: wait a bit between submissions to avoid overwhelming the scheduler
            sleep 1
        done
    done
done

# Clean up temporary job scripts
rm -f temp_job_heat_mpi_*.sh

echo "All jobs submitted. Results will be in ${RESULTS_DIR}/"