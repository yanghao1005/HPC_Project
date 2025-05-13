import os
import re
import pandas as pd

def parse_log_files():
    # Directory containing the log files
    log_dir = "/mnt/d/UDL/HPC/HPC_Project/MPI/results"
    
    # Dictionary to store execution times: {processes: {(size, steps): time}}
    data = {}
    
    # Regular expression to extract information from filenames
    filename_pattern = re.compile(r'heat_mpi_(\d+)_(\d+)_(\d+)p\.log')
    
    # Regular expression to extract execution time from log content
    time_pattern = re.compile(r'The Execution Time = (\d+\.\d+) seconds with a matrix size of (\d+)x\d+ and (\d+) steps')
    
    # Process all log files in the directory
    for filename in os.listdir(log_dir):
        if filename.endswith('.log'):
            # Extract size, steps, and process count from the filename
            match = filename_pattern.match(filename)
            if match:
                size, steps, processes = match.groups()
                size, steps, processes = int(size), int(steps), int(processes)
                
                # Convert processes to string to use as dictionary key
                processes_key = str(processes)
                
                # Create entry for this process count if it doesn't exist
                if processes_key not in data:
                    data[processes_key] = {}
                
                # Read the log file and extract the execution time
                filepath = os.path.join(log_dir, filename)
                try:
                    with open(filepath, 'r') as f:
                        content = f.read().strip()
                        time_match = time_pattern.search(content)
                        if time_match:
                            exec_time = float(time_match.group(1))
                            data[processes_key][(size, steps)] = exec_time
                except Exception as e:
                    print(f"Error reading {filepath}: {e}")
    
    return data

def create_tables(data):
    # Directory for saving the results
    results_dir = "/mnt/d/UDL/HPC/HPC_Project/MPI/results"
    
    # Generate a table for each process count (sort keys numerically but keep as strings)
    for processes in sorted(data.keys(), key=int):
        print(f"\n## Execution Times for {processes} Process{'es' if int(processes) > 1 else ''}")
        
        # Find all unique sizes and steps
        sizes = sorted(set(size for (size, _) in data[processes].keys()))
        steps = sorted(set(steps for (_, steps) in data[processes].keys()))
        
        # Create a DataFrame with sizes as rows and steps as columns
        df = pd.DataFrame(index=sizes, columns=steps)
        df.index.name = "Size\\Steps"
        
        # Fill the DataFrame with execution times
        for (size, step), time in data[processes].items():
            df.loc[size, step] = time
        
        # Display the table with proper formatting
        pd.set_option('display.float_format', '{:.6f}'.format)
        print(df.to_string())
        
        # Save to CSV in the results directory
        csv_filename = os.path.join(results_dir, f"execution_times_{processes}_processes.csv")
        df.to_csv(csv_filename)
        print(f"Table saved to {csv_filename}")
    
    # Create a summary table for speedup comparison across process counts
    if len(data) > 1:
        print("\n## Generating speedup comparison tables...")
        
        # Get all size/step combinations across all processes
        all_sizes = sorted(set(size for p in data.values() for (size, _) in p.keys()))
        all_steps = sorted(set(steps for p in data.values() for (_, steps) in p.keys()))
        
        for size in all_sizes:
            for steps in all_steps:
                # Create DataFrame for this size/step combination
                process_counts = sorted(data.keys(), key=int)
                df_speedup = pd.DataFrame(index=['Execution Time (s)', 'Speedup'], 
                                         columns=[int(p) for p in process_counts])
                
                # Get baseline time (usually from lowest process count)
                baseline_proc = min(process_counts, key=int)
                if (size, steps) in data[baseline_proc]:
                    baseline_time = data[baseline_proc][(size, steps)]
                    
                    # Fill DataFrame with times and speedups
                    for proc in process_counts:
                        if (size, steps) in data[proc]:
                            time = data[proc][(size, steps)]
                            speedup = baseline_time / time
                            df_speedup[int(proc)]['Execution Time (s)'] = time
                            df_speedup[int(proc)]['Speedup'] = speedup
                    
                    # Only print/save if we have data for this combination
                    if not df_speedup.isnull().all().all():
                        # Print and save speedup table
                        print(f"\nSpeedup for Grid Size {size}x{size} with {steps} Steps:")
                        print(df_speedup.to_string())
                        
                        csv_filename = os.path.join(results_dir, f"speedup_{size}_{steps}.csv")
                        df_speedup.to_csv(csv_filename)
                        print(f"Speedup table saved to {csv_filename}")

if __name__ == "__main__":
    print("Analyzing MPI heat diffusion execution times...")
    data = parse_log_files()
    
    if not data:
        print("No log files found or could not parse any execution times.")
    else:
        create_tables(data)
        print("\nAnalysis complete.")