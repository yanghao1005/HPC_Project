import os
import re
import pandas as pd

def parse_log_files():
    # Directory containing the log files
    log_dir = "/mnt/d/UDL/HPC/HPC_Project/OpenMP/results"
    
    # Dictionary to store execution times: {threads: {(size, steps): time}}
    data = {}
    
    # Regular expression to extract information from filenames
    filename_pattern = re.compile(r'heat_omp_(\d+)_(\d+)_(\d+)t\.log')
    
    # Regular expression to extract execution time from log content
    time_pattern = re.compile(r'The Execution Time = (\d+\.\d+) seconds with a matrix size of (\d+)x\d+ and (\d+) steps')
    
    # Process all log files in the directory
    for filename in os.listdir(log_dir):
        if filename.endswith('.log'):
            # Extract size, steps, and thread count from the filename
            match = filename_pattern.match(filename)
            if match:
                size, steps, threads = match.groups()
                size, steps, threads = int(size), int(steps), int(threads)
                
                # Create entry for this thread count if it doesn't exist
                if threads not in data:
                    data[threads] = {}
                
                # Read the log file and extract the execution time
                filepath = os.path.join(log_dir, filename)
                try:
                    with open(filepath, 'r') as f:
                        content = f.read().strip()
                        time_match = time_pattern.search(content)
                        if time_match:
                            exec_time = float(time_match.group(1))
                            data[threads][(size, steps)] = exec_time
                except Exception as e:
                    print(f"Error reading {filepath}: {e}")
    
    return data

def create_tables(data):
    # Directory for saving the results
    results_dir = "/mnt/d/UDL/HPC/HPC_Project/OpenMP/results"
    
    # Generate a table for each thread count
    for threads in sorted(data.keys()):
        print(f"\n## Execution Times for {threads} Thread{'s' if threads > 1 else ''}")
        
        # Find all unique sizes and steps
        sizes = sorted(set(size for (size, _) in data[threads].keys()))
        steps = sorted(set(steps for (_, steps) in data[threads].keys()))
        
        # Create a DataFrame with sizes as rows and steps as columns
        df = pd.DataFrame(index=sizes, columns=steps)
        df.index.name = "Size\\Steps"
        
        # Fill the DataFrame with execution times
        for (size, step), time in data[threads].items():
            df.loc[size, step] = time
        
        # Display the table with proper formatting
        pd.set_option('display.float_format', '{:.6f}'.format)
        print(df.to_string())
        
        # Save to CSV in the results directory
        csv_filename = os.path.join(results_dir, f"execution_times_{threads}_thread.csv")
        df.to_csv(csv_filename)
        print(f"Table saved to {csv_filename}")

if __name__ == "__main__":
    print("Analyzing OpenMP heat diffusion execution times...")
    data = parse_log_files()
    
    if not data:
        print("No log files found or could not parse any execution times.")
    else:
        create_tables(data)
        print("\nAnalysis complete.")