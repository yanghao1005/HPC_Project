import os
import re
import csv

def parse_log_files():
    # Directory containing the log files
    log_dir = "/mnt/d/UDL/HPC/HPC_Project/Serial/results"
    
    # Dictionary to store execution times: {(size, steps): time}
    data = {}
    
    # Regular expression to extract information from filenames
    filename_pattern = re.compile(r'heat_(\d+)_(\d+)\.log')
    
    # Regular expression to extract execution time from log content
    time_pattern = re.compile(r'The Execution Time=(\d+\.\d+)s with a matrix size of (\d+)x\d+ and (\d+) steps')
    
    # Process all log files in the directory
    for filename in os.listdir(log_dir):
        if filename.endswith('.log'):
            # Extract size and steps from the filename
            match = filename_pattern.match(filename)
            if match:
                size, steps = match.groups()
                size, steps = int(size), int(steps)
                
                # Read the log file and extract the execution time
                filepath = os.path.join(log_dir, filename)
                try:
                    with open(filepath, 'r') as f:
                        content = f.read().strip()
                        time_match = time_pattern.search(content)
                        if time_match:
                            exec_time = float(time_match.group(1))
                            data[(size, steps)] = exec_time
                except Exception as e:
                    print(f"Error reading {filepath}: {e}")
    
    return data

def create_table(data):
    # Directory for saving the results
    results_dir = "/mnt/d/UDL/HPC/HPC_Project/Serial/results"
    
    # Find all unique sizes and steps
    sizes = sorted(set(size for (size, _) in data.keys()))
    steps = sorted(set(steps for (_, steps) in data.keys()))
    
    print("\n## Serial Execution Times")
    
    # Create and print a simple text table
    print("Size\\Steps", end="\t")
    for step in steps:
        print(f"{step}", end="\t")
    print()
    
    for size in sizes:
        print(f"{size}", end="\t")
        for step in steps:
            time = data.get((size, step), "-")
            print(f"{time:.6f}" if time != "-" else "-", end="\t")
        print()
    
    # Create CSV file
    csv_filename = os.path.join(results_dir, "execution_times_serial.csv")
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        header = ["Size\\Steps"] + steps
        writer.writerow(header)
        
        # Write data rows
        for size in sizes:
            row = [size]
            for step in steps:
                time = data.get((size, step), "")
                row.append(time)
            writer.writerow(row)
    
    print(f"\nTable saved to {csv_filename}")

if __name__ == "__main__":
    print("Analyzing serial heat diffusion execution times...")
    data = parse_log_files()
    
    if not data:
        print("No log files found or could not parse any execution times.")
    else:
        create_table(data)
        print("\nAnalysis complete.")