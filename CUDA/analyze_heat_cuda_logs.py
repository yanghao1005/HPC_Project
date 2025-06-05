import os
import re
import pandas as pd

def parse_block_log_files():
    """Parse log files and organize by block configuration"""
    # Directory containing the block result folders
    base_dir = "/mnt/d/UDL/HPC/HPC_Project/CUDA/results"
    
    # Dictionary to store data: {block_config: {(size, steps): time}}
    block_data = {}
    
    # Multiple regex patterns to try (based on your actual log format)
    filename_patterns = [
        re.compile(r'heat_cuda_(\d+)x(\d+)_(\d+)_(\d+)\.log'),
        re.compile(r'heat_block_(\d+)x(\d+)_(\d+)_(\d+)\.log'),
        re.compile(r'heat_.*_(\d+)_(\d+)\.log')
    ]
    
    time_patterns = [
        re.compile(r'The Execution Time = ([\d.]+) seconds with a matrix size of (\d+)x\d+ and (\d+) steps'),
        re.compile(r'Execution Time: ([\d.]+) seconds'),
        re.compile(r'execution time: ([\d.]+)'),
        re.compile(r'Time: ([\d.]+)s')
    ]
    
    # Process all subdirectories (block configurations)
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path):
            print(f"Processing directory: {item}")
            
            # Initialize data for this block configuration
            if item not in block_data:
                block_data[item] = {}
            
            # Process all log files in this directory
            log_files = [f for f in os.listdir(item_path) if f.endswith('.log')]
            print(f"  Found {len(log_files)} log files")
            
            for filename in log_files:
                print(f"  Processing: {filename}")
                
                # Extract info from filename first
                size, steps, block_config = None, None, None
                
                # Try to extract from filename
                for pattern in filename_patterns:
                    match = pattern.match(filename)
                    if match:
                        if len(match.groups()) == 4:  # heat_cuda_8x8_100_100.log format
                            block_x, block_y, size, steps = match.groups()
                            block_config = f"{block_x}x{block_y}"
                            size, steps = int(size), int(steps)
                        elif len(match.groups()) == 2:  # simpler format
                            size, steps = match.groups()
                            size, steps = int(size), int(steps)
                            # Extract block config from directory name
                            if 'block_' in item:
                                block_config = item.replace('block_', '')
                        break
                
                # If we couldn't parse filename, try to extract from directory name
                if not block_config and 'block_' in item:
                    block_config = item.replace('block_', '')
                
                # Read the log file and extract execution time and other info
                filepath = os.path.join(item_path, filename)
                try:
                    with open(filepath, 'r') as f:
                        content = f.read()
                    
                    print(f"    File content preview: {content[:200]}...")
                    
                    # Try to extract execution time
                    exec_time = None
                    for pattern in time_patterns:
                        time_match = pattern.search(content)
                        if time_match:
                            exec_time = float(time_match.group(1))
                            print(f"    Found execution time: {exec_time}s")
                            break
                    
                    # Try to extract size and steps from content if not from filename
                    if size is None or steps is None:
                        size_match = re.search(r'Grid size: (\d+)x(\d+)', content)
                        steps_match = re.search(r'Time steps: (\d+)', content)
                        
                        if size_match:
                            size = int(size_match.group(1))
                            print(f"    Found size from content: {size}")
                        if steps_match:
                            steps = int(steps_match.group(1))
                            print(f"    Found steps from content: {steps}")
                    
                    # Try to extract block config from content if not available
                    if not block_config:
                        block_match = re.search(r'Block size: (\d+)x(\d+)', content)
                        if block_match:
                            block_config = f"{block_match.group(1)}x{block_match.group(2)}"
                            print(f"    Found block config from content: {block_config}")
                    
                    # Store data if we have all required info
                    if exec_time is not None and size is not None and steps is not None:
                        key = (size, steps)
                        block_data[item][key] = exec_time
                        print(f"    ✓ Stored: {item} -> {key} = {exec_time}s")
                    else:
                        print(f"    ✗ Missing data: time={exec_time}, size={size}, steps={steps}")
                        
                except Exception as e:
                    print(f"    Error reading {filepath}: {e}")
    
    return block_data

def create_simple_csv_tables(block_data):
    """Create simple CSV tables for each block configuration"""
    
    if not block_data:
        print("No data found!")
        return
    
    # Directory for saving results
    results_dir = "/mnt/d/UDL/HPC/HPC_Project/CUDA/results"
    
    print(f"\nCreating CSV files...")
    
    # Create individual CSV files for each block configuration
    for block_config, data in block_data.items():
        if not data:
            print(f"No data for {block_config}")
            continue
            
        print(f"\n## Block Configuration: {block_config}")
        print(f"Found {len(data)} data points")
        
        # Find all unique sizes and steps for this block config
        sizes = sorted(set(size for (size, _) in data.keys()))
        steps = sorted(set(steps for (_, steps) in data.keys()))
        
        print(f"Sizes: {sizes}")
        print(f"Steps: {steps}")
        
        # Create DataFrame with sizes as rows and steps as columns
        df = pd.DataFrame(index=sizes, columns=steps)
        df.index.name = "Size\\Steps"
        
        # Fill the DataFrame with execution times
        for (size, step), time in data.items():
            df.loc[size, step] = time
        
        # Display the table
        print(df.to_string(float_format='{:.6f}'.format))
        
        # Save individual CSV
        csv_filename = os.path.join(results_dir, f"execution_times_{block_config}.csv")
        df.to_csv(csv_filename)
        print(f"Saved: {csv_filename}")

def create_comparison_table(block_data):
    """Create a comparison table showing all configurations side by side"""
    
    if not block_data:
        print("No data for comparison!")
        return
    
    results_dir = "/mnt/d/UDL/HPC/HPC_Project/CUDA/results"
    
    # Collect all unique (size, steps) combinations
    all_combinations = set()
    for data in block_data.values():
        all_combinations.update(data.keys())
    
    all_combinations = sorted(all_combinations)
    
    if not all_combinations:
        print("No combinations found!")
        return
    
    # Create comparison data
    comparison_data = []
    
    for size, steps in all_combinations:
        row = {'Size': size, 'Steps': steps}
        
        # Get execution time for each block configuration
        for block_config, data in block_data.items():
            if (size, steps) in data:
                row[f'{block_config}'] = data[(size, steps)]
            else:
                row[f'{block_config}'] = None
        
        # Find best configuration for this combination
        times = {k: v for k, v in row.items() if k not in ['Size', 'Steps'] and v is not None}
        if times:
            best_config = min(times.keys(), key=lambda k: times[k])
            best_time = times[best_config]
            row['Best_Config'] = best_config
            row['Best_Time'] = best_time
        
        comparison_data.append(row)
    
    # Create DataFrame
    comparison_df = pd.DataFrame(comparison_data)
    
    # Save comparison table
    comparison_filename = os.path.join(results_dir, "block_comparison.csv")
    comparison_df.to_csv(comparison_filename, index=False)
    
    print(f"\n## Block Configuration Comparison")
    print(comparison_df.to_string(index=False, float_format='{:.6f}'.format))
    print(f"\nComparison table saved: {comparison_filename}")

def print_performance_analysis(block_data):
    """Print performance analysis summary"""
    
    print(f"\n{'='*60}")
    print("CUDA BLOCK SIZE PERFORMANCE ANALYSIS")
    print(f"{'='*60}")
    
    if not block_data:
        print("No data to analyze!")
        return
    
    block_stats = {}
    
    for block_config, data in block_data.items():
        if not data:
            continue
            
        times = list(data.values())
        avg_time = sum(times) / len(times)
        min_time = min(times)
        max_time = max(times)
        
        block_stats[block_config] = {
            'avg': avg_time,
            'min': min_time,
            'max': max_time,
            'count': len(times)
        }
        
        print(f"\nBlock {block_config}:")
        print(f"  Average Time: {avg_time:.6f} seconds")
        print(f"  Best Time:    {min_time:.6f} seconds")
        print(f"  Worst Time:   {max_time:.6f} seconds")
        print(f"  Tests Run:    {len(times)}")
    
    # Find overall best configuration
    if block_stats:
        best_avg_config = min(block_stats.keys(), key=lambda k: block_stats[k]['avg'])
        best_min_config = min(block_stats.keys(), key=lambda k: block_stats[k]['min'])
        
        print(f"\n{'='*40}")
        print("SUMMARY:")
        print(f"Best Average Performance: {best_avg_config} ({block_stats[best_avg_config]['avg']:.6f}s)")
        print(f"Best Single Performance:  {best_min_config} ({block_stats[best_min_config]['min']:.6f}s)")
        print(f"{'='*40}")

if __name__ == "__main__":
    print("Analyzing CUDA heat diffusion execution times by block configuration...")
    
    # Parse all log files organized by block configuration
    block_data = parse_block_log_files()
    
    print(f"\nSummary of parsed data:")
    for block_config, data in block_data.items():
        print(f"  {block_config}: {len(data)} data points")
    
    if not any(block_data.values()):
        print("\nNo log files found or could not parse any execution times.")
        print("Please check:")
        print("1. Directory path: /mnt/d/UDL/HPC/HPC_Project/CUDA/results")
        print("2. Log file format and content")
        print("3. File permissions")
    else:
        # Create simple CSV tables (avoid Excel issues)
        create_simple_csv_tables(block_data)
        
        # Create comparison table
        create_comparison_table(block_data)
        
        # Print performance analysis
        print_performance_analysis(block_data)
        
        print("\nAnalysis complete!")
        print("\nGenerated files:")
        print("1. execution_times_*.csv - Individual CSV for each block config")
        print("2. block_comparison.csv - Side-by-side comparison")