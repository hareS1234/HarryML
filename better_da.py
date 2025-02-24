import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# Function to extract the protein name dynamically
def extract_protein_name(path_file: Path) -> str:
    match = re.match(r"density_(.+?)_\d+\.dat", path_file.name)
    
    return match.group(1) if match else "Unknown"


def extract_temperature(path_file: Path) -> int:

    match = re.match(r".*_(\d+)\.dat$", path_file.name)
    
    return int(match.group(1)) if match else None
    
    
# Function to read and process density data
def get_density(path_file: Path):
    """"""
    dense = {}  # Using a dictionary to handle unknown bin sizes
    
    with open(path_file) as file:
        for line in file:
            line = line.strip()
            
            if not line or len(line.split()) == 3 or line.startswith("#"):
                continue
            try:
                parts = line.split()
                bin_index = int(parts[0]) - 1  # Convert bin to zero-based index
                density_value = float(parts[3])  # Extract density value
    
                if bin_index not in dense:
                    dense[bin_index] = []
    
                dense[bin_index].append(density_value)
    
            except (IndexError, ValueError):
                print(f"Skipping malformed line in {path_file}: {line}")
                continue
    
    # Convert dictionary to a sorted list of lists
    max_bin = max(dense.keys(), default=0)  # Get the maximum bin index
    dense_list = [dense.get(i, []) for i in range(max_bin + 1)]  # Fill missing bins with empty lists
    density_profile = [np.mean(d) if d else 0 for d in dense_list]
    stddev_profile = [np.std(d) if d else 0 for d in dense_list]
    
    return density_profile, stddev_profile


# Function to plot density profiles
def plot_density(density_profile, protein_name):
    
    z_position = np.linspace(0.5 / len(density_profile), 1.0, len(density_profile))
    
    plt.figure(figsize=(10, 5))
    plt.plot(z_position, density_profile, linestyle="-", color="blue", label=protein_name)
    plt.ylabel("Density (g/cm³)")
    plt.title(f"Density Profile for {protein_name}")
    plt.xlabel("Z Position")
    plt.legend()
    #plt.savefig(f"{file_name}_density_profile.png")
    plt.close()


# Function to calculate dilute and dense phase concentrations
def get_dilute_and_dense_concentrations(
    density, 
    manual_override=True, 
    min_dense_threshold=0.3, 
    max_dilute_threshold=0.1
):
    """"""
    z_step = 0.02
    slope = [(density[i+1] - density[i]) / z_step for i in range(len(density) - 1)]
    
    if manual_override:
        dense_values = [value for value in density if value >= min_dense_threshold]
        dilute_values = [value for value in density if value <= max_dilute_threshold]
        dense_avg = np.mean(dense_values) if dense_values else 0
        dilute_avg = np.mean(dilute_values) if dilute_values else 0
        dense_std = np.std(dense_values) if dense_values else 0
        dilute_std = np.std(dilute_values) if dilute_values else 0
    
        return dense_avg, dilute_avg, dense_std, dilute_std
    
    return 0, 0, 0, 0  # Placeholder if manual_override is False


# Process a single density file
def process_density_file(path_file: Path):
    
    protein_name = extract_protein_name(path_file)
    temperature = extract_temperature(path_file)
    
    density, stddevs = get_density(path_file)
    
    plot_density(density, protein_name)
    
    dense, dilute, dense_stdev, dilute_stdev = get_dilute_and_dense_concentrations(
        density, manual_override=True
    )
    
    return protein_name, dense, dilute, dense_stdev, dilute_stdev, temperature


# Main function to run processing
def main():
    
    if len(sys.argv) != 2:
        print("Usage: python process_density.py <density_file>")
        return
    
    file_name = sys.argv[1]
    protein_name, dense, dilute, dense_stdev, dilute_stdev, temp = process_density_file(file_name)
    
    print(f"Protein: {protein_name}")
    print(f"Dense phase concentration: {dense:.4f} ± {dense_stdev:.4f}")
    print(f"Dilute phase concentration: {dilute:.4f} ± {dilute_stdev:.4f}")
    # Save results to a file
    with open(f"{file_name}_concentrations.txt", "w") as f:
        f.write(f"Protein: {protein_name}\n")
        f.write(f"Dense phase concentration: {dense:.4f} ± {dense_stdev:.4f}\n")
        f.write(f"Dilute phase concentration: {dilute:.4f} ± {dilute_stdev:.4f}\n")
        

if __name__ == "__main__":
    main()
