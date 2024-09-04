import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.optimize import curve_fit

# Function to read and process density data
def get_density(fileName):
    file = open(fileName)
    dense = [[] for _ in range(50)]

    for line in file:
        line = line.rstrip()
        if len(line) == 0 or len(line.split()) == 3:
            continue
        elif line[0] == "#":
            continue
        else:
            try:
                dense[int(line.split()[0])-1].append(float(line.split()[3]))
            except IndexError:
                print(f"Skipping line due to IndexError: {line}")
            except ValueError:
                print(f"Skipping line due to ValueError: {line}")

    file.close()

    split_value = len(dense[0]) // 2
    dense = [np.array(chunk[split_value:]) for chunk in dense if len(chunk) > 0]

    if len(dense) == 0:
        raise ValueError("No valid data found in the file")

    densityProfile = np.array([np.mean(chunk) if len(chunk) > 0 else 0 for chunk in dense])
    stddevProfile = np.array([np.std(chunk) if len(chunk) > 0 else 0 for chunk in dense])

    return densityProfile, stddevProfile

# Function to plot density profiles
def plot_density(densityProfile, file_name, dense, dilute, dense_stdev, dilute_stdev):
    NBINS = 50
    zPosition = []
    start = 0.5 / NBINS 
    for i in range(NBINS): 
        zPosition.append(start)
        start += 1.0 / NBINS

    plt.figure(figsize=(10, 5))
    plt.plot(zPosition, densityProfile, linestyle='-', color='blue', label=file_name)
    plt.ylabel("Density (g/cm3)")
    plt.title(f"Density Profile for {file_name}")
    plt.xlabel("Z Position")
    plt.legend()
    plt.text(0.02, 0.98, f"Dense: {dense:.4f} +/- {dense_stdev:.4f}\nDilute: {dilute:.4f} +/- {dilute_stdev:.4f}",
             horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)
    plt.savefig(f"{file_name}_density_profile.png")
    plt.close()

# Function to fit a super Gaussian and find interfaces
def super_gaussian_fd_method(coords, avg_profile):
    # Fit super Gaussian
    super_gaussian = lambda x, A, x0, sigma, p: A * np.exp(-2 * ((x - x0) / sigma) ** (2 * np.round(p)))
    initial_guess = [np.percentile(avg_profile, 80), np.mean(coords), np.std(coords), 2]
    popt, _ = curve_fit(super_gaussian, coords, avg_profile, p0=initial_guess, maxfev=20000)
    fine_coords = np.linspace(min(coords), max(coords), num=1000)
    fitted_profile = super_gaussian(fine_coords, *popt)

    # Compute first derivative
    first_derivative = np.gradient(fitted_profile)
    threshold = 1e-5
    binarized = [abs(x) / max(fitted_profile) <= threshold for x in first_derivative]

    # Find interfaces
    num_clusters = 1
    true_cluster_boundaries = []
    curr_cluster = [0] if binarized[0] else []
    for idx in range(1, len(binarized)):
        if binarized[idx] != binarized[idx - 1]:
            num_clusters += 1
            if binarized[idx]:
                curr_cluster.append(fine_coords[idx])
            else:
                curr_cluster.append(fine_coords[idx - 1])
                true_cluster_boundaries.append(curr_cluster)
                curr_cluster = []
    if len(curr_cluster) == 1:
        curr_cluster.append(fine_coords[-1])
        true_cluster_boundaries.append(curr_cluster)
    num_clusters += 1
    num_clusters /= 2
    if binarized[0] == True and binarized[-1] == True and num_clusters == 3:
        left_interface_coord = fine_coords[binarized.index(False)]
        right_interface_coord = fine_coords[len(binarized) - binarized[::-1].index(False) - 1]
    else:
        left_interface_coord, right_interface_coord = np.percentile(coords, 15), np.percentile(coords, 85)
    
    return left_interface_coord, right_interface_coord

# Function to calculate dilute and dense phase concentrations
def get_dilute_and_dense_concentrations(density, coords):
    left_interface_coord, right_interface_coord = super_gaussian_fd_method(coords, density)
    
    dilute_profile = np.array(density)[np.where((np.array(coords) < left_interface_coord) | (np.array(coords) > right_interface_coord))[0]]
    dense_profile = np.array(density)[np.where((np.array(coords) >= left_interface_coord) & (np.array(coords) <= right_interface_coord))[0]]
    
    mean_dense = np.mean(dense_profile)
    std_dense = np.std(dense_profile)
    mean_dilute = np.mean(dilute_profile)
    std_dilute = np.std(dilute_profile)
    
    return mean_dense, mean_dilute, std_dense, std_dilute

# Process a single density chunk file
def process_density_file(file_name):
    density, stddevs = get_density(file_name)
    NBINS = 50
    zPosition = []
    start = 0.5 / NBINS 
    for i in range(NBINS): 
        zPosition.append(start)
        start += 1.0 / NBINS

    dense, dilute, dense_stdev, dilute_stdev = get_dilute_and_dense_concentrations(density, zPosition)
    plot_density(density, file_name, dense, dilute, dense_stdev, dilute_stdev)
    return dense, dilute, dense_stdev, dilute_stdev

# Main function to run the processing
def main():
    if len(sys.argv) != 2:
        print("Usage: python process_density.py <density_chunk_file>")
        return

    file_name = sys.argv[1]
    try:
        dense, dilute, dense_stdev, dilute_stdev = process_density_file(file_name)

        print(f"Dense phase concentration: {dense} +/- {dense_stdev}")
        print(f"Dilute phase concentration: {dilute} +/- {dilute_stdev}")

        # Save results to a file
        with open(f"{file_name}_concentrations.txt", 'w') as f:
            f.write(f"Dense phase concentration: {dense} +/- {dense_stdev}\n")
            f.write(f"Dilute phase concentration: {dilute} +/- {dilute_stdev}\n")
    except ValueError as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    main()