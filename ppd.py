import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import re
from scipy.optimize import curve_fit


def read_density_file(file_name):
    
    with open(file_name, 'r') as file:
        data_lines = file.readlines()
        
        print(f"File {file_name} read successfully with {len(data_lines)} lines.")
        print("Original data lines:")
    
        for line in data_lines:
            print(line.strip())
        
        # Skip lines until we find the relevant data
        data_lines = [line for line in data_lines if "phase concentration:" in line.lower()]
        
        print(f"Filtered data lines: {len(data_lines)}")
        for line in data_lines:
            print(line.strip())
        
        if len(data_lines) < 2:
            print(f"Error: Expected at least 2 data lines in {file_name}, but got {len(data_lines)}.")
            return None, None, None, None
        
        # Extract values
        dense_line = data_lines[0].split()
        dilute_line = data_lines[1].split()

        dense = float(dense_line[3])
        dense_stdev = float(dense_line[5])
        dilute = float(dilute_line[3])
        dilute_stdev = float(dilute_line[5])
        
    return dense, dense_stdev, dilute, dilute_stdev


def super_gaussian_fd_method(coords, avg_profile, temp, percentile_low=15, percentile_high=85, min_dilute_size_multiplier=0.3):
    
    super_gaussian = lambda x, A, x0, sigma, p: A*np.exp(-2*((x-x0)/sigma)**(2*np.round(p)))
    initial_guess = [np.percentile(avg_profile, 80), np.mean(coords), np.std(coords), 2]
    popt, pcov = curve_fit(super_gaussian, coords, avg_profile, p0=initial_guess, maxfev=10000)
    A, x0, sigma, p = popt
    P = 2*np.round(p)
    fine_coords = np.linspace(min(coords), max(coords), num=1000)
    fitted_profile = super_gaussian(fine_coords, *popt)
    first_derivative = np.gradient(fitted_profile)
    threshold = 1e-5
    binarized = [abs(x)/max(fitted_profile) <= threshold for x in first_derivative]
    num_clusters = 1
    true_cluster_boundaries = []
    curr_cluster = [0] if binarized[0] else []
    
    for idx in range(1, len(binarized)):
        if binarized[idx] != binarized[idx-1]:
            num_clusters += 1
            if binarized[idx]:
                curr_cluster.append(fine_coords[idx])
            else:
                curr_cluster.append(fine_coords[idx-1])
                true_cluster_boundaries.append(curr_cluster)
                curr_cluster = []
    
    if len(curr_cluster) == 1:
        curr_cluster.append(fine_coords[-1])
        true_cluster_boundaries.append(curr_cluster)
    
    num_clusters += 1
    num_clusters /= 2
    
    if binarized[0]==True and binarized[-1]==True and num_clusters == 3:
        left_interface_coord = fine_coords[binarized.index(False)]
        right_interface_coord = fine_coords[len(binarized) - binarized[::-1].index(False) - 1]
    
        if (left_interface_coord + (max(fine_coords)-right_interface_coord)) <= min_dilute_size_multiplier * right_interface_coord - left_interface_coord:
            left_interface_coord = np.percentile(coords, percentile_low)
            right_interface_coord = np.percentile(coords, percentile_high)
    else:
        left_interface_coord = np.percentile(coords, percentile_low)
        right_interface_coord = np.percentile(coords, percentile_high)
    
    return left_interface_coord, right_interface_coord, sigma


def plot_phase_diagram(dense, dilute, temps, dense_stdev, dilute_stdev, title="Phase Diagram"):
    Tc, d = find_crit_temp(dense, dilute, temps)
    Pc, A = find_crit_dense(dense, dilute, temps, Tc)

    temp_lin = np.linspace(np.min(temps) - 1, Tc, 1000)
    p_high = []
    p_low = []
    
    for i in range(len(temp_lin)):
        a = np.array([[1, -1], [1, 1]])
        b = np.array([(d * (1 - (temp_lin[i] / Tc))) ** (1 / 3.06), 2 * Pc + 2 * A * (temp_lin[i] - Tc)])
        pH, pL = np.linalg.solve(a, b)
        p_high.append(pH)
        p_low.append(pL)

    dense_plus_dilute = [i for i in dense]
    temps_dense_plus_dilute = [i for i in temps]
    
    for i in range(len(dilute) - 1, -1, -1):
        dense_plus_dilute.append(dilute[i])
        temps_dense_plus_dilute.append(temps[i])

    plt.figure(figsize=(8, 8))
    plt.errorbar(dense, temps, xerr=dense_stdev, fmt='o', color="blue", label="Dense Phase Concentration")
    plt.errorbar(dilute, temps, xerr=dilute_stdev, fmt='o', color="red", label="Dilute Phase Concentration")
    plt.scatter(Pc, Tc, marker='o', facecolor="none", edgecolor="black", label="Critical Point")
    plt.plot(p_high, temp_lin, linestyle="-", color="blue")
    plt.plot(p_low, temp_lin, linestyle="-", color="red")
    plt.ylabel("Temperature (K)", size=11)
    plt.title(title, size=11)
    plt.xlabel("Density (g/cm3)", size=11)
    plt.legend(fontsize=11)
    plt.xticks(size=10)
    plt.yticks(size=10)
    plt.savefig("PhaseDiagram.pdf", dpi=300)
    
    return Tc


def find_crit_temp(dense, dilute, temps):
    dense_sub_dilute = [(i - j) ** 3.06 for i, j in zip(dense, dilute)]
    slope, intercept = np.polyfit(temps, dense_sub_dilute, 1)

    Tc = (-1) * (slope / intercept) ** (-1)

    return Tc, intercept


def find_crit_dense(dense, dilute, temps, Tc):
    dense_sub_dilute = [(i + j) for i, j in zip(dense, dilute)]
    slope, intercept = np.polyfit(temps, dense_sub_dilute, 1)

    Pc = (intercept + slope * Tc) / 2

    return Pc, slope / 2


def main():
    temperature_files = [f for f in os.listdir('.') if re.match(r'.*temp\d+\.txt', f)]
    temperatures = sorted([int(re.findall(r'\d+', f)[0]) for f in temperature_files])
    
    dense_phase_concentrations = []
    dilute_phase_concentrations = []
    dense_stdev = []
    dilute_stdev = []

    for temp in temperatures:
        file_name = f'densities_chunked_temp{temp}.txt'
        if os.path.exists(file_name):
            dense, dense_stdev_temp, dilute, dilute_stdev_temp = read_density_file(file_name)
            
            if dense is not None and dilute is not None:
                dense_phase_concentrations.append(dense)
                dilute_phase_concentrations.append(dilute)
                dense_stdev.append(dense_stdev_temp)
                dilute_stdev.append(dilute_stdev_temp)
        else:
            print(f'File {file_name} not found. Skipping this temperature.')

    print("Temperatures:", temperatures)
    print("Dense phase concentrations:", dense_phase_concentrations)
    print("Dilute phase concentrations:", dilute_phase_concentrations)
    print("Dense std dev:", dense_stdev)
    print("Dilute std dev:", dilute_stdev)

    if dense_phase_concentrations and dilute_phase_concentrations:
        Tc = plot_phase_diagram(dense_phase_concentrations, dilute_phase_concentrations, temperatures, dense_stdev, dilute_stdev, title="Phase Diagram")
        
        # Write the critical solution temperature to a text file
        with open("critical_solution_temperature.txt", "w") as file:
            file.write(f"Critical Solution Temperature: {Tc:.2f} K\n")
            print(f"Critical Solution Temperature: {Tc:.2f} K written to critical_solution_temperature.txt")


if __name__ == "__main__":
    main()