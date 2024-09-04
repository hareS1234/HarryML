import matplotlib.pyplot as plt
import numpy as np
import os

def read_density_file(file_name):
    with open(file_name, 'r') as file:
        data_lines = file.readlines()
        
        # Skip lines until we find the relevant data
        data_lines = [line for line in data_lines if "phase concentration:" in line]
        
        # Extract values
        dense_line = data_lines[0].split()
        dilute_line = data_lines[1].split()

        dense = float(dense_line[3])
        dense_stdev = float(dense_line[5])
        dilute = float(dilute_line[3])
        dilute_stdev = float(dilute_line[5])
        
    return dense, dense_stdev, dilute, dilute_stdev

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
    temperatures = [250, 270, 290, 300, 307, 314, 321, 328, 335, 342, 349, 356, 363, 370]
    dense_phase_concentrations = []
    dilute_phase_concentrations = []
    dense_stdev = []
    dilute_stdev = []

    for temp in temperatures:
        file_name = f'densities_chunked_TIA_temp{temp}.txt'
        if os.path.exists(file_name):
            dense, dense_stdev_temp, dilute, dilute_stdev_temp = read_density_file(file_name)
            dense_phase_concentrations.append(dense)
            dilute_phase_concentrations.append(dilute)
            dense_stdev.append(dense_stdev_temp)
            dilute_stdev.append(dilute_stdev_temp)
        else:
            print(f'File {file_name} not found. Skipping this temperature.')

    plot_phase_diagram(dense_phase_concentrations, dilute_phase_concentrations, temperatures, dense_stdev, dilute_stdev, title="Phase Diagram")

if __name__ == "__main__":
    main()
