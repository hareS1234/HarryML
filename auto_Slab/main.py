
import numpy as np
from pathlib import Path
import os
import sys



def prepare_NPT(file_path, N=64): 
    """
    Prepares a simulation box with N (default 64) copies for compression.

    Parameters:
    file_path (str): Name of the protein sequence file.
    N (int, optional): Number of copies. Defaults to 64.

    Returns:
    parsed sequence name
    the sequence length
    the number copies
    """
    # Initialize an empty list to store the letter sequences from the input file.
    letter_sequences = []
    # Open the input file in read mode and iterate over each line.
    with open(file_path, 'r') as file:
        for line in file:
            # Skip lines that start with '#' or '>' (comments and sequence headers).
            if not (line.startswith('#') or line.startswith('>')):
                # Append the line to the list of letter sequences.
                letter_sequences.append(line.strip())
    # Join all the letter sequences into a single string.
    sequence = ''.join(letter_sequences)
    # Define a dictionary that maps each amino acid to its corresponding numerical ID.
    seq2num = {
        'M':  1,   'G':  2,   'K':  3,   'T':  4,   'R':  5,   'A':  6,
        'D':  7,   'E':  8,   'Y':  9,   'V': 10,   'L': 11,   'Q': 12,
        'W': 13,   'F': 14,   'S': 15,   'H': 16,   'N': 17,   'P': 18,
        'C': 19,   'I': 20
    }
    # Use a list comprehension to map each amino acid in the sequence to its numerical ID.
    mapped_sequence = [seq2num.get(aa, 0) for aa in sequence]
    
    # Calculate the total number of atoms and bonds based on the length of the mapped sequence.
    numAtoms = len(mapped_sequence)
    numBonds = numAtoms - 1
    # Calculate the overall maximum size of the simulation box based on the number of atoms.
    overallMax = numAtoms * 3.81 * 0.6
    # Extract the name of the input file sequence from its path.
    seq_name = Path(file_path).stem
    # Open a configuration file in write mode and write the necessary information to it.
    config = open('config_NPT_' + seq_name + '.dat', "w")
    # Write the header of the LAMMPS script, including the number of atoms and bonds.
    config.write("LAMMPS Script for IDP " + seq_name + " \n")
    config.write("\n")
    config.write("{} atoms\n".format(numAtoms))
    config.write("{} bonds\n".format(numBonds))
    # Write the types of atom and bond, as well as the box dimensions.
    config.write("20 atom types\n")
    config.write("1 bond types\n")
    # config.write("1 angle types\n")
    config.write("\n")
    config.write("{} {} xlo xhi\n".format(-overallMax, overallMax))
    config.write("{} {} ylo yhi\n".format(-overallMax, overallMax))
    config.write("{} {} zlo zhi\n".format(-overallMax, overallMax))
    # Write the information about atoms in the simulation.
    config.write("\n\nAtoms\n\n")
    k = 0
    for i in range(numAtoms):
        k += 1
        # Write the atom type and position coordinates to the configuration file.
        config.write(f"{k} {1} {mapped_sequence[i]} {0} {i * 3.81 - overallMax * 0.5:.3f} {0:.3f} {0:.3f}\n")
    # Write the information about bonds in the simulation.
    config.write("\n\nBonds\n\n")
    k = 0
    for i in range(numAtoms - 1):
        k += 1
        # Write the bond type and atom indices to the configuration file.
        config.write(f"{k} {1} {k} {k + 1}\n")
    
    # Create a new directory to store the simulation configuration files.
    sim_dir = 'NPT_' + seq_name
    os.system('rm -r ' + sim_dir)
    os.system('mkdir ' + sim_dir)
    os.system('mv ' + 'config_NPT_' + seq_name + '.dat' + ' ' + sim_dir + '/')
    os.system('cp potential.param' + ' ' + sim_dir + '/')

    # Return the name, sequence length, and number of copies
    return seq_name, len(mapped_sequence), N


def compress_NPT(seq_name, Lammps, T = 100, n = 30000):
    """
    Simulates NPT (Nose-Hoover-Parrinello-Teter) compression of a protein sequence.

    Parameters:
    seq_name (str): Name of the protein sequence.
    LAMMPS (str): Path to the LAMMPS executable.
    T (float, optional): Temperature in Kelvin. Defaults to 100 K.
    n (int, optional): Number of compressed steps. Defaults to 30000.

    Returns:
    None
    """

    # Create a new LAMMPS input file for NPT compression
    lmp_file = 'NPT_compress_' + seq_name + '.in'
    
    # Copy the template input file and overwrite it with the new contents
    os.system('cp NPT_compress_template.in ' + lmp_file)

    # Open the new LAMMPS input file in read-write mode
    with open(lmp_file, 'r+') as f:
        # Read the entire content of the file
        content = f.read()
        # Move the cursor back to the beginning of the file for rewriting
        f.seek(0, 0)
        # Write new variables and settings to the file
        f.write(f'variable myTemp equal {T:.0f}\n')  # Set temperature to T K
        f.write(f'variable compressedStep equal {n:.0f}\n')  # Set number of compressed steps
        f.write(f'variable protein string ' + seq_name + '\n')  # Define the protein name
        
        # Append the original content to the rewritten file
        f.write('\n' + content)

    # Create a directory for the NPT simulation results
    sim_dir = 'NPT_' + seq_name
    
    # Move the LAMMPS input file into the new directory
    os.system('mv ' + lmp_file + ' ' + sim_dir + '/')

    # Change into the new directory for running the simulation
    original_cwd = os.getcwd()
    os.chdir(sim_dir + '/')
    
    # Run the NPT compression simulation using LAMMPS
    os.system('srun ' + Lammps + ' -in ' + lmp_file + '> tmp.log')

    # The log file name for the simulation results
    log_file = seq_name + '_' + str(T) + '_compressed.log'

    # Extract the density from the log file
    density = 0
    with open(log_file, 'r') as f:
        for line in f:
            data = line.split()
            if len(data) == 5 and data[0] == str(n):
                # print('found line!')
                density = float(data[-1])

    # Print the simulation results
    print(f'NPT compression at {T:.0f}K is done! Density is {density:.5f} g/cm^3')

    # Change back to the original directory
    os.chdir(original_cwd)

    
def relax_NPT(seq_name, Lammps, T, n1, n2=150000):
    """
    Simulates a NPT relaxation at zero pressure simulation using LAMMPS.
    
    Parameters:
    seq_name (str): Name of the protein sequence
    LAMMPS (str): Path to the LAMMPS executable
    T (float): Temperature for the simulation in Kelvin
    n1 (int): Number of steps for compressed phase
    n2 (int, optional): Number of steps for relaxation phase. Defaults to 150000.
    
    Returns:
    density (float): Average tail density obtained from the log file
    """

    # Create a new LAMMPS input file for NPT relaxation simulation
    lmp_file = 'NPT_relax_' + seq_name + '.in'
    
    # Copy the template LAMMPS input file to the new one
    os.system('cp NPT_relax_template.in ' + lmp_file)
    
    # Open the new LAMMPS input file for modification
    with open(lmp_file, 'r+') as f:
        # Read the content of the file
        content = f.read()
        
        # Move the cursor to the beginning of the file
        f.seek(0, 0)
        
        # Write the modified variables and parameters to the file
        f.write(f'variable myTemp equal {T:.0f}\n')
        f.write(f'variable compressedStep equal {n1:.0f}\n')
        f.write(f'variable relaxStep equal {n2:.0f}\n')
        f.write(f'variable protein string ' + seq_name + '\n')
        
        # Append the original content to the file
        f.write('\n' + content)
    
    # directory for the NPT simulation output
    sim_dir = 'NPT_' + seq_name
    
    # Move the modified LAMMPS input file to the directory
    os.system('mv ' + lmp_file + ' ' + sim_dir + '/')
    
    # Change into the directory
    original_cwd = os.getcwd()
    os.chdir(sim_dir + '/')
    
    # Run the LAMMPS simulation using the modified input file
    os.system('srun ' + Lammps + ' -in ' + lmp_file + '> tmp.log')
    
    # Get the name of the log file based on the sequence name and temperature
    log_file = seq_name + '_' + f'{T:.0f}' + '.log'
    
    # Initialize an empty list to store the density values from the log file
    density_list = []
    
    # Open the log file for reading
    with open(log_file, 'r') as f:
        # Iterate over each line in the log file
        for line in f:
            # Check if the line contains the "Step" keyword
            if line.startswith('   Step'):
                # Get the density value from the next line and append it to the list
                for i in range(int(n2/100)):
                    data = next(f).split()
                    density_list.append(float(data[-1]))
    
    # Calculate the average density from the last 100 lines of the log file
    density = np.mean(density_list[-100:-1])
    
    # Print a message indicating that the simulation is done and the average density is obtained
    print(f'NPT relaxation at {T:.1f}K is done! density is {density:.5f} g/cm^3')
    
    # Change back into the original directory
    os.chdir(original_cwd)
    
    # Return the average density value
    return density



def prepare_slab(seq_name, seq_len, N=64):
    """
    Prepare a LAMMPS simulation config  for an IDP
    
    Parameters:
        seq_name (str): Name of the sequence (e.g., "hnRNPA1")
        seq_len (int): Length of the sequence
        N (int, optional): Number of copies. Defaults to 64.
    
    Returns:
        None
    """

    # directory name for the simulation based on the sequence name
    sim_dir = 'NPT_' + seq_name

    # Save the current working directory so we can go back to it later
    original_cwd = os.getcwd()
    
    # Change into the NPT simulation directory
    os.chdir(sim_dir + '/')

    # Extract the last frame from the compressed trajectory file and save it as "lastFrame.dat"
    # The number of lines to extract is seq_len * 64 (num atoms) + 9 (header lines)
    os.system('tail -n ' + str(seq_len*N+9) + ' compressed_' + seq_name + '_100.lammpstrj > lastFrame.dat')

    # Open the "lastFrame.dat" file and read its contents
    frame = open('lastFrame.dat')
    
    # Read lines 5-7 (xlo, xhi, ylo, yhi, zlo, zhi) to determine the simulation box dimensions
    content = frame.readlines()
    
    # Extract the box dimensions from the file content
    xlo = float(content[5].split()[0])  # xlo coordinate
    xhi = float(content[5].split()[1])  # xhi coordinate
    ylo = float(content[6].split()[0])  # ylo coordinate
    yhi = float(content[6].split()[1])  # yhi coordinate
    zlo = float(content[7].split()[0])  # zlo coordinate
    zhi = float(content[7].split()[1])  # zhi coordinate

    # Calculate the midpoints of the box dimensions
    mid_x = (xlo + xhi) / 2
    mid_y = (ylo + yhi) / 2
    mid_z = (zlo + zhi) / 2

    # Calculate the absolute differences between the box coordinates (box dimensions)
    dim_x = abs(xlo - xhi)
    dim_y = abs(ylo - yhi)
    dim_z = abs(zlo - zhi)

    # Extract only the last seq_len * 64 lines from the compressed trajectory file and save it as "lastFrame.dat"
    os.system('tail -n ' + str(seq_len*N) + ' compressed_' + seq_name + '_100.lammpstrj > lastFrame.dat')

    # Load the atomic positions from the "lastFrame.dat" file into a NumPy array
    pos_data = np.loadtxt('lastFrame.dat')

    # Change back to the original working directory
    os.chdir(original_cwd)

    # Open a new file for writing the LAMMPS simulation configuration
    config = open('slab_config_' + seq_name + '.dat', "w")

    # Define some constants (num atoms, num bonds)
    numAtoms = seq_len * 64
    
    # Write some header lines to the configuration file
    config.write("LAMMPS Script for IDP " + seq_name + " \n")
    config.write("\n")
    
    # Write the number of atoms and bonds to the configuration file
    config.write("{} atoms\n".format(numAtoms))
    config.write("{} bonds\n".format((seq_len-1)*N))

    # Define some type counts (atom types, bond types, angle types)
    config.write("\n")
    config.write("20 atom types\n")
    config.write("1 bond types\n")
    # config.write("1 angle types\n")

    # Write the box dimensions to the configuration file
    config.write("\n")
    
    # only extend the simulation box in x direction
    config.write("{} {} xlo xhi\n".format(-3*dim_x, 3*dim_x))
    config.write("{} {} ylo yhi\n".format(-dim_y/2, dim_y/2))
    config.write("{} {} zlo zhi\n".format(-dim_z/2, dim_z/2))

    # Write atom section to configuration file
    config.write("\n\nAtoms\n\n")

    # Initialize counter for atomic number (k)
    k = 0

    # Iterate over positions of all atoms
    for i in range(numAtoms):
        # Increment k by 1 for each atom
        k += 1
        
        # Write position and ID of current atom to configuration file
        config.write(f"{k} {pos_data[i,1]:.0f} {pos_data[i,2]:.0f} 0 ")  
        
        # Calculate and write positions relative to mid_x, mid_y, mid_z
        config.write(f"{pos_data[i,-3]-mid_x:.3f} {pos_data[i,-2]-mid_y:.3f} {pos_data[i,-1]-mid_z:.3f}\n")  # x, y, z coordinates centered around (mid_x, mid_y, mid_z)

    # Write bond section to configuration file
    config.write("\n\nBonds\n\n")

    # Initialize counter for bond ID (b_id)
    b_id = 0

    # Iterate over residue indices and sequence length
    for i in range(N):
        # Iterate over residues within current residue index
        for k in range(seq_len-1):
            b_id += 1
            
            # Write bond information to configuration file
            config.write(f"{b_id} {1} {i*seq_len+k+1} {i*seq_len+k+2}\n")  # ID, type (assumed to be 1), and two atom IDs involved in the bond

    # Close configuration file
    config.close()



def run_NVT(seq_name, seq_len, T_list, N = 64):
    """
    This function submit multiple NVT slab simulations at different temperatures
    for a given protein sequence.

    Parameters:
        seq_name (str): Name of the protein sequence
        seq_len (int): Length of the protein sequence
        T_list (list): List of temperatures in Kelvin to run simulations at
        N (int, optional): Number of copies. Defaults to 64.

    Returns:
        None (simulations are run on a remote cluster)
    """

    # Create the base directory for NVT simulations
    sim_dir = 'NVT_' + seq_name
    original_cwd = os.getcwd()  # Save the current working directory

    # Remove any existing directories with this name and create a new one
    os.system('rm -r ' + sim_dir)  # Remove any existing directory
    os.system('mkdir ' + sim_dir)  # Create a new directory

    # Loop over each temperature in the list
    for T in T_list:
        # Create a subdirectory for this temperature
        sim_dir_T = sim_dir + f'/{T:.0f}K'
        os.system('mkdir ' + sim_dir_T)  # Create a new directory for this temperature

        # Copy the LAMMPS input file template to the current directory
        lmp_file = 'NVT_slab.in'
        os.system('cp NVT_slab_template.in ' + lmp_file)

        # Open the LMP file and modify it with the correct parameters
        with open(lmp_file, 'r+') as f:
            content = f.read()  # Read the current contents of the file

            # Move the cursor to the beginning of the file
            f.seek(0, 0)

            # Write the new parameters (temperature and number of atoms) to the top of the file
            f.write(f'variable myTemp equal {T:.0f}\n')  # Set the temperature variable
            f.write(f'variable numAtom equal {seq_len*N:.0f}\n')  # Set the number of atoms variable
            f.write(f'variable IDP string ' + seq_name + '\n')  # Set the protein sequence variable

            # Write the rest of the contents back to the file
            f.write('\n' + content)

        # Move the modified LMP file to the correct directory
        os.system('mv ' + lmp_file + ' ' + sim_dir_T)

        # Copy other necessary files to the current directory
        os.system('cp potential.param ' + sim_dir_T)  # Copy the potential parameters file
        os.system('cp ' + 'slab_config_' + seq_name + '.dat ' + sim_dir_T)  # Copy the slab configuration file

        # Create a new job script for this temperature
        sh_file = 'job.sh'
        os.system('cp job_slab_NVT_template.sh ' + sh_file)

        # Open the job script and modify it with the correct parameters
        with open(sh_file, 'r+') as f:
            content = f.read()  # Read the current contents of the file

            # Move the cursor to the beginning of the file
            f.seek(0, 0)

            # Write the new parameters (job name and temperature) to the top of the file
            f.write('#!/bin/bash\n')  # Set the shebang line
            f.write(f'\n#SBATCH --job-name={T:.0f}-{seq_name}\n')  # Set the job name variable

            # Write the rest of the contents back to the file
            f.write('\n' + content)

        # Move the modified job script to the correct directory
        os.system('mv ' + sh_file + ' ' + sim_dir_T)

        # Run the LAMMPS simulation on the remote cluster
        os.chdir(sim_dir_T)  # Change to the current directory
        os.system('sbatch job.sh')  # Run the job

        # Go back to the original working directory
        os.chdir(original_cwd)

#----------------------------------------------------------------------------
# Main function: this script performs a series of simulations using LAMMPS

# Get the input sequence file from the command line arguments
filename = sys.argv[1]  # sequence file from the terminal

# Specify the path to the LAMMPS executable (assuming it's installed in ~/.local/bin)
# Lammps = '~/.local/bin/lmp_della_test_full'  # LAMMPS executable
Lammps = sys.argv[2]  # LAMMPS executable

# Prepare the NPT simulation parameters by reading the input sequence file
seq_name, seq_len, n_copies = prepare_NPT(filename)

# Set the number of compression steps (30000)
n_compress = 30000

# Run the compress_NPT simulation to create an initial compressed state
compress_NPT(seq_name, Lammps, 100, 30000)  # NPT simulation at 100 K and n_compress steps

# Define the temperature range for the subsequent relaxation simulations
T_low = 100  # Lower bound of the temperature range (K)
T_hi = 500   # Upper bound of the temperature range (K)
delta_T = T_hi - T_low  # Initial temperature difference

# Perform the relax_NPT simulation at the lower and upper temperatures to obtain initial density estimates
density_l = relax_NPT(seq_name, Lammps, T_low, n_compress)  # Relaxation at low temperature (T_low K)
density_h = relax_NPT(seq_name, Lammps, T_hi, n_compress)  # Relaxation at high temperature (T_hi K)

# Perform a bisection search to estimate the critical temperature for the subsequent simulations
while delta_T > 2:
    # Calculate the new midpoint temperature
    newT = int((T_hi + T_low) / 2)
    
    # Run the relax_NPT simulation at the new temperature and obtain a new density estimate
    newDensity = relax_NPT(seq_name, Lammps, newT, n_compress)
    
    # Update the lower or upper bound of the temperature range based on the new density estimate
    if newDensity > 0.30:
        T_low = newT  # Lower bound updated to the new midpoint temperature
        density_l = newDensity  # New density estimate at low temperature
    else:
        T_hi = newT  # Upper bound updated to the new midpoint temperature
        density_h = newDensity  # New density estimate at high temperature
    
    # Update the temperature difference for the next iteration
    delta_T = T_hi - T_low

# Prepare the input files and parameters for the slab simulation
prepare_slab(seq_name, seq_len)

# Run the NVT simulations at 4 temperatures below the critical temperature found above
run_NVT(seq_name, seq_len, [T_low-5, T_low-10,T_low-15,T_low-20])


    
