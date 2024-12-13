LAMMPS = '~/.local/bin/lmp_della_double'

# copy files
#! Edit the email and excutables in the template file as well
for idp in IDP_list:
    folder = idp+'/'
    os.system('cp job_slab_NVT_template.sh '+ folder)
    os.system('cp job_start_template.sh '+ folder)
    os.system('cp potential.param '+ folder)
    os.system('cp NPT_compress_template.in '+ folder)
    os.system('cp NPT_relax_template.in '+ folder)
    os.system('cp NVT_slab_template.in '+ folder)
    os.system('cp main.py '+ folder)
    # os.system('cp main.py '+ folder)

n = 1
for idp in IDP_list:
    original_cwd = os.getcwd()
    os.chdir(idp + '/')
    # Create a new LAMMPS input file for NPT relaxation simulation
    job_file = 'job_' + idp + '.sh'
    
    # Copy the template LAMMPS input file to the new one
    os.system('cp job_start_template.sh ' + job_file)
    # Open the new LAMMPS input file for modification
    with open(job_file, 'r+') as f:
        # Read the content of the file
        content = f.read()
        # Move the cursor to the beginning of the file
        f.seek(0, 0)
        
        # Write the modified variables and parameters to the file
        f.write(f'#!/bin/bash\n')
        f.write(f'#SBATCH --job-name={n:.0f}-{idp:s}')
        # Append the original content to the file
        f.write('\n' + content)
        f.write('\n\npython3 main.py '+idp+'.fasta '+LAMMPS)
        n += 1
    os.system('sbatch job_start_template.sh')
    os.chdir(original_cwd)