# Auto Slab Simulations
The files in this folder for automated slab simulationsfor phase diagram predictions


The work flow is described as below:
## Email and LAMMPS
1. Change the email addresses in the job bash script templates [job_slab_NVT_template.sh](./job_slab_NVT_template.sh) and [job_start_template.sh](./job_start_template.sh).
2. Change the LAMMPS excutable name in [job_slab_NVT_template.sh](./job_slab_NVT_template.sh) and [submit_jobs.py](./submit_jobs.py)
## Prepare simulation
3. [prep_sequence.py](./prep_sequence.py) is the file that prepares the files. It samples 150 sequences from the data base. Then it creates  three folders to have 50 sequences for each. In each subfolder, there will be 50 folders with their `.fasta` sequence file in it. After that it copies the necessary files for simulations to the three main subfolders. Finally, it will add the sequence list to the header of the file [submit_jobs.py](./submit_jobs.py) while copiying.
### **This is the part need to be changed correspondly.** 
4. In each main subfolder, there will be a `submit_jobs.py` file. In it, there will be the list of the sequences, similar to following
```py
import os
# all sequences in the current folder
IDP_list = ['P0CB33', 'Q5T681', 'Q9Y297', 'Q5JTW2', 'Q3ZCQ2', 'P49257', 'Q9NS39', 'Q96Q04', 'Q15424', 'P78385', 'O14964', 'P04156', 'P52952', 'Q8IY17', 'Q7Z6I8', 'Q9Y5W7', 'Q9NQW5', 'Q6UB98', 'Q9NQR7', 'Q96T17', 'Q14511', 'Q53EL6', 'Q8IVB4', 'Q12955', 'Q6ZMS4', 'Q8ND25', 'Q9NRB3', 'Q9H0H9', 'Q92889', 'Q6UWB1', 'Q8WWZ4', 'Q6ZWB6', 'O95263', 'Q8NA47', 'Q30KQ4', 'Q9NRC6', 'P08581', 'P54727', 'Q9H1A4', 'Q15052', 'Q14204', 'Q9P2C4', 'O75592', 'Q86W26', 'P43489', 'A8MQ03', 'Q96LC9', 'Q96N64', 'Q499Y3', 'Q6ZVM7', ]

LAMMPS = '~/.local/bin/lmp_della_double'
```
5. Finally, in the subfolder run the following command, all simulations will be submitted accordingly.
```sh
user@della: /selected_1/$ python3 submit_jobs_1.py
```