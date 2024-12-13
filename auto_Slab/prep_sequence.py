import pandas as pd
import os

df = pd.read_csv("IDRome_DB.csv")

# select the sequences based on the lengths
df1 = df[df['N']> 100]
df2 = df1[df1['N']< 200]


# sample 150 sequences
df3 = df2.sample(n=150)

# save the seuqeunces
df3.to_csv('./selected_sequences.csv', index=True)


# create the simulation folder the first person
folder1 = 'selected_1/'
os.system('rm -r '+folder1)
os.system('mkdir '+folder1)

# create the simulation folder the second person
folder2 = 'selected_2/'
os.system('rm -r '+folder2)
os.system('mkdir '+folder2)

# create the simulation folder the third person
folder3 = 'selected_3/'
os.system('rm -r '+folder3)
os.system('mkdir '+folder3)


# for the first person
for index, row in df3[:50].iterrows():
    prot_id = row['UniProt_ID']
    new_dir = folder1+prot_id+'/'

    # create the subfolder "/selected_1/uniprot_id/
    os.system('mkdir '+new_dir)
    sequence = row['fasta']

    # # create the .fasta sequence file "/selected_1/uniprot_id/uniprit_id.fasta
    f = open(new_dir+prot_id+'.fasta', "w")
    f.write(sequence + "\n")
    f.close()


# for the second person
for index, row in df3[50:100].iterrows():
    prot_id = row['UniProt_ID']
    new_dir = folder2+prot_id+'/'

    # create the subfolder "/selected_2/uniprot_id/
    os.system('mkdir '+new_dir)
    sequence = row['fasta']

    # # create the .fasta sequence file "/selected_2/uniprot_id/uniprit_id.fasta
    f = open(new_dir+prot_id+'.fasta', "w")
    f.write(sequence + "\n")
    f.close()

# for the third person
for index, row in df3[100:].iterrows():
    prot_id = row['UniProt_ID']
    new_dir = folder3+prot_id+'/'

    # create the subfolder "/selected_2/uniprot_id/
    os.system('mkdir '+new_dir)
    sequence = row['fasta']

    # # create the .fasta sequence file "/selected_3/uniprot_id/uniprit_id.fasta
    f = open(new_dir+prot_id+'.fasta', "w")
    f.write(sequence + "\n")
    f.close()


# Copy the template files to the subfolders
for i in range(1,4):
    folder = 'selected_'+str(i)+'/'
    os.system('cp job_slab_NVT_template.sh '+ folder)
    os.system('cp job_start_template.sh '+ folder)
    os.system('cp potential.param '+ folder)
    os.system('cp NPT_compress_template.in '+ folder)
    os.system('cp NPT_relax_template.in '+ folder)
    os.system('cp NVT_slab_template.in '+ folder)
    os.system('cp main.py '+ folder)

    # the final file to run in te cluster
    final_file = 'submit_jobs_'+str(i)+'.py'
    
    # Copy the template LAMMPS input file to the new one
    os.system('cp submit_jobs.py ' + final_file)
    # Open the new LAMMPS input file for modification
    with open(final_file, 'r+') as f:
        # Read the content of the file
        content = f.read()
        # Move the cursor to the beginning of the file

        # add the sequence list
        f.seek(0, 0)
        f.write('import os\n')
        f.write("# all sequences in the current folder\n")
        f.write("IDP_list = [")
        for index, row in df3[(i-1)*50:i*50].iterrows():
            f.write("\'"+row['UniProt_ID']+"\', ")
        f.write("]\n")
        f.write('\n'+content)
    os.system('mv ' + final_file+' '+folder)