#SBATCH --nodes=1                # node count
#SBATCH --ntasks=8               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=200MB         # memory per cpu-core (4G is default)
#SBATCH --time=23:50:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all        # send email when job begins
##SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=da9666@princeton.edu
#SBATCH --constraint=cascade,skylake

module purge
module load intel/19.1.1.217
module load intel-mpi/intel/2019.7
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

srun $HOME/.local/bin/lmp_della_double -in NVT_slab.in
