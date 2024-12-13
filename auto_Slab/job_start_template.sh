#SBATCH --nodes=1                # node count
#SBATCH --ntasks=8               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1G         # memory per cpu-core (4G is default)
#SBATCH --time=03:62:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all        # send email when job begins
##SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=da9666@princeton.edu
#SBATCH --constraint=cascade,skylake

module purge
module load intel/2022.2.0
module load intel-mpi/intel/2021.7.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
