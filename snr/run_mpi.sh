#!/bin/bash
#SBATCH --job-name=xray          # create a short name for your job
#SBATCH --output=slurm-%A.out # stdout file
#SBATCH --error=slurm-%A.err  # stderr file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=12              # total number of tasks across all nodes
##SBATCH --exclusive
#SBATCH --mem=2000G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=06:00:00          # total run time limit (HH:MM:SS)

module purge
module load openmpi/gcc/4.1.0 
module load anaconda3/2021.11
conda activate snr-xray

cd $HOME/snr/snr/

# srun python xray.py snrnet_h05_t16_128_k0
# srun python xray.py snrnet_h05_t16_256_k0
srun python xraynew.py $1
