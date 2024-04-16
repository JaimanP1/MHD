#!/bin/bash -l
# The above line must be first and must include the "-l"

#SBATCH --job-name=mhd_merge
#SBATCH --output=Out_files/mhd_merge%j.out # %j expands to slurm JobID
#SBATCH --error=Err_files/mhd_merge%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --qos=high_wangj
#SBATCH --partition=general

# Use "sinfo" to see what partitions are available to you
#SBATCH --partition=general

# Memory required; lower amount gets scheduling priority
#SBATCH --mem-per-cpu=4G 

# Time required in d-hh:mm:ss format; lower time gets scheduling priority
#SBATCH --time=72:00:00

# Purge and load the correct modules
module purge > /dev/null 2>&1
module load wulver
module load foss/2022b
module use /opt/site/easybuild/modules/all/Core
make

# Run the mpi program
srun MergeCode

