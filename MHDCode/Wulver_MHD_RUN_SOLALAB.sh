#!/bin/bash -l
# The above line must be first and must include the "-l"

#SBATCH --job-name=mhd_run
#SBATCH --output=mhd_run%j.out # %j expands to slurm JobID
#SBATCH --error=mhd_run%j.err
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=128
#SBATCH --qos=high_wangj

# Use "sinfo" to see what partitions are available to you
#SBATCH --partition=general

# Memory required; lower amount gets scheduling priority
#SBATCH --mem-per-cpu=2G 

# Time required in d-hh:mm:ss format; lower time gets scheduling priority
#SBATCH --time=72:00:00

# Purge and load the correct modules
module purge > /dev/null 2>&1
module use /opt/site/easybuild/modules/all/Core
module load wulver
module load intel/2022b
make

# Run the mpi program
srun mhd_run

