#!/bin/bash -l
# The above line must be first and must include the "-l"

#SBATCH --job-name=pi_test.impi
#SBATCH --output=pi_test.impi.%j.out # %j expands to slurm JobID
#SBATCH --error=pi_test.impi.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --qos=high_wangj

# Use "sinfo" to see what partitions are available to you
#SBATCH --partition=general

# Memory required; lower amount gets scheduling priority
#SBATCH --mem-per-cpu=4000M 

# Time required in d-hh:mm:ss format; lower time gets scheduling priority
#SBATCH --time=01:00:00

# Purge and load the correct modules
module purge > /dev/null 2>&1
module load wulver
module load intel/2022b

mpiifort -o pi.impi pi.f

# Run the mpi program
srun pi.impi

