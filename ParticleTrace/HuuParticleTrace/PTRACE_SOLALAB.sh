#!/bin/bash -l
#SBATCH --job-name="o-point_p-trace"
#SBATCH --output=OutFiles/o-point_p-trace%j.out
#SBATCH --error=ErrFiles/o-point_p-trace%j.err
#SBATCH --partition=bigmem
#SBATCH --account=wangj
#SBATCH --qos=high_wangj
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=4G

# Purge and load the correct modules
module purge > /dev/null 2>&1
module load wulver
module load foss/2022b
module load GCC/12.2.0
module use /opt/site/easybuild/modules/all/Core
make

# Run the mpi program
srun ./bline

