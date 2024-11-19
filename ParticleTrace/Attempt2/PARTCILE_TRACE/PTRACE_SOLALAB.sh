#!/bin/bash
#SBATCH --job-name="NLFFF MHD ANA PTRACE"
#SBATCH --partition=general
#SBATCH --account=wangj
#SBATCH --qos=high_wangj
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=4G

module load foss/2022b
module load GCC/12.2.0
make

echo "Starting job at: "
date
srun ./bline
echo "Finished"
date

