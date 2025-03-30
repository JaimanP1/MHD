#!/bin/bash
#SBATCH -p solarlab
#SBATCH --ntasks-per-node=1
# Number of nodes
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --nodelist=node[819]
#SBATCH --mem-per-cpu=1G
# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1
module use /opt/site/easybuild/modules/all/Core
module load GCC/9.3.0 OpenMPI/4.0.3
# Load the module environment suitable for the job
echo "Starting job at: "
date
srun ./bline
echo "Finished"
date
