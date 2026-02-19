#!/bin/bash -l                                                                   
#SBATCH --job-name=twist_mhd
#SBATCH --output=Out_files/%x.%j.out                                                       
#SBATCH --error=Err_files/%x.%j.err                                                        
#SBATCH --partition=general                                                      
#SBATCH --qos=high_wangj                                                        
#SBATCH --account=si22                                                          
#SBATCH --nodes=1                                                                
#SBATCH --ntasks-per-node=11                                                     
#SBATCH --time=10:00:00  # D-HH:MM:SS                                            
#SBATCH --mem-per-cpu=4000M  # (4000M = 4 GB)                                    

# Load the modules required for the job                                          
module load foss/2025a
module load GCC/14.2.0

make
srun ./twist_mhd
