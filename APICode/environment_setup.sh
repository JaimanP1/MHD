#DO NOT RUN THIS SCRIPT YET
#!/bin/bash

#compute node, try to get GPU
srun -p general -n 4 --qos=standard --account=si22 --time=1:00:00 --pty bash

module load Miniforge3

module load foss/2021b FFmpeg/4.3.2

module load foss/2021b ImageMagick/7.1.0-4

conda create --prefix=/your/project/directory

conda activate ./your/environment

echo `mamba --version`

mamba install -c conda-forge -c ncar-vapor vapor
