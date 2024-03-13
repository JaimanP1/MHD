#!/bin/bash
#PBS -q NODEX1
#PBS -l elapstim_req=48:00:00
#PBS -lcpunum_job=16
#PBS -T intmpi
#PBS -j o

mpirun ${NQSII_MPIOPTS} -n 16 /cidashome/sc/c0075inoue/MAIN_PROGRAM_2019/MHD_FULL/MHD_VER0/mhd_run



