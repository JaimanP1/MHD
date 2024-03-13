#!/bin/sh
#PJM -L "rscgrp=cx2-middle"
#PJM -L "vnode=100"
#PJM -L "elapse=10:00:00"
#PJM -j

mpiexec -n 100 ./sdo_mhd_1_ar11012
