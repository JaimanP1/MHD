# MHD
MHD code from work with Prof. Satoshi Inoue, CSTR, NJIT using the HPC Wulver.

This code outputs both binary and png data to visualize magnetic flux ropes in quadrupolar photospheric topologies for coronal mass ejections. The visualization process is done through NCAR VAPOR API. 

This code is parallelized to run on up to 144 cores on 2 nodes. 

Recent updates include finishing the implementation of the API so data does not have to be downloaded locally for visualization.

Current work revolves around investigating the roles of MHD instability and magnetic reconnection in the Breakout Model.

Future work will likely include an effort to create a more elastic parallelization scheme and implement a more robust PDE solver.

Will also make a wiki page when I get to it.
