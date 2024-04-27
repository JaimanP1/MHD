# MHD
MHD code from work with Prof Satoshi Inoue, CSTR, NJIT

This code outputs both binary and vdc data to generate Magnetic Flux ropes in quadrupolar photosphere topologies for Coronal Mass Ejections. The vdc data can be visualized in VAPOR. 

This code is parallelized to run on up to 144 cores on 2 nodes. 

Refer to Wiki page for more information (in progress)

Recent updates include finishing the implementation of a function to stop the twisting motion at a specified time.

Current work revolves around investigating why halting the twisting motion yields (plasmoid/breakout) reconnection.

Future work will likely include an effort to create a more elastic parallelization scheme, using Galerkin/Finite-Element methods as opposed to Runge-Kutta methods, and comparing converging motions with twisting motions.
