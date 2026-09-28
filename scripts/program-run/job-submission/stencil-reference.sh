#!/bin/bash
#PJM -L "node=64"
#PJM -L "rscgrp=small"
#PJM -L "elapse=4:00:00"
#PJM -g hp260450

#PJM --mpi "shape=64"
#PJM --mpi "max-proc-per-node=32"

#PJM -s

# setting command line parameters
plate=plate-204800-2048
jacobi=25000
padding=4
write_output=0

# run the MPI program
mpiexec ./stencil.o $plate $jacobi $padding $write_output
