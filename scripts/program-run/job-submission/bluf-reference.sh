#!/bin/bash
#PJM -L "node=200"
#PJM -L "rscgrp=small"
#PJM -L "elapse=1:00:00"
#PJM -g hp260450

#PJM --mpi "shape=200"
#PJM --mpi "max-proc-per-node=32"

#PJM -s

# setting command line parameters
block_size=12
write_output=0
input_file=a-76800

# run the MPI program
mpiexec ./bluf.o $block_size $input_file $write_output
