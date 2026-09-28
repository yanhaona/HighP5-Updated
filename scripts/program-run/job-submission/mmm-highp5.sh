#!/bin/bash
#PJM -L "node=32"
#PJM -L "rscgrp=small"
#PJM -L "elapse=0:15:00"
#PJM -g hp250392

#PJM --mpi "shape=32"
#PJM --mpi "max-proc-per-node=1"

# setting up the thread count configuration
export PARALLEL=32
export OMP_NUM_THREADS=${PARALLEL}

#PJM -s

# setting command line parameters
cache_block=128
A=a-10240
B=a-10240

# run the MPI program
mpiexec ./mmult.o input_file_1=$A input_file_2=$B k=$cache_block l=$cache_block q=$cache_block
