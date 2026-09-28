#!/bin/bash
#PJM -L "node=64"
#PJM -L "rscgrp=small"
#PJM -L "elapse=2:00:00"
#PJM -g hp260450

#PJM --mpi "shape=64"
#PJM --mpi "max-proc-per-node=1"

# setting up the thread count configuration
export PARALLEL=48
export OMP_NUM_THREADS=${PARALLEL}

#PJM -s

# setting command line parameters
cache_block=64
block_part=256
A=a-49152

# run the MPI program
mpiexec ./bluf-64-mod0.o input=$A block_size=$cache_block part_size=$block_part

# clean extra log files
rm -v segment_1*
rm -v segment_2*
rm -v segment_3*
rm -v segment_4*
rm -v segment_5*
rm -v segment_6*
