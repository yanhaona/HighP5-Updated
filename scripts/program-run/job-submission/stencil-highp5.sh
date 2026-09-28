#!/bin/bash
#PJM -L "node=64"
#PJM -L "rscgrp=small"
#PJM -L "elapse=2:00:00"
#PJM -g hp260450

#PJM --mpi "shape=64"
#PJM --mpi "max-proc-per-node=1"
#PJM -s

# setting up the thread count configuration
export PARALLEL=48
export OMP_NUM_THREADS=${PARALLEL}

# setting command line parameters
input_file=plate-409600-2048
upper_block=64
lower_block=48
p1=32
p2=2
iterations=1250


# run the MPI program
mpiexec ./stencil-64-orig.o input_file=$input_file k=$upper_block l=1 m=$lower_block n=1 p1=$p1 p2=$p2 iterations=$iterations

# remove extra log files
rm -v segment_1*
rm -v segment_2*
rm -v segment_3*
rm -v segment_4*
rm -v segment_5*
rm -v segment_6*
