#!/bin/bash
#PJM -L "node=4"
#PJM -L "rscgrp=small"
#PJM -L "elapse=0:10:00"
#PJM -g hp250392

#PJM --mpi "shape=4"
#PJM --mpi "max-proc-per-node=1"

# setting up the thread count configuration
export PARALLEL=48
export OMP_NUM_THREADS=${PARALLEL}

#PJM -s

# setting command line parameters
precision=0.0001
max_rounds=100
cell_length=10
grid_dim=100
points_per_cell=1000
block_size=50

# run the MPI program
mpiexec ./monte.o cell_length=$cell_length grid_dim=$grid_dim points_per_cell=$points_per_cell b=$block_size max_rounds=$max_rounds
