#!/bin/bash
#PJM -L "node=6"
#PJM -L "rscgrp=small"
#PJM -L "elapse=0:20:00"
#PJM -g hp250392

#PJM --mpi "shape=6"
#PJM --mpi "max-proc-per-node=32"

#PJM -s

# setting command line parameters
block_size=64
write_output=0
input_file_1=a-10240
input_file_2=a-10240

# run the MPI program
mpiexec ./mmult.o $block_size $input_file_1 $input_file_2 $write_output
