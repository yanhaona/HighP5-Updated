#!/bin/bash
#PJM -L "node=100"
#PJM -L "rscgrp=small"
#PJM -L "elapse=1:00:00"
#PJM -g hp250392

#PJM --mpi "shape=100"
#PJM --mpi "max-proc-per-node=8"

#PJM -s

# setting command line parameters
values_file=values-51200-90-percent
rows_file=rows-51200-90-percent
cols_file=columns-51200-90-percent
known_vector_file=b_51200
pred_vector_file=x_51200
iterations=51200
writing=0

# run the MPI program
mpiexec ./conj-grad.o $values_file $cols_file $rows_file $known_vector_file $pred_vector_file $iterations $writing 
