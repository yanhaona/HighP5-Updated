#!/bin/bash

cell_length=100 
grid_dim=1024
points_per_cell=1000
process_per_machine=12

# run the MPI program
mpirun -n $process_per_machine --oversubscribe mpi-monte.o $cell_length $grid_dim $points_per_cell > output.txt
cat output.txt

