#!/bin/bash

matrix_dim=$1
iterations=$2
process_per_machine=$3
sparsity=90
write_output=0
array_generator=../../../../tools/binary-array-generator
matrix_generator=../../../../tools/sparse-matrix-generator

# generate the arrays correspond to sparse matrix
$matrix_generator $matrix_dim $matrix_dim $sparsity 3 1

# generate known and prediction vectors
$array_generator 3 2 known $matrix_dim $matrix_dim
$array_generator 3 2 pred $matrix_dim $matrix_dim

# run the MPI program
mpirun -n $process_per_machine --oversubscribe mpi-conjgrad.o \
       values columns rows known pred $iterations $write_output > output.txt
cat output.txt

# delete the input files
rm values columns rows known pred
