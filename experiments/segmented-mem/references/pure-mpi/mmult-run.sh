#!/bin/bash

matrix_dim=$1
process_per_machine=$2
block_size=64
write_output=0
input_file_1=a.bin
input_file_2=b.bin
array_generator=../../../../tools/binary-array-generator

# generate input matrix
$array_generator 3 2 $input_file_1 $matrix_dim $matrix_dim
$array_generator 3 2 $input_file_2 $matrix_dim $matrix_dim

# run the MPI program
mpirun -n $process_per_machine --oversubscribe mpi-mmult.o $block_size $input_file_1 $input_file_2 $write_output > output.txt
cat output.txt

# delete the input file
rm $input_file_1
rm $input_file_2
