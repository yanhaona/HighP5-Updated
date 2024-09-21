#!/bin/bash

echo "This script requires the following command line inputs:"
echo "	1. number of rows in the input plate"
echo "	2. number of columns in the input plate"
echo "	3. number of jaqobi iterations"
echo "	4. amount of padding between successive processes"
echo "	5. number of MPI processes"

plate_rows=$1
plate_cols=$2
iterations=$3
padding=$4
process_per_machine=$5
write_output=0
input_file=a.bin
array_generator=../../../../tools/binary-array-generator

# generate input matrix
$array_generator 3 2 $input_file $plate_rows $plate_cols

# run the MPI program
mpirun -n $process_per_machine --oversubscribe mpi-stencil.o $input_file $iterations $padding $write_output > output.txt
cat output.txt

# delete the input file
rm $input_file
