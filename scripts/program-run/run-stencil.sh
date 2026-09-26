#!/bin/bash

# parameter settings
rows=$1
columns=$2
iterations=$3
p1=$4
p2=$5
processes=4
threads=2

echo "removing previous versions of the input file"
rm plate > /dev/null 2>&1 
echo "Generating the input plate"
../tools/binary-array-generator 3 2 plate $rows $columns

echo "Running the stencil program"
mpirun --oversubscribe -n=$processes ./stencil.o input_file=plate k=$processes l=1 m=$threads n=1 p2=$p1 p1=$p2 iterations=$iterations

echo "removing the input files"
rm plate > /dev/null 2>&1

echo "removing the log files"
rm *.log
