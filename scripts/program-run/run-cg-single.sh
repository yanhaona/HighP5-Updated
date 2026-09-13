#!/bin/bash

# parameter settings
sparsity=90
dimLength=$1
iterations=$2
point_sparsity=$(( sparsity / 100 ))
processes=4
threads=2
blocks=$(( processes * threads ))

echo "removing previous versions of the input files"
rm b x rows columns values > /dev/null 2>&1 
echo "Generating the known vector"
../tools/binary-array-generator 3 1 b $dimLength
echo "Generating the prediction vector"
../tools/binary-array-generator 3 1 x $dimLength
echo "Generating the sparse matrix"
../tools/sparse-matrix-generator $dimLength $dimLength $sparsity 3 1

echo "Running the conjugate gradient program"
mpirun --oversubscribe -n=$processes ./cg-single.o \
	arg_matrix_rows=rows arg_matrix_cols=columns arg_matrix_values=values \
	known_vector=b prediction_vector=x \
	maxIterations=$iterations \
	r=$blocks


echo "removing the input files"
rm b x rows columns values > /dev/null 2>&1
echo "removing the log files"
rm *.log
