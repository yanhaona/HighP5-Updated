#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# adjust these two variables based on the mode of the experiment
exe_class=pinned
mode=../executables/thread-pinned

# parallel run counts and problem size configurations
parallel_runs=2
matrix_size=1024
sparsity=90
maxIterations=10
block_size=64
precision=0.0001

# array generator executables
sparse_matrix_generator=../../../../tools/sparse-matrix-generator
array_generator=../../../../tools/binary-array-generator

cg_dir=${mode}/ConjGrad
parallelism=( 20 12 8 6 2 1 )

# make data directory for experiments
mkdir -p ../data/${exe_class}/ConjGrad
root_data_dir=../data/${exe_class}/ConjGrad

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

echo "generating input data"
# generate sparse matrix and move associated files in the root data directory
$sparse_matrix_generator $matrix_size $matrix_size $sparsity 3 1
cp values ${root_data_dir}/values
cp columns ${root_data_dir}/columns
cp rows ${root_data_dir}/rows
rm values columns rows

# generate predication and known vectors
$array_generator 3 1 ${root_data_dir}/known $matrix_size
$array_generator 3 1 ${root_data_dir}/pred $matrix_size


for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=${cg_dir}/conj-grad-${exe_class}-${version}-way.o
	echo "--------------------------------------------------------------------------------------------------------"
	echo "running $version-way parallel experiments"
	echo "Executable $executable"
	
	# create a subdirectory for the current executable
	mkdir ${root_data_dir}/version-${version}
	exper_dir=${root_data_dir}/version-${version}

	# copy the executable in the created directory
	cp $executable $exper_dir/executable.o

	# go to the experiment directory
	cd $exper_dir

	# complete the parallel runs
	for (( run=1; run<=parallel_runs; run++ ))
	do
		echo ""
		echo ""
		echo "Experiment run $run"

		# create experiment sub-directory
		mkdir run-${run}
		# go inside that directory
		cd run-${run}

		let upper_block=1
		let lower_block=$version
		
		# execute the program with proper command line input
		( time ../executable.o \
			arg_matrix_cols=../../columns arg_matrix_rows=../../rows \
			arg_matrix_values=../../values known_vector=../../known \
			prediction_vector=../../pred b=$block_size r=$block_size \
			maxIterations=$maxIterations precision=$precision ) > output.txt 2>&1
		echo "matrix of $matrix_size by $matrix_size" >> output.txt
		echo "matrix sparsity $sparsity" >> output.txt
		echo "gradient descent iterations $maxIterations" >> output.txt
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
