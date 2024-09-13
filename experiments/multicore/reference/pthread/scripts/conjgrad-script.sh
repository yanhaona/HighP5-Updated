#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# adjust these two variables based on the mode of the experiment
ref_target=brac-node
machine=brac
parallelism=( 20 12 8 6 2 1 )
executable=../executables/pthread-conjgrad.o

# parallel runs count
parallel_runs=1

# sampling parameters
matrix_size=10240
sparsity=90
maxIterations=1000

# make data directory for experiments
mkdir -p ../../../$ref_target/data/pthread/ConjGrad
root_data_dir=../../../$ref_target/data/pthread/ConjGrad

# array generator executables
sparse_matrix_generator=../../../../../tools/sparse-matrix-generator
array_generator=../../../../../tools/binary-array-generator

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# copy the executable in the root data directory
cp $executable $root_data_dir/pthread_exec.o

echo "generating input data"
# generate sparse matrix and move associated files in the root data directory
$sparse_matrix_generator $matrix_size $matrix_size $sparsity 3 1
# copy files to experiment directory
cp values ${root_data_dir}/values
cp columns ${root_data_dir}/columns
cp rows ${root_data_dir}/rows
# remove versions stored in current directory
rm values columns rows

# generate predication and known vectors in experiment root directory
$array_generator 3 1 ${root_data_dir}/known $matrix_size
$array_generator 3 1 ${root_data_dir}/pred $matrix_size


for version in "${parallelism[@]}"
do
	echo "--------------------------------------------------------------------------------------------------------"
	echo "running pthread $version-way parallel experiments"
	echo "Executable $executable"
	
	# create a subdirectory for the current degree of parallelism
	mkdir ${root_data_dir}/version-${version}
	exper_dir=${root_data_dir}/version-${version}

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

		# execute the program using the input file piping command line inputs
		(time ../../pthread_exec.o ../../values ../../columns ../../rows \
			../../known ../../pred $maxIterations $version $machine) > output.txt 2>&1
		echo "Matrix Dimension $matrix_size by $matrix_size" >> output.txt
		echo "Sparsity $sparsity" >> output.txt
		echo "Refinement iterations $maxIterations" >> output.txt
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
