#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# adjust these two variables based on the mode of the experiment
ref_target=brac-node
machine=brac
parallelism=( 20 12 8 6 2 1 )
executable=../executables/pthread-conjgrad.o

# parallel runs count
parallel_runs=2

# sampling parameters
matrix_size=10240
sparsity=90
maxIterations=100
precision=0.0001

# make data directory for experiments
mkdir -p ../../../$ref_target/data/pthread/ConjGrad
root_data_dir=../../../$ref_target/data/pthread/ConjGrad

# array generator executables
sparse_matrix_generator=../../../../../tools/sparse-matrix-generator
array_generator=../../../../../tools/array-generator


# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# copy the executable in the root data directory
cp $executable $root_data_dir/pthread_exec.o

echo "generating input data"
# generate sparse matrix and move associated files in the root data directory
echo "values" > tmp.txt
echo "columns" >> tmp.txt
echo "rows" >> tmp.txt
$sparse_matrix_generator $matrix_size $matrix_size $sparsity 3 0 < tmp.txt
cp values ${root_data_dir}/values
cp columns ${root_data_dir}/columns
cp rows ${root_data_dir}/rows
rm values columns rows
rm tmp.txt

# generate predication and known vectors
$array_generator $matrix_size double 0.25 100 ${root_data_dir}/known $matrix_size
$array_generator $matrix_size double 0.25 100 ${root_data_dir}/pred $matrix_size


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

		# generate an input file for the program
		echo "../../columns" > input.txt
		echo "../../rows" >> input.txt
		echo "../../values" >> input.txt
		echo "../../known" >> input.txt
		echo "../../pred" >> input.txt

		# execute the program using the input file piping command line inputs
		(time ../../pthread_exec.o $maxIterations $precision $version $machine < input.txt) > output.txt 2>&1
		echo "Matrix Dimension $matrix_size by $matrix_size" >> output.txt
		echo "Sparsity $sparsity" >> output.txt
		echo "Refinement iterations $maxIterations" >> output.txt
		echo "Refinement precision $precision" >> output.txt
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
