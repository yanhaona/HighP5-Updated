#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

machine=tozammel
parallelism=( 128 64 32 16 8 )
executable=../executables/pthread-cg.o

# parallel run counts and input size configuration
matrix_size=$1
iterations=$2
parallel_runs=$3
sparsity=90


# make data directory for experiments
mkdir -p ../data/CG/pthread-high
root_data_dir=../data/CG/pthread-high

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# copy the executable in the root data directory
cp $executable $root_data_dir/pthread_exec.o

# tools references
array_generator=../tools/binary-array-generator
matrix_generator=../tools/sparse-matrix-generator

echo "generating input data"
# generate sparse matrix and move associated files in the root data directory
$matrix_generator $matrix_size $matrix_size $sparsity 3 1 > /dev/null
cp values ${root_data_dir}/values
cp columns ${root_data_dir}/columns
cp rows ${root_data_dir}/rows
rm values columns rows

# generate predication and known vectors
$array_generator 3 1 ${root_data_dir}/known $matrix_size > /dev/null
$array_generator 3 1 ${root_data_dir}/pred $matrix_size > /dev/null

echo "Experiment Description:"
echo "	sparse matrix of $matrix_size by $matrix_size size"
echo "	matrix sparsity $sparsity"
echo "	gradient descent iterations $iterations"
for version in "${parallelism[@]}"
do
	echo "--------------------------------------------------------------------------------------------------------"
	echo "running pthread $version-way parallel experiments"
	
	# create a subdirectory for the current degree of parallelism
	mkdir ${root_data_dir}/version-${version}
	exper_dir=${root_data_dir}/version-${version}

	# go to the experiment directory
	cd $exper_dir

	# complete the parallel runs
	for (( run=1; run<=parallel_runs; run++ ))
	do
		# create experiment sub-directory
		mkdir run-${run}
		# go inside that directory
		cd run-${run}

		# execute the program with proper command line input
                ( time ../../pthread_exec.o ../../values ../../columns ../../rows \
                        ../../known ../../pred $iterations $version $machine ) > output.txt 2>&1
		echo "sparse matrix of $matrix_size by $matrix_size size" >> output.txt
		echo "matrix sparsity $sparsity" >> output.txt
		echo "gradient descent iterations $iterations" >> output.txt

		# check if experiment run correctly
                run_checker=`cat output.txt | grep 'Execution Time' | wc -l`
                if [[ $run_checker == 1 ]]; then
                        echo "experiment number $run of version $version-way parallel code was successfull"
			execution_time=`cat output.txt | grep 'Execution Time'`
                        echo $execution_time
			computation_time=`cat output.txt | grep 'Computation Time'`
                        echo $computation_time
                else
                        echo "experiment number $run of version $version-way parallel code failed"
                fi

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done

# remove the bulky plate file from data directory
cd ${root_data_dir}
rm columns values rows known pred

cd $curr_dir
