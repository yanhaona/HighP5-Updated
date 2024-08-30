#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# adjust these two variables based on the mode of the experiment
ref_target=brac-node
parallelism=( 20 12 8 6 2 1 )
executable=../executables/pthread-monte-simple.o

# parallel runs count
parallel_runs=2

# sampling parameters
cell_length=1000
grid_dim=1024
points_per_cell=100

# make data directory for experiments
mkdir -p ../../../$ref_target/data/pthread/Monte
root_data_dir=../../../$ref_target/data/pthread/Monte

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# copy the executable in the root data directory
cp $executable $root_data_dir/pthread_exec.o

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
		(time ../../pthread_exec.o $cell_length $grid_dim $points_per_cell $version) > output.txt 2>&1
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
