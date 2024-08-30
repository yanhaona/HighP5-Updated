#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# adjust these two variables based on the mode of the experiment
ref_target=brac-node
parallelism=( 20 12 8 6 2 1 )
executable=../executables/pthread-bluf.o
array_generator=../../../../../tools/array-generator

# parallel run counts and input size configuration
parallel_runs=2
input_size=1024
cache_block=64

# make data directory for experiments
mkdir -p ../../../$ref_target/data/pthread/Bluf
root_data_dir=../../../$ref_target/data/pthread/Bluf

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# copy the executable in the root data directory
cp $executable $root_data_dir/pthread_exec.o

# generate the reference input data file using the array generator
./${array_generator} $input_size*$input_size double 0.25 1000 $root_data_dir/a.txt

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
		../../pthread_exec.o $cache_block $version $input_size > output.txt
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
