#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# adjust these two variables based on the mode of the experiment
ref_target=brac-node
machine=brac
parallelism=( 20 12 8 6 2 1 )
executable=../executables/pthread-mmult.o
array_generator=../../../../../tools/array-generator

# parallel run counts and input size configuration
parallel_runs=2
input_size=1024
cache_block=64

# make data directory for experiments
mkdir -p ../../../$ref_target/data/pthread/MMult
root_data_dir=../../../$ref_target/data/pthread/MMult

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# copy the executable in the root data directory
cp $executable $root_data_dir/pthread_exec.o

# generate reference input data files using the array generator
./${array_generator} $input_size*$input_size double 0.25 1000 $root_data_dir/m1.txt
./${array_generator} $input_size*$input_size double 0.25 1000 $root_data_dir/m2.txt

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

		# create an input file to pass user inputs
		echo ../../m1.txt > input.txt
		echo ../../m2.txt >> input.txt

		# execute the program using the input file piping command line inputs
		(time ../../pthread_exec.o m1 m2 $cache_block $version $machine < input.txt) > output.txt 2>&1
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
