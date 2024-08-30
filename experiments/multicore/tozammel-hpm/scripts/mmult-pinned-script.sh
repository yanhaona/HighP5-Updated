#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# adjust these two variables based on the mode of the experiment
exe_class=pinned
mode=../executables/thread-pinned

# parallel run counts and input size configuration
parallel_runs=2
input_size=1024
cache_block=64

mmult_dir=${mode}/MMult
parallelism=( 1 2 64 128 )

# make data directory for experiments
mkdir -p ../data/${exe_class}/MMult
root_data_dir=../data/${exe_class}/MMult

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=${mmult_dir}/mmult-${exe_class}-${version}-way.o
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

		# create data file input
		# random initialization for two input matrices
		echo N > input.txt
		let lastIndex=$input_size-1
		echo 0:$lastIndex*0:$lastIndex >> input.txt
		echo N >> input.txt
		echo 0:$lastIndex*0:$lastIndex >> input.txt
		# three partition parameters
		echo $cache_block >> input.txt
		echo $cache_block >> input.txt
		echo $cache_block >> input.txt
		# disable storing of output and input matrices
		echo N >> input.txt
		echo N >> input.txt
		echo N >> input.txt

		# execute the program using the input file piping command line inputs
		../executable.o < input.txt > output.txt
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
