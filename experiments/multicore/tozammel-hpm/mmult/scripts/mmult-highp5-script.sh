#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# parallel run counts and input size configuration
input_size=$1
parallel_runs=$2
cache_block=64

parallelism=( 128 64 32 16 8 4 2 1 )

# make data directory for experiments
mkdir -p ../data/MMult/highP5
root_data_dir=../data/MMult/HighP5

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=../executables/mmult-toza-${version}-way.o
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
		echo "matrix dimension $input_size by $input_size" >> output.txt
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
