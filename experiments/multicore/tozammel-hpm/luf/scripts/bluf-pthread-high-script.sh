#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

machine=tozammel
parallelism=( 128 64 32 16 8 )
executable=../executables/pthread-bluf.o

# parallel run counts and input size configuration
input_size=$1
parallel_runs=$2
cache_block=64

# make data directory for experiments
mkdir -p ../data/Bluf/pthread-high
root_data_dir=../data/Bluf/pthread-high

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# copy the executable in the root data directory
cp $executable $root_data_dir/pthread_exec.o

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

		# execute the program using the input file piping command line inputs
                (time ../../pthread_exec.o $cache_block $version $machine $input_size) > output.txt 2>&1

		echo "Matrix dimension $input_size by $input_size" >> output.txt
		echo "$version way parallel pthread version" >> output.txt

		# check if experiment run correctly
                run_checker=`cat output.txt | grep 'total time' | wc -l`
                if [[ $run_checker == 1 ]]; then
                        echo "experiment number $run of version $version-way parallel code was successfull"
			execution_time=`cat output.txt | grep 'computation time'`
                        echo $execution_time
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
