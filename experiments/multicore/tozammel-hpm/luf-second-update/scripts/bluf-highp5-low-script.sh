#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# parallel run counts and input size configuration
input_size=$1
parallel_runs=$2
cache_block=64
exec_type=pinned
exec_sub_dir=pinned-versions

parallelism=( 4 2 1 )

# make data directory for experiments
mkdir -p ../data/Bluf/highp5-low
root_data_dir=../data/Bluf/highp5-low

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=../executables/$exec_sub_dir/bluf-toza-$exec_type-${version}-way.o
	echo "--------------------------------------------------------------------------------------------------------"
	echo "running $version-way parallel experiments"
	
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
		# create experiment sub-directory
		mkdir run-${run}
		# go inside that directory
		cd run-${run}

		# create data file input
                # random initialization for the input matrix
                echo N > input.txt
                let lastIndex=$input_size-1
                echo 0:$lastIndex*0:$lastIndex >> input.txt
                # provide the block partition parameter of block-stride partition
                echo $cache_block >> input.txt
                # disable storing of output and input matrices
                echo N >> input.txt
                echo N >> input.txt
                echo N >> input.txt
                # disable storing of the pivot array
                echo N >> input.txt

                # execute the program using the input file piping command line inputs
               	(time ../executable.o < input.txt) > output.txt 2>&1
		echo "matrix dimension $input_size by $input_size" >> output.txt
		
		# check if experiment run correctly
		run_checker=`cat output.txt | grep 'Execution Time' | wc -l`
		if [[ $run_checker == 1 ]]; then
			echo "experiment number $run of version $version parallel code was successfull"
			execution_time=`cat output.txt | grep 'Execution Time'`
                        echo $execution_time
		else
			echo "experiment number $run of version $version parallel code failed"
		fi	

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
