#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

machine=tozammel
parallelism=( 4 2 1 )
executable=../executables/pthread-monte-simple.o

# parallel run counts and input size configuration
cell_length=100
grid_dim=$1
points_per_cell=$2
parallel_runs=$3


# make data directory for experiments
mkdir -p ../data/Monte/pthread-low
root_data_dir=../data/Monte/pthread-low

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
                (time ../../pthread_exec.o $cell_length $grid_dim $points_per_cell $version $machine) > output.txt 2>&1
                echo "cell length: $cell_length" >> output.txt
                echo "grid dim: $grid_dim by $grid_dim" >> output.txt
                echo "samples per cell: $points_per_cell" >> output.txt

		# check if experiment run correctly
                run_checker=`cat output.txt | grep 'Execution Time' | wc -l`
                if [[ $run_checker == 1 ]]; then
                        echo "experiment number $run of version $version-way parallel code was successfull"
			execution_time=`cat output.txt | grep 'Execution Time'`
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
