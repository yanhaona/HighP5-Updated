#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# parallel run counts and input size configuration
cell_length=100
grid_dim=$1
points_per_cell=$2
parallel_runs=$3

exec_type=pinned
exec_sub_dir=pinned-versions

parallelism=( 4 2 1 )

# make data directory for experiments
mkdir -p ../data/Monte/highp5-low
root_data_dir=../data/Monte/highp5-low

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=../executables/$exec_sub_dir/monte-simple-toza-$exec_type-${version}-way.o
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

                # execute the program with proper parameters
                (time mpirun -n 1 ../executable.o \
                        cell_length=$cell_length grid_dim=$grid_dim \
                        points_per_cell=$points_per_cell b=$version) > output.txt 2>&1
		echo "Grid dim $grid_dim by $grid_dim" >> output.txt
		echo "Samples per grid cell $points_per_cell" >> output.txt
		echo "Grid cell dimension $cell_length by $cell_length" >> output.txt

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
