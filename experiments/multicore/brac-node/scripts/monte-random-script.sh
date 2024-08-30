#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# adjust these two variables based on the mode of the experiment
exe_class=random
mode=../executables/random-affinity

# parallel run counts
parallel_runs=2

# input parameter configuration
precision=0.000001
cell_length=1000
points_per_cell=10
grid_dim=1024
max_rounds=10

monte_dir=${mode}/Monte
parallelism=( 20 12 8 6 2 1 )

# make data directory for experiments
mkdir -p ../data/${exe_class}/Monte
root_data_dir=../data/${exe_class}/Monte

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=${monte_dir}/monte-${exe_class}-${version}-way.o
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

		# execute the program with proper parameters
		(time mpirun -n 1 ../executable.o \
			precision=$precision max_rounds=$max_rounds \
			cell_length=$cell_length grid_dim=$grid_dim \
			points_per_cell=$points_per_cell) > output.txt 2>&1
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
