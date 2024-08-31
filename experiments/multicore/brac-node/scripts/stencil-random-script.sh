#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# adjust these two variables based on the mode of the experiment
exe_class=random
mode=../executables/random-affinity

# parallel run counts and problem size configurations
parallel_runs=2
input_size=1024
iterations=1000
upper_padding=4
lower_padding=1

# array generator executable
array_generator=../../../../tools/array-generator

stencil_dir=${mode}/Stencil
parallelism=( 20 12 8 6 2 1 )

# make data directory for experiments
mkdir -p ../data/${exe_class}/Stencil
root_data_dir=../data/${exe_class}/Stencil

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=${stencil_dir}/stencil-${exe_class}-${version}-way.o
	echo "--------------------------------------------------------------------------------------------------------"
	echo "running $version-way parallel experiments"
	echo "Executable $executable"
	
	# create a subdirectory for the current executable
	mkdir ${root_data_dir}/version-${version}
	exper_dir=${root_data_dir}/version-${version}

	# copy the executable in the created directory
	cp $executable $exper_dir/executable.o

	# create a plate object input file
	./$array_generator $input_size*$input_size double 0.25 100  $exper_dir/plate.txt

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

		let upper_block=1
		let lower_block=$version
		
		# execute the program with proper command line input
		(time ../executable.o input_file=../plate.txt output_file=plate-final.txt \
			k=$upper_block l=$upper_block m=$lower_block n=$lower_block \
			p1=$upper_padding p2=$lower_padding \
			iterations=$iterations) > output.txt 2>&1
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
