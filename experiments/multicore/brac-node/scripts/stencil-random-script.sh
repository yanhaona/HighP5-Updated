#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# adjust these two variables based on the mode of the experiment
exe_class=random
mode=../executables/random-affinity

# parallel run counts and problem size configurations
parallel_runs=1
plate_rows=4096
row_column_ratio=1
iterations=100
upper_padding=1
lower_padding=1

# array generator executable
array_generator=../../../../tools/binary-array-generator

stencil_dir=${mode}/Stencil
parallelism=( 8 )

# make data directory for experiments
mkdir -p ../data/${exe_class}/Stencil
root_data_dir=../data/${exe_class}/Stencil

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# generate the plate object in the root data directory
let plate_cols=$plate_rows/$row_column_ratio 
let updated_rows=$plate_rows*$row_column_ratio
$array_generator 3 2 ${root_data_dir}/plate $updated_rows $plate_cols

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
		let refinements=$lower_padding*$iterations
		
		# execute the program with proper command line input
		(time ../executable.o input_file=../../plate \
			k=$upper_block l=$upper_block m=$lower_block n=1 \
			p1=$upper_padding p2=$lower_padding \
			iterations=$iterations) > output.txt 2>&1

		# add some information about the input
		echo "plate dimension is $updated_rows by $plate_cols" >> output.txt
		echo "refinement iterations $refinements" >> output.txt

		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
