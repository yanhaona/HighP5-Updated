#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# adjust these two variables based on the mode of the experiment
ref_target=brac-node
machine=brac
parallelism=( 20 12 8 6 2 1 )
executable=../executables/pthread-stencil.o

# parallel runs count
parallel_runs=2

# sampling parameters
plate_rows=2048
plate_cols=512
iterations=1000
padding=1

# make data directory for experiments
mkdir -p ../../../$ref_target/data/pthread/Stencil
root_data_dir=../../../$ref_target/data/pthread/Stencil

# get the reference of binary array generator
array_generator=../../../../../tools/binary-array-generator

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# copy the executable in the root data directory
cp $executable $root_data_dir/pthread_exec.o

# generate the input plate in the root data directory
$array_generator 3 2 $root_data_dir/plate $plate_rows $plate_cols

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

		# execute the program using the input file piping command line inputs
		(time ../../pthread_exec.o ../../plate $iterations $padding $version $machine) > output.txt 2>&1
		let refinements=$iterations*$padding
		echo "Plate Dimension $plate_rows by $plate_cols" >> output.txt
		echo "Refinement iterations $refinements" >> output.txt
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
