#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

machine=tozammel
parallelism=( 4 2 1 )
executable=../executables/pthread-stencil.o

# parallel run counts and input size configuration
plate_dim=$1
parallel_runs=$2
iterations=500
padding=4


# make data directory for experiments
mkdir -p ../data/Stencil/pthread-low
root_data_dir=../data/Stencil/pthread-low

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# copy the executable in the root data directory
cp $executable $root_data_dir/pthread_exec.o

# array generator executable
array_generator=../tools/binary-array-generator

# generate the plate object in the root data directory
$array_generator 3 2 ${root_data_dir}/plate $plate_dim $plate_dim

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
                (time ../../pthread_exec.o ../../plate $iterations $padding $version $machine) > output.txt 2>&1
                let refinements=$iterations*$padding
                echo "Plate Dimension $plate_rows by $plate_cols" >> output.txt
                echo "Refinement iterations $refinements" >> output.txt
                echo "Padding rows $padding" >> output.txt

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

# remove the bulky plate file from data directory
rm -f  ${root_data_dir}/plate
