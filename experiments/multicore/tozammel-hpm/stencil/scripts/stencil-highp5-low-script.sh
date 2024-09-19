#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# parallel run counts and input size configuration
plate_dim=$1
parallel_runs=$2
iterations=500
upper_padding=4
lower_padding=4

exec_type=pinned
exec_sub_dir=pinned-versions

parallelism=( 4 2 1 )

# make data directory for experiments
mkdir -p ../data/Stencil/highp5-low
root_data_dir=../data/Stencil/highp5-low

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# array generator executable
array_generator=../tools/binary-array-generator

# generate the plate object in the root data directory
$array_generator 3 2 ${root_data_dir}/plate $plate_dim $plate_dim

for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=../executables/$exec_sub_dir/stencil-toza-$exec_type-${version}-way.o
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

		let upper_block=1
                let lower_block=$version
                let refinements=$lower_padding*$iterations

                # execute the program with proper command line input
                (time mpirun -n 1 ../executable.o input_file=../../plate \
                        k=$upper_block l=$upper_block m=$lower_block n=1 \
                        p1=$upper_padding p2=$lower_padding \
                        iterations=$iterations) > output.txt 2>&1

                # add some information about the input
                echo "plate dimension is $updated_rows by $plate_cols" >> output.txt
                echo "refinement iterations $refinements" >> output.txt
                echo "Upper blocking $upper_block-by-$upper_block" >> output.txt
                echo "Lower blocking $lower_block-by-1" >> output.txt
                echo "Upper padding $upper_padding" >> output.txt
                echo "Lower padding $lower_padding" >> output.txt

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

# remove the bulky plate file from data directory
rm -f  ${root_data_dir}/plate
