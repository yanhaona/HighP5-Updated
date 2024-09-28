#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# parallel run counts and input size configuration
plate_rows=$1
plate_cols=$2
iterations=$3
padding=$4
parallel_runs=$5
type=$6
threads=6

# then check whether proper number of parameters have been provided
if [[ $# < 6 ]]; then
        echo "You have to provide three cmd arguments to run the experiments:"
        echo "	First, the number of rows in the plate to do heat propagation"
        echo "	Second, the number of columns in the plate"
	echo "	Third, the number of iterations of heat propagation calculation"
        echo "	Fourth, the amound of padding rows between consecutive MPI prococesses"
        echo "	Fifth, the number of times each program version should run"
        echo "	Sixth, type of threads in MPI processes 'random' or 'pinned'"
        exit -1
fi

parallelism=( 8 7 6 5 4 3 2 1 )

# make data directory for experiments
mkdir -p ./data/Stencil/highP5
root_data_dir=./data/Stencil/HighP5

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

echo "In the experiments, threads to CPU core assignment is $type"
echo ""

echo "creating input data files"
# array generator executable
array_generator=../../../../tools/binary-array-generator
# generate the plate object in the root data directory
$array_generator 3 2 ${root_data_dir}/plate $plate_rows $plate_cols

for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=./executables/${type}/stencil-${version}Nodes-${threads}Threads.o
	echo "--------------------------------------------------------------------------------------------------------"
	echo "running $version-nodes and $threads per node parallel experiments"
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

		let upper_block=$version
                let lower_block=$threads
                let refinements=$padding*$iterations

		# execute the program with proper command line input
                (time mpirun -n $version --oversubscribe ../executable.o input_file=../../plate \
                        k=$upper_block l=1 m=$lower_block n=1 \
                        p1=$padding p2=1 \
                        iterations=$iterations) > output.txt 2>&1

                # add some information about the input
                echo "plate dimension is $updated_rows by $plate_cols" >> output.txt
                echo "refinement iterations: $refinements" >> output.txt
                echo "Upper blocking: $upper_block-by-1" >> output.txt
                echo "Lower blocking: $lower_block-by-1" >> output.txt
                echo "Padding between MPI processes: $padding" >> output.txt
                echo "Padding between threads inside an MPI process: 1" >> output.txt

		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done

# delete the bulky input plate file
rm ${root_data_dir}/plate
