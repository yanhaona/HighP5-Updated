#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# parallel run counts and input size configuration
input_size=$1
parallel_runs=$2
type=$3
threads=$4
cache_block=64

# then check whether proper number of parameters have been provided
if [[ $# < 4 ]]; then
        echo "You have to provide three cmd arguments to run the experiments:"
        echo "  First, the dimension length of the square input matrices"
        echo "  Second, the number of times each program version should run"
        echo "  Third, type of threads in MPI processes 'random' or 'pinned'"
        echo "  Foruthm, number of threads per MPI process: options are 1 and 6"
        exit -1
fi

parallelism=( 8 7 6 5 4 3 2 1 )

# make data directory for experiments
mkdir -p ./data/Bluf/highP5
root_data_dir=./data/Bluf/HighP5

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

echo "In the experiments, threads to CPU core assignment is $type"
echo ""

echo "creating input data files"
# array generator executable
array_generator=../../../../tools/binary-array-generator
# generate the input matrix in the root data directory
$array_generator 3 2 ${root_data_dir}/a $input_size $input_size

for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=./executables/${type}/bluf-${version}Nodes-${threads}Threads.o
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

		# execute the program using the input file piping command line inputs
		( time mpirun -n $version --oversubscribe ../executable.o \
			input_matrix_file=../../a block_size=$cache_block ) > output.txt 2>&1
		echo "matrix dimension $input_size by $input_size" >> output.txt
		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done

# delete the bulky input matrix file
rm ${root_data_dir}/a
