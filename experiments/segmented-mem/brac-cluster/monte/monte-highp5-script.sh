#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# parallel run counts and input size configuration
cell_length=100
grid_dim=$1
points_per_cell=$2
parallel_runs=$3
type=$4
threads=6

# then check whether proper number of parameters have been provided
if [[ $# < 4 ]]; then
        echo "You have to provide three cmd arguments to run the experiments:"
        echo "  First, the dimension of the square grid within which Monte Carlo area estimation to be done"
        echo "  Second, the number of random samples per grid cell"
        echo "  Third, the number of times each program version should run"
        echo "  Fourth, type of threads in MPI processes 'random' or 'pinned'"
        exit -1
fi

parallelism=( 8 7 6 5 4 3 2 1 )

# make data directory for experiments
mkdir -p ./data/Monte/highP5
root_data_dir=./data/Monte/HighP5

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

echo "In the experiments, threads to CPU core assignment is $type"
echo ""

for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=./executables/${type}/monte-simple-${version}Nodes-${threads}Threads.o
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

		let partitions=$version*$threads
		
		# execute the program with proper parameters
                (time mpirun -n $version --oversubscribe ../executable.o \
                        cell_length=$cell_length grid_dim=$grid_dim \
                        points_per_cell=$points_per_cell b=$partitions) > output.txt 2>&1
                echo "Grid dim $grid_dim by $grid_dim" >> output.txt
                echo "Samples per grid cell $points_per_cell" >> output.txt
                echo "Grid cell dimension $cell_length by $cell_length" >> output.txt	

		cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done
