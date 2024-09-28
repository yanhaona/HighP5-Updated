#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

# parallel run counts and input size configuration
matrix_size=$1
sparsity=90
maxIterations=$2
parallel_runs=$3
type=$4
threads=1
precision=0.0001

# then check whether proper number of parameters have been provided
if [[ $# < 4 ]]; then
        echo "You have to provide three cmd arguments to run the experiments:"
	echo "  First, the dimension length of the sparse matrix (default sparsity is 90%)"
        echo "  Second, the number of conjugate gradient iterations"
        echo "  Third, the number of times each experiment should run"
        echo "  Fourth, type of threads in MPI processes 'random' or 'pinned'"
        exit -1
fi

parallelism=( 8 7 6 5 4 3 2 1 )

# make data directory for experiments
mkdir -p ./data/CG/highP5
root_data_dir=./data/CG/HighP5

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*

# tool executables
matrix_generator=../../../../tools/sparse-matrix-generator
array_generator=../../../../tools/binary-array-generator

echo "In the experiments, threads to CPU core assignment is $type"
echo ""

echo "creating input data files"
# generate sparse matrix and move associated files in the root data directory
$matrix_generator $matrix_size $matrix_size $sparsity 3 1
cp values ${root_data_dir}/values
cp columns ${root_data_dir}/columns
cp rows ${root_data_dir}/rows
rm values columns rows

# generate predication and known vectors
$array_generator 3 1 ${root_data_dir}/known $matrix_size
$array_generator 3 1 ${root_data_dir}/pred $matrix_size

for version in "${parallelism[@]}"
do
	# find the location of the executable file
	executable=./executables/${type}/cg-multitask-${version}Nodes-${threads}Threads.o
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

                # execute the program with proper command line input
                ( time mpirun -n $version --oversubscribe ../executable.o \
                        arg_matrix_cols=../../columns arg_matrix_rows=../../rows \
                        arg_matrix_values=../../values known_vector=../../known \
                        prediction_vector=../../pred r=$version b=$version \
                        maxIterations=$maxIterations precision=$precision ) > output.txt 2>&1
                echo "matrix of $matrix_size by $matrix_size" >> output.txt
                echo "matrix sparsity $sparsity" >> output.txt
                echo "gradient descent iterations $maxIterations" >> output.txt
                cat output.txt

		# come back to the parent directory
		cd ..
	done	
	
	echo ""

	# come back to the current directory
	cd $curr_dir
done

# delete the bulky files
cd ${root_data_dir}
rm values columns rows known pred
