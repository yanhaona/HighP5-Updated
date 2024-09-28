#!/bin/bash

# keep track of the current directory
curr_dir=`pwd`

exec_type=pinned
exec_sub_dir=pinned-versions

# parallel run counts and input size configuration
matrix_size=$1
iterations=$2
parallel_runs=$3
sparsity=90

parallelism=( 128 64 32 16 8 )

# make data directory for experiments
mkdir -p ../data/CG/highp5-high
root_data_dir=../data/CG/highp5-high

# clean previous content of the data directory if exists
rm -rf ${root_data_dir}/*


# tools references
array_generator=../tools/array-generator
matrix_generator=../tools/sparse-matrix-generator

echo "generating input data"
# generate sparse matrix and move associated files in the root data directory
echo "values" > input.txt
echo "columns" >> input.txt
echo "rows" >> input.txt
$matrix_generator $matrix_size $matrix_size $sparsity 3 0 < input.txt > /dev/null
cp values ${root_data_dir}/values
cp columns ${root_data_dir}/columns
cp rows ${root_data_dir}/rows
rm values columns rows input.txt

# generate predication and known vectors
$array_generator $matrix_size double 0.25 100 ${root_data_dir}/known > /dev/null
$array_generator $matrix_size double 0.25 100 ${root_data_dir}/pred > /dev/null

echo "Experiment Description:"
echo "	sparse matrix of $matrix_size by $matrix_size size"
echo "	matrix sparsity $sparsity"
echo "	gradient descent iterations $iterations"
for version in "${parallelism[@]}"
do
	# find the location of the executable file
        executable=../executables/$exec_sub_dir/cg-$exec_type-${version}-way.o
	echo "--------------------------------------------------------------------------------------------------------"
	echo "running pthread $version-way parallel experiments"
	
	# create a subdirectory for the current degree of parallelism
	mkdir ${root_data_dir}/version-${version}
	exper_dir=${root_data_dir}/version-${version}

	# copy the executable in that directory
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

		# generate an input file for the program
		echo 'Y' > input.txt
		echo 'Y' >> input.txt
		echo 'Y' >> input.txt
		echo 'Y' >> input.txt
		echo 'Y' >> input.txt

		( time ../executable.o \
                        arg_matrix_cols=../../columns arg_matrix_rows=../../rows \
                        arg_matrix_values=../../values known_vector=../../known \
                        prediction_vector=../../pred r=$version \
                        iterations=$iterations < input.txt ) > output.txt 2>&1
                echo "matrix of $matrix_size by $matrix_size" >> output.txt
                echo "matrix sparsity $sparsity" >> output.txt
                echo "gradient descent iterations $iterations" >> output.txt
		rm input.txt

		# check if experiment run correctly
                run_checker=`cat output.txt | grep 'Execution Time' | wc -l`
                if [[ $run_checker == 1 ]]; then
                        echo "experiment number $run of version $version-way parallel code was successfull"
			execution_time=`cat output.txt | grep 'Execution Time'`
			mem_init_time=`cat it-program.log | grep 'Memory Preparation Time'`
			echo $mem_init_time
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
cd ${root_data_dir}
rm columns values rows known pred

cd $curr_dir
