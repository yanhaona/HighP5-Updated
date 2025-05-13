#!/bin/bash

# first check whether the autorun script is called from the proper directory
bluf_tracker_found=`cat bluf-dir-track | wc -l`
if [[ $bluf_tracker_found == 0 ]]; then
	echo "You must run the autorun.sh script from the directory containing that file"
	exit -1
fi

# then check whether proper number of parameters have been provided
if [[ $# < 4 ]]; then
	echo "You have to provide three cmd arguments to run the experiments:"
	echo "	The dimension of the square matrix for doing LU factorization"
	echo "	The dimension of the square sparse matrix for doing conjugate gradient "
	echo "	The number of time each parallel version should repeat"
	echo "	Threading type 'random' or 'pinned'"
	exit -1
fi

# try to remove the data directory from the previous experiment
rm -rf data

# check if there is already a data directory then ask the user to remove it
if test -d data; then
	echo "It seems there is an existing data directory in this folder."
	echo "If you need the data that exist there then save the directory content elsewhere."
	echo "Then delete the data directory by running rm -rf data"
	echo "Then rerun this script"
	exit -1
fi

luf_matrix_size=$1
cg_matrix_size=$2
program_repeats=$3
thread_type=$4

# make the necessary scripts and binaries executable
curr_dir=`pwd`
chmod a+x luf-run.sh
chmod a+x conj-run.sh

echo "----------------------------------------------------------------------------------------"
echo "LU Decomposition"
echo "----------------------------------------------------------------------------------------"
./luf-run.sh $luf_matrix_size $program_repeats $thread_type
echo ""
echo ""
echo ""
echo "----------------------------------------------------------------------------------------"
echo "Conjugate Gradient"
echo "----------------------------------------------------------------------------------------"
./conj-run.sh $cg_matrix_size $program_repeats $thread_type



