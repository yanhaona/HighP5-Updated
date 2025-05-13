#!/bin/bash

# first check whether the autorun script is called from the proper directory
cg_tracker_found=`cat cg-dir-track | wc -l`
if [[ $cg_tracker_found == 0 ]]; then
	echo "You must run the autorun.sh script from the directory containing that file"
	exit -1
fi

# then check whether proper number of parameters have been provided
if [[ $# < 3 ]]; then
	echo "You have to provide three cmd arguments to run the experiments:"
	echo "	The dimension of the sparse matrix for conjugate gradient calculation"
	echo "	The number of time each high order parallel version should repeat"
	echo "	Threading type 'random' or 'pinned'"
	exit -1
fi

# make the necessary scripts and binaries executable
curr_dir=`pwd`
cd scripts
chmod a+x *
cd $curr_dir
cd executables
chmod a+x *
cd conj/pinned-versions
chmod a+x *
cd ../random-versions
chmod a+x *
cd $curr_dir
cd tools
chmod a+x *
cd $curr_dir

matrix_size=$1
high_repeat_count=$2
thread_type=$3
iterations=1000
export matrix_size=$matrix_size
export high_repeat_count=$high_repeat_count
export thread_type=$thread_type
export iterations=$iterations

echo ""
echo ""
echo "HighP5 128, 64, 32, 16 way parallel codes will repeat for $high_repeat_count times"
echo ""
echo ""

# go to the scripts directory
cd scripts

echo "Going to run the HighP5 program versions"
./cg-highp5-high-script.sh $matrix_size $iterations $high_repeat_count $thread_type
