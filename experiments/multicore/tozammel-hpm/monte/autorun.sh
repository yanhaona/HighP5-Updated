#!/bin/bash

# first check whether the autorun script is called from the proper directory
bluf_tracker_found=`cat monte-dir-track | wc -l`
if [[ $bluf_tracker_found == 0 ]]; then
	echo "You must run the autorun.sh script from the directory containing that file"
	exit -1
fi

# then check whether proper number of parameters have been provided
if [[ $# < 3 ]]; then
	echo "You have to provide three cmd arguments to run the experiments:"
	echo "	The number of Monte Carlo sampling experiments"
	echo "	The number of time each high order parallel version should repeat"
	echo "	The number of time ech low order parallel version should repeat"
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

# make the necessary scripts and binaries executable
curr_dir=`pwd`
cd scripts
chmod a+x *
cd $curr_dir
cd executables
chmod a+x *
cd pinned-versions
chmod a+x *
cd ../random-versions
chmod a+x *
cd $curr_dir

grid_dim=1024
monte_samples=$1
high_repeat_count=$2
low_repeat_count=$3
export grid_dim=$grid_dim
export monte_samples=$monte_samples
export high_repeat_count=$high_repeat_count
export low_repeat_count=$low_repeat_count

echo ""
echo ""
echo "HighP5 and handwritten Pthread 128, 64, 32, 16, 8 way parallel codes will repeat for $high_repeat_count times"
echo "HighP5 and handwritten Pthread 4, 2, 1 way parallel codes will repeat for $low_repeat_count times"
echo ""
echo ""

# go to the scripts directory
cd scripts

echo "Going to run the HighP5 program versions"
./monte-highp5-high-script.sh $grid_dim $monte_samples $high_repeat_count
./monte-highp5-low-script.sh $grid_dim $monte_samples $low_repeat_count
echo ""
echo ""
echo ""
echo "Going to run the Pthreads program versions"
./monte-pthread-high-script.sh $grid_dim $monte_samples $high_repeat_count
./monte-pthread-low-script.sh $grid_dim $monte_samples $low_repeat_count




