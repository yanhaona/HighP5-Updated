!#/bin/bash

curr_dir=`pwd`
cd scripts
chmod a+x *
cd $curr_dir
cd executables
chmod a+x *
cd $curr_dir
cd tools
chmod a+x *
cd $curr_dir
matrix_size=$1
repeat_count=$2
export matrix_size=$matrix_size
export repeat_count=$repeat_count

echo "Going to run the HighP5 program versions"
cd scripts
./mmult-highp5-script.sh $matrix_size $repeat_count
echo "HighP5 experiment runs are complete"
echo ""
echo ""
echo ""
echo "Going to run the Pthreads program versions"
./mmult-pthread-script.sh $matrix_size $repeat_count
echo "Pthread experiment runs are complete"


