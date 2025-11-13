#!/bin/bash

START=0
END=47
fileName=fugaku.cn
touch $fileName
truncate -s 0 $fileName

for i in $(seq $START $END); do
	echo "processor	: $i	physical id	: 0	core id		: $i" >> $fileName
done

