#!/bin/sh

for i in `seq 1 114`
do
	echo $i
	sbatch run_single_realdata.scr $i
done