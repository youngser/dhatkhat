#!/bin/sh

for p in `seq 1 6`
do
	for i in `seq 1 100`
	do
		echo $p $i
		sbatch run_single_simulation.scr $p $i
	done 
done