#!/bin/bash

for i in {1..4..1}
do 
	j=$( echo "scale = 1; $i/2" | bc)
	mkdir -p k-$j/
	cd k-$j/
	cp ../* .
	sed -i "s/input-coupling/$j/g" potential_mD.f90
	make 
        sbatch run_pool.sh
	cd ../
done
