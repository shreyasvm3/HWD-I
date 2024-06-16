#!/bin/bash

for i in 2  
do 
     	j=$( echo "scale = 1; $i/2" | bc)
	mkdir -p k-$j/
	cd k-$j/
        for l in {7..8..1}
 	do
		mkdir -p ntraj-10-$l/
		cd ntraj-10-$l/
		cp ../../dyn.x ../../input_mD ../../run_multinode.sh .
		sed -i "s/input-coupling/$j/g" input_mD
		sed -i "s/input-ntraj/$l/g" input_mD
	        sbatch run_multinode.sh
		cd ../
	done
	cd ../
done
