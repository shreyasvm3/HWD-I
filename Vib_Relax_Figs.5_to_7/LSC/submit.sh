#!/bin/bash

for t in 100
do 
	mkdir -p T-$t/
	cd T-$t/
	for i in {5..5..10} 
	do 
		j=$( echo "scale = 2; $i/100" | bc)
        	#k=$((15000/$i))
		mkdir -p eta-$j/
		cd eta-$j/
        	for k in {20..20..5}
		do 
			mkdir -p Nb-$k/
			cd Nb-$k/
			for l in {2..4..1}
			do 
				mkdir -p ntraj-10-$l/
				cd ntraj-10-$l/
		  		cp ../../../../dyn.x ../../../../run_multinode.sh ../../../../input_mD .
				sed -i "s/input-t/$t/g" input_mD
				sed -i "s/input-eta/$j/g" input_mD
				sed -i "s/input-nb/$k/g" input_mD
				sed -i "s/input-ntraj/$l/g" input_mD
		        	sbatch run_multinode.sh
				cd ../
			done
	        	cd ../
		done
		cd ../
	done
	cd ../
done
