#!/bin/bash

for t in 1 100 300
do 
	for i in {5..25..10} 
	do 
		j=$( echo "scale = 2; $i/100" | bc)
		tail -2 T-$t/eta-$j/Nb-20/ntraj-10-7/VR_HWD.o* >> broken_info.out  
	done
done
