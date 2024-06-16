#!/bin/bash

for t in 100 
do 
	for i in {5..5..10} 
	do
		j=$( echo "scale = 2; $i/100" | bc)
		for l in {2..8..1}
		do 
			cp T-$t/eta-$j/Nb-20/ntraj-10-$l/fort.192 collected_data/HWD_nt"$l"_t192.out
			cp T-$t/eta-$j/Nb-20/ntraj-10-$l/fort.640 collected_data/HWD_nt"$l"_t640.out
			cp T-$t/eta-$j/Nb-20/ntraj-10-$l/fort.1600 collected_data/HWD_nt"$l"_t1600.out
		done
	done
done
