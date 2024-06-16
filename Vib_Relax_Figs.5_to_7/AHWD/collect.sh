#!/bin/bash

for t in 1 100 300
do 
	for i in {5..25..10} 
	do 
		j=$( echo "scale = 2; $i/100" | bc)
		cp T-$t/eta-$j/Nb-20/ntraj-10-8/fort.192 collected_data/HWD_T_"$t"_eta_"$j"_t192.out
		cp T-$t/eta-$j/Nb-20/ntraj-10-8/fort.640 collected_data/HWD_T_"$t"_eta_"$j"_t640.out
		cp T-$t/eta-$j/Nb-20/ntraj-10-8/fort.1600 collected_data/HWD_T_"$t"_eta_"$j"_t1600.out
	done
done
