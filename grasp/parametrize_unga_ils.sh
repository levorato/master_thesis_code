#!/bin/bash
export LD_LIBRARY_PATH=/home/levorato/lib:$LD_LIBRARY_PATH
cd build
alpha_array=("0.4" "1.0")
iter_array=("5" "10")
iterils_array=("2" "5")
pert_array=("3" "30")
r_array=("1" "2")
for i in "${iter_array[@]}"; do   # The quotes are necessary here
	for j in "${alpha_array[@]}"; do 
		for neig in "${r_array[@]}"; do
			for p in "${pert_array[@]}"; do			# five repetitions of each combination
				cmd="./graspcc -l $neig --iter=$i --alpha=$j --perturbationLevelMax=$p --input-file-dir '../tests/Instances/UNGA/UNGA_Weighted_SignedGraphs/paramet' --output-folder '../output/parametrize/unga-ils/l$neig-a$j-i$i-iterils5-p$p' --gain-function-type=0 --time-limit 7200 --strategy ILS --rcc false --jobid=parametric-unga-ils-l$neig-a$j-i$i-iterils5-p$p"
				echo $cmd
				eval "$cmd"
			done
		done
	done
done
