#!/bin/bash
export LD_LIBRARY_PATH=/home/levorato/lib:$LD_LIBRARY_PATH
cd build
alpha_array=("0.0" "0.4" "0.8" "1.0")
iter_array=("100" "200" "300" "400")
r_array=("1" "2")
for i in "${iter_array[@]}"; do   # The quotes are necessary here
	for j in "${alpha_array[@]}"; do 
		for r in "${r_array[@]}"; do
			#for x in {1..5}; do			# five repetitions of each combination
				cmd="./graspcc -l $r --iter=$i --alpha=$j --input-file-dir '../tests/Instances/UNGA/UNGA_Weighted_SignedGraphs/paramet' --output-folder '../output/parametrize/unga-grasp/l$r-a$j-i$i' --gain-function-type=0 --time-limit 7200 --strategy GRASP --rcc false --jobid=parametric-unga-grasp-l$r-a$j-i$i"
				echo $cmd
				eval "$cmd"
			#done
		done
	done
done
