export LD_LIBRARY_PATH=/home/levorato/lib:$LD_LIBRARY_PATH
cd build

for run in $(seq 30); do
        cmd='./graspcc -l 1 --iter=10 --alpha=1.0 --input-file-dir "/home/levorato/mario/git/cuda/mestrado/grasp/tests/Instances/Random_SG_Bigger" --output-folder "../output/cuda-ils-sbpo/rbig-seq-l1-i10-1h" --gain-function-type=0 --time-limit=3600 --rcc=false --strategy ILS'
        echo $cmd
        eval "$cmd"
done
