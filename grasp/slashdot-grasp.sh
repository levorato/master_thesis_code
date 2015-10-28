export LD_LIBRARY_PATH=/home/levorato/lib:$LD_LIBRARY_PATH
cd build

for run in $(seq 30); do
        cmd='./graspcc -l 1 --iter=400 --alpha=1.0 --input-file-dir "/home/levorato/mario/git/cuda/mestrado/grasp/tests/Instances/Social Media/slashdot-undirected/part0" --output-folder "../output/cuda-grasp-sbpo/sdot-seq-l1-i400-1h" --gain-function-type=0 --time-limit=3600 --rcc=false --strategy GRASP'
        echo $cmd
        eval "$cmd"
done

