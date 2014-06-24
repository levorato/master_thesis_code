python generate_random_graph_chinese.py -c 4 -n 16 -k 16 -pin 0.2 -p_minus 0 -p_plus 0
python generate_random_graph_chinese.py -c 4 -n 16 -k 16 -pin 0.7 -p_minus 0 -p_plus 0
python generate_random_graph_chinese.py -c 4 -n 32 -k 20 -pin 0.7 -p_minus 0.6 -p_plus 0.6
python generate_random_graph_chinese.py -c 4 -n 32 -k 24 -pin 0.7 -p_minus 0.2 -p_plus 0.2
python generate_random_graph_chinese.py -c 4 -n 64 -k 24 -pin 0.7 -p_minus 0.2 -p_plus 0.2
python generate_random_graph_chinese.py -c 25 -n 30 -k 20 -pin 0.6 -p_minus 0.3 -p_plus 0.3
python generate_random_graph_chinese.py -c 4 -n 64 -k 32 -pin 0.8 -p_minus 0 -p_plus 0.1
python generate_random_graph_chinese.py -c 4 -n 96 -k 24 -pin 0.7 -p_minus 0.2 -p_plus 0.2


array=("0.0" "0.2" "0.4" "0.6" "0.8" "1.0")
for i in "${array[@]}"; do   # The quotes are necessary here
    for j in "${array[@]}"; do 
        cmd="python generate_random_graph_chinese.py -c 4 -n 64 -k 32 -pin 0.8 -p_minus $i -p_plus $j"
        echo $cmd
        eval "$cmd"
    done
done


