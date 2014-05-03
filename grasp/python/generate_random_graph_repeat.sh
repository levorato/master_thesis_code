array=("0.0" "0.2" "0.4" "0.6" "0.8" "1.0")
for i in "${array[@]}"; do   # The quotes are necessary here
    for j in "${array[@]}"; do 
        cmd="python generate_random_graph.py -c 4 -n 64 -k 32 -pin 0.8 -p_minus $i -p_plus $j"
        echo $cmd
        eval "$cmd"
    done
done

