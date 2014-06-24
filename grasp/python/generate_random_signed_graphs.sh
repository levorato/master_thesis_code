array1=("100" "150" "200" "300" "400" "600")
array2=("0.1" "0.2" "0.5" "0.8")
array3=("0.2" "0.5" "0.8")
for n in "${array1[@]}"; do   # The quotes are necessary here
  for i in "${array2[@]}"; do
    for j in "${array3[@]}"; do 
        cmd="python generate_random_signed_graph.py -n $n -d $i -nd $j"
        echo $cmd
        eval "$cmd"
    done
  done
done


