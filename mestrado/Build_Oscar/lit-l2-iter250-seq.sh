#Isso e um comentario, INICIO DO SCRIPT
#o script precisa entrar na pasta onde o programa esta
export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH
cd /home_nfs/mlevorato/mestrado/Build_Oscar
#agora a primeira instancia sera executada com n=100
#para que a o script nao fique bloqueado esperando a instancia
#terminar antes de chamar a seguinte, o comando e executado em
#background com o caracter & no fim do comando
mpdboot -n 1
mpiexec -n 1 ./MestradoMario -l 2 --iter=250 --alpha=0.6 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Literatura/Formato\ PAJEK" --output-folder "./output/lit-l2-iter250-seq-1" &
mpiexec -n 1 ./MestradoMario -l 2 --iter=250 --alpha=0.6 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Literatura/Formato\ PAJEK" --output-folder "./output/lit-l2-iter250-seq-2" &
mpiexec -n 1 ./MestradoMario -l 2 --iter=250 --alpha=0.6 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Literatura/Formato\ PAJEK" --output-folder "./output/lit-l2-iter250-seq-3" &
mpiexec -n 1 ./MestradoMario -l 2 --iter=250 --alpha=0.6 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Literatura/Formato\ PAJEK" --output-folder "./output/lit-l2-iter250-seq-4" &
mpiexec -n 1 ./MestradoMario -l 2 --iter=250 --alpha=0.6 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Literatura/Formato\ PAJEK" --output-folder "./output/lit-l2-iter250-seq-5" &
#simultaneamente, colocamos a segunda instancia para rodar
#mpiexec -n 1 ./MestradoMario -l 2 --iter=250 --alpha=0.6 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Random"  --output-folder "./output/random-l2-iter250-seq-2" &
#agora a terceira instancia
#mpiexec -n 1 ./MestradoMario -l 2 --iter=250 --alpha=0.6 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Random"  --output-folder "./output/random-l2-iter250-seq-3" &
#agora a quarta instancia
#mpiexec -n 1 ./MestradoMario -l 2 --iter=250 --alpha=0.6 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Random"  --output-folder "./output/random-l2-iter250-seq-4" &
#e finalmente a ultima
#mpiexec -n 1 ./MestradoMario -l 2 --iter=250 --alpha=0.6 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Random"  --output-folder "./output/random-l2-iter250-seq-5" &
#e preciso fazer o script esperar pelo termino de todas as instancias com o comando
#wait. Se o script terminar e deixar alguma instancia rodando em background
#para tras o PBS pode matar a instancia antes dela terminar.
wait
mpdallexit
#FIM DO SCRIPT

