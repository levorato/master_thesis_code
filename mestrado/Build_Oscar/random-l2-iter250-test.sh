#Isso e um comentario, INICIO DO SCRIPT
#o script precisa entrar na pasta onde o programa esta
export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH
cd /home_nfs/mlevorato/mestrado/Build_Oscar
#agora a primeira instancia sera executada com n=100
#para que a o script nao fique bloqueado esperando a instancia
#terminar antes de chamar a seguinte, o comando e executado em
#background com o caracter & no fim do comando
mpdboot -n 1
mpiexec -n 1 ./MestradoMario -l 2 --iter=250 --alpha=0.6 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Random"  --output-folder "./output/random-l2-iter250-seq-1"
mpdallexit
#FIM DO SCRIPT

