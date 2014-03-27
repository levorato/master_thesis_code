#INICIO DO SCRIPT
#o script precisa entrar na pasta onde o programa esta
export LD_LIBRARY_PATH=$HOME/lib:/opt/intel/mpi/3.2.2.006/lib64:$LD_LIBRARY_PATH
cd /home_nfs/mlevorato/graspcc-1.0/build
#agora a primeira instancia sera executada com n=100
#para que a o script nao fique bloqueado esperando a instancia
#terminar antes de chamar a seguinte, o comando e executado em
#background com o caracter & no fim do comando
mpdboot -n 1
mpiexec -n 1 ./graspcc -l 2 --iter=400 --alpha=0.8 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Random" --output-folder "./output/rnd-l2-i400-0.8-fII" --gain-function-type=4
mpdallexit
#FIM DO SCRIPT

