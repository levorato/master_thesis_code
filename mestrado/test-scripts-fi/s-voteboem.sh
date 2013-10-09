#INICIO DO SCRIPT
#o script precisa entrar na pasta onde o programa esta
export LD_LIBRARY_PATH=$HOME/lib:/opt/intel/mpi/3.2.2.006/lib64:$LD_LIBRARY_PATH
cd /home_nfs/mlevorato/graspcc-1.0/build
#agora a primeira instancia sera executada com n=100
#para que a o script nao fique bloqueado esperando a instancia
#terminar antes de chamar a seguinte, o comando e executado em
#background com o caracter & no fim do comando
mpdboot -n 1
mpiexec -n 1 ./graspcc -l 1 --iter=-1 --alpha=1.0 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Social Media/slashdot-undirected/part0" --output-folder "./output/sdot-voteboem" --gain-function-type=0 --firstImprovementOnOneNeig=false --time-limit=3600
mpdallexit
#FIM DO SCRIPT

