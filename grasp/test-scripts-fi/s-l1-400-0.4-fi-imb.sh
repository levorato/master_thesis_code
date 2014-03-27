#INICIO DO SCRIPT
#o script precisa entrar na pasta onde o programa esta
export LD_LIBRARY_PATH=$HOME/lib:/opt/intel/mpi/3.2.2.006/lib64:$LD_LIBRARY_PATH
cd /home_nfs/mlevorato/graspcc-1.0/build
mpdboot -n 1
mpiexec -n 1 ./graspcc -l 1 --iter=400 --alpha=0.4 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Social Media/slashdot-undirected/part0" --output-folder "./output/sdot-l1-i400-0.4-fi-imb" --gain-function-type=0 --time-limit=3600
mpdallexit
#FIM DO SCRIPT

